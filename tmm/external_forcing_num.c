/*
 * NUM coupled to the Transport Matrix Method driver.
 *
 * Spike scope: generalists-only, which has no sinking species, so the sparse
 * upwind sinking operator that simulateGlobal builds in matlab is not needed
 * yet. Temperature and light are the two forcings NUM needs; light is derived
 * from the TMM insolation routine with exponential attenuation, matching
 * p.kw in parametersGlobal.
 *
 * NUM supplies tendencies via f_calcderivatives, so PETSc owns the time
 * stepping. That is deliberate: it is what makes the Newton-Krylov spinup
 * usable, which an integrate-in-place coupling would give up.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "petscmat.h"
#include "tmm_petsc_matvec_utils.h"
#include "tmm_share.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_timer.h"
#include "tmm.h"

/* NUM library entry points (Fortran, bind(c)) */
extern void f_setupgeneralistsonly(const int n, int *errorio, char *errorstr);
extern void f_getnumgrid(int *n, int *idxB);
extern void f_calcderivatives(const double u[], const double L, const double T,
                              const double dt, double dudt[]);

extern void insolation_(PetscInt *N, PetscScalar *myTime, PetscScalar *locLatitude,
                        PetscScalar *daysPerYear, PetscScalar *locSWRad, PetscScalar *tau);

typedef struct {
  PetscScalar **localTR;   /* tracer state, [tracer][local point]  */
  PetscScalar **localJTR;  /* tendency TMM adds, same layout       */
  PetscScalar *localTs;    /* temperature, per local point         */
  PetscScalar *locallatitude, *localswrad, *localtau;
  PetscScalar *localdz;    /* layer thickness, per local point     */
  PetscScalar *u, *dudt;   /* one cell's state vector              */
  PetscInt nGrid;
} NUMContext;

static NUMContext numctx;

static PetscScalar kw = 0.1;          /* light attenuation, m^-1 */
static PetscScalar PARfrac = 0.4;     /* fraction of shortwave available as PAR */
static PetscScalar EinConv = 4.57;    /* W m^-2 -> umol photons s^-1 m^-2 */
static PetscScalar daysPerYear = 365.0;
static PetscScalar Tconst = 10.0;     /* spike: uniform temperature until Theta forcing is wired */
static PetscInt nSizeGroups = 10;

PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx)
{
  PetscErrorCode ierr;
  PetscInt numTracers, ip, kl, nzloc, k;
  int nGrid, idxB, errorio = 0;
  char errorstr[80];

  numTracers = state->numTracers;

  f_setupgeneralistsonly(nSizeGroups, &errorio, errorstr);
  if (errorio) SETERRQ(PETSC_COMM_WORLD, 1, "NUM setup failed: %s", errorstr);

  f_getnumgrid(&nGrid, &idxB);
  numctx.nGrid = (PetscInt)nGrid;
  if (numctx.nGrid != numTracers)
    SETERRQ(PETSC_COMM_WORLD, 1,
            "NUM has %d state variables but TMM was given %d tracers",
            nGrid, (int)numTracers);

  ierr = VecGetArrays(state->c, numTracers, &numctx.localTR);CHKERRQ(ierr);
  ierr = VecGetArrays(state->qef, numTracers, &numctx.localJTR);CHKERRQ(ierr);

  ierr = PetscMalloc(numctx.nGrid*sizeof(PetscScalar), &numctx.u);CHKERRQ(ierr);
  ierr = PetscMalloc(numctx.nGrid*sizeof(PetscScalar), &numctx.dudt);CHKERRQ(ierr);
  ierr = PetscMalloc(lSize*sizeof(PetscScalar), &numctx.localTs);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar), &numctx.locallatitude);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar), &numctx.localswrad);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar), &numctx.localtau);CHKERRQ(ierr);

  ierr = PetscMalloc(lSize*sizeof(PetscScalar), &numctx.localdz);CHKERRQ(ierr);
  ierr = VecLoadVecIntoArray(state->c[0], "dz.petsc", numctx.localdz);CHKERRQ(ierr);

  for (k = 0; k < lSize; k++) numctx.localTs[k] = Tconst;

  ierr = readProfileSurfaceScalarData("latitude.bin", numctx.locallatitude, 1);CHKERRQ(ierr);

  return 0;
}

PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop,
                                   TMMState state, void *ctx)
{
  PetscInt ip, kl, nzloc, k, itr;
  PetscScalar myTime, L, depth;

  myTime = 86400.0*tc;
  insolation_(&lNumProfiles, &myTime, &numctx.locallatitude[0], &daysPerYear,
              &numctx.localswrad[0], &numctx.localtau[0]);

  for (ip = 0; ip < lNumProfiles; ip++) {
    nzloc = lProfileLength[ip];
    kl = lStartIndices[ip];
    depth = 0.0;

    for (k = 0; k < nzloc; k++) {
      /* light at the middle of this layer */
      depth += 0.5*numctx.localdz[kl+k];
      L = numctx.localswrad[ip]*PARfrac*EinConv*exp(-kw*depth);
      depth += 0.5*numctx.localdz[kl+k];

      for (itr = 0; itr < numctx.nGrid; itr++)
        numctx.u[itr] = numctx.localTR[itr][kl+k];

      f_calcderivatives(numctx.u, L, numctx.localTs[kl+k], 0.0, numctx.dudt);

      for (itr = 0; itr < numctx.nGrid; itr++)
        numctx.localJTR[itr][kl+k] = numctx.dudt[itr];
    }
  }

  return 0;
}

PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop,
                                    TMMState state, void *ctx)
{
  return 0;
}

PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop,
                                           TMMState state, void *ctx)
{
  return 0;
}

PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx)
{
  PetscErrorCode ierr;

  ierr = VecRestoreArrays(state->c, state->numTracers, &numctx.localTR);CHKERRQ(ierr);
  ierr = VecRestoreArrays(state->qef, state->numTracers, &numctx.localJTR);CHKERRQ(ierr);
  ierr = PetscFree(numctx.u);CHKERRQ(ierr);
  ierr = PetscFree(numctx.dudt);CHKERRQ(ierr);
  ierr = PetscFree(numctx.localTs);CHKERRQ(ierr);
  ierr = PetscFree(numctx.localdz);CHKERRQ(ierr);
  ierr = PetscFree(numctx.locallatitude);CHKERRQ(ierr);
  ierr = PetscFree(numctx.localswrad);CHKERRQ(ierr);
  ierr = PetscFree(numctx.localtau);CHKERRQ(ierr);

  return 0;
}
