module generalists
  use globals
  use sizespectrum
  implicit none

  real(dp), parameter:: rhoCN = 5.68
  real(dp), parameter:: epsilonL = 0.9 ! Light uptake efficiency
  real(dp), parameter:: epsilonFF = 0.8 ! Assimilation efficiency
  real(dp), parameter:: cLeakage = 0.00015 ! passive leakage of C and N
  real(dp), parameter:: c = 0.0015 ! Parameter for cell wall fraction of mass.
            !The constant is increased a bit to limit the lower cell size
  real(dp), parameter:: alphaN = 0.00012 !0.00004 % 0.000162 % Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g)^(-1/3) / 1.5 (g/cm); Andersen et al 2015
  real(dp), parameter:: cN = 0.1
  real(dp), parameter:: alphaL = 0.000914 ! if using Andys shading formula for non-diatoms
  real(dp), parameter:: cL = 21 ! if using Andys shading formula for non-diatoms
  real(dp), parameter:: alphaF = 0.018 !  Fits to TK data for protists
  real(dp), parameter:: cF = 0.6 ! Just a guess
  real(dp), parameter:: alphaJ = 1.5 ! Constant for jmax.  per day
  real(dp), parameter:: cR = 0.1
  real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to N and DOC
  real(dp), parameter:: remin2 = 1.d0 ! fraction of virulysis remineralized to N and DOC
  
  real(dp),  dimension(:), allocatable:: AN(:), AL(:), Jmax(:),  JlossPassive(:)
  real(dp),  dimension(:), allocatable:: nu(:), mort(:)
  real(dp),  dimension(:), allocatable:: JN(:), JL(:), Jresp(:), JFreal(:)
  real(dp):: mort2

contains

  function initGeneralists(n, ixStart) result(this)
      type(typeSpectrum):: this
      integer, intent(in):: n, ixStart

      this = initSpectrum(n, ixStart, 10.d0**(-8.5), 1.d0)

      allocate(AN(n))
      allocate(AL(n))
      allocate(Jresp(n))
      allocate(JlossPassive(n))
      allocate(nu(n))
      allocate(mort(n))

      allocate(JN(n))
      allocate(JL(n))
      allocate(JFreal(n))

      AN = alphaN*this%m**onethird / (1 + cN*this%m**(-onethird))
      AL = alphaL*this%m**twothirds * (1-exp(-cL*this%m**onethird ))  ! shading formula
      this%AF = alphaF*this%m
      JlossPassive = cLeakage * this%m**twothirds ! in units of C
      this%JFmax = 0.6*this%m**twothirds
      nu = c * this%m**(-onethird)
      Jmax = alphaJ * this%m * (1.d0-nu) ! mugC/day
      Jresp = cR*alphaJ*this%m
      mort = 0*0.005*(Jmax/this%m) * this%m**(-0.25);
      mort2 = 0.0002*n
      !
      ! Export grid to globals:
      !
      m(this%ixStart : this%ixEnd) = this%m
      beta(this%ixStart : this%ixEnd) = 500.
      sigma(this%ixStart : this%ixEnd) = 1.3
      !
      ! Export feeding parameters:
      !
      AF(this%ixStart : this%ixEnd) = this%AF
      JFmax(this%ixStart : this%ixEnd) = this%JFmax
      epsilonF(this%ixStart : this%ixEnd) = epsilonFF
      
    end function initGeneralists

    subroutine calcRatesGeneralists(this, u, rates, L, gammaN, gammaDOC)
      type(typeSpectrum), intent(in):: this
      real(dp), intent(in):: u(:), gammaN, gammaDOC
      type(typeRates), intent(inout):: rates
      real(dp), intent(in):: L
      real(dp):: N, DOC, f
      integer:: ix, i 
      !
      ! Make sure values are positive
      !
      N = max(0., u(idxN))
      DOC = max(0., u(idxDOC))

      do ix = this%ixStart, this%ixEnd
         i = ix-this%ixStart+1
         !
         ! Uptakes
         !
         rates%JN(ix) =   gammaN * AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
         rates%JDOC(ix) = gammaDOC * AN(i)*DOC ! Diffusive DOC uptake, units of C/time
         rates%JL(ix) =   epsilonL * AL(i)*L  ! Photoharvesting
         ! Total nitrogen uptake:
         rates%JNtot(ix) = rates%JN(ix)+rates%JF(ix)-Jlosspassive(i) ! In units of C
         ! Total carbon uptake
         rates%JCtot(ix) = rates%JL(ix)+rates%JF(ix)+rates%JDOC(ix)-Jresp(i)-JlossPassive(i)
         ! Liebig + synthesis limitation:
         rates%Jtot(ix) = min( rates%JNtot(ix), rates%JCtot(ix) )
         f = rates%Jtot(ix)/(rates%Jtot(ix) + Jmax(i))
         ! If synthesis-limited then down-regulate feeding:
         if (rates%Jtot(ix) .gt. 0) then
            JFreal(i) = max(0.d0, rates%JF(ix) - (rates%Jtot(ix)-f*Jmax(i)))
         else
            JFreal(i) = max(0.d0, rates%JF(ix))
         end if
         rates%Jtot(ix) = f * Jmax(i)
         rates%JLreal(ix) = rates%JL(ix) - max( 0.d0, &
              min((rates%JCtot(ix) - (rates%JF(ix)-JFreal(i))-rates%Jtot(ix)), rates%JL(ix)))
         
         ! Actual uptakes:
         rates%JCtot(ix) = &
              + rates%JLreal(ix)  &
              + rates%JDOC(ix)  &
              + JFreal(i)  &
              - Jresp(i)  &
              - JlossPassive(i)
         rates%JNtot(ix) = &
              rates%JN(ix) + &
              JFreal(i) - &
              JlossPassive(i)
         !
         ! Losses:
         !
         rates%JCloss_feeding(ix) = (1.-epsilonFF)/epsilonFF*JFreal(i) ! Incomplete feeding (units of carbon per time)
         rates%JCloss_photouptake(ix) = (1.-epsilonL)/epsilonL * rates%JLreal(ix)
         rates%JNlossLiebig(ix) = max( 0.d0, rates%JNtot(ix)-rates%Jtot(ix))  ! In units of C
         rates%JClossLiebig(ix) = max( 0.d0, rates%JCtot(ix)-rates%Jtot(ix)) ! C losses from Liebig, not counting losses from photoharvesting

         rates%JNloss(ix) = &
              rates%JCloss_feeding(ix) + &
              rates%JNlossLiebig(ix) +&
              JlossPassive(i) ! In units of C
         rates%JCloss(ix) = &
              rates%JCloss_feeding(ix) + &
              rates%JCloss_photouptake(ix) + &
              rates%JClossLiebig(ix) +&
              JlossPassive(i)
         rates%JF(ix) = JFreal(i)
      end do
    end subroutine calcRatesGeneralists
    
    subroutine calcDerivativesGeneralists(this, u, rates)
      type(typeSpectrum), intent(in):: this
      type(typeRates), intent(inout):: rates
      real(dp), intent(in):: u(:)
      real(dp):: mortloss
      
      integer:: i, ix
      
      do ix = this%ixStart, this%ixEnd
         i = ix - this%ixStart+1
         mortloss = u(ix)*((1.d0-remin2)*mort2*u(ix) + rates%mortHTL(ix))
         !
         ! Update nitrogen:
         !
         rates%dudt(idxN) = rates%dudt(idxN)  &
              + (-rates%JN(ix) &
              + rates%JNloss(ix))*u(ix)/this%m(i) &
              + (remin2*mort2*u(ix)*u(ix) &
              + remin*mortloss)/rhoCN
         !
         ! Update DOC:
         !
         rates%dudt(idxDOC) = rates%dudt(idxDOC) &
            + (-rates%JDOC(ix) &
            + rates%JCloss(ix))*u(ix)/this%m(i) &
            + remin2*mort2*u(ix)*u(ix) &
            + remin*mortloss
         !
         ! Update the generalists:
         !
         rates%dudt(ix) = (rates%Jtot(ix)/this%m(i)  &
            - mort(i) &
            - rates%mortpred(ix) &
            - mort2*u(ix) &
            - rates%mortHTL(ix))*u(ix)
      end do
    end subroutine calcDerivativesGeneralists
    
  end module generalists
