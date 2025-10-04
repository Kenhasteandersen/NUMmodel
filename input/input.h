!                                                                       !           
! ======================================================================!
! Input parameters for NUM model                                        !
! ======================================================================!
! This file contains the input parameters for the NUM Model.            !
! The parameters are read in in the different initialization            !
! routines with a call to read_input(filename).                         !   
! Comments can be added with an exclamation mark (!)                    !
! A descrition of each parameter is provided alongside the parameter    !
! along with the units in square brackets ([unit]). A glossary is found !
! at the end of the file where standard values for each parameter is    !
! given in curly brackets ({standard value}).                           !
!                                                                       !
!***********************************************************************!
!
!***********************************************************************
! GENERAL PARAMETERS
!***********************************************************************
!
  rhoCN = 5.68                  ! C:N mass ratio of cell [gC/gN]

! Variables for HTL mortalities:
!-------------------------------
  fracHTL_to_N = 0.5            ! Half becomes urine that is routed back to N
  fracHTL_to_POM = 0.5          ! Another half is fecal pellets that are routed back to the largest POM size class

!
!***********************************************************************
! GENERALISTS SIMPLE INPUT PARAMETERS
!***********************************************************************
!
  mMinGeneralist = 1.1623d-9    ! Smallest cell size [mug C]
  mMaxGeneralist = 1.0d0	! Largest cell size [mug C]

! Light uptake:
!--------------
  epsilonL = 0.8d0              ! Light uptake efficiency []
  alphaL = 0.3d0                ! Light affinity coef. [1/(uE/m2/s) 1/day um]
  rLstar = 7.5d0                ! Light affinity cross-over size [um]

! Dissolved nutrient and DOC uptake:
!-----------------------------------
  alphaN = 0.972d0                ! Diffusive affinity coefficient [L/d/mugC/um^2] 
  rNstar = 0.4d0                  ! Diffusive affinity cross-over size [um]

! Phagotrophy:
!-------------
  epsilonF = 0.8d0                ! Food assimilation efficiency [-]
  alphaF = 0.018d0                ! Clearance rate [L/d/ug C]
  cF = 30.0d0                      ! Max phagotrophy coefficient [um/day]
  beta = 500.d0                 ! Preferred predator-prey mass ratio
  sigma = 1.3d0                 ! Preferred predator-prey mass range

! Metabolism:
!------------
  cLeakage = 0.03d0             ! Passive leakage of C and N
  delta = 0.05d0                 ! Thickness of cell wall [um]
  alphaJ = 1.5d0                 ! Constant for jMax [day-1]
  cR = 0.1d0                     ! Basal metabolism relative to jMax [-]

! Biogeo:
!--------
  remin2 = 0.5d0                ! Fraction of viral lysis remineralized to DOC
  reminF = 0.1d0                ! Fraction of feeding losses remineralized



!***********************************************************************
! GENERALISTS INPUT PARAMETERS
! -generalists with explict metabolic costs
!***********************************************************************
  mMinGeneralist = 1.0d-9   ! Description [mug C]
  mMaxGeneralist = 1.0d0	! Description [mug C]
  rho = 0.4d-6

! Light uptake:
!--------------
  epsilonL = 0.8d0                ! Light uptake efficiency []
  alphaL = 0.3d0                  ! Scaling factor for light [unit]
  mUpperAlphaL = 3.003d0          ! Upper size for phototrophy
  rLstar = 7.5d0                  ! Description [unit]
  bL = 0.08d0                     ! Cost of light harvesting [mugC/mugC]

! Dissolved nutrient and DOC uptake:
!-----------------------------------
  alphaN = 0.972d0                ! Description [L/d/mugC/mum^2] 
  rNstar = 0.4d0                  ! Description [mum]
  bN = 0.3d0                     	! cost of N uptake [mugC/mugN]
  bDOC = 0.3d0                   	! cost of DOC uptake [mugC/mugC]

! Phagotrophy:
!-------------
  epsilonF = 0.8d0                ! Food Assimilation efficiency [unit]
  alphaF = 0.018d0                ! Food affinity scaling factor [L mug C-1 d-1]
  cF = 30.0d0
  beta = 500.d0
  sigma = 1.3d0
  bF = 0.3d0                     	! Cost of food uptake [mugC/mugC]

! Metabolism:
!------------
  cLeakage = 0.03d0               ! Passive leakage of C and N
  delta = 0.05d0                  ! Thickness of cell wall [mum]
  alphaJ = 1.5d0                  ! Constant for jmax. [day-1]
  cR = 0.03d0       
  bg = 0.2d0                      ! Cost of biosynthsesis

! Biogeo:
!--------
  remin2 = 0.5d0                ! Fraction of virulysis remineralized to DOC
  reminF = 0.1d0

!***********************************************************************
! DIATOMS SIMPLE INPUT PARAMETERS
!***********************************************************************
  mMin = 1.1623d-8            	! Smallest size [mug C]
  rhoCSi = 3.4                  ! Carbon:Si mass ratio
  v  = 0.6                      ! Vacuole volume fraction

! Light uptake:
!--------------
  epsilonL = 0.8                ! Light uptake efficiency []
  alphaL = 0.3                  ! Scaling factor for light [unit]
  rLstar = 7.5                  ! Description [unit]

! Dissolved nutrient uptake:
!---------------------------
  alphaN = 0.972                ! alphaN/(1-v)[L/d/mugC/mum^2] 
  rNstar = 0.4                  ! Description [mum]
  bN = 0.0                      ! cost of N uptake [mugC/mugN]

! Silicate uptake:
!-----------------
  bSi = 0.0                     ! Cost of food uptake [gC/gSi]

! Metabolism:
!------------
  cLeakage = 0.03               ! Passive leakage of C and N
  delta = 0.05                  ! Thickness of cell wall [mum]
  alphaJ = 1.5                  ! Constant for jmax. [day-1]
  cR = 0.03       

! Biogeo:
!--------
  remin2 = 0.5d0                ! Fraction of virulysis remineralized to DOC

! Vulnerability to predation:
!----------------------------
  palatability = 0.5
  
  
!***********************************************************************
! DIATOMS INPUT PARAMETERS
!***********************************************************************
  mMinDiatom = 1.d-6            ! Description [mug C]
  mMaxDiatom = 0.01		          ! Description [mug C]
  rhoCSi = 3.4                  ! Carbon:Si mass ratio
  v  = 0.8                      ! Vacuole volume fraction

! Light uptake:
!--------------
  epsilonL = 0.8                ! Light uptake efficiency []
  alphaL = 0.3                  ! Scaling factor for light [unit]
  rLstar = 7.5                  ! Description [unit]
  bL = 0.08                     ! Cost of light harvesting [mugC/mugC]

! Dissolved nutrient uptake:
!---------------------------
  alphaN = 0.972                ! Description [L/d/mugC/mum^2] 
  rNstar = 0.4                  ! Description [mum]
  bN = 0.3                      ! cost of N uptake [mugC/mugN]
  bDOC = 0.3                    ! cost of DOC uptake [mugC/mugC]

! Silicate uptake:
!-----------------
  bSi = 0.3                    	! Cost of food uptake [gC/gSi]

! Metabolism:
!------------
  cLeakage = 0.03              	! Passive leakage of C and N
  delta = 0.05                  ! Thickness of cell wall [mum]
  alphaJ = 1.5                  ! Constant for jmax. [day-1]
  cR = 0.03       
  bg = 0.2                      ! Cost of biosynthsesis

! Biogeo:
!--------
! remin = 0.0                   ! Fraction of mortality losses reminerilized to DOC
  remin2 = 0.5d0                ! Fraction of virulysis remineralized to DOC

! Vulnerability to predation:
!----------------------------
  palatability = 0.5
  
!***********************************************************************
! COPEPODS ACTIVE INPUT PARAMETERS
! - Values taken from Serra-Pompei et al (2020) for actively feeding copepods,
! but adjusted to a reference temperature of 10 degrees
!***********************************************************************

  epsilonF = 0.67               ! Assimilation efficiency
  epsilonR = 0.25               ! Reproductive efficiency
  beta = 10000.                 ! Preferred predator-prey mass ratio
  sigma = 1.5                   ! Preferred predator-prey mass range
  alphaF = 0.011                ! Clearance rate coefficient
  q = 0.75                      ! Exponent of clearance rate
  h = 0.97                      ! Coefficient for maximum ingestion rate
  hExponent = 0.75              ! Exponent for maximum ingestion rate

! kBasal  is a factor for basal metabolism {0.006}. This value represents basal metabolism at 
! starvation. Following Kiørboe (1985) the starvation metabolism is approximatly
! 0.2*0.18=0.036 times the maximum metabolism (kSDA). Increased to 0.01 to avoid 
! too long transients. 
!---------------------
  kBasal = 0.01  
       
! kSDA = Factor for SDA metabolism (Serra-Pompei 2020) {0.16}. This value assumes that the 
! data in Kiørboe and Hirst (2014) are for fully fed copepods.   
!------------------------------------------------------------       
  kSDA = 0.16                   
   
  AdultOffspring = 100.         ! Adult:offspring mass ratio [-]
  vulnerability = 1.            ! Active copepods have full risk of predation
  DiatomsPreference = 1.0       ! Feeding preference on diatoms

!***********************************************************************
! COPEPODS PASSIVE INPUT PARAMETERS
! - Values taken from Serra-Pompei et al (2020) for passively feeding copepods
! but adjusted to reference temperature of 10 degrees
!***********************************************************************

   epsilonF = 0.67               ! Assimilation efficiency
   epsilonR = 0.25               ! Reproductive efficiency
   beta = 100.                   ! Preferred predator-prey mass ratio
   sigma = 1.                    ! Preferred predator-prey mass range
   alphaF = 0.0052               ! Clearance rate coefficient
   q = 0.75                      ! Exponent of clearance rate
   h = 0.29                      ! Coefficient for maximum ingestion rate
   hExponent = 0.75              ! Exponent for maximum ingestion rate

! kBasal  is a factor for basal metabolism {0.006}. This value represents basal metabolism at 
! starvation. Following Kiørboe (1985) the starvation metabolism is approximatly
! 0.2*0.18=0.036 times the maximum metabolism (kSDA). Increased to 0.01 to avoid 
! too long transients. 
!---------------------
  kBasal = 0.01  
       
! kSDA = Factor for SDA metabolism (Serra-Pompei 2020) {0.16}. This value assumes that the 
! data in Kiørboe and Hirst (2014) are for fully fed copepods.   
!------------------------------------------------------------       
  kSDA = 0.16                   
   
  AdultOffspring = 100.         ! Adult:offspring mass ratio [-]
  vulnerability = 0.2           ! Passive copepods have reduced risk of predation
  DiatomsPreference = 0.2       ! Feeding preference on diatoms (lowered feeding on diatoms)

!***********************************************************************
! PARTICULATE ORGANIC MATTER (POM) INPUT PARAMETERS
!*********************************************************************** 
  mMin = 1.d-9                  ! Smallest POM mass
  remin = 0.07d0                ! remineralisation rate (1/day) (Serra-Pompei (2022)) @10 degrees
  palatability = 0.1d0          ! Preference of other groups for eating POM

  
