module read_input_module
  use globals
  implicit none
  
  
 
  ! CN mass ratio:
  real(dp) :: rhoCN 
  real(dp) :: rhoCSi
  real(dp) :: v                 ! Vacuole fraction

  !
  ! Light uptake:
  !
  real(dp) :: epsilonL          ! Light uptake efficiency
  real(dp) :: alphaL  
  real(dp) :: rLstar  
  real(dp) :: bL                ! cost of light harvesting mugC(mugC)^-1
  !
  ! Dissolved nutrient uptake:
  !
  real(dp) :: alphaN            ! L/d/mugC/mum^2
  real(dp) :: rNstar            ! mum
  real(dp) :: bN                ! cost of N uptake mugC(mugN)^-1
  
  
  real(dp) :: bDOC              ! cost of DOC uptake mugC(mugN)^-1
  real(dp) :: bSi               ! cost of Si uptake mugC(mugSi)^-1

  !
  ! Phagotrophy:
  !
  real(dp) :: epsilonF          ! Assimilation efficiency
  real(dp) :: alphaF            ! Clearance rate coefficient
  real(dp) :: cF 
  real(dp) :: beta 
  real(dp) :: sigma 
  real(dp) :: bF                ! cost of food uptake mugC(mugSi)^-1
  real(dp) :: q                 ! Exponent of clerance rate
  real(dp) :: h                 ! Factor for maximum ingestion rate
  real(dp) :: hExponent         ! Exponent for maximum ingestions rate
  real(dp) :: kBasal            ! 0.006 ! Factor for basal metabolism.This value represents basal
  real(dp) :: kSDA              ! Factor for SDA metabolism (Serra-Pompei 2020). This value assumes that the
  real(dp) :: AdultOffspring    ! Adult-offspring mass ratio
  real(dp) :: DiatomsPreference
  !
  ! Metabolism
  !
  real(dp) :: cLeakage          ! passive leakage of C and N
  real(dp) :: delta             ! Thickness of cell wall in mum
  real(dp) :: alphaJ            ! Constant for jmax.  per day
  real(dp) :: cR 
  real(dp) :: bg                ! cost of biosynthsesis -- parameter from literature pending
  
  !
  ! Reproduction
  !
  real(dp) :: epsilonR          ! Reproductive efficiency

  !
  ! Biogeo:
  !
  real(dp) :: remin             ! remineralisation rate of POM (1/day) (Serra-Pompei (2022)) @10 degrees
  real(dp) :: remin2            ! fraction of virulysis remineralized to N and DOC
  real(dp) :: reminF            ! fraction of feeding losses to DOC
  !
  ! Max and min sizes
  !
  real(dp) :: mMinGeneralist
  real(dp) :: mMaxGeneralist
  real(dp) :: mMinDiatom
  real(dp) :: mMaxDiatom
  real(dp) :: mMin
  
  real(dp) :: fracHTL_to_N     ! Half becomes urine that is routed back to N
  real(dp) :: fracHTL_to_POM   ! Another half is fecal pellets that are routed back to the largest POM size class
  
  !
  ! Predation risk:
  !
  real(dp) :: palatability
  real(dp) :: vulnerability     ! Passed to "palatability" in the parent spectrum class
  


  contains
        
  subroutine read_input(filename,listname)
    character(*), intent(in) :: filename, listname
    character(len=100) :: line, keyword, val, thislist
    character(len=100) :: str_General
    character(len=100) :: str_Generalists_simple
    character(len=100) :: str_Generalists
    character(len=100) :: str_Diatoms_simple
    character(len=100) :: str_Diatoms
    character(len=100) :: str_Copepods_passive
    character(len=100) :: str_Copepods_active
    character(len=100) :: str_POM
    integer:: ios, n, l=1
  
    !------------------------------------------------------------------
    !  Open file and find keywords and parameter values
    !------------------------------------------------------------------
    print*, 'Loading parameter for ', listname, ' from file ', filename,':' 
    str_General='! GENERAL PARAMETERS'
    str_Generalists_simple='! GENERALISTS SIMPLE INPUT PARAMETERS'
    str_Generalists='! GENERALISTS INPUT PARAMETERS'
    str_Diatoms_simple='! DIATOMS SIMPLE INPUT PARAMETERS'
    str_Diatoms='! DIATOMS INPUT PARAMETERS'
    str_Copepods_passive='! COPEPODS PASSIVE INPUT PARAMETERS'
    str_Copepods_active='! COPEPODS ACTIVE INPUT PARAMETERS'
    str_POM='! PARTICULATE ORGANIC MATTER (POM) INPUT PARAMETERS'
    thislist='no list defined yet'
    
    OPEN(1, file = filename)
      DO WHILE(.TRUE.)
        read(1, '(A)', iostat=ios) line
        IF (ios /= 0) exit
    !
    ! Check if we are on a group define statement
    !
          if (line(1:len(str_General)).eq.str_General) then
            thislist='general'
          else if (line(1:len(str_Generalists_simple)).eq.str_Generalists_simple) then
            thislist='generalists_simple'
          else if (line(1:len(str_Generalists)).eq.str_Generalists) then
            thislist='generalists'
          else if (line(1:len(str_Diatoms_simple)).eq.str_Diatoms_simple) then
            thislist='diatoms_simple'
          else if (line(1:len(str_Diatoms)).eq.str_Diatoms) then
            thislist='diatoms'
          else if (line(1:len(str_Copepods_passive)).eq.str_Copepods_passive) then
            thislist='copepods_passive'
          else if (line(1:len(str_Copepods_active)).eq.str_Copepods_active) then
            thislist='copepods_active'
          else if (line(1:len(str_POM)).eq.str_POM) then
            thislist='POM'
          end if   
                    
    !
    ! Remove leading blank spaces
    !
          line=TRIM(ADJUSTL(line(1:LEN(line))))
    !
    ! If the line is not a comment and it is on the right list then:
    !
           if ((line(1:1).ne.'!').and.(LEN_TRIM(line).gt.0).and.(TRIM(thislist).eq.TRIM(listname))) THEN
    !
    ! remove comments, if any
    !
            n=index(line,'!')
            IF (n /=0) THEN
              line=TRIM(ADJUSTL(line(1:n-1)))
            END IF
    !
    ! Find equal sign
    !
            n=index(line,'=')
    !
    ! Isolate keywords and values
    !
            keyword=TRIM(ADJUSTL(line(1:n-1)))
            val=TRIM(ADJUSTL(line(n+1:LEN(line))))

            print*, '   ', TRIM(keyword), ' = ', val
    !------------------------------------------------------------------
    !  Assign keyword values to parameters
    !------------------------------------------------------------------
            SELECT CASE (TRIM(keyword))
              CASE ('rhoCN')
                read(val,*) rhoCN
              CASE ('rhoCSi')
                read(val,*) rhoCSi
              CASE ('v')
                read(val,*) v
              CASE ('fracHTL_to_N')
                read(val,*) fracHTL_to_N
              CASE ('fracHTL_to_POM')
                read(val,*) fracHTL_to_POM
              CASE ('epsilonL')
                read(val,*) epsilonL
              CASE ('alphaL')
                read(val,*) alphaL
              CASE ('rLstar')
                read(val,*) rLstar
              CASE ('bL')
                read(val,*) bL
              CASE ('alphaN')
                read(val,*) alphaN
              CASE ('rNstar')
                read(val,*) rNstar
              CASE ('bN')
                read(val,*) bN
              CASE ('bDOC')
                read(val,*) bDOC
              CASE ('bSi')
                read(val,*) bSi
              CASE ('epsilonF')
                read(val,*) epsilonF
              CASE ('alphaF')
                read(val,*) alphaF
              CASE ('cF')
                read(val,*) cF                
              CASE ('beta')
                read(val,*) beta                
              CASE ('sigma')
                read(val,*) sigma                
              CASE ('bF')
                read(val,*) bF                
              CASE ('q')
                read(val,*) q                
              CASE ('h')
                read(val,*) h                  
              CASE ('hExponent')
                read(val,*) hExponent                  
              CASE ('kBasal')
                read(val,*) kBasal                  
              CASE ('kSDA')
                read(val,*) kSDA    
              CASE ('AdultOffspring')
                read(val,*) AdultOffspring              
              CASE ('DiatomsPreference')
                read(val,*) DiatomsPreference   
              CASE ('cLeakage')
                read(val,*) cLeakage                
              CASE ('delta')
                read(val,*) delta                  
              CASE ('alphaJ')
                read(val,*) alphaJ                  
              CASE ('cR')
                read(val,*) cR                  
              CASE ('bg')
                read(val,*) bg    
              CASE ('epsilonR')
                read(val,*) epsilonR 
              CASE ('remin')
                read(val,*) remin             
              CASE ('remin2')
                read(val,*) remin2  
              CASE ('reminF')
                read(val,*) reminF   
              CASE ('palatability')
                read(val,*) palatability                 
              CASE ('vulnerability')
                read(val,*) vulnerability  
              CASE ('mMinGeneralist')
                read(val,*) mMinGeneralist  
              CASE ('mMaxGeneralist')
                read(val,*) mMaxGeneralist  
              CASE ('mMinDiatom')
                read(val,*) mMinDiatom  
              CASE ('mMaxDiatom')
                read(val,*) mMaxDiatom  
              CASE ('mMin')
                read(val,*) mMin  
            END SELECT
                   
        END IF
        l=l+1
      END DO
    CLOSE(1)
    
    20    continue
     
  end subroutine read_input
                        

end module read_input_module


