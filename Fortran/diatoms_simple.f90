!
! Module to handle diatoms
! Maximal simplicity; based on generalists with the only difference:
!  - a fixed-sized vacuole
!  - lack of ability to do phagotrophy
!  - reliance on silicate
!  - lower predation risk due to silicate shell
!  - inability to take up DOC
!
module diatoms_simple
    use globals
    use spectrum
    use read_input_module
    implicit none

    private

    type, extends(spectrumUnicellular) :: spectrumDiatoms_simple
      real(dp), dimension(:), allocatable:: JSi

      contains
      procedure, pass :: initDiatoms_simple
      procedure :: calcRates => calcRatesDiatoms_simple
      procedure :: calcDerivativesDiatoms_simple
      procedure :: printRates => printRatesDiatoms_simple
      procedure :: getNbalance
      procedure :: getCbalance
      procedure :: getSiBalance
    end type spectrumDiatoms_simple

    public spectrumDiatoms_simple, initDiatoms_simple, calcRatesDiatoms_simple
    public calcDerivativesDiatoms_simple, printRatesDiatoms_simple
  contains
      
    subroutine initDiatoms_simple(this, n, mMax)
      class(spectrumDiatoms_simple):: this
      real(dp), intent(in):: mMax
      integer, intent(in):: n
      real(dp), parameter:: rho = 0.4*1d6*1d-12
      call read_input(inputfile,'diatoms_simple')
      call this%initUnicellular(n, mMin, mMax)
      allocate(this%JSi(this%n))
      !
      ! Radius:
      !
      this%r = (threequarters/pi * this%m/rho/(1-v))**onethird  ! Andy's approximation
      this%nu = 6**twothirds*pi**onethird*delta * (this%m/rho)**(-onethird) * &
        (v**twothirds + (1.+v)**twothirds)

      ! Affinities are the same as for the generalists divided by a factor (1-v):
      this%AN = alphaN * this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m / (1-v)
      this%AL = alphaL/this%r * (1-exp(-this%r/rLstar)) * this%m * (1.d0-this%nu) / (1-v)
      this%AF = 0.d0
      this%JFmax = 0.d0

      this%JlossPassive = cLeakage/this%r * this%m ! in units of C
  
      this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
      this%Jresp = cR*alphaJ*this%m
      this%JDOCreal = 0.d0
  
      this%beta = 0.d0 ! No feeding
      this%palatability = palatability ! Lower risk of predation
    end subroutine initDiatoms_simple
 
    subroutine calcRatesDiatoms_simple(this, L, N, Si, gammaN, gammaSi)
      class(spectrumDiatoms_simple), intent(inout):: this
      real(dp), intent(in):: gammaN, gammaSi
      real(dp), intent(in):: L, N, Si
      real(dp):: JmaxT
      integer:: i
  
      do i = 1, this%n
         !
         ! Uptakes
         !
         this%JN(i) =   fTemp15 * gammaN * this%AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
         !this%JDOC(i) = gammaDOC * this%AN(i)*DOC ! Diffusive DOC uptake, units of C/time
         this%JSi(i) = fTemp15 * gammaSi * this%AN(i)*Si*rhoCSi! Diffusive Si uptake, units of C/time
         this%JL(i) =   epsilonL * this%AL(i)*L  ! Photoharvesting
         !
         ! Estimate the limiting growth nutrient (Liebig):
         !
         JmaxT = max(0.d0, fTemp2 * this%Jmax(i))
         this%Jtot(i) = min( JmaxT, this%JN(i), this%JL(i)-ftemp2*this%Jresp(i), this%JSi(i) )
         !    
         ! Account for possible carbon limitation due to carbon costs of uptakes:
         !
         this%Jtot(i) = min( this%Jtot(i), &
            this%JL(i)-ftemp2*this%Jresp(i) - (bN/rhoCN + bSi/rhoCSi)*this%Jtot(i) )
         
         ! Jtot synthesis limitation:
        if (this%Jtot(i) .gt. 0.) then
              this%Jtot(i) = JmaxT*this%Jtot(i) / ( this%Jtot(i) + JmaxT )
        end if

        this%JLreal(i) = this%JL(i)
        this%Jresptot(i) = fTemp2*this%Jresp(i)
      end do
    end subroutine calcRatesDiatoms_simple
  
    subroutine calcDerivativesDiatoms_simple(this, u, dNdt, dSidt, dudt)
      class(spectrumDiatoms_simple), intent(inout):: this
      real(dp), intent(in):: u(this%n)
      real(dp), intent(inout):: dNdt, dSidt, dudt( this%n )
      real(dp):: mortloss
      integer:: i
  
      this%mort2 = this%mort2constant*u
      this%jPOM = (1-remin2)*this%mort2 ! non-remineralized mort2 => POM

      do i = 1, this%n
        mortloss = remin2 * this%mort2(i)*u(i) !+reminHTL*this%mortHTL(i))
        !
        ! Update nitrogen:
        !
        dNdt = dNdt  &
            - ((this%Jtot(i))*u(i)/this%m(i) &
            + mortloss)/rhoCN
        !
        ! Update DOC:
        !
        !this%dudt(idxDOC) = 0.d0 !this%dudt(idxDOC) &
             !+ (-this%JDOC(i) &
             !+ this%JCloss(i))*u(i)/this%m(i) &
             !+ mortloss
        !
        ! Update Si:
        !
        dSidt = dSidt &
             + ((-this%Jtot(i))*u(i)/this%m(i) & ! Uptakes of silicate
             + mortloss  & ! all mortality due to remineralisation, HTL and virulysis is returned
             + this%mortpred(i) )/rhoCSi ! Silicate relased during predation by generalists is assumed to be remineralized immidiately
        !
        ! Update the diatoms:
        !
        dudt(i) = (this%Jtot(i)/this%m(i)  &
             - this%mortpred(i) &
             - this%mort2(i) &
             - this%mortHTL(i))*u(i)

     end do
   end subroutine calcDerivativesDiatoms_simple

  subroutine printRatesDiatoms_simple(this)
     class(spectrumDiatoms_simple), intent(in):: this

     write(*,*) "Diatoms_simple with ", this%n, " size classes:"
     call this%printRatesUnicellular()

     99 format (a10, 20f10.6)
 
     write(*,99) "jSi:", this%JSi / this%m
   end subroutine printRatesDiatoms_simple

  end module diatoms_simple
