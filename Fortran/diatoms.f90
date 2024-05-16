!
! Module to handle diatoms
! AMALIA latest version
! with Down-regulation 
!
module diatoms
     use globals
     use spectrum
     use read_input_module
     implicit none
 
     private
     real(dp) :: bN, bL, bDOC, bSi, bg, epsilonL
     real(dp) :: rhoCSi, remin2
     
     type, extends(spectrumUnicellular) :: spectrumDiatoms
     real(dp), dimension(:), allocatable:: JSi,JSireal, jNet

     contains
       procedure, pass :: initDiatoms
       procedure :: calcRates => calcRatesDiatoms
       procedure :: calcDerivativesDiatoms
       procedure :: printRates => printRatesDiatoms
       procedure :: getProdNet => getProdNetDiatoms
     end type spectrumDiatoms

     public  initDiatoms, spectrumDiatoms, calcRatesDiatoms, calcDerivativesDiatoms
     public printRatesDiatoms

   contains
       
     subroutine initDiatoms(this, n,errorio,errorstr)
       use iso_c_binding, only: c_char
       class(spectrumDiatoms):: this
       integer, intent(in):: n
       logical(1), intent(out):: errorio 
       character(c_char), dimension(*), intent(out) :: errorstr
       integer:: i
       real(dp), parameter:: mMin = 3.1623d-9
       real(dp), parameter:: rho = 0.4*1d-6
       real(dp) :: mMinDiatom, mMaxDiatom, v, alphaL, rLstar,  alphaN
       real(dp) :: rNstar, cLeakage, delta, alphaJ, cR, palatability
       !real(dp) :: fl
       
       ! no errors to begin with
       errorio=.false.
       
       print*, 'Loading parameter for diatoms from ', inputfile, ':'
       call read_input(inputfile,'diatoms','mMinDiatom',mMinDiatom,errorio,errorstr)
       call read_input(inputfile,'diatoms','mMaxDiatom',mMaxDiatom,errorio,errorstr)
       call this%initUnicellular(n, mMinDiatom, mMaxDiatom)
       call read_input(inputfile,'diatoms','rhoCSi',rhoCSi,errorio,errorstr)
       call read_input(inputfile,'diatoms','v',v,errorio,errorstr)
       call read_input(inputfile,'diatoms','epsilonL',epsilonL,errorio,errorstr)
       call read_input(inputfile,'diatoms','alphaL',alphaL,errorio,errorstr)
       call read_input(inputfile,'diatoms','rLstar',rLstar,errorio,errorstr)
       call read_input(inputfile,'diatoms','bL',bL,errorio,errorstr)
       call read_input(inputfile,'diatoms','alphaN',alphaN,errorio,errorstr)
       call read_input(inputfile,'diatoms','rNstar',rNstar,errorio,errorstr)
       
       call read_input(inputfile,'diatoms','bN',bN,errorio,errorstr)
       call read_input(inputfile,'diatoms','bDOC',bDOC,errorio,errorstr)
       call read_input(inputfile,'diatoms','bSi',bSi,errorio,errorstr)
       call read_input(inputfile,'diatoms','cLeakage',cLeakage,errorio,errorstr)
       call read_input(inputfile,'diatoms','delta',delta,errorio,errorstr)
       call read_input(inputfile,'diatoms','alphaJ',alphaJ,errorio,errorstr)
       call read_input(inputfile,'diatoms','cR',cR,errorio,errorstr)
       call read_input(inputfile,'diatoms','bg',bg,errorio,errorstr)
       call read_input(inputfile,'diatoms','remin2',remin2,errorio,errorstr)
       call read_input(inputfile,'diatoms','palatability',palatability,errorio,errorstr)
       
       allocate(this%JSi(this%n))
       allocate(this%JSireal(this%n))
       allocate(this%jNet(this%n))
       !
       ! Radius:
       !
       this%r = (threequarters/pi * this%m/rho/(1-v))**onethird  ! Andy's approximation
       this%nu = 6**twothirds*pi**onethird*delta * (this%m/rho)**(-onethird) * &
              (v**twothirds + (1.+v)**twothirds)
       do i = 1,this%n
         this%nu(i) = min(1.d0, this%nu(i))
       enddo
       
       ! Affinities are the same as for the generalists divided by a factor (1-v):
       this%AN = alphaN * this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m / (1-v)
       this%AL = alphaL/this%r * (1-exp(-this%r/rLstar)) * this%m * (1.d0-this%nu) / (1-v)
       this%AF = 0.d0 ! No feeding
       this%JFmax = 0.d0
 
       this%JlossPassive = cLeakage/this%r * this%m ! in units of C
 
       this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
       this%Jresp = cR*alphaJ*this%m 
   
       this%beta = 0.d0 ! No feeding
       this%palatability = palatability ! Lower risk of predation
     end subroutine initDiatoms
  
     subroutine calcRatesDiatoms(this, L, N,DOC, Si, gammaN, gammaDOC, gammaSi)
       class(spectrumDiatoms), intent(inout):: this
       real(dp), intent(in):: gammaN, gammaSi, gammaDOC
       real(dp), intent(in):: L, N, Si, DOC
       real(dp):: JmaxT, tmp
       real(dp):: dN(this%n), dSi(this%n), Jnetp(this%n), Jnet(this%n)
       integer:: i
   
       do i = 1, this%n
          !
          ! Uptakes
          !
          this%JN(i) =  ftemp15* gammaN * this%AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
          this%JDOC(i) = ftemp15*gammaDOC * this%AN(i)*DOC ! Diffusive DOC uptake, units of C/time
          this%JSi(i) = ftemp15* gammaSi * this%AN(i)*Si*rhoCSi! Diffusive Si uptake, units of C/time
          this%JL(i) =  epsilonL * this%AL(i)*L  ! Photoharvesting
          JmaxT = fTemp2*this%Jmax(i)
          ! 
          ! Calculate downregulation of nutrient uptakes:
          !
          Jnetp(i) = this%JL(i)*(1-bL) + this%JDOC(i)*(1-bDOC) - ftemp2*this%Jresp(i)
          tmp = Jnetp(i) / (1.d0 + bg + bN + bSi)
          ! N:
          if (this%JN(i) .eq. 0.d0) then
            dN(i) = 1.d0
          else
            dN(i) = min( 1.d0, tmp/this%JN(i) )
          end if
          ! Silicate:
          if (this%JSi(i) .eq. 0.d0) then
            dSi(i) = 1.d0
          else
            dSi(i) = min( 1.d0, tmp/this%JSi(i) )
          end if
          Jnet(i) = 1.d0/(1+bg) * ( this%JDOC(i)*(1-bDOC) + this%JL(i)*(1-bL) &
             - ftemp2*this%Jresp(i) - bN*dN(i)*this%JN(i) - bSi*dSi(i)*this%JSi(i) )
          !
          ! Adjust Jnet according to limitation:
          !
          if ( JNetp(i) .lt. 0.d0) then  ! Check for severe light limitation
            Jnet(i) = -ftemp2*this%Jresp(i) + this%JDOC(i)*(1-bDOC) + this%JL(i)*(1-bL)
          else  
            ! Check for nutrient limitation:
            if ( (dSi(i).ge.1.d0) .or. (dN(i).ge.1.d0) ) then
              if (this%JSi(i) .gt. this%JN(i) ) then
                JNet(i) = this%JN(i) ! N limitation
              else  
                JNet(i) = this%JSi(i) ! Si limitation
              end if
            end if
          end if 
          !
          ! Synthesis limitation:
          !
          if ( Jnet(i) .gt. this%JlossPassive(i) ) then ! Apply FR only if net growth is positive
            Jnet(i) = JmaxT * Jnet(i) / ( Jnet(i) + JmaxT )
          endif
          this%Jtot(i) = Jnet(i) - this%JlossPassive(i)
          ! Calculate synthesis-limited uptakes:
          this%JNreal(i) = max( 0.d0, JNet(i) )
          this%JSireal(i) = max( 0.d0, JNet(i) )
          !
          ! Regulate carbon uptakes for growth + respiration towards lowered jNet.
          !

          ! Prioritize DOC:
          tmp = Jnet(i) + bg*max(0.d0, Jnet(i)) + bN*this%JNreal(i) + bSi*this%JSireal(i) + ftemp2*this%Jresp(i)
          this%jDOCreal(i) = min( this%JDOC(i), tmp/(1-bDOC) )
          ! And then light:
          this%JLreal(i) = min( this%JL(i), (tmp - this%jDOCreal(i)*(1-bDOC))/(1-bL) )
          !
          ! Losses:
          !
          this%JNlossLiebig(i) = max( 0.d0, -Jnet(i) )  ! Exude surplus N and Si
          this%JClossLiebig(i) = 0.d0 ! There are never surplus C uptakes
          this%JCloss_photouptake(i) = (1.-epsilonL)/epsilonL * this%JLreal(i)
          this%Jresptot(i)= &
            fTemp2*this%Jresp(i) + &
            bDOC*this%JDOCreal(i) + &
            bL*this%JLreal(i) + &
            bN*this%JNreal(i) + &
            bSi*this%JSireal(i) + &
            max(0.d0, bg*Jnet(i))
          
          this%f(i) = 0.d0
         end do

        ! Needed to get the right rates with getRates:
        this%jN = this%jNreal  
        this%JSi = this%JSireal
        this%jDOC = this%jDOCreal
        this%JNloss = this%JNlossLiebig
     end subroutine calcRatesDiatoms
   
     subroutine calcDerivativesDiatoms(this, u, dNdt, dDOCdt, dSidt, dudt)
       class(spectrumDiatoms), intent(inout):: this
       real(dp), intent(in):: u(this%n)
       real(dp), intent(inout):: dNdt, dDOCdt,dSidt, dudt( this%n )
       integer:: i
   
       this%mort2 = this%mort2constant*u
       this%jPOM = (1-remin2)*this%mort2 ! non-remineralized mort2 => POM

       do i = 1, this%n
        ! mortloss = u(i)*(remin2*this%mort2(i) +reminHTL* this%mortHTL(i))
         !
         ! Update nitrogen:
         !
         dNdt = dNdt  &
         + ((-this%JNreal(i) &!+ (this%JN(i)- this%JlossPassive(i) -this%Jtot(i)) &
         + this%JlossPassive(i) &
         + this%JNlossLiebig(i))/this%m(i) & ! Losses when growth is negative
         + remin2*this%mort2(i) &
           )* u(i)/rhoCN
         !
         ! Update DOC:
         !
         dDOCdt = dDOCdt &
         + ((-this%JDOCreal(i) &
         +  this%JlossPassive(i) &
         +  this%JCloss_photouptake(i))/this%m(i) &
         +  remin2*this%mort2(i) &
           ) * u(i)
         ! update Si
         dSidt = dSidt &
         + ((-this%JSireal(i) &
         +   this%JlossPassive(i) &
         +   this%JNlossLiebig(i) )/this%m(i) & ! Losses when growth is negative
         +  remin2*this%mort2(i) &
           ) * u(i)/rhoCSi
         !
         ! Update the diatoms:
         !
         dudt(i) = &
           (this%Jtot(i)/this%m(i)  &
         - this%mortpred(i) &
         - this%mort2(i) &
         - this%mortHTL(i) &
           )*u(i)
     end do
    end subroutine calcDerivativesDiatoms
      
    subroutine printRatesDiatoms(this)
      class(spectrumDiatoms), intent(in):: this
 
      write(*,*) "Diatoms with ", this%n, " size classes:"
      call this%printRatesUnicellular()
 
      99 format (a10, 20f10.6)
  
      write(*,99) "jSi:", this%JSi / this%m
      write(*,99) "jSireal:", this%JSireal / this%m
      write(*,99) "jResptot:", this%Jresptot / this%m
    end subroutine printRatesDiatoms
    
!
  ! Returns the net primary production calculated as the total amount of carbon fixed
  ! by photsynthesis minus the respiration due to basal respiration,
  ! photosynthesis, nutrint uptake, and growth. Units: mugC/day/m3
  ! (See Andersen and Visser (2023) table 5)
  !
    function getProdNetDiatoms(this, u) result(ProdNet)
      real(dp):: ProdNet
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: u(this%n)
      integer:: i
      real(dp):: resp, tmp
    
      ProdNet = 0.d0

      do i = 1, this%n
        if ( (this%JLreal(i) + this%JDOCreal(i)) .ne. 0.d0 ) then
          tmp = this%JLreal(i) / (this%JLreal(i) + this%JDOCreal(i))
        else
          tmp = 0.d0
        endif

        resp = &
          fTemp2*this%Jresp(i) + & ! Basal metabolism
          bL*this%JLreal(i) + &    ! Light uptake metabolism
          bN*this%JNreal(i) * tmp + &  ! The fraction of N uptake that is not associated to DOC uptake  
          bg*this%Jnet(i) * tmp ! The fraction of growth not associated with DOC
        ProdNet = ProdNet + max( 0.d0, (this%JLreal(i) - resp) * u(i)/this%m(i) )

      end do
    end function getProdNetDiatoms

    function getProdBactDiatoms(this, u) result(ProdBact)
      real(dp):: ProdBact
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: u(this%n)
      integer:: i
  
      ProdBact = 0.d0
      do i = 1, this%n
        ProdBact = ProdBact + max(0.d0, this%JDOC(i) - ftemp2*this%Jresp(i))*u(i)/this%m(i)
      enddo
  
    end function getProdBactDiatoms
   end module diatoms
