!
! Module to handle diatoms
! AMALIA latest version
! with Down-regulation 
!
module diatoms
     use globals
     use spectrum
     use input
     implicit none
 
     private
     !
     ! Stoichiometry:
     !
     !real(dp), parameter:: rhoCN = 5.68
     real(dp) :: rhoCSi
     !cell properties
     !real(dp), parameter:: cm = 6d-7 ! rhoCmem 6*d.-7 ! Carbon content in cell membrane mugC(mum^-3)
     !real(dp), parameter:: cb = 1.25d-7 ! rhoCcyt 1.25*d.-7 ! Carbon content in cell cytoplasm mugC(mum^-3)  
     real(dp) :: v ! Vacuole fraction
     !
     ! Light uptake:
     !
     real(dp) :: epsilonL ! Light uptake efficiency
     real(dp) :: alphaL
     real(dp) :: rLstar
     real(dp) :: bL ! cost of light harvesting mugC(mugC)^-1

     !real(dp), parameter :: cL=0.18 ! photosynthetic affinity constant
     !
     ! Dissolved nutrient uptake:
     !
     real(dp) :: alphaN ! L/d/mugC/mum^2
     real(dp) :: rNstar ! mum
     real(dp) :: bN ! cost of N uptake mugC(mugN)^-1
     real(dp) :: bDOC ! cost of DOC uptake mugC(mugN)^-1
     real(dp) :: bSi ! cost of Si uptake mugC(mugSi)^-1
     !
     ! Metabolism
     !
     real(dp) :: cLeakage ! passive leakage of C and N
     !real(dp) :: c ! Parameter for cell wall fraction of mass.
     real(dp) :: delta ! Thickness of cell wall in mum
     real(dp) :: alphaJ ! Constant for jmax.  per day
     real(dp) :: cR
     real(dp) :: bg ! cost of biosynthsesis -- parameter from literature pending

     !
     ! Bio-geo:
     !
     real(dp) :: remin2 ! fraction of virulysis remineralized to N and DOC
     !real(dp) :: reminHTL ! fraction of HTL mortality remineralized
     real(dp) :: mMinDiatom
     real(dp) :: mMaxDiatom
     
     type, extends(spectrumUnicellular) :: spectrumDiatoms
     real(dp), dimension(:), allocatable:: JSi,JSireal

     contains
       procedure, pass :: initDiatoms
       procedure :: calcRates => calcRatesDiatoms
       procedure :: calcDerivativesDiatoms
       procedure :: printRates => printRatesDiatoms
       procedure :: getNbalanceDiatoms
       procedure :: getCbalanceDiatoms
       procedure :: getSibalanceDiatoms
     end type spectrumDiatoms

     public  initDiatoms, spectrumDiatoms, calcRatesDiatoms, calcDerivativesDiatoms
     public printRatesDiatoms, getNbalanceDiatoms, getCbalanceDiatoms, getSiBalanceDiatoms

   contains
       
     subroutine read_namelist()
        integer :: file_unit,io_err

        namelist /input_diatoms / &
             & RhoCSi, v, &
             & epsilonL, alphaL, rLstar, bL, &
             & alphaN,rNstar, bN, bDOC, &
             & bSi, &
             & cLeakage, delta, alphaJ, cR, bg, &
             & remin2,mMinDiatom, mMaxDiatom


        call open_inputfile(file_unit, io_err)
        read(file_unit, nml=input_diatoms, iostat=io_err)
        call close_inputfile(file_unit, io_err)
     end subroutine read_namelist

     subroutine initDiatoms(this, n)
       class(spectrumDiatoms):: this
       integer, intent(in):: n
       integer:: i
       real(dp), parameter:: mMin = 3.1623d-9
       real(dp), parameter:: rho = 0.57*1d-6
       !real(dp) :: fl

       call read_namelist()
       call this%initUnicellular(n, mMinDiatom, mMaxDiatom)
       allocate(this%JSi(this%n))
       allocate(this%JSireal(this%n))

       !
       ! Radius:
       !
       this%r = (threequarters/pi * this%m/rho/(1-v))**onethird  ! Andy's approximation
      ! this%nu = 3*delta/this%r
      this%nu = 6**twothirds*pi**onethird*delta * (this%m/rho)**(-onethird) * &
        (v**twothirds + (1.+v)**twothirds)
       do i = 1,this%n
        this%nu(i) = min(1.d0, this%nu(i))
       enddo
       
       this%AN = alphaN * this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m
       this%AL = alphaL/this%r * (1-exp(-this%r/rLstar)) * this%m * (1.d0-this%nu)
       this%AF = 0.d0
       this%JFmax = 0.d0
 
       this%JlossPassive = cLeakage/this%r * this%m ! in units of C
 
       this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
       this%Jresp = 0.3*cR*alphaJ*this%m ! decrease suggested by Ken
   
       this%beta = 0.d0 ! No feeding
       this%palatability = 0.5d0 ! Lower risk of predation
     end subroutine initDiatoms
  
     subroutine calcRatesDiatoms(this, L, N,DOC, Si, gammaN, gammaDOC, gammaSi)
       class(spectrumDiatoms), intent(inout):: this
       real(dp), intent(in):: gammaN, gammaSi, gammaDOC
       real(dp), intent(in):: L, N, Si, DOC
       real(dp):: f, JmaxT
       real(dp):: dL(this%n), dN(this%n), dDOC(this%n), dSi(this%n), Jnetp(this%n), Jnet(this%n),Jlim(this%n)
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

         ! write(*,*) this%JL(i)
          jlim(i)=0.0
          !
          ! Potential net uptake
          Jnetp(i)=this%JL(i)*(1-bL)+this%JDOC(i)*(1-bDOC)-ftemp2*this%Jresp(i)
          !Jnet(i)  = max(0.,1./(1+bg)*(jnetp(i)-(bN*dN(i)*this%jN(i)+bSi*dSi(i)*this%JSi(i))))
          !
          ! Calculation of down-regulation factors for
          ! N-uptake
          if (this%JN(i).gt.0) then
            dN(i) = min(1., 1./this%JN(i)*Jnetp(i)/(1+bg +bN+bSi))
          else 
            dN(i)=1.
          end if  
          ! Si-uptake
          if (this%JSi(i).gt.0.) then
            dSi(i) = min(1., 1./this%JSi(i)*Jnetp(i)/(1+bg +bN+bSi))
          else 
            dSi(i)=1.
          end if 
          dDOC(i)=1.
          !------------ up to this point all good -----------------
          !
          !write(*,*) dN(i),dSi(i)

           if (dN(i)==1 .and. dSi(i)==1) then
             if ( this%JN(i)<this%JSi(i) ) then
                 jlim(i)=this%JN(i)
                if (this%JSi(i).gt.0.) then
                 dSi(i)= jlim(i)/this%JSi(i)
                end if
                ! write(*,*) i, 'if (dN(i)==1 .and. dSi(i)==1)'
             else if ( this%JSi(i)<=this%JN(i) ) then
                 jlim(i)= this%JSi(i)
                if (this%JN(i).gt.0.) then
                 dN(i)= jlim(i)/this%JN(i) 
                end if 
               !  write(*,*) i,'( this%JSi(i)<this%JN(i) )'
             end if
           end if   
           
           if  ( dSi(i)==1 .and. dN(i)<1 ) then       ! If Si is the limiting resource
              jlim(i)=this%JSi(i)
              if (this%JN(i).gt.0.) then
               dN(i)= jlim(i)/this%jN(i)
              end if
              !write(*,*) i,'( dSi(i)==1 .and. dN(i)<1 ) '
           else if ( dN(i)==1 .and. dSi(i)<1 ) then   ! If N is the limiting resource
              jlim(i)=this%JN(i)
              if (this%JSi(i).gt.0.) then
               dSi(i)= jlim(i)/this%JSi(i)
              end if
           end if
           !
           if  ( dSi(i)<1 .and. dN(i)<1 ) then  ! this check might be redundant, but just to make sure
              dL(i)=1
             Jnetp(i) = dL(i)*this%jL(i)*(1-bL)+dDOC(i)*this%JDOC(i)*(1-bDOC)-ftemp2*this%Jresp(i)
             Jnet(i)  = max(0.,1./(1+bg)*(jnetp(i)-(bN*dN(i)*this%jN(i)+bSi*dSi(i)*this%JSi(i))))
            jlim(i)=Jnet(i)
           end if
           if (this%JL(i)>0.) then
             dL(i) = min(1., 1./(this%jL(i)*(1-bL))*(jlim(i)*(1+bg+bSi+bN)-(1-bDOC)*this%JDOC(i)+ftemp2*this%Jresp(i)) )


             !write(*,*) 1./(this%jL(i)*(1-bL))*(jlim(i)*(1+bg+bSi+bN)-(1-bDOC)*this%JDOC(i)+ftemp2*this%Jresp(i)) 
           
           else 
             dL(i)=-1.
           end if
           !write(*,*) jlim(i)/this%m(i),dL(i)
           !
           if (dL(i)<0.) then ! It means that there's is too much C taken up 
             dL(i) = 0       ! We set light uptake to 0 and downregulate DOC uptake
             dDOC(i) = min(1., 1./(this%JDOC(i)*(1-bDOC))*(Jlim(i)*(1+bg+bSi+bN) + ftemp2*this%Jresp(i)))  
           end if
             Jnetp(i) = dL(i)*this%jL(i)*(1-bL)+dDOC(i)*this%JDOC(i)*(1-bDOC)-ftemp2*this%Jresp(i)
             !Jnet(i)  = max(0.,1./(1+bg)*(jnetp(i)-(bN*dN(i)*this%jN(i)+bSi*dSi(i)*this%JSi(i))))
             Jnet(i)  = 1./(1+bg)*(jnetp(i)-(bN*dN(i)*this%jN(i)+bSi*dSi(i)*this%JSi(i)))

             !write(*,*) i, Jnetp(i)/this%m(i),dL(i)
             if ( dSi(i)<1 .and. dN(i)<1 ) then
              jlim(i)=Jnet(i)
              if (this%JN(i).gt.0.) then
              dN(i)=jlim(i)/this%JN(i)
              end if
              if (this%JSi(i).gt.0.) then
              dSi(i)=jlim(i)/this%JSi(i)
              end if

             end if
          !
          ! Down-regulated fluxes: the uptake flux JX is multiplied by the reduction factor dX , X= N, DOC, Si or L
          !
          ! update jnetp with delta's
          !
          !Jnetp  = this%JL(i)*(1-bL)+dDOC(i)*this%jDOC(i)*(1-bDOC)-this%Jresp(i)   
          !Jnet = max(0., 1./(1+bg)*(Jnetp - (bN*dN(i)*this%JN(i)+bSi*dSi(i)*this%JSi(i)))) 
          !
          ! Saturation of net growth
          !
          f = (Jnet(i))/(Jnet(i)+ JmaxT)
          if ((Jnet(i) + JmaxT).eq.0) then
           f=0.
          end if

          this%JNreal(i) = (1-f)*dN(i) * this%JN(i) 
          this%JDOCreal(i) = (1-f)*dDOC(i) * this%JDOC(i)
          this%JLreal(i) = (1-f)*dL(i) * this%JL(i)
          this%JSireal(i) = (1-f)*dSi(i) * this%JSi(i)
          this%Jtot(i)= f * JmaxT - (1-f)*( this%JlossPassive(i) )!+ fTemp2*this%Jresp(i))


         this%JCtot(i) = & 
        !
        ! (1-f)*(dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i) )&
        ! - (1-f)*fTemp2*this%Jresp(i) &
        ! - ( (1-f)*(bDOC*dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i)*bL+bN*dN(i)*this%JN(i)+bN*dSi(i)*this%JSi(i))&
        ! + (1-f)*bg*Jnet(i) )
          (1-f)*(dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i)- fTemp2*this%Jresp(i) &
         - bN*dN(i)*this%JN(i)-bSi*dSi(i)*this%JSi(i) -bg*Jnet(i) -bDOC*dDOC(i)*this%JDOC(i)-dL(i)*this%JL(i)*bL)
             
          this%JCloss_photouptake(i) = (1.-epsilonL)/epsilonL * this%JLreal(i)
          this%Jresptot(i)= (1-f)*fTemp2*this%Jresp(i) + & !
               (1-f)*(bDOC*dDOC(i)*this%JDOC(i) + &
                      bL*dL(i)*this%JL(i) + &
                      bN*dN(i)*this%JN(i) + &
                      bSi*dSi(i)*this%JSi(i) + &
                      bg*Jnet(i))
          !
          !write(*,*) jlim(i),dN(i),dSi(i)

          !write(*,*) 'N budget', i,':',(this%JNreal(i)-(1-f)*this%JlossPassive(i) &
          !- this%Jtot(i))/this%m(i)
          !
          !write(*,*) 'C budget', i,':',(this%JCtot(i) -(1-f)*this%JlossPassive(i)&
          !- this%Jtot(i))/this%m(i) !this works only if we take the negative values of jnet

         !write(*,*) 'Si budget', i,':',(this%JSi(i)-this%JlossPassive(i) - (this%JSi(i)-this%JlossPassive(i)-this%Jtot(i)) &
         !- this%Jtot(i))/this%m(i)
         this%f(i)=f
         !write(*,*) this%JCloss_photouptake(i), this%JLreal(i),dL(i)
        end do
        this%jN = this%jNreal  ! Needed to get the the actual uptakes out with "getRates"
        this%jDOC = this%jDOCreal
        this%JSi = this%jSireal
     end subroutine calcRatesDiatoms
   
     subroutine calcDerivativesDiatoms(this, u, dNdt, dDOCdt, dSidt, dudt)
       class(spectrumDiatoms), intent(inout):: this
       real(dp), intent(in):: u(this%n)
       real(dp), intent(inout):: dNdt, dDOCdt,dSidt, dudt( this%n )
       real(dp):: mortloss
       integer:: i
   
       this%mort2 = this%mort2constant*u
       this%jPOM = 0*(1-remin2)*this%mort2 ! non-remineralized mort2 => POM

       do i = 1, this%n
        ! mortloss = u(i)*(remin2*this%mort2(i) +reminHTL* this%mortHTL(i))
         !
         ! Update nitrogen:
         !
         dNdt = dNdt  &
         + ((-this%JNreal(i) &!+ (this%JN(i)- this%JlossPassive(i) -this%Jtot(i)) &
         + (1-this%f(i))*this%JlossPassive(i))/this%m(i) &
         + remin2*this%mort2(i) &
         !+ reminHTL*this%mortHTL(i) &
         )* u(i)/rhoCN
         !
         ! Update DOC:
         !
         dDOCdt = dDOCdt &
         + ((-this%JDOCreal(i) &
         +  (1-this%f(i))*this%JlossPassive(i) &
         +  this%JCloss_photouptake(i))/this%m(i) &
         +  remin2*this%mort2(i) &
         !+  reminHTL*this%mortHTL(i) &
         ) * u(i)
         ! update Si
         dSidt = dSidt &
         + ((-this%JSireal(i) &
         +   (1-this%f(i))*this%JlossPassive(i) )/this%m(i) &
         +  remin2*this%mort2(i) &
         !+  reminHTL*this%mortHTL(i) &
         ) * u(i)/rhoCSi
         !
         ! Update the diatoms:
         !
         dudt(i) = (this%Jtot(i)/this%m(i)  &
         !- mort(i) &
         - this%mortpred(i) &
         - this%mort2(i) &
         - this%mortHTL(i))*u(i)

         ! write(*,*) this%mortpred(i)*u(i) 


         !write(*,*) 'Nbalance =',i, ':', (dNdt + sum( dudt &
         !+ (1-reminHTL)*this%mortHTL*u &
         !+ (1-1)*this%mort2*u & ! full N remineralization of viral mortality
         !   )/rhoCN) 

            !write(*,*) 'Nbalance =',i, ':', (dNdt + sum( dudt(i:this%n) &
            !+ (1-reminHTL)*this%mortHTL(i:this%n)*u(i:this%n) &
            !+ (1-1)*this%mort2(i:this%n)*u(i:this%n) & ! full N remineralization of viral mortality
             !  )/rhoCN) 
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

    function getNbalanceDiatoms(this, N, dNdt, u, dudt) result(Nbalance)
      real(dp):: Nbalance
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: N,dNdt, u(this%n), dudt(this%n)
  
      Nbalance = (dNdt + sum( dudt &
      + (1-fracHTL_to_N)*this%mortHTL*u &
      + (1-remin2)*this%mort2*u & ! full N remineralization of viral mortality
         )/rhoCN)/N 
    end function getNbalanceDiatoms 

    function getSibalanceDiatoms(this, Si, dSidt, u, dudt) result(Sibalance)
      real(dp):: Sibalance
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: Si,dSidt, u(this%n), dudt(this%n)
  
      Sibalance = (dSidt + sum( dudt &
      + this%mortHTL*u &
      + (1-remin2)*this%mort2*u & ! full Si remineralization of viral mortality
      )/rhoCSi)/Si 
    end function getSibalanceDiatoms 
  
    function getCbalanceDiatoms(this, DOC, dDOCdt, u, dudt) result(Cbalance)
      real(dp):: Cbalance
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: DOC, dDOCdt, u(this%n), dudt(this%n)
  
      Cbalance = (dDOCdt + sum(dudt &
      + this%mortHTL*u &
      + (1-remin2)*this%mort2*u &
      - this%JLreal*u/this%m & ! ??
      - this%JCloss_photouptake*u/this%m &
      + this%Jresptot*u/this%m & !plus uptake costs
      )) / DOC
    end function getCbalanceDiatoms 
    
  
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