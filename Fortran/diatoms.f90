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
      ! this%nu = 3*delta/this%r
      this%nu = 6**twothirds*pi**onethird*delta * (this%m/rho)**(-onethird) * &
        (v**twothirds + (1.+v)**twothirds)
       do i = 1,this%n
        this%nu(i) = min(1.d0, this%nu(i))
       enddo
       
       ! Affinities are the same as for the generalists divided by a factor (1-v):
       this%AN = alphaN * this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m / (1-v)
       this%AL = alphaL/this%r * (1-exp(-this%r/rLstar)) * this%m * (1.d0-this%nu) / (1-v)
       this%AF = 0.d0
       this%JFmax = 0.d0
 
       this%JlossPassive = cLeakage/this%r * this%m ! in units of C
 
       this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
       this%Jresp = cR*alphaJ*this%m ! decrease suggested by Ken
   
       this%beta = 0.d0 ! No feeding
       this%palatability = palatability ! Lower risk of predation
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
          Jnetp(i) = this%JL(i)*(1-bL)+this%JDOC(i)*(1-bDOC)-ftemp2*this%Jresp(i)
          if (Jnetp(i) .lt. 0.d0) then ! There is not enough carbon to satisfy standard metabolism. Then only take up carbon and no resources
             dN(i) = 0.
             dSi(i) = 0.
             dDOC(i) = 1.
             dL(i) = 1.
             Jnet(i) = 0.
          else
           !Jnet(i)  = max(0.,1./(1+bg)*(jnetp(i)-(bN*dN(i)*this%jN(i)+bSi*dSi(i)*this%JSi(i))))
          !
          ! Calculation of down-regulation factors for
          ! N-uptake
          if (this%JN(i).gt.0) then
            dN(i) = min(1.d0, 1.d0/this%JN(i)*Jnetp(i)/(1.d0+bg +bN+bSi))
          else 
            dN(i)=1.
          end if  
          ! Si-uptake
          if (this%JSi(i).gt.0.) then
            dSi(i) = min(1.d0, 1.d0/this%JSi(i)*Jnetp(i)/(1.d0+bg +bN+bSi))
          else 
            dSi(i)=1.
          end if 
          dDOC(i)=1.
          !------------ up to this point all good -----------------
          !
          !write(*,*) dN(i),dSi(i), Jnetp(i), jlim(i)

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
             Jnet(i)  = max(0.d0,1.d0/(1.d0+bg)*(jnetp(i)-(bN*dN(i)*this%jN(i)+bSi*dSi(i)*this%JSi(i))))
            jlim(i)=Jnet(i)
           end if
           if (this%JL(i)>0.) then
             dL(i) = min(1.d0, 1.d0/(this%jL(i)*(1.d0-bL))*(jlim(i)*(1+bg+bSi+bN) &
             -(1-bDOC)*this%JDOC(i)+ftemp2*this%Jresp(i)) )


             !write(*,*) 1./(this%jL(i)*(1-bL))*(jlim(i)*(1+bg+bSi+bN)-(1-bDOC)*this%JDOC(i)+ftemp2*this%Jresp(i)) 
           
           else 
             dL(i)=-1.
           end if
           !write(*,*) jlim(i)/this%m(i),dL(i)
           !
           if (dL(i)<0.) then ! It means that there's is too much C taken up 
             dL(i) = 0       ! We set light uptake to 0 and downregulate DOC uptake
             dDOC(i) = min(1.d0, 1.d0/(this%JDOC(i)*(1-bDOC))*(Jlim(i)*(1.d0+bg+bSi+bN) + ftemp2*this%Jresp(i)))  
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
            end if   
          !
          ! Saturation of net growth
          !
          f = (Jnet(i))/(Jnet(i)+ JmaxT)
          if ((Jnet(i) + JmaxT).eq.0) then
           f=0.
          end if
      

          !write(*,*) i, dN(i),dDOC(i),dL(i),dSi(i)

          this%JNreal(i) =   (1-f)*dN(i)   * this%JN(i) 
          this%JDOCreal(i) = (1-f)*dDOC(i) * this%JDOC(i)
          this%JLreal(i) =   (1-f)*dL(i)   * this%JL(i)
          this%JSireal(i) =  (1-f)*dSi(i)  * this%JSi(i)
          this%Jtot(i)= f*JmaxT - (1-f)*this%JlossPassive(i)!+ fTemp2*this%Jresp(i))
          this%jNet(i) = Jnet(i)
        ! this%JCtot(i) = & 
        !
        ! (1-f)*(dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i) )&
        ! - (1-f)*fTemp2*this%Jresp(i) &
        ! - ( (1-f)*(bDOC*dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i)*bL+bN*dN(i)*this%JN(i)+bN*dSi(i)*this%JSi(i))&
        ! + (1-f)*bg*Jnet(i) )
        !  (1-f)*(dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i)- fTemp2*this%Jresp(i) &
        ! - bN*dN(i)*this%JN(i)-bSi*dSi(i)*this%JSi(i) -bg*Jnet(i) &
        ! -bDOC*dDOC(i)*this%JDOC(i)-dL(i)*this%JL(i)*bL)
             
          this%JCloss_photouptake(i) = (1.-epsilonL)/epsilonL * this%JLreal(i)
          this%Jresptot(i) = &
              fTemp2*this%Jresp(i) + & !
              bDOC*this%JDOCreal(i) + &
              bL*this%JLreal(i) + &
              bN*this%JNreal(i) + &
              bSi*this%JSireal(i) + &
              bg*Jnet(i)
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
