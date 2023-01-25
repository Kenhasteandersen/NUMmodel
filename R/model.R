# --------------------------------------------------
# Core logic for the model
# --------------------------------------------------


parameters <- function(n=25, mmin10=-8.5, mmax10=0) {
  p = list()
  
  p$n = as.integer(n); # No of groups
  p$m = 10^seq(mmin10,mmax10,length.out = p$n)  # Mass bins in mugC
  rho = 0.4*1e6*1e-12 # mug/cm3 (Andersen et al 2016; rho = m/V = 0.3*d^3/(4/3*pi*(d/2)^3) )
  p$r = calcRadius(p$m)
  p$Delta = p$m[2]/p$m[1]
  
  p$rhoCN = 5.68 # C:N mass ratio
  p$epsilonL = 0.8 # Light uptake efficiency
  p$epsilonF = 0.8 # Assimilation efficiency
  p$cLeakage = 0.03 # passive leakage of C and N
  #
  # Cell wall fraction of mass:
  #
  p$c = 0.0015 # the constant is increased a bit to limit the lower cell size
  p$delta = 0.05 # Thickness of cell wall in mum
  p$nu = 3*p$delta/p$r 
  #p$nu=p$c * p$m^(-1/3)  
  #
  # Clearance rates:
  #
  #factor = (1e-6)^(1/3)/1.5 
  p$aN = 0.972 #0.682 #  L/d/mugC/mum^2.   0.00012 #0.00004 # 0.000162 # Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g)^(-1/3) / 1.5 (g/cm); Andersen et al 2015
  p$rNstar = .4 # mum
  #p$AL = 0.0019 # if using Al propto m^(2/3) for non-diatoms
  #p$AL = 0.0012 # if using my shading formula for non-diatoms
  #p$cL = 0.021 # if using my shading formula for non-diatoms
  #p$AL = 0.000914 # if using Andys shading formula for non-diatoms
  #p$cL = 21 # if using Andys shading formula for non-diatoms
  p$alphaL = 0.3 # 0.206
  p$rLstar = 7.5
  p$aF = 0.018  #  Fits to TK data for protists
  p$cF = 30# Just a guess
  #
  # Calc rates as a function of m:
  #
  #p$ANm = p$AN*p$m^(1/3) / (1 + p$cN*p$m^(-1/3))
  #p$ANm = 0.00012*p$m^(1/3) / (1 + 0.1*p$m^(-1/3))
  #p$ALm = p$AL*p$m^(2/3)*(1-nu)
  #p$ALm = p$cL*p$m * p$AL*p$m^(2/3) / ( p$cL*p$m + p$AL*p$m^(2/3) )  # shading formula
  #p$ALm = p$AL*p$m^(2/3) * (1-exp(- p$cL*p$m^(1/3) ))  # shading formula
  p$aNm = p$aN*p$r^(-2) / (1 + (p$r/p$rNstar)^(-2))
  p$aLm = p$alphaL/p$r * (1-exp(-p$r/p$rLstar))
  p$aFm = p$aF*p$m/p$m
  p$jFmaxm = p$cF/p$r

  p$Jloss_passive_m = p$cLeakage / p$r *p$m # in units of C
  #
  # Prey encounter
  #
  p$beta = 500
  p$sigma = 1.3
  p$theta = matrix(nrow=p$n, ncol=p$n)
  for (i in 1:p$n) {
    #p$theta[i,] = phi(p$m[i]/p$m, p$beta, p$sigma)   # (predator, prey)
    for (j in 1:p$n)
      p$theta[i,j] = Phi(p$m[i]/p$m[j], p)
  }
  #
  # Metabolism:
  #
  p$alphaJ = 1.5 # per day
  p$jMax = pmax(0, p$alphaJ * (1-p$nu)) # mugC/day
  p$cR = 0.1
  p$jR = p$cR*p$alphaJ
  #
  # Losses:
  #
  p$mort = 0*0.005*(p$jMax) * p$m^(-1/4)
  p$mort2 = 0.004/log(p$Delta) # 
  p$mortHTL = 0.1
  p$mHTL = max(p$m)/p$beta^1.5 # Bins affected by HTL mortality
  p$mortHTLm = 1/(1+(p$m/p$mHTL)^(-2))
  #
  # Production of DOC:
  #
  p$remin = 0.0 # fraction of mortality losses reminerilized to N and DOC
  p$remin2 = 0.5 # fraction of virulysis being available to N and DOC
  p$reminF = 0.1 # Fraction of feeding losses being directly available
  p$reminHTL = 0
  
  p$T = 10
  p$latitude=0
  p$L = 100
  #
  # Initial conditions:
  #
  p$N0 = 150
  p$DOC0 = 0
  p$B0 = rep(1,p$n)
  
  return(p)
}

calcRadius = function(m) {
  rho = 0.4*1e6*1e-12 # mug/cm3 Menden-Deuer (2000). If we assume that m propto V we get approximately rho = 0.4e-6 ugC/um^3.
  return( (3/(4*pi)*m/rho)^(1/3) )
}

calcMass = function(r) {
  rho = 0.4*1e6*1e-12 # mug/cm3 (Menden-Deuer (2000). If we assume that m propto V we get approximately rho = 0.4e-6 ugC/um^3.
  return( r^3*rho*4*pi/3)
}

#
# Prey size preference function:
#
phi = function(z, beta, sigma) {
  exp( -(log(z/beta))^2/(2*sigma^2) )
}
#
# Prey size function integrated over size groups:
#
Phi = function(z, p=parameters(), Delta=p$m[2]/p$m[1]) {
  m = 10^seq(-12,3,length.out = 1000)
  dm = diff(m)
  m = m[2:1000]
  
  fPrey = function(m, w0, Delta) {
    integrate( function(logw) phi(m/exp(logw),beta=p$beta,sigma=p$sigma), log(w0/sqrt(Delta)), log(w0*sqrt(Delta)))$value
  }
  
  fTot = function(m0,w0, Delta) {
    ix = m>m0/sqrt(Delta) & m<m0*sqrt(Delta)
    sum(
      as.numeric(
        lapply( m[ix], function(m) fPrey(m,w0,Delta)/m) 
      )*dm[ix] / log(Delta)^2
    )
  }

  return(fTot(1,1/z,Delta))
}

#
# Q10 temperature function:
#
fTemp = function(Q10, T, Tref=10) {
  return(Q10^((T-Tref)/10))
}
#
# Convert to ESD:
#
calcESD = function(m) {
  10000 * 1.5 * (m*1e-6)^(1/3)
}
#
# Version 3: heuristic down-regulation again
#
calcRates = function(Light,N,DOC,B,p) {
  with(p, {
    B = pmax(0,B)
    N = max(0,N)
    DOC = max(0, DOC)
    
    #print(c(N,DOC,B))
    #if (sum(is.nan(c(B,N,DOC))))
    #  browser()
    #
    # Temperature corrections:
    #
    #ANmT = ANm*fTemp(1.5,p$T)
    ANmT = aNm*m*fTemp(1.5,p$T)
    f2 = fTemp(2,p$T)
    JmaxT = jMax*f2*m
    JR = jR*f2*m
    JFmaxmT = jFmaxm*m*f2
    #
    # Uptakes
    #
    JN =   ANmT*N*rhoCN # Diffusive nutrient uptake in units of C/time
    # 
    JDOC = ANmT*DOC # Diffusive DOC uptake, units of C/time
    
    JL =   epsilonL * aLm*m*Light  # Photoharvesting

    # Light acclimation:
    #JLreal = pmax( 0, JL - pmax(0,(JL+JDOC-JR - JN)) )
    
    # Feeding as type II with surface limitation:
    F = theta %*% B
    JF = epsilonF * JFmaxmT * aFm*p$m*F / (aFm*p$m*F + JFmaxmT) #        # Feeding
    
    # Passive losses:
    #Jloss_passive = p$cLeakage * m^(2/3) # in units of C
    
    # Total nitrogen uptake:    
    JNtot = JN+JF-Jloss_passive_m # In units of C
    # Total carbon uptake
    JCtot = JL+JF+JDOC-JR-Jloss_passive_m
    
    # Liebig + synthesis limitation:
    #Jtot = pmin( JNtot, JCtot, JmaxT )
    Jtot = pmin( JNtot, JCtot ) 
    f = Jtot/(Jtot + JmaxT)
    
    # If synthesis-limited then down-regulate feeding:
    #JFreal = pmax(0, JF - pmax( 0,  pmax(0, Jtot-JmaxT) ))
    JFreal = pmax(0, JF - (Jtot-f*JmaxT)*(Jtot>0))
    Jtot = f*JmaxT
    JLreal = JL-pmax(0, pmin((JCtot - (JF-JFreal)-Jtot), JL))
    
    # Actual uptakes:
    JCtot = JLreal + JDOC + JFreal - JR - Jloss_passive_m
    JNtot = JN + JFreal - Jloss_passive_m
    # 
    # Losses:
    #
    JCloss_feeding = (1-epsilonF)/epsilonF*JFreal # Incomplete feeding (units of carbon per time)
    JCloss_photouptake = (1-epsilonL)/epsilonL*JLreal
    JNlossLiebig = pmax(0,JNtot-Jtot)  # In units of C
    JClossLiebig = pmax(0,JCtot-Jtot) # C losses from Liebig, not counting losses from photoharvesting

    JNloss = JCloss_feeding + JNlossLiebig + Jloss_passive_m # In units of C
    JCloss = JCloss_feeding + JCloss_photouptake + JClossLiebig + Jloss_passive_m
    #
    # Check:
    #
    #print( ((JLreal + JFreal + JDOC - JR - Jloss_passive - JClossLiebig) - Jtot )/m )
    #print( ((JN + JFreal - Jloss_passive - JNlossLiebig) - Jtot )/m )
    
    #if (sum(c(JNloss,JCloss,B)<0))
    #  browser()
    #
    # Mortality:
    #
    mortpred =  t(theta) %*% (JFreal/epsilonF*B/m/(F+1e-20))

    #if (sum(is.nan(Jtot)))
    #  browser()
    
    
    return(list( 
      jN=JN/m, 
      jDOC=JDOC/m, 
      jL=JL/m, 
      jF=JF/m,
      jFreal = JFreal/m,
      jLreal = JLreal/m, 
      jR=rep(jR,n),
      jFmax=JFmaxmT/m,
      jMax = JmaxT/m,
      jNlossLiebig=JNlossLiebig/m, 
      jClossLiebig=JClossLiebig/m,
      jCloss_photouptake=JCloss_photouptake/m,
      jLossPassive = Jloss_passive_m/m,
      jCloss_feeding=JCloss_feeding/m, 
      jCloss=JCloss/m, 
      jNloss=JNloss/m,
      jTot=Jtot/m, 
      F=F, 
      jCtot = JCtot/m, 
      jNtot=JNtot/m,
      mortpred=mortpred, 
      mort=mort,
      mort2=mort2*B,
      mortHTL = mortHTL*mortHTLm,
      totKilled = sum(JF/epsilonF*B/m), 
      totEaten = sum(mortpred*B), 
      totGrowth=sum(Jtot*B/m)
      ))
  })
}
#
# Version 2: regulation only through functional responses:
#
calcRatesFR = function(t,L,N,DOC,B,p) {
  with(p, {
    B = pmax(0,B)
    N = max(0,N)
    DOC = max(0, DOC)
    #
    # Temperature corrections:
    #
    ANmT = ANm*fTemp(1.5,p$T)
    JmaxT = jmax*fTemp(2,p$T)
    JR = Jresp*fTemp(2,p$T)
    #
    # Uptakes
    #
    JN =   JmaxT/p$rhoCN * ANmT*N / (JmaxT/p$rhoCN + ANmT*N) # Diffusive nutrient uptake
                                                        # in units of N/time
    JDOC = JmaxT * ANmT*DOC / (JmaxT + ANmT*DOC) # Diffusive DOC uptake, units of C/time
    
    JL =   epsilonL * JmaxT * ALm*L / (JmaxT + ALm*L)  # Photoharvesting
    
    F = theta %*% B
    JF = epsilonF * JmaxT * AFm*F / (JmaxT + AFm*F)        # Feeding

    # Total nitrogen uptake:    
    JNtot = JN+JF/rhoCN # In units of N
    
    # Down-regulation of light uptake:
    JLreal = pmin(JL, pmax(0, JNtot*rhoCN - (JDOC-JR)))
                     
    JCtot = JLreal+JF+JDOC-JR # Total carbon uptake

    Jlim = pmin( JCtot, JNtot*rhoCN )  # Liebigs law; units of C
    # 
    # Losses:
    #
    JCloss_feeding = (1-epsilonF)/epsilonF*JF # Incomplete feeding (units of carbon per time)
    JCloss_photouptake = (1-epsilonL)/epsilonL*JLreal
    JNlossLiebig = pmax(0, JNtot*rhoCN-JCtot)/rhoCN  # N losses from Liebig
    JClossLiebig = pmax(0, JCtot-JNtot*rhoCN) # C losses from Liebig, not counting losses from photoharvesting
    #JClossLiebig = pmax(0, Jtot - JNtot*rhoCN) # C losses from Liebig, not counting losses from photoharvesting
    #JClossLiebig = pmin(JClossLiebig, JDOC) # However, light surplus is not leaked but is downregulated

    Jloss_passive = p$cLeakage * m^(2/3) # in units of C
    
    JNloss = JCloss_feeding/rhoCN + JNlossLiebig + Jloss_passive/rhoCN
    JCloss = JCloss_feeding + JCloss_photouptake + JClossLiebig + Jloss_passive
    
    #if (sum(c(JNloss,JCloss,B)<0))
    #  browser()
    #
    # Mortality:
    #
    mortpred =  t(theta) %*% (JF/epsilonF*B/m/F)

    Jtot = Jlim - Jloss_passive
    
    return(list( 
      JN=JN, JDOC=JDOC, JL=JL, JF=JF,
      JLreal = JLreal, JR=JR,
      JNlossLiebig=JNlossLiebig, JClossLiebig=JClossLiebig,
      JCloss_photouptake=JCloss_photouptake,
      Jloss_passive = Jloss_passive,
      JCloss_feeding=JCloss_feeding, JCloss=JCloss, JNloss=JNloss,
      Jtot=Jtot, F=F, JCtot = JCtot, JNtot=JNtot,
      mortpred=mortpred, mort=mort,
      mort2=mort2*B,
      totKilled = sum(JF/epsilonF*B/m), 
      totEaten = sum(mortpred*B), 
      totGrowth=sum(Jtot*B/m)))  
  })
}
#
# Version 1. Heuristic down-regulation; first of feeding,
#
calcRatesOld = function(t,N,DOC,B,p) {
  with(p, {
    #
    # Potential uptakes:
    #
    JN =   ANm*N  # Diffusive nutrient uptake
    JDOC = ANm*DOC # Diffusive DOC uptake
    JL =   epsilonL*ALm*L*(1+amplitudeL*(-cos(t*2*pi/365)))  # Photoharvesting
    
    F = theta %*% B
    JF = epsilonF*AFm*F        # Feeding
    
    JCtot = JL+JF-Jresp+JDOC # Total carbon intake
    JNtot = JN+JF/rhoCN # In units of N
    Jtot = pmin( JCtot, JNtot*rhoCN )  # Liebigs law; units of C
    
    f = Jtot / (Jtot + Jmax) # feeding level
    #
    # Actual uptakes:
    #
    JFreal = pmax(0, JF - (Jtot-f*Jmax))
    JLreal = JL-pmin((JCtot - (JF-JFreal)-f*Jmax), JL)
    JDOCreal = pmin(JDOC, Jresp + f*Jmax - JLreal - JFreal) # The min is only needed to reduce round-off errors
    JNreal = pmax(0, (f*Jmax - JFreal))/rhoCN # In units of N
    # 
    # Losses:
    #
    JNloss_piss = -pmin(0, (f*Jmax - JFreal))/rhoCN  # Exudation of surplus nutrients by heterotrophs
    JNloss = (1-epsilonF)/epsilonF*JFreal/rhoCN + JNloss_piss
    JCloss_feeding = (1-epsilonF)/epsilonF*JFreal
    JCloss = JCloss_feeding + (1-epsilonL)/epsilonL*JLreal
    
    # These two should be equal to zero:
    #JCcheck = f*Jmax - JLreal - JFreal - JDOCreal + Jresp
    #JNcheck = f*Jmax/rhoCN - JNreal - JFreal/rhoCN + JNloss_feeding
    #print(JNcheck)
    #
    # Mortality:
    #
    mortpred =  t(theta) %*% (JFreal/epsilonF*B/m/F)
    #
    # System:
    #
    dBdt = (Jmax*f/m  - (mort + mortpred + mort2*B + mortHTL*(m>=mHTL)))*B
    dNdt   =  d*(N0-N)  - sum(JNreal*B/m)   + sum(JNloss*B/m)
    dDOCdt =  d*(0-DOC) - sum(JDOCreal*B/m) + sum(JCloss*B/m)

    # Check of nutrient conservation; should be close to zero
    #Nin = d*(N0-N)
    #Nout = (1-remin) * ( sum(mort2*B*B) + sum(mortHTL*(m>=mHTL)*B) ) / rhoCN
    #NcheckSystem = Nin - Nout - sum(dBdt)/rhoCN - dNdt
    #print(NcheckSystem)
    
    return(list( 
      dNdt=dNdt, dDOCdt=dDOCdt, dBdt=dBdt, 
      JN=JN, JDOC=JDOC, JL=JL, JF=JF,
      JNreal=JNreal, JDOCreal=JDOCreal, JLreal=JLreal, JFreal=JFreal, 
      JNloss_piss=JNloss_piss, JNloss=JNloss, 
      JCloss_feeding=JCloss_feeding, JCloss=JCloss,
      Jtot=Jtot, f=f, F=F, JCtot = JCtot, JNtot=JNtot,
      mortpred=mortpred, mort=mort,
      mort2=mort2*B,
      totKilled = sum(JFreal/epsilonF*B/m), totEaten = sum(mortpred*B), totGrowth=sum(Jmax*f*B/m)))  
  })
}
