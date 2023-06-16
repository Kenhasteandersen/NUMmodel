source("model.R")
source("basetools.R")
#library("sundialr")
library("deSolve")
#library(tictoc)

alpha = 0.25
colOsmo = rgb(165/256,42/256,42/256,alpha=alpha)
colPhoto = rgb(0,1,0,alpha=alpha)
colN = rgb(0,0,1,alpha=alpha)
colMixo = rgb(1,0.5,0.5,alpha=alpha)
colHetero = rgb(1,0,0,alpha=alpha)

parametersChemostat = function(p=parameters()) {
  #
  # Biogeochemical model:
  #
  p$d = 0.05  # diffusion rate, m/day
  p$M = 20   # Thickness of the mixed layer, m
  p$T = 10   # Temperature
  p$N0 = 50 # Deep nutrient levels
  
  #
  # Light:
  #
  p$L = 30  # PAR, mu E/m2/s
  p$latitude = 0 # amplitude of seasonal light variation in fractions of L
  
  p$tEnd = 365 # Simulation length (days)
  
  return(p)
}

#
# Seasonal variation in exchange rate as a function of latitude (degrees)
# and time (days)
#
# SeasonalExchange = function(latitude, t) {
#   t = t %% 365
#   
#   dmax = 0.05*(1+tanh(0.05*(latitude-40)))
#   dsummer = 0.01
#   tspring = 180 * latitude/120
#   tautumn = 200 + 180 *(90-latitude)/90
#   widthautumn = 1
#   
#   summer = 1-0.5*(1+tanh(8*((t-tspring)/365*2*pi)))
#   winter = 1-0.5*(1 - tanh(widthautumn*((t-tautumn)/365*2*pi)))
#   spring = 1-0.5*(1 - tanh(widthautumn*(((t+365)-tautumn)/365*2*pi)))
#   summer = pmin(summer, spring)
#   
#   d = dsummer + dmax*(winter + summer)
# }
#
# Seasonal variation in light. Roughly taken from Evans and Parslows. 
# M is the depth of the mixed layer.
#
#SeasonalLight = function(p,t) {
#  p$L*exp(-0.025*p$M)*(1 - 0.8*sin(pi*p$latitude/180)*cos(2*pi*t/365))
#}

# derivative = function(t,y,p) {
#   N = y[1]
#   DOC = y[2]
#   B = y[3:(2+p$n)]
#   
#   if (p$latitude > 0)
#     L = SeasonalLight(p,t)
#   else
#     L = p$L
#   
#   rates = calcRates(L,N,DOC,B,p)
#   #
#   # System:
#   #
#   if (p$latitude>0)
#     diff = SeasonalExchange(p$latitude, t)
#   else
#     diff = p$d
#   
#   dBdt = diff*(0-B) + (rates$jTot  - 
#                          (  rates$mort+ 
#                               rates$mortpred + 
#                               rates$mort2 + 
#                               p$mortHTL*p$mortHTLm))*B
#   #dBdt[(B<1e-3) & (dBdt<0)] = 0 # Impose a minimum concentration even if it means loss of mass balance
#   
#   #mortloss = sum(B*(rates$mort2 + p$mortHTL*(p$m>=p$mHTL)))
#   mortloss = sum(B*((1-p$remin2)*rates$mort2 + p$mortHTL*p$mortHTLm))
#   dNdt   =  diff*(p$N0-N) -
#     sum(rates$jN/p$rhoCN*B) +
#     sum(rates$jNloss/p$rhoCN*B) +
#     p$remin2*sum(rates$mort2*B)/p$rhoCN +
#     p$remin*mortloss/p$rhoCN
#   dDOCdt =  diff*(0-DOC) -
#     sum(rates$jDOC*B) +
#     sum(rates$jCloss*B) +
#     p$remin2*sum(rates$mort2*B) +
#     p$remin*mortloss
#   
#   if (TRUE==FALSE) {
#     # Check of nutrient conservation; should be close to zero
#     Nin = diff*(p$N0-N) + diff*sum(0-B)/p$rhoCN
#     Nout = (1-p$remin) * mortloss / p$rhoCN 
#     NcheckSystem = Nin - Nout - sum(dBdt)/p$rhoCN - dNdt
#     print(c(NcheckSystem, Nin, Nout))
#     
#     # Check carbon balance; should be close to zero:
#     Cin = diff*(0-DOC) + diff*sum(0-B) + sum(rates$JLreal/p$m*B/p$epsilonL)
#     Cout = (1-p$remin) * mortloss + sum(p$Jresp/p$m*B)
#     CcheckSystem = Cin - Cout - sum(dBdt) - dDOCdt
#     #print(CcheckSystem)
#   }
#   #print(sum(rates$JF/p$epsilonF*B/p$m) - sum(rates$mortpred*B))
#   
#   # Check. Expensive to evaluate, so commented out  
#   #  if ( sum(c( is.nan(unlist(rates)), is.infinite(unlist(rates)), rates$N<0, rates$B<0))>0)
#   #    browser()
#   
#   return(c(dNdt, dDOCdt, dBdt))
# }
loadNUMmodel = function() {
  sys=Sys.info()['sysname']
  
  if (sys=='Darwin') 
    sLibname = '../lib/libNUMmodel_R.dylib'
  if (sys=='Linux') 
    sLibname = '../lib/libNUMmodel_linux_R.so'
  if (sys=='Windows')
    sLibname = '../lib/libNUMmodel_R.dll'

  dyn.load(sLibname)
}

derivativeF = function(t,y,p) {
  derivF = .Fortran("f_calcderivatives", 
                    #nGrid=as.integer(length(y)),
                    u = as.numeric(y),
                    L=as.numeric(p$L),
                    T=as.numeric(p$T),
                    dt = as.numeric(0),
                    dudt=as.numeric(y))
  # Do the chemostat:
  derivF$dudt = derivF$dudt - p$d*y
  derivF$dudt[1] = derivF$dudt[1] + p$d*p$N0
  #derivF$dudt[1] = derivF$dudt[1] + p$d*(p$N0-y[1])
  #derivF$dudt[2] = derivF$dudt[2] - p$d*y[2]
  
  return(derivF$dudt)
}

# compareFandRmodel = function(p=parameters(),N=p$N0,DOC=p$DOC0,B=p$B0) {
#   y = c(N,DOC,B)
#   dudtR = derivative(0,y,p)
#   
#   # Load library
#   loadNUMmodel()
#   # Set parameters
#   dummy = .Fortran("f_setupgeneralistsonly", as.integer(p$n))
#   
#   dudtF = derivativeF(0,y,p)
#   
#   plot(p$m, dudtR[3:(p$n+2)], type="l", log="x", col="red")
#   points(p$m, dudtF[3:(p$n+2)], col="blue", lty=dashed)
#   print(dudtR[1:2])
#   print(dudtF[1:2])
#   
#   return( list(calcRates(p$L,N,DOC,B,p), getFrates(p)) )
# }

simulateChemostatEuler = function(p=parametersChemostat(), bLosses=TRUE) {
  # Load library
  loadNUMmodel()
  # Set parameters
  dummy = .Fortran("f_setupgeneralistsonly", as.integer(p$n))
  
  u0 = c(0.1*p$N0, p$DOC0, p$B0)
  uDeep = c(p$N0, 0)
  out = .Fortran("f_simulatechemostateuler",
                 u = as.numeric(u0),
                 L=as.numeric(p$L), 
                 T=as.numeric(p$T),
                 nNutrients = as.integer(length(uDeep)),
                 uDeep = as.numeric(uDeep),
                 diff = as.numeric(p$d),
                 tend = as.numeric(p$tEnd),
                 dt = as.numeric(0.01),
                 bLosses = as.logical(bLosses))
  
  out = out$u
  ixB = 3:(p$n+2)
  
  B = out[ixB]
  
  result = list(
    p = p,
    t = 0,
    y = out,
    
    N = mean(out[1]),
    DOC = mean(out[2]),
    B = B,
    
    Bmin = B,
    Bmax = B)
  
  rates = getFrates(p)
  
  result$rates = rates
  return(result)
}

simulateChemostat = function(p=parametersChemostat(), 
                             useC=FALSE, useF=TRUE,
                             sSetup="GeneralistsSimpleOnly") {
  
  # Get the version of sundialr:
  #pkg = installed.packages(fields = "Built")
  
  # if (useC) {
  #   # Load library
  #   dyn.load("../Cpp/model.so")
  #   # Set parameters
  #   dummy = .C("setParameters", as.integer(p$n), 
  #              p$m, p$rhoCN, p$epsilonL, p$epsilonF,
  #              p$aNm*p$m, p$aLm*p$m, p$aFm*p$m, p$Jmax, p$jFmaxm*p$m,
  #              p$Jresp, p$Jloss_passive_m,
  #              p$theta,
  #              p$mort, p$mort2, p$mortHTL*p$mortHTLm, p$remin,
  #              p$remin2, p$cLeakage);
  #   
  #   dudt = assign("dudt", rep(0,p$n+2), envir = .GlobalEnv) # Need a static global for speed
  #   if (pkg[pkg[,1]=="sundialr"][3]=="0.1.2")
  #     out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
  #                 IC = c(0.1*p$N0, p$DOC0, p$B0),
  #                 input_function = function(t,y, dummy) derivativeC(t,y,p),
  #                 reltolerance = 1e-6,
  #                 abstolerance = 1e-10+1e-6*c(0.1*p$N0, p$DOC0, p$B0))
  #   else
  #     out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
  #                 IC = c(0.1*p$N0, p$DOC0, p$B0),
  #                 input_function = function(t,y, dummy) derivativeC(t,y,p),
  #                 reltolerance = 1e-6,
  #                 Parameters = 0,
  #                 abstolerance = 1e-10+1e-6*c(0.1*p$N0, p$DOC0, p$B0))
  # }
  # else if (useF) {
    # Load library
    loadNUMmodel()
    # Choose setup:
    if (sSetup=="GeneralistsSimpleOnly") { 
      dummy = .Fortran("f_setupgeneralistssimpleonly", as.integer(p$n))
      dummy = .Fortran("f_sethtl", as.double(p$mHTL), as.double(p$mortHTL), 
                       as.logical(FALSE), as.logical(FALSE))
    } 
    else if (sSetup=="GeneralistsOnly") { 
      dummy = .Fortran("f_setupgeneralistsonly", as.integer(p$n))
      dummy = .Fortran("f_sethtl", as.double(p$mHTL), as.double(p$mortHTL), 
                       as.logical(FALSE), as.logical(FALSE))
    } 
    else if (sSetup=="Generic") {
      dummy = .Fortran("f_setupgeneric", as.integer(4), c(1,10,100,1000) )
      p$n = p$n + 4*10
      p$B0 = c(p$B0, rep(1,4*10))
      #dummy = .Fortran("f_sethtl", as.double(p$mHTL), as.double(p$mortHTL), 
      #                 as.logical(FALSE), as.logical(FALSE))
    } else
    {
      cat('Setup ',sSetup,' no known\n')
      keyboard()
    }
    
    dudt = assign("dudt", rep(as.double(0),p$n+2), envir = .GlobalEnv) # Need a static global for speed
    
    out = ode(y=c(0.1*p$N0, p$DOC0, p$B0),
              times=seq(0, p$tEnd, length.out = p$tEnd+1),
              func = function(t,y,parms)  list(derivativeF(t,y,parms)),
              parms = p)
    
    # if (pkg[pkg[,1]=="sundialr"][3]=="0.1.2")
    #   out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
    #               IC = c(0.1*p$N0, p$DOC0, p$B0),
    #               input_function = function(t,y, dummy) derivativeF(t,y,p),
    #               reltolerance = 1e-4,
    #               abstolerance = 1e-4+1e-4*c(0.1*p$N0, p$DOC0, p$B0))
    # else {
    #   out = t(as.matrix(c(0,0.1*p$N0, p$DOC0, p$B0)))
    # 
    #  # try(
    # #    out = ode(c(0.1*p$N0, p$DOC0, p$B0),
    # #              seq(0, p$tEnd, length.out = p$tEnd),
    # #              function(t,y, dummy) list(derivativeF(t,y,p)))
    #     out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
    #                 IC = c(0.1*p$N0, p$DOC0, p$B0),
    #                input_function = function(t,y, dummy) derivativeF(t,y,p),
    #                 reltolerance = 1e-4,
    #                 Parameters = 0,
    #                 abstolerance = 1e-8+1e-4*c(0.1*p$N0, p$DOC0, p$B0))
    #    ,TRUE)
    #}
  # }
  # else
  # {
  #   if (pkg[pkg[,1]=="sundialr"][3]=="0.1.2")
  #     out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
  #                 IC = c(0.1*p$N0, p$DOC0, p$B0),
  #                 input_function = function(t,y, dummy) derivative(t,y,p),
  #                 reltolerance = 1e-6,
  #                 abstolerance = 1e-10+1e-6*c(0.1*p$N0, p$DOC0, p$B0))
  #   else
  #     out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
  #                 IC = c(0.1*p$N0, p$DOC0, p$B0),
  #                 input_function = function(t,y, dummy) derivative(t,y,p),
  #                 reltolerance = 1e-6,
  #                 Parameters = 0, 
  #                 abstolerance = 1e-10+1e-6*c(0.1*p$N0, p$DOC0, p$B0))
  # }
  nSave = dim(out)[1]
  # Assemble results:
  if (nSave==1)
    ix = 1
  else
    ix = seq(floor(nSave/2),nSave)
  ixB = 4:(p$n+3)
  
  Bmin = 0*p$m
  Bmax = 0*p$m
  for (i in 1:p$n) {
    Bmin[i] = max(1e-20, min(out[ix,ixB[i]]))
    Bmax[i] = max(1e-20, out[ix,ixB[i]])
  }
  
  if (nSave==1)
    B = out[ix, ixB]
  else
    B = colMeans(out[ix,ixB])
  
  result = list(
    p = p,
    t = out[,1],
    y = out[,2:(p$n+3)],
    
    N = mean(out[ix,2]),
    DOC = mean(out[ix,3]),
    B = B,
    
    Bmin = Bmin,
    Bmax = Bmax)
  
  if (useF) 
    rates = getFrates(p)
  else
    rates = calcRates(p$L, result$N, result$DOC, result$B,p)
  
  result$rates = rates
  return(result)
}

getFrates = function(p) {
  nGrid = p$n
  rates = .Fortran("f_getrates", 
                   jN = vector(length=nGrid, mode="numeric"), # Units of carbon
                   jDOC = vector(length=nGrid, mode="numeric"),
                   jL = vector(length=nGrid, mode="numeric"),
                   jSi = vector(length=nGrid, mode="numeric"),
                   jF = vector(length=nGrid, mode="numeric"),
                   jFreal = vector(length=nGrid, mode="numeric"),
                   f = vector(length=nGrid, mode="numeric"),
                   jTot = vector(length=nGrid, mode="numeric"),
                   jMax = vector(length=nGrid, mode="numeric"),
                   jFmax = vector(length=nGrid, mode="numeric"),
                   jR = vector(length=nGrid, mode="numeric"),
                   jRespTot = vector(length=nGrid, mode="numeric"),
                   jLossPassive = vector(length=nGrid, mode="numeric"),
                   jNloss = vector(length=nGrid, mode="numeric"),
                   jLreal = vector(length=nGrid, mode="numeric"),
                   jPOM = vector(length=nGrid, mode="numeric"),
                   mortpred = vector(length=nGrid, mode="numeric"),
                   mortHTL = vector(length=nGrid, mode="numeric"),
                   mort2 = vector(length=nGrid, mode="numeric"),
                   mort = vector(length=nGrid, mode="numeric")
  )
  
  rates$jLossPassive = p$Jloss_passive_m/p$m # NOTE HARDCODED
  rates$jCloss_photouptake = (1-p$epsilonL)/p$epsilonL*rates$jLreal
  rates$jCloss_feeding = (1-p$epsilonF)/p$epsilonF*rates$jFreal
  rates$jNtot = rates$jN + rates$jFreal - rates$jLossPassive
  
  rates$jDOCprodPassive = rates$jLossPassive
  rates$jDOCprodPhoto = rates$jCloss_photouptake
  rates$jDOCprodFeeding = p$reminF * rates$jCloss_feeding 
  rates$jDOCprodVirulysis = p$remin2 * rates$mort2
  
  return(rates)
}

calcFunctionsChemostat = function(p,r,L,N,B) {

  m = p$m
  conversion  = 365*p$M*1e-6*1000 # Convert to gC/yr/m2
  #
  # New production calculated as the flux of nitrogen into the system. Units of carbon:
  #
  prodNew = conversion * p$d*(p$N0-N) * p$rhoCN  
  #
  # Primary production (carbon fixed)
  #
  prodCgross = conversion * sum(r$jLreal*B)/p$epsilonL
  prodCnet = conversion * sum( pmax(0, r$jLreal-p$jR)*B )
  #prodCnet = conversion * sum( r$JLreal*(1 - r$JR/(r$JCtot+r$JR))*B/m )
  #
  # Loss to HTL:
  #
  prodHTL = conversion * sum( p$mortHTL*(m>=p$mHTL)*B )
  #
  # Loss to depth:
  #
  prodSeq = 0#conversion * r$POMgeneration * (1-epsilonPOM)
  #
  # Respiration
  #
  resp = conversion * sum( p$jR*B )
  #
  # Bacterial production:
  #
  jBact_m = pmin(pmax(0,r$jDOC-p$jR-r$jLossPassive), r$jTot)
  prodBact_m = conversion *  jBact_m * B
  prodBact = sum(prodBact_m)
  #prodBact = conversion * sum( r$JDOC*(1-r$JR/(r$JCtot+r$JR))*B/m )
  #
  # Efficiencies:
  #
  effHTL = prodHTL/prodNew # CHECK: CORRECT uNitS?
  effBact = prodBact / prodCnet
  #if (effBact>1)
  #  browser()
  
  #
  # Losses
  #
  lossPassive = conversion * sum( r$jLossPassive*B ) 
  lossPhotouptake = conversion * sum( r$jCloss_photouptake*B )
  lossFeeding = conversion * sum( r$jCloss_feeding*B )
  lossFeedingHTL = sum( (1-p$epsilonF)*prodHTL )
  lossTotalC = lossPassive+lossPhotouptake+lossFeeding+lossFeedingHTL
  lossTotalN = (lossPassive+lossFeeding+lossFeedingHTL)/p$rhoCN
  #
  # Biomasses:
  #
  conversion = p$M*1000*1e-6  # Convert to gC/m2
  d = 2*p$r
  Bpico = conversion * sum( B[d < 2] )
  Bnano = conversion * sum( B[d>=2 & d <20])
  Bmicro = conversion * sum( B[d>=20])
  #
  # Production rate (ratio between net C production and biomass):
  # (in units of day^-1)
  prodRate = prodCnet/(Bpico+Bnano+Bmicro)/365
  
  #
  # Chl-a : Use rough conversion from Edwards et al (2015) that Chl-a propto alpha
  #
  Chl_per_l = sum( r$jLreal/(p$epsilonL*L)*B ) # mugChl/l
  Chl_per_m2 = conversion * Chl_per_l # gC/m2
  
  return(list(
    prodNew = prodNew,
    prodCgross = prodCgross,
    prodCnet = prodCnet,
    prodHTL = prodHTL,
    jBact_m = jBact_m,
    prodBact_m = prodBact_m,
    prodBact = prodBact,
    resp = resp,
    effHTL = effHTL,
    effBact = effBact,
    lossPassive = lossPassive,
    lossPhotouptake = lossPhotouptake,
    lossFeeding = lossFeeding,
    lossFeedingHTL = lossFeedingHTL,
    lossTotalC = lossTotalC,
    lossTotalN = lossTotalN,
    Bpico=Bpico, Bnano=Bnano, Bmicro=Bmicro,
    prodRate=prodRate,
    Chl_per_m2 = Chl_per_m2,
    Chl_per_l = Chl_per_l
  ))
}

plotTimeline = function(sim, time=max(sim$t)) {
  p = sim$p
  t = sim$t
  
  par(cex.axis=cex,
      cex.lab=cex,
      mar=c(4, 5, 6, 2) + 0.1)
  
  y = sim$y
  y[y <= 0] = 1e-30
  
  ylim = c(max(1e-5, min(sim$y)), max(sim$y))
  if (p$latitude==0) {
    xlim = range(t)  
  } else {
    xlim = c(0,365)
    t = 365+t-max(t)
  }
  
  plot(t, y[,1], log="y", type="l", col="blue", 
       ylim=ylim, xlim=xlim, lwd=2,
       xlab="Time (day)", ylab=TeX("Biomass ($\\mu$gC/l)"))
  lines(t, y[,2], col="magenta", lwd=2)
  lines(t, y[,3], col="orange", lwd=2)
  for (i in 1:p$n)
    lines(t, y[,i+2], lwd=i/p$n*3, col="black")
  
  if (p$latitude>0) {
    lines(time*c(1,1), ylim, lty=dotted)
    lines(t, SeasonalLight(p,t),
          col="orange", lwd=2)
  }
}

# plotSeasonal = function(p,time) {
#   defaultplot()
#   defaultpanel(xlim=c(0,365),
#                ylim=c(0,1),
#                xlab="Time (days)",
#                ylab="d/max(d) / L/max(L)")
#   
#   t = seq(0,365,length.out = 100)
#   lines(t, SeasonalExchange(p$latitude, t)/max(SeasonalExchange(p$latitude, t)),col="red", lwd=3)
#   lines(t, SeasonalLight(p, t)/p$L, col="green", lwd=3)
#   vline(x=time)
# }
# 
# plotSeasonalTimeline = function(sim) {
#   require(lattice)
#   require(latticeExtra)
#   
#   p = sim$p
#   ix = sim$t>max(sim$t-365)
#   t = sim$t[ix]
#   
#   B = log10(sim$y[ix,3:(p$n+2)])
#   B[B<(-1)] = -1
#   B[B>3] = 3
#   
#   levelplot(
#     B, 
#     column.values=(p$m), 
#     row.values = t, 
#     aspect="fill",
#     scales = list(y = list(log = 10)),
#     yscale.components = yscale.components.log10ticks,
#     xlab = "Time (days)",
#     ylab = TeX("Carbon mass ($\\mu$gC)"),
#     col.regions = terrain.colors(100),
#     ylim=c(p$m[1]/3, max(p$m)*3)
#   )
# }

#
# Plot functions:
#
plotFunctionsChemostat <- function(sim) {
  #
  # Get functions
  # p,r,L,N,B
  func = calcFunctionsChemostat(sim$p, sim$rates, sim$p$L, sim$N, sim$B)
  # Get the func value from the previous call:
  oldfunc = attr(plotFunctionsChemostat, "oldfunc")
  if (is.null(oldfunc))
    oldfunc = func
  attr(plotFunctionsChemostat, "oldfunc") <<- func
  
  par(mfcol=c(2,1), mar=c(5,12,1,2))
  #
  # Fluxes:
  #
  flux = c(func$prodNew, func$prodCgross, func$prodCnet, func$prodHTL)
  oldflux = c(oldfunc$prodNew, oldfunc$prodCgross, oldfunc$prodCnet, oldfunc$prodHTL)
  
  heights = matrix(c(flux, oldflux), nrow=2, byrow = TRUE)
  barplot(height=heights,
          names.arg = c("New production", "Gross PP", "Net PP", "HTL"),
          xlab = TeX("Production (gC/m$^2$/yr)"),
          beside=TRUE, col=c("black","grey"),
          horiz=TRUE, las=1,
          border=NA)
  legend("topright",
         c("This simulation","Previous simulation"),
         fill=c("black","grey"),
         bty="n")
  #
  # Efficiencies:
  #
  PP = func$prodCnet
  eff = c(func$lossPassive/PP, func$lossPhotouptake/PP,
          func$lossFeeding/PP, func$lossTotalC/PP)
  PP = oldfunc$prodCnet
  oldeff = c(oldfunc$lossPassive/PP, oldfunc$lossPhotouptake/PP,
             oldfunc$lossFeeding/PP, oldfunc$lossTotalC/PP)
  heights = matrix(c(eff, oldeff), nrow=2, byrow=TRUE)
  barplot(height=heights,
          names.arg = c("ePassiveloss", "ePhotoloss", "eFeedingloss", "eTotalloss"),
          xlab = TeX("Fraction of net PP"),
          beside=TRUE, col=c("black","grey"),
          horiz=TRUE, las=1,
          border=NA)
}
#
# Returns the trophic strategy as one of: osmoheterotroph, light or nutrient-limited
# phototroph, mixotroph, or heterotroph.
#
calcStrategy = function(sim, bPlot=FALSE, ylim=c(0.1,1000)) {
  p = sim$p
  r = sim$rates
  
  strategy = rep('Unknown', p$n)
  jNencounter = p$aNm * sim$N
  jLencounter = p$aLm * p$L
  strategy[jNencounter>jLencounter] = "Light limited"
  strategy[jLencounter>=jNencounter] = "Nutrient limited"
  strategy[r$jDOC > r$jLreal] = "Osmoheterotroph"
  strategy[(r$jNloss>1e-5) & (r$jN<r$jF)] = "Heterotroph"
  strategy[((jLencounter/r$jFreal > 1)  & (r$jF>r$jN)) ] = "Mixotroph"  
  
  # Plot shadings:
  if (bPlot) {
    fac = sqrt(p$m[2]/p$m[1])
    
    #strategy = calcStrategy(p,r)
    
    if (ylim[1] <= 0) {
      ys = c(ylim[1]-1, ylim[1]-1, 2*ylim[2], 2*ylim[2])
    } else {
      ys = c(ylim[1]*0.5, ylim[1]*0.5, 2*ylim[2], 2*ylim[2])
    }
    
    for (i in 1:p$n) {
      if (strategy[i]=="Heterotroph")
        col = colHetero
      if (strategy[i]=="Mixotroph")
        col = colMixo
      if (strategy[i]=="Light limited")
        col = colPhoto
      if (strategy[i]=="Nutrient limited")
        col = colN
      if (strategy[i]=="Osmoheterotroph")
        col = colOsmo
      
      polygon(p$m[i]*c(1/fac,fac,fac,1/fac), 
              ys,
              col=col,
              lty=0)
    }
  }
  
  return(strategy)
}

plotSpectrum <- function(sim, t=max(sim$t), bPlot=TRUE, 
                         bXlabel=TRUE, bLegend=TRUE) {
  p = sim$p
  m = p$m
  r = sim$rates
  
  #  par(cex.axis=cex,
  #      cex.lab=cex,
  #      mar=c(4, 5, 6, 2) + 0.1)
  

  ylim = c(0.3,200)

  ixt = which(floor(sim$t)==t+365)[1]
  if (p$latitude==0) {
    B = sim$B / log(p$Delta)
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    B = sim$y[ixt, 3:(p$n+2)] / log(p$Delta)
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(p$L, N, DOC, B, sim$p)
  }
  
  if (bPlot)
    defaultplot(mar=c(2.1,2.3,2.1,0))
  if (bXlabel)
    xlab = "Carbon mass ($\\mu$gC)"
  else
    xlab = ""
  loglogpanel(xlim=p$m, ylim=ylim, xaxis = bXlabel,
              xlab=xlab,
              ylab="Sheldon biomass ($\\mu$gC/l)")
  
  lines(m, B, lwd=8)
  if (p$n<15)
    points(m,B)
  #
  # Theoretical prediction:
  #
  kappa = calcSheldonKappa(sim$p)
  lines(m, rep(kappa, p$n), lty=dotted, col=grey(0.5))
  text(x = 0.7e-8, y = 0.8*kappa,
       labels='Theory',
       cex=0.75*cex, pos=3,col=grey(0.5))
  #     mar=c(4,5,8,2)+0.1)
  #
  # Add gray-scale variation
  #
  if (p$latitude==0)par()
  polygon(c(m, m[seq(p$n,1,by = -1)]), c(sim$Bmin, sim$Bmax[seq(p$n,1,by = -1)])/ log(p$Delta), 
          col=rgb(0.5,0.5,0.5,alpha=0.25), border=NA)
  
  # Determine limiting process:
  strategy = calcStrategy(sim,bPlot=TRUE)
  #
  # Add extra size labels
  #
  d = 10^seq(-6,-1,by=1)
  axis(side=3, line=0,
       at=0.3e6*d^3,
       labels = c("0.01","0.1","1","10","100","1e3"))
  mtext(TeX("Diameter ($\\mu$m)"), side=3, line=1.25, at=1e-3, adj=1,cex=cex)
  #
  # Legend:
  #
  if (bLegend)
  legend(x="topright", bty="n", cex=cex,
         legend=c("Osmoheterotrophs", "Light limited phototrophs","N limited phototrophs","Mixotrophs","Heterotrophs"),
         fill=c(colOsmo, colPhoto,colN,colMixo,colHetero,"transparent"),
         border=c("black","black","black","black","black","transparent"),
         lwd = c(0,0,0,0,0,3),
         col=c(NA,NA,NA,NA,NA,1))
  #
  # Summary state variables: 
  #
  y0 = ylim[1]
  text(x=m[1], y=y0+0.3, labels=TeX(sprintf("DIN: %2.2f $\\mu$mol N/l", N/14)) , cex=cex, pos=4, col=grey(0.5))
  text(x=m[1], y=y0+0.1, labels=TeX(sprintf("DOC: %2.2f $\\mu$mol C/l", DOC/12)), cex=cex, pos=4, col=grey(0.5))
  
  func = calcFunctionsChemostat(sim$p, sim$rates, sim$p$L, sim$N, sim$B)
  text(x=1, y0+1., 
       labels=TeX(sprintf("Chl-a: %1.3f ${\\mu}$gC/l", func$Chl_per_l)),
       cex=cex, pos=2, col=grey(0.5))
  text(x=1, y0+0.6, 
       labels=TeX(sprintf("Picoplankton: %1.2f $gC/m$^2$", func$Bpico)),
       cex=cex, pos=2, col=grey(0.5))
  text(x=1, y0+0.3, 
       labels=TeX(sprintf("Nanoplankton: %1.2f $gC/m$^2$", func$Bnano)),
       cex=cex, pos=2, col=grey(0.5))
  text(x=1, y=y0+0.1, 
       labels=TeX(sprintf("Microplankton: %1.2f $gC/m$^2$", func$Bmicro)),
       cex=cex, pos=2, col=grey(0.5))
  
  box()
}

plotRates = function(sim, p=sim$p, 
                     B=sim$B, N=sim$N, DOC=sim$DOC,
                     t=max(sim$t), 
                     bPlot=TRUE, bLosses=TRUE, bXaxis=TRUE) {
  mm = 10^seq(-8,2, length.out = 100)  
  
  L = p$L
  if (p$latitude!=0) {
    ixt = which(floor(sim$t)==t+365)[1]
    L = SeasonalLight(p,t)
    B = sim$y[ixt, 3:(p$n+2)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
  }
  r = sim$rates
  
  if (bPlot)
    defaultplot()
  ylim = c(-1.3*bLosses-0.1,1.6)
  semilogxpanel(xlim=p$m, ylim=ylim,
                xlab="Carbon mass ($\\mu$gC)",
                ylab="Rates (1/day)",
                xaxis=bXaxis)
  calcStrategy(sim,bPlot=TRUE, ylim=ylim)
  #
  # Gains
  #
  lines(p$m, r$jTot, lwd=10, type="l", col="black")# log="x", xlim=range(p$m),
  lines(p$m, r$jMax, lty=3)
  
  #lines(mm, p$AL*mm^(2/3)*input$L/mm, lty=3, lwd=1, col="green")
  #JLreal = r$Jtot - r$JF+p$Jresp-r$JDOC
  #lines(p$m, p$ALm*p$L/p$m, lty=dotted, lwd=1, col="green")
  lines(p$m, r$jL , lty=dotted, lwd=1, col="green")
  lines(p$m, r$jLreal, lwd=4, col="green")
  
  #lines(mm, p$AN*mm^(1/3)*p$rhoCN*sim$N/mm, lwd=1, lty=3, col="blue")
  #lines(p$m, p$Jmax * p$ANm*N / (p$Jmax/p$rhoCN + p$ANm*N)/p$m, lwd=1, lty=dotted, col="blue")
  #lines(p$m, r$JN/p$m*p$rhoCN, lwd=4, col="blue")
  #lines(p$m, p$JN/p$m, lwd=1, lty=dotted, col="blue")
  lines(p$m, r$jN, lwd=4, col="blue")
  
  #lines(mm, p$AN*mm^(1/3)*DOC/mm, lwd=1, lty=3, col="brown")
  lines(p$m, r$jDOC, lwd=4, col="brown")
  
  lines(p$m, r$jFreal,lwd=4,col="red")
  lines(p$m, p$epsilonF * r$jFmax ,col="red", lty=dotted)
  
  legend(x="topright", cex=cex,
         legend=c("Gains:","Light harvesting","Nutrient uptake","DOC uptake","Food consumption","Division rate"),
         col=c("white","green","blue","brown","red","black"),
         lwd=c(0,4,4,4,4,4),
         bty="n")
  #
  # Losses
  #
  if (bLosses) {
    #polygon(c(1e-9,10,10,1e-9), c(-1.5,-1.5,0,0), 
    #        col=rgb(1,0,0,alpha=0.25), border=NA)
    
    JNexude = r$JNloss 
    lines(p$m, -(r$mortpred + p$mortHTL + r$mort2 + p$mort), lwd=10)
    lines(p$m, -r$mortpred, col="red", lwd=4)
    lines(p$m, -r$mortHTL, col="magenta", lwd=4)
    lines(p$m, -r$mort2, col="orange", lwd=4)
    lines(p$m, -r$jR, col="grey", lwd=4)
    lines(p$m, -r$jLossPassive, col="darkgreen", lwd=4)
    
    #BSheldon =exp(mean(log(B)))
    #delta = (p$m[2]-p$m[1]) / sqrt(p$m[2]*p$m[1])
    #mortPredTheoretical = BSheldon * (1-0.6) * p$AF *sqrt(2*pi)*p$sigma / delta
    #lines(range(p$m), -mortPredTheoretical*c(1,1), lty=dotted, col="red")
    
    legend(x="bottomright", cex=cex,
           legend=c("Losses:", "Predation", "Virulysis", 
                    "Higher trophic levels","Respiration","Passive"),
           col=c(NA,"red", "orange", "magenta","grey","darkgreen"),
           lwd=c(0,4,4,4,4,4), bty="n")
    
    lines(p$m, 0*p$m, col="white", lwd=4)
    lines(p$m, 0*p$m, lty=dashed, lwd=2)
  }
  return(r)
}

plotLeaks = function(sim, t=max(sim$t)) {
  p = sim$p
  
  m = p$m
  ixt = which(floor(sim$t)==t+365)[1]
  if (p$latitude==0) {
    L = p$L
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    L = SeasonalLight(p,t)
    B = sim$y[ixt, 3:(p$n+2)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(L,N, DOC, B, p)
  }
  
  defaultplot()
  semilogxpanel(xlim=m, ylim=c(0,0.4),
                xlab="Carbon mass ($\\mu$gC)",
                ylab="Loss rates (1/day)")
  
  lines(m, r$jNloss-r$jCloss_feeding, col="blue", lwd=4)
  lines(m, r$jCloss_feeding, col="red", lwd=4)
  lines(m, r$jCloss_photouptake, col="green", lwd=4)
  lines(m, r$jLossPassive, col="darkgreen", lwd=4)
  
  legend(x="topright", cex=cex,
         legend=c("Leaks:","N exudation", "C exudation",
                  "N+C sloppy feeding","Passive exudation"),
         col=c("white","blue","green","red","darkgreen"),
         lwd=4, bty="n")
}

plotDeltas = function(sim) {
  rates = sim$rates
  deltaL = rates$jLreal/((1-rates$f)*rates$jL)
  deltaF = rates$jFreal/((1-rates$f)*rates$jF)
  
  defaultplot()
  semilogxpanel(xlim=sim$p$m, ylim=c(0,1),
                xlab="Cell mass ($\\mu$gC)",
                ylab="$\\delta$")
  
  lines(sim$p$m, deltaL, col="green", lwd=4)
  lines(sim$p$m, deltaF, col="red", lwd=4)
}

plotComplexRates = function(sim, t=max(sim$t)) {
  p = sim$p
  
  ixt = which(floor(sim$t)==t+365)[1]
  if (p$latitude==0) {
    L = p$L
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    L = SeasonalLight(p,t)
    B = sim$y[ixt, 3:(p$n+2)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(L, N, DOC, B, p)
  }
  
  par(cex.axis=cex,
      cex.lab=cex,
      mar=c(4, 5, 6, 2) + 0.1)
  
  m = p$m
  plot(m, p$rhoCN*r$JN/m, log="x", type="l", col="blue", 
       ylim=c(-1,2.5), xlim=c(1e-8, 1000),
       lwd=3,
       xlab="Carbon mass (mu gC)",
       ylab="Rates (1/day)")
  lines(m, r$JL/m, col="green", lwd=3)
  lines(m, r$JF/m, col="red", lwd=3)
  lines(m, r$Jtot/m, lty=2)
  lines(m, p$Jmax/m*r$f)
  lines(m, r$dBdt/B,lwd=2)
  lines(m, 0*m,lty=3)
  lines(m, -r$mortpred, col="red")
  lines(m, -r$mort, col="red", lty=2)
  lines(m, -r$JR/m, col="magenta")
  lines(m, -p$mort2*B, col="blue")
  lines(m, -p$mortHTL*(m>=p$mHTL), col="orange")
  legend(x="topright",
         legend=c("rhoCN*JN/m","JL/m","JF/m","Jtot/m","Jmax/m","r","0","mort_pred","mort","resp","mort2","mortHTL"),
         col=c("blue", "green", "red", "black","black","black","black","red", "red", "magenta", "blue","orange"),
         lty=c(1,1,1,2,1,1,3,1,2,1,1,1),
         lwd=c(3,3,3,1,1,2,1,1,1,1,1,1))
}
#
# Theoreitical kappa
#
calcSheldonKappa = function(p) {
  return(
    p$d*p$N0*p$rhoCN / 
      #  (p$mortHTL*log(p$beta) + 0.1*log(p$m[length(p$m)]/p$m[1])))
      (p$mortHTL*log(p$beta) + (p$jR+p$d)*log(p$m[length(p$m)]/p$m[1])))
}
#
# Fitted kappa
#
calcSheldonFit = function(sim, bPlot=TRUE, bXaxis=TRUE) {
  m = sim$p$m
  Delta = m[2]/m[1]
  B = sim$B/Delta
  Bmin = 0.0001*max(B)
  B = pmax(Bmin, sim$B)
  
  #ix = B>Bmin
  #m = m[ix]
  #B = B[ix]
  #
  # Fit slope:
  #
  #f = function(kappa,lambda,mMin,mMax) 
  #  pmax(Bmin, kappa*m^-lambda * (1 - (m/mMin)^-2 - (m/mMax)^2))
  #theta0 = c(2,0,-7,-1)
  
  
  f = function(kappa,lambda,mMin,mMax) {
    pmax(Bmin, kappa*m^-lambda * (1 - (m/mMin)^-2 - (m/mMax)^2))
  }
  
  obj = function(theta) sum( (log(B) - log(f(10^theta[1], theta[2], 10^theta[3], 10^theta[4])))^2)
  fit = nlm(obj, c(-1,0,-9,-1), print.level=0 )
  theta = fit$estimate
  mMin = 10^theta[3]
  mMax = 10^theta[4]
  Bmean = exp(mean(log(B[m>mMin & m<mMax & B>Bmin])))
  
  if (bPlot) {
    #defaultplot()
    if (bXaxis)
      sXlab = "Carbon mass (${\\mu}$gC/l)"
    else
      sXlab = ""
    
    loglogpanel(xlim=m, ylim=c(1e-3,400),
                xlab=sXlab,
                ylab="Sheldon biomass ($\\mu$gC/l)", xaxis=bXaxis)
    B[B==Bmin]=1e-10
    lines(sim$p$m, B, lwd=3)
    #lines(m, rep(Bmean,sim$p$n),lwd=1)
    m = 10^seq(-9,1,length.out=500)
    Bfit = f(10^theta[1], theta[2], 10^theta[3], 10^theta[4])
    Bfit[Bfit==Bmin] = 1e-10
    lines(m, Bfit, lty=dotted)
  }
  
  return(list(Bmean=Bmean, mMin=mMin, mMax=mMax, 
              kappa = 10^theta[1], lambda = theta[2],
              kappaTheo = calcSheldonKappa(sim$p)))
}

testSheldon = function(sim) {
  p=sim$p
  
  JN = p$d*p$N0*p$rhoCN
  Bhtl = sum(sim$B*p$mortHTLm)
  Delta = p$m[2]/p$m[1]
  kappa = calcSheldonKappa(p)
  Bpredict = kappa*Delta
  
  cat('Flux in ', JN, 'mugC/L/day\n')
  cat('Bhtl ', Bhtl, 'mugC/L\n')
  cat('average B', mean(sim$B), 'mugC/L\n')
  cat('HTL flux ', Bhtl*p$mortHTL, 'mugC/L/day\n')
  cat('------\n')
  cat('Predicted Bhtl: ', JN/p$mortHTL,'mugC/L\n')
  cat('Predicted B:', Bpredict, 'mugC/L\n')
  #cat('Predicted Bhtl from B:', JN/(p$mortHTL*p$beta)*Delta, 'mugC/L\n')
  cat('------\n')
  cat('Encounter ',p$aF*Bpredict, ' need ', p$mortHTL,'\n')  
}

baserunChemostat = function(p = parametersChemostat(), useC=FALSE, useF=TRUE) {
  sim = simulateChemostat(p, useC, useF)

  defaultplot(c(2,1))
  plotSpectrum(sim, bPlot=FALSE)
  plotRates(sim, bPlot=FALSE)
  
  return(sim)
}


