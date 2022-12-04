# by Jinsen Zheng @ FFPRI 
# 202201127
# The following R code describes a new numerical model framework built on the 
# structure of the published one (NtraceNitrite by Rutting and Muller, 2008) 
# for simulating the soil nitrite dynamics. 

# Users can solve the differential equations using the ode function from the 
# deSolve package (by Soetaert et al.).

mod.NO2Trace <- function(time,      # time step; a numeric vector
                         stocks,    # N pools; a named numeric vector
                         parms      # parameter sets; a named numeric vector
                         ) { 
  with(as.list(c(stocks, parms)), {

    # define the fluxes
    
    # 14N fluxes
    fMin_Nrec_14N = kmin_Nrec*(sNrec_14N/(sNrec_14N+sNrec_15N))               # mineralization of recailcitrant organic N (Nrec)
    fiNH4_Nrec_14N = kiNH4_Nrec*(sNH4_14N/(sNH4_14N+sNH4_15N))                # NH4 immobilized to Nrec
    fMin_Nlab_14N = kmin_Nlab*sNlab_14N                                       # mineralization of labile organic N (Nlab)
    fiNH4_Nlab_14N = kiNH4_Nlab*sNH4_14N                                      # NH4 immobilized to Nlab
    fONrec_NO3_14N = kONrec*(sNrec_14N/(sNrec_14N+sNrec_15N))                 # Nrec oxidized to NO3nit
    fiNO3_Nrec_14N = kiNO3_Nrec*(sNO3nit_14N/(sNO3nit_14N+sNO3nit_15N))       # NO3nit immobilized to Nrec
    fONlab_NO3_14N = kONlab*sNlab_14N                                         # Nlab oxidized to NO3nit
    fiNO3_Nlab_14N = kiNO3_Nlab*sNO3nit_14N                                   # NO3nit immobilized to Nlab
    fONH4_14N = kONH4*sNH4_14N                                                # NH4 oxidation
    fDNRA_14N = kDNRA*sNO2den_14N                                             # DNRA sourced from sNO2den (NO3den may be an alternative)
    fiads_NH4_14N = kiads_NH4*sNH4_14N                                        # adsorption of NH4
    fads_NH4_14N = kads_NH4*sNH4ads_14N                                       # release of adsorped NH4
    fONO2nit_14N = kONO2nit*sNO2nit_14N                                       # NO2nit oxidation
    fhetnit_14N = khet*(sNrec_14N/(sNrec_14N+sNrec_15N))                      # Nrec oxidized to NO2het
    fRNO3nit_14N = kRNO3nit*sNO3nit_14N                                       # reduction of NO3nit
    fRNO2nit_14N = kRNO2nit*sNO2nit_14N                                       # reduction of NO2nit
    fRNO2den_14N = kRNO2den*sNO2den_14N                                       # reduction of NO2den
    fRNO2het_14N = kRNO2het*sNO2het_14N                                       # reduction of NO2het
    fcomammox_14N = kcomammox*sNH4_14N                                        # complete ammonia oxidation 
    fRNO3den_14N = kRNO3den*sNO3den_14N                                       # reduction of NO3den
    fRNO2cnd_14N = kRNO2cnd*sNO2cnd_14N                                       # reduction of NO2cnd
    fONO2den_14N = kONO2den*sNO2den_14N                                       # NO2den oxidation

    # 15N fluxes
    fMin_Nrec_15N = kmin_Nrec*(sNrec_15N/(sNrec_14N+sNrec_15N))
    fiNH4_Nrec_15N = kiNH4_Nrec*(sNH4_15N/(sNH4_14N+sNH4_15N))
    fMin_Nlab_15N = kmin_Nlab*sNlab_15N
    fiNH4_Nlab_15N = kiNH4_Nlab*sNH4_15N
    fONrec_NO3_15N = kONrec*(sNrec_15N/(sNrec_14N+sNrec_15N))
    fiNO3_Nrec_15N = kiNO3_Nrec*(sNO3nit_15N/(sNO3nit_14N+sNO3nit_15N))
    fONlab_NO3_15N = kONlab*sNlab_15N
    fiNO3_Nlab_15N = kiNO3_Nlab*sNO3nit_15N
    fONH4_15N = kONH4*sNH4_15N
    fDNRA_15N = kDNRA*sNO2den_15N
    fiads_NH4_15N = kiads_NH4*sNH4_15N
    fads_NH4_15N = kads_NH4*sNH4ads_15N
    fONO2nit_15N = kONO2nit*sNO2nit_15N
    fhetnit_15N = khet*(sNrec_15N/(sNrec_14N+sNrec_15N))
    fRNO3nit_15N = kRNO3nit*sNO3nit_15N
    fRNO2nit_15N = kRNO2nit*sNO2nit_15N
    fRNO2den_15N = kRNO2den*sNO2den_15N
    fRNO2het_15N = kRNO2het*sNO2het_15N
    fcomammox_15N = kcomammox*sNH4_15N
    fRNO3den_15N = kRNO3den*sNO3den_15N
    fRNO2cnd_15N = kRNO2cnd*sNO2cnd_15N
    fONO2den_15N = kONO2den*sNO2den_15N
    
    
    
    # define the derivatives 
    # NH4 pool
    ddt_NH4_14N <- fMin_Nrec_14N + fMin_Nlab_14N - fiNH4_Nrec_14N - 
      fONH4_14N + fads_NH4_14N - fiNH4_Nlab_14N - fiads_NH4_14N + 
      fDNRA_14N - fcomammox_14N
    ddt_NH4_15N <- fMin_Nrec_15N + fMin_Nlab_15N - fiNH4_Nrec_15N - 
      fONH4_15N + fads_NH4_15N - fiNH4_Nlab_15N - fiads_NH4_15N +
      fDNRA_15N - fcomammox_15N
    # NO2nit pool
    ddt_NO2nit_14N <- fONH4_14N - fONO2nit_14N - fRNO2nit_14N
    ddt_NO2nit_15N <- fONH4_15N - fONO2nit_15N - fRNO2nit_15N
    # NO2den pool
    ddt_NO2den_14N <- fRNO3den_14N - fRNO2den_14N - fDNRA_14N - fONO2den_14N 
    ddt_NO2den_15N <- fRNO3den_15N - fRNO2den_15N - fDNRA_15N - fONO2den_15N 
    # NO2het pool
    ddt_NO2het_14N <- fhetnit_14N - fRNO2het_14N
    ddt_NO2het_15N <- fhetnit_15N - fRNO2het_15N
    # NO2cnd pool
    ddt_NO2cnd_14N <- fRNO3nit_14N - fRNO2cnd_14N
    ddt_NO2cnd_15N <- fRNO3nit_15N - fRNO2cnd_15N
    # NO3 pool
    ddt_NO3nit_14N <- fONO2nit_14N - fRNO3nit_14N + fONrec_NO3_14N - fiNO3_Nrec_14N + 
      fcomammox_14N + fONlab_NO3_14N - fiNO3_Nlab_14N + fONO2den_14N  
    ddt_NO3nit_15N <- fONO2nit_15N - fRNO3nit_15N + fONrec_NO3_15N - fiNO3_Nrec_15N +
      fcomammox_15N + fONlab_NO3_15N - fiNO3_Nlab_15N + fONO2den_15N   
    # NO3den pool
    ddt_NO3den_14N <- -fRNO3den_14N
    ddt_NO3den_15N <- -fRNO3den_15N
    # Nlab pool
    ddt_Nlab_14N <- fiNH4_Nlab_14N - fONlab_NO3_14N + -fMin_Nlab_14N + 
      fiNO3_Nlab_14N 
    ddt_Nlab_15N <- fiNH4_Nlab_15N - fONlab_NO3_15N + -fMin_Nlab_15N + 
      fiNO3_Nlab_15N 
    # Nrec pool
    ddt_Nrec_14N <- -fMin_Nrec_14N - fhetnit_14N  + fiNH4_Nrec_14N + fiNO3_Nrec_14N - fONrec_NO3_14N 
    ddt_Nrec_15N <- -fMin_Nrec_15N - fhetnit_15N  + fiNH4_Nrec_15N + fiNO3_Nrec_15N - fONrec_NO3_15N 
    # NH4ads pool
    ddt_NH4ads_14N <- -fads_NH4_14N + fiads_NH4_14N
    ddt_NH4ads_15N <- -fads_NH4_15N + fiads_NH4_15N
    # Ngas pool
    ddt_Ngas_14N <- fRNO2nit_14N + fRNO2den_14N + fRNO2het_14N + fRNO2cnd_14N #
    ddt_Ngas_15N <- fRNO2nit_15N + fRNO2den_15N + fRNO2het_15N + fRNO2cnd_15N #
    
    # calculate other outputs
    sNH4 = sNH4_14N + sNH4_15N
    sNO2nit = sNO2nit_14N+sNO2nit_15N
    sNO2den = sNO2den_14N+sNO2den_15N
    sNO2het = sNO2het_14N+sNO2het_15N
    sNO2cnd = sNO2cnd_14N+sNO2cnd_15N
    sNO2_15N = sNO2nit_15N+sNO2den_15N+sNO2het_15N + sNO2cnd_15N
    sNO2 = sNO2nit + sNO2den + sNO2het + sNO2cnd
    sNO3nit = sNO3nit_14N + sNO3nit_15N
    sNO3den = sNO3den_14N + sNO3den_15N
    sNO3_15N = sNO3nit_15N+sNO3den_15N
    sNO3 = sNO3nit+sNO3den
    
    atomper_NH4 = sNH4_15N/sNH4
    atomper_NO2 = sNO2_15N/sNO2
    atomper_NO3 = sNO3_15N/sNO3
    
    # output the flux for gross rates
    fMin_Nrec = fMin_Nrec_14N + fMin_Nrec_15N
    fiNH4_Nrec = fiNH4_Nrec_14N + fiNH4_Nrec_15N
    fMin_Nlab = fMin_Nlab_14N + fMin_Nlab_15N
    fiNH4_Nlab = fiNH4_Nlab_14N + fiNH4_Nlab_15N
    fONrec_NO3 = fONrec_NO3_14N + fONrec_NO3_15N
    fiNO3_Nrec = fiNO3_Nrec_14N + fiNO3_Nrec_15N
    fONlab_NO3 = fONlab_NO3_14N + fONlab_NO3_15N
    fiNO3_Nlab = fiNO3_Nlab_14N + fiNO3_Nlab_15N
    fONH4 = fONH4_14N + fONH4_15N
    fDNRA = fDNRA_14N + fDNRA_15N
    fNH4ads = fads_NH4_14N + fads_NH4_15N
    fiNH4ads = fiads_NH4_14N + fiads_NH4_15N
    fONO2nit = fONO2nit_14N + fONO2nit_15N
    fhetnit = fhetnit_14N + fhetnit_15N
    fRNO3nit = fRNO3nit_14N + fRNO3nit_15N
    fRNO2nit = fRNO2nit_14N + fRNO2nit_15N
    fRNO2den = fRNO2den_14N + fRNO2den_15N
    fRNO2het = fRNO2het_14N + fRNO2het_15N
    fcomammox = fcomammox_14N + fcomammox_15N
    fRNO3den = fRNO3den_14N + fRNO3den_15N
    fRNO3 = fRNO3nit + fRNO3den
    fRNO2cnd = fRNO2cnd_14N + fRNO2cnd_15N
    fONO2den = fONO2den_14N + fONO2den_15N
    
    # calculate gross rates
    grossMin = fMin_Nlab + fMin_Nrec  # only consdier the fast the slow ammonification from Nrec here    
    grossImb = fiNH4_Nlab + fiNH4_Nrec # only consider the immobilization of NH4 here
    
    grossP_NH4 = fMin_Nlab +fMin_Nrec + fNH4ads + fDNRA   
    grossC_NH4 = fiNH4_Nlab  + fiNH4ads + fONH4 + fiNH4_Nrec + fcomammox  
    
    grossP_NO2 = fONH4 + fRNO3nit + fRNO3den + fhetnit
    grossC_NO2 = fRNO2nit + fONO2nit + fRNO2den + fRNO2het + fDNRA + fRNO2cnd + fONO2den 
    
    grossP_NO3 = fONlab_NO3 + fONrec_NO3 + fONO2nit + fcomammox + fONO2den 
    grossC_NO3 = fRNO3nit + fRNO3den + fiNO3_Nrec + fiNO3_Nlab #+ fDNRA # note that DNRA sourced from NO3den may offer better fit. 
    

    return(list(c(ddt_NH4_14N,
                  ddt_NH4_15N,
                  ddt_NO2nit_14N,
                  ddt_NO2nit_15N,
                  ddt_NO2den_14N,
                  ddt_NO2den_15N,
                  ddt_NO2het_14N,
                  ddt_NO2het_15N,
                  ddt_NO2cnd_14N,
                  ddt_NO2cnd_15N,
                  ddt_NO3nit_14N,
                  ddt_NO3nit_15N,
                  ddt_NO3den_14N,
                  ddt_NO3den_15N,
                  ddt_Nlab_14N,
                  ddt_Nlab_15N,
                  ddt_Nrec_14N,
                  ddt_Nrec_15N,
                  ddt_NH4ads_14N,
                  ddt_NH4ads_15N,
                  ddt_Ngas_14N,
                  ddt_Ngas_15N),
                sNH4 = sNH4,
                sNO2 = sNO2,
                sNO3 = sNO3,
                aNH4 = atomper_NH4,
                aNO2 = atomper_NO2,
                aNO3 = atomper_NO3,
                
                fMin_Nrec = fMin_Nrec,
                fiNH4_Nrec = fiNH4_Nrec,
                fMin_Nlab = fMin_Nlab,
                fiNH4_Nlab = fiNH4_Nlab,
                fONrec_NO3 = fONrec_NO3,
                fiNO3_Nrec = fiNO3_Nrec,
                fONlab_NO3 = fONlab_NO3,
                fiNO3_Nlab = fiNO3_Nlab,
                fONH4 = fONH4,
                fDNRA = fDNRA,
                fNH4ads = fNH4ads,
                fiNH4ads = fiNH4ads,
                fONO2nit = fONO2nit,
                fhetnit = fhetnit,
                fRNO3nit = fRNO3nit,
                fRNO3den = fRNO3den,
                fRNO3 = fRNO3,
                fRNO2nit = fRNO2nit,
                fRNO2den = fRNO2den,
                fRNO2het = fRNO2het,
                fRNO2cnd = fRNO2cnd,
                fcomammox = fcomammox,
                fONO2den = fONO2den,
                grossMin = grossMin,
                grossImb = grossImb,
                grossP_NH4 = grossP_NH4,
                grossC_NH4 = grossC_NH4,
                grossP_NO2 = grossP_NO2,
                grossC_NO2 = grossC_NO2,
                grossP_NO3 = grossP_NO3,
                grossC_NO3 = grossC_NO3
    ))
  })
}
