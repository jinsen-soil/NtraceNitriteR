# by Jinsen Zheng @ FFPRI 
# 20220212
# The following R code describes a new numerical model framework built on the 
# structure of the published one (NtraceNitrite by Rutting and Muller, 2008) 
# for simulating the soil nitrite dynamics. 

# Users can solve the differential equations using the ode function from the 
# deSolve package (by Soetaert et al.).




# time is the time step (defined as a numeric vector)
# stocks refer to N pools (defined as a named numeric vector, including, e.g., sNrec, sNlab, sNH4 below)
# parms for parameter set (defined as a named numeric vector, including, e.g., kmin_zero, kminfast, kONH4 below) 
NtraceNitriteR <- function(time, stocks, parms) { 
  with(as.list(c(stocks, parms)), {

    # define the fluxes
    # 14N
    fMinzero_14N = kmin_zero*(sNrec_14N/(sNrec_14N+sNrec_15N)) # Mineralization from recalcitrant organic N pool
    fMinfast_14N = sNlab_14N*kminfast                      # Mineralization from labile organic N pool
    fONH4_14N = sNH4_14N*kONH4     # NH4 oxidation
    fONO2_14N = sNO2nit_14N*kONO2  # NO2nit oxidation
    fiNH4_14N = kiNH4_zero*(sNH4_14N/(sNH4_14N+sNH4_15N)) # Immobilization of NH4 to labil organic N pool
    fads_NH4_14N = sNH4ads_14N*kads_NH4                   # Release of adsorbed NH4
    fiNO3_14N = kiNO3_zero*(sNO3nit_14N/(sNO3nit_14N+sNO3nit_15N)) # Immobilization of NO3nit 
    fRNO3nit_14N = sNO3nit_14N*kRNO3nit   # Reduction of NO3nit 
    fRNO3den_14N = sNO3den_14N*kRNO3den   # Reduction of NO3den, which is dedicated to partial denitrification
    fhetnit_14N = khet_zero*(sNrec_14N/(sNrec_14N+sNrec_15N)) # heteotrophic nitrification from sNrec
    fRNO2nit_14N = sNO2nit_14N*kRNO2nit # NO2nit -> NOx & N2
    fRNO2den_14N = sNO2den_14N*kRNO2den # NO2den -> NOx & N2
    fRNO2het_14N = sNO2het_14N*kRNO2het # NO2het -> NOx & N2
    fRNO2cnd_14N = sNO2cnd_14N*kRNO2cnd # NO2cnd -> NOx & N2
    fONO2den_14N = kONO2den*sNO2den_14N # Oxidation of NO2den
    fcomammox_14N = kcomammox*sNH4_14N  # Comammox of NH4
    fiNH4_Nrec_14N = kiNH4_Nrec*(sNH4_14N/(sNH4_14N+sNH4_15N)) # Immobilization of NH4 to Nrec
    fDNRA_14N = kDNRA*sNO2den_14N # DNRA (from NO2den to NH4)

    # 15N
    fMinzero_15N = kmin_zero*(sNrec_15N/(sNrec_14N+sNrec_15N))
    fMinfast_15N = sNlab_15N*kminfast
    fONH4_15N = sNH4_15N*kONH4
    fONO2_15N = sNO2nit_15N*kONO2
    fiNH4_15N = kiNH4_zero*(sNH4_15N/(sNH4_14N+sNH4_15N))
    fads_NH4_15N = sNH4ads_15N*kads_NH4
    fiNO3_15N = kiNO3_zero*(sNO3nit_15N/(sNO3nit_14N+sNO3nit_15N))
    fRNO3nit_15N = sNO3nit_15N*kRNO3nit    
    fRNO3den_15N = sNO3den_15N*kRNO3den
    fhetnit_15N = khet_zero*(sNrec_15N/(sNrec_14N+sNrec_15N))
    fRNO2nit_15N = sNO2nit_15N*kRNO2nit 
    fRNO2den_15N = sNO2den_15N*kRNO2den
    fRNO2het_15N = sNO2het_15N*kRNO2het
    fRNO2cnd_15N = sNO2cnd_15N*kRNO2cnd
    fONO2den_15N = kONO2den*sNO2den_15N
    fcomammox_15N = kcomammox*sNH4_15N
    fiNH4_Nrec_15N = kiNH4_Nrec*(sNH4_15N/(sNH4_14N+sNH4_15N))
    fDNRA_15N = kDNRA*sNO2den_15N
 
       
    # define the derivatives 

    # NH4 pool
    ddt_NH4_14N  <- fMinzero_14N + fMinfast_14N - fONH4_14N  + 
      fads_NH4_14N - fiNH4_14N  - fcomammox_14N - fiNH4_Nrec_14N + fDNRA_14N
    ddt_NH4_15N  <- fMinzero_15N + fMinfast_15N - fONH4_15N  + 
      fads_NH4_15N - fiNH4_15N  - fcomammox_15N - fiNH4_Nrec_15N + fDNRA_15N

    # NO2nit pool
    ddt_NO2nit_14N <- fONH4_14N - fONO2_14N - fRNO2nit_14N
    ddt_NO2nit_15N <- fONH4_15N - fONO2_15N - fRNO2nit_15N 
    
    # NO2den pool
    ddt_NO2den_14N <- fRNO3den_14N - fRNO2den_14N - fONO2den_14N - fDNRA_14N
    ddt_NO2den_15N <- fRNO3den_15N - fRNO2den_15N - fONO2den_15N - fDNRA_15N
    
    # NO2het pool
    ddt_NO2het_14N <- fhetnit_14N - fRNO2het_14N
    ddt_NO2het_15N <- fhetnit_15N - fRNO2het_15N
    
    # NO2cnd pool
    ddt_NO2cnd_14N <- fRNO3nit_14N - fRNO2cnd_14N
    ddt_NO2cnd_15N <- fRNO3nit_15N - fRNO2cnd_15N
    
    # NO3nit pool
    ddt_NO3nit_14N <- fONO2_14N - fiNO3_14N - fRNO3nit_14N + fONO2den_14N + fcomammox_14N 
    ddt_NO3nit_15N <- fONO2_15N - fiNO3_15N - fRNO3nit_15N + fONO2den_15N + fcomammox_15N 
    
    # NO3den pool
    ddt_NO3den_14N <- -fRNO3den_14N
    ddt_NO3den_15N <- -fRNO3den_15N
    
    # Nlab pool
    ddt_Nlab_14N <- -fMinfast_14N + fiNH4_14N 
    ddt_Nlab_15N <- -fMinfast_15N + fiNH4_15N 
    
    # Nrec pool
    ddt_Nrec_14N <- -fMinzero_14N + fiNO3_14N - fhetnit_14N  + fiNH4_Nrec_14N  
    ddt_Nrec_15N <- -fMinzero_15N + fiNO3_15N - fhetnit_15N  + fiNH4_Nrec_15N  
    
    # NH4ads pool
    ddt_NH4ads_14N <- -fads_NH4_14N 
    ddt_NH4ads_15N <- -fads_NH4_15N 
    
    # Ngas pool
    ddt_Ngas_14N <- fRNO2nit_14N + fRNO2den_14N + fRNO2het_14N + fRNO2cnd_14N
    ddt_Ngas_15N <- fRNO2nit_15N + fRNO2den_15N + fRNO2het_15N + fRNO2cnd_15N
    
    
    # calculate other outputs for the objective function
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
    
    
    # output the fluxes for gross rates
    fMinfast = fMinfast_14N + fMinfast_15N
    fMinzero = fMinzero_14N + fMinzero_15N
    fNH4ads = fads_NH4_14N + fads_NH4_15N
    fiNH4 = fiNH4_14N + fiNH4_15N
    fONH4 = fONH4_14N+fONH4_15N
    fhetnit = fhetnit_14N + fhetnit_15N
    fRNO3nit = fRNO3nit_14N + fRNO3nit_15N
    fRNO3den = fRNO3den_14N+fRNO3den_15N
    fRNO3 = fRNO3nit+fRNO3den
    fONO2nit = fONO2_14N + fONO2_15N
    fiNO3 = fiNO3_14N + fiNO3_15N
    fONO2den = fONO2den_14N + fONO2den_15N
    fiNH4_Nrec = fiNH4_Nrec_14N + fiNH4_Nrec_15N
    fcomammox = fcomammox_14N + fcomammox_15N
    fDNRA = fDNRA_14N + fDNRA_15N
    fRNO2den = fRNO2den_14N + fRNO2den_15N
    fRNO2het = fRNO2het_14N + fRNO2het_15N
    fRNO2nit = fRNO2nit_14N + fRNO2nit_15N
    fRNO2cnd = fRNO2cnd_14N + fRNO2cnd_15N

    # calculate gross rates
    grossMin = fMinfast + fMinzero + fhetnit # gross mineralization rate
    grossImb = fiNH4  + fiNH4_Nrec + fiNO3   # gross immobilization rate
    
    grossP_NH4 = fMinfast + fMinzero + fNH4ads + fDNRA   # gross NH4 production rate
    grossC_NH4 = fiNH4  + fONH4 + fiNH4_Nrec + fcomammox # gross NH4 consumption rate
    
    grossP_NO2 = fONH4 + fRNO3nit + fRNO3den + fhetnit                                      # gross NO2 production rate
    grossC_NO2 = fONO2nit + fRNO2nit + fRNO2den + fRNO2het  + fRNO2cnd  + fONO2den + fDNRA  # gross NO2 consumption rate
    
    grossP_NO3 = fONO2nit + fcomammox  + fONO2den  # gross NO3 production rate
    grossC_NO3 = fRNO3nit + fRNO3den + fiNO3       # gross NO3 consumption rate
    
    
    # final returns
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
                fMinfast = fMinfast,
                fMinzero = fMinzero,
                fiNH4 = fiNH4,
                fiNH4_Nrec = fiNH4_Nrec,
                fNH4ads = fNH4ads,
                fONH4 = fONH4,
                fhetnit = fhetnit,
                fRNO3nit = fRNO3nit,
                fRNO3den = fRNO3den,
                fRNO3 = fRNO3,
                fONO2nit = fONO2nit,
                fiNO3 = fiNO3,
                fONO2den = fONO2den,
                fcomammox = fcomammox,
                fDNRA = fDNRA,
                fRNO2den = fRNO2den,
                fRNO2het = fRNO2het,
                fRNO2nit = fRNO2nit,
                fRNO2cnd = fRNO2cnd,
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
