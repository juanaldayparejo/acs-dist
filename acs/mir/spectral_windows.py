#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# spectral_windows - Dictionary for each spectral window to be used in the retrievals
#
# Copyright (C) 2026 Juan Alday
#
# ACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

process_orders_info = {

    ###########################################################################################
    # POSITION 6 - DERIVATION OF 18O/16O in CO2 
    ###########################################################################################
    
    "alday_pos06_237_win1": {
        "Position": 6,                  #ACS secondary grating position
        "Diffor": 237,                  #Diffraction order
        "WaveMin": 3974.8,              #Minimum wavenumber of each spectral window
        "WaveMax": 3988,                #Maximum wavenumber of each spectral window
        "MinTrans": 0.1,                #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,                 #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 100.,               #Maximum altitude to include in the retrievals (km)
        "IDact": [2,2,1],               #Radtran ID of the active gases in atmosphere
        "ISOact": [1,3,1],              #Radtran isotope ID of the active gases in atmosphere
        "split_CO_iso": False,          #Flag indicating whether CO must be separated into its 4 main isotopes
        "split_H2O_iso": True,          #Flag indicating whether H2O must be separated into its 4 main isotopes
        "split_CO2_iso": True,          #Flag indicating whether CO2 must be separated into its 4 main isotopes
        "Flag_Press": False,            #Flag to retrieve pressure
        "Htan": 30.0,                   #If Flag_Press=True, this is the altitude at which pressure is retrieved
        "Flag_Hcorr": False,            #Flag to retrieve tangent height correction
        "Flag_Temp": False,             #Flag to retrieve temperature (based on hydrostatic equilibrium)             
        "Flag_Temp_Analytic": False,    #Flag to retrieve temeprature (based on rotational temperature)
        "Flag_Gas": True,               #Flag to retrieve gas vmrs
        "IDret": [2,2,1],               #Radtran ID of the gases that are meant to be retrieved
        "ISOret": [1,3,1],              #Radtran ID of the gases that are meant to be retrieved
        "Flag_Baseline": True,          #Flag to retrieve the baseline
        "Baseline_Degree": 2,           #Degree of the polynomial to fit the baseline
        "DELDG_apr": -0.14,             #A priori separation between two Gaussians
        "FWHM_apr": 0.09,               #A priori FWHM of the Gaussian
        "AMP1_apr": 0.3,                #A priori amplitude of the second Gaussian at first wavenumber
        "AMP2_apr": 0.3,                #A priori amplitude of the second Gaussian at last wavenumber
    },
    
    "alday_pos06_239_win1": {
        "Position": 6,                  #ACS secondary grating position
        "Diffor": 239,                  #Diffraction order
        "WaveMin": 4012.,               #Minimum wavenumber of each spectral window
        "WaveMax": 4025.,               #Maximum wavenumber of each spectral window
        "MinTrans": 0.1,                #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,                 #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 100.,               #Maximum altitude to include in the retrievals (km)
        "IDact": [2,2,1],               #Radtran ID of the active gases in atmosphere
        "ISOact": [1,3,1],              #Radtran isotope ID of the active gases in atmosphere
        "split_CO_iso": False,          #Flag indicating whether CO must be separated into its 4 main isotopes
        "split_H2O_iso": True,          #Flag indicating whether H2O must be separated into its 4 main isotopes
        "split_CO2_iso": True,          #Flag indicating whether CO2 must be separated into its 4 main isotopes
        "Flag_Press": False,            #Flag to retrieve pressure
        "Htan": 30.0,                   #If Flag_Press=True, this is the altitude at which pressure is retrieved
        "Flag_Hcorr": False,            #Flag to retrieve tangent height correction
        "Flag_Temp": False,             #Flag to retrieve temperature (based on hydrostatic equilibrium)             
        "Flag_Temp_Analytic": False,    #Flag to retrieve temeprature (based on rotational temperature)
        "Flag_Gas": True,               #Flag to retrieve gas vmrs
        "IDret": [2,2,1],               #Radtran ID of the gases that are meant to be retrieved
        "ISOret": [1,3,1],              #Radtran ID of the gases that are meant to be retrieved
        "Flag_Baseline": True,          #Flag to retrieve the baseline
        "Baseline_Degree": 2,           #Degree of the polynomial to fit the baseline
        "DELDG_apr": -0.14,             #A priori separation between two Gaussians
        "FWHM_apr": 0.09,               #A priori FWHM of the Gaussian
        "AMP1_apr": 0.3,                #A priori amplitude of the second Gaussian at first wavenumber
        "AMP2_apr": 0.3,                #A priori amplitude of the second Gaussian at last wavenumber
    },




    ###########################################################################################
    # POSITION 11 - DERIVATION OF D/H ratio
    ###########################################################################################

    "alday_pos11_159_win1": {
        "Position": 11,                 #ACS secondary grating position
        "Diffor": 159,                  #Diffraction order
        "WaveMin": 2662.3,              #Minimum wavenumber of each spectral window
        "WaveMax": 2669.5,              #Maximum wavenumber of each spectral window
        "MinTrans": 0.05,               #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,                 #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 100.,               #Maximum altitude to include in the retrievals (km)
        "IDact": [2,1],                 #Radtran ID of the active gases in atmosphere
        "ISOact": [0,4],                #Radtran isotope ID of the active gases in atmosphere
        "split_CO_iso": False,          #Flag indicating whether CO must be separated into its 4 main isotopes
        "split_H2O_iso": True,          #Flag indicating whether H2O must be separated into its 4 main isotopes
        "split_CO2_iso": False,         #Flag indicating whether CO2 must be separated into its 4 main isotopes

        "Flag_Press": False,            #Flag to retrieve pressure
        "Flag_Hcorr": False,            #Flag to retrieve tangent height correction
        "Flag_Temp": False,             #Flag to retrieve temperature (based on hydrostatic equilibrium)             
        "Flag_Temp_Analytic": False,    #Flag to retrieve temeprature (based on rotational temperature)

        "Flag_Gas": True,               #Flag to retrieve gas vmrs
        "IDret": [2,1],                 #Radtran ID of the gases that are meant to be retrieved
        "ISOret": [0,4],                #Radtran ID of the gases that are meant to be retrieved
        "Gas_Err": [0.2,0.5],           #A priori relative error in the gas vmr profiles 
        "Gas_CLEN": [1.5,1.5],          #Correlation length in the gas vmr profiles 

        "Flag_Baseline": True,          #Flag to retrieve the baseline
        "Baseline_Degree": 2,           #Degree of the polynomial to fit the baseline

        "DELDG_apr": -0.09,             #A priori separation between two Gaussians   (optimised)
        "FWHM_apr": 0.08,               #A priori FWHM of the Gaussian (optimised)
        "AMP1_apr": 0.15,               #A priori amplitude of the second Gaussian at first wavenumber (optimised)
        "AMP2_apr": 0.25,               #A priori amplitude of the second Gaussian at last wavenumber (optimised)
    },

    "alday_pos11_162_win1": {
        "Position": 11,                 #ACS secondary grating position
        "Diffor": 162,                  #Diffraction order
        "WaveMin": 2719.8,              #Minimum wavenumber of each spectral window
        "WaveMax": 2723.7,              #Maximum wavenumber of each spectral window
        "MinTrans": 0.05,               #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,                 #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 100.,               #Maximum altitude to include in the retrievals (km)
        "IDact": [2,1],                 #Radtran ID of the active gases in atmosphere
        "ISOact": [0,4],                #Radtran isotope ID of the active gases in atmosphere
        "split_CO_iso": False,          #Flag indicating whether CO must be separated into its 4 main isotopes
        "split_H2O_iso": True,          #Flag indicating whether H2O must be separated into its 4 main isotopes
        "split_CO2_iso": False,         #Flag indicating whether CO2 must be separated into its 4 main isotopes

        "Flag_Press": False,            #Flag to retrieve pressure
        "Flag_Hcorr": False,            #Flag to retrieve tangent height correction
        "Flag_Temp": False,             #Flag to retrieve temperature (based on hydrostatic equilibrium)             
        "Flag_Temp_Analytic": False,    #Flag to retrieve temeprature (based on rotational temperature)

        "Flag_Gas": True,               #Flag to retrieve gas vmrs
        "IDret": [2,1],                 #Radtran ID of the gases that are meant to be retrieved
        "ISOret": [0,4],                #Radtran ID of the gases that are meant to be retrieved
        "Gas_Err": [0.2,0.5],           #A priori relative error in the gas vmr profiles 
        "Gas_CLEN": [1.5,1.5],          #Correlation length in the gas vmr profiles

        "Flag_Baseline": True,          #Flag to retrieve the baseline
        "Baseline_Degree": 2,           #Degree of the polynomial to fit the baseline
        "DELDG_apr": -0.13,             #A priori separation between two Gaussians (optimised)
        "FWHM_apr": 0.11,               #A priori FWHM of the Gaussian (optimised)
        "AMP1_apr": 1.0,                #A priori amplitude of the second Gaussian at first wavenumber (optimised)
        "AMP2_apr": 1.0,                #A priori amplitude of the second Gaussian at last wavenumber (optimised)
    },

    "alday_pos11_175_win1": {
        "Position": 11,                 #ACS secondary grating position
        "Diffor": 175,                  #Diffraction order
        "WaveMin": 2932.0,              #Minimum wavenumber of each spectral window
        "WaveMax": 2939.0,              #Maximum wavenumber of each spectral window
        "MinTrans": 0.05,               #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,                 #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 100.,               #Maximum altitude to include in the retrievals (km)
        "IDact": [2,1,1],               #Radtran ID of the active gases in atmosphere
        "ISOact": [0,1,4],              #Radtran isotope ID of the active gases in atmosphere
        "split_CO_iso": False,          #Flag indicating whether CO must be separated into its 4 main isotopes
        "split_H2O_iso": True,          #Flag indicating whether H2O must be separated into its 4 main isotopes
        "split_CO2_iso": False,         #Flag indicating whether CO2 must be separated into its 4 main isotopes

        "Flag_Press": False,            #Flag to retrieve pressure
        "Flag_Hcorr": False,            #Flag to retrieve tangent height correction
        "Flag_Temp": False,             #Flag to retrieve temperature (based on hydrostatic equilibrium)             
        "Flag_Temp_Analytic": False,    #Flag to retrieve temeprature (based on rotational temperature)       

        "Flag_Gas": True,               #Flag to retrieve gas vmrs
        "IDret": [1],                   #Radtran ID of the gases that are meant to be retrieved
        "ISOret": [1],                  #Radtran ID of the gases that are meant to be retrieved
        "Gas_Err": [0.5],               #A priori relative error in the gas vmr profiles 
        "Gas_CLEN": [1.5],              #Correlation length in the gas vmr profiles

        "Flag_Baseline": True,          #Flag to retrieve the baseline
        "Baseline_Degree": 2,           #Degree of the polynomial to fit the baseline

        "DELDG_apr": 0.125,             #A priori separation between two Gaussians (optimised)
        "FWHM_apr": 0.125,              #A priori FWHM of the Gaussian (optimised)
        "AMP1_apr": 0.05,               #A priori amplitude of the second Gaussian at first wavenumber (optimised)
        "AMP2_apr": 0.20,               #A priori amplitude of the second Gaussian at last wavenumber (optimised)
    },

}





###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

process_atmosphere_info = {
    
    "3": {
        "Diffor": 206,      #Diffraction order from which to estimate the baseline level
        "WaveMin": 3453.0,         #Diffraction order from which to estimate the baseline level
        "WaveMax": 3461.0,         #Diffraction order from which to estimate the baseline level
        "MinTrans": 0.05,           #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,            #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 110.,          #Maximum altitude to include in the retrievals (km)
    },
    
    "13": {
        "Diffor": 206,      #Diffraction order from which to estimate the baseline level
        "WaveMin": 3453.0,         #Diffraction order from which to estimate the baseline level
        "WaveMax": 3461.0,         #Diffraction order from which to estimate the baseline level
        "MinTrans": 0.05,           #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,            #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 180.,          #Maximum altitude to include in the retrievals (km)
    },
    
    "4": {
        "Diffor": 224,      #Diffraction order from which to estimate the baseline level
        "WaveMin": 3760.0,         #Diffraction order from which to estimate the baseline level
        "WaveMax": 3765.0,         #Diffraction order from which to estimate the baseline level
        "MinTrans": 0.2,           #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,            #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 140.,          #Maximum altitude to include in the retrievals (km)
    },
    
    "5": {
        "Diffor": 227,      #Diffraction order from which to estimate the baseline level
        "WaveMin": 3810.0,         #Diffraction order from which to estimate the baseline level
        "WaveMax": 3815.0,         #Diffraction order from which to estimate the baseline level
        "MinTrans": 0.1,           #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,            #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 100.,          #Maximum altitude to include in the retrievals (km)
    },
    
    "11": {
        "Diffor": 175,        #Diffraction order from which to estimate the baseline level
        "WaveMin": 2937.0,    #Diffraction order from which to estimate the baseline level
        "WaveMax": 2940.0,    #Diffraction order from which to estimate the baseline level
        "MinTrans": 0.05,     #Cut-off of the minimum transmission to include in retrievals (defined if MinTanhe is not defined)
        "MinTanhe": 0.,       #Minimum altitude to include in the retrievals (km)
        "MaxTanhe": 110.,     #Maximum altitude to include in the retrievals (km)
    },
    
}
