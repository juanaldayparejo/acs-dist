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
        "Flag_OH": False,               #Flag to retrieve OH
        "Flag_HO2": False,              #Flag to retrieve HO2
        "Flag_N2O": False,              #Flag to retrieve N2O
        "Flag_CO2": False,              #Flag to retrieve CO2
        "Flag_CO2_iso1": True,          #Flag to retrieve (12C)(16O)2
        "Flag_CO2_iso2": False,         #Flag to retrieve (13C)(16O)2
        "Flag_CO2_iso3": True,          #Flag to retrieve (18O)(12C)(16O)
        "Flag_CO2_iso4": False,         #Flag to retrieve (17O)(12C)(16O)
        "Flag_H2O_iso1": False,         #Flag to retrieve H2(16O)
        "Flag_H2O_iso2": False,         #Flag to retrieve H2(18O)
        "Flag_H2O_iso3": False,         #Flag to retrieve H2(17O)
        "Flag_H2O_iso4": True,          #Flag to retrieve HD(16O)
        "Flag_CO_iso1": False,          #Flag to retrieve (12C)(16O)
        "Flag_CO_iso2": False,          #Flag to retrieve (13C)(16O)
        "Flag_CO_iso3": False,          #Flag to retrieve (12C)(18O)
        "Flag_CO_iso4": False,          #Flag to retrieve (12C)(17O)
        "Flag_ILS": False,              #Flag to retrieve the ILS
        "Flag_Baseline": True,          #Flag to retrieve the baseline
        "Baseline_Degree": 2,           #Degree of the polynomial to fit the baseline
        "DELDG_apr": -0.14,              #A priori separation between two Gaussians
        "FWHM_apr": 0.09,               #A priori FWHM of the Gaussian
        "AMP1_apr": 0.3,                #A priori amplitude of the second Gaussian at first wavenumber
        "AMP2_apr": 0.3,                #A priori amplitude of the second Gaussian at last wavenumber
    },
    
    "alday_pos06_239_win1": {
        "Position": 6,                  #ACS secondary grating position
        "Diffor": 239,                  #Diffraction order
        "WaveMin": 4012.,              #Minimum wavenumber of each spectral window
        "WaveMax": 4025.,                #Maximum wavenumber of each spectral window
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
        "Flag_OH": False,               #Flag to retrieve OH
        "Flag_HO2": False,              #Flag to retrieve HO2
        "Flag_N2O": False,              #Flag to retrieve N2O
        "Flag_CO2": False,              #Flag to retrieve CO2
        "Flag_CO2_iso1": True,          #Flag to retrieve (12C)(16O)2
        "Flag_CO2_iso2": False,         #Flag to retrieve (13C)(16O)2
        "Flag_CO2_iso3": True,          #Flag to retrieve (18O)(12C)(16O)
        "Flag_CO2_iso4": False,         #Flag to retrieve (17O)(12C)(16O)
        "Flag_H2O_iso1": False,         #Flag to retrieve H2(16O)
        "Flag_H2O_iso2": False,         #Flag to retrieve H2(18O)
        "Flag_H2O_iso3": False,         #Flag to retrieve H2(17O)
        "Flag_H2O_iso4": True,          #Flag to retrieve HD(16O)
        "Flag_CO_iso1": False,          #Flag to retrieve (12C)(16O)
        "Flag_CO_iso2": False,          #Flag to retrieve (13C)(16O)
        "Flag_CO_iso3": False,          #Flag to retrieve (12C)(18O)
        "Flag_CO_iso4": False,          #Flag to retrieve (12C)(17O)
        "Flag_ILS": False,              #Flag to retrieve the ILS
        "Flag_Baseline": True,          #Flag to retrieve the baseline
        "Baseline_Degree": 2,           #Degree of the polynomial to fit the baseline
        "DELDG_apr": -0.14,              #A priori separation between two Gaussians
        "FWHM_apr": 0.09,               #A priori FWHM of the Gaussian
        "AMP1_apr": 0.3,                #A priori amplitude of the second Gaussian at first wavenumber
        "AMP2_apr": 0.3,                #A priori amplitude of the second Gaussian at last wavenumber
    },

}