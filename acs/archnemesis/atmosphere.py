#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# archnemesis.atmosphere - Functions to create the input atmosphere for an archnemesis simulation
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

import numpy as np
import archnemesis as ans
import acs

###############################################################################################

def create_mcd_atmosphere_class(lat,lon,Ls,LST,h):

    """
    FUNCTION NAME : create_mcd_atmosphere_class()

    DESCRIPTION : Function to create the archNEMESIS Atmosphere class using profiles from the
                    Mars Climate Database

    INPUTS : 

        lat :: Latitude of the tangent point (degrees)
        lon :: Longitude of the tangent point (degrees)
        Ls :: Solar longitude (degrees)
        LST :: Local time (degrees)
        h(np) :: Altitude array (km)

    OPTIONAL INPUTS:
    
        mcd_year :: If None, it uses the climatology scenario of the MCD. 

    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = create_mcd_atmosphere_class(lat,lon,Ls,LST,tanhe)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    gasID = np.array([1,2,3,5,7,22,39,45,48,76]) #RADTRAN ID for each gas
    press,temp,rho,vmr = acs.get_mcd_profile(h,gasID,lat,lon,Ls,LST,scenario=1,zkey=2)

    #Calculating the molecular weight
    k_B = 1.380649e-23
    amu = 1.660539e-27
    molwt = (rho*k_B*temp/press) / amu / 1.0e3

    Atmosphere_MCD = ans.Atmosphere_0()
    Atmosphere_MCD.NVMR = len(gasID)
    Atmosphere_MCD.ID = np.zeros(Atmosphere_MCD.NVMR,dtype='int32')
    Atmosphere_MCD.ISO = np.zeros(Atmosphere_MCD.NVMR,dtype='int32')
    Atmosphere_MCD.ID[:] = gasID
    Atmosphere_MCD.ISO[:] = np.zeros(Atmosphere_MCD.NVMR,dtype='int32')
    Atmosphere_MCD.NP = len(h)
    Atmosphere_MCD.H = h*1.0e3 #m
    Atmosphere_MCD.edit_P(press)       #Pa
    Atmosphere_MCD.edit_T(temp) 
    Atmosphere_MCD.edit_VMR(vmr) 
    Atmosphere_MCD.MOLWT = molwt  #kg mol-1

    Atmosphere_MCD.LATITUDE = lat
    Atmosphere_MCD.LONGITUDE = lon
    Atmosphere_MCD.AMFORM = 0
    Atmosphere_MCD.IPLANET = 4

    return Atmosphere_MCD

###############################################################################################

def split_H2O_isotopes(Atmosphere,dhratio=5.,o18ratio=1.,o17ratio=1.):

    """
    FUNCTION NAME : split_water_isotopes()

    DESCRIPTION : Function to split H2O into its 4 main isotopes

    INPUTS : 

        Atmosphere :: archNEMESIS atmosphere class

    OPTIONAL INPUTS:
    
        dhratio :: Value of D/H (VSMOW)
        o18ratio :: Value of 18O/16O (VSMOW)
        o17ratio :: Value of 17O/16O (VSMOW)

    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = split_H2O_isotopes(Atmosphere,dhratio=5.,o18ratio=1.,o17ratio=1.)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Splitting H2O into 4 isotopes
    ih2o = np.where( (Atmosphere.ID==1) & (Atmosphere.ISO==0) )[0][0]

    vmr_h2o = Atmosphere.VMR[:,ih2o]

    Atmosphere.remove_gas(1,0)

    #Adding the isotopes of H2O
    Atmosphere.add_gas(1,1,vmr_h2o)
    Atmosphere.add_gas(1,2,vmr_h2o*2005.2e-6*o18ratio)      #VSMOW
    Atmosphere.add_gas(1,3,vmr_h2o*379.9e-6*o17ratio)       #VSMOW
    Atmosphere.add_gas(1,4,vmr_h2o*(155.76e-6*2.)*dhratio)  #VSMOW

    return Atmosphere

###############################################################################################

def split_CO2_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.):

    """
    FUNCTION NAME : split_water_isotopes()

    DESCRIPTION : Function to split H2O into its 4 main isotopes

    INPUTS : 

        Atmosphere :: archNEMESIS atmosphere class

    OPTIONAL INPUTS:
    
        c13ratio :: Value of 13C/12C (VPDB)
        o18ratio :: Value of 18O/16O (VSMOW)
        o17ratio :: Value of 17O/16O (VSMOW)

    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = split_CO2_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Splitting CO2 into 4 isotopes
    ico2 = np.where( (Atmosphere.ID==2) & (Atmosphere.ISO==0) )[0][0]

    vmr_co2 = Atmosphere.VMR[:,ico2]

    Atmosphere.remove_gas(2,0)

    #Adding the isotopes of CO2 isotopes (following HITRAN fractionation)
    Atmosphere.add_gas(2,1,vmr_co2*0.984204)
    Atmosphere.add_gas(2,2,vmr_co2*0.011057*c13ratio)       #VPDB
    Atmosphere.add_gas(2,3,vmr_co2*0.003947*o18ratio)       #VSMOW
    Atmosphere.add_gas(2,4,vmr_co2*7.339890e-4*o18ratio)    #VSMOW

    return Atmosphere

###############################################################################################

def split_CO_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.):

    """
    FUNCTION NAME : split_CO_isotopes()

    DESCRIPTION : Function to split H2O into its 4 main isotopes

    INPUTS : 

        Atmosphere :: archNEMESIS atmosphere class

    OPTIONAL INPUTS:
    
        c13ratio :: Value of 13C/12C (VPDB)
        o18ratio :: Value of 18O/16O (VSMOW)
        o17ratio :: Value of 17O/16O (VSMOW)

    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = split_CO2_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Splitting CO into 4 isotopes
    ico = np.where( (Atmosphere.ID==5) & (Atmosphere.ISO==0) )[0][0]

    vmr_co = Atmosphere.VMR[:,ico]

    Atmosphere.remove_gas(5,0)

    #Adding the isotopes of CO
    Atmosphere.add_gas(5,1,vmr_co)
    Atmosphere.add_gas(5,2,vmr_co*0.0112372*c13ratio)   #VPDB
    Atmosphere.add_gas(5,3,vmr_co*2005.2e-6*o18ratio)   #VSMOW
    Atmosphere.add_gas(5,4,vmr_co*379.9e-6*o17ratio)    #VSMOW

    return Atmosphere