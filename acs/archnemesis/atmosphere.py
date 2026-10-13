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
import matplotlib.pyplot as plt

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
    press,temp,rho,vmr = acs.mcd.get_mcd_profile(h,gasID,lat,lon,Ls,LST,scenario=1,zkey=2)

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

def update_pt_atmosphere_with_profile(Atmosphere,h_nir,p_nir,t_nir,extend_downward=False,shift_nir=0.0,MakePlot=False):

    """
    FUNCTION NAME : update_pt_atmosphere_with_profile()

    DESCRIPTION : Update an Atmosphere class from MIR with the NIR pressure and temperature profiles.
    
    INPUTS : 

        Atmosphere :: Atmosphere class created for the MIR profiles
        h_nir :: Altitude array for the NIR profiles (km)
        p_nir :: Pressure array for the NIR profiles
        t_nir :: Temperature array for the NIR profiles

    OPTIONAL INPUTS:
    
        extend_downward :: Flag to indicate whether the NIR profiles should be extended downward to the MIR altitudes
        shift_nir :: Shift to apply to the NIR altitudes (km)

    OUTPUTS : 
 
        Atmosphere :: Atmosphere class updated with the NIR pressure and temperature profiles

    CALLING SEQUENCE:

        Atmosphere = update_pt_atmosphere_with_profile(Atmosphere,h_nir,p_nir,t_nir,extend_downward=False)

    MODIFICATION HISTORY : Juan Alday (22/08/2026)

    """

    h_nir += shift_nir

    #Cropping the maximum altitude for NIR since we cannot extrapolate upward there
    itop = np.where( (Atmosphere.H/1.0e3<=h_nir.max()) & (Atmosphere.T>0.0))[0]

    h_mir = Atmosphere.H[itop]/1.0e3
    t_mir = Atmosphere.T[itop]
    p_mir = Atmosphere.P[itop]
    vmr_mir = Atmosphere.VMR[itop,:]

    Atmosphere.NP = len(itop)
    Atmosphere.edit_H(h_mir*1.0e3)
    Atmosphere.edit_T(t_mir)
    Atmosphere.edit_P(p_mir)
    Atmosphere.edit_VMR(vmr_mir)
    
    if Atmosphere.AMFORM==0:
        molwt = Atmosphere.MOLWT[itop]
        Atmosphere.MOLWT = molwt


    #Interpolating NIR profiles into MIR altitude
    if extend_downward==True:

        z_b = h_nir.min() + 10.
        iextend = np.where(Atmosphere.H/1.0e3<z_b)[0]
        ikeep = np.where(Atmosphere.H/1.0e3>=z_b)[0]

        t_mir_new = np.zeros(t_mir.shape)
        t_mir_new[ikeep] = np.interp(h_mir[ikeep],h_nir,t_nir)

        p_mir_new = np.zeros(p_mir.shape)
        p_mir_new[ikeep] = np.exp(np.interp(h_mir[ikeep],h_nir,np.log(p_nir)))

        if len(iextend)>0:
            
            #Extending the NIR temperature profile downward to the MIR altitudes
            t_nir_b = np.interp(z_b,h_nir,t_nir)
            t_mir_b = np.interp(z_b,h_mir,t_mir)

            deltat = t_nir_b - t_mir_b

            clen = 10.0
            t_mir_new[iextend] = t_mir[iextend] + deltat * np.exp( (h_mir[iextend] - z_b) / clen )

            #Extending the MIR pressure profile using hydrostatic equilibrium and the NIR temperature profile
            
            Atmosphere.calc_grav()
            R_const = 8.31446261815324
            scaleh = R_const * t_mir_new / (Atmosphere.MOLWT * Atmosphere.GRAV)   #scale height (m)

            for i in range(len(iextend)):
                ialt = iextend[len(iextend)-1-i]

                sh =  0.5*(scaleh[ialt]+scaleh[ialt+1])
                delh = (h_mir[ialt]-h_mir[ialt+1])*1.0e3
                p_mir_new[ialt] = p_mir_new[ialt+1]*np.exp(-delh/sh)

        Atmosphere.edit_T(t_mir_new)
        Atmosphere.edit_P(p_mir_new)

    else:

        iuse = np.where(Atmosphere.H/1.0e3>=h_nir.min())[0]
        t_mir_new = np.interp(h_mir[iuse],h_nir,t_nir)
        p_mir_new = np.exp(np.interp(h_mir[iuse],h_nir,np.log(p_nir)))

        Atmosphere.NP = len(iuse)
        Atmosphere.edit_H(h_mir[iuse]*1.0e3)
        Atmosphere.edit_T(t_mir_new)
        Atmosphere.edit_P(p_mir_new)
        Atmosphere.edit_VMR(vmr_mir[iuse,:])

        if Atmosphere.AMFORM==0:
            molwt = Atmosphere.MOLWT[iuse]
            Atmosphere.MOLWT = molwt

    Atmosphere.assess()

    if MakePlot:

        fig,ax = plt.subplots(1,2,figsize=(10,5))
        ax[0].plot(t_nir,h_nir,'r-',label='NIR')
        ax[0].plot(t_mir,h_mir,'b-',label='MIR (initial)')
        ax[0].plot(Atmosphere.T,Atmosphere.H/1.0e3,'g-',label='MIR (updated)')
        ax[0].set_xlabel('Temperature (K)')
        ax[0].set_ylabel('Altitude (km)')
        ax[0].grid()
        ax[1].semilogx(p_nir,h_nir,'r-',label='NIR')
        ax[1].semilogx(p_mir,h_mir,'b-',label='MIR (initial)')
        ax[1].semilogx(Atmosphere.P,Atmosphere.H/1.0e3,'g-',label='MIR (updated)')
        ax[1].set_xlabel('Pressure (Pa)')
        ax[1].set_ylabel('Altitude (km)')
        ax[1].grid()

        plt.legend()
        plt.tight_layout()

    return Atmosphere


###############################################################################################

def update_h2o_atmosphere_with_profile(Atmosphere,h_nir,h2o_nir,extend_downward=False,shift_nir=0.0,MakePlot=False):

    """
    FUNCTION NAME : update_h2o_atmosphere_with_profile()

    DESCRIPTION : Update an Atmosphere class from MIR with the NIR pressure and temperature profiles.
    
    INPUTS : 

        Atmosphere :: Atmosphere class created for the MIR profiles
        h_nir :: Altitude array for the NIR profiles (km)
        h2o_nir :: Volume mixing ratio of H2O 

    OPTIONAL INPUTS:
    
        extend_downward :: Flag to indicate whether the NIR profiles should be extended downward to the MIR altitudes
        shift_nir :: Shift to apply to the NIR altitudes (km)

    OUTPUTS : 
 
        Atmosphere :: Atmosphere class updated with the NIR H2O profiles

    CALLING SEQUENCE:

        Atmosphere = update_h2o_atmosphere_with_profile(Atmosphere,h_nir,h2o_nir,extend_downward=False)

    MODIFICATION HISTORY : Juan Alday (22/08/2026)

    """

    h_nir += shift_nir

    #Locating the index of H2O
    ih2o = Atmosphere.locate_gas(1,0)
    if ih2o is None:
        ih2o = Atmosphere.locate_gas(1,1)
        if ih2o is None:
            raise ValueError("H2O is not present in the reference atmosphere, so we cannot update it with NIR")

    #Cropping the maximum altitude for NIR since we cannot extrapolate upward there
    itop = np.where( (Atmosphere.H/1.0e3<=h_nir.max()) & (Atmosphere.T>0.0))[0]

    Atmosphere.NP = len(itop)
    h_mir = Atmosphere.H[itop]/1.0e3
    t_mir = Atmosphere.T[itop]
    p_mir = Atmosphere.P[itop]
    vmr_mir = Atmosphere.VMR[itop,:]
    h2o_mir = Atmosphere.VMR[itop,ih2o]

    Atmosphere.NP = len(itop)
    Atmosphere.edit_H(h_mir*1.0e3)
    Atmosphere.edit_T(t_mir)
    Atmosphere.edit_P(p_mir)
    Atmosphere.edit_VMR(vmr_mir)

    if Atmosphere.AMFORM==0:
        molwt = Atmosphere.MOLWT[itop]
        Atmosphere.MOLWT = molwt

    #Interpolating NIR profiles into MIR altitude
    if extend_downward==True:

        z_b = h_nir.min() + 0.
        iextend = np.where(Atmosphere.H/1.0e3<z_b)[0]
        ikeep = np.where(Atmosphere.H/1.0e3>=z_b)[0]

        h2o_mir_new = np.zeros(h2o_mir.shape)
        h2o_mir_new[ikeep] = np.interp(h_mir[ikeep],h_nir,h2o_nir)

        if len(iextend)>0:
            
            #Extending the NIR temperature profile downward to the MIR altitudes
            h2o_nir_b = np.interp(z_b,h_nir,h2o_nir)
            h2o_mir_b = np.interp(z_b,h_mir,h2o_mir)

            deltat = h2o_nir_b - h2o_mir_b

            clen = 10.0
            h2o_mir_new[iextend] = h2o_mir[iextend] + deltat * np.exp( (h_mir[iextend] - z_b) / clen )

        Atmosphere.VMR[:,ih2o] = h2o_mir_new

    else:

        iuse = np.where(Atmosphere.H/1.0e3>=h_nir.min())[0]
        h2o_mir_new = np.interp(h_mir[iuse],h_nir,h2o_nir)

        vmr_mir[iuse,ih2o] = h2o_mir_new

        Atmosphere.NP = len(iuse)
        Atmosphere.edit_H(h_mir[iuse]*1.0e3)
        Atmosphere.edit_T(t_mir[iuse])
        Atmosphere.edit_P(p_mir[iuse])
        Atmosphere.edit_VMR(vmr_mir[iuse,:])

        if Atmosphere.AMFORM==0:
            molwt = Atmosphere.MOLWT[iuse]
            Atmosphere.MOLWT = molwt

    Atmosphere.assess()

    if MakePlot:

        fig,ax = plt.subplots(1,1,figsize=(5,5))
        ax.plot(h2o_nir,h_nir,'r-',label='NIR')
        ax.plot(h2o_mir,h_mir,'b-',label='MIR (initial)')
        ax.plot(Atmosphere.VMR[:,ih2o],Atmosphere.H/1.0e3,'g-',label='MIR (updated)')
        ax.set_xlabel('H$_2$O mixing ratio')
        ax.set_ylabel('Altitude (km)')
        ax.grid()

        plt.legend()
        plt.tight_layout()

    return Atmosphere


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

