#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# files - Functions to extract transmission spectra from the NIR files processed by Anna Fedorova
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

import sys,os
import acs
import numpy as np
import matplotlib.pyplot as plt

#####################################################################################################################################

def prepare_nir_inputs_pt_v1(orbit,nflag,Datadir_NIR='/exomars/data/external/tgo/acs/nir_retrievals/2025_H2O_TP/'):
    """

    FUNCTION NAME : prepare_nir_inputs_pt_v1()

    DESCRIPTION : Read the retrieved profiles from ACS NIR. This is the most updated version of the data up to July 2025.

    INPUTS : 

        orbit :: ACS orbit number
        nflag :: ACS NIR n-flag

    OPTIONAL INPUTS: None
    
        Datadir_NIR :: Path to the dataset published by Fedorova
            
    OUTPUTS : 
 
        h(nalt) :: Altitude above the areoid (km)
        p(nalt) :: Pressure above the areoid (Pa)
        t(nalt) :: Temperature above the areoid (K)
        terr(nalt) :: Temperature above the areoid (K)

    CALLING SEQUENCE:

        exists,h,p,t,terr = prepare_nir_inputs_pt_v1(Datadir_NIR_S,orbit,nflag)

    MODIFICATION HISTORY : Juan Alday (29/05/2024)

    """
    
    Datadir_T = Datadir_NIR+'TP/'
    filename_t = 'orb'+str(orbit).zfill(6)+'_N'+str(nflag)+'_TP_verMay2024.res'
    
    exists = False
    
    if( (os.path.isfile(Datadir_T+filename_t)) ):
        exists = True
    
        #Reading files
        h_nirx,p_nirx,t_nirx,terr_nirx = read_nir_t_v1(Datadir_T+filename_t)

        if len(h_nirx)==0:
            exists=False
            return exists,0.0,0.0,0.0,0.0,0.0,0.0

    else:
        
        h_nirx = 0.
        p_nirx = 0.
        t_nirx = 0.
        terr_nirx = 0.

    return exists,h_nirx,p_nirx,t_nirx,terr_nirx

#####################################################################################################################################

def prepare_nir_inputs_h2o_v1(orbit,nflag,Datadir_NIR='/exomars/data/external/tgo/acs/nir_retrievals/2025_H2O_TP/'):
    """

    FUNCTION NAME : prepare_nir_inputs_h2o_v1()

    DESCRIPTION : Read the retrieved profiles from ACS NIR. This is the most updated version of the data up to July 2025.

    INPUTS : 

        orbit :: ACS orbit number
        nflag :: ACS NIR n-flag

    OPTIONAL INPUTS: None
    
        Datadir_NIR :: Path to the dataset published by Fedorova
            
    OUTPUTS : 
 
        h(nalt) :: Altitude above the areoid (km)
        p(nalt) :: Pressure above the areoid (Pa)
        t(nalt) :: Temperature above the areoid (K)
        terr(nalt) :: Temperature above the areoid (K)
        h2o(nalt) :: Water vapour mixing ratio
        h2oerr(nalt) :: Uncertainty in water vapour mixing ratio

    CALLING SEQUENCE:

        exists,h,p,t,terr,h2o,h2oerr = prepare_nir_inputs_h2o(Datadir_NIR_S,orbit,nflag)

    MODIFICATION HISTORY : Juan Alday (29/05/2024)

    """
    
    Datadir_T = Datadir_NIR+'TP/'
    Datadir_H2O = Datadir_NIR+'H2O/'
    
    filename_h2o = 'orb'+str(orbit).zfill(6)+'_N'+str(nflag)+'_H2O_verMay2024.res'
    filename_t = 'orb'+str(orbit).zfill(6)+'_N'+str(nflag)+'_TP_verMay2024.res'
    
    exists = False
    
    if( (os.path.isfile(Datadir_H2O+filename_h2o)) & (os.path.isfile(Datadir_T+filename_t)) ):
        exists = True
    
        #Reading files
        h_h2o_nir,h2o_nir,h2oerr_nir = read_nir_h2o_v1(Datadir_H2O+filename_h2o)
        h_nir,p_nir,t_nir,terr_nir = read_nir_t_v1(Datadir_T+filename_t)

        if len(h_h2o_nir)==0 or len(h_nir)==0:
            exists=False
            return exists,0.0,0.0,0.0,0.0,0.0,0.0

        #Interpolate T to the water grid
        t_nirx = np.interp(h_h2o_nir,h_nir,t_nir)
        terr_nirx = np.interp(h_h2o_nir,h_nir,terr_nir)
        p_nirx = np.exp(np.interp(h_h2o_nir,h_nir,np.log(p_nir)))
        
    else:
        
        h_h2o_nir = 0.
        press = 0.
        p_nirx = 0.
        t_nirx = 0.
        terr_nirx = 0.
        h2o_nir = 0.
        h2oerr_nir = 0.

    return exists,h_h2o_nir,p_nirx,t_nirx,terr_nirx,h2o_nir,h2oerr_nir

#####################################################################################################################################

def read_nir_h2o_v1(filename):
    """

    FUNCTION NAME : read_nir_h2o_v1()

    DESCRIPTION : Read the retrieved profiles of H2O from ACS NIR. 
                  This is the most updated version of the data up to July 2025.

    INPUTS : 

        filename :: Name of the file

    OPTIONAL INPUTS: None
            
    OUTPUTS : 
 
        h(nalt) :: Altitude above the areoid (km)
        h2o(nalt) :: Water vapour mixing ratio
        h2oerr(nalt) :: Uncertainty in water vapour mixing ratio

    CALLING SEQUENCE:

        h,h2o,h2oerr = read_nir_h2o_v1(filename)

    MODIFICATION HISTORY : Juan Alday (29/05/2024)

    """
    
    from archnemesis.Files import file_lines

    #Reading NIR file
    f = open(filename,'r')
    nalt_h2o = file_lines(filename) - 1
    h_h2o_nir = np.zeros(nalt_h2o)
    h2o_nir = np.zeros(nalt_h2o)
    h2oerr_nir = np.zeros(nalt_h2o)
    header = f.readline()
    for i in range(nalt_h2o):
        s = f.readline().split()
        h_h2o_nir[i] = float(s[0])
        h2o_nir[i] = float(s[5])
        h2oerr_nir[i] = float(s[6])
    f.close()

    #Sorting the arrays
    isort = np.argsort(h_h2o_nir)
    h_h2o_nir = h_h2o_nir[isort]
    h2o_nir = h2o_nir[isort]
    h2oerr_nir = h2oerr_nir[isort]
    
    #Removing negative values
    h2o_nir[h2o_nir<0.0] = 1.0e-10

    #Removing any NaN values
    isel = np.where( (np.isnan(h2o_nir)==False) )[0]
    h_h2o_nir = h_h2o_nir[isel]
    h2o_nir = h2o_nir[isel] * 1.0e-6
    h2oerr_nir = h2oerr_nir[isel] * 1.0e-6

    return h_h2o_nir,h2o_nir,h2oerr_nir

#####################################################################################################################################

def read_nir_t_v1(filename):
    """

    FUNCTION NAME : read_nir_t_s()

    DESCRIPTION : Read the retrieved profiles of temperature published in Fedorova et al. (2022)

    INPUTS : 

        filename :: Name of the file

    OPTIONAL INPUTS: None
            
    OUTPUTS : 
 
        h(nalt) :: Altitude above the areoid (km)
        p(nalt) :: Pressure (Pa)
        t(nalt) :: Temperature (K)
        terr(nalt) :: Temperature uncertainty (K)

    CALLING SEQUENCE:

        h,t,terr = read_nir_t_s(filename)

    MODIFICATION HISTORY : Juan Alday (29/05/2024)

    """
    
    from archnemesis.Files import file_lines

    #Reading NIR file
    f = open(filename,'r')
    nalt = file_lines(filename) - 1
    h_nir = np.zeros(nalt)
    t_nir = np.zeros(nalt)
    terr_nir = np.zeros(nalt)
    p_nir = np.zeros(nalt)
    perr_nir = np.zeros(nalt)
    header = f.readline()
    for i in range(nalt):
        s = f.readline().split()
        h_nir[i] = float(s[0]) 
        t_nir[i] = float(s[5])
        terr_nir[i] = float(s[6])
        p_nir[i] = float(s[7])   #bar
        perr_nir[i] = float(s[8]) #bar
    f.close()

    #Sorting the arrays
    isort = np.argsort(h_nir)
    h_nir = h_nir[isort]
    t_nir = t_nir[isort]
    terr_nir = terr_nir[isort]
    p_nir = p_nir[isort]
    perr_nir = perr_nir[isort]

    #Removing any NaN values
    isel = np.where( (np.isnan(t_nir)==False) )[0]
    h_nir = h_nir[isel]
    t_nir = t_nir[isel]
    terr_nir = terr_nir[isel]
    p_nir = p_nir[isel] * 1.0e5   #Converting bar to Pa

    return h_nir,p_nir,t_nir,terr_nir

