#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# mcd - Functions to extract vertical profiles from the Mars Climate Database
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


#MARS CLIMATE DATABASE
#################################################################################

# MCD scenarios:
# 1=Climatology Scenario, solarEUVaverageconditions 
# 2=Climatology Scenario, solarEUVminimumconditions 
# 3=Climatology Scenario, solarEUVmaximumconditions 
# 4=dust storm τ=5, solarminimumconditions 
# 5=dust storm τ=5, solaraveragedconditions 
# 6=dust storm τ=5, solarmaximumconditions 
# 7=warm scenario: dusty atmosphere, solarmax 
# 8=cold scenario: low-dust conditions, solarmin. 
# 24=Mars Year 24, with associated solar EUV conditions. 
# 25=Mars Year 25, with associated solar EUV conditions. 
# 26=Mars Year 26, with associated solar EUV conditions. 
# 27=Mars Year 27, with associated solar EUV conditions. 
# 28=Mars Year 28, with associated solar EUV conditions. 
# 29=Mars Year 29, with associated solar EUV conditions. 
# 30=Mars Year 30, with associated solar EUV conditions. 
# 31=Mars Year 31, with associated solar EUV conditions. 
# 32=Mars Year 32, with associated solar EUV conditions. 
# 33=Mars Year 33, with associated solar EUV conditions. 
# 34=Mars Year 34, with associated solar EUV conditions. 
# 35=Mars Year 35, with associated solar EUV conditions.

from acs.mcd.fmcd import mcd
import numpy as np
import matplotlib.pyplot as plt
import sys,os

#Path for Mars Climate Database
mcddir = '/srv/workspace/data/mcd/MCD_6.1/data/'

########################################################################################################################

def get_mcd_profile(h,gasID,lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,zkey=2,hrkey=1):
    '''
    Function to get the vertical profiles of the volume mixing ratios 

    Inputs
    ------

    h(nh) :: Altitude above the Martian areoid (km)
    gasID(ngas) :: Gas ID
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    zkey :: Altitudes defined above Martian areoid (2) or surface (3)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    vmr(nh,ngas) :: Volume mixing ratios of the different gases
    '''

    #MCD inputs
    ##########################################

    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    nh = len(h)
    ngas = len(gasID)

    vmr = np.zeros((nh,ngas))
    press = np.zeros(nh)
    temp = np.zeros(nh)
    rho = np.zeros(nh)
    for i in range(nh):
        (press[i], rho[i], temp[i], zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,h[i]*1000.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100,dtype='int32'))


        for igas in range(ngas):

            if gasID[igas]==1: #H2O
                vmr[i,igas] = extvar[48-1]
            elif gasID[igas]==2: #CO2
                vmr[i,igas] = extvar[64-1]
            elif gasID[igas]==22: #N2
                vmr[i,igas] = extvar[65-1]
            elif gasID[igas]==76: #Ar
                vmr[i,igas] = extvar[66-1]
            elif gasID[igas]==5: #CO
                vmr[i,igas] = extvar[67-1]
            elif gasID[igas]==45: #O
                vmr[i,igas] = extvar[68-1]
            elif gasID[igas]==7: #O2
                vmr[i,igas] = extvar[69-1]
            elif gasID[igas]==3: #O3
                vmr[i,igas] = extvar[70-1]
            elif gasID[igas]==48: #H
                vmr[i,igas] = extvar[71-1]
            elif gasID[igas]==39: #H2
                vmr[i,igas] = extvar[72-1]
            elif gasID[igas]==40: #He
                vmr[i,igas] = extvar[73-1]
            else:
                sys.exit('error in get_mcd_vmr_profile :: gasID is not included in the MCD')

    return press,temp,rho,vmr