#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# archnemesis.spectroscopy - Functions to create the input spectroscopy for an archnemesis simulation
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

location_lbl_tables_hitran24 = "/exomars/retrievals/nemesis/spectroscopy/LBLtables/ACSMIR/HITRAN24/"

###############################################################################################

def create_spectroscopy_runtime_class(
    gas_ids,
    iso_ids,
    wavemin,
    delwave,
    nwave,
    line_database,
    wn_calc_window=25.,
    wn_approx_window=75.,
    iproc=0):

    """
    FUNCTION NAME : create_spectroscopy_runtime_class()

    DESCRIPTION : Function to create the archNEMESIS Spectroscopy class 

    INPUTS : 

        gas_ids(ngas) :: List of RADTRAN IDs to include as active gases
        iso_ids(ngas) :: List of RADTRAN isotope IDs to include as active gases
        wavemin :: Minimum wavenumber (cm-1)
        delwave :: Wavenumber step (cm-1)
        nwave :: Number of spectral points
        line_database :: archnemesis spectroscopic database to use
        wn_calc_window :: Wavenumber window for calculation of lineshape (default = 25 cm-1)
        wn_approx_window :: Wavenumber window for approximation of lineshape (default = 75 cm-1)
        iproc :: Lineshape to use (default = 0, Voigt profile)
        

    OUTPUTS : 
 
        Spectroscopy :: archNEMESIS Spectroscopy class

    CALLING SEQUENCE:

        Spectroscopy = create_spectroscopy_runtime_class(gas_ids,iso_ids,wavemin,delwave,nwave,line_database)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Initialising spectroscopy class
    Spectroscopy  = ans.Spectroscopy_0(ILBL=1)
    Spectroscopy.NGAS = 0

    #Calculating spectral points
    wavemax = wavemin + delwave * (nwave - 1)

    waves = np.arange( wavemin-1., wavemax+1., delwave )

    #Defining line data parameters
    line_data_params = ans.MolLineDataParams(
        lineshape=iproc,
        wn_calc_window=wn_calc_window,
        wn_approx_window=wn_approx_window,
        include_pressure_shift=True,
        s_min=1.0e-50,
        s_floor=0.0,
        amb_gas=[ans.enum.AmbientGasEnum.AIR],
    )

    for igas in range(len(gas_ids)):
        #Editing class
        Spectroscopy.add_line_by_line_runtime(
                mol_id=gas_ids[igas],
                iso_id=iso_ids[igas],
                waves=waves,
                fpath_ld=line_database, 
                wave_unit=0,  #wavenumber
                mol_line_data_params=line_data_params,
        )

    return Spectroscopy

###############################################################################################

def create_spectroscopy_lookup_class(
    lbl_tables,
    ):

    """
    FUNCTION NAME : create_spectroscopy_lookup_class()

    DESCRIPTION : Function to create the archNEMESIS Spectroscopy class 

    INPUTS : 

        lbl_tables :: List of strings including the paths to the pre-computed look-up tables
        

    OUTPUTS : 
 
        Spectroscopy :: archNEMESIS Spectroscopy class

    CALLING SEQUENCE:

        Spectroscopy = create_spectroscopy_runtime_class(lbl_tables)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Initialising spectroscopy class
    Spectroscopy  = ans.Spectroscopy_0(ILBL=2)
    Spectroscopy.NGAS = len(lbl_tables)

    Spectroscopy.ONLINE = True
    Spectroscopy.LOCATION = lbl_tables
    Spectroscopy.assess()

    return Spectroscopy

###############################################################################################

def location_lls_acsmir_hitran24(ngas,gasID,isoID,Datadir=location_lbl_tables_hitran24,online_k=True):

    """

    FUNCTION NAME : location_lls_acsmir_hitran24()

    DESCRIPTION : Create array with the location of the line-by-line tables for ACS MIR retrievals

    INPUTS : 

        ngas :: Number of active gases 
        gasID(ngas) :: RADTRAN ID of each gas
        isoID(ngas) :: RADTRAN isotopologue ID (0 for all isotopes)

    OPTIONAL INPUTS: None
            
    OUTPUTS : 
 
        location(ngas) :: Path to the line-by-line tables

    CALLING SEQUENCE:

        location = location_lls_acsmir(runname,ngas,gasID,isoID)

    MODIFICATION HISTORY : Juan Alday (29/08/2019)

    """

    location = ['']*ngas
    for i in range(ngas):
        if gasID[i]==1:
            strgas = 'H2O'
        elif gasID[i]==2:
            strgas = 'CO2'
        elif gasID[i]==3:
            strgas = 'O3' 
        elif gasID[i]==4:
            strgas = 'N2O'
        elif gasID[i]==5:
            strgas = 'CO'
        elif gasID[i]==8:
            strgas = 'NO'
        elif gasID[i]==11:
            strgas = 'NH3'
        elif gasID[i]==13:
            strgas = 'OH'
        elif gasID[i]==14:
            strgas = 'HF'
        elif gasID[i]==15:
            strgas = 'HCl'
        elif gasID[i]==16:
            strgas = 'HBr'
        elif gasID[i]==18:
            strgas = 'ClO'
        elif gasID[i]==22:
            strgas = 'N2'
        elif gasID[i]==44:
            strgas = 'HO2'
        elif gasID[i]==131:
            strgas = ' Cl2'
        elif gasID[i]==132:
            strgas = ' ClO2'
        else:
            sys.exit('error in write_lls_acsmir :: include gas name in this function')

        if isoID[i]==0:
            filename='ACSMIR_WN_'+strgas+'_HITRAN24'
        else:
            filename='ACSMIR_WN_'+strgas+'_iso'+str(isoID[i])+'_HITRAN24'


        if ((gasID[i]==1) & (isoID[i]==1)):
            filename='ACSMIR_WN_H2O_iso1_SR2022'

        if ((gasID[i]==1) & (isoID[i]==2)):
            filename='ACSMIR_WN_H2O_iso2_SR2022'

        if ((gasID[i]==1) & (isoID[i]==3)):
            filename='ACSMIR_WN_H2O_iso3_SR2022'

        if ((gasID[i]==1) & (isoID[i]==4)):
            filename='ACSMIR_WN_H2O_iso4_SR2022'

        if ((gasID[i]==1) & (isoID[i]==5)):
            filename='ACSMIR_WN_H2O_iso5_SR2022'

        #Defining the extension of the tables
        if online_k==True:
            extension = '.h5'
        else:
            extension = '.lta'
        filename = filename+extension

        location[i] = Datadir+filename
        
    return location

###############################################################################################

def create_spectroscopy_lookup_hitran24(
    ngas,gasID,isoID,Datadir=location_lbl_tables_hitran24,online_k=True
    ):

    """
    FUNCTION NAME : create_spectroscopy_lookup_hitran24()

    DESCRIPTION : Function to create the archNEMESIS Spectroscopy class based on the 
                    HITRAN24 format

    INPUTS : 

        ngas :: Number of active gases 
        gasID(ngas) :: RADTRAN ID of each gas
        isoID(ngas) :: RADTRAN isotopologue ID (0 for all isotopes)

    OPTIONAL INPUTS: None
            
    OUTPUTS : 
 
        Spectroscopy :: Spectroscopy class

    CALLING SEQUENCE:

        Spectroscopy = create_spectroscopy_lookup_hitran24(runname,ngas,gasID,isoID)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    lbl_tables = location_lls_acsmir_hitran24(ngas,gasID,isoID,Datadir=Datadir,online_k=online_k)

    Spectroscopy = create_spectroscopy_lookup_class(lbl_tables)

    return Spectroscopy