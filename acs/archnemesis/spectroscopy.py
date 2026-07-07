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