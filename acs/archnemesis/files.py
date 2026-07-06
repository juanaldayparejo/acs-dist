#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# archnemesis.files - Functions to create the input files for an archnemesis simulation
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

def create_files_forward_model(filename,Atmosphere,Measurement,Spectroscopy):

    """
    FUNCTION NAME : create_files_forward_model()

    DESCRIPTION : Function to create the archNEMESIS input file for a forwad model simulation

    INPUTS : 

        Atmosphere :: archnemesis Atmosphere class
        Measurement :: archnemesis Measurement class
        Spectroscopy :: archnemesis Spectroscopy class

    OUTPUTS : 
 
        archnemesis "filename.h5" file

    CALLING SEQUENCE:

        create_files_forward_model(filename,Atmosphere,Measurement,Spectroscopy)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Creating extra unnecessary classes

    #Scatter class
    Scatter = ans.Scatter_0(NDUST=0,ISPACE=0,ISCAT=0)
    Scatter.NWAVE = 7
    Scatter.WAVE = np.linspace(Measurement.VCONV.min()-2.,Measurement.VCONV.max()+2.,Scatter.NWAVE)
    Scatter.KEXT = np.ones((Scatter.NWAVE,Scatter.NDUST))
    Scatter.SGLALB = np.ones((Scatter.NWAVE,Scatter.NDUST))
    Scatter.KSCA = np.ones((Scatter.NWAVE,Scatter.NDUST))
    Scatter.assess()

    #Surface class
    Surface = ans.Surface_0(GALB=0.0,LOWBC=0,NLOCATIONS=1)
    Surface.TSURF = 220.
    Surface.LATITUDE = Atmosphere.LATITUDE
    Surface.LONGITUDE = Atmosphere.LONGITUDE
    Surface.NEM = 2
    Surface.VEM = np.linspace(Measurement.VCONV.min()-2.,Measurement.VCONV.max()+2.,Surface.NEM)
    Surface.EMISSIVITY = np.ones(Surface.NEM)
    Surface.assess()

    #Stellar class
    Stellar = ans.Stellar_0(SOLEXIST=False,ISPACE=0)
    Stellar.assess()

    #Defining Layer class
    Atmosphere.calc_radius()
    Layer = ans.Layer_0(LAYTYP=5,LAYINT=1,LAYHT=0.0,RADIUS=Atmosphere.RADIUS)
    Layer.NLAY = Atmosphere.NP - 1
    Layer.H_base = Atmosphere.H[0:Atmosphere.NP-1]
    Layer.assess()

    #Retrieval class
    Retrieval = ans.OptimalEstimation_0(IRET=0)
    Retrieval.NITER = -1       #Number of iterations
    Retrieval.PHILIMIT = 0.1   #Convergence criterion
    Retrieval.NCORES = 1       #Number of available cores
    Retrieval.assess_input()

    #Writing a dummy .apr file with a single parameter so that we can initisialise the Variables class
    f = open(filename+'.apr','w')
    f.write('#dummy file \n')
    f.write('1 \n')
    f.write(str(0)+' '+str(0)+' 2 \n')
    f.write('1.0 0.1 \n')
    f.close()

    #Writing input archnemesis file
    Atmosphere.write_hdf5(filename)
    Measurement.write_hdf5(filename)
    Spectroscopy.write_hdf5(filename)
    Layer.write_hdf5(filename)
    Surface.write_hdf5(filename)
    Scatter.write_hdf5(filename)
    Stellar.write_hdf5(filename)
    Retrieval.write_input_hdf5(filename)


