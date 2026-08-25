#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# archnemesis.measurement - Functions to create the input measurement for an archnemesis simulation
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

def create_measurement_forward_model(
    lat,
    lon,
    tanhe,
    wavemin,
    wavemax,
    resolving_power=30000.,
    nwave_per_fwhm = 5.,
    ):
    """
    FUNCTION NAME : create_measurement_forward_model()

    DESCRIPTION : Function to create the archNEMESIS Measurement class for a forward model case
                    (i.e., not including real ACS data)

    INPUTS : 

        lat :: Latitude of the tangent point
        lon :: Longitude of the tangent point
        tanhe(ngeom) :: Tangent altitude above the surface (km)
        wavemin :: Minimum wavenumber (cm-1)
        wavemax :: Maximum wavenumber (cm-1)
        resolving_power :: Resolving power (default = 30,000)

    OUTPUTS : 
 
        Measurement :: archNEMESIS Measurement class

    CALLING SEQUENCE:

        Spectroscopy = create_spectroscopy_runtime_class(gas_ids,iso_ids,wavemin,delwave,nwave,line_database)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    Measurement = ans.Measurement_0(ISPACE=0)
    
    midwave = (wavemin+wavemax)/2.
    fwhm = midwave / resolving_power
    nconv = int((wavemax-wavemin)/fwhm * nwave_per_fwhm)
    vconv = np.linspace(wavemin,wavemax,nconv)

    

    Measurement.FWHM = fwhm
    Measurement.ISHAPE = 2   #Gaussian ILS
    Measurement.NGEOM = len(tanhe)
    Measurement.LATITUDE = lat
    Measurement.LONGITUDE = lon
    Measurement.NCONV = np.zeros(Measurement.NGEOM,dtype='int32') + nconv

    ngeom = Measurement.NGEOM
    vconvx = np.zeros((nconv,ngeom))
    for i in range(ngeom):
        vconvx[:,i] = vconv

    Measurement.edit_VCONV(vconvx)
    Measurement.edit_MEAS(np.ones((nconv,ngeom)))
    Measurement.edit_ERRMEAS(np.ones((nconv,ngeom)))
    Measurement.NAV = np.ones(Measurement.NGEOM,dtype='int32')
    Measurement.edit_FLAT(np.zeros((Measurement.NGEOM,1))+lat)
    Measurement.edit_FLON(np.zeros((Measurement.NGEOM,1))+lon)
    Measurement.edit_WGEOM(np.zeros((Measurement.NGEOM,1))+1.0)
    Measurement.edit_EMISS_ANG(np.zeros((Measurement.NGEOM,1))-1.0)  #Negative emission angle to indicate limb-viewing observation
    Measurement.edit_SOL_ANG(np.zeros((Measurement.NGEOM,1))+90.) 
    Measurement.edit_AZI_ANG(np.zeros((Measurement.NGEOM,1))+0.)
    TANHEp = np.zeros((Measurement.NGEOM,1))
    TANHEp[:,0] = tanhe
    Measurement.edit_TANHE(TANHEp)

    return Measurement

#####################################################################################################################################################

def create_measurement_mir(
    lat,lon,
    tanhe,
    vconv,meas,errmeas,
    deldg,fwhm,amp1,amp2,
    ):
    """
    FUNCTION NAME : create_measurement_forward_model()

    DESCRIPTION : Function to create the archNEMESIS Measurement class for a forward model case
                    (i.e., not including real ACS data)

    INPUTS : 

        lat :: Latitude of the tangent point
        lon :: Longitude of the tangent point
        tanhe(ngeom) :: Tangent altitude above the surface (km)
        vconv(nconv) :: Convolution wavenumbers (cm-1)
        meas(nconv,ngeom) :: Measured transmission spectra
        errmeas(nconv,ngeom) :: Uncertainties in measured transmission spectra
        deldg :: Offset between the two Gaussians defining the ILS (according to model 230 in NEMESIS)
        fwhm :: FWHM of the Gaussians defining the ILS (according to model 230 in NEMESIS)
        amp1 :: Amplitude of the second Gaussian at the first wavenumber defining the ILS (according to model 230 in NEMESIS)
        amp2 :: Amplitude of the second Gaussian at the last wavenumber defining the ILS (according to model 230 in NEMESIS)

    OUTPUTS : 
 
        Measurement :: archNEMESIS Measurement class

    CALLING SEQUENCE:

        Measurement = create_measurement_mir(
                        lat,lon,
                        tanhe,
                        vconv,meas,errmeas,
                        deldg,fwhm,amp1,amp2,
                        )

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    Measurement = ans.Measurement_0()
    Measurement.ISPACE = 0      #Wavenumbers 
    Measurement.FWHM = -0.1 
    Measurement.NGEOM = tanhe.shape[0]
    Measurement.LATITUDE = lat
    Measurement.LONGITUDE = lon
    Measurement.NCONV = np.zeros(Measurement.NGEOM,dtype='int32') + vconv.shape[0]
    
    vconvx = np.zeros((vconv.shape[0],Measurement.NGEOM))
    vconvx[:,:] = vconv[:,np.newaxis]

    Measurement.edit_VCONV(vconvx)
    Measurement.edit_MEAS(meas)
    Measurement.edit_ERRMEAS(errmeas)
    Measurement.NAV = np.ones(Measurement.NGEOM,dtype='int32')
    Measurement.edit_FLAT(np.zeros((Measurement.NGEOM,1))+lat)
    Measurement.edit_FLON(np.zeros((Measurement.NGEOM,1))+lon)
    Measurement.edit_WGEOM(np.zeros((Measurement.NGEOM,1))+1.0)
    Measurement.edit_EMISS_ANG(np.zeros((Measurement.NGEOM,1))-1.0)  #Negative emission angle to indicate limb-viewing observation
    Measurement.edit_SOL_ANG(np.zeros((Measurement.NGEOM,1))+90.) 
    Measurement.edit_AZI_ANG(np.zeros((Measurement.NGEOM,1))+0.)
    TANHEp = np.zeros((Measurement.NGEOM,1))
    TANHEp[:,0] = tanhe
    Measurement.edit_TANHE(TANHEp)


    #Building the ILS model
    NWindows = 1
    PAR_FIL = np.zeros((7,1))
    PAR_FIL[0,0] = 1.0e-8
    PAR_FIL[1,0] = 1.0e-8
    PAR_FIL[2,0] = 1.0e-8
    PAR_FIL[3,0] = deldg
    PAR_FIL[4,0] = fwhm
    PAR_FIL[5,0] = amp1
    PAR_FIL[6,0] = amp2
    Measurement = ans.Models[230].calculate(Measurement,NWindows,[Measurement.VCONV[0,0]],[Measurement.VCONV[-1,0]],PAR_FIL)


    return Measurement