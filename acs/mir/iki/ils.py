#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# iki.ils - Functions to fit the ILS in an observation processed by IKI
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
import acs
import matplotlib.pyplot as plt
from pathlib import Path
import archnemesis as ans
import scipy
import sys,os
import glob

##############################################################################################
##############################################################################################
# DEFINING FUNCTIONS FOR BASELINE SUBSTRACTION
##############################################################################################
##############################################################################################

#@ray.remote
def fit_ils_ave_rows(ACS_datadir,Observation,Window,Rows,Caldir,mintrans=0.1,hsel=None,ratiosel=0.95):
    """
        FUNCTION NAME : fit_ils_ave_rows()
        
        DESCRIPTION : Read the ACS MIR observation and fit the instrument lineshape for
                      one single row on the detector.
        
        INPUTS :
        
            ACS_datadir :: Path to the directory where the IKI ACS MIR data is stored
            Observation :: Name of the ACS MIR observation
            Window      :: Name of the spectral window to be used for the fit
            Rows        :: List of rows to be averaged and used for the fit
            Caldir      :: Path to the directory where the calibration files will be stored  
            
        OPTIONAL INPUTS:
        
            mintrans :: Minimum transmission level of the measurement to be used for the fit (default = 0.2)
            hsel :: Selected tangent height for the fit (default = None)
            ratiosel :: Approximate depth of the absorption lines to be used for the fit (default = 0.9)
        
        OUTPUTS :
        
            ILS calibration file
        
        CALLING SEQUENCE:
        
            fit_ils_ave_rows(ACS_datadir,Observation,Window,Rows,Caldir,mintrans=0.1,hsel=None,ratiosel=0.95)
        
        MODIFICATION HISTORY : Juan Alday (15/06/2024)
    """
    
    #Finding corresponding file

    curr = os.getcwd()
    os.chdir(ACS_datadir)
    filenames = glob.glob("*"+Observation+"*")
    os.chdir(curr)

    if len(filenames) > 1:
        print(filenames)
        raise ValueError("error :: several files were found corresponding to the required observation")
    if len(filenames) == 0:
        raise ValueError("error :: no files were found corresponding to the required observation")

    filename = filenames[0]

    #Reading some parameters specific for this window
    ########################################################################################
    
    Position = acs.mir.spectral_windows.process_orders_info[Window]["Position"]
    DifforSel = acs.mir.spectral_windows.process_orders_info[Window]["Diffor"]
    lowinX = acs.mir.spectral_windows.process_orders_info[Window]["WaveMin"]
    hiwinX = acs.mir.spectral_windows.process_orders_info[Window]["WaveMax"]
    MinTrans = acs.mir.spectral_windows.process_orders_info[Window]["MinTrans"]
    MinTanhe = acs.mir.spectral_windows.process_orders_info[Window]["MinTanhe"]
    MaxTanhe = acs.mir.spectral_windows.process_orders_info[Window]["MaxTanhe"]

    gasIDact = acs.mir.spectral_windows.process_orders_info[Window]["IDact"]
    isoIDact = acs.mir.spectral_windows.process_orders_info[Window]["ISOact"]
    split_CO_iso = acs.mir.spectral_windows.process_orders_info[Window]["split_CO_iso"]
    split_H2O_iso = acs.mir.spectral_windows.process_orders_info[Window]["split_H2O_iso"]
    split_CO2_iso = acs.mir.spectral_windows.process_orders_info[Window]["split_CO2_iso"]   

    DELDG_apr = acs.mir.spectral_windows.process_orders_info[Window]["DELDG_apr"]
    FWHM_apr = acs.mir.spectral_windows.process_orders_info[Window]["FWHM_apr"]
    AMP1_apr = acs.mir.spectral_windows.process_orders_info[Window]["AMP1_apr"]
    AMP2_apr = acs.mir.spectral_windows.process_orders_info[Window]["AMP2_apr"]


    #Getting the ACS MIR data (whole diffraction order) from the FITS files
    ##########################################################################################

    IRows = np.array(Rows)
    refalt = 180.
    
    latX,lonX,LsX,LoctX,VCONVX,MEASX,ERRMEASX,Tanhe_AreoidX = acs.mir.iki.files.extract_order(ACS_datadir+filename,DifforSel,IRows)

    Tanhe_Areoid_AveRowX = np.mean(Tanhe_AreoidX,axis=1) 
    MEAS_AveRowX = np.mean(MEASX,axis=2)
    ERRMEAS_AveRowX = np.mean(ERRMEASX,axis=2)

    #Creating Atmosphere class from the MCD
    hmin = np.ceil(Tanhe_Areoid_AveRowX.min()) - 1
    hmax = np.floor(Tanhe_Areoid_AveRowX.max()).astype(int) + 1
    h = np.arange(hmin,hmax+1,1)
    Atmosphere = acs.archnemesis.atmosphere.create_mcd_atmosphere_class(latX,lonX,LsX,LoctX,h)


    #Preparing Atmosphere class accordingly
    ##########################################################################################

    #Splitting atmosphere into different isotopes if required
    if split_CO_iso==True:
        
        #Splitting CO into 4 isotopes
        ico = np.where( (Atmosphere.ID==5) & (Atmosphere.ISO==0) )[0][0]

        vmr_co = Atmosphere.VMR[:,ico]

        Atmosphere.remove_gas(5,0)

        #Adding the isotopes of CO
        Atmosphere.add_gas(5,1,vmr_co)
        Atmosphere.add_gas(5,2,vmr_co*0.0112372)   #VPDB
        Atmosphere.add_gas(5,3,vmr_co*2005.2e-6)   #VSMOW
        Atmosphere.add_gas(5,4,vmr_co*379.9e-6)    #VSMOW


    if split_CO2_iso==True:

        #Splitting CO2 into 4 isotopes
        ico2 = np.where( (Atmosphere.ID==2) & (Atmosphere.ISO==0) )[0][0]

        vmr_co2 = Atmosphere.VMR[:,ico2]

        Atmosphere.remove_gas(2,0)

        #Adding the isotopes of CO2 isotopes (following HITRAN fractionation)
        Atmosphere.add_gas(2,1,vmr_co2*0.984204)
        Atmosphere.add_gas(2,2,vmr_co2*0.011057)       #VPDB
        Atmosphere.add_gas(2,3,vmr_co2*0.003947)       #VSMOW
        Atmosphere.add_gas(2,4,vmr_co2*7.339890e-4)    #VSMOW


    if split_H2O_iso==True:

        #Splitting H2O into 4 isotopes
        ih2o = np.where( (Atmosphere.ID==1) & (Atmosphere.ISO==0) )[0][0]

        vmr_h2o = Atmosphere.VMR[:,ih2o]

        Atmosphere.remove_gas(1,0)

        #Adding the isotopes of H2O
        Atmosphere.add_gas(1,1,vmr_h2o)
        Atmosphere.add_gas(1,2,vmr_h2o*2005.2e-6)            #VSMOW
        Atmosphere.add_gas(1,3,vmr_h2o*379.9e-6)             #VSMOW
        Atmosphere.add_gas(1,4,vmr_h2o*(155.76e-6*2.)*5.)    #VSMOW
        Atmosphere.add_gas(1,5,vmr_h2o*(155.76e-6*2.)*5.*2005.2e-6)    #VSMOW

    #Defining the dust abundance to 0
    Atmosphere.NDUST = 1
    dust = np.zeros([Atmosphere.NP,Atmosphere.NDUST])
    Atmosphere.edit_DUST(dust)
        
    #Limiting the altitude to only >0 and <200 km

    iin = np.where( (Atmosphere.H>0.0) & (Atmosphere.P>0.0) & (Atmosphere.H/1.0e3<200.) )[0]

    Atmosphere.NP = len(iin)
    try:
        Atmosphere.edit_H(Atmosphere.H[iin])
    except Exception as e:
        
        diffs = Atmosphere.H[iin][1:] - Atmosphere.H[iin][:-1]
        bad_indices = np.where(diffs <= 0)[0]  # These are the *start* indices of non-increasing pairs
        print(Atmosphere.H[iin])

        assert bad_indices.size == 0, (
            f"Array is not strictly increasing at indices: {bad_indices}, "
            f"values: {[(Atmosphere.H[iin][i], Atmosphere.H[iin][i+1]) for i in bad_indices]}"
        )
        raise ValueError("Error editing H in Atmosphere class ("+Observation+"). Check the input data.")
    
    Atmosphere.edit_VMR(Atmosphere.VMR[iin,:])
    Atmosphere.edit_P(Atmosphere.P[iin])
    Atmosphere.edit_T(Atmosphere.T[iin])
    Atmosphere.edit_DUST(Atmosphere.DUST[iin,:])
    Atmosphere.MOLWT = Atmosphere.MOLWT[iin]

    
    #Preparing Measurement class and selecting tangent height for the fit
    ##########################################################################################
    
    #Calculating the measurement in the spectral window
    iwin = np.where( (VCONVX>=lowinX) & (VCONVX<=hiwinX) )[0]

    NCONV_COMB = len(iwin)
    VCONV_COMB = VCONVX[iwin]
    NGEOMSEL = len(IRows)
    MEAS_COMBx = MEAS_AveRowX[iwin,:]
    ERRMEAS_COMBx = ERRMEAS_AveRowX[iwin,:] * 5.0
    TANHESEL_AREOIDx = Tanhe_Areoid_AveRowX[:]

    LATSELx = latX
    LONSELx = lonX

    #Filtering the resulst with a certain level of transmission
    iGEOM = np.where( (np.max(MEAS_COMBx,axis=0)>=mintrans) & (TANHESEL_AREOIDx>Atmosphere.H.min()/1.0e3) )[0]
    MEAS_COMBx = MEAS_COMBx[:,iGEOM]
    ERRMEAS_COMBx = ERRMEAS_COMBx[:,iGEOM]
    TANHESEL_AREOIDx = TANHESEL_AREOIDx[iGEOM]
    
    #Removing the baseline from the spectra
    baseline = np.zeros(MEAS_COMBx.shape)
    for i in range(len(TANHESEL_AREOIDx)):
        #baseline[:,i] = als_baseline(MEAS_COMBx[:,i], asymmetry_param = 0.999,smoothness_param=1e3, max_iters=50, conv_thresh=1e-5, verbose=False)
        baseline[:,i] = als_baseline(MEAS_COMBx[:,i], asymmetry_param = 0.9,smoothness_param=1e1, max_iters=50, conv_thresh=1e-5, verbose=False)

    meas_norm = MEAS_COMBx / baseline

    #Calculating the tangent height at which to fit the ILS
    maxtrans = np.max(meas_norm,axis=0)
    mintrans = np.min(meas_norm,axis=0)

    ratiotrans = mintrans / maxtrans

    if hsel is not None:
        iHsel = np.argmin( np.abs(TANHESEL_AREOIDx-hsel) )
    else:
        iHsel = np.argmin( np.abs(ratiotrans - ratiosel) )

    print('Selected tangent height = ',TANHESEL_AREOIDx[iHsel],'km')

    #Writing information into class
    Measurement = ans.Measurement_0(ISPACE=0)
    Measurement.FWHM = -0.1
    Measurement.NGEOM = 1
    Measurement.LATITUDE = Atmosphere.LATITUDE
    Measurement.LONGITUDE = Atmosphere.LONGITUDE
    Measurement.NCONV = np.zeros(Measurement.NGEOM,dtype='int32') + NCONV_COMB
    VCONV_COMBx = np.zeros((NCONV_COMB,Measurement.NGEOM))
    VCONV_COMBx[:,:] = VCONV_COMB[:,None]
    Measurement.edit_VCONV(VCONV_COMBx)
    Measurement.edit_MEAS(MEAS_COMBx[:,iHsel,None])
    Measurement.edit_ERRMEAS(ERRMEAS_COMBx[:,iHsel,None])
    Measurement.NAV = np.ones(Measurement.NGEOM,dtype='int32')
    Measurement.edit_FLAT(np.zeros((Measurement.NGEOM,1))+Atmosphere.LATITUDE)
    Measurement.edit_FLON(np.zeros((Measurement.NGEOM,1))+Atmosphere.LONGITUDE)
    Measurement.edit_WGEOM(np.zeros((Measurement.NGEOM,1))+1.0)
    Measurement.edit_EMISS_ANG(np.zeros((Measurement.NGEOM,1))-1.0)  #Negative emission angle to indicate limb-viewing observation
    Measurement.edit_SOL_ANG(np.zeros((Measurement.NGEOM,1))+90.) 
    Measurement.edit_AZI_ANG(np.zeros((Measurement.NGEOM,1))+0.)
    TANHEp = np.zeros((Measurement.NGEOM,1))
    TANHEp[:,0] = TANHESEL_AREOIDx[iHsel]
    Measurement.edit_TANHE(TANHEp)

    #Dfining the lineshape
    NWindows = 1
    PAR_FIL = np.zeros((7,1))
    PAR_FIL[0,0] = 1.0e-8
    PAR_FIL[1,0] = 1.0e-8
    PAR_FIL[2,0] = 1.0e-8
    PAR_FIL[3,0] = DELDG_apr
    PAR_FIL[4,0] = FWHM_apr
    PAR_FIL[5,0] = AMP1_apr
    PAR_FIL[6,0] = AMP2_apr

    lowin = Measurement.VCONV.min()-0.001
    hiwin = Measurement.VCONV.max()+0.001
    Measurement = ans.Models[230].calculate(Measurement,NWindows,[lowin],[hiwin],PAR_FIL)
    
    #Creating Spectroscopy class
    ###############################################################################
    
    #Defining the Spectroscopy
    Spectroscopy = ans.Spectroscopy_0()
    Spectroscopy.NGAS = len(gasIDact)
    Spectroscopy.ILBL = 2

    llslocation = acs.archnemesis.spectroscopy.location_lls_acsmir_hitran24(len(gasIDact),gasIDact,isoIDact)

    Spectroscopy.LOCATION = llslocation
    Spectroscopy.ONLINE = True
    Spectroscopy.read_header()

    vmin, vmax = Measurement.calc_wave_range()

    #Adding some more extra point to the wavelength
    extra_wave = 4.
    vmin -= extra_wave
    vmax += extra_wave

    #Now, reading k-tables or lbl-tables for the spectral range of interest
    Spectroscopy.read_tables(wavemin=vmin,wavemax=vmax)
    
    #Calculating the optical depth of each gas from the reference atmosphere
    ##########################################################################################

    WAVE,TAUGAS_PATH = calc_tau_gas(Atmosphere,Measurement,Spectroscopy)

    TRANS = np.exp(-np.sum(TAUGAS_PATH,axis=1))

    #Performing the first initial wavenumber calibration
    ##########################################################################################
    
    vary_wavecal = True
    extra_wave = 4.0
    if vary_wavecal is True:
    
        #Defining Measurement class for first wavenumber calibration
        step = (FWHM_apr)/2./5.
        VCONVc = np.arange(Measurement.VCONV[:,0].min()-(extra_wave-1.),Measurement.VCONV[:,0].max()+(extra_wave-1.),step)
        Measurementc = ans.Measurement_0()
        Measurementc.NGEOM=1
        Measurementc.NCONV = np.array([len(VCONVc)])
        Measurementc.VCONV = np.zeros((len(VCONVc),1))
        Measurementc.VCONV[:,0] = VCONVc[:]
        
        Measurementc = modelILS(Measurementc,DELDG_apr,FWHM_apr,AMP1_apr,AMP2_apr)
        TRANSCONVX = ans.lblconv_fil(Spectroscopy.NWAVE,Spectroscopy.WAVE,TRANS,Measurementc.NCONV[0],Measurementc.VCONV[:,0],Measurementc.NFIL,Measurementc.VFIL,Measurementc.AFIL)

        #Specifying a priori values for wavenumber calibration
        pixx = np.arange(0,Measurement.NCONV[0],1)
        v0 = Measurement.VCONV[0,0]
        dV = Measurement.VCONV[:,0] - Measurement.VCONV[0,0]
        
        px = np.polyfit(pixx,dV,2)
        C0 = px[2]
        C1 = px[1]
        C2 = px[0]
        dVm = construct_dVm(pixx,C0,C1,C2)
        newconv = construct_vconv(pixx,v0+0.0,0.0,C1,C2)
        
        #Fitting polynomial function to both measured and modelled spectra to normalise them
        ndeg = 2
        px = np.polyfit(pixx,Measurement.MEAS[:,0],ndeg)
        pol = px[2] + px[1]*pixx + px[0]*(pixx)**2.
        
        #meas = Measurement.MEAS[:,0] - pol
        meas = Measurement.MEAS[:,0]
        
        px = np.polyfit(Measurementc.VCONV[:,0],TRANSCONVX,ndeg)
        polm = px[2] + px[1]*Measurementc.VCONV[:,0] + px[0]*(Measurementc.VCONV[:,0])**2.
        
        model = TRANSCONVX - polm
        model = TRANSCONVX
        smodel = scipy.interpolate.interp1d(Measurementc.VCONV[:,0],model)
        
        v0range = 0.5
        v0step = 0.02
        C1range = 10.
        C1step = 0.5
        C2range = 20.
        C2step = 5.
        
        v0x = np.arange(v0-v0range,v0+v0range+v0step,v0step)
        C1x = C1 * (1.0 + np.arange(-C1range,C1range+C1step,C1step)/100.)
        C2x = C2 * (1.0 + np.arange(-C2range,C2range+C2step,C2step)/100.)
        
        #Performing cross correlation for all cases
        corr = np.zeros((len(v0x),len(C1x),len(C2x)))
        for idV in range(len(v0x)):
            for iC1 in range(len(C1x)):
                for iC2 in range(len(C2x)):
                    corr[idV,iC1,iC2] = correlate_spectra(pixx,v0x[idV],C1x[iC1],C2x[iC2],meas,smodel,normalise=True)
        
        #Getting the point that maximises the correlation between measured and modelled spectra
        imax = np.unravel_index(corr.argmax(), corr.shape)

        v0fit = v0x[imax[0]]
        C1fit = C1x[imax[1]]
        C2fit = C2x[imax[2]]
        
        dVfit = construct_dVm(pixx,0.0,C1fit,C2fit) 
        vconvfit = construct_vconv(pixx,v0fit,0.0,C1fit,C2fit)

    else:
        
        #Specifying a priori values for wavenumber calibration
        pixx = np.arange(0,Measurement.NCONV[0],1)
        v0 = Measurement.VCONV[0,0]
        dV = Measurement.VCONV[:,0] - Measurement.VCONV[0,0]
        
        px = np.polyfit(pixx,dV,2)
        v0fit = px[2] + Measurement.VCONV[0,0]
        C1fit = px[1]
        C2fit = px[0]

        dVfit = construct_dVm(pixx,0.0,C1fit,C2fit) 
        vconvfit = construct_vconv(pixx,v0fit,0.0,C1fit,C2fit)


    if np.where(np.isnan(TAUGAS_PATH)==True)[0].size>0:
        print(Observation,Window)
        raise ValueError('There are NaN values in the computed optical depths')

    if np.where(np.isnan(Measurement.MEAS)==True)[0].size>0:
        print(Observation,Window)
        raise ValueError('There are NaN values in the measurement')

    #Building model for the fitting of the ILS
    ##########################################################################################
    
    from lmfit import Parameters, Model, Minimizer

    params = Parameters()

    vary_DELDG = True
    vary_FWHM = True

    # Parameters for ILS model
    params.add('DELDG', value=DELDG_apr,vary=vary_DELDG,min=-0.2,max=0.2)
    params.add('FWHM', value=FWHM_apr,vary=vary_FWHM,min=0.05,max=0.25)
    params.add('AAMP1', value=AMP1_apr,vary=True,min=0.0)
    params.add('AAMP2', value=AMP2_apr,vary=True,min=0.0)

    # Parameters for wavenumber calibration
    C0ini = v0fit
    C1ini = C1fit
    C2ini = C2fit

    params.add('C0',vary=True,min=C0ini-0.5,max=C0ini+0.5,value=C0ini)
    params.add('C1',vary=True,min=C1ini*0.30,max=C1ini*1.6,value=C1ini)
    params.add('C2',vary=True,min=C2ini*0.30,max=C2ini*1.6,value=C2ini)

    # Parameters for the gas abundance fitting in each row
    n_geoms = Measurement.NGEOM
    ngases = Spectroscopy.NGAS

    for i in range(n_geoms):
        for g in range(ngases):
            params.add(f'SC_GAS{g+1}_{i+1}', value=1.0, min=0.0, vary=True)

    Measurement.NCONV[:] = len(vconvfit)
    Measurement.VCONV[:,:] = vconvfit[:,None]


    #Fitting the ILS model
    ##########################################################################################

    #Define the measurement vector
    ymeas = np.zeros(np.sum(Measurement.NCONV))
    ix = 0
    for i in range(Measurement.NGEOM):
        ymeas[ix:ix+Measurement.NCONV[i]] = Measurement.MEAS[:,i]
        ix = ix + Measurement.NCONV[i]

    # Define a residual function to compare model output to measured data
    def residual(params, Measurement, WAVE, TAUGAS_PATH, ymeas):
        ymodel = fm_calibrate(Measurement, WAVE, TAUGAS_PATH, params)
        return ymeas - ymodel

    # Setup Minimizer (assume `params` is already created and populated)
    minner = Minimizer(residual,
                    params,
                    fcn_args=(Measurement, Spectroscopy.WAVE, TAUGAS_PATH, ymeas))

    # Run minimisation
    result = minner.minimize(method='leastsq', max_nfev=(Spectroscopy.NGAS*Measurement.NGEOM + 9)*10)

    # Best-fit model
    yfit = fm_calibrate(Measurement, Spectroscopy.WAVE, TAUGAS_PATH, result.params)

    #for name, param in result.params.items():
    #    print(f"{name}: {param.value} ± {param.stderr}")

    DELDG_fit = result.params['DELDG'].value
    FWHM_fit = result.params['FWHM'].value
    AAMP1_fit = result.params['AAMP1'].value
    AAMP2_fit = result.params['AAMP2'].value
    C0_fit = result.params['C0'].value
    C1_fit = result.params['C1'].value
    C2_fit = result.params['C2'].value
    
    #Calculating the error in the retrieved parameters
    ##########################################################################################
    
    DELDG_fit_err = result.params['DELDG'].stderr
    FWHM_fit_err = result.params['FWHM'].stderr
    AAMP1_fit_err = result.params['AAMP1'].stderr
    AAMP2_fit_err = result.params['AAMP2'].stderr

    #Making a summary figure
    ##########################################################################################

    fig = plt.figure(figsize=(12,8))
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
    ax2 = plt.subplot2grid((3, 3), (2, 0))
    ax3 = plt.subplot2grid((3, 3), (2, 1))
    ax4 = plt.subplot2grid((3, 3), (2, 2))

    #Updating spectral registration
    pixx = np.arange(0,Measurement.NCONV[0],1)
    for i in range(Measurement.NGEOM):
        Measurement.VCONV[:,i] = construct_vconv(pixx,C0_fit+0.0,0.0,C1_fit,C2_fit)

    ix = 0
    SPECONV = np.zeros(Measurement.MEAS.shape)
    for i in range(Measurement.NGEOM):
        SPECONV[:,i] = yfit[ix:ix+Measurement.NCONV[i]]
        ix += Measurement.NCONV[i]

    offset = 0.
    for i in range(Measurement.NGEOM):
        ax1.plot(Measurement.VCONV[:,i],Measurement.MEAS[:,i]/SPECONV[:,i].max()+offset,c='black',linewidth=1.)
        offset += (1. - (Measurement.MEAS[:,i]/Measurement.MEAS[:,i].max()).min()) * 0.5

    offset = 0.
    for i in range(Measurement.NGEOM):
        ax1.plot(Measurement.VCONV[:,i],SPECONV[:,i]/SPECONV[:,i].max()+offset,linewidth=1.)
        offset += (1. - (Measurement.MEAS[:,i]/Measurement.MEAS[:,i].max()).min()) * 0.5

    for i in range(Measurement.NGEOM):
        offset = ((Measurement.MEAS[:,0]/Measurement.MEAS[:,0].max()).min())
        offset1 = (1. - (Measurement.MEAS[:,i]/Measurement.MEAS[:,i].max()).min()) * 0.5
        ax1.plot(Measurement.VCONV[:,i],(Measurement.MEAS[:,i]/SPECONV[:,i].max()-SPECONV[:,i]/SPECONV[:,i].max()) + offset - offset1 ,linewidth=1.)

    for i in range(Measurement.NGEOM):

        AMP1_fit = AAMP1_fit
        AMP2_fit = AAMP2_fit
        Measurement = modelILS(Measurement,DELDG_fit,FWHM_fit,AMP1_fit,AMP2_fit)
        Measurement.AFIL[np.where(Measurement.AFIL<0.0)] = 0.0

        ic = 0
        imax = int(Measurement.NFIL[ic]/2)
        ax2.plot(Measurement.VFIL[0:Measurement.NFIL[ic],ic]-Measurement.VFIL[imax,ic],Measurement.AFIL[0:Measurement.NFIL[ic],ic]/Measurement.AFIL[0:Measurement.NFIL[ic],ic].max(),linewidth=1.)
        ic = int(Measurement.NCONV[i]/2)
        imax = int(Measurement.NFIL[ic]/2)
        ax3.plot(Measurement.VFIL[0:Measurement.NFIL[ic],ic]-Measurement.VFIL[imax,ic],Measurement.AFIL[0:Measurement.NFIL[ic],ic]/Measurement.AFIL[0:Measurement.NFIL[ic],ic].max(),linewidth=1.)
        ic = Measurement.NCONV[i]-1
        imax = int(Measurement.NFIL[ic]/2)
        ax4.plot(Measurement.VFIL[0:Measurement.NFIL[ic],ic]-Measurement.VFIL[imax,ic],Measurement.AFIL[0:Measurement.NFIL[ic],ic]/Measurement.AFIL[0:Measurement.NFIL[ic],ic].max(),linewidth=1.)

    ax2.grid()
    ax3.grid()
    ax4.grid()
    ax2.set_facecolor('lightgray')
    ax3.set_facecolor('lightgray')
    ax4.set_facecolor('lightgray')
    ax2.set_xlabel('Wavenumber (cm$^{-1}$)')
    ax2.set_ylabel('Amplitude of ILS')
    ax3.set_xlabel('Wavenumber (cm$^{-1}$)')
    ax3.set_ylabel('Amplitude of ILS')
    ax4.set_xlabel('Wavenumber (cm$^{-1}$)')
    ax4.set_ylabel('Amplitude of ILS')

    ax2.set_xlim(-1.0,1.0)
    ax3.set_xlim(-1.0,1.0)
    ax4.set_xlim(-1.0,1.0)

    ax1.set_facecolor('lightgray')
    ax1.set_xlabel('Wavenumber (cm$^{-1}$)')
    ax1.set_ylabel('Transmission')
    plt.tight_layout()

    if os.path.exists(Caldir+'Summary_Plots/'+Window)==False:
        os.mkdir(Caldir+'Summary_Plots/'+Window)

    fig.savefig(Caldir+'Summary_Plots/'+Window+'/'+Observation+'_'+Window+'.png',dpi=300)
    
    #fig,ax1 = plt.subplots(1,1,figsize=(10,3))
    #for i in range(len(TANHESEL_AREOIDx)):
    #    ax1.plot(VCONV_COMB, meas_norm[:,i],linewidth=0.5)
    #fig.savefig(Caldir+'Summary_Plots/'+Window+'/'+Observation+'_'+Window+'_meas_norm.png',dpi=300)


    #Writing the results into a file for each row
    ##########################################################################################
    
    #Creating folder if it does not exist
    if os.path.exists(Caldir+'Calibration_Files/'+Window)==False:
        os.mkdir(Caldir+'Calibration_Files/'+Window)

    print('saving file for Row '+str(-1))

    #Creating folder if it does not exist
    if os.path.exists(Caldir+'Calibration_Files/'+Window+'/row'+str(-1))==False:
        os.mkdir(Caldir+'Calibration_Files/'+Window+'/row'+str(-1))
    
    filen = Caldir+'Calibration_Files/'+Window+'/row'+str(-1)+'/'+Observation+'.h5'

    AMP1_fit = AAMP1_fit
    AMP2_fit = AAMP2_fit
    Measurement = modelILS(Measurement,DELDG_fit,FWHM_fit,AMP1_fit,AMP2_fit)
    Measurement.AFIL[np.where(Measurement.AFIL<0.0)] = 0.0

    write_calfile(filen,Measurement.VCONV[:,i],DELDG_fit,FWHM_fit,AMP1_fit,AMP2_fit,DELDG_err=DELDG_fit_err,FWHM_err=FWHM_fit_err,AMP1_err=AAMP1_fit_err,AMP2_err=AAMP2_fit_err)

    finished = True
    return finished



##############################################################################################
##############################################################################################
# DEFINING FUNCTIONS FOR ILS MODELLING
##############################################################################################
##############################################################################################

def modelILS(Measurement,par4,par5,par6,par7):
    """
        FUNCTION NAME : modelILS()
        
        DESCRIPTION : Model de ILS from ACS MIR with a double Gaussian function.
        
        INPUTS :
        
            Measurement :: NEMESIS class for the Measurement, including the wavelength grid
            par4 :: Wavenumber offset of the second Gaussian with respect to the first one
            par5 :: FWHM of the first Gaussian
            par6 :: Amplitude of the second Gaussian with respect to the first one at the first wavelength
            par7 :: Amplitude of the second Gaussian with respect to the first one at the last wavelength
            
        OPTIONAL INPUTS: None
        
        OUTPUTS :
        
            Measurement :: NEMESIS class for the Measurement, including the ILS parameters
        
        CALLING SEQUENCE:
        
            Measurement = modelILS(Measurement,DELDG,FWHM,AMP1,AMP2)
        
        MODIFICATION HISTORY : Juan Alday (15/06/2024)
    """

    from archnemesis.helpers.maths_helper import ngauss

    par1 = 0.0
    par2 = 0.0
    par3 = 0.0
    #par4 = -par5

    #Calculating the parameters for each spectral point
    nconv = Measurement.NCONV[0]
    vconv1 = Measurement.VCONV[0:nconv,0]
    ng = 2

    # 1. Wavenumber offset of the two gaussians
    #    We divide it in two sections with linear polynomials     
    iconvmid = int(nconv/2.)
    wavemax = vconv1[nconv-1]
    wavemin = vconv1[0]
    wavemid = vconv1[iconvmid]
    offgrad1 = (par2 - par1)/(wavemid-wavemin)
    offgrad2 = (par2 - par3)/(wavemid-wavemax)
    offset = np.zeros([nconv,ng])
    for i in range(iconvmid):
        offset[i,0] = (vconv1[i] - wavemin) * offgrad1 + par1
        offset[i,1] = offset[i,0] + par4
    for i in range(nconv-iconvmid):
        offset[i+iconvmid,0] = (vconv1[i+iconvmid] - wavemax) * offgrad2 + par3
        offset[i+iconvmid,1] = offset[i+iconvmid,0] + par4

    # 2. FWHM for the two gaussians (assumed to be constant in wavelength, not in wavenumber)
    fwhm = np.zeros([nconv,ng])
    fwhml = par5 / wavemin**2.0
    for i in range(nconv):
        fwhm[i,0] = fwhml * (vconv1[i])**2.
        fwhm[i,1] = fwhm[i,0]

    # 3. Amplitde of the second gaussian with respect to the main one
    amp = np.zeros([nconv,ng])
    ampgrad = (par7 - par6)/(wavemax-wavemin)
    for i in range(nconv):
        amp[i,0] = 1.0
        amp[i,1] = (vconv1[i] - wavemin) * ampgrad + par6

    #Running for each spectral point
    nfil = np.zeros(nconv,dtype='int32')
    mfil1 = 200
    vfil1 = np.zeros([mfil1,nconv])
    afil1 = np.zeros([mfil1,nconv])
    for i in range(nconv):

        #determining the lowest and highest wavenumbers to calculate
        xlim = 0.0
        xdist = 5.0 
        for j in range(ng):
            xcen = offset[i,j]
            xmin = abs(xcen - xdist*fwhm[i,j]/2.)
            if xmin > xlim:
                xlim = xmin
            xmax = abs(xcen + xdist*fwhm[i,j]/2.)
            if xmax > xlim:
                xlim = xmax

        #determining the wavenumber spacing we need to sample properly the gaussians
        xsamp = 7.0   #number of points we require to sample one HWHM 
        xhwhm = 10000.0
        for j in range(ng):
            xhwhmx = fwhm[i,j]/2. 
            if xhwhmx < xhwhm:
                xhwhm = xhwhmx
        deltawave = xhwhm/xsamp
        np1 = 2.0 * xlim / deltawave
        npx = int(np1) + 1

        #Calculating the ILS in this spectral point
        iamp = np.zeros([ng])
        imean = np.zeros([ng])
        ifwhm = np.zeros([ng])
        fun = np.zeros([npx])
        xwave = np.linspace(vconv1[i]-deltawave*(npx-1)/2.,vconv1[i]+deltawave*(npx-1)/2.,npx)        
        for j in range(ng):
            iamp[j] = amp[i,j]
            imean[j] = offset[i,j] + vconv1[i]
            ifwhm[j] = fwhm[i,j]

        fun = ngauss(npx,xwave,ng,iamp,imean,ifwhm)  
        nfil[i] = npx
        vfil1[0:nfil[i],i] = xwave[:]
        afil1[0:nfil[i],i] = fun[:]

    mfil = nfil.max()
    vfil = np.zeros([mfil,nconv])
    afil = np.zeros([mfil,nconv])
    for i in range(nconv):
        vfil[0:nfil[i],i] = vfil1[0:nfil[i],i]
        afil[0:nfil[i],i] = afil1[0:nfil[i],i]
    
    Measurement.NFIL = nfil
    Measurement.VFIL = vfil
    Measurement.AFIL = afil

    return Measurement

##############################################################################################

def fm_calibrate_multi(Measurement, WAVE, TAUGAS_PATH, params):
    """
        FUNCTION NAME : fm_calibrate_multi()
        
        DESCRIPTION : Forward model to fit the instrument lineshape for several rows on the ACS MIR detector.
                        The ILS is modelled with a double Gaussian function. The amplitude of the second Gaussian
                        is a linear function of wavenumber. The amplitude of the second Gaussian is also a linear
                        function of the row number.
        
        INPUTS :
        
            Measurement :: NEMESIS class for the Measurement, including the wavelength grid
            WAVE        :: Wavelength grid for the spectroscopy (i.e., TAUGAS_PATH)
            TAUGAS_PATH :: Optical depth for each gas in the path
            params      :: Dictionary with the parameters for the ILS model:
                        DELDG  : Wavenumber offset of the second Gaussian with respect to the first one
                        FWHM   : FWHM of the first Gaussian
                        AAMP1  : Amplitude of the second Gaussian with respect to the first one at the first wavelength
                        AAMP2  : Amplitude of the second Gaussian with respect to the first one at the last wavelength
                        BAMP1  : Gradient of the amplitude of the second Gaussian as a function of the row number for the first wavelength
                        BAMP2  : Gradient of the amplitude of the second Gaussian as a function of the row number for the last wavelength
                        C0     : Wavenumber calibration parameter
                        C1     : Wavenumber calibration parameter
                        C2     : Wavenumber calibration parameter
                        SC_GASg_i : Scaling factor for each gas (g) in each row (i)
            
        OPTIONAL INPUTS: None
        
        OUTPUTS :
        
            Y :: Modelled spectra for each row on the detector
        
        CALLING SEQUENCE:
        
            Y = fm_calibrate_multi(Measurement, WAVE, TAUGAS_PATH, params)
        
        MODIFICATION HISTORY : Juan Alday (15/04/2025)
    """
    
    from archnemesis.Measurement_0 import lblconv_fil
    import numpy as np

    # Shared parameters
    DELDG = params['DELDG']
    FWHM  = params['FWHM']
    AAMP1  = params['AAMP1']
    AAMP2  = params['AAMP2']

    BAMP1  = params['BAMP1']
    BAMP2  = params['BAMP2']

    # Wavenumber calibration parameters for this geometry
    C0 = params['C0']
    C1 = params['C1']
    C2 = params['C2']

    modeled_spectra = []

    for i in range(Measurement.NGEOM):
        ngas = TAUGAS_PATH.shape[1]

        # Get scaling factors per gas for this geometry
        TAUTOT = np.zeros_like(TAUGAS_PATH[:, 0])
        for g in range(ngas):
            scale = params[f'SC_GAS{g+1}_{i+1}']
            TAUTOT += TAUGAS_PATH[:, g] * scale

        TRANS = np.exp(-TAUTOT)

        pixx = np.arange(0, Measurement.NCONV[i], 1)
        Measurement.VCONV[:, i] = construct_vconv(pixx, C0, 0.0, C1, C2)

        AMP1 = AAMP1 + BAMP1 * i  #Calculating the ILS for this specific row
        AMP2 = AAMP2 + BAMP2 * i  #Calculating the ILS for this specific row
        Measurement = modelILS(Measurement, DELDG, FWHM, AMP1, AMP2)

        TRANSCONV = lblconv_fil(
            len(WAVE),
            WAVE,
            TRANS,
            Measurement.NCONV[i],
            Measurement.VCONV[:, i],
            Measurement.NFIL,
            Measurement.VFIL,
            Measurement.AFIL,
        )

        # Baseline correction with 3rd-degree polynomial
        BASEMEAS = Measurement.MEAS[:, i] / TRANSCONV
        PCOEFF = np.polyfit(Measurement.VCONV[:, i], BASEMEAS, 3)
        TRANS_BASE = np.polyval(PCOEFF, Measurement.VCONV[:, i])
        TRANSCONV *= TRANS_BASE

        modeled_spectra.append(TRANSCONV)

    Y = np.zeros(np.sum(Measurement.NCONV))
    ix = 0
    for i in range(Measurement.NGEOM):
        Y[ix:ix+Measurement.NCONV[i]] = modeled_spectra[i]
        ix = ix + Measurement.NCONV[i]
    
    if np.where(np.isnan(Y)==True)[0].size>0:
        print(params['C0'].value,params['C1'].value,params['C2'].value)
        print(params['DELDG'].value,params['FWHM'].value,params['AAMP1'].value,params['AAMP2'].value)
        print(params['BAMP1'].value,params['BAMP2'].value)
        raise ValueError('There are NaN values in the modelled spectra')
    
    return Y


##############################################################################################

def fm_calibrate(Measurement, WAVE, TAUGAS_PATH, params):
    """
        FUNCTION NAME : fm_calibrate()
        
        DESCRIPTION : Forward model to fit the instrument lineshape for several rows on the ACS MIR detector.
                        The ILS is modelled with a double Gaussian function. The amplitude of the second Gaussian
                        is a linear function of wavenumber.
        
        INPUTS :
        
            Measurement :: NEMESIS class for the Measurement, including the wavelength grid
            WAVE        :: Wavelength grid for the spectroscopy (i.e., TAUGAS_PATH)
            TAUGAS_PATH :: Optical depth for each gas in the path
            params      :: Dictionary with the parameters for the ILS model:
                        DELDG  : Wavenumber offset of the second Gaussian with respect to the first one
                        FWHM   : FWHM of the first Gaussian
                        AAMP1  : Amplitude of the second Gaussian with respect to the first one at the first wavelength
                        AAMP2  : Amplitude of the second Gaussian with respect to the first one at the last wavelength
                        C0     : Wavenumber calibration parameter
                        C1     : Wavenumber calibration parameter
                        C2     : Wavenumber calibration parameter
                        SC_GASg_i : Scaling factor for each gas (g) in each row (i)
            
        OPTIONAL INPUTS: None
        
        OUTPUTS :
        
            Y :: Modelled spectra for each row on the detector
        
        CALLING SEQUENCE:
        
            Y = fm_calibrate_multi(Measurement, WAVE, TAUGAS_PATH, params)
        
        MODIFICATION HISTORY : Juan Alday (15/04/2025)
    """
    
    from archnemesis.Measurement_0 import lblconv_fil
    import numpy as np

    # Shared parameters
    DELDG = params['DELDG']
    FWHM  = params['FWHM']
    AAMP1  = params['AAMP1']
    AAMP2  = params['AAMP2']

    # Wavenumber calibration parameters for this geometry
    C0 = params['C0']
    C1 = params['C1']
    C2 = params['C2']

    modeled_spectra = []

    for i in range(Measurement.NGEOM):
        ngas = TAUGAS_PATH.shape[1]

        # Get scaling factors per gas for this geometry
        TAUTOT = np.zeros_like(TAUGAS_PATH[:, 0])
        for g in range(ngas):
            scale = params[f'SC_GAS{g+1}_{i+1}']
            TAUTOT += TAUGAS_PATH[:, g] * scale

        TRANS = np.exp(-TAUTOT)

        pixx = np.arange(0, Measurement.NCONV[i], 1)
        Measurement.VCONV[:, i] = construct_vconv(pixx, C0, 0.0, C1, C2)

        AMP1 = AAMP1
        AMP2 = AAMP2
        Measurement = modelILS(Measurement, DELDG, FWHM, AMP1, AMP2)

        TRANSCONV = lblconv_fil(
            len(WAVE),
            WAVE,
            TRANS,
            Measurement.NCONV[i],
            Measurement.VCONV[:, i],
            Measurement.NFIL,
            Measurement.VFIL,
            Measurement.AFIL,
        )

        # Baseline correction with 3rd-degree polynomial
        BASEMEAS = Measurement.MEAS[:, i] / TRANSCONV
        PCOEFF = np.polyfit(Measurement.VCONV[:, i], BASEMEAS, 3)
        TRANS_BASE = np.polyval(PCOEFF, Measurement.VCONV[:, i])
        TRANSCONV *= TRANS_BASE

        modeled_spectra.append(TRANSCONV)

    Y = np.zeros(np.sum(Measurement.NCONV))
    ix = 0
    for i in range(Measurement.NGEOM):
        Y[ix:ix+Measurement.NCONV[i]] = modeled_spectra[i]
        ix = ix + Measurement.NCONV[i]
    
    if np.where(np.isnan(Y)==True)[0].size>0:
        print(params['C0'].value,params['C1'].value,params['C2'].value)
        print(params['DELDG'].value,params['FWHM'].value,params['AAMP1'].value,params['AAMP2'].value)
        raise ValueError('There are NaN values in the modelled spectra')
    
    return Y

####################################################################################################

def calc_tau_gas(Atmosphere,Measurement,Spectroscopy):
    """
        FUNCTION NAME : fm_calibrate_multi()
        
        DESCRIPTION : Calculate the optical depth for each gas in the path
        
        INPUTS :
        
            Atmosphere :: NEMESIS class for the Atmosphere, including the gas concentrations
            Measurement :: NEMESIS class for the Measurement, including the wavelength grid
            Spectroscopy :: NEMESIS class for the Spectroscopy, including the spectroscopy data
            
        OPTIONAL INPUTS: None
        
        OUTPUTS :
        
            WAVE(NWAVE) :: Wavelength grid for the spectroscopy (i.e., TAUGAS_PATH)
            TAUGAS_PATH(NWAVE,NGAS) :: Optical depth for each gas in the path
        
        CALLING SEQUENCE:
        
            WAVE,TAUGAS_PATH = calc_tau_gas(Atmosphere,Measurement,Spectroscopy)
        
        MODIFICATION HISTORY : Juan Alday (15/04/2025)
    """

    #Defining other classes
    ##########################################################################################
    
    #Defining Scatter class
    Scatter = ans.Scatter_0(NDUST=1,ISPACE=0,ISCAT=0)
    Scatter.NWAVE = 2
    Scatter.WAVE = np.linspace(0.,10000.,Scatter.NWAVE)
    Scatter.KEXT = np.ones((Scatter.NWAVE,Scatter.NDUST))
    Scatter.SGLALB = np.ones((Scatter.NWAVE,Scatter.NDUST))

    #Defining Layer class
    Layer = ans.Layer_0(LAYTYP=5,LAYINT=1,LAYHT=0.0)
    Layer.NLAY = Atmosphere.NP - 1
    Layer.H_base = Atmosphere.H[0:Atmosphere.NP-1]
    Atmosphere.calc_grav()
    Atmosphere.calc_radius()
    Layer.RADIUS = Atmosphere.RADIUS

    #Defining Stellar class
    Stellar = ans.Stellar_0(SOLEXIST=False,ISPACE=0)

    #Defining Surface class
    Surface = ans.Surface_0(GALB=0.0,LOWBC=0)
    Surface.NEM = 2
    Scatter.VEM = np.linspace(0.,10000.,Scatter.NWAVE)
    Surface.EMISSIVITY = np.ones(Surface.NEM)

    #Defining the retrieval parameters
    Retrieval = ans.OptimalEstimation_0(PHILIMIT=0.1,NITER=10,NCORES=1,IRET=0)
    
    #Calculating the optical depth for each gas in the path
    FM = ans.ForwardModel_0(Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Scatter=Scatter,CIA=None,Layer=Layer)

    FM.calc_path_SO()
    
    if np.where(np.isnan(FM.LayerX.AMOUNT)==True)[0].size>0:
        print('There are NaN values in the amounts')
        raise ValueError('There are NaN values in the amounts')

    if np.where(np.isnan(FM.PathX.SCALE)==True)[0].size>0:
        print('There are NaN values in the path scale')
        raise ValueError('There are NaN values in the path scale')

    TAUGAS = np.zeros((FM.SpectroscopyX.NWAVE,FM.SpectroscopyX.NG,FM.LayerX.NLAY,FM.SpectroscopyX.NGAS))  #Vertical opacity of each gas in each layer

    #Calculating the cross sections for each gas in each layer
    k = FM.SpectroscopyX.calc_klbl(FM.LayerX.NLAY,FM.LayerX.PRESS/101325.,FM.LayerX.TEMP)

    for i in range(FM.SpectroscopyX.NGAS):
        IGAS = np.where( (FM.AtmosphereX.ID==FM.SpectroscopyX.ID[i]) & (FM.AtmosphereX.ISO==FM.SpectroscopyX.ISO[i]) )[0]

        #Calculating vertical column density in each layer
        VLOSDENS = FM.LayerX.AMOUNT[:,IGAS].T * 1.0e-4   #cm-2

        #Calculating vertical opacity for each gas in each layer
        TAUGAS[:,0,:,i] = k[:,:,i] * VLOSDENS   #(NWAVE,NG,NLAY,NGAS)

    TAUGAS = np.transpose(TAUGAS,axes=(0,3,1,2)) #(NWAVE,NGAS,NG,NLAY)

    #Calculating the line-of-sight opacities
    TAUGAS_LAYINC = TAUGAS[:,:,:,FM.PathX.LAYINC[:,:]] * FM.PathX.SCALE[:,:]  #(NWAVE,NGAS,NG,NLAYIN,NPATH)

    #Calculating the total opacity over the path
    TAUGAS_PATH = np.sum(TAUGAS_LAYINC,3) #(NWAVE,NGAS,NG,NPATH)
    TAUGAS_PATH = np.mean(TAUGAS_PATH,axis=3)
    TAUGAS_PATH_ini = TAUGAS_PATH[:,:,0]
    TRANS_ini = np.exp(-np.sum(TAUGAS_PATH,axis=1))
    
    return FM.SpectroscopyX.WAVE,TAUGAS_PATH_ini
    
##############################################################################################

def write_calfile(filen,waveconv,DELDG,FWHM,AMP1,AMP2,DELDG_err=None,FWHM_err=None,AMP1_err=None,AMP2_err=None):
    '''
    Write the calibration results into HDF5 file
    '''

    import h5py
    
    hf = h5py.File(filen, 'w') 
    
    #Writing the wavenumber array
    dset = hf.create_dataset('VCONV_CALIBRATED',data=waveconv)
    dset.attrs['title'] = "Calibrated wavenumber array for the spectral window"

    #Writing the ILS parameters
    dset = hf.create_dataset('DELDG',data=DELDG)
    dset.attrs['title'] = "Distance of the first and second Gaussian (cm-1)"

    if DELDG_err is not None:
        dset.attrs['error'] = DELDG_err
        dset.attrs['error_title'] = "Error in the distance of the first and second Gaussian (cm-1)"

    dset = hf.create_dataset('FWHM',data=FWHM)
    dset.attrs['title'] = "Full width at half maximum for the main Gaussian (cm-1)"

    if FWHM_err is not None:
        dset.attrs['error'] = FWHM_err
        dset.attrs['error_title'] = "Error in the full width at half maximum for the main Gaussian (cm-1)"

    dset = hf.create_dataset('AMP1',data=AMP1)
    dset.attrs['title'] = "Amplitude of second Gaussian wrt the main one at first wavenumber"

    if AMP1_err is not None:
        dset.attrs['error'] = AMP1_err
        dset.attrs['error_title'] = "Error in the amplitude of second Gaussian wrt the main one at first wavenumber"

    dset = hf.create_dataset('AMP2',data=AMP2)
    dset.attrs['title'] = "Amplitude of second Gaussian wrt the main one at last wavenumber"

    if AMP2_err is not None:
        dset.attrs['error'] = AMP2_err
        dset.attrs['error_title'] = "Error in the amplitude of second Gaussian wrt the main one at last wavenumber"

    hf.close()
    
##############################################################################################

def read_calfile(filen):
    '''
    Read the calibration results from HDF5 file
    '''

    import h5py
    
    hf = h5py.File(filen, 'r') 
    WAVECONV = np.array(hf.get('VCONV_CALIBRATED'))
    DELDG = np.array(hf.get('DELDG'))
    FWHM = np.array(hf.get('FWHM'))
    AMP1 = np.array(hf.get('AMP1'))
    AMP2 = np.array(hf.get('AMP2'))
    hf.close()

    return WAVECONV,DELDG,FWHM,AMP1,AMP2

##############################################################################################

def construct_dVm(pixx,C0,C1,C2):
    
    dVm = C0 + C1*pixx + C2*(pixx)**2.
    
    return dVm
        
##############################################################################################

def construct_vconv(pixx,v0,C0,C1,C2):
    
    iCENTER = np.argmin(np.abs(pixx))
    
    dVm = construct_dVm(pixx,C0,C1,C2)
    
    vconv = np.zeros(len(pixx))
    vconv[iCENTER] = v0
    vconv[:] = v0 + dVm
    
    return vconv

##############################################################################################

def correlate_spectra(pixx,v0,C1,C2,meas,smodel,normalise=True):
    
    newconv = construct_vconv(pixx,v0,0.0,C1,C2)
    model = smodel(newconv)
    
    #Fitting polynomial function to both measured and modelled spectra to normalise them
    if normalise==True:
        ndeg = 2
        px = np.polyfit(pixx,meas,ndeg)
        pol = px[2] + px[1]*pixx + px[0]*(pixx)**2.

        px = np.polyfit(pixx,model,ndeg)
        pol2 = px[2] + px[1]*pixx + px[0]*(pixx)**2.
    
    corr = np.correlate(meas-pol,model-pol2)[0]
    
    return corr

##############################################################################################
##############################################################################################
# DEFINING FUNCTIONS FOR BASELINE SUBSTRACTION
##############################################################################################
##############################################################################################

from scipy.linalg import solveh_banded

def als_baseline(intensities, asymmetry_param=0.05, smoothness_param=1e6,
                 max_iters=10, conv_thresh=1e-5, verbose=False):
    '''Computes the asymmetric least squares baseline.
    * http://www.science.uva.nl/~hboelens/publications/draftpub/Eilers_2005.pdf
    smoothness_param: Relative importance of smoothness of the predicted response.
    asymmetry_param (p): if y > z, w = p, otherwise w = 1-p.
                         Setting p=1 is effectively a hinge loss.
    '''
    smoother = WhittakerSmoother(intensities, smoothness_param, deriv_order=2)
    # Rename p for concision.
    p = asymmetry_param
    # Initialize weights.
    w = np.ones(intensities.shape[0])
    for i in range(max_iters):
        z = smoother.smooth(w)
        mask = intensities > z
        new_w = p * mask + (1-p) * (~mask)
        conv = np.linalg.norm(new_w - w)
        if verbose:
            print(i + 1, conv)
        if conv < conv_thresh:
            break
        w = new_w
    else:
        print('ALS did not converge in %d iterations' % max_iters)
    return z

class WhittakerSmoother(object):
    def __init__(self, signal, smoothness_param, deriv_order=1):
        self.y = signal
        assert deriv_order > 0, 'deriv_order must be an int > 0'
        # Compute the fixed derivative of identity (D).
        d = np.zeros(deriv_order * 2 + 1, dtype=int)
        d[deriv_order] = 1
        d = np.diff(d, n=deriv_order)
        n = self.y.shape[0]
        k = len(d)
        s = float(smoothness_param)

        # Here be dragons: essentially we're faking a big banded matrix D,
        # doing s * D.T.dot(D) with it, then taking the upper triangular bands.
        diag_sums = np.vstack([
            np.pad(s * np.cumsum(d[-i:] * d[:i]), ((k-i, 0),), 'constant')
            for i in range(1, k + 1)])
        upper_bands = np.tile(diag_sums[:, -1:], n)
        upper_bands[:, :k] = diag_sums
        for i, ds in enumerate(diag_sums):
            upper_bands[i, -i-1:] = ds[::-1][:i+1]
        self.upper_bands = upper_bands

    def smooth(self, w):
        foo = self.upper_bands.copy()
        foo[-1] += w  # last row is the diagonal
        return solveh_banded(foo, w * self.y, overwrite_ab=True, overwrite_b=True)

