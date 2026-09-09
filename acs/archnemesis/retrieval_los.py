#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# retrievals_los - Set of functions to run an ACS retrieval of the line-of-sight densities
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
from struct import *
import acs
import archnemesis as ans
from copy import deepcopy
import h5py, scipy

###############################################################################################

def run_retrieval(outname,
                  Atmosphere,Measurement,Spectroscopy,ID_ret,ISO_ret,
                  ndegree=2,
                  fourierdegree=None,
                  fit_baseline_automatic=False,
                  maxh=None,
                  use_errmeas=True):
    
    """
    FUNCTION NAME : retrieval_LOS_all()

    DESCRIPTION : Perform the retrieval of solar occultation spectra using the modified method. 
                  In this retrieval we analyse each spectrum independently and fit only the line-of-sight density
                  of some gases (together with some polynomial for the baseline and the ILS).
                  
                  In this version, we call several times retrieval_LOS() to perform the retrieval independently
                  in each tangent height.

    INPUTS : 

        Atmosphere :: archNEMESIS Atmosphere class
        Measurement :: archNEMESIS Measurement class
        Spectroscopy :: archNEMESIS Spectroscopy class
        ID_ret(NGAS) :: Gas ID of the species to retrieve 
        ISO_ret(NGAS) :: Isotope ID of the species to retrieve

    OPTIONAL INPUTS:
    
 
        ndegree :: Degree of the polynomial to fit for the baseline
        maxh :: Maximum altitude to which perform the retrievals (If None, all geometries are used)
        fit_baseline_automatic :: The baseline is fitted in each iteration of the forward model using ndegree and fourierdegree
            
    OUTPUTS : 

    CALLING SEQUENCE:

        run_retrieval_LOS_all(dir_atm,runname_atm,dir_meas,runname_meas)

    MODIFICATION HISTORY : Juan Alday (09/08/2023)

    """

    ID_ret = np.array(ID_ret)
    ISO_ret = np.array(ISO_ret)

    #Reading measurement file to get array sizes
    ##################################################################
    
    if maxh is None:
        NGEOM = Measurement.NGEOM
        iin = np.linspace(0,NGEOM-1,NGEOM,dtype='int32')
    else:
        iin = np.where(Measurement.TANHE[:,0]<=maxh)[0]
        NGEOM = len(iin)
        
        
    NCONV = Measurement.NCONV[0]
    SPECMOD = np.zeros((NCONV,NGEOM))
    BASELINE = np.zeros((NCONV,NGEOM))
    LOSDENS_TOT = np.zeros(NGEOM)
    VMR_ret = np.zeros((NGEOM,Spectroscopy.NGAS))
    VMR_reterr = np.zeros((NGEOM,Spectroscopy.NGAS))

    #Checking whether the measurement uncertainty has NaN values
    ##################################################################
    
    if use_errmeas is True:
        if len(np.where(np.isnan(Measurement.ERRMEAS)==True)[0])>0:
            print('warning :: use_errmeas set to False because there are NaN values in ERRMEAS')
            use_errmeas=False

    #Performing the retrieval in each tangent height
    ##################################################################
    
    SPECMOD,SPECMOD_sigma,BASELINE,LOSDENS_TOT,VMR_ret,VMR_reterr = retrieval_LOS_all(Atmosphere,Measurement,Spectroscopy,ID_ret,ISO_ret,ndegree=ndegree,fourierdegree=fourierdegree,fit_baseline_automatic=fit_baseline_automatic,use_errmeas=use_errmeas)
    VMR_ret = VMR_ret.T
    VMR_reterr = VMR_reterr.T


    #Writing output file
    ##################################################################
    
    f = h5py.File(outname+'.h5','w')
    
    #Writing the main dimensions
    dset = f.create_dataset('NGEOM',data=NGEOM)
    dset.attrs['title'] = "Number of tangent heights"
    
    dset = f.create_dataset('NCONV',data=Measurement.NCONV[0])
    dset.attrs['title'] = "Number of convolution wavelengths"
    
    dset = f.create_dataset('TANHE',data=Measurement.TANHE[iin,0])
    dset.attrs['title'] = "Tangent height / km"
    
    dset = f.create_dataset('VCONV',data=Measurement.VCONV[0:Measurement.NCONV[0],0])
    dset.attrs['title'] = "Convolution wavelengths"
    
    dset = f.create_dataset('MEAS',data=Measurement.MEAS[0:Measurement.NCONV[0],iin])
    dset.attrs['title'] = "Measured spectra"
    
    dset = f.create_dataset('ERRMEAS',data=Measurement.ERRMEAS[0:Measurement.NCONV[0],iin])
    dset.attrs['title'] = "Uncertainty in measured spectra"
    
    dset = f.create_dataset('SPECMOD_trace',data=SPECMOD)
    dset.attrs['title'] = "Modelled spectra with trace gases"
    
    dset = f.create_dataset('SPECMOD_sigma',data=SPECMOD_sigma)
    dset.attrs['title'] = "Modelled spectra with 1,3 and 5-sigma contribution from trace gas"
    
    dset = f.create_dataset('SPECMOD_background',data=SPECMOD)
    dset.attrs['title'] = "Modelled spectra with background gases (w/o trace gases)"
    
    dset = f.create_dataset('BASELINE',data=BASELINE)
    dset.attrs['title'] = "Modelled baseline"
    
    
    #Writing the retrieved parameters
    dset = f.create_dataset('ID',data=np.array(ID_ret))
    dset.attrs['title'] = "Gas ID of fitted species"
    
    dset = f.create_dataset('ISO',data=np.array(ISO_ret))
    dset.attrs['title'] = "Isotope ID of fitted species"
    
    dset = f.create_dataset('LOSDENS_TOT',data=LOSDENS_TOT)
    dset.attrs['title'] = "Line-of-sight density at each tangent height (m-2)"
    
    dset = f.create_dataset('VMR_ret',data=VMR_ret)
    dset.attrs['title'] = "Volume mixing ratio of fitted species"
    
    dset = f.create_dataset('VMR_reterr',data=VMR_reterr)
    dset.attrs['title'] = "Uncertainty in volume mixing ratio of fitted species"
    
    f.close()


###############################################################################################

def retrieval_LOS_all(Atmosphere,Measurement,Spectroscopy,
                      ID_ret,ISO_ret,
                      ndegree=2,
                      fourierdegree=None,
                      fit_baseline_automatic=False,
                      use_errmeas=True):
    
    """
    FUNCTION NAME : retrieval_LOS_all()

    DESCRIPTION : Perform the retrieval of solar occultation spectra using the modified method. 
                  In this retrieval we analyse each spectrum independently and fit only the line-of-sight density
                  of some gases (together with some polynomial for the baseline)

    INPUTS : 

        RunName :: Name of the NEMESIS file

    OPTIONAL INPUTS:
            
    OUTPUTS : 

    CALLING SEQUENCE:

        retrieval_LOS(RunName)

    MODIFICATION HISTORY : Juan Alday (09/08/2023)

    """
    
    from lmfit import Parameters,Model,Minimizer

    ######################################################
    ######################################################
    #    READING INPUT FILES AND SETTING UP VARIABLES
    ######################################################
    ######################################################

    #Checking that retrieved gas exists in Spectroscopy and Atmosphere
    ###############################################################

    #Reading the header information
    Spectroscopy.read_header()

    #Adding the gas to be searched for
    NGAS_ret = len(ID_ret)
    vary_gas = [False] * Spectroscopy.NGAS
    for igas in range(NGAS_ret):

        igasx = np.where( (Spectroscopy.ID==ID_ret[igas]) & (Spectroscopy.ISO==ISO_ret[igas]) )[0]
        if len(igasx)==0:
            print(ID_ret[igas],ISO_ret[igas])
            print(Specotroscopy.ID)
            raise ValueError("error :: retrieved gas does not exist in Spectroscopy class")
        if len(igasx)>1:
            print(ID_ret[igas],ISO_ret[igas])
            print(Specotroscopy.ID)
            raise ValueError("error :: there are multiple indices in Spectroscopy class corresponding to one retrieved gas")
        if len(igasx)==1:
            vary_gas[igasx[0]] = True

        igasx = np.where( (Atmosphere.ID==ID_ret[igas]) & (Atmosphere.ISO==ISO_ret[igas]) )[0]
        if len(igasx)==0:
            print(ID_ret[igas],ISO_ret[igas])
            print(Atmosphere.ID,Atmosphere.ISO)
            raise ValueError("error :: retrieved gas does not exist in Atmosphere class")
        if len(igasx)>1:
            print(ID_ret[igas],ISO_ret[igas])
            print(Atmosphere.ID,Atmosphere.ISO)
            raise ValueError("error :: there are multiple indices in Atmosphere class corresponding to one retrieved gas")


    #Initialise Measurement class and reading LBL-tables file
    ###############################################################

    Measurement.calc_MeasurementVector()

    #Calculating the 'calculation wavelengths'
    vmin,vmax = Measurement.calc_wave_range()

    #Now, reading k-tables or lbl-tables for the spectral range of interest
    Spectroscopy.read_tables(wavemin=vmin,wavemax=vmax)

    #Initialise other classes
    ###############################################################

    #Scatter class
    Scatter = None

    #Layer class
    Layer = ans.Layer_0(LAYTYP=5,LAYINT=1,LAYHT=Atmosphere.H.min(),RADIUS=Atmosphere.RADIUS)
    Layer.NLAY = Atmosphere.NP - 1
    Layer.H_base = Atmosphere.H[0:Atmosphere.NP-1]
    Layer.assess()

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

    #CIA class
    CIA = None

    #Retrieval class
    Retrieval = ans.OptimalEstimation_0(IRET=0)
    Retrieval.NITER = -1       #Number of iterations
    Retrieval.PHILIMIT = 0.1   #Convergence criterion
    Retrieval.NCORES = 1       #Number of available cores
    Retrieval.assess_input()

    Variables = None

    ######################################################
    ######################################################
    #           CALCULATING LOS DENSITY
    ######################################################
    ######################################################

    #Firstly, we calculate the line-of-sight densities
    LOSDENS_PATH = calc_los_density_all(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval)  #(NGAS,NGEOM)
    
    LOSDENS_ACT_PATH = np.zeros((Spectroscopy.NGAS,Measurement.NGEOM))
    for igas in range(Spectroscopy.NGAS):
        igasx = np.where( (Atmosphere.ID==Spectroscopy.ID[igas]) & (Atmosphere.ISO==Spectroscopy.ISO[igas]) )[0][0]
        LOSDENS_ACT_PATH[igas,:] = LOSDENS_PATH[igasx,:]
    
    LOSDENS_TOT_PATH = np.sum(LOSDENS_PATH,axis=0)  #(NGEOM)
    VMR_ACT_PATH = LOSDENS_ACT_PATH / LOSDENS_TOT_PATH  #(NGAS,NGEOM)
    
    #Secondly, we calculate the high-spectral resolution optical depth in the whole spectral range
    TAUGAS_PATH = calc_los_tau_all(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval)  #(NWAVE,NGAS,NGEOM)

    ######################################################
    ######################################################
    #       PREPARING MODEL AND PERFORMING RETRIEVAL
    ######################################################
    ######################################################
    
    igas_vary = np.flatnonzero(vary_gas)
    ngas_vary = len(igas_vary)


    SPECMOD = np.zeros(Measurement.MEAS.shape)
    SPECMODx = np.zeros(Measurement.MEAS.shape)
    BASELINEx = np.zeros(Measurement.MEAS.shape)
    SPECMOD_sigma = np.full((*SPECMODx.shape, ngas_vary, 3), np.nan)
    
    VMR_ret = np.full((ngas_vary, Measurement.NGEOM), np.nan)
    VMR_reterr = np.full((ngas_vary, Measurement.NGEOM), np.nan)
    
    for iGEOM in range(Measurement.NGEOM):
        
        print('Tangent height ',iGEOM+1,'of',Measurement.NGEOM)
    
        #First we perform retrieval without last gas
        #################################################
        
        lmodel = Model(forward_model,independent_vars=['Measurement','WAVE','TAUGAS_PATH','ndegree','fourierdegree'])
        
        for i in range(Spectroscopy.NGAS):
            lmodel.set_param_hint(
                f"SC_GAS{i + 1}",
                value=1.0,
                vary=vary_gas[i],
            )
        
        if fit_baseline_automatic==True:
            lmodel.set_param_hint('A0',vary=False)
            lmodel.set_param_hint('A1',vary=False)
            lmodel.set_param_hint('A2',vary=False)
            lmodel.set_param_hint('A3',vary=False)
            A0ini = 0.0 ; A1ini = 0.0 ; A2ini = 0.0 ;A3ini = 0.0
        else:

            for i in range(ndegree + 1):
                lmodel.set_param_hint(
                    f"A{i}",
                    value=1.0 if i == 0 else 0.0,
                    vary=True,
                )

            params = lmodel.make_params()

        if use_errmeas==True:
            result = lmodel.fit(Measurement.MEAS[:,iGEOM],params=params,Measurement=Measurement,WAVE=Spectroscopy.WAVE,TAUGAS_PATH=TAUGAS_PATH[:,:,iGEOM],ndegree=ndegree,fourierdegree=fourierdegree,weights=1./Measurement.ERRMEAS[:,0]**2.)
        else:
            result = lmodel.fit(Measurement.MEAS[:,iGEOM],params=params,Measurement=Measurement,WAVE=Spectroscopy.WAVE,TAUGAS_PATH=TAUGAS_PATH[:,:,iGEOM],ndegree=ndegree,fourierdegree=fourierdegree)
        

        # Extract retrieved parameters for varied gases only
        gas_params = [
            result.params[f"SC_GAS{i + 1}"]
            for i in igas_vary
        ]

        SC_GAS_ret = np.array([p.value for p in gas_params])

        SC_GAS_reterr = np.array([
            p.stderr if p.stderr is not None else np.nan
            for p in gas_params
        ])

        # Best-fit spectrum
        SPECMODx[:, iGEOM] = result.best_fit

        # Retrieved abundances and uncertainties for varied gases only
        vmr_reference = VMR_ACT_PATH[igas_vary, iGEOM]

        VMR_ret[:, iGEOM] = SC_GAS_ret * vmr_reference
        VMR_reterr[:, iGEOM] = SC_GAS_reterr * np.abs(vmr_reference)

        # Baseline coefficients
        A_ret = np.array([
            result.params[f"A{i}"].value
            for i in range(ndegree + 1)
        ])

        # Calculate spectra at different sigma levels
        tau = TAUGAS_PATH[:, :, iGEOM]
        ngas = tau.shape[1]

        model_inputs = dict(
            Measurement=Measurement,
            WAVE=Spectroscopy.WAVE,
            TAUGAS_PATH=tau,
            ndegree=ndegree,
            fourierdegree=fourierdegree,
        )

        best_parameters = result.params.valuesdict()

        # Perturb only varied gases; keep all other parameters at best-fit values
        for j, igas in enumerate(igas_vary):
            name = f"SC_GAS{igas + 1}"
            uncertainty = result.params[name].stderr

            if uncertainty is None or not np.isfinite(uncertainty):
                continue

            for isigma, nsigma in enumerate([1, 3, 5]):
                perturbed = best_parameters.copy()
                perturbed[name] += nsigma * uncertainty

                SPECMOD_sigma[:, iGEOM, j, isigma] = forward_model(
                    **model_inputs, **perturbed
                )

        # Baseline only: zero all gases, including fixed gases
        baseline_parameters = best_parameters.copy()
        for igas in range(ngas):
            baseline_parameters[f"SC_GAS{igas + 1}"] = 0.0

        BASELINEx[:, iGEOM] = forward_model(
            **model_inputs, **baseline_parameters
        )

    return SPECMODx,SPECMOD_sigma,BASELINEx,LOSDENS_TOT_PATH,VMR_ret,VMR_reterr


###############################################################################################

def calc_los_density(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval):
    
    """
    FUNCTION NAME : calc_los_density()

    DESCRIPTION : Calculate the line-of-sight density of all active gases 

    INPUTS : 

        NEMESIS python classes

    OPTIONAL INPUTS:
            
    OUTPUTS : 
    
        los_dens_gas(ngas) :: Line-of-sight density of all gases in atmosphere (m-2)

    CALLING SEQUENCE:

        los_dens_gas = calc_los_density(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval)

    MODIFICATION HISTORY : Juan Alday (09/08/2023)

    """
    
    #Initialise forward model class to perform calculations
    FM = ans.ForwardModel_0(runname='retrieval_LOS', Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables)

    #Initialising the X classes just for storing one location
    FM.MeasurementX = deepcopy(FM.Measurement)
    FM.AtmosphereX = deepcopy(FM.Atmosphere)
    FM.SurfaceX = deepcopy(FM.Surface)
    FM.ScatterX = deepcopy(FM.Scatter)
    FM.Scatter = None
    FM.StellarX = deepcopy(FM.Stellar)
    FM.Stellar = None
    FM.SpectroscopyX = deepcopy(FM.Spectroscopy)
    FM.Spectroscopy = None
    FM.LayerX = deepcopy(FM.Layer)
    FM.CIAX = deepcopy(FM.CIA)
    FM.CIA = None
    flagh2p = False
        
    #Calculating the path
    FM.calc_path_SO()
    
    BASEH_TANHE = np.zeros(FM.PathX.NPATH)
    for i in range(FM.PathX.NPATH):
        BASEH_TANHE[i] = FM.LayerX.BASEH[FM.PathX.LAYINC[int(FM.PathX.NLAYIN[i]/2),i]]/1.0e3
    
    #Calculating the line-of-sight densities in each layer
    LOSDENS_LAYINC = FM.LayerX.AMOUNT[FM.PathX.LAYINC[:,:],:].T * FM.PathX.SCALE[:,:].T  #(NGAS,NPATH,NLAYIN)
        
    #Calculating the line-of-sight over the path
    LOSDENS_PATH = np.sum(LOSDENS_LAYINC,2) #(NGAS,NPATH)
    
    #Calculating the LOS density exactly at our tangent height
    if FM.PathX.NPATH!=2:
        sys.exit('error :: Number of calculated paths in calc_LOS_density() must be 2 (only one tangent height in Measurement)')
    LOSDENS_PATHX = np.zeros(FM.AtmosphereX.NVMR)
    LOSDENS_PATHX[:] = LOSDENS_PATH[:,0] + (LOSDENS_PATH[:,1]-LOSDENS_PATH[:,0])/(BASEH_TANHE[1]-BASEH_TANHE[0]) * (FM.MeasurementX.TANHE[0,0]-BASEH_TANHE[0]) 
    
    return LOSDENS_PATHX


###############################################################################################

def calc_los_density_all(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval):
    
    """
    FUNCTION NAME : calc_los_density()

    DESCRIPTION : Calculate the line-of-sight density of all active gases 

    INPUTS : 

        NEMESIS python classes

    OPTIONAL INPUTS:
            
    OUTPUTS : 
    
        los_dens_gas(ngas) :: Line-of-sight density of all gases in atmosphere (m-2)

    CALLING SEQUENCE:

        los_dens_gas = calc_los_density(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval)

    MODIFICATION HISTORY : Juan Alday (09/08/2023)

    """
    
    #Initialise forward model class to perform calculations
    FM = ans.ForwardModel_0(runname='retrieval_LOS', Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables)

    #Initialising the X classes just for storing one location
    FM.MeasurementX = deepcopy(FM.Measurement)
    FM.AtmosphereX = deepcopy(FM.Atmosphere)
    FM.SurfaceX = deepcopy(FM.Surface)
    FM.ScatterX = deepcopy(FM.Scatter)
    FM.Scatter = None
    FM.StellarX = deepcopy(FM.Stellar)
    FM.Stellar = None
    FM.SpectroscopyX = deepcopy(FM.Spectroscopy)
    FM.Spectroscopy = None
    FM.LayerX = deepcopy(FM.Layer)
    FM.CIAX = deepcopy(FM.CIA)
    FM.CIA = None
    flagh2p = False
        
    #Calculating the path
    FM.calc_path_SO()
    
    #Calculating the base height of each layer (at which the paths are calculated)
    BASEH_TANHE = np.zeros(FM.PathX.NPATH)
    for i in range(FM.PathX.NPATH):
        BASEH_TANHE[i] = FM.LayerX.BASEH[FM.PathX.LAYINC[int(FM.PathX.NLAYIN[i]/2),i]]/1.0e3
    
    #Calculating the line-of-sight densities in each layer
    LOSDENS_LAYINC = FM.LayerX.AMOUNT[FM.PathX.LAYINC[:,:],:].T * FM.PathX.SCALE[:,:].T  #(NGAS,NPATH,NLAYIN)
        
    #Calculating the line-of-sight over the path
    LOSDENS_PATH = np.sum(LOSDENS_LAYINC,2) #(NGAS,NPATH)
    
    #Calculating the LOS density exactly at our tangent heights
    LOSDENS_PATHX = np.zeros((FM.AtmosphereX.NVMR,FM.MeasurementX.NGEOM))
    for i in range(FM.MeasurementX.NGEOM):
        
        #Find altitudes above and below the actual tangent height
        ibase = np.argmin(np.abs(BASEH_TANHE-FM.MeasurementX.TANHE[i]))
        base0 = BASEH_TANHE[ibase]
        if base0<=FM.MeasurementX.TANHE[i]:
            ibasel = ibase
            ibaseh = ibase + 1
        else:
            ibasel = ibase - 1
            ibaseh = ibase
            
        if ibaseh>FM.PathX.NPATH-1:
            LOSDENS_PATHX[:,i] = LOSDENS_PATH[:,ibasel]
        else:
            fhl = (FM.MeasurementX.TANHE[i]-BASEH_TANHE[ibasel])/(BASEH_TANHE[ibaseh]-BASEH_TANHE[ibasel])
            fhh = (BASEH_TANHE[ibaseh]-FM.MeasurementX.TANHE[i])/(BASEH_TANHE[ibaseh]-BASEH_TANHE[ibasel])
        
            LOSDENS_PATHX[:,i] = LOSDENS_PATH[:,ibasel] * (1.-fhl) + LOSDENS_PATH[:,ibaseh] * (1.-fhh)
    
    return LOSDENS_PATHX

###############################################################################################

def calc_los_tau(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval):
    
    """
    FUNCTION NAME : calc_los_tau()

    DESCRIPTION : Calculate the optical depth of all active gases at each altitude

    INPUTS : 

        NEMESIS python classes

    OPTIONAL INPUTS:
            
    OUTPUTS : 
    
        TAUTOT(NWAVE,NGAS,NPATH) :: Optical depth of each gas at each altitude level

    CALLING SEQUENCE:

        TAUTOT = calc_LOS_density(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval)

    MODIFICATION HISTORY : Juan Alday (09/08/2023)

    """
    
    #Initialise forward model class to perform calculations
    FM = ans.ForwardModel_0(runname='retrieval_LOS', Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables)

    #Initialising the X classes just for storing one location
    FM.MeasurementX = deepcopy(FM.Measurement)
    FM.AtmosphereX = deepcopy(FM.Atmosphere)
    FM.SurfaceX = deepcopy(FM.Surface)
    FM.ScatterX = deepcopy(FM.Scatter)
    FM.Scatter = None
    FM.StellarX = deepcopy(FM.Stellar)
    FM.Stellar = None
    FM.SpectroscopyX = deepcopy(FM.Spectroscopy)
    FM.Spectroscopy = None
    FM.LayerX = deepcopy(FM.Layer)
    FM.CIAX = deepcopy(FM.CIA)
    FM.CIA = None
    flagh2p = False
        
    #Calculating the path
    FM.calc_path_SO()
    
    BASEH_TANHE = np.zeros(FM.PathX.NPATH)
    for i in range(FM.PathX.NPATH):
        BASEH_TANHE[i] = FM.LayerX.BASEH[FM.PathX.LAYINC[int(FM.PathX.NLAYIN[i]/2),i]]/1.0e3
    
    #Calculating the cross sections at the required pressures and temperatures
    if FM.SpectroscopyX.ILBL == 2:
        k = FM.SpectroscopyX.calc_klbl(FM.LayerX.NLAY,FM.LayerX.PRESS/101325.,FM.LayerX.TEMP) #(NWAVE,NLAY,NGAS)
    elif FM.SpectroscopyX.ILBL == 1:
        k = FM.SpectroscopyX.calc_klbl_online(FM.LayerX.NLAY,FM.LayerX.PRESS/101325.,FM.LayerX.TEMP) #(NWAVE,NLAY,NGAS)

    #Calculating the vertical column optical depth
    TAUGAS = np.zeros((FM.SpectroscopyX.NWAVE,FM.SpectroscopyX.NGAS,FM.LayerX.NLAY))  #Vertical opacity of each gas in each layer
    for i in range(FM.SpectroscopyX.NGAS):
        IGAS = np.where( (FM.AtmosphereX.ID==FM.SpectroscopyX.ID[i]) & (FM.AtmosphereX.ISO==FM.SpectroscopyX.ISO[i]) )
        IGAS = IGAS[0]

        #Calculating vertical column density in each layer
        VLOSDENS = FM.LayerX.AMOUNT[:,IGAS].T * 1.0e-4   #cm-2

        #Calculating vertical opacity for each gas in each layer
        TAUGAS[:,i,:] = k[:,:,i] * VLOSDENS
        
    #Calculating the line-of-sight opacities
    TAUGAS_LAYINC = TAUGAS[:,:,FM.PathX.LAYINC[:,:]] * FM.PathX.SCALE[:,:]  #(NWAVE,NGAS,NLAYIN,NPATH)
        
    #Calculating the total opacity over the path
    TAUGAS_PATH = np.sum(TAUGAS_LAYINC,2) #(NWAVE,NGAS,NPATH)
    
    #Calculating the LOS opacity exactly at our tangent height
    if FM.PathX.NPATH!=2:
        sys.exit('error :: Number of calculated paths in calc_LOS_density() must be 2 (only one tangent height in Measurement)')
    TAUGAS_PATHX = np.zeros((FM.SpectroscopyX.NWAVE,FM.SpectroscopyX.NGAS))
    TAUGAS_PATHX[:,:] = TAUGAS_PATH[:,:,0] + (TAUGAS_PATH[:,:,1]-TAUGAS_PATH[:,:,0])/(BASEH_TANHE[1]-BASEH_TANHE[0]) * (FM.MeasurementX.TANHE[0,0]-BASEH_TANHE[0]) 
    
    return TAUGAS_PATHX

###############################################################################################

def calc_los_tau_all(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval):
    
    """
    FUNCTION NAME : calc_los_tau_all()

    DESCRIPTION : Calculate the optical depth of all active gases at each altitude

    INPUTS : 

        NEMESIS python classes

    OPTIONAL INPUTS:
            
    OUTPUTS : 
    
        TAUGAS_TANHE(NWAVE,NGAS,NGEOM) :: Optical depth of each gas at each tangent height

    CALLING SEQUENCE:

        TAUGAS_TANHE = calc_los_tau_all(Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval)

    MODIFICATION HISTORY : Juan Alday (09/08/2023)

    """
    
    #Initialise forward model class to perform calculations
    FM = ans.ForwardModel_0(runname='retrieval_LOS', Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables)

    #Initialising the X classes just for storing one location
    FM.MeasurementX = deepcopy(FM.Measurement)
    FM.AtmosphereX = deepcopy(FM.Atmosphere)
    FM.SurfaceX = deepcopy(FM.Surface)
    FM.ScatterX = deepcopy(FM.Scatter)
    FM.Scatter = None
    FM.StellarX = deepcopy(FM.Stellar)
    FM.Stellar = None
    FM.SpectroscopyX = deepcopy(FM.Spectroscopy)
    FM.Spectroscopy = None
    FM.LayerX = deepcopy(FM.Layer)
    FM.CIAX = deepcopy(FM.CIA)
    FM.CIA = None
    flagh2p = False
        
    #Calculating the path
    FM.calc_path_SO()
    
    BASEH_TANHE = np.zeros(FM.PathX.NPATH)
    for i in range(FM.PathX.NPATH):
        BASEH_TANHE[i] = FM.LayerX.BASEH[FM.PathX.LAYINC[int(FM.PathX.NLAYIN[i]/2),i]]/1.0e3
    
    #Calculating the cross sections at the required pressures and temperatures
    if FM.SpectroscopyX.ILBL == 2:
        k = FM.SpectroscopyX.calc_klbl(FM.LayerX.NLAY,FM.LayerX.PRESS/101325.,FM.LayerX.TEMP) #(NWAVE,NLAY,NGAS)
    elif FM.SpectroscopyX.ILBL == 1:
        k = FM.SpectroscopyX.calc_klbl_online(FM.LayerX.NLAY,FM.LayerX.PRESS/101325.,FM.LayerX.TEMP) #(NWAVE,NLAY,NGAS)

    #Calculating the vertical column optical depth
    TAUGAS = np.zeros((FM.SpectroscopyX.NWAVE,FM.SpectroscopyX.NGAS,FM.LayerX.NLAY))  #Vertical opacity of each gas in each layer
    for i in range(FM.SpectroscopyX.NGAS):
        IGAS = np.where( (FM.AtmosphereX.ID==FM.SpectroscopyX.ID[i]) & (FM.AtmosphereX.ISO==FM.SpectroscopyX.ISO[i]) )
        IGAS = IGAS[0]

        #Calculating vertical column density in each layer
        VLOSDENS = FM.LayerX.AMOUNT[:,IGAS].T * 1.0e-4   #cm-2

        #Calculating vertical opacity for each gas in each layer
        TAUGAS[:,i,:] = k[:,:,i] * VLOSDENS
        
    #Calculating the line-of-sight opacities
    TAUGAS_LAYINC = TAUGAS[:,:,FM.PathX.LAYINC[:,:]] * FM.PathX.SCALE[:,:]  #(NWAVE,NGAS,NLAYIN,NPATH)
        
    #Calculating the total opacity over the path
    TAUGAS_PATH = np.sum(TAUGAS_LAYINC,2) #(NWAVE,NGAS,NPATH)
    
    #Interpolating the optical depths to the correct tangent height
    TAUGAS_TANHE = np.zeros((Spectroscopy.NWAVE,Spectroscopy.NGAS,Measurement.NGEOM))

    for i in range(FM.MeasurementX.NGEOM):

        #Find altitudes above and below the actual tangent height
        ibase = np.argmin(np.abs(BASEH_TANHE-FM.MeasurementX.TANHE[i]))
        base0 = BASEH_TANHE[ibase]
        if base0<=FM.MeasurementX.TANHE[i]:
            ibasel = ibase
            ibaseh = ibase + 1
        else:
            ibasel = ibase - 1
            ibaseh = ibase

        if ibaseh>FM.PathX.NPATH-1:
            TAUGAS_TANHE[:,:,i] = TAUGAS_PATH[:,:,ibasel]
        else:
            fhl = (FM.MeasurementX.TANHE[i]-BASEH_TANHE[ibasel])/(BASEH_TANHE[ibaseh]-BASEH_TANHE[ibasel])
            fhh = (BASEH_TANHE[ibaseh]-FM.MeasurementX.TANHE[i])/(BASEH_TANHE[ibaseh]-BASEH_TANHE[ibasel])

            TAUGAS_TANHE[:,:,i] = TAUGAS_PATH[:,:,ibasel]*(1.-fhl) + TAUGAS_PATH[:,:,ibaseh]*(1.-fhh)
    
    return TAUGAS_TANHE

#Functions to add a fourier series to the baseline
#
# Define the Fourier series equation
def fourier_series(x, a0, *a):
    series = a0 / 2
    n = len(a) // 2
    for i in range(n):
        series += a[i] * np.cos((i + 1) * x) + a[n+i] * np.sin((i + 1) * x)
    return series

# Define the loss function to be minimized
def loss_function(params, x, y):
    a0 = params[0]
    a = params[1:]
    y_pred = fourier_series(x, a0, *a)
    return np.sum((y - y_pred) ** 2)
 
 
#Function defining the forward model
#   
import matplotlib.pyplot as plt

def forward_model(Measurement,WAVE,TAUGAS_PATH,ndegree,fourierdegree,**parameters):
    
    """
    FUNCTION NAME : forward_model()

    DESCRIPTION : Calculate the forward model with three gases as inputs.
                  In addition, we include a second order polynomial for the baseline.

    INPUTS : 

        Measurement :: NEMESIS measurement class including the measured spectrum and ILS
        TAUGAS_PATH(NWAVE,NGAS) :: Optical depth of each gas at each wavelength

    OPTIONAL INPUTS:
            
    OUTPUTS : 
    
        SPECMOD(NWAVE) :: Modelled spectrum 

    CALLING SEQUENCE:

        SPECMOD = forward_model(Measurement,TAUGAS_PATH,SC_GAS1,SC_GAS2,SC_GAS3,A0,A1,A2,A3)

    MODIFICATION HISTORY : Juan Alday (09/01/2024)

    """
    
    from scipy.optimize import minimize
    
    ngas = TAUGAS_PATH.shape[1]

    # Gas transmission
    scales = np.array([
        parameters[f"SC_GAS{i + 1}"]
        for i in range(ngas)
    ])

    TAUTOT = TAUGAS_PATH @ scales
    TRANS = np.exp(-TAUTOT)

    # Polynomial baseline: A0 + A1*WAVE + ... + Andegree*WAVE**ndegree
    coefficients = np.array([
        parameters[f"A{i}"]
        for i in range(ndegree + 1)
    ])

    TRANS_BASE = np.polynomial.polynomial.polyval(WAVE, coefficients)
    TRANS *= TRANS_BASE

    #Convolving the transmission spectrum with the ILS
    if Measurement.FWHM>0.0:
        TRANSCONV = ans.lblconv(len(WAVE),WAVE,TRANS,Measurement.NCONV[0],Measurement.VCONV[:,0],Measurement.ISHAPE,Measurement.FWHM)
    elif Measurement.FWHM<0.0:
        TRANSCONV = ans.lblconv_fil(len(WAVE),WAVE,TRANS,Measurement.NCONV[0],Measurement.VCONV[:,0],Measurement.NFIL,Measurement.VFIL,Measurement.AFIL)
    else:
        TRANSCONV = np.interp(Measurement.VCONV[:,0],WAVE,TRANS)

    return TRANSCONV


###############################################################################################################

def read_output(filename):
    
    """
    FUNCTION NAME : read_output()

    DESCRIPTION : Read the output file generated by retrieval_LOS_all()

    INPUTS : 

        filename :: Name of the HDF5 file

    OPTIONAL INPUTS:
            
    OUTPUTS : 
    
        NGEOM :: Number of geometries
        NCONV :: Number of convolution wavelengths
        VCONV(NCONV) :: Convolution wavelengths
        MEAS(NCONV,NGEOM) :: Measured spectra
        ERRMEAS(NCONV,NGEOM) :: Uncertainty in measured spectra
        SPECMOD(NCONV,NGEOM) :: Modelled spectra
        BASELINE(NCONV,NGEOM) :: Modelled baseline
        
        ID :: Gas ID of the fitted species
        ISO :: Isotope ID of the fitted species
        LOSDENS(NGEOM) :: Line-of-sight density at each tangent height (m-2)
        VMR_ret(NGEOM,NGAS) :: Retrieved VMR
        VMR_reterr(NGEOM,NGAS) :: Uncertainty in retrieved VMR

    CALLING SEQUENCE:

        NGEOM,NCONV,TANHE,VCONV,MEAS,ERRMEAS,SPECMOD_background,SPECMOD_trace,SPECMOD_sigma,BASELINE,ID,ISO,LOSDENS,VMR_ret,VMR_reterr = read_output(filename)

    MODIFICATION HISTORY : Juan Alday (09/01/2024)

    """
    
    dataset_names = [
        "NGEOM", "NCONV", "TANHE",
        "VCONV", "MEAS", "ERRMEAS",
        "SPECMOD_trace", "SPECMOD_sigma", "SPECMOD_background",
        "BASELINE", "ID", "ISO",
        "VMR_ret", "VMR_reterr",
    ]

    with h5py.File(filename + ".h5", "r") as f:
        data = {name: f[name][()] for name in dataset_names}

        # Store under your preferred name
        data["LOSDENS"] = f["LOSDENS_TOT"][()]

    
    return data


