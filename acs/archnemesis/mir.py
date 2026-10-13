#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# archnemesis.mir - Functions to create the input classes for MIR retrievals
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
import sys,os

location_lbl_tables_hitran24 = "/exomars/retrievals/nemesis/spectroscopy/LBLtables/ACSMIR/HITRAN24/"

#########################################################################################################################################

def create_atmosphere(ACS_datadir,
                        Observation,
                        Windows,
                        latmos=True,
                        iki_geometry=True,
                        use_pt_nir=False,use_h2o_nir=False,extend_downward=False,
                        delh=1.,
                    ):
    """
    FUNCTION NAME : create_atmosphere()

    DESCRIPTION : Create the Atmosphere class for a MIR solar occultation observation

    INPUTS : 

        Observation :: Name of the Observation (e.g., "MIR_2B_ORB001941_N1_I_P1_00B00_F")
        Windows :: List of the spectral windows that will be included in the retrieval (see acs.mir.spectral_windows)

    OPTIONAL INPUTS:

        latmos :: If True, it will process the LATMOS FITS files
        iki_geometry :: If True, it will use the IKI geometry files
        use_pt_nir :: If True, it will update the Atmosphere using p-T profiles from ACS NIR (see acs.nir.products)
        use_h2o_nir :: If True, it will update the Atmosphere using H2O profiles from ACS NIR (see acs.nir.products)
        extend_downward :: If True, it will extend the p-T and H2O profiles from ACS NIR downward (using the MCD)
        delh :: Vertical sampling of the profiles (default = 1 km)
        
    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = create_atmosphere(Observation,Windows)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    if latmos is False:
        raise ValueError("The IKI files are not yet implemented for the Atmosphere class. Please use latmos=True")

    #Checking that some relevant parameters are shared between the selected spectral windows
    split_H2O_iso = False
    split_CO2_iso = False
    split_CO_iso = False
    position = 0
    ids = []
    isos = []
    
    ix = 0
    for Window in Windows:
        
        flag_h2o_iso = acs.mir.spectral_windows.process_orders_info[Window]["split_H2O_iso"]
        if ix > 0:
            if flag_h2o_iso != split_H2O_iso:
                raise ValueError("the spectral windows must all use the same split_H2O_iso flag")
        else:
            split_H2O_iso = flag_h2o_iso

        flag_co2_iso = acs.mir.spectral_windows.process_orders_info[Window]["split_CO2_iso"]
        if ix > 0:
            if flag_co2_iso != split_CO2_iso:
                raise ValueError("the spectral windows must all use the same split_CO2_iso flag")
        else:
            split_CO2_iso = flag_co2_iso

        flag_co_iso = acs.mir.spectral_windows.process_orders_info[Window]["split_CO_iso"]
        if ix > 0:
            if flag_co_iso != split_CO_iso:
                raise ValueError("the spectral windows must all use the same split_CO_iso flag")
        else:
            split_CO_iso = flag_co_iso


        positionx = acs.mir.spectral_windows.process_orders_info[Window]["Position"]
        if ix > 0:
            if positionx != position:
                raise ValueError("the spectral windows must all use the same position")
        else:
            position = positionx


        idx = acs.mir.spectral_windows.process_orders_info[Window]["IDact"]
        for j in range(len(idx)):
            ids.append(idx[j])

        isox = acs.mir.spectral_windows.process_orders_info[Window]["ISOact"]
        for j in range(len(isox)):
            isos.append(isox[j])
        
        ix += 1

    pairs = np.column_stack((ids, isos))
    unique_pairs = np.unique(pairs, axis=0)
    
    #Extracting information about occultation
    if latmos is True:
        orbit,num,ieflag,fpflag,position = acs.mir.latmos.files.get_info_from_name(Observation)

    #Searching record in mir database
    mir_database = acs.mir.database.read_database()
    record = acs.mir.database.search_records(mir_database,orbit,num)

    #Extracting the information about the position (generally)
    position_str = str(position)
    DifforSel = acs.mir.spectral_windows.process_atmosphere_info[position_str]["Diffor"]
    lowinX = acs.mir.spectral_windows.process_atmosphere_info[position_str]["WaveMin"]
    hiwinX = acs.mir.spectral_windows.process_atmosphere_info[position_str]["WaveMax"]
    MinTrans = acs.mir.spectral_windows.process_atmosphere_info[position_str]["MinTrans"]
    MinTanhe = acs.mir.spectral_windows.process_atmosphere_info[position_str]["MinTanhe"]
    MaxTanhe = acs.mir.spectral_windows.process_atmosphere_info[position_str]["MaxTanhe"]

    #Extracting the spectra for the selected order
    if latmos is True:
        acs_file = ACS_datadir+Observation+'.fits.fz'

    irows_template = [5,6,7,8,9,10,11]
    lat,lon,Ls,LST,VCONVX,MEASX,ERRMEASX,TANHEX = acs.mir.latmos.files.extract_order(acs_file,DifforSel,irows_template,refalt=180.,iki_geometry=iki_geometry)

    iwin = np.where( (VCONVX>=lowinX) & (VCONVX<=hiwinX) )[0]
    VCONVX = VCONVX[iwin]
    MEASX = MEASX[iwin,:,:]
    ERRMEASX = ERRMEASX[iwin,:,:]

    #Removing the rows that have no NaN values in the measurement
    irow_remove = []
    for i in range(len(irows_template)):
        if np.isnan(MEASX[:,:,i]).any():
            irow_remove.append(i)
    MEASX = np.delete(MEASX,irow_remove,axis=2)
    ERRMEASX = np.delete(ERRMEASX,irow_remove,axis=2)
    TANHEX = np.delete(TANHEX,irow_remove,axis=1)
    irows_template = np.delete(irows_template,irow_remove,axis=0)

    #Calculating average spectra from all rows
    TANHE_AVEX = np.mean(TANHEX,axis=1) 
    MEAS_AVEX = np.mean(MEASX,axis=2)
    ERRMEAS_AVEX = np.mean(ERRMEASX,axis=2)

    #Selecting the altitude range to consider in the observation
    meantrans = np.zeros(MEASX.shape[1])
    for j in range(MEASX.shape[1]):
        meantrans[j] = np.mean(MEAS_AVEX[:,j])
        
    iGeom = np.where( (TANHE_AVEX<MaxTanhe) & (TANHE_AVEX>0.) & (TANHE_AVEX>=MinTanhe) & (meantrans>=MinTrans) )[0]

    min_tanhe = np.floor(TANHE_AVEX[iGeom].min())
    max_tanhe = np.ceil(TANHE_AVEX[iGeom].max())

    #Creating altitude array
    h = np.arange(min_tanhe,max_tanhe+delh,delh)

    #Creating Atmosphere from MCD
    Atmosphere = acs.archnemesis.atmosphere.create_mcd_atmosphere_class(lat,lon,Ls,LST,h)

    #Updating atmosphere with NIR profile if required
    if use_pt_nir is True:

        #Searching record in mir database
        nir_database = acs.nir.database.read_database()
        nir_record = acs.nir.database.search_records_from_utc(nir_database, record["date"].iloc[0], record["time"].iloc[0], threshold_minutes=20.)
        
        #Reading NIR p-T files
        exists,h_nir,p_nir,t_nir,terr_nir = acs.nir.products.prepare_nir_inputs_pt_v1(nir_record["orbit"].iloc[0],nir_record["num"].iloc[0])

        #Updating Atmosphere based on NIR temperature and pressure profiles
        Atmosphere = acs.archnemesis.atmosphere.update_pt_atmosphere_with_profile(Atmosphere,h_nir,p_nir,t_nir,extend_downward=extend_downward,shift_nir=0.0,MakePlot=False)

    if use_h2o_nir is True:

        #Searching record in mir database
        nir_database = acs.nir.database.read_database()
        nir_record = acs.nir.database.search_records_from_utc(nir_database, record["date"].iloc[0], record["time"].iloc[0], threshold_minutes=20.)
        
        #Reading H2O profiles
        exists,h_h2o_nir,p_nir,t_nir,terr_nir,h2o_nir,h2oerr_nir = acs.nir.products.prepare_nir_inputs_h2o_v1(nir_record["orbit"].iloc[0],nir_record["num"].iloc[0])

        #Updating Atmosphere based on NIR H2O profiles
        Atmosphere = acs.archnemesis.atmosphere.update_h2o_atmosphere_with_profile(Atmosphere,h_h2o_nir,h2o_nir,extend_downward=extend_downward,shift_nir=0.0,MakePlot=True)


    #Splitting isotopes if required
    if split_H2O_iso is True:
        Atmosphere = acs.archnemesis.atmosphere.split_H2O_isotopes(Atmosphere,dhratio=5.,o18ratio=1.,o17ratio=1.)
    if split_CO2_iso is True:
        Atmosphere = acs.archnemesis.atmosphere.split_CO2_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.)
    if split_CO_iso is True:
        Atmosphere = acs.archnemesis.atmosphere.split_CO_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.)

    #Adding active species that are not present in the atmosphere by default
    for i in range(len(unique_pairs)):
        pair = unique_pairs[i,:]
        idx = pair[0] ; isox = pair[1]
        iindex = Atmosphere.locate_gas(idx,isox)

        if iindex is None:
            #Adding gas with 1 ppbv
            Atmosphere.add_gas(idx,isox,np.ones(Atmosphere.NP)*1.0e-9)
    
    return Atmosphere

#########################################################################################################################################

def create_measurement(ACS_datadir,
                           Observation,
                           Window,
                           irows,
                           Calfiles=None,
                           latmos=True,
                           iki_geometry=True,
                           ave_rows=False
                    ):
    """
    FUNCTION NAME : create_measurement()

    DESCRIPTION : Create the Measurement class for a MIR solar occultation observation

    INPUTS : 

        Observation :: Name of the Observation (e.g., "MIR_2B_ORB001941_N1_I_P1_00B00_F")
        Window :: The spectral window that will be included in the retrieval (see acs.mir.spectral_windows)
        irows :: List of detector rows to include in the retrieval (e.g., [1,2,3,4,5,6,7,8])

    OPTIONAL INPUTS:

        latmos :: If True, it will process the LATMOS FITS files
        iki_geometry :: If True, it will use the IKI geometry files
        ave_rows :: If True, it will average the selected rows to create a single measurement (default = False)
                       Otherwise, it will create a measurement for each row (default = False)
        Calfiles :: List of calibration files to use for the ILS (default = None, in which case the default ILS will be used)
                    Note that the default ILS will not be correct
        
    OUTPUTS : 
 
        Measurement :: archNEMESIS Measurement class

    CALLING SEQUENCE:

        Measurement = create_measurement(ACS_datadir,
                                             Observation,
                                             Window,
                                             latmos=True,
                                             iki_geometry=True,
                                             irows=[1,2,3,4,5,6,7,8],
                                             ave_rows=False
                                          )

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    if latmos is False:
        raise ValueError("The IKI files are not yet implemented for the Measurement class. Please use latmos=True")

    if Calfiles is not None:
        if ave_rows is True:
            #Checking that only 1 file exists
            check_calfile = np.array(Calfiles)
            if len(check_calfile.shape)!=1:
                raise ValueError("If ave_rows is True, only 1 calibration file should be provided. Please check the input Calfiles")
        else:
            #Checking that the number of calibration files is equal to the number of rows
            if len(Calfiles)!=len(irows):
                raise ValueError("The number of calibration files should be equal to the number of rows. Please check the input Calfiles and irows")

    #Extracting the information about the Window
    Position = acs.mir.spectral_windows.process_orders_info[Window]["Position"]
    DifforSel = acs.mir.spectral_windows.process_orders_info[Window]["Diffor"]
    lowinX = acs.mir.spectral_windows.process_orders_info[Window]["WaveMin"]
    hiwinX = acs.mir.spectral_windows.process_orders_info[Window]["WaveMax"]
    MinTrans = acs.mir.spectral_windows.process_orders_info[Window]["MinTrans"]
    MinTanhe = acs.mir.spectral_windows.process_orders_info[Window]["MinTanhe"]
    MaxTanhe = acs.mir.spectral_windows.process_orders_info[Window]["MaxTanhe"]

    DELDG_apr = acs.mir.spectral_windows.process_orders_info[Window]["DELDG_apr"]
    FWHM_apr = acs.mir.spectral_windows.process_orders_info[Window]["FWHM_apr"]
    AMP1_apr = acs.mir.spectral_windows.process_orders_info[Window]["AMP1_apr"]
    AMP2_apr = acs.mir.spectral_windows.process_orders_info[Window]["AMP2_apr"]

    #Extracting the spectra for the selected order
    if latmos is True:
        acs_file = ACS_datadir+Observation+'.fits.fz'

    lat,lon,Ls,LST,VCONVX,MEASX,ERRMEASX,TANHEX = acs.mir.latmos.files.extract_order(acs_file,DifforSel,irows,refalt=180.,iki_geometry=iki_geometry)

    iwin = np.where( (VCONVX>=lowinX) & (VCONVX<=hiwinX) )[0]
    VCONVX = VCONVX[iwin]
    MEASX = MEASX[iwin,:,:]
    ERRMEASX = ERRMEASX[iwin,:,:]

    #Removing the rows that have no NaN values in the measurement
    irow_remove = []
    for i in range(len(irows)):
        if np.isnan(MEASX[:,:,i]).any():
            print("Row "+str(irows[i])+" has NaN values in the measurement. It will be removed from the retrieval")
            irow_remove.append(i)
    MEASX = np.delete(MEASX,irow_remove,axis=2)
    ERRMEASX = np.delete(ERRMEASX,irow_remove,axis=2)
    TANHEX = np.delete(TANHEX,irow_remove,axis=1)
    irows = np.delete(irows,irow_remove,axis=0)

    #Averaging the rows if requested
    if ave_rows is True:
        nrows = 1
        #Averaging the rows
        TANHEX = np.mean(TANHEX,axis=1, keepdims=True) 
        MEASX = np.mean(MEASX,axis=2, keepdims=True)
        ERRMEASX = np.mean(ERRMEASX,axis=2, keepdims=True)
    else:
        nrows = len(irows)

    Measurements = []
    for irow in range(nrows):

        #Checking if the calibration file is provided
        if Calfiles is not None:
            if ave_rows is True:
                calfile = Calfiles[0]
            else:
                calfile = Calfiles[irow]

            WAVECONV,DELDG,FWHM,AMP1,AMP2 = acs.mir.latmos.ils.read_calfile(calfile)

        else:

            WAVECONV = VCONVX
            DELDG = DELDG_apr
            FWHM = FWHM_apr
            AMP1 = AMP1_apr
            AMP2 = AMP2_apr

        #Selecting tangent heights
        isel = np.where( (TANHEX[:,irow]>=MinTanhe) & (TANHEX[:,irow]<=MaxTanhe) )[0]

        #Creating the Measurement class for each row
        Measurement = create_measurement_general(
            lat,lon,
            TANHEX[isel,irow],
            WAVECONV,MEASX[:,isel,irow],ERRMEASX[:,isel,irow],
            DELDG,FWHM,AMP1,AMP2,
        )

        Measurements.append(Measurement)
    
    return Measurements


########################################################################################################################################

def create_measurement_general(
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

        Measurement = create_measurement_general(
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

#########################################################################################################################################

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

def create_spectroscopy(Window,
                        Datadir=location_lbl_tables_hitran24,online_k=True):

    """
    FUNCTION NAME : create_spectroscopy_lookup_hitran24()

    DESCRIPTION : Function to create the archNEMESIS Spectroscopy class based on the 
                    HITRAN24 format

    INPUTS : 

        Window :: The spectral window that will be included in the retrieval (see acs.mir.spectral_windows)

    OPTIONAL INPUTS:

        Datadir :: Path to the directory where the line-by-line tables are stored (default = location_lbl_tables_hitran24)
        online_k :: If True, it will use the online k-tables (default = True)
            
    OUTPUTS : 
 
        Spectroscopy :: Spectroscopy class

    CALLING SEQUENCE:

        Spectroscopy = create_spectroscopy(Window)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Extracting the information about the Window
    ids = acs.mir.spectral_windows.process_orders_info[Window]["IDact"]
    isos = acs.mir.spectral_windows.process_orders_info[Window]["ISOact"]
    ngas = len(ids)

    lbl_tables = location_lls_acsmir_hitran24(ngas,ids,isos,Datadir=Datadir,online_k=online_k)
    Spectroscopy = acs.archnemesis.spectroscopy.create_spectroscopy_lookup_class(lbl_tables)

    return Spectroscopy

###############################################################################################

def create_other_classes(Atmosphere):

    """
    FUNCTION NAME : create_other_classes()

    DESCRIPTION : Function to create the archNEMESIS classes that are standard for all solar
                    occultation observations

    INPUTS : 

        Atmosphere :: archNEMESIS Atmosphere class

    OPTIONAL INPUTS: None
            
    OUTPUTS : 
 
        

    CALLING SEQUENCE:

        Spectroscopy = create_spectroscopy(Window)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Defining Scatter class
    Scatter = ans.Scatter_0(NDUST=Atmosphere.NDUST,ISPACE=0,ISCAT=0)
    Scatter.NWAVE = 2
    Scatter.WAVE = np.linspace(1000.,5000.,Scatter.NWAVE)
    Scatter.KEXT = np.ones((Scatter.NWAVE,Scatter.NDUST))
    Scatter.SGLALB = np.ones((Scatter.NWAVE,Scatter.NDUST))
    
    #Defining Layer class
    Layer = ans.Layer_0(LAYTYP=5,LAYINT=1,LAYHT=0.0)
    Layer.NLAY = Atmosphere.NP - 1
    Layer.H_base = Atmosphere.H[0:Atmosphere.NP-1]
    
    #Defining Stellar class
    Stellar = ans.Stellar_0(SOLEXIST=False,ISPACE=0)
    
    #Defining Surface class
    Surface = ans.Surface_0(GALB=0.0,LOWBC=0)
    Surface.NEM = 2
    Scatter.VEM = np.linspace(1000.,5000.,Scatter.NWAVE)
    Surface.EMISSIVITY = np.ones(Surface.NEM)
    
    #Defining the retrieval parameters
    niter = 10
    Retrieval = ans.OptimalEstimation_0(PHILIMIT=0.1,NITER=niter,NCORES=1,IRET=0)

    return Scatter, Layer, Stellar, Surface, Retrieval

###############################################################################################

def create_variables(runname,Window,Atmosphere,Measurement,outdir="./"):

    """
    FUNCTION NAME : create_variables()

    DESCRIPTION : Create .apr file

    INPUTS : 

        runname :: Name of the file
        Window :: The spectral window that will be included in the retrieval (see acs.mir.spectral_windows)
        Atmosphere :: archNEMESIS Atmosphere class
        Measurement :: archNEMESIS Measurement class

    OPTIONAL INPUTS: None
            
    OUTPUTS : 
 
        archNEMESIS .apr file

    CALLING SEQUENCE:

        acs.archnemesis.mir.create_variables(Observation,Window,Atmosphere,Measurement)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Reading the information regarding the retrieval parameters for the selected window
    flag_press = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Press"]
    htan = None
    if flag_press:
        htan = acs.mir.spectral_windows.process_orders_info[Window]["Htan"]

    flag_hcorr = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Hcorr"]

    flag_temp = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Temp"]
    flag_temp_analytic = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Temp_Analytic"]
    temp_err = None
    temp_clen = None
    if flag_temp | flag_temp_analytic:
        temp_err = acs.mir.spectral_windows.process_orders_info[Window]["Temp_Err"]
        temp_clen = acs.mir.spectral_windows.process_orders_info[Window]["Temp_CLEN"]

    flag_gas = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Gas"]
    idret = None
    isoret = None
    cont_err = None
    cont_clen = None
    if flag_gas:
        idret = acs.mir.spectral_windows.process_orders_info[Window]["IDret"]
        isoret = acs.mir.spectral_windows.process_orders_info[Window]["ISOret"]
        cont_err = acs.mir.spectral_windows.process_orders_info[Window]["Gas_Err"]
        cont_clen = acs.mir.spectral_windows.process_orders_info[Window]["Gas_CLEN"]

    flag_baseline = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Baseline"]
    baseline_degree = None
    if flag_baseline:
        baseline_degree = acs.mir.spectral_windows.process_orders_info[Window]["Baseline_Degree"]

    #Creating .apr file
    curr = os.getcwd()
    os.chdir(outdir)
    acs.archnemesis.variables.create_apr_file(runname,
                                                Atmosphere, 
                                                Measurement,
                                                retrieve_temp=flag_temp,
                                                retrieve_press=flag_press,
                                                retrieve_continuous_prof=flag_gas,
                                                retrieve_baseline=flag_baseline,
                                                flag_temp_analytic=flag_temp_analytic,
                                                varID1=idret,varID2=isoret,
                                                cont_clen=cont_clen,cont_err=cont_err,
                                                temp_clen=temp_clen,temp_err=temp_err,
                                                htan=htan,ptan=None,ptanerr=0.1,
                                                baseline_degree=baseline_degree)
    os.chdir(curr)