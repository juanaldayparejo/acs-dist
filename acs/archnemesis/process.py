#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# process - Set of functions to create the input files for an ACS MIR retrieval
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

import acs
import archnemesis as ans
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os

###############################################################################################

def process_observation_mir(data_path,observation,windows,irows,
                            out_path_meas,out_path_atm,
                            latmos=True,iki=False,
                            ave_rows=True,
                            h_atm_res=1.,
                            ils_calfile=None,
                            hitran24_tables=True,):
    """
    FUNCTION NAME : process_observation_mir()

    DESCRIPTION : Function to extract the data from the ACS MIR file and create the archnemesis input files
    
    INPUTS : 

        data_path :: Directory where the ACS MIR files are stored
        observation :: Name of the input ACS MIR file
        windows(nwin) :: Name of the ACS MIR spectral windows to process
        irows(nrows) :: Index of the rows to use (from the slit edge on IKI version and from the usable ones in LATMOS version) 

    OPTIONAL INPUTS:
    
        latmos :: If True, it indicates the observation is one of the LATMOS FITS files
        iki :: If True, it indicates the observation is one of the IKI binary files
        ave_rows :: If True, the spectral from all the rows will be averaged
        ils_calfile :: If it exists, it is meant to be a list with the paths to the calibration file created with the ils fitting procedure of the package
                        the shape of the list is meant to be (nwin,nrow) if not averaging rows and (nwin,1) if averaging rows
        hitran24_tables :: If True, it indicates that the lbl tables are stored in the main HITRAN24 database

    OUTPUTS : 
 
        Output HDF5 file

    CALLING SEQUENCE:

        process_observation_acsmir_latmos(filename)

    MODIFICATION HISTORY : Juan Alday (02/08/2024)
    """

    #Doing some checks
    ####################################################################################

    if os.path.exists(data_path) is False:
        print("data_path = ",data_path)
        raise ValueError("error :: data_path does not exist")

    if os.path.exists(out_path_meas) is False:
        print("out_path_meas = ",out_path_meas)
        raise ValueError("error :: out_path_meas does not exist")

    if os.path.exists(out_path_atm) is False:
        print("out_path_atm = ",out_path_atm)
        raise ValueError("error :: out_path_atm does not exist")

    #Finding the file corresponding to the observation
    ####################################################################################

    curr = os.getcwd()
    os.chdir(data_path)
    filenames = glob.glob("*"+observation+"*")
    os.chdir(curr)

    if len(filenames) > 1:
        print(filenames)
        raise ValueError("error :: several files were found corresponding to the required observation")
    if len(filenames) == 0:
        raise ValueError("error :: no files were found corresponding to the required observation")

    filename = filenames[0]

    #Creating measurement class
    ####################################################################################

    nwin = len(windows)
    nrow = len(irows)

    hmin = 1000.
    hmax = -1000.

    trans_all = []
    transerr_all = []
    tanhe_all = []

    split_CO_iso = False
    split_H2O_iso = False
    split_CO2_iso = False

    Measurements_all = []
    Spectroscopy_all = []
    for iwin in range(nwin):

        Window = windows[iwin]

        #Reading the information about the spectral window
        Position = acs.mir.spectral_windows.process_orders_info[Window]["Position"]
        DifforSel = acs.mir.spectral_windows.process_orders_info[Window]["Diffor"]
        lowinX = acs.mir.spectral_windows.process_orders_info[Window]["WaveMin"]
        hiwinX = acs.mir.spectral_windows.process_orders_info[Window]["WaveMax"]
        MinTrans = acs.mir.spectral_windows.process_orders_info[Window]["MinTrans"]
        MinTanhe = acs.mir.spectral_windows.process_orders_info[Window]["MinTanhe"]
        MaxTanhe = acs.mir.spectral_windows.process_orders_info[Window]["MaxTanhe"]

        split_CO_isox = acs.mir.spectral_windows.process_orders_info[Window]["split_CO_iso"]
        split_H2O_isox = acs.mir.spectral_windows.process_orders_info[Window]["split_H2O_iso"]
        split_CO2_isox = acs.mir.spectral_windows.process_orders_info[Window]["split_CO2_iso"]   

        DELDG_apr = acs.mir.spectral_windows.process_orders_info[Window]["DELDG_apr"]
        FWHM_apr = acs.mir.spectral_windows.process_orders_info[Window]["FWHM_apr"]
        AMP1_apr = acs.mir.spectral_windows.process_orders_info[Window]["AMP1_apr"]
        AMP2_apr = acs.mir.spectral_windows.process_orders_info[Window]["AMP2_apr"]


        lat,lon,Ls,LST,waven,trans,transerr,tanhe = acs.mir.iki.files.extract_order(data_path+filename,DifforSel,irows)

        #Filtering the spectral range of the window
        iwave = np.where((waven>=lowinX) & (waven<=hiwinX))[0]
        waven = waven[iwave]
        trans = trans[iwave,:,:]
        transerr = transerr[iwave,:,:]

        #Updating the ILS and spectral calibration if required
        if ils_calfile is not None:
            raise ValueError("error :: ils update has not been implemented yet")
        else:
            deldg = DELDG_apr
            fwhm = FWHM_apr
            amp1 = AMP1_apr
            amp2 = AMP2_apr

        #Averaging rows if required
        if ave_rows is True:
            trans = np.mean(trans,axis=2, keepdims=True)
            transerr = np.sqrt(np.sum(transerr**2, axis=2, keepdims=True)) / nrow
            tanhe = np.mean(tanhe,axis=1, keepdims=True)
            nrow = 1

        #Calculating the mean transmission level in the window to assess the range of altitudes
        trans_level = np.mean(trans,axis=(0,2))
        tanhe_level = np.mean(tanhe,axis=1)

        #Calculating the tangent heights complying with the transmission level requirement
        for i in range(trans_level.shape[0]):
            
            if ((trans_level[i+1]>=MinTrans) & (trans_level[i]<MinTrans)):
                hmin_trans = tanhe_level[i+1]
                break

        #Calculating the range of altitudes to process
        hminx = np.max([hmin_trans,MinTanhe])
        hmaxx = MaxTanhe
        if hminx < hmin:
            hmin = hminx
        if hmaxx > hmax:
            hmax = hmaxx

        #Filtering the spectra
        iuse = np.where( (tanhe_level>=hminx) & (tanhe_level<=hmaxx) )[0]
        trans = trans[:,iuse,:]
        transerr = transerr[:,iuse,:]
        tanhe = tanhe[iuse,:]
        tanhe_level = tanhe_level[iuse]
        ngeom = len(iuse)

        #Calculating an optimised tangent height array in case there are too many altitude points (delH < 1 km)
        if (tanhe_level[ngeom-1]-tanhe_level[0])/ngeom>1.:
            tanhe_integer = np.round(tanhe_level,0)
            tanhe_good,itangood = np.unique(tanhe_integer, return_index=True)

            trans = trans[:,itangood,:]
            transerr = transerr[:,itangood,:]
            tanhe = tanhe[itangood,:]
            tanhe_level = tanhe_level[itangood]
            ngeom = len(itangood)


        #Creating measurement class 
        Measurements = []
        for irow in range(nrow):
            Measurement = acs.archnemesis.measurement.create_measurement_mir(lat,lon,tanhe[:,irow],waven,trans[:,:,irow],transerr[:,:,irow],deldg,fwhm,amp1,amp2)
            Measurements.append(Measurement)
        Measurements_all.append(Measurements)

        #Defining some common things for the Atmosphere
        if split_CO_isox is True:
            split_CO_iso = True

        if split_H2O_isox is True:
            split_H2O_iso = True

        if split_CO2_isox is True:
            split_CO2_iso = True

    Measurements_all = np.asarray(Measurements_all, dtype=object)


    #Creating Spectrocopy class
    ####################################################################################

    Spectroscopy_all = []
    for iwin in range(nwin):

        Window = windows[iwin]

        id_act = acs.mir.spectral_windows.process_orders_info[Window]["IDact"]
        iso_act = acs.mir.spectral_windows.process_orders_info[Window]["ISOact"]
        ngas = len(id_act)
        
        Spectroscopy = acs.archnemesis.spectroscopy.create_spectroscopy_lookup_hitran24(ngas,id_act,iso_act)

        Spectroscopy_all.append(Spectroscopy)


    #Creating common atmosphere for all windows
    ####################################################################################

    hmin_atm = np.ceil(hmin) - 2
    hmax_atm = np.floor(hmax).astype(int) + 2

    if hmin_atm < 0:
        hmin_atm = 0

    h_atm = np.arange(hmin,hmax+h_atm_res,h_atm_res)

    Atmosphere = acs.archnemesis.atmosphere.create_mcd_atmosphere_class(lat,lon,Ls,LST,h_atm)

    #Splitting isotopes if required
    if split_H2O_iso:
        Atmosphere = acs.archnemesis.atmosphere.split_H2O_isotopes(Atmosphere,dhratio=5.,o18ratio=1.,o17ratio=1.)
    if split_CO_iso:
        Atmosphere = acs.archnemesis.atmosphere.split_CO_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.)
    if split_CO2_iso:
        Atmosphere = acs.archnemesis.atmosphere.split_CO2_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.)

    #Writing some of the other classes
    #################################################################################################################
    
    Scatter = ans.Scatter_0(NDUST=1,ISPACE=0,ISCAT=0)
    Scatter.NWAVE = 2
    Scatter.WAVE = np.linspace(3000.,100000.,Scatter.NWAVE)
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
    Scatter.VEM = np.linspace(3000.,100000.,Surface.NEM)
    Surface.EMISSIVITY = np.ones(Surface.NEM)
    
    #Defining the retrieval parameters
    niter = 10
    Retrieval = ans.OptimalEstimation_0(PHILIMIT=0.1,NITER=niter,NCORES=1,IRET=0)

    #Saving the common classes into files
    ####################################################################################

    #Writing the HDF5 file
    os.chdir(out_path_atm)
    Atmosphere.write_hdf5(observation)
    Scatter.write_hdf5(observation)
    Layer.write_hdf5(observation)
    Stellar.write_hdf5(observation)
    Retrieval.write_input_hdf5(observation)
    os.chdir(curr)

    #Saving the Measurement dependent classes into files
    ####################################################################################

    for iwin in range(nwin):

        os.chdir(out_path_meas)

        #Changing folder to window level
        winpath = str(windows[iwin])
        if not os.path.exists(winpath):
            os.makedirs(winpath)
        os.chdir(winpath) 

        for irow in range(nrow):

            if ave_rows is True:
                rowpath = str("row-1")
            else:
                rowpath = "row"+str(int(irow))

            if not os.path.exists(rowpath):
                os.makedirs(rowpath)
            os.chdir(rowpath) 

            #Writing the HDF5 files
            Measurements_all[iwin,irow].write_hdf5(observation)
            Spectroscopy_all[iwin].write_hdf5(observation)

            os.chdir("../")

        os.chdir(curr)


    #Creating .apr files
    ####################################################################################

    for iwin in range(nwin):

        flag_temp = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Temp"]
        flag_temp_analytic = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Temp_Analytic"]  
        flag_press = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Press"]
        flag_hcorr = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Hcorr"]
        flag_gas = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Gas"]
        flag_baseline = acs.mir.spectral_windows.process_orders_info[Window]["Flag_Baseline"]

        if flag_press is True:
            htan = acs.mir.spectral_windows.process_orders_info[Window]["Htan"]
            ptan = None   #calculated automatically from the Atmosphere class at the selected tangent height
        else:
            htan = None
            ptan = None

        if flag_baseline is True:
            baseline_degree = acs.mir.spectral_windows.process_orders_info[Window]["Baseline_Degree"]
        else:
            baseline_degree = None

        if flag_gas is True:
            id_ret = acs.mir.spectral_windows.process_orders_info[Window]["IDret"]
            iso_ret = acs.mir.spectral_windows.process_orders_info[Window]["ISOret"]
            cont_err = np.ones(len(id_ret)) * 0.2
            cont_clen = np.ones(len(id_ret)) * 1.5
        else:
            id_ret = None
            iso_ret = None
            cont_err = None
            cont_clen = None

        if flag_temp is True or flag_temp_analytic is True:
            temp_err = 5.0
            temp_clen = 1.5
        else:
            temp_err = None
            temp_clen = None
            

        for irow in range(nrow):

            if ave_rows is True:
                rowpath = str("row-1")
            else:
                rowpath = "row"+str(int(irow))
            winpath = str(windows[iwin])

            os.chdir(out_path_meas+"/"+winpath+"/"+rowpath)

            acs.archnemesis.variables.create_apr_file(observation,
                                                    Atmosphere, 
                                                    Measurement,
                                                    retrieve_temp=flag_temp,
                                                    retrieve_press=flag_press,
                                                    retrieve_continuous_prof=flag_gas,
                                                    retrieve_baseline=flag_baseline,
                                                    flag_temp_analytic=flag_temp_analytic,
                                                    varID1=id_ret,varID2=iso_ret,
                                                    cont_clen=cont_clen,cont_err=cont_err,
                                                    temp_clen=temp_clen,temp_err=temp_err,
                                                    htan=htan,ptan=ptan,ptanerr=0.1,
                                                    baseline_degree=baseline_degree)

            os.chdir(curr)

