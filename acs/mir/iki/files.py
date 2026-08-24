#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# iki.files - Functions to extract transmission spectra from the IKI binary files
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

acs_dir = Path(acs.__file__).resolve().parent

iki_geometry_dir = "/srv/workspace/data/tgo/acs/mir/iki/geometry/"

###############################################################################################

def extract_order(filename,ordersel,irows,iki_geometry_dir=iki_geometry_dir):

    """
    FUNCTION NAME : extract_order()

    DESCRIPTION : Function to extract a set of detector rows within a selected diffraction order 
                   from the IKI binary files

    INPUTS : 

        filename :: Name of the file
        ordersel :: Diffraction order
        irows(nrows) :: Index of row above the slit edge
        geomdir :: Base directory for the IKI geometry files

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        lat :: Latitude 
        lon :: Longitude
        Ls :: Solar longitude
        Loct :: Local time
        waven(nwave) :: Wavenumber array for the selected order (cm-1)
        trans(nwave,ngeom,nrows) :: Transmission array (3D) for the selected order
        transerr(nwave,ngeom,nrows) :: Uncertainty in transmission array
        tanhe(ngeom,nrows) :: Tangent altitude of each spectrum (km)

    CALLING SEQUENCE:

        lat,lon,Ls,Loct,Atmosphere,waven,trans,transerr,tanhe = process_observation_acsmir_latmos_order(filename,ordersel,irows,refalt=160.,iki_geometry=True)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    start = 'MIR_2A_'
    end = '_v7.'
    Observation = start+filename[filename.find(start)+len(start):filename.rfind(end)]

    #Reading the position number
    Position = int(Observation[25:27])

    #Reading the IKI binary file
    ##############################

    PointsCount,SeriesCount0A,RowsCount0A,ColumnsCount0A,OrdersCount0A,SeriesNumber0A,diffor_array,stripe_center,sun_ref,wave_array, \
         trans_array,transerr_array = readdata_acsmir(filename)

    diffor_array = diffor_array[0]
    stripe_center = stripe_center[0]

    ndiffor = len(diffor_array)

    iordersel = np.where( diffor_array == ordersel )[0]
    if len(iordersel) == 0:
        print("Observation = ",filename)
        print("Selected order = ",ordersel)
        print("Available orders = ",diffor_array)
        raise ValueError("error :: order not found")
    idifforsel = iordersel[0]

    #Calculating the slit size for each diffraction order
    ##########################################################

    npix = 512
    if RowsCount0A<npix:
        raise ValueError('error :: need to adapt the calculation of the slit size for RowsCount0A<512')
        
    #Reading the look-up table
    slitsizep = np.zeros(npix)
    for ipix in range(npix):
        slitsizep[ipix] = calc_slitsize_acsmir(Position,ipix)

    #Placing the look-up table in the middle of the frame
    slit_size_corr = np.zeros(RowsCount0A)
    row_corr = np.linspace(0.,RowsCount0A-1,RowsCount0A)

    ilo = int(int(RowsCount0A/2)-1-npix/2)
    ihi = int(int(RowsCount0A/2)-1+npix/2)

    slit_size_corr[ilo:ihi] = slitsizep[:]

    #Fitting a polynomial function to capture the row/slitsize dependence
    deg = 2
    pol = np.polyfit(row_corr[ilo:ihi],slit_size_corr[ilo:ihi],deg=2)

    slit_size_corr = np.zeros(RowsCount0A)
    for i in range(deg+1):
        slit_size_corr[:] = slit_size_corr[:] + pol[deg-i] * (row_corr[:])**i
        
    #Calculating the size of the slit at the row corresponding to the stripe centre
    slit_size = np.zeros(ndiffor)
    for idiffor in range(ndiffor):
        slit_size[idiffor] = slit_size_corr[stripe_center[idiffor]]


    #Calculating the position of the slit edge for each diffraction order
    ########################################################################
    
    #NEW APPROACH BASED ON THE LIMIT OF THE WAVENUMBER ARRAY

    #Calculating the row edge
    row_edge = np.zeros(ndiffor,dtype='int32')
    for i in range(ndiffor):
        
        jmax=30
        for j in range(jmax):
            row = stripe_center[i] - j
            if np.mean(wave_array[0,row,:])==0:
                row_edge[i] = row + 1
                break
                
    NRow = np.zeros(ndiffor,dtype='int32')
    NRow[:] = stripe_center-row_edge   


    #Assigning an altitude to the slit edge and calculating the row-to-row altitude variations
    ###########################################################################################

    #Reading the main observational parameters from the database
    orbit = int(Observation[10:16])
    ssflag = Observation[17:19]
    Nflag = int(Observation[18:19])
    IEflag = Observation[20:21]

    database = acs.mir.database.read_database()
    records = acs.mir.database.search_records(database,orbit,Nflag)

    mtp = records["mtp"].iloc[0]
    stp = records["stp"].iloc[0]

    #Reading the geometry files
    filename_geom1 = 'GEO_1_MIR_0B_'+Observation[7:29]+'_SWP.txt'    #Slit centre
    filename_geom2 = 'GEO_2_MIR_0B_'+Observation[7:29]+'_SWP.txt'    #Slit edge (upper part - not visible in images)
    filename_geom3 = 'GEO_3_MIR_0B_'+Observation[7:29]+'_SWP.txt'    #Slit edge (lower part)


    Geomdir = iki_geometry_dir+mtp+'/'+stp+"/"

    nacq1,lat_tgo1,lon_tgo1,lat_obs1,lon_obs1,lat_subsolar1,lon_subsolar1,alt_tgo1,\
        tanhe_areoid1,tanhe_ellipsoid1,alt_topo_areoid1,lst_obs1 = acs.mir.iki.read_geometry(Geomdir+filename_geom1)

    nacq2,lat_tgo2,lon_tgo2,lat_obs2,lon_obs2,lat_subsolar2,lon_subsolar2,alt_tgo2,\
        tanhe_areoid2,tanhe_ellipsoid2,alt_topo_areoid2,lst_obs2 = acs.mir.iki.read_geometry(Geomdir+filename_geom2)

    nacq3,lat_tgo3,lon_tgo3,lat_obs3,lon_obs3,lat_subsolar3,lon_subsolar3,alt_tgo3,\
        tanhe_areoid3,tanhe_ellipsoid3,alt_topo_areoid3,lst_obs3 = acs.mir.iki.read_geometry(Geomdir+filename_geom3)

    #Defining geometry for the slit edge, which is given by GEO3
    nacq_slit_edge = nacq3
    lat_tgo_slit_edge = lat_tgo3
    lon_tgo_slit_edge = lon_tgo3
    lat_obs_slit_edge = lat_obs3
    lon_obs_slit_edge = lon_obs3
    tanhe_areoid_slit_edge = tanhe_areoid3
    tanhe_ellipsoid_slit_edge = tanhe_ellipsoid3
    alt_topo_areoid_slit_edge = alt_topo_areoid3


    #Calculating the pixel-to-pixel altitude variation
    alt_res = np.zeros((ndiffor,nacq3))
    for idiffor in range(ndiffor):
        alt_res[idiffor,:] = (tanhe_areoid2[:]-tanhe_areoid3[:])/(slit_size[idiffor])
        
    
    #Assigning an altitude to each row in the diffraction order
    #############################################################
    
    #Defining geometry for the stripe centre
    nacq_stripe_centre = nacq3
    lat_tgo_stripe_centre = np.zeros((ndiffor,nacq3))
    lon_tgo_stripe_centre = np.zeros((ndiffor,nacq3))
    lat_obs_stripe_centre = np.zeros((ndiffor,nacq3))
    lon_obs_stripe_centre = np.zeros((ndiffor,nacq3))
    tanhe_areoid_stripe_centre = np.zeros((ndiffor,nacq3))
    tanhe_ellipsoid_stripe_centre = np.zeros((ndiffor,nacq3))
    alt_topo_areoid_stripe_centre = np.zeros((ndiffor,nacq3))
    for idiffor in range(ndiffor):
        lat_tgo_stripe_centre[idiffor,:] = lat_tgo3[:] + (stripe_center[idiffor]-row_edge[idiffor]) * (lat_tgo2[:]-lat_tgo3[:])/(slit_size[idiffor])
        lon_tgo_stripe_centre[idiffor,:] = lon_tgo3[:] + (stripe_center[idiffor]-row_edge[idiffor]) * (lon_tgo2[:]-lon_tgo3[:])/(slit_size[idiffor])
        lat_obs_stripe_centre[idiffor,:] = lat_obs3[:] + (stripe_center[idiffor]-row_edge[idiffor]) * (lat_obs2[:]-lat_obs3[:])/(slit_size[idiffor])
        lon_obs_stripe_centre[idiffor,:] = lon_obs3[:] + (stripe_center[idiffor]-row_edge[idiffor]) * (lon_obs2[:]-lon_obs3[:])/(slit_size[idiffor])
        tanhe_areoid_stripe_centre[idiffor,:] = tanhe_areoid3[:] + (stripe_center[idiffor]-row_edge[idiffor]) * (tanhe_areoid2[:]-tanhe_areoid3[:])/(slit_size[idiffor])
        tanhe_ellipsoid_stripe_centre[idiffor,:] = tanhe_ellipsoid3[:] + (stripe_center[idiffor]-row_edge[idiffor]) * (tanhe_ellipsoid2[:]-tanhe_ellipsoid3[:])/(slit_size[idiffor])


    #Defining the geometry for the average between slit edge and stripe centre
    ##############################################################################

    nacq_v1 = nacq_stripe_centre
    lat_tgo_v1 = (lat_tgo_slit_edge+lat_tgo_stripe_centre)/2.
    lon_tgo_v1 = (lon_tgo_slit_edge+lon_tgo_stripe_centre)/2.
    lat_obs_v1 = (lat_obs_slit_edge+lat_obs_stripe_centre)/2.
    lon_obs_v1 = (lon_obs_slit_edge+lon_obs_stripe_centre)/2.
    tanhe_areoid_v1 = (tanhe_areoid_slit_edge+tanhe_areoid_stripe_centre)/2.
    tanhe_ellipsoid_v1 = (tanhe_ellipsoid_slit_edge+tanhe_ellipsoid_stripe_centre)/2.
    alt_topo_areoid_v1 = (alt_topo_areoid_slit_edge+alt_topo_areoid_stripe_centre)/2.


    NDiffor = OrdersCount0A
    NCONV = ColumnsCount0A
    NRow = np.zeros(NDiffor,dtype='int32')
    NRow[:] = stripe_center-row_edge
    NRow = NRow[idifforsel]
    NAcq = SeriesCount0A

    #Determining the geometry for the acquisition included in the observation file (SeriesNumber0A) for the Average Row
    #####################################################################################################################

    Tanhe_Areoid_AveRow1 = np.zeros(NAcq)
    Tanhe_Ellipsoid_AveRow1 = np.zeros(NAcq)
    Latobs_AveRow1 = np.zeros(NAcq)
    Lonobs_AveRow1 = np.zeros(NAcq)
    for i in range(NAcq):
        Tanhe_Areoid_AveRow1[i] = tanhe_areoid_v1[idifforsel,SeriesNumber0A[i]]
        Tanhe_Ellipsoid_AveRow1[i] = tanhe_ellipsoid_v1[idifforsel,SeriesNumber0A[i]]
        Latobs_AveRow1[i] = lat_obs_v1[idifforsel,SeriesNumber0A[i]]
        Lonobs_AveRow1[i] = lon_obs_v1[idifforsel,SeriesNumber0A[i]]

    #Extracting selected spectra
    ################################################

    if np.max(irows) > NRow - 1:
        print("selected rows = ",irows)
        print("maximum number of rows for selected diffraction order = ",NRow)
        raise ValueError("error")
    nrowsel = len(irows)

    Tanhe_Areoid1 = np.zeros((NAcq,nrowsel))
    Tanhe_Ellipsoid1 = np.zeros((NAcq,nrowsel))
    Latobs1 = np.zeros((NAcq,nrowsel))
    Lonobs1 = np.zeros((NAcq,nrowsel)) 

    VCONV_AveRow1 = np.zeros((NCONV,NAcq))
    MEAS_AveRow1 = np.zeros((NCONV,NAcq))
    ERRMEAS_AveRow1 = np.zeros((NCONV,NAcq))

    VCONV1 = np.zeros((NCONV,NAcq,nrowsel))
    MEAS1 = np.zeros((NCONV,NAcq,nrowsel))
    ERRMEAS1 = np.zeros((NCONV,NAcq,nrowsel))


    irowx = 0
    for IRow in irows:

        #Getting the geometry for each row
        for i in range(NAcq):
            Tanhe_Areoid1[i,irowx] = tanhe_areoid_slit_edge[SeriesNumber0A[i]] + IRow * (tanhe_areoid_stripe_centre[idifforsel,SeriesNumber0A[i]]-tanhe_areoid_slit_edge[SeriesNumber0A[i]])/(stripe_center[idifforsel]-row_edge[idifforsel]) 
            Tanhe_Ellipsoid1[i,irowx] = tanhe_ellipsoid_slit_edge[SeriesNumber0A[i]] + IRow * (tanhe_ellipsoid_stripe_centre[idifforsel,SeriesNumber0A[i]]-tanhe_ellipsoid_slit_edge[SeriesNumber0A[i]])/(stripe_center[idifforsel]-row_edge[idifforsel])
            Latobs1[i,irowx] = lat_obs_slit_edge[SeriesNumber0A[i]] + IRow * (lat_obs_stripe_centre[idifforsel,SeriesNumber0A[i]]-lat_obs_slit_edge[SeriesNumber0A[i]])/(stripe_center[idifforsel]-row_edge[idifforsel])
            Lonobs1[i,irowx] = lon_obs_slit_edge[SeriesNumber0A[i]] + IRow * (lon_obs_stripe_centre[idifforsel,SeriesNumber0A[i]]-lon_obs_slit_edge[SeriesNumber0A[i]])/(stripe_center[idifforsel]-row_edge[idifforsel])

            VCONV_AveRow1[0:NCONV,i] = VCONV_AveRow1[0:NCONV,i] + wave_array[0,row_edge[idifforsel]+IRow,0:NCONV]
            MEAS_AveRow1[0:NCONV,i] = MEAS_AveRow1[0:NCONV,i] + trans_array[i,0,row_edge[idifforsel]+IRow,0:NCONV]
            ERRMEAS_AveRow1[0:NCONV,i] = ERRMEAS_AveRow1[0:NCONV,i] + (transerr_array[i,0,row_edge[idifforsel]+IRow,0:NCONV])**2.

            VCONV1[0:NCONV,i,irowx] = wave_array[0,row_edge[idifforsel]+IRow,0:NCONV]
            MEAS1[0:NCONV,i,irowx] = trans_array[i,0,row_edge[idifforsel]+IRow,0:NCONV]
            ERRMEAS1[0:NCONV,i,irowx] = transerr_array[i,0,row_edge[idifforsel]+IRow,0:NCONV]

        irowx += 1

    VCONV_AveRow1[0:NCONV,:] = VCONV_AveRow1[0:NCONV,:] / nrowsel
    MEAS_AveRow1[0:NCONV,:] = MEAS_AveRow1[0:NCONV,:] / nrowsel
    ERRMEAS_AveRow1[0:NCONV,:] = np.sqrt(ERRMEAS_AveRow1[0:NCONV,:] / nrowsel)


    #Sorting the acquisitions to be in ascending order from bottom to top of atmosphere
    #####################################################################################

    isort = np.argsort(Tanhe_Areoid_AveRow1)

    Tanhe_Areoid_AveRow = Tanhe_Areoid_AveRow1[isort]
    Tanhe_Ellipsoid_AveRow = Tanhe_Ellipsoid_AveRow1[isort]
    Latobs_AveRow = Latobs_AveRow1[isort]
    Lonobs_AveRow = Lonobs_AveRow1[isort]    

    Tanhe_Areoid = np.zeros((NAcq,nrowsel))
    Tanhe_Ellipsoid = np.zeros((NAcq,nrowsel))
    Latobs = np.zeros((NAcq,nrowsel))
    Lonobs = np.zeros((NAcq,nrowsel))

    VCONV_AveRow = np.zeros((NCONV,NAcq))
    MEAS_AveRow = np.zeros((NCONV,NAcq))
    ERRMEAS_AveRow = np.zeros((NCONV,NAcq))

    VCONV = np.zeros((NCONV,NAcq,nrowsel))
    MEAS = np.zeros((NCONV,NAcq,nrowsel))
    ERRMEAS = np.zeros((NCONV,NAcq,nrowsel))

    for IRow in range(nrowsel):
        for i in range(NAcq):
            Tanhe_Areoid[i,IRow] = Tanhe_Areoid1[isort[i],IRow]
            Tanhe_Ellipsoid[i,IRow] = Tanhe_Ellipsoid1[isort[i],IRow]
            Latobs[i,IRow] = Latobs1[isort[i],IRow]
            Lonobs[i,IRow] = Lonobs1[isort[i],IRow]          

            VCONV_AveRow[0:NCONV,i] = VCONV_AveRow1[0:NCONV,isort[i]]
            MEAS_AveRow[0:NCONV,i] = MEAS_AveRow1[0:NCONV,isort[i]]
            ERRMEAS_AveRow[0:NCONV,i] = ERRMEAS_AveRow1[0:NCONV,isort[i]]

            VCONV[0:NCONV,i,:] = VCONV1[0:NCONV,isort[i],:]
            MEAS[0:NCONV,i,:] = MEAS1[0:NCONV,isort[i],:]
            ERRMEAS[0:NCONV,i,:] = ERRMEAS1[0:NCONV,isort[i],:]

    Tanhe_Areoid_AveRowX = np.mean(Tanhe_Areoid,axis=1)
    Tanhe_Ellipsoid_AveRowX = np.mean(Tanhe_Ellipsoid,axis=1)
    Latobs_AveRowX = np.mean(Latobs,axis=1)
    Lonobs_AveRowX = np.mean(Lonobs,axis=1)

    lat = database["Latitude"]
    lon = database["Longitude"]
    Ls = database["Ls"]
    Loct = database["LST"]

    return lat,lon,Ls,Loct,VCONV_AveRow[:,0],MEAS,ERRMEAS,Tanhe_Areoid


###############################################################################################

def readdata_acsmir(filename):

    """
    FUNCTION NAME : readdata_acsmir()

    DESCRIPTION : This function reads the 2A data files from an ACS MIR observation in a format
                  compatible with PDS4

    INPUTS : 

        filename :: Name of the file

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        PointsCount :: Number of secondary grating positions used in the observation
        SeriesCount0A :: Number of acquisitions
        RowsCount0A :: Number of rows in detector frame (spatial direction)
        ColumnsCount0A :: Number of columns in detector frame (spectral direction)
        OrdersCount0A :: Number of diffraction orders
        SeriesNumber0a(SeriesCount0A) :: Frame number associated in the geometry files
        diffor_array(PointsCount,OrdersCount0A) :: Diffraction orders
        stripe_center(PointsCount,OrdersCount0A) :: Row number of center of diffraction order
        sun_ref(PointsCount,RowsCount0A,ColumnsCount0A) :: Solar reference 
        wave_array(PointsCount,RowsCount0A, ColumnsCount0A) :: Wavenumber of each pixel (cm-1)
        trans_array(SeriesCount0A, PointsCount, RowsCount0A, ColumnsCount0A) :: Transmission
        transerr_array(SeriesCount0A, PointsCount, RowsCount0A, ColumnsCount0A) :: Uncertainty in transmission

    CALLING SEQUENCE:

        PointsCount,SeriesCount0A,RowsCount0A,ColumnsCount0A,OrdersCount0A,SeriesNumber0A,diffor_array,stripe_center,sun_ref,wave_array, \
         trans_array,transerr_array = readdata_acsmir(filename)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    f = open(filename,'r')

    #Reading version
    version = np.fromfile(f,dtype='int32',count=1)

    #Reading command
    s_command = np.fromfile(f,dtype='int32',count=4)

    #Reading config
    s_config = np.fromfile(f,dtype='int32',count=4)

    #Reading control
    s_control = np.fromfile(f,dtype='int32',count=3)

    #Reading timing
    s_timing = np.fromfile(f,dtype='int32',count=2)

    #Reading detector
    s1_detector = np.fromfile(f,dtype='int32',count=1)
    s2_detector = np.fromfile(f,dtype='float64',count=3)
    s3_detector = np.fromfile(f,dtype='int32',count=2)
    s4_detector = np.fromfile(f,dtype='float64',count=1)

    #Reading frame
    s_begin1 = np.zeros([5],dtype='int32')
    s_end1 = np.zeros([5],dtype='int32')
    for i in range(5):
        s_temp = np.fromfile(f,dtype='int32',count=1)
        s_begin1[i] = s_temp
        s_temp = np.fromfile(f,dtype='int32',count=1)
        s_end1[i] = s_temp

    #Reading position
    s_position =  np.zeros([5,10],dtype='int32')
    for i in range(10):
        s_temp = np.fromfile(f,dtype='int32',count=5)
        s_position[:,i] = s_temp[:]

    #Reading cooling
    s_cooling = np.fromfile(f,dtype='int32',count=1)


    #Reading vByte
    s_vByte = np.fromfile(f,dtype='int32',count=1)

    #Reading number of series
    SeriesCount0A1 = np.fromfile(f,dtype='int32',count=1)
    SeriesCount0A = SeriesCount0A1[0]

    #Reading number of points (position)
    PointsCount1 = np.fromfile(f,dtype='int32',count=1)
    PointsCount = PointsCount1[0]

    #Reading number of rows (spatial dimension)
    RowsCount0A1 = np.fromfile(f,dtype='int32',count=1)
    OrdersCount0A1 = np.fromfile(f,dtype='int32',count=1)
    RowsCount0A = RowsCount0A1[0]
    OrdersCount0A = OrdersCount0A1[0]

    #Reading number of columns (spectral dimension)
    s_temp = np.fromfile(f,dtype='int32',count=2)
    pix_left = s_temp[0]
    ColumnsCount0A = s_temp[1]

    #Reading the difference between board and local time
    s_diff = np.fromfile(f,dtype='float64',count=1)

    trans_array = np.zeros([SeriesCount0A,PointsCount,RowsCount0A,ColumnsCount0A])
    transerr_array = np.zeros([SeriesCount0A,PointsCount,RowsCount0A,ColumnsCount0A])
    wave_array = np.zeros([PointsCount,RowsCount0A,ColumnsCount0A])
    diffor_array = np.zeros([PointsCount,OrdersCount0A],dtype='int32')
    stripe_center = np.zeros([PointsCount,OrdersCount0A],dtype='int32')
    sun_ref = np.zeros([PointsCount,RowsCount0A,ColumnsCount0A])
    SeriesNumber0A = np.zeros([SeriesCount0A],dtype='int32')

    for i in range(PointsCount):
        for j in range(RowsCount0A):
            s_temp = np.fromfile(f,dtype='int32',count=ColumnsCount0A)
            sun_ref[i,j,:] = s_temp[:]

    for i in range(SeriesCount0A):
        for j in range(PointsCount):
            #Reading local time
            s_loctime = np.fromfile(f,dtype='float64',count=1)

            #Reading frame number
            s_framenumber = np.fromfile(f,dtype='int32',count=2)
            SeriesNumber0A[i] = int(s_framenumber[1])

            #Reading config
            s_config_frame = np.fromfile(f,dtype='int32',count=4)

            #Reading control
            s_control_frame = np.fromfile(f,dtype='int32',count=3)

            #Reading detector
            s_video_mode = np.fromfile(f,dtype='int32',count=1)
            s_gain = np.fromfile(f,dtype='int32',count=1)
            s_offset = np.fromfile(f,dtype='float64',count=1)
            s_exposition = np.fromfile(f,dtype='float64',count=1)
            s_accum = np.fromfile(f,dtype='int32',count=1)
            s_Temp_D = np.fromfile(f,dtype='float64',count=1)

            #Reading position
            s_position = np.fromfile(f,dtype='int32',count=6)


    for i in range(SeriesCount0A):
        for j in range(PointsCount):
                for k in range(RowsCount0A):
                    s_temp = np.fromfile(f,dtype='float64',count=ColumnsCount0A)
                    trans_array[i,j,k,:] = s_temp[:]
                    s_temp = np.fromfile(f,dtype='float64',count=ColumnsCount0A)
                    transerr_array[i,j,k,:] = s_temp[:]

    for j in range(PointsCount):
        for k in range(RowsCount0A):
            s_temp = np.fromfile(f,dtype='float64',count=ColumnsCount0A)
            wave_array[j,k,:] = s_temp[:]

    for j in range(PointsCount):
        for k in range(OrdersCount0A):
            s_temp = np.fromfile(f,dtype='int32',count=1)
            diffor_array[j,k] = int(s_temp[0])
            s_temp = np.fromfile(f,dtype='int32',count=1)
            stripe_center[j,k] = int(s_temp[0])

    f.close()

    return PointsCount,SeriesCount0A,RowsCount0A,ColumnsCount0A,OrdersCount0A,SeriesNumber0A,diffor_array,stripe_center,\
           sun_ref,wave_array,trans_array,transerr_array

###############################################################################################

def read_geometry(filename,correct_boresight=True):

    """
    FUNCTION NAME : read_geometry()

    DESCRIPTION : This function reads the IKI geometry file for an ACS MIR observation given
                  by the format in the MIR_geometry_pr_tc/

    INPUTS : 

        filename :: Name of the file

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
          
        nacq :: Number of acquisitions made in the observation
        lat_tgo(nacq) :: Latitude of the sub-TGO point at each acquisition (degrees)
        lon_tgo(nacq) :: Longitude of the sub-TGO point at each acquisition (degrees)
        lat_obs(nacq) :: Latitude of the tangent point at each acquisition (degrees)
        lon_obs(nacq) :: Longitude of the tangent point at each acquisition (degrees)
        lat_subsolar(nacq) :: Latitude of the subsolar point at each acquisition (degrees)
        lon_subsolar(nacq) :: Longitude of the subsolar point at each acquisition (degrees)
        alt_tgo(nacq) :: Altitude of TGO above the martian ellipsoid (km)
        tanhe_aeroid(nacq) :: Tangent height of each acquisition above the Martian aeroid (km)
        tanhe_ellipsoid(nacq) :: Tangent height of each acquisition above the Martian ellipsoid (km)
        alt_topo_aeroid(nacq) :: Altitude of the local surface above the Martian aeroid (km)

    CALLING SEQUENCE:

        nacq,lat_tgo,lon_tgo,lat_obs,lon_obs,lat_subsolar,lon_subsolar,alt_tgo,tanhe_areoid,\
         tanhe_ellipsoid,alt_topo_areoid = read_geometry_pr_tc(filename)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    #Getting number of lines in file
    with open(filename, "r", encoding="utf-8") as file:
        line_count = sum(1 for _ in file)

    nacq = line_count - 3

    f = open(filename,'r')
    dummy1 = f.readline().split()
    dummy2 = f.readline().split()
    dummy3 = f.readline().split()

    lat_tgo = np.zeros([nacq])
    lon_tgo = np.zeros([nacq])
    lat_obs = np.zeros([nacq])
    lon_obs = np.zeros([nacq])
    lst_obs = np.zeros([nacq])
    lat_subsolar = np.zeros([nacq])
    lon_subsolar = np.zeros([nacq])
    alt_tgo = np.zeros([nacq])
    alt_topo_aeroid = np.zeros([nacq])
    r_aeroid = np.zeros([nacq])
    r_ellipsoid = np.zeros([nacq])
    tanhe_aeroid = np.zeros([nacq])
    tanhe_ellipsoid = np.zeros([nacq])
    occ_distance = np.zeros([nacq])

    for i in range(nacq):
        tmp = f.readline().split()
        lat_tgo[i] = tmp[5]
        lon_tgo[i] = tmp[6]
        lat_obs[i] = tmp[12]
        lon_obs[i] = tmp[13]
        lat_subsolar[i] = tmp[18]
        lon_subsolar[i] = tmp[19]
        alt_tgo[i] = tmp[20]
        alt_topo_aeroid[i] = tmp[15]
        occ_distance[i] = tmp[16]
        #r_aeroid[i] = tmp[23]
        #r_ellipsoid[i] = tmp[24]
        tanhe_aeroid[i] = tmp[20]
        tanhe_ellipsoid[i] = tmp[21]
        
        Loct1 = tmp[8]
        hour = float(Loct1[0:2])
        minute = float(Loct1[3:5])
        second = float(Loct1[6:8])
        lst_obs[i] = hour + minute/60. + second/3600.
        
    f.close()
    
    if correct_boresight is True:
        
        #Boresight angle correction
        boresight_correction = 0.045  #angle that the boresight is off in degrees
        
        delta_z = occ_distance * np.tan( np.radians(boresight_correction) )  #km
        tanhe_aeroid += delta_z
        tanhe_ellipsoid += delta_z
        
    return nacq,lat_tgo,lon_tgo,lat_obs,lon_obs,lat_subsolar,lon_subsolar,alt_tgo,tanhe_aeroid,tanhe_ellipsoid,alt_topo_aeroid,lst_obs


###############################################################################################

def calc_slitedge_acsmir(filename,MakePlot=False):

    """
    FUNCTION NAME : calc_slitedge_acsmir()

    DESCRIPTION : Calculate the row of the slit edge of ACS MIR observations

    INPUTS : 

        filename :: Name of the file

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        ndiffor :: Number of diffraction orders
        row_edge(ndiffor) :: Row number associated with the slit edge of each diffraction order


    CALLING SEQUENCE:

        ndiffor,row_edge = calc_slitedge_acsmir(filename)

    MODIFICATION HISTORY : Juan Alday (15/07/2021)

    """

    from scipy.signal import find_peaks_cwt
 
    PointsCount,SeriesCount0A,RowsCount0A,ColumnsCount0A,OrdersCount0A,SeriesNumber0A,diffor_array,stripe_center,sun_ref,wave_array, \
         trans_array,transerr_array = readdata_acsmir(filename)

    diffor_array = diffor_array[0]
    stripe_center = stripe_center[0]
    ndiffor = len(diffor_array)

    icol = 250   #Using a column at the centre of the detector frame
    irowmid = 4   #Row above the slit edge to include in the plots

    #Calculating the mean transmission of all series to check carefully the results
    meantrans = np.zeros(SeriesCount0A)
    for i in range(SeriesCount0A):
        meantrans[i] = np.mean(trans_array[i,0,stripe_center[ndiffor-1],:])

    #tanhe0,iseries = find_nearest(meantrans,0.6)    #Finding the series where the mean transmission is 0.6
    diff = np.abs(meantrans - 0.6)
    iseries = np.argmin(diff)

    #Calculating where the slit edge is for each stripe     
    vertical_cut = np.zeros([RowsCount0A])
    vertical_cut = sun_ref[0,:,icol]
    grad = np.gradient(vertical_cut)
    peaks2 = find_peaks_cwt(grad,np.arange(1,15))   #Finding the maxima in the gradient
    peaks = np.zeros([ndiffor],dtype='int32')
    for i in range(ndiffor):
        ipeak1 = np.where(peaks2<stripe_center[i])
        ipeak = ipeak1[0]
        peaks[i] = peaks2[ipeak[len(ipeak)-1]]


    if MakePlot==True:
        fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,10))
        ax1.imshow(sun_ref[0,:,:],origin='lower',interpolation='nearest', aspect='auto')    
 
        ax2.plot(vertical_cut,range(RowsCount0A),c='green')
        ax2.scatter(vertical_cut[stripe_center],stripe_center,c='blue')
        ax2.scatter(vertical_cut[peaks],peaks,c='orange')
        #ax2.scatter(vertical_cut[peaks+irowmid],peaks+irowmid,c='red')
        ax2.plot(grad,range(RowsCount0A),c='orange')

        ax3.imshow(sun_ref[0,:,:],origin='lower',interpolation='nearest', aspect='auto')
        ax3.axvline(icol,color='green')
        for i in range(ndiffor):
            ax3.axhline(stripe_center[i],color='blue')
            ax3.axhline(peaks[i],color='orange')
            #ax3.axhline(peaks[i]+irowmid,color='red')

        plt.tight_layout()
        plt.show()

        for i in range(ndiffor):
            plt.plot(wave_array[0,stripe_center[i],:],trans_array[iseries,0,stripe_center[i],:],c='blue')
            plt.plot(wave_array[0,peaks[i],:],trans_array[iseries,0,peaks[i],:],c='orange')
            plt.plot(wave_array[0,peaks[i]+irowmid,:],trans_array[iseries,0,peaks[i]+irowmid,:],c='red')

        plt.show()

    row_edge = peaks

    return ndiffor,row_edge

###############################################################################################

def calc_slitsize_acsmir(position,row):
    
    """
        FUNCTION NAME : calc_slitsize_mir()
        
        DESCRIPTION : Calculate the slit width for a given secondary grating position 
                      and a given detector row
        
        INPUTS :
        
            position :: Secondary grating position
            row :: Detector row
        
        OPTIONAL INPUTS: None
        
        OUTPUTS :

            slit_width :: Size of the slit (pixels)
        
        CALLING SEQUENCE:
        
            create_database(directory)
        
        MODIFICATION HISTORY : Juan Alday (15/07/2021)
        
    """

    from astropy.io import fits

    filename = (acs_dir.parent / "data" / "ACS_MIR_slit_size.fits").resolve()

    f = fits.open(filename)
    table = f[0].data

    pos_array = np.array([3,4,5,6,7,9,10,11,12,13])

    ipos = np.where(pos_array==position)[0]

    slitsize = float(table[ipos,row][0])

    f.close()

    return slitsize