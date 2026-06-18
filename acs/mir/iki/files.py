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


###############################################################################################

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

def read_geometry(filename,correct_boresight=False):

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

    nacq = file_lines(filename) - 3

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
