#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# latmos.files - Functions to extract transmission spectra from the LATMOS FITS files
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

###############################################################################################

def extract_order(filename,ordersel,irows,refalt=180.,iki_geometry=False,shift_row_edge=0):

    """
    FUNCTION NAME : extract_order()

    DESCRIPTION : Function to extract a set of detector rows within a selected diffraction order 
                   from the LATMOS FITS files
    
    INPUTS : 

        filename :: Name of the input FITS file
        ordersel :: Diffraction order
        irows(nrows) :: Index of row above the usable ones

    OPTIONAL INPUTS:
    
        refalt :: Reference above which to calculate the uncertainty for the topmost spectra (km)
        iki_geometry :: If True, the geometry is calculated using the IKI method. If False, the geometry is calculated using the LATMOS method.
        shift_row_edge :: Number of rows to shift the bottom edge of the slit for the computation of the geometry.
            
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

    MODIFICATION HISTORY : Juan Alday (02/08/2024)

    """
    
    from astropy.io import fits
    
    if filename[-5::]=='.fits':
        hdul = fits.open(filename)
        
    elif filename[-8::]=='.fits.fz':
        hdul = fits.open(filename)

    Observation = filename[len(filename)-40:len(filename)-8]
    IEflag = Observation[20]
    FPflag = Observation[-1]
    irows = np.array(irows,dtype='int32')

    #Reading File
    ##########################################################
    
    # Access the primary HDU (Header Data Unit)
    trans_hdu = hdul['COMPRESSED_IMAGE'] 
    
    ncols = trans_hdu.header['NAXIS1']   #Number of columns in matrix (i.e. number of wavelengths)
    nrows = trans_hdu.header['NAXIS2']   #Number of rows in detector matrix
    nacq = trans_hdu.header['NAXIS3']    #Number of acquisitions
    trans = trans_hdu.data               #Transmission spectra (nacq,nrows,ncols)

    date_utc_start = hdul['COMPRESSED_IMAGE'].header['DATE-BEG']   #Date and time at start of occultation
    date_utc_end = hdul['COMPRESSED_IMAGE'].header['DATE-END']     #Date and time at end of occultation
    mtp = hdul['COMPRESSED_IMAGE'].header['MTP']                   #MTP of the observation
    stp = hdul['COMPRESSED_IMAGE'].header['STP']                   #STP of the observation
    orbit = hdul['COMPRESSED_IMAGE'].header['ORBIT']               #Orbit number
    sequence = hdul['COMPRESSED_IMAGE'].header['SEQUENCE']         #Sequence (N1 or N2)
    position = hdul['COMPRESSED_IMAGE'].header['POS_1']                 #Secondary grating position
    
    #Reading CALIB
    calib_hdu = hdul['CALIB']
    
    calib_data = calib_hdu.data
    orders = calib_data['N'][0]           #Name of the orders 
    norders = len(orders)                 #Number of diffraction orders
    irow_slitbottom = calib_data['Y'][0]  #Row index of the bottom end of the slit
    irow_slitbottom += shift_row_edge
    irow_low = calib_data['Y_LOW'][0]     #Suggested nominal lower row index
    irow_high = calib_data['Y_HIGH'][0]   #Suggested nominal higher row index
    wavel = calib_data['WL'][0]           #(norders,ncols) wavelength in nm
    waven = 1./ (wavel/1.0e3) * 1.0e4     #(norders,ncols) wavenumber in cm-1

    iorder = -1
    for i in range(norders):
        if orders[i]==ordersel:
            iorder = i
    if iorder==-1:
        print('available orders :: ',orders)
        raise ValueError('error :: order '+str(ordersel)+' does not exist in this observation')

    #Reading GEO_BOTTOM
    geob_hdu = hdul['GEO_BOTTOM']

    sza_bottom = geob_hdu.data['SZA'][0]   #degrees
    alt_bottom = geob_hdu.data['ALT'][0]   #tangent altitude w/r to MOLA topography (km)
    lat_bottom = geob_hdu.data['LAT_PCE'][0]   #degrees
    lon_bottom = geob_hdu.data['LON_PCE'][0]   #degrees
    lst_bottom = geob_hdu.data['LST'][0]   #hours

    ellips_bottom = geob_hdu.data['SUB_ELPS'][0]   #ellipsoid radius at sub-tangent point (km)
    areoid_bottom = geob_hdu.data['SUB_ARE'][0]   #areoid radius at sub-tangent point (km)
    mola_bottom = geob_hdu.data['SUB_MOLA'][0]   #MOLA radius at sub-tangent point (km)

    tanhe_ellipsoid_bottom = (alt_bottom +  mola_bottom) - ellips_bottom  #tangent altitude wrt ellipsoid (km)
    tanhe_areoid_bottom = (alt_bottom +  mola_bottom) - areoid_bottom  #tangent altitude wrt ellipsoid (km)

    #Reading GEO_GEN
    geoc_hdu = hdul['GEO_CENTER']

    sza_cen = geoc_hdu.data['SZA'][0]   #degrees
    alt_cen = geoc_hdu.data['ALT'][0]   #tangent altitude w/r to MOLA topography (km)
    lat_cen = geoc_hdu.data['LAT_PCE'][0]   #degrees
    lon_cen = geoc_hdu.data['LON_PCE'][0]   #degrees
    lst_cen = geoc_hdu.data['LST'][0]   #hours

    ellips_cen = geoc_hdu.data['SUB_ELPS'][0]   #ellipsoid radius at sub-tangent point (km)
    areoid_cen = geoc_hdu.data['SUB_ARE'][0]   #areoid radius at sub-tangent point (km)
    mola_cen = geoc_hdu.data['SUB_MOLA'][0]   #MOLA radius at sub-tangent point (km)

    tanhe_ellipsoid_cen = (alt_cen +  mola_cen) - ellips_cen  #tangent altitude wrt ellipsoid (km)
    tanhe_areoid_cen = (alt_cen +  mola_cen) - areoid_cen  #tangent altitude wrt areoid (km)

    #Reading GEO_SAT
    SUN_DIST = hdul['GEO_SAT'].header['SUN_DIST']   #Sun-Mars distance (km)
    Ls = hdul['GEO_SAT'].header['L_S']        #Solar longitude (degrees)
    
    # Close the FITS file
    hdul.close()
    
    #Reading the IKI geometry to compare or use
    if iki_geometry is True:
        
        #filename_geom1 = 'GEO_1_MIR_0B_'+Observation[7:29]+'_SWP.txt'    #Slit centre
        #filename_geom2 = 'GEO_2_MIR_0B_'+Observation[7:29]+'_SWP.txt'    #Slit edge (upper part - not visible in images)
        #filename_geom3 = 'GEO_3_MIR_0B_'+Observation[7:29]+'_SWP.txt'    #Slit edge (lower part)

        geomdir='/exomars/data/external/tgo/acs/Geometry/'
        mtpdir = 'MTP'+str(mtp).zfill(3)
        stpdir = 'STP'+str(stp).zfill(3)
        Geomdir = geomdir+mtpdir+'/'+stpdir+'/MIR_geometry_pr_tc/'
        
        read_files = True
        if os.path.isfile(Geomdir+'GEO_1_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P1_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'):
            file_slit_centre = 'GEO_1_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P1_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'
            file_slit_bottom = 'GEO_3_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P1_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'
        elif os.path.isfile(Geomdir+'GEO_1_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P2_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'):
            file_slit_centre = 'GEO_1_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P2_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'
            file_slit_bottom = 'GEO_3_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P2_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'
        elif os.path.isfile(Geomdir+'GEO_1_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P4_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'):
            file_slit_centre = 'GEO_1_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P4_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'
            file_slit_bottom = 'GEO_3_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P4_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt'
        else:
            print(Geomdir+'GEO_1_MIR_0B_ORB'+str(orbit).zfill(6)+'_N'+str(sequence)+'_'+IEflag+'_P2_'+str(position).zfill(2)+'_'+FPflag+'_SWP.txt', 'does not exist')
            read_files = False
            
        if read_files is True:
            correct_boresight = True   #We apply the correction to the boresight vector
            print('Using IKI geometry')
            nacq1,lat_tgo1,lon_tgo1,lat_cen,lon_cen,lat_subsolar_cen,lon_subsolar_cen,alt_tgo_cen,\
                tanhe_areoid_cen,tanhe_ellipsoid_cen,alt_topo_areoid_cen,lst_cen = ps.acs.geometry.read_geometry_pr_tc(Geomdir+file_slit_centre,correct_boresight=correct_boresight)

            nacq3,lat_tgo3,lon_tgo3,lat_bottom,lon_bottom,lat_subsolar_bottom,lon_subsolar_bottom,alt_tgo_bottom,\
                tanhe_areoid_bottom,tanhe_ellipsoid_bottom,alt_topo_areoid_bottom,lst_bottom = ps.acs.geometry.read_geometry_pr_tc(Geomdir+file_slit_bottom,correct_boresight=correct_boresight)

            if( (nacq1==nacq+1) or (nacq1==nacq)):
                lat_cen = lat_cen[0:nacq]
                lon_cen = lon_cen[0:nacq]
                lst_cen = lst_cen[0:nacq]
                lat_bottom = lat_bottom[0:nacq]
                lon_bottom = lon_bottom[0:nacq]
                lst_bottom = lst_bottom[0:nacq]
                tanhe_areoid_cen = tanhe_areoid_cen[0:nacq]
                tanhe_ellipsoid_cen = tanhe_ellipsoid_cen[0:nacq]
                tanhe_areoid_bottom = tanhe_areoid_bottom[0:nacq]
                tanhe_ellipsoid_bottom = tanhe_ellipsoid_bottom[0:nacq]
            else:
                print('IKI acquisitions :: ',nacq1)
                print('LATMOS acquisitions ::',nacq)
                raise ValueError('error :: number of acquisitions in LATMOS and IKI geometry files should be the same')


    #Re-shaping the transmission data and calculating geometry for each row
    #######################################################################################
    
    #In this version of the ACS MIR data the size of the slit is fixed to 34.3 rows
    #Therefore, knowing the position of the slit bottom we can easily calculate the 
    #geometry at any arbitrary row
    slit_size = 34.3 #pixels
    
    tanhe_ellipsoid_per_row = (tanhe_ellipsoid_cen-tanhe_ellipsoid_bottom)/(slit_size/2.)
    tanhe_areoid_per_row = (tanhe_areoid_cen-tanhe_areoid_bottom)/(slit_size/2.)

    lat_per_row = (lat_cen-lat_bottom)/(slit_size/2.)
    lon_per_row = (lon_cen-lon_bottom)/(slit_size/2.)
    lst_per_row = (lst_cen-lst_bottom)/(slit_size/2.)
    
    #Re-shaping the data
    nrows_use = irow_high[iorder] - irow_low[iorder]   #Usable number of rows per diffraction order
    
    if np.max(irows) > nrows_use-1:
        print('warning :: some of the desired rows are not included in the order', filename,':: order', ordersel)
        print(np.where(irows<nrows_use)[0])
        irows = irows[np.where(irows<nrows_use)[0]]
        if len(irows)==0:
            print('irows :: ',irows)
            print('number of available rows :: ',nrows_use)
            raise ValueError('error')
    
    nrowssel = len(irows)
    
    trans_use = np.zeros((ncols,nacq,nrowssel))
    tanhe_areoid_use = np.zeros((nacq,nrowssel))
    tanhe_ellipsoid_use = np.zeros((nacq,nrowssel))
    lat_use = np.zeros((nacq,nrowssel))
    lon_use = np.zeros((nacq,nrowssel))
    lst_use = np.zeros((nacq,nrowssel))
    
    if nrowssel>0:

        for irow in range(nrowssel):
            trans_use[:,:,irow] = trans[:,irow_low[iorder]+irows[irow],:].T

            tanhe_ellipsoid_use[:,irow] = tanhe_ellipsoid_bottom[:] + tanhe_ellipsoid_per_row[:] * (irow_low[iorder] + irows[irow] - irow_slitbottom[iorder])
            tanhe_areoid_use[:,irow] = tanhe_areoid_bottom[:] + tanhe_areoid_per_row[:] * (irow_low[iorder] + irows[irow] - irow_slitbottom[iorder])
            lat_use[:,irow] = lat_bottom[:] + lat_per_row[:] * (irow_low[iorder] + irows[irow] - irow_slitbottom[iorder])
            lon_use[:,irow] = lon_bottom[:] + lon_per_row[:] * (irow_low[iorder] + irows[irow] - irow_slitbottom[iorder])
            lst_use[:,irow] = lst_bottom[:] + lst_per_row[:] * (irow_low[iorder] + irows[irow] - irow_slitbottom[iorder])
    
    else:
        
        raise ValueError('error :: There are no usable rows for the desired order')

    #Sorting the array to be in ascending order from bottom to top of atmosphere
    ###############################################################################

    isort = np.argsort(tanhe_areoid_use[:,0])

    trans_use = trans_use[:,isort,:]
    tanhe_areoid_use = tanhe_areoid_use[isort,:]
    tanhe_ellipsoid_use = tanhe_ellipsoid_use[isort,:]
    lat_use = lat_use[isort,:]
    lon_use = lon_use[isort,:]
    lst_use = lst_use[isort,:]

    #Removing the altitudes with nan frames
    ################################################################################

    #Removing the frames with nan values
    irow = int(nrowssel/2)
    inotnan = np.where(np.isnan(trans_use[int(waven.shape[1]/2),:,irow])==False)[0]

    trans_use = trans_use[:,inotnan,:]
    tanhe_areoid_use = tanhe_areoid_use[inotnan,:]
    tanhe_ellipsoid_use = tanhe_ellipsoid_use[inotnan,:]
    lat_use = lat_use[inotnan,:]
    lon_use = lon_use[inotnan,:]
    lst_use = lst_use[inotnan,:]
    
    inotnan = np.where(np.isnan(tanhe_areoid_use[:,irow])==False)[0]
    if len(inotnan)==0:
        print('error :: no usable rows in the order', ordersel, 'for observation', filename)
        raise ValueError('error :: no usable rows in the order '+str(ordersel)+' for observation '+filename)
    
    trans_use = trans_use[:,inotnan,:]
    tanhe_areoid_use = tanhe_areoid_use[inotnan,:]
    tanhe_ellipsoid_use = tanhe_ellipsoid_use[inotnan,:]
    lat_use = lat_use[inotnan,:]
    lon_use = lon_use[inotnan,:]
    lst_use = lst_use[inotnan,:]    

    #Calculating the average geometry for each acquisition
    ############################################################################
    
    tanhe_areoid_averow = np.zeros(trans_use.shape[1])
    tanhe_ellipsoid_averow = np.zeros(trans_use.shape[1])
    lat_averow = np.zeros(trans_use.shape[1])
    lon_averow = np.zeros(trans_use.shape[1])
    lst_averow = np.zeros(trans_use.shape[1])
    
    for iacq in range(trans_use.shape[1]):
        vals_areoid = []
        vals_ellipsoid = []
        vals_lat = []
        vals_lon = []
        vals_lst = []
        for irow in range(nrowssel): 
            vals_areoid.append(tanhe_areoid_use[iacq,irow])
            vals_ellipsoid.append(tanhe_ellipsoid_use[iacq,irow])
            vals_lat.append(lat_use[iacq,irow])
            vals_lon.append(lon_use[iacq,irow])
            vals_lst.append(lst_use[iacq,irow])
        
        tanhe_areoid_averow[iacq] = np.mean(vals_areoid)
        tanhe_ellipsoid_averow[iacq] = np.mean(vals_ellipsoid)
        lat_averow[iacq] = np.mean(vals_lat)
        lon_averow[iacq] = np.mean(vals_lon)
        lst_averow[iacq] = np.mean(vals_lst)
        
    #Estimating the uncertainty from the topmost altitudes
    ############################################################################

    transerr_top = np.ones((ncols,nrowssel)) * np.nan
        
    for irow in range(nrowssel):

        iin = np.where( (tanhe_ellipsoid_use[:,irow]>=refalt) & (tanhe_ellipsoid_use[:,irow]<=270.))[0]
        
        ntop = len(iin)
        vals = np.ones((ncols,ntop))
        
        for i in range(ntop):
        
            itan = iin[i]
            
            #Smoothing the transmission with a 10-pixel box kernel
            kernel_size = 30
            kernel = np.ones(kernel_size) / kernel_size
            smoothed = np.convolve(trans_use[:,itan,irow], kernel, mode='same')
        
            #Computing the uncertainty as the difference between the smoothed and non-smoothed
            error = np.abs(smoothed-trans_use[:,itan,irow])
        
            #For each value of the 
            vals[:,i] = trans_use[:,itan,irow]/smoothed - 1.

        mu = np.mean(vals,axis=1)
        std = np.std(vals,axis=1)
        error = np.max(np.abs(vals),axis=1)
        
        kernel_size = 15
        kernel = np.ones(kernel_size) / kernel_size
        smoothed = np.convolve(error, kernel, mode='same')
        
        inotnan_error = np.where(np.isnan(smoothed)==False)[0]
        inotnan_wave = np.where(np.isnan(trans_use[:,itan,irow])==False)[0]
        
        transerr_top[inotnan_wave,irow] = np.interp(waven[iorder,inotnan_wave],waven[iorder,inotnan_error],smoothed[inotnan_error])


    #Propagating the errors to the rest of the altitudes
    ###################################################################################

    errmeas_v1 = np.zeros(trans_use.shape)
    errmeas_v2 = np.zeros(trans_use.shape)
    transerr = np.zeros(trans_use.shape)
    frac = 0.5 #fraction to be applied between the v1 and v2 methods for propagation of the errors in the transmission level

    for irow in range(nrowssel):

        for itan in range(trans_use.shape[1]):
        
            errmeas_v1[:,itan,irow] = np.abs(transerr_top[:,irow] * np.sqrt( (1.+trans_use[:,itan,irow]**2.)/2. ))
            errmeas_v2[:,itan,irow] = np.abs(transerr_top[:,irow] * np.sqrt( trans_use[:,itan,irow]*(1.+trans_use[:,itan,irow]**2.)/2. ))
            transerr[:,itan,irow] = frac*errmeas_v1[:,itan,irow] + (1.-frac)*errmeas_v2[:,itan,irow]
    
    
    #Calculating the lat,lon and Loct in the level where the tangent height is 30
    ival = np.argmin(np.abs(tanhe_areoid_averow-30.))
    lat = lat_averow[ival]
    lon = lon_averow[ival]
    Loct = lst_averow[ival]
    
    return lat,lon,Ls,Loct,waven[iorder,:],trans_use,transerr,tanhe_areoid_use


###############################################################################################

def check_rows_order(filename,ordersel,irows):

    """
    FUNCTION NAME : check_rows_order()

    DESCRIPTION : Check that the desired rows in the ACS MIR file and order are available
    
    INPUTS : 

        filename :: Name of the input binary file
        order :: Diffraction order
        rows :: Index of row above the usable ones

    OPTIONAL INPUTS: None
            
    OUTPUTS : 
 
        exist :: Flag indicating whether the rows exist
        nrows :: Number of rows that exist

    CALLING SEQUENCE:

        exist,nrows = check_rows_acsmir_latmos_order(filename,ordersel,irows)

    MODIFICATION HISTORY : Juan Alday (02/08/2024)

    """
    
    from astropy.io import fits
    
    if filename[-5::]=='.fits':
        hdul = fits.open(filename)
        
    elif filename[-8::]=='.fits.fz':
        hdul = fits.open(filename)

    Observation = filename[len(filename)-40:len(filename)-8]
    IEflag = Observation[20]
    FPflag = Observation[-1]
    irows = np.array(irows,dtype='int32')

    #Reading File
    ##########################################################
    
    # Access the primary HDU (Header Data Unit)
    trans_hdu = hdul['COMPRESSED_IMAGE'] 
    
    ncols = trans_hdu.header['NAXIS1']   #Number of columns in matrix (i.e. number of wavelengths)
    nrows = trans_hdu.header['NAXIS2']   #Number of rows in detector matrix
    nacq = trans_hdu.header['NAXIS3']    #Number of acquisitions
    trans = trans_hdu.data               #Transmission spectra (nacq,nrows,ncols)

    date_utc_start = hdul['COMPRESSED_IMAGE'].header['DATE-BEG']   #Date and time at start of occultation
    date_utc_end = hdul['COMPRESSED_IMAGE'].header['DATE-END']     #Date and time at end of occultation
    mtp = hdul['COMPRESSED_IMAGE'].header['MTP']                   #MTP of the observation
    stp = hdul['COMPRESSED_IMAGE'].header['STP']                   #STP of the observation
    orbit = hdul['COMPRESSED_IMAGE'].header['ORBIT']               #Orbit number
    sequence = hdul['COMPRESSED_IMAGE'].header['SEQUENCE']         #Sequence (N1 or N2)
    position = hdul['COMPRESSED_IMAGE'].header['POS_1']                 #Secondary grating position
    
    #Reading CALIB
    calib_hdu = hdul['CALIB']
    
    calib_data = calib_hdu.data
    orders = calib_data['N'][0]           #Name of the orders 
    norders = len(orders)                 #Number of diffraction orders
    irow_slitbottom = calib_data['Y'][0]  #Row index of the bottom end of the slit
    irow_low = calib_data['Y_LOW'][0]     #Suggested nominal lower row index
    irow_high = calib_data['Y_HIGH'][0]   #Suggested nominal higher row index
    wavel = calib_data['WL'][0]           #(norders,ncols) wavelength in nm
    waven = 1./ (wavel/1.0e3) * 1.0e4     #(norders,ncols) wavenumber in cm-1

    # Close the FITS file
    hdul.close()

    iorder = -1
    for i in range(norders):
        if orders[i]==ordersel:
            iorder = i
    if iorder==-1:
        print('available orders :: ',orders)
        raise ValueError('error :: order '+str(ordersel)+' does not exist in observation '+filename)

    if irow_high[iorder]==-1:
        print('orders :: ',orders)
        print('irow_high :: ',irow_high)
        print('error :: order '+str(ordersel)+' does not exist in observation '+filename+' (likely partial frame)')
        exist=False
        return exist,0

    #Checking the number of rows that are usable
    nrows_use = irow_high[iorder] - irow_low[iorder]   #Usable number of rows per diffraction order

    if nrows_use==0:
        exist=False
        return exist,nrows_use
    
    if irows.max()>nrows_use-1:
        exist=False
        return exist,nrows_use
    
    exist=True
    return exist,nrows_use
        
