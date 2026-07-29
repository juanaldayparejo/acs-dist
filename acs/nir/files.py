#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# files - Functions to extract transmission spectra from the NIR files processed by Anna Fedorova
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
import struct

###############################################################################################

def read_clb(filename):
    """
    FUNCTION NAME : read_clb()

    DESCRIPTION : Function to read the ACS NIR calibrated files (.clb)
    
    INPUTS : 

        filename :: Name of the input calibrated file

    OPTIONAL INPUTS:
    
    OUTPUTS : 
 
        ngeom :: Number of geometries
        norders :: Number of diffraction orders
        orders(norders) :: Diffraction orders
        geo :: Dictionary with geometry parameters
            geo["latitude"] = (ngeom,norders) :: Latitude of each spectrum
            geo["longitude"] = (ngeom,norders) :: Longitude of each spectrum
            geo["Ls"] = (ngeom,norders) :: Ls of each spectrum
            geo["LST"] = (ngeom,norders) :: Local time of each spectrum
            geo["TANHE"] = (ngeom,norders) :: Tangent height (km)
            geo["TANHE_AREOID"] = (ngeom,norders) :: Tangent height above the areoid (km)

        data :: Dictionary with data parameters
            data["wlm2"] = (640,ngeom,norders) :: Wavenumber array for order -2
            data["wlm1"] = (640,ngeom,norders) :: Wavenumber array for order -1
            data["wl0"] = (640,ngeom,norders) :: Wavenumber array for main order 0
            data["wlp1"] = (640,ngeom,norders) :: Wavenumber array for order +1
            data["wlp2"] = (640,ngeom,norders) :: Wavenumber array for order +2
            data["trans1"] = (640,ngeom,norders) :: Transmission spectra (initial version of calibration)
            data["trans2"] = (640,ngeom,norders) :: Transmission spectra (second version of calibration)
            data["trans3"] = (640,ngeom,norders) :: Transmission spectra (third version of calibration)
            data["trans_error"] = (640,ngeom,norders) :: Uncertainty in transmission spectra

    CALLING SEQUENCE:

        ngeom, norders, orders, geo, data = read_clb(filename)

    MODIFICATION HISTORY : Juan Alday (13/07/2026)

    """
    
    with open(filename, 'rb') as fid:

        SeriesCount, PointsCount = struct.unpack('<hh', fid.read(4))

        ordermy = np.zeros(PointsCount, dtype=np.int16)

        geo = {
            'latitude': np.zeros((SeriesCount, PointsCount), dtype=np.float32),
            'longitude': np.zeros((SeriesCount, PointsCount), dtype=np.float32),
            'ls': np.zeros((SeriesCount, PointsCount), dtype=np.float32),
            'lst': np.zeros((SeriesCount, PointsCount), dtype=np.float32),
            'tangent_height': np.zeros((SeriesCount, PointsCount), dtype=np.float32),
            'tangent_height_areoid': np.zeros((SeriesCount, PointsCount), dtype=np.float32),
            'wlfr_final': np.zeros((SeriesCount, PointsCount), dtype=np.float32),
        }

        data = {key: np.zeros((640, SeriesCount, PointsCount), dtype=np.float32)
                for key in ['wlm2', 'wlm1', 'wl0', 'wlp1', 'wlp2', \
                            'trans1', 'trans2', 'trans3', 'trans_error']}

        for i in range(SeriesCount):
            for j in range(PointsCount):
                # Read i1, j1 (int16)
                fid.read(4)

                # Read ordermy (int16)
                ordermy[j], = struct.unpack('<h', fid.read(2))

                # Read 7 floats (float32)
                ah = struct.unpack('<7f', fid.read(28))
                geo['latitude'][i, j], geo['longitude'][i, j], geo['ls'][i, j], geo['lst'][i, j], \
                geo['tangent_height'][i, j], geo['tangent_height_areoid'][i, j], geo['wlfr_final'][i, j] = ah

                for key in ['wlm2', 'wlm1', 'wl0', 'wlp1', 'wlp2',
                            'trans1', 'trans2', 'trans3', 'trans_error']:
                    buf = fid.read(640 * 4)
                    arr = struct.unpack('<640f', buf)
                    data[key][:, i, j] = arr

    return SeriesCount, PointsCount, ordermy, geo, data

###############################################################################################

def extract_order_clb(filename,order):
    """
    FUNCTION NAME : extract_order_clb()

    DESCRIPTION : Extract information about one diffraction order from the ACS NIR calibrated files (.clb)
    
    INPUTS : 

        filename :: Name of the input calibrated file

    OPTIONAL INPUTS:
    
    OUTPUTS : 
 
        ngeom :: Number of geometries
        geo :: Dictionary with geometry parameters
            geo["latitude"] = (ngeom) :: Latitude of each spectrum
            geo["longitude"] = (ngeom) :: Longitude of each spectrum
            geo["Ls"] = (ngeom) :: Ls of each spectrum
            geo["LST"] = (ngeom) :: Local time of each spectrum
            geo["TANHE"] = (ngeom) :: Tangent height (km)
            geo["TANHE_AREOID"] = (ngeom) :: Tangent height above the areoid (km)

        data :: Dictionary with data parameters
            data["wlm2"] = (640,ngeom) :: Wavenumber array for order -2
            data["wlm1"] = (640,ngeom) :: Wavenumber array for order -1
            data["wl0"] = (640,ngeom) :: Wavenumber array for main order 0
            data["wlp1"] = (640,ngeom) :: Wavenumber array for order +1
            data["wlp2"] = (640,ngeom) :: Wavenumber array for order +2
            data["trans1"] = (640,ngeom) :: Transmission spectra (initial version of calibration)
            data["trans2"] = (640,ngeom) :: Transmission spectra (second version of calibration)
            data["trans3"] = (640,ngeom) :: Transmission spectra (third version of calibration)
            data["trans_error"] = (640,ngeom) :: Uncertainty in transmission spectra

    CALLING SEQUENCE:

        ngeom geo, data = read_acsnir_clb(filename,order)

    MODIFICATION HISTORY : Juan Alday (13/07/2026)

    """
    
    with open(filename, 'rb') as fid:

        SeriesCount, PointsCount = struct.unpack('<hh', fid.read(4))

        geo = {
            'latitude': np.zeros((SeriesCount), dtype=np.float32),
            'longitude': np.zeros((SeriesCount), dtype=np.float32),
            'ls': np.zeros((SeriesCount), dtype=np.float32),
            'lst': np.zeros((SeriesCount), dtype=np.float32),
            'tangent_height': np.zeros((SeriesCount), dtype=np.float32),
            'tangent_height_areoid': np.zeros((SeriesCount), dtype=np.float32),
            'wlfr_final': np.zeros((SeriesCount), dtype=np.float32),
        }

        data = {key: np.zeros((640, SeriesCount), dtype=np.float32)
                for key in ['wlm2', 'wlm1', 'wl0', 'wlp1', 'wlp2', \
                            'trans1', 'trans2', 'trans3', 'trans_error']}

        for i in range(SeriesCount):
            for j in range(PointsCount):
                # Read i1, j1 (int16)
                fid.read(4)

                # Read ordermy (int16)
                ordermy, = struct.unpack('<h', fid.read(2))

                # Read 7 floats (float32)
                ah = struct.unpack('<7f', fid.read(28))
                if ordermy == order:
                    geo['latitude'][i], geo['longitude'][i], geo['ls'][i], geo['lst'][i], \
                    geo['tangent_height'][i], geo['tangent_height_areoid'][i], geo['wlfr_final'][i] = ah

                for key in ['wlm2', 'wlm1', 'wl0', 'wlp1', 'wlp2',
                            'trans1', 'trans2', 'trans3', 'trans_error']:
                    buf = fid.read(640 * 4)
                    arr = struct.unpack('<640f', buf)
                    if ordermy == order:
                        data[key][:, i] = arr

    return SeriesCount, geo, data