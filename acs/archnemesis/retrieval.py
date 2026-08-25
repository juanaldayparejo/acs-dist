#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# retrieval - Set of functions to perform an archnemesis retrieval
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

import archnemesis as ans
import matplotlib.pyplot as plt
import numpy as np
import time
import h5py
from copy import copy,deepcopy
from scipy.optimize import minimize
import sys,os

#####################################################################################
#####################################################################################
#                                        nemesisACS
#####################################################################################
#####################################################################################

# Version of Nemesis for doing retrievals in ACS solar occultation observations

#####################################################################################
#####################################################################################

def run_retrieval(dir_atm,runname_atm,dir_meas,runname_meas):

    ######################################################
    ######################################################
    #    READING INPUT FILES AND SETTING UP VARIABLES
    ######################################################
    ######################################################

    NCores = 20

    smooth_baseline = None
    normalise_spectra=False
    poldegree=None
    fourierdegree=None

    remove_window = False
    remove_vmin = 3764.5
    remove_vmax = 3789.72

    start = time.time()

    #Initialise Atmosphere class and read file
    ##############################################################

    Atmosphere = ans.Atmosphere_0()

    #Read gaseous atmosphere
    Atmosphere.read_hdf5(dir_atm+'/'+runname_atm)

    #Initialise Layer class and read file
    ###############################################################

    Layer = ans.Layer_0(Atmosphere.RADIUS)
    Layer.read_hdf5(dir_atm+'/'+runname_atm)

    #Initialise Surface class and read file
    ###############################################################

    isurf = ans.Data.planet_data.planet_info[str(Atmosphere.IPLANET)]["isurf"]
    Surface = ans.Surface_0()
    if isurf==1:
        Surface.read_hdf5(dir_atm+'/'+runname_atm)
    else:
        Surface.GASGIANT=True

    #Initialise Scatter class and read file
    ###############################################################

    Scatter = ans.Scatter_0()
    Scatter.read_hdf5(dir_atm+'/'+runname_atm)

    #Initialise CIA class and read files (.cia)  - NOT FROM HDF5 YET
    ##############################################################

    CIA = None

    #Initialise Spectroscopy class and read file
    ###############################################################

    f = h5py.File(dir_meas+'/'+runname_meas+'.h5','r')
    #Checking if Spectroscopy exists
    e = "/Spectroscopy" in f
    f.close()

    if e is True:
        Spectroscopy = ans.Spectroscopy_0()
        Spectroscopy.read_hdf5(dir_meas+'/'+runname_meas)
    else:
        Spectroscopy = None



    #Initialise Measurement class and read file
    ###############################################################

    Measurement = ans.Measurement_0()
    Measurement.read_hdf5(dir_meas+'/'+runname_meas)

    if normalise_spectra==True:
        Measurement.MEAS = (Measurement.MEAS[:,:].T / Measurement.MEAS[:,Measurement.NGEOM-1]).T
        
    if remove_window==True:
        
        #Removing the specified pixels
        irem = np.where(~( (Measurement.VCONV[:,0]>=remove_vmin) & (Measurement.VCONV[:,0]<=remove_vmax) ))[0]

        Measurement.NCONV[:] = len(irem)
        Measurement.edit_VCONV(Measurement.VCONV[irem,:])
        Measurement.edit_MEAS(Measurement.MEAS[irem,:])
        Measurement.edit_ERRMEAS(Measurement.ERRMEAS[irem,:])

        Measurement.NFIL = Measurement.NFIL[irem]
        Measurement.VFIL = Measurement.VFIL[:,irem]
        Measurement.AFIL = Measurement.AFIL[:,irem]

    Measurement.calc_MeasurementVector()

    if Spectroscopy is not None:

        #Calculating the 'calculation wavelengths'
        vmin, vmax = Measurement.calc_wave_range()

        #Now, reading k-tables or lbl-tables for the spectral range of interest
        Spectroscopy.read_tables(wavemin=vmin,wavemax=vmax)
        
    else:
        
        vmin, vmax = Measurement.calc_wave_range()
        wavex = np.arange(vmin,vmax,0.05)
        
        #Creating dummy Spectroscopy file if it does not exist
        Spectroscopy = ans.Spectroscopy_0()
        Spectroscopy.NWAVE = len(wavex)
        Spectroscopy.WAVE = wavex
        Spectroscopy.NG = 1
        Spectroscopy.ILBL = 2
        Spectroscopy.NGAS = 1
        Spectroscopy.ID = np.array([Atmosphere.ID[0]],dtype='int32')
        Spectroscopy.ISO = np.array([Atmosphere.ISO[0]],dtype='int32')
        Spectroscopy.NP = 2
        Spectroscopy.NT = 2
        Spectroscopy.PRESS = np.array([Atmosphere.P.min()/101325.,Atmosphere.P.max()/101325.])
        Spectroscopy.TEMP = np.array([Atmosphere.T.min(),Atmosphere.T.max()])
        Spectroscopy.K = np.zeros([Spectroscopy.NWAVE,Spectroscopy.NP,Spectroscopy.NT,Spectroscopy.NGAS])
        Spectroscopy.DELG = np.array([1])

    #Reading Stellar class
    ################################################################

    Stellar = ans.Stellar_0()
    Stellar.read_hdf5(dir_atm+'/'+runname_atm)

    #Reading .apr file and Variables Class
    #################################################################

    os.chdir(dir_meas)
    Variables = ans.Variables_0()
    Variables.read_apr(runname_meas,Atmosphere.NP,nlocations=1,ngas=Atmosphere.NVMR,ndust=Atmosphere.NDUST)


    Variables.XN = deepcopy(Variables.XA)
    Variables.SX = deepcopy(Variables.SA)


    #Reading retrieval setup
    #################################################################

    Retrieval = ans.OptimalEstimation_0()
    Retrieval.read_hdf5(dir_atm+'/'+runname_atm)
    #Retrieval.NITER = 30

    ######################################################
    ######################################################
    #      RUN THE RETRIEVAL USING ANY APPROACH
    ######################################################
    ######################################################

    if Retrieval.IRET==0: #(0) Optimal Estimation
        OptimalEstimation = coreretOE_ACS(runname_meas,Variables,Measurement,Atmosphere,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,\
                                        NITER=Retrieval.NITER,PHILIMIT=Retrieval.PHILIMIT,nemesisSO=True,NCores=NCores,write_itr=False,\
                                        poldegree=poldegree,fourierdegree=fourierdegree,smooth_baseline=smooth_baseline)
        Retrieval = OptimalEstimation
        Retrieval.IRET = 0

    else:
        sys.exit('error in nemesisSO :: Retrieval scheme has not been implemented yet')

    ######################################################
    ######################################################
    #                WRITE OUTPUT FILES
    ######################################################
    ######################################################

    if Retrieval.IRET==0:
        Retrieval.write_output_hdf5(dir_meas+'/'+runname_meas,Variables,write_cov=False)

    #Finishing pogram
    end = time.time()
    print('Model run OK')
    print(' Elapsed time (s) = '+str(end-start))


##################################################################################################
# MAIN CONVERGENCE LOOP
##################################################################################################

def coreretOE_ACS(runname,Variables,Measurement,Atmosphere,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,\
                 NITER=10,PHILIMIT=0.1,NCores=1,nemesisSO=True,write_itr=True,poldegree=None,fourierdegree=None,smooth_baseline=None):


    """
        FUNCTION NAME : coreretOE()
        
        DESCRIPTION : 

            This subroutine runs the Optimal Estimation iterarive algorithm to solve the inverse
            problem and find the set of parameters that fit the spectral measurements and are closest
            to the a priori estimates of the parameters.

        INPUTS :
       
            runname :: Name of the Nemesis run
            Variables :: Python class defining the parameterisations and state vector
            Measurement :: Python class defining the measurements 
            Atmosphere :: Python class defining the reference atmosphere
            Spectroscopy :: Python class defining the spectroscopic parameters of gaseous species
            Scatter :: Python class defining the parameters required for scattering calculations
            Stellar :: Python class defining the stellar spectrum
            Surface :: Python class defining the surface
            CIA :: Python class defining the Collision-Induced-Absorption cross-sections
            Layer :: Python class defining the layering scheme to be applied in the calculations

        OPTIONAL INPUTS:

            NITER :: Number of iterations in retrieval
            PHILIMIT :: Percentage convergence limit. If the percentage reduction of the cost function PHI
                        is less than philimit then the retrieval is deemed to have converged.

            nemesisSO :: If True, the retrieval uses the function jacobian_nemesisSO(), adapated specifically
                         for solar occultation observations, rather than the more general jacobian_nemesis() function.
            poldegree :: If not None, the forward model and jacobian matrix are normalised with a polynomial function

        OUTPUTS :

            OptimalEstimation :: Python class defining all the variables required as input or output
                                 from the Optimal Estimation retrieval
 
        CALLING SEQUENCE:
        
            OptimalEstimation = coreretOE(runname,Variables,Measurement,Atmosphere,Spectroscopy,Scatter,Stellar,Surface,Layer)
 
        MODIFICATION HISTORY : Juan Alday (06/08/2021)

    """

    #Creating class and including inputs
    #############################################

    OptimalEstimation = ans.OptimalEstimation_0()

    OptimalEstimation.NITER = NITER
    OptimalEstimation.PHILIMIT = PHILIMIT
    OptimalEstimation.NX = Variables.NX
    OptimalEstimation.NY = Measurement.NY
    OptimalEstimation.edit_XA(Variables.XA)
    OptimalEstimation.edit_XN(Variables.XN)
    OptimalEstimation.edit_SA(Variables.SA)
    OptimalEstimation.edit_Y(Measurement.Y)
    OptimalEstimation.edit_SE(Measurement.SE)

    #Opening .itr file
    #################################################################

    if OptimalEstimation.NITER>0:
        if write_itr==True:
            fitr = open(runname+'.itr','w')
            fitr.write("\t %i \t %i \t %i\n" % (OptimalEstimation.NX,OptimalEstimation.NY,OptimalEstimation.NITER))

    #Calculate the first measurement vector and jacobian matrix
    #################################################################

    ForwardModel = ans.ForwardModel_0(runname=runname, Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables)
    print('nemesis :: Calculating Jacobian matrix KK')
    YN,KK = ForwardModel.jacobian_nemesis(NCores=NCores,nemesisSO=nemesisSO)
    
    #Correcting the measurement vector and jacobian matrix with polynomial degree
    ##################################################################################
    
    if((poldegree is not None) & (smooth_baseline is not None)):
        raise ValueError('error :: baseline must be either with poldegree or smooth_baseline, but not both')
    
    if poldegree is not None:
        
        ix = 0
        YN_BASELINE = np.zeros(YN.shape)
        for iGEOM in range(Measurement.NGEOM):
            
            SPECMOD = np.zeros(Measurement.NCONV[iGEOM])
            SPECMOD[:] = YN[ix:ix+Measurement.NCONV[iGEOM]]

            #MEAS = MOD * BASELINE
            #Calculating the baseline from the measurement
            polmeas = Measurement.MEAS[0:Measurement.NCONV[iGEOM],iGEOM] / SPECMOD[:]
            
            #Fitting the baseline with a polynomial function
            pcoeff = np.polyfit(Measurement.VCONV[:,iGEOM],polmeas,poldegree)
            polmod1 = np.polyval(pcoeff,Measurement.VCONV[:,iGEOM])
            
            if fourierdegree is not None:
                
                #Calculating the new baseline after removing polynomial
                polmeas = Measurement.MEAS[0:Measurement.NCONV[iGEOM],iGEOM] / (SPECMOD[:] * polmod1[:])
                
                #Fitting the baseline with Fourier series
                initial_guess = np.ones(fourierdegree)
                
                # Minimize the loss function to obtain Fourier coefficients
                result = minimize(loss_function, initial_guess, args=(Measurement.VCONV[0:Measurement.NCONV[iGEOM],iGEOM], polmeas))
                
                # Get optimized parameters
                opt_params = result.x
        
                polmod = fourier_series(Measurement.VCONV[0:Measurement.NCONV[iGEOM],iGEOM], *opt_params)  
                
                polmodf = polmod * polmod1
                
            else:
                
                polmodf = polmod1   
            
            YN_BASELINE[ix:ix+Measurement.NCONV[iGEOM]] = polmodf[:]
    
            ix = ix + Measurement.NCONV[iGEOM]
    
        #Applying baseline correction to measurement vector
        YN[:] = YN[:] * YN_BASELINE[:]
        for ix in range(Variables.NX):
            KK[:,ix] = KK[:,ix] * YN_BASELINE[:]
            
    if smooth_baseline is not None:
        
        ix = 0
        YN_BASELINE = np.zeros(YN.shape)
        for iGEOM in range(Measurement.NGEOM):
            
            SPECMOD = np.zeros(Measurement.NCONV[iGEOM])
            SPECMOD[:] = YN[ix:ix+Measurement.NCONV[iGEOM]]

            #MEAS = MOD * BASELINE
            #Calculating the baseline from the measurement
            polmeas = Measurement.MEAS[0:Measurement.NCONV[iGEOM],iGEOM] / SPECMOD[:]
            
            #Smoothing the baseline to avoid the presence of absorption lines
            kernel_size = smooth_baseline
            kernel = np.ones(kernel_size) / kernel_size
            polmeas = np.convolve(polmeas, kernel, mode='same')
            
            #At the edges we get some spurious features due to the boundary conditions
            #We find those points and extrapolate
            s = interp1d(Measurement.VCONV[int(kernel_size/2):-int(kernel_size/2),iGEOM],polmeas[int(kernel_size/2):-int(kernel_size/2)],fill_value='extrapolate')
            polmeas = s(Measurement.VCONV[:,iGEOM])
            
            #fig,ax1 = plt.subplots(1,1,figsize=(10,3))
            #ax1.plot(Measurement.VCONV[:,iGEOM],polmeas)
            #ax1.plot(Measurement.VCONV[:,iGEOM],Measurement.MEAS[0:Measurement.NCONV[iGEOM],iGEOM])
            #ax1.grid()
            
            YN_BASELINE[ix:ix+Measurement.NCONV[iGEOM]] = polmeas[:]
    
            ix = ix + Measurement.NCONV[iGEOM]
            
        #Applying baseline correction to measurement vector
        YN[:] = YN[:] * YN_BASELINE[:]
        for ix in range(Variables.NX):
            KK[:,ix] = KK[:,ix] * YN_BASELINE[:]
    
    OptimalEstimation.edit_YN(YN)
    OptimalEstimation.edit_KK(KK)

    #Calculate gain matrix and average kernels
    #################################################################

    print('nemesis :: Calculating gain matrix')
    OptimalEstimation.calc_gain_matrix()

    #Calculate initial value of cost function phi
    #################################################################

    print('nemesis :: Calculating cost function')
    OptimalEstimation.calc_phiret()

    OPHI = OptimalEstimation.PHI
    print('chisq/ny = '+str(OptimalEstimation.CHISQ))

    #Assessing whether retrieval is going to be OK
    #################################################################

    OptimalEstimation.assess()

    #Run retrieval for each iteration
    #################################################################

    #Initializing some variables
    alambda = 1.0   #Marquardt-Levenberg-type 'braking parameter'
    NX11 = np.zeros(OptimalEstimation.NX)
    XN1 = deepcopy(OptimalEstimation.XN)
    NY1 = np.zeros(OptimalEstimation.NY)
    YN1 = deepcopy(OptimalEstimation.YN)

    for it in range(OptimalEstimation.NITER):

        print('nemesis :: Iteration '+str(it)+'/'+str(OptimalEstimation.NITER))

        if write_itr==True:
            
        #Writing into .itr file
        ####################################

            fitr.write('%10.5f %10.5f \n' % (OptimalEstimation.CHISQ,OptimalEstimation.PHI))
            for i in range(OptimalEstimation.NX):fitr.write('%10.5f \n' % (XN1[i]))
            for i in range(OptimalEstimation.NX):fitr.write('%10.5f \n' % (OptimalEstimation.XA[i]))
            for i in range(OptimalEstimation.NY):fitr.write('%10.5f \n' % (OptimalEstimation.Y[i]))
            for i in range(OptimalEstimation.NY):fitr.write('%10.5f \n' % (OptimalEstimation.SE[i,i]))
            for i in range(OptimalEstimation.NY):fitr.write('%10.5f \n' % (YN1[i]))
            for i in range(OptimalEstimation.NY):fitr.write('%10.5f \n' % (OptimalEstimation.YN[i]))
            for i in range(OptimalEstimation.NX):
                for j in range(OptimalEstimation.NY):fitr.write('%10.5f \n' % (OptimalEstimation.KK[j,i]))


        #Calculating next state vector
        #######################################

        print('nemesis :: Calculating next iterated state vector')
        X_OUT = OptimalEstimation.calc_next_xn()
        #  x_out(nx) is the next iterated value of xn using classical N-L
        #  optimal estimation. However, we want to apply a braking parameter
        #  alambda to stop the new trial vector xn1 being too far from the
        #  last 'best-fit' value xn

        IBRAKE = 0
        while IBRAKE==0: #We continue in this while loop until we do not find problems with the state vector
    
            for j in range(OptimalEstimation.NX):
                XN1[j] = OptimalEstimation.XN[j] + (X_OUT[j]-OptimalEstimation.XN[j])/(1.0+alambda)
                
                #Check to see if log numbers have gone out of range
                if Variables.LX[j]==1:
                    if((XN1[j]>85.) or (XN1[j]<-85.)):
                        print('nemesis :: log(number gone out of range) --- increasing brake')
                        alambda = alambda * 10.
                        IBRAKE = 0
                        if alambda>1.e30:
                            sys.exit('error in nemesis :: Death spiral in braking parameters - stopping')
                        break
                    else:
                        IBRAKE = 1
                else:
                    IBRAKE = 1
                    pass
                        
            if IBRAKE==0:
                continue
                        
            #Check to see if any VMRs or other parameters have gone negative.
            Variables1 = deepcopy(Variables)
            Variables1.XN = XN1

            ForwardModel1 = ans.ForwardModel_0(runname=runname, Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables1)
            #Variables1 = copy(Variables)
            #Variables1.XN = XN1
            #Measurement1 = copy(Measurement)
            #Atmosphere1 = copy(Atmosphere)
            #Scatter1 = copy(Scatter)
            #Stellar1 = copy(Stellar)
            #Surface1 = copy(Surface)
            #Spectroscopy1 = copy(Spectroscopy)
            #Layer1 = copy(Layer)
            #flagh2p = False
            #xmap = subprofretg(runname,Variables1,Measurement1,Atmosphere1,Spectroscopy1,Scatter1,Stellar1,Surface1,Layer1,flagh2p)
            ForwardModel1.subprofretg()

            #if(len(np.where(Atmosphere1.VMR<0.0))>0):
            #    print('nemesisSO :: VMR has gone negative --- increasing brake')
            #    alambda = alambda * 10.
            #    IBRAKE = 0
            #    continue
            
            #iwhere = np.where(Atmosphere1.T<0.0)
            iwhere = np.where(ForwardModel1.AtmosphereX.T<0.0)
            if(len(iwhere[0])>0):
                print('nemesis :: Temperature has gone negative --- increasing brake')
                alambda = alambda * 10.
                IBRAKE = 0
                continue


        #Calculate test spectrum using trial state vector xn1. 
        #Put output spectrum into temporary spectrum yn1 with
        #temporary kernel matrix kk1. Does it improve the fit? 
        Variables.edit_XN(XN1)
        print('nemesis :: Calculating Jacobian matrix KK')

        ForwardModel = ans.ForwardModel_0(runname=runname, Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables)
        YN1,KK1 = ForwardModel.jacobian_nemesis(NCores=NCores,nemesisSO=nemesisSO)


        #Correcting the measurement vector and jacobian matrix with polynomial degree
        ##################################################################################
        
        if poldegree is not None:
            
            ix = 0
            YN_BASELINE = np.zeros(YN.shape)
            for iGEOM in range(Measurement.NGEOM):
                
                SPECMOD = np.zeros(Measurement.NCONV[iGEOM])
                SPECMOD[:] = YN1[ix:ix+Measurement.NCONV[iGEOM]]

                #MEAS = MOD * BASELINE
                #Calculating the baseline from the measurement
                polmeas = Measurement.MEAS[0:Measurement.NCONV[iGEOM],iGEOM] / SPECMOD[:]
                
                #Fitting the baseline with a polynomial function
                pcoeff = np.polyfit(Measurement.VCONV[:,iGEOM],polmeas,poldegree)
                polmod1 = np.polyval(pcoeff,Measurement.VCONV[:,iGEOM])
                
                if fourierdegree is not None:
                    
                    #Calculating the new baseline after removing polynomial
                    polmeas = Measurement.MEAS[0:Measurement.NCONV[iGEOM],iGEOM] / ( SPECMOD[:] * polmod1)
                    
                    #Fitting the baseline with Fourier series
                    initial_guess = np.ones(fourierdegree)
                    
                    # Minimize the loss function to obtain Fourier coefficients
                    result = minimize(loss_function, initial_guess, args=(Measurement.VCONV[0:Measurement.NCONV[iGEOM],iGEOM], polmeas))
                    
                    # Get optimized parameters
                    opt_params = result.x
            
                    polmod = fourier_series(Measurement.VCONV[0:Measurement.NCONV[iGEOM],iGEOM], *opt_params)    
                    
                    polmodf = polmod * polmod1
                    
                else:
                    
                    polmodf = polmod1        
                
                YN_BASELINE[ix:ix+Measurement.NCONV[iGEOM]] = polmodf[:]
        
                ix = ix + Measurement.NCONV[iGEOM]
        
            #Applying baseline correction to measurement vector
            YN1[:] = YN1[:] * YN_BASELINE[:]
            for ix in range(Variables.NX):
                KK1[:,ix] = KK1[:,ix] * YN_BASELINE[:]
    
        if smooth_baseline is not None:
            
            ix = 0
            YN_BASELINE = np.zeros(YN.shape)
            for iGEOM in range(Measurement.NGEOM):
                
                SPECMOD = np.zeros(Measurement.NCONV[iGEOM])
                SPECMOD[:] = YN1[ix:ix+Measurement.NCONV[iGEOM]]

                #MEAS = MOD * BASELINE
                #Calculating the baseline from the measurement
                polmeas = Measurement.MEAS[0:Measurement.NCONV[iGEOM],iGEOM] / SPECMOD[:]
                
                #Smoothing the baseline to avoid the presence of absorption lines
                kernel_size = smooth_baseline
                kernel = np.ones(kernel_size) / kernel_size
                polmeas = np.convolve(polmeas, kernel, mode='same')
                
                #At the edges we get some spurious features due to the boundary conditions
                #We find those points and extrapolate
                s = interp1d(Measurement.VCONV[int(kernel_size/2):-int(kernel_size/2),iGEOM],polmeas[int(kernel_size/2):-int(kernel_size/2)],fill_value='extrapolate')
                polmeas = s(Measurement.VCONV[:,iGEOM])
                
                YN_BASELINE[ix:ix+Measurement.NCONV[iGEOM]] = polmeas[:]
        
                ix = ix + Measurement.NCONV[iGEOM]
                
            #Applying baseline correction to measurement vector
            YN1[:] = YN1[:] * YN_BASELINE[:]
            for ix in range(Variables.NX):
                KK1[:,ix] = KK1[:,ix] * YN_BASELINE[:]
    
        OptimalEstimation1 = deepcopy(OptimalEstimation)
        OptimalEstimation1.edit_YN(YN1)
        OptimalEstimation1.edit_XN(XN1)
        OptimalEstimation1.edit_KK(KK1)
        OptimalEstimation1.calc_phiret()
        print('chisq/ny = '+str(OptimalEstimation1.CHISQ))

        #Does the trial solution fit the data better?
        if (OptimalEstimation1.PHI <= OPHI):
            print('Successful iteration. Updating xn,yn and kk')
            OptimalEstimation.edit_XN(XN1)
            OptimalEstimation.edit_YN(YN1)
            OptimalEstimation.edit_KK(KK1)
            Variables.edit_XN(XN1)

            #Now calculate the gain matrix and averaging kernels
            OptimalEstimation.calc_gain_matrix()

            #Updating the cost function
            OptimalEstimation.calc_phiret()

            #Has the solution converged?
            tphi = 100.0*(OPHI-OptimalEstimation.PHI)/OPHI
            if (tphi>=0.0 and tphi<=OptimalEstimation.PHILIMIT and alambda<1.0):
                print('phi, phlimit : '+str(tphi)+','+str(OptimalEstimation.PHILIMIT))
                print('Phi has converged')
                print('Terminating retrieval')
                break
            else:
                OPHI=OptimalEstimation.PHI
                alambda = alambda*0.3  #reduce Marquardt brake

        else:
            #Leave xn and kk alone and try again with more braking
            alambda = alambda*10.0  #increase Marquardt brake


    #Calculating output parameters
    ######################################################

    #Calculating retrieved covariance matrices
    OptimalEstimation.calc_serr()

    #Make sure errors stay as a priori for kiter < 0
    if OptimalEstimation.NITER<0:
        OptimalEstimation.ST = deepcopy(OptimalEstimation.SA)

    #Closing .itr file
    if write_itr==True:
        if OptimalEstimation.NITER>0:
            fitr.close()

    #Writing the contribution of each gas to .gcn file
    #if nemesisSO==True:
    #    calc_gascn(runname,Variables,Measurement,Atmosphere,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer)

    return OptimalEstimation


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







