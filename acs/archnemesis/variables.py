#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# archnemesis.variables - Functions to initialise the model parameterisations that will be used in the retrieval
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



###############################################################################################

def write_apr_prof(filename,Atmosphere,varID,var2ID,error,clen=1.5,ScaleFactor=True,xprof=None):

    """

    FUNCTION NAME : write_apr_prof()

    DESCRIPTION : Write the file for retrieving a continuous profile of gas VMR, temperature or dust density

    INPUTS : 

        filename :: Name of the file (must be the same as indicated in .apr file)
        Atmosphere :: Python class defining the atmosphere
        varID :: Variable to be retrieved
                  - varID = 0 indicates that temperature is to be retrieved
                  - varID > 0 indicates that a gas vmr profile is to be retrieved. varID
                              must then be equal to the gasID of the gas that wants to be retrieved
                  - varID < 0 indicates that an aerosol density profile is to be retrieved.
                              varID must then indicate the aerosol population number that
                              wants to be retrieved (e.g. -2 indicates that the second 
                              population of aerosols in the aerosol.ref wants to be retrieved)

        var2ID :: Isotopologue to be retrieved. If varID = 0 or varID < 0 then this variables,
                  must be set to 0. If a gas vmr profile is to be retrieved, then var2ID must
                  be equal to the isotopologue ID of the gas that wants to be retrieved. 
        err :: A priori uncertainty on the profile to be retrieved, expressed as an absolute
               magnitude in the same units as the profile is defined. 

    OPTIONAL INPUTS:
    
        clen :: Correlation length, expressed as a factor of scale height
        ScaleFactor :: If True, then error represents a fractional error rather than absolute
        xprof(npro) :: If we do not want to use the profile in the Atmosphere as the initial file
                        we can define xprof to use it as the a priori profile (default is None)
            
    OUTPUTS : 
 
        A priori file

    CALLING SEQUENCE:

        write_apr_prof(filename,Atmosphere,gasID,isoID)

    MODIFICATION HISTORY : Juan Alday (06/05/2023)

    """
    
    npro = Atmosphere.NP
    
    #Temperature
    if varID==0:
        xref = np.zeros(Atmosphere.NP)
        if xprof is None:
            xref[:] = Atmosphere.T
        else:
            if len(xprof)!=npro:
                raise ValueError('error :: xprof must have the same length as NPRO')
            else:
                xref[:] = xprof
    
    #Dust    
    elif varID<0:
        caero = abs(varID)
        if caero>Atmosphere.NDUST:
            raise ValueError('error :: The aerosol population that wants to be retrieved does not exist in aerosol.ref')
        
        xref = np.zeros(npro)
        if xprof is None:
            xref[0:npro] = Atmosphere.DUST[0:npro,caero-1]
        else:
            if len(xprof)!=npro:
                raise ValueError('error :: xprof must have the same length as NPRO')
            else:
                xref[:] = xprof
                
    #Gas abundance      
    elif varID>0:
        cgas1 = np.where((Atmosphere.ID==varID) & (Atmosphere.ISO==var2ID))
        cgas = cgas1[0]

        xref = np.zeros(npro)
        if xprof is None:
            xref[0:npro] = Atmosphere.VMR[0:npro,cgas[0]]
        else:
            if len(xprof)!=npro:
                raise ValueError('error :: xprof must have the same length as NPRO')
            else:
                xref[:] = xprof
                
    #Creating error array
    errarr = np.zeros(npro)
    errarr[0:npro] = error

    if ScaleFactor==True:
        errarr = xref*errarr

    #Write file
    fref = open(filename,'w')
    fref.write('\t %i \t %10.3f \n' % (npro,clen))
    for i in range(npro):
        fref.write('\t %10.6e \t %10.6e \t %10.6e \n' % (Atmosphere.P[i],xref[i],errarr[i]))
    fref.close()


###############################################################################################

def create_apr_file(runname,
        Atmosphere, 
        Measurement,
        retrieve_temp=False,
        retrieve_press=False,
        retrieve_continuous_prof=False,
        retrieve_baseline=False,
        flag_temp_analytic=True,
        varID1=None,varID2=None,
        cont_clen=None,cont_err=None,
        temp_clen=1.5,temp_err=None,
        htan=40.,ptan=None,ptanerr=0.1,
        baseline_degree=2):

    """
    FUNCTION NAME : create_apr_file()

    DESCRIPTION : Function to create the archNEMESIS .apr file to define the retrieved parameters

    INPUTS : 

        runname :: Name of the archNEMESIS run
        Atmosphere :: Atmosphere class
        Measurement :: Measurement class

    OPTIONAL INPUTS:
    
        retrieve_temp :: Flag to retrieve temperature profile with numerical calculation of Jacobian (i.e., hydrostatic approach)
            flag_temp_analytic :: Flag indicating whether the temperature retrieval approach must be from hydrostatic (False) or spectroscopy (True)
            temp_clen :: Correlation length expressed as a fraction of the scale height (required if retrieve_temp=True)
            temp :: A priori uncertainty in the temperature (K) (required if retrieve_temp=True)
        retrieve_press :: Flag to retrieve the pressure at a given tangent height
            htan :: Tangent height at which the pressure must be retrieved in km (required if retrieve_press=True)
            ptan :: Pressure at the given tangent height (optional if retrieve_press=True. If None then it will be taken from the Atmosphere class)
            ptanerr :: Fractional error in the a priori pressure (required if retrieve_press=True)
        retrieve_continuous_prof :: Flag to retrieve a continuous profile of gas VMR or dust number density
            varID1(ngas) :: ID of the gases or aerosols to be retrieved (see archnemesis documentation. required if retrieve_continuous_prof=True)
            varID2(ngas) :: ID of the isotopes to be retrieved (see archnemesis documentation. required if retrieve_continuous_prof=True)
            cont_clen :: Correlation length expressed as a fraction of the scale height (required if retrieve_continuous_prof=True)
            cont_err(ngas) :: Fractional error in the gas vmrs or dust densities (required if retrieve_continuous_prof=True)
        retrieve_baseline :: Flag to retrieve the baseline at each tangent height with a polynomial function
            baseline_degree :: Degree of the polynomial that must be used to fit the baseline (required if retrieve_baseline=True)

    OUTPUTS : 
 
        archNEMESIS .apr file

    MODIFICATION HISTORY : Juan Alday (29/07/2026)

    """

    #Counting number of retrieved parameters
    ####################################################################################################

    nvar = 0

    #Counting number of variables from continuous profiles
    if retrieve_continuous_prof is True:
        nvar += len(varID1)

    #Counting other variables
    nvar += int(retrieve_temp) + int(retrieve_press) + int(retrieve_baseline)

    #Writing the file
    ####################################################################################################

    fapr = open(runname+'.apr','w')
    fapr.write('#ACS retrieval \n')
    fapr.write('\t'+str(nvar)+' \n')

    #Going through the continuous profiles
    if retrieve_continuous_prof is True:

        ngas = len(varID1)
        if len(varID1) != len(varID2):
            raise ValueError("error while writing continuous profiles in .apr file :: varID1 and varID2 must be of the same length")

        if cont_err is None:
            raise ValueError("error while writing continuous profiles in .apr file :: cont_err must be defined for each of the retrieved species")

        if cont_clen is None:
            cont_clen = np.ones(ngas) * 1.5 

        for igas in range(ngas):

            gasid = varID1[igas] ; isoid = varID2[igas]

            fapr.write('\t %i \t %i \t %i \n' % (gasid,isoid,0))
            filename = 'contprof'+str(igas)+'.dat'
            fapr.write(filename+' \n')

            write_apr_prof(filename,Atmosphere,gasid,isoid,cont_err[igas],clen=cont_clen[igas],ScaleFactor=True)

    #Going through the retrieval of the baseline
    if retrieve_baseline is True:
        raise ValueError("error while writing baseline retrieval in .apr file :: This functionality has not yet been implemented")

    #Going through the temperature retrieval
    if retrieve_temp is True:

        if flag_temp_analytic is True:
            fapr.write('\t %i \t %i \t %i \n' % (0,0,0))
        else:
            fapr.write('\t %i \t %i \t %i \n' % (0,-1,0))
        fapr.write('tempapr.dat \n')

        filename = 'tempapr.dat'
        write_apr_prof(filename,Atmosphere,0,0,temp_err,clen=temp_clen,ScaleFactor=False)

    
    #Going through the pressure retrieval
    if retrieve_press is True:

        if htan is None:
            raise ValueError("error while writing the pressure retrieval in the .apr file :: htan must be defined if we want to retrieve pressure")


        fapr.write('\t %i \t %i \t %i \n' % (666,0,666))
        ipress0 = np.argmin(np.abs(Atmosphere.H/1.0e3-htan))
        tanhe0 = Atmosphere.H[ipress0]/1.0e3

        fapr.write('\t %10.7f \n' % (tanhe0))

        if ptan is None:
            ptan = Atmosphere.P[ipress0]/101325.

        perr = ptan * ptanerr
        fapr.write('\t %10.7e \t %10.7e \n' % (ptan,perr))


    fapr.close()
