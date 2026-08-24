#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
import numpy as np
from datetime import datetime, timedelta, timezone

###########################################################################################################################
###########################################################################################################################
#                                                MARS-EARTH TIME CONVERSIONS
###########################################################################################################################
###########################################################################################################################

def julian_to_utc(julian_date, seconds_of_day):
    """
    Convert Julian Date and seconds of the day to UTC datetime.
    """
    
    JD_UNIX_EPOCH = 2440587.5  # JD corresponding to 1970-01-01 00:00:00 UTC
    unix_time = (julian_date - JD_UNIX_EPOCH) * 86400 + seconds_of_day  # Convert to seconds since epoch
    return datetime.utcfromtimestamp(unix_time)

###########################################################################################################################

def utc_to_julian(dt_list):
    """
    Convert a list of UTC datetimes or strings to (Julian Date, seconds of the day).

    Parameters:
        dt_list (list of datetime or str): A list of datetime objects in UTC or strings 
                                           in "%Y-%m-%d %H:%M:%S" format.

    Returns:
        A list of tuples [(julian_date, seconds_of_day), ...]
    """
    
    # Ensure the input is a list
    singleinput = False
    if isinstance(dt_list, str) or isinstance(dt_list, datetime):
        singleinput = True
        dt_list = [dt_list]

    JD_UNIX_EPOCH = 2440587.5  # JD corresponding to 1970-01-01 00:00:00 UTC
    
    jd = []
    sod = []

    for dt in dt_list:
        
        # Convert string to datetime if needed
        if isinstance(dt, str):
            dt = datetime.strptime(dt, "%Y-%m-%d %H:%M:%S")
        
        # Ensure UTC timezone
        dt = dt.astimezone(timezone.utc)

        # Convert datetime to Unix timestamp
        unix_time = dt.timestamp()

        # Convert Unix time to Julian Date
        jd_full = JD_UNIX_EPOCH + (unix_time / 86400.0)

        # Extract integer Julian Date and remainder
        julian_date = int(jd_full)
        seconds_of_day = (jd_full - julian_date) * 86400.0
        
        jd.append(julian_date)
        sod.append(seconds_of_day)
        
    jd = np.array(jd)
    sod = np.array(sod)

    if singleinput is True:
        return jd[0],sod[0]
    else:
        return jd,sod

###########################################################################################################################

def time_since_j2000(utc_datetimes, unit="seconds"):
    """
    Calculate time since J2000 (2000-01-01 12:00:00 UTC).
    
    Parameters:
    - utc_datetimes (numpy array of datetime or str): Array of UTC datetime values.
    - unit (str, optional): "days" (default) or "seconds" for output units.
    
    Returns:
    - numpy array of time since J2000 in chosen units.
    """
    
    #Check if input is a single value
    if isinstance(utc_datetimes, datetime):
        utc_datetimes = np.array([utc_datetimes])
    elif isinstance(utc_datetimes, str):
        utc_datetimes = np.array([datetime.strptime(utc_datetimes, "%Y-%m-%d %H:%M:%S")])
    
    #Check if the input is a string
    if isinstance(utc_datetimes[0], str):
        utc_datetimes = [datetime.strptime(dt, "%Y-%m-%d %H:%M:%S") for dt in utc_datetimes]
    
    # Define J2000 epoch
    J2000_EPOCH = datetime(2000, 1, 1, 12, 0, 0)

    # Compute time differences
    time_deltas = np.array([(dt - J2000_EPOCH) for dt in utc_datetimes])

    # Convert to desired units
    if unit == "days":
        result = np.array([td.total_seconds() / 86400 for td in time_deltas])
    elif unit == "seconds":
        result = np.array([td.total_seconds() for td in time_deltas])
    else:
        raise ValueError("Unit must be 'days' or 'seconds'")
    
    #Check if output is a single value
    if len(result)==1:
        return result[0]
    else:    
        return result

###########################################################################################################################

def getMarsParams(j2000):
    '''
    Calculate Mars Parameters of Date (Allison, M., and M. McEwen 2000; equations 16-20)
    '''

    Coefs = np.array(
    [[0.0071,2.2353,49.409],
    [0.0057,2.7543,168.173],
    [0.0039,1.1177,191.837],
    [0.0037,15.7866,21.736],
    [0.0021,2.1354,15.704],
    [0.0020,2.4694,95.528],
    [0.0018,32.8493,49.095]])

    dims = Coefs.shape
    #Mars mean anomaly:
    M = 19.3870 + 0.52402075 * j2000

    #angle of Fiction Mean Sun
    alpha = 270.3863 + 0.52403840*j2000

    #Perturbers
    PBS = 0
    for i in range(dims[0]):
        PBS += Coefs[i,0]*np.cos(((0.985626* j2000 / Coefs[i,1]) + Coefs[i,2])*np.pi/180.)

    #Equation of Center
    vMinusM = ((10.691 + 3.0e-7 *j2000)*np.sin(M*np.pi/180.) + 0.623*np.sin(2*M*np.pi/180.) + 
    0.050*np.sin(3*M*np.pi/180.) + 0.005*np.sin(4*M*np.pi/180.) + 0.0005*np.sin(5*M*np.pi/180.) + PBS)

    return M, alpha, PBS, vMinusM

###########################################################################################################################

def utc_to_marsdate(utc_datetimes):
    """
    Convert UTC datetime to Mars Solar Longitude (Ls) using the Mars Year 24 model.
    
    Parameters:
    - utc_datetimes (numpy array of datetime): Array of UTC datetime values.
    
    Returns:
    - numpy array of Mars Years (MY) corresponding to the input UTC datetimes
    - numpy array of Mars Solar Longitudes (Ls) corresponding to the input UTC datetimes.
    - numpy array of sol number in the year corresponding to the input UTC datetimes.
    
    """
    
    #Calculating the Mars year
    DPY = 686.9713  #days in a Mars year
    #refTime = [1955,4,11,10,56,0] #Mars year 1
    rDate,rSOD = utc_to_julian(['1955-04-11 10:56:00'])
    thisDate,thisSOD = utc_to_julian(utc_datetimes)
    rJD = rDate + rSOD/86400.
    thisJD = thisDate + thisSOD/86400.
    Myear = np.floor((thisJD - rJD)/DPY)+1
    Myear = np.array(Myear,dtype='int32')
    
    #Calculating the sol
    solperyear = 669.6
    sol = ( (thisJD - rJD)/DPY + 1  - Myear ) * solperyear
    
    #Check if input is a single value
    if isinstance(utc_datetimes, datetime):
        utc_datetimes = np.array([utc_datetimes])
    elif isinstance(utc_datetimes, str):
        utc_datetimes = np.array([datetime.strptime(utc_datetimes, "%Y-%m-%d %H:%M:%S")])
    
    #Check if the input is a string
    if isinstance(utc_datetimes[0], str):
        utc_datetimes = [datetime.strptime(dt, "%Y-%m-%d %H:%M:%S") for dt in utc_datetimes]
    
    # Convert UTC datetime to time since J2000 in seconds
    time_j2000 = time_since_j2000(utc_datetimes, unit="seconds")
    
    # Mars Year 24 model parameters (Allison & McEwen 2000)
    M, alpha, PBS, vMinusM = getMarsParams(time_j2000 / 86400.)
    
    # Mars Solar Longitude (Ls) in degrees
    Ls = alpha + vMinusM
    
    # Wrap Ls to [0, 360) degrees
    Ls = np.mod(Ls, 360)
    
    return Myear, Ls, sol


###########################################################################################################################

def marsdate_to_utc_sgl(MY,Ls):
    """
    Convert Mars Year (MY) and Solar Longitude (Ls) to UTC datetime.

    Parameters:
    - MY (int): Mars Year
    - Ls (float): Solar Longitude in degrees

    Returns:
    - datetime: UTC datetime corresponding to the input Mars date
    """
    
    #Calculating the Mars year
    DPY = 686.9713  #days in a Mars year
    rDate,rSOD = utc_to_julian('1955-04-11 10:56:00')
    rJD = rDate + rSOD/86400.
    
    rJD_MY = rJD+(MY-1)*DPY
    rSOD_MY = (rJD_MY - int(rJD_MY)) * 86400.
    
    utc_ls = julian_to_utc(int(rJD_MY),rSOD_MY) #LS 0 of given mars year
    
    factor = 1 #do we increment up or down?
    iter = 0  #Iteration number
    dt = 12  #hours.  This will get smaller as we get closer
    counter = 0
    olddiff = 1000.
    diff = 100
    error = 1.0e-7
    while diff > error:
		
        utc_ls = utc_ls+factor*timedelta(hours=dt)
        cMyear, cLs, csol = utc_to_marsdate(utc_ls)

        diff = np.abs(cLs - Ls)

        if diff > olddiff:
            factor = -1*factor 
            counter += 1
            if counter > 1:
                dt = dt/60.

        olddiff = diff
        #print(f'iter {iter} :: MY {cMyear} Ls {cLs:.4f} diff {diff:.6f} hours {dt:.4f}')
        iter += 1

        if iter > 10000:
            raise ValueError('error :: reached maximum number of iterations to calculate UTC from Mars date')
    
    return utc_ls

###########################################################################################################################

def marsdate_to_utc(MY,Ls):
    """
    Convert Mars Year (MY) and Solar Longitude (Ls) to UTC datetime.

    Parameters:
    - MY (int): Mars Year
    - Ls (float): Solar Longitude in degrees

    Returns:
    - datetime: UTC datetime corresponding to the input Mars date
    """
    
    #Check if input is a single value
    singlevalue = False
    if isinstance(MY, (int, np.integer)):  # np.integer is the base class for np.int32, np.int64, etc.
        singlevalue = True
        MY = np.array([MY])
        Ls = np.array([Ls])
    
    utcs = []
    for i in range(len(MY)):
        utcs.append(marsdate_to_utc_sgl(MY[i],Ls[i]))

    if singlevalue is True:
        return utcs[0]
    else:
        return utcs
