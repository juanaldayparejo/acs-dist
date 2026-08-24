#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# acs - Python package to process observations from TGO/ACS
# database - Functions to make a database of the ACS MIR observations
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
import pandas as pd
from pathlib import Path
import acs

acs_dir = Path(acs.__file__).resolve().parent
database_file = (acs_dir.parent / "data" / "mir_database.csv").resolve()

###########################################################################################################################

def initialise_database():
    """
        FUNCTION NAME : initialise_database()
        
        DESCRIPTION : Initialise the ACS MIR database
        
        INPUTS : None
        
        OPTIONAL INPUTS: None
        
        OUTPUTS :
        
            database :: Pandas DataFrame containing the database
        
        CALLING SEQUENCE:
        
            database = initialise_database()
        
        MODIFICATION HISTORY : Juan Alday (25/02/2025)
        
    """
    
    database = pd.DataFrame(columns=["orbit", "num", "date", "time", "position", "mtp", "stp", "IEflag", "FPflag", "Ls", "LST", "Latitude", "Longitude", "MY"])  # Empty DataFrame

    return database


###########################################################################################################################

def read_database(csvname=database_file):
    """
        FUNCTION NAME : read_database()
        
        DESCRIPTION : Reads the ACS MIR database from a CSV file
        
        INPUTS : None
        
        OPTIONAL INPUTS:
        
            csvname :: Name of the CSV file containing the database (default: "database.csv")
        
        OUTPUTS :
        
            database :: Pandas DataFrame containing the database
        
        CALLING SEQUENCE:
        
            database = read_database()
        
        MODIFICATION HISTORY : Juan Alday (25/02/2025)
        
    """

    try:
        df = pd.read_csv(csvname)  # Load from file if it exists
        return df
    except FileNotFoundError:
        print("File not found")
    
###########################################################################################################################
    
# Function to validate IEflag and FPflag
def validate_flags(IEflag, FPflag):
    """
        FUNCTION NAME : validate_flags()
        
        DESCRIPTION : Ensure flags are consistent
        
        INPUTS : 
        
            IEflag :: Ingress/Egress flag (I or E)
            FPflag :: Full/Partial frame flag (F or P)
        
        OPTIONAL INPUTS: None
        
        OUTPUTS :
        
            database :: Pandas DataFrame containing the database
        
        CALLING SEQUENCE:
        
            validate_flags(IEflag, FPflag)
        
        MODIFICATION HISTORY : Juan Alday (25/02/2025)
    """
    if IEflag not in {"I", "E", "C"}:
        raise ValueError(f"Invalid IEflag '{IEflag}'. Must be 'I' or 'E' or 'C'.")
    if FPflag not in {"F", "P", "D"}:
        raise ValueError(f"Invalid FPflag '{FPflag}'. Must be 'F' or 'P' or 'D'.")

###########################################################################################################################

def add_record(database, orbit, num, date, time, mtp, stp, position, IEflag, FPflag, Ls, LST, Latitude, Longitude, MY):
    """
        FUNCTION NAME : add_record()
        
        DESCRIPTION : Add record to the database
        
        INPUTS : 
        
            database :: Pandas DataFrame containing the database
            orbit :: Orbit number
            num :: Observation number within orbit
            date :: Observation date
            time :: Observation time
            mtp :: MTP number
            stp :: STP number
            position :: Secondary grating position
            IEflag :: Ingress/Egress flag (I or E)
            FPflag :: Full/Partial frame flag (F or P)
            Ls :: Solar longitude
            LST :: Local solar time
            Latitude :: Latitude
            Longitude :: Longitude
            MY :: Mars Year
            
        OPTIONAL INPUTS: None
        
        OUTPUTS :
        
            database :: Updated Pandas DataFrame containing the database
        
        CALLING SEQUENCE:
        
            database = add_record(database, orbit, num, date, time, mtp, stp, position, IEflag, FPflag, Ls, LST, Latitude, Longitude, MY)
        
        MODIFICATION HISTORY : Juan Alday (25/02/2025)
    """

    # Validate flags
    #validate_flags(IEflag, FPflag)

    # Convert to integer for correct ordering
    orbit = int(orbit)
    if isinstance(num, str):
        num = int(num[1])

    # New record dictionary
    new_record = {
        "orbit": orbit, "num": num, "date": date, "time": time, "mtp": mtp, "stp": stp,
        "position": position, "IEflag": IEflag, "FPflag": FPflag, "Ls": Ls, "LST": LST,
        "Latitude": Latitude, "Longitude": Longitude, "MY": MY
    }

    # Check if (orbit, num) already exists
    mask = (database["orbit"] == orbit) & (database["num"] == num)

    if mask.any():
        # Update existing record
        database.loc[mask, new_record.keys()] = new_record.values()
        print(f"Updated existing record for orbit {orbit}, num {num}.")
    else:
        # Append new record
        database = pd.concat([database, pd.DataFrame([new_record])], ignore_index=True)
        print(f"Added new record for orbit {orbit}, num {num}.")

    # Sort by orbit and num
    database = database.sort_values(by=["orbit", "num"]).reset_index(drop=True)

    return database

###########################################################################################################################

def write_database(database,csvname=database_file):
    """
        FUNCTION NAME : read_database()
        
        DESCRIPTION : Reads the ACS MIR database from a CSV file
        
        INPUTS :
        
            database :: Pandas DataFrame containing the database
        
        OPTIONAL INPUTS:
        
            csvname :: Name of the CSV file (default: "database.csv")
        
        OUTPUTS :
        
            CSV file containing the database
        
        CALLING SEQUENCE:
        
            write_database()
        
        MODIFICATION HISTORY : Juan Alday (25/02/2025)
        
    """

    database.to_csv(csvname, index=False)
    
###########################################################################################################################

# Function to search for a record based on orbit and num
def search_records(database, orbits, nums):
    """
        FUNCTION NAME : search_record()
        
        DESCRIPTION : Search for a record on database based on orbit and num
        
        INPUTS :
        
            database :: Pandas DataFrame containing the database
            orbits :: Orbit numbers
            nums :: Observation numbers within orbit
        
        OPTIONAL INPUTS: None
        
        OUTPUTS :
        
            record :: Dictionary containing the record
        
        CALLING SEQUENCE:
        
            record = search_record(database, orbit, num)
        
        MODIFICATION HISTORY : Juan Alday (25/02/2025)
        
    """

    # Ensure input is a list (even if it's a single value)
    if isinstance(orbits, int):
        orbits = [orbits]
    if isinstance(nums, int):
        nums = [nums]
    
    # Ensure the two lists have the same length
    if len(orbits) != len(nums):
        raise ValueError("The 'orbits' and 'nums' lists must have the same length.")

    # Create a mask to filter the DataFrame
    mask = database.apply(lambda row: (row["orbit"], row["num"]) in zip(orbits, nums), axis=1)

    # Get matching records
    matching_records = database[mask]

    if matching_records.empty:
        print("No matching records found.")
        return None
    else:
        #print(f"Found {len(matching_records)} matching records.")
        return matching_records

###########################################################################################################################

def read_geometry_essentials(filename,altx=30.):

    """
    FUNCTION NAME : read_geometry_essentials()

    DESCRIPTION : This function reads the essential geometry parameters to add to the database

    INPUTS : 

        filename :: Name of the file

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
          
        orb :: Orbit
        nflag:: Observation number within orbit
        date :: Date of observation
        time :: Time of observation
        IEobs :: Ingress/Egress flag
        pos :: Secondary grating position
        FPobs :: Full/Partial frame flag
        Ls :: Solar longitude
        LST :: Local solar time
        lat :: Latitude 
        lon :: Longitude

    CALLING SEQUENCE:

        orb,nflag,date,time,IEobs,pos,FPobs,Ls,Loct,lat,lon = read_geometry_essentials(filename)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    file1 = filename
    orb = file1[16:22]
    nflag = file1[23:25]
    IEobs = file1[26:27]
    pos = file1[31:33]
    FPobs = file1[34:35]
    
    #Getting number of lines in file
    with open(file1, "r", encoding="utf-8") as file:
        line_count = sum(1 for _ in file)

    #Reading the file
    nacq = line_count - 3
    f = open(file1,'r')
    dummy1 = f.readline().split()
    dummy2 = f.readline().split()
    dummy3 = f.readline().split()
    
    Lsx = np.zeros(nacq)
    Loctx = np.zeros(nacq)
    latx = np.zeros(nacq)
    lonx = np.zeros(nacq)
    tanheax = np.zeros(nacq)
    tanheex = np.zeros(nacq)
    datex = ['']*nacq
    timex = ['']*nacq
    for i in range(nacq):
        tmp = f.readline().split()
        datetime = tmp[0]
        datex[i] = datetime[0:10]
        timex[i] = datetime[11:19]
        Lsx[i] = float(tmp[2])
        Loct1 = tmp[8]
        hour = float(Loct1[0:2])
        minute = float(Loct1[3:5])
        second = float(Loct1[6:8])
        Loctx[i] = hour + minute/60. + second/3600.
        latx[i] = float(tmp[12])
        lonx[i] = float(tmp[13])
        tanheax[i] = float(tmp[20])
        tanheex[i] = float(tmp[21])
    
    f.close()
    
    #Finding the value that corresponds to a given altitude
    ialt = np.argmin( np.abs(tanheax-altx) )
    
    Ls = Lsx[ialt]
    Loct = Loctx[ialt]
    lat = latx[ialt]
    lon = lonx[ialt]
    date = datex[ialt]
    time = timex[ialt]
                    
    return orb,nflag,date,time,IEobs,pos,FPobs,Ls,Loct,lat,lon

###########################################################################################################################