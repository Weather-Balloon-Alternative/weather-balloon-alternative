#%%
import json
import glob
import math
import os.path
from math import radians, degrees, sin, cos, atan2, sqrt, pi
import numpy as np
import subprocess
import shutil
import itertools
from time import sleep

import multiprocessing as mp

#%%
def get_sonde_file_list(folder="."):
    """ Use glob to recurse through our sonde data store and return a list of all sondes files """
    return glob.glob(os.path.join(folder,"*/*/*.json"))

def load_launch_sites(filename='sites.json'):
    """
    Load in the launch sites dataset and rearrange it a bit to be useful later
    Updates to work with the new sites API structure.
    """
    _f = open(filename,'r')
    _data = _f.read()
    _f.close()

    data = json.loads(_data)

    for _station in data.keys():
        data[_station]['lat'] = float(data[_station]['position'][1])
        data[_station]['lon'] = float(data[_station]['position'][0])

    return data



#%%
def download_sonde_day(year, month, day):
    cmd = f"aws s3 cp --recursive --no-sign-request s3://sondehub-history/date/{year}/{month:02}/{day:02} sondes_{year}/{month:02}/{day:02}"
    p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT)
    
    while p.poll() is None:
        sleep(1)

    print("finished")

def download_sonde_month(year, month):
    cmd = f"aws s3 cp --recursive --no-sign-request s3://sondehub-history/date/{year}/{month:02}/ sondes_{year}/{month:02}/"
    p = subprocess.run(cmd)

#%%
def process_month(year, month):

    sondes = glob.glob(os.path.join(f"sondes_{year}/{month:02}","*/*.json"))
    for sonde in sondes:
        with open(sonde, 'r') as f:
            dat = json.load(f)
            if "launch_site" in dat[0]:
                if dat[0]["launch_site"] == "06260":
                    shutil.copy(sonde, "nl_data")

def download_process_month(year, month):
    download_sonde_month(year, month)
    process_month(year, month)



    #%%

    

# def main():
#     sites = load_launch_sites()
#     sondes = get_sonde_file_list()

# if __name__ == "__main__":
#     main()

# %%
def main():

    year_list = [ 2021]
    months = range(1,13)
    days = range(1,32)
    combinations = list(itertools.product(year_list, months, days))
    print(combinations)

    # process_month(2022, 1)
    with mp.Pool(5) as pool:
        pool.starmap(download_sonde_day, combinations)
   
        # pool.starmap(process_month, list(itertools.product(year_list, months)))

if __name__ == "__main__":
    main()
