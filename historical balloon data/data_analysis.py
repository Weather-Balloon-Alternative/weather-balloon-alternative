# %%
import matplotlib.pyplot as plt
import geopy.distance
import glob
import os
import json
import math

import numpy as np

# %%
def calc_bearing(pos1, pos2):
    lat1, long1 = pos1
    lat2, long2 = pos2
    # Convert latitude and longitude to radians
    lat1 = math.radians(lat1)
    long1 = math.radians(long1)
    lat2 = math.radians(lat2)
    long2 = math.radians(long2)

    # Calculate the bearing
    bearing = math.atan2(
        math.sin(long2 - long1) * math.cos(lat2),
        math.cos(lat1) * math.sin(lat2) - math.sin(lat1) *
        math.cos(lat2) * math.cos(long2 - long1)
    )

    # Convert the bearing to degrees
    bearing = math.degrees(bearing)

    # Make sure the bearing is positive
    bearing = (bearing + 360) % 360

    return bearing

# %%
def analyze_files():
    sonde_files = glob.glob(os.path.join("nl_data", "*.json"))
    launch_position = (52.0989, 5.1797)

    burst_distance = []
    final_distance = []
    # bear

    for sonde in sonde_files:
        with open(sonde, 'r') as f:
            dat = json.load(f)
            if len(dat) != 3:
                print("skipped")
                continue
            final_pos = (dat[-1]['lat'], dat[-1]['lon'])
            burst_pos = (dat[1]['lat'], dat[1]['lon'])

            dist_burst = geopy.distance.distance(launch_position, burst_pos)
            dist_final = geopy.distance.distance(launch_position, final_pos)
            bear = calc_bearing(launch_position, final_pos)
            burst_distance.append(dist_burst.km)
            final_distance.append(dist_final.km)
    return burst_distance,final_distance

def main():
    pass
    

if __name__ == "__main__":
    main()

# %%
