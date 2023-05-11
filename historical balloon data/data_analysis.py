# %%
import matplotlib.pyplot as plt
import geopy.distance
import glob
import os
import json
import math
import pandas as pd
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


def analyze_stats(dat):
    data_out = dict()
    launch_position = (52.0989, 5.1797)  # De bilt launch location
    data_out["launch_position_lat"], data_out["launch_position_lon"] = launch_position
    burst_position = (float(dat[1]['lat']), float(dat[1]['lon']))
    data_out["burst_position_lat"], data_out["burst_position_lon"] = burst_position
    land_position = (float(dat[2]['lat']), float(dat[2]['lon']))
    data_out["land_position_lat"], data_out["land_position_lon"] = land_position
    data_out["launch_time"] = pd.to_datetime(dat[0]["datetime"]).tz_localize(None)
    data_out["burst_time"] = pd.to_datetime(dat[1]["datetime"]).tz_localize(None)
    data_out["landing_time"] = pd.to_datetime(dat[2]["datetime"]).tz_localize(None)
    data_out["time_to_burst"] = data_out["burst_time"] - \
        data_out["launch_time"]
    data_out["time_to_burst"] = data_out["time_to_burst"].total_seconds()
    data_out["time_to_land"] = data_out["landing_time"] - \
        data_out["launch_time"]
    data_out["time_to_land"] = data_out["time_to_land"].total_seconds()

    data_out["burst_atlitude"] = dat[1]["alt"]

    data_out["burst_distance"] = geopy.distance.distance(
        launch_position, burst_position).km
    data_out["land_distance"] = geopy.distance.distance(
        launch_position, land_position).km
    data_out["average_drift_velocity_ascend"] = data_out["burst_distance"]*1000 / \
        data_out["time_to_burst"]
    data_out["bearing_burst"] = calc_bearing(launch_position, land_position)
    df = pd.DataFrame(data_out, columns=data_out.keys(), index=[0])
 
    return df

# %%
# Define a function to convert datetime columns to timezone-unaware datetime objects
def convert_datetime(column):
    if column.dtype == 'datetime64[ns]':
        dt = pd.to_datetime(column)
        return dt.tz_localize(None)
    else:
        return column

def analyze_files():
    sonde_files = glob.glob(os.path.join("nl_data", "*.json"))
    launch_position = (52.0989, 5.1797)

    burst_distance = []
    final_distance = []

    df = pd.DataFrame(columns=['launch_position_lat', 'launch_position_lon', 'burst_position_lat', 'burst_position_lon', 'land_position_lat', 'land_position_lon',
                      'launch_time', 'burst_time', 'landing_time', 'time_to_burst', 'time_to_land', 'burst_atlitude', 'burst_distance', 'land_distance', 'average_drift_velocity_ascend'])
    # bear

    for sonde in sonde_files:
        with open(sonde, 'r') as f:
            dat = json.load(f)
            if len(dat) != 3:
                print("skipped")
                continue
            df = pd.concat([df,analyze_stats(dat)], ignore_index=True)
            # burst_distance.append(dist_burst.km)
            # final_distance.append(dist_final.km)

        

    # Apply the function to all columns in the DataFrame
    df = df.apply(convert_datetime)
    return df

#%%
def main():
    df = analyze_files()
    df.to_excel("Bilt_balloon_data.xlsx")


if __name__ == "__main__":
    main()

# %%
