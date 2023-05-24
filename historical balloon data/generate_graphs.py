# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Ellipse

from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull

# %%


def distance_histogram(df):

    # Generate some random data
    data1 = df["burst_distance"]
    data2 = df["land_distance"]

    # Create a figure with two rows of histograms that share the same x-axis
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(6, 6))

    # Create histograms for the first row
    bin_edges = np.arange(0, 500, 25)

    axs[0].hist(data1, bins=bin_edges, color='blue', alpha=0.5)
    axs[1].hist(data2, bins=bin_edges, color='red', alpha=0.5)

    # # Create histograms for the second row
    # axs[1].hist(data1, bins=20, color='green', alpha=0.5)
    # axs[1].hist(data2, bins=20, color='orange', alpha=0.5)

    # Add titles and labels for the subplots and the overall figure
    axs[0].set_title('Distance traveled until burst')
    axs[1].set_title('Distance traveled to landing')
    axs[1].set_xlabel('Distance [km]')
    axs[0].set_ylabel('Frequency')
    axs[1].set_ylabel('Frequency')

    # Adjust the spacing between the subplots
    fig.subplots_adjust(hspace=0.4)
    fig.tight_layout()

    # Display the figure
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig("figs/distance_traveled.png")

# %%


def direction_polar_plot(df):
    # Convert directions from degrees to radians
    angles = np.deg2rad(df["bearing_burst"])

    # Create a polar histogram of the directions
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.hist(angles, bins=36)

    # Set the direction of the radial axis to East
    ax.set_theta_zero_location('N')

    # Set the direction of the theta axis to clockwise
    ax.set_theta_direction(-1)

    # Set the title of the plot
    ax.set_title('Polar Histogram of Drift Directions')

    # Display the plot
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig("figs/polar_drift_direction.png")

# %%


def plot_locations_map(df):

    # Set the projection of the map
    projection = ccrs.PlateCarree()

    # Create a figure and axes for the map
    fig, ax = plt.subplots(subplot_kw=dict(projection=projection))

    # Set the extent of the map to the Netherlands
    extent = [3.2, 9.6, 50.4, 53.7]
    ax.set_extent(extent, crs=projection)

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
    ax.add_feature(cfeature.COASTLINE)

    # Plot the coordinates on the map
    lons = df["burst_position_lon"]
    lats = df["burst_position_lat"]
    ax.scatter(lons, lats, color='red', transform=projection, alpha=0.5, s=7)

    d = df.iloc[0]
    ax.plot(d["launch_position_lon"], d["launch_position_lat"],
            marker='o', markersize=10, color='green', transform=projection)

    # Add a title to the map
    ax.set_title('plotted Burst locations from balloons from the bilt')

    plt.tight_layout()

    # Display the map
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig("figs/burst_locations.png")

    # Create a figure and axes for the map
    fig, ax = plt.subplots(subplot_kw=dict(projection=projection))

    # Set the extent of the map to the Netherlands
    extent = [3.2, 9.6, 50.4, 53.7]
    ax.set_extent(extent, crs=projection)

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
    ax.add_feature(cfeature.COASTLINE)

    # Plot the coordinates on the map
    lons = df["land_position_lon"]
    lats = df["land_position_lat"]
    ax.scatter(lons, lats, color='red', transform=projection, alpha=0.5, s=7)

    d = df.iloc[0]
    ax.plot(d["launch_position_lon"], d["launch_position_lat"],
            marker='o', markersize=10, color='green', transform=projection)
    # Add a title to the map
    ax.set_title('plotted Burst locations from balloons from the bilt')

    plt.tight_layout()

    # Display the map
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig("figs/land_locations.png")

    # plt.show()

# %%


def calculate_elipse(lats, lons, proj, n_std=1):
    # Convert the coordinates to map projection coordinates
    x, y, _ = proj.transform_points(ccrs.Geodetic(), lons, lats).T

    # Calculate the average location of the points
    average_location = np.mean(np.column_stack((x, y)), axis=0)

    # Calculate the distances between the average location and all the points
    distances_to_average = np.sqrt(
        np.sum((np.column_stack((x, y)) - average_location)**2, axis=1))

    # Sort the distances in ascending order
    sorted_distances = np.sort(distances_to_average)

    # Calculate the covariance matrix of the coordinates
    covariance_matrix = np.cov(np.column_stack((x, y)).T)

    # Get the eigenvalues and eigenvectors of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)

    # Sort the eigenvalues in descending order
    sorted_indices = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[sorted_indices]
    eigenvectors = eigenvectors[:, sorted_indices]

    # Calculate the semi-axes of the ellipse
    semi_major_axis = np.sqrt(eigenvalues[0])*n_std
    semi_minor_axis = np.sqrt(eigenvalues[1])*n_std

    # Calculate the angle of rotation of the ellipse
    angle = np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])

    return semi_major_axis, semi_minor_axis, angle, average_location


def plot_average_elipse(df, n_std=2):
    # Set the projection of the map
    projection = ccrs.PlateCarree()

    # Create a figure and axes for the map
    fig, ax = plt.subplots(subplot_kw=dict(projection=projection))

    # Set the extent of the map to the Netherlands
    extent = [3.2, 9.6, 50.4, 53.7]
    ax.set_extent(extent, crs=projection)

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
    ax.add_feature(cfeature.COASTLINE)

    # Plot the coordinates on the map
    lons = df["burst_position_lon"]
    lats = df["burst_position_lat"]
    ax.scatter(lons, lats, color='red', transform=projection, alpha=0.2, s=5)

    # actaully calculate the elipse

    semi_major_axis, semi_minor_axis, angle, average_location = calculate_elipse(
        lats, lons, projection, n_std=n_std)

    # Draw the ellipse on the map
    ellipse = Ellipse(xy=average_location, width=2*semi_major_axis, height=2*semi_minor_axis,
                      angle=angle*180/np.pi, facecolor='none', edgecolor='green', transform=projection)
    ax.add_artist(ellipse)

    d = df.iloc[0]
    ax.plot(d["launch_position_lon"], d["launch_position_lat"],
            marker='o', markersize=10, color='green', transform=projection)

    ax.set_title(f"burst locations with {n_std}std confidence elipse")

    plt.tight_layout()

    plt.rcParams['savefig.dpi'] = 300
    plt.savefig(f"figs/{n_std}_elipse.png")

    # PLOT the landing location
    # Set the projection of the map
    projection = ccrs.PlateCarree()

    # Create a figure and axes for the map
    fig, ax = plt.subplots(subplot_kw=dict(projection=projection))

    # Set the extent of the map to the Netherlands
    extent = [3.2, 9.6, 50.4, 53.7]
    ax.set_extent(extent, crs=projection)

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
    ax.add_feature(cfeature.COASTLINE)

    # Plot the coordinates on the map
    lons = df["land_position_lon"]
    lats = df["land_position_lat"]
    ax.scatter(lons, lats, color='red', transform=projection, alpha=0.2, s=5)

    # actaully calculate the elipse

    semi_major_axis, semi_minor_axis, angle, average_location = calculate_elipse(
        lats, lons, projection, n_std=n_std)

    # Draw the ellipse on the map
    ellipse = Ellipse(xy=average_location, width=2*semi_major_axis, height=2*semi_minor_axis,
                      angle=angle*180/np.pi, facecolor='none', edgecolor='green', transform=projection)
    ax.add_artist(ellipse)

    d = df.iloc[0]
    ax.plot(d["launch_position_lon"], d["launch_position_lat"],
            marker='o', markersize=10, color='green', transform=projection)

    ax.set_title(f"landing locations with {n_std}std confidence elipse")

    plt.tight_layout()

    plt.rcParams['savefig.dpi'] = 300
    plt.savefig(f"figs/{n_std}_elipse_landing.png")
    # plt.show()

# %%


def plot_elipse_reduced_drift(df, percentile=95, n_std=2):
    # Set the projection of the map
    projection = ccrs.PlateCarree()

    # Create a figure and axes for the map
    fig, ax = plt.subplots(subplot_kw=dict(projection=projection))

    # Set the extent of the map to the Netherlands
    extent = [3.2, 9.6, 50.4, 53.7]
    ax.set_extent(extent, crs=projection)

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
    ax.add_feature(cfeature.COASTLINE)

    # Plot the coordinates on the map
    percentile_drift = np.percentile(
        df["average_drift_velocity_ascend"], percentile)
    reduced_df = df[df["average_drift_velocity_ascend"] < percentile_drift]
    lons = reduced_df["burst_position_lon"]
    lats = reduced_df["burst_position_lat"]

    ax.scatter(lons, lats, color='red', transform=projection, alpha=0.2, s=5)

    # actaully calculate the elipse

    semi_major_axis, semi_minor_axis, angle, average_location = calculate_elipse(
        lats, lons, projection, n_std=n_std)

    # Draw the ellipse on the map
    ellipse = Ellipse(xy=average_location, width=2*semi_major_axis, height=2*semi_minor_axis,
                      angle=angle*180/np.pi, facecolor='none', edgecolor='green', transform=projection)
    ax.add_artist(ellipse)

    d = df.iloc[0]
    ax.plot(d["launch_position_lon"], d["launch_position_lat"],
            marker='o', markersize=10, color='green', transform=projection)

    ax.set_title(
        f"{n_std} std when ommiting drift velocities above {percentile} percentile")

    plt.tight_layout()

    plt.rcParams['savefig.dpi'] = 300
    plt.savefig(f"figs/{n_std}_{percentile}_omitted_elipse.png")
    # plt.show()

# %%


def plot_average_drift_velocity(df):
    # Convert directions from degrees to radians
    vel = df["average_drift_velocity_ascend"]

    # Create a polar histogram of the directions
    fig, ax = plt.subplots()
    ax.hist(vel, bins="auto")

    # Set the direction of the radial axis to East
    # ax.set_theta_zero_location('N')

    # Set the direction of the theta axis to clockwise
    # ax.set_theta_direction(-1)
    ax.set_xlabel("velocity [m/s]")
    ax.set_ylabel("frequency")

    # Set the title of the plot
    ax.set_title('Histogram of Drift velocity on ascend')

    plt.tight_layout()
    # Display the plot
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig("figs/hist_drift_velocity.png")


# %%
def load_data():
    return pd.read_excel("Bilt_balloon_data.xlsx")

# %%


def main():
    df = load_data()
    distance_histogram(df)
    direction_polar_plot(df)
    plot_locations_map(df)
    plot_average_elipse(df)
    plot_elipse_reduced_drift(df)
    plot_average_drift_velocity(df)

    # plt.show()


if __name__ == "__main__":
    main()
