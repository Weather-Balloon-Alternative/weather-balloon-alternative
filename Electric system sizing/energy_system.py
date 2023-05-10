import numpy as np
import matplotlib.pyplot as plt

# class battery_cell():
#     mass =

#     def wh_A(self, A):
#         # watthour at discharge current

#     def __init__(self) -> None:
#         pass


@np.vectorize
def liion_battery(average_power, time, max_power=None):
    """generic lithium ion battery sizing based on power requirement

    Args:
        Power (_type_): average power in watts
        time (_type_): time in seconds

    Returns:
        mass[kg]+-10%, volume[L]+-10%
    """
    wh_kg = 210  # Wh/kg, averaged from DOI: 10.1149/2.0281814jes
    wh_l = 580  # Wh/L, averaged from DOI: 10.1149/2.0281814jes
    w_kg = 600  # W/kg
    wh_eur = 3.0  # can actully be a lot cheaper

    max_power = max_power if not None else average_power
    energy = average_power*time/3600  # Wh
    mass_energy = energy/wh_kg
    mass_power = max_power/w_kg

    # np.max([mass_energy,mass_power])
    mass_actual = mass_energy if mass_energy > mass_power else mass_power

    volume = energy/wh_l

    max_energy = energy if mass_energy > mass_power else mass_power*wh_kg

    cost = max_energy/wh_eur

    return mass_energy, mass_power, mass_actual, volume, cost


# def lifepo_battery(average_power, time, max_power=None):
#     """lifepo4 lithium ion battery sizing based on power requirement

#     Args:
#         Power (_type_): average power in watts
#         time (_type_): time in seconds

#     Returns:
#         mass[kg]+-10%, volume[L]+-10%
#     """
#     max_power = max_power if not None else average_power
#     energy = average_power*time/3600  # Wh
#     mass_energy = energy/120  # Wh/kg, based on average from lifepo cells
#     mass_power = max_power/1800  # W/kg based on average from lifepo cells

#     # np.max([mass_energy,mass_power])
#     mass_actual = mass_energy if mass_energy > mass_power else mass_power

#     volume = energy/350  # Wh/L, averaged from DOI: 10.1149/2.0281814jes

#     return mass_energy, mass_power, mass_actual, volume


def first_order_energy_system_sizing():
    pass


def main():
    average_power = np.linspace(100, 4000, 100)
    mission_time = np.linspace(3600/4, 3600*5, 100)

    x = average_power
    y = mission_time
    X, Y = np.meshgrid(x, y)

    Z1, Z2, Z3, Z4, __ = liion_battery(X, Y, 4000)
    # print(Z)
    # Create a figure with 3 subplots arranged in a 1x3 grid
    fig, axs = plt.subplots(1, 4, figsize=(12, 4))

    # Plot the first subplot
    a = axs[0].pcolor(X, Y, Z1)
    fig.colorbar(a, ax=axs[0])
    # axs[0].imshow(Z1, cmap='coolwarm', extent=[x[0], x[-1], y[-1], y[0]])
    axs[0].set_title('mass for energy')
    axs[0].set_ylabel('average Power')
    axs[0].set_xlabel('flight Time')

    # Plot the second subplot
    a = axs[1].pcolor(X, Y, Z2)
    fig.colorbar(a, ax=axs[1])
    # axs[1].imshow(Z2, cmap='viridis', extent=[x[0], x[-1], y[-1], y[0]])
    axs[1].set_title('mass for power')
    axs[1].set_ylabel('average Power')
    axs[1].set_xlabel('flight Time')

    # Plot the third subplot
    a = axs[2].pcolor(X, Y, Z3)
    fig.colorbar(a, ax=axs[2])
    # axs[2].imshow(Z3, cmap='plasma', extent=[x[0], x[-1], y[-1], y[0]])
    axs[2].set_title('mass actual')
    axs[2].set_ylabel('average Power')
    axs[2].set_xlabel('flight Time')

    a = axs[3].pcolor(X, Y, Z4)
    fig.colorbar(a, ax=axs[3])
    # axs[2].imshow(Z3, cmap='plasma', extent=[x[0], x[-1], y[-1], y[0]])
    axs[3].set_title('volume')
    axs[3].set_ylabel('average Power')
    axs[3].set_xlabel('flight Time')

    # Add a colorbar for each subplot
    # for ax in axs:
    #     cbar = fig.colorbar(ax, ax=ax)
    #     cbar.ax.set_ylabel('Intensity')

    # Add a title for the entire figure
    fig.suptitle('Three 2D Color Maps')
    plt.tight_layout()
    # Display the figure
    plt.show()


if __name__ == "__main__":
    main()
