import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from source import simulation as sim


def plot_trajectory(motion_sequence):
    """
    Plots the trajectory and current position of objects, from a dataframe containing positions at multiple times
    :param motion_sequence: object data returned from simulate
    """
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    planets = motion_sequence.index.unique('name')

    for p in planets:
        xdata = motion_sequence.loc[p, :]['x']
        ydata = motion_sequence.loc[p, :]['y']
        zdata = motion_sequence.loc[p, :]['z']

        ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='viridis')


def plot_sun_distances(sequence, preload=False):
    """
    Plots the distance of planets from the sun as a function of time
    :param sequence: object data returned from simulate
    :param preload: flag indicating if the distances are preloaded and given in the sequence
    """
    if preload:
        planets = sequence.columns.unique()

        for i, planet in enumerate(planets):
            distance = sequence[planet]
            plt.plot(distance, label=planet)
    else:
        planets = sequence.index.unique('name').drop('sun')

        for i, planet in enumerate(planets):
            distance = sim.get_sun_distances(sequence[['x', 'y', 'z']], planet)
            plt.plot(distance, label=planet)

    plt.legend(loc="upper left")


def plot_kepler_third_law(sequence, preload=False):
    """
    Plots the (orbital period)^2 of planet's motion as a function of (mean distance to the sun)^3
    :param sequence: object data returned from simulate
    :param preload: flag indicating if the distances are preloaded and given in the sequence
    """
    planets = []
    mean_distances = []
    distance_sequence = []

    # Preloading or generating distances data:
    if preload:
        planets = sequence.columns.unique()
        distance_sequences = sequence

    else:
        planets = sequence.index.unique('name').drop('sun')
        distances = []

        for i, planet in enumerate(planets):
            distance = sim.get_sun_distances(sequence[['x', 'y', 'z']], planet)
            distances.extend([distance])

        distance_sequences = pd.concat(distances, axis=1)

    # Calculating (distance)^3 for each planet
    for planet in planets:
        dist_seq = distance_sequences[planet].to_numpy()
        distance = pd.DataFrame({'name': [planet],
                                 'distance': [np.average(dist_seq) ** 3]})
        mean_distances.append(distance)

    mean_distances = pd.concat(mean_distances)
    mean_distances.set_index(['name'], inplace=True)

    # Calculating (period)^2 for each planet
    periods = sim.get_orbit_periods(distance_sequences)['period'] ** 2

    plt.scatter(mean_distances, periods)


def plot_energy(motion_sequence):
    """
    Plots kinetic, potential and total energy of the system as a function of time
    :param motion_sequence: object data returned from simulate
    """
    energy = sim.get_total_energy(motion_sequence)

    plt.plot(energy, label='total energy')
    plt.legend(loc="upper left")


def show():
    plt.show(block=True)


if __name__ == "__main__":
    pass

