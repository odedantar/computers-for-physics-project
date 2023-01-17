import numpy as np
import pandas as pd

from source import parse_data
from source import physics_engine as engine


def simulation_step(motion_data, time_delta):
    """
    Computes the positions of objects after one time_delta
    :param motion_data: dataframe containing stellar objects' details at a certain time
    :param time_delta: time in seconds to advance positions
    :return: dataframe containing updated data for the objects
    """
    step = motion_data['step'] + 1
    masses = motion_data['mass']
    initial_location = motion_data.loc[:, ['x', 'y', 'z']]
    initial_velocity = motion_data.loc[:, ['vx', 'vy', 'vz']]

    forces = engine.get_forces(motion_data)
    acceleration = engine.get_accelerations(forces, masses)
    velocities = engine.get_velocities(initial_velocity, acceleration, time_delta)
    locations = engine.get_locations(initial_location, velocities, time_delta)

    update = pd.concat([locations, velocities, step, masses], axis=1)

    return update


def simulate(steps, time_delta, initial_data=None, path='./data/data.csv'):
    """
    Performs a specified amount of manager steps with time_delta seconds between steps.
    :param steps: number of steps to simulate
    :param time_delta: time in seconds to advance positions between steps
    :param initial_data: dataframe containing stellar objects' details at a certain time
    :param path: string containing data's file path - in case the function is required to load
    :return: a dataframe containing all objects' motion data at all times
    """
    if initial_data is None:
        initial_data = parse_data.open_data_csv(path)

    motion_list = []
    motion_data = simulation_step(initial_data, time_delta)
    motion_list.extend([motion_data])

    for s in range(2, steps + 1):
        motion_data = simulation_step(motion_data, time_delta)
        motion_list.extend([motion_data])

    motion_sequence = pd.concat(motion_list)
    motion_sequence.set_index(['step'], append=True, inplace=True)
    return motion_sequence


def get_kinetic_energy(motion_sequence):
    """
    Computes the kinetic energy of the system as a function of manager step
    :param motion_sequence: a DataFrame containing planet data for many steps
    :return: the kinetic energy at each step (dataframe)
    """
    energy_sequence = engine.get_kinetic_sequence(motion_sequence['mass'].xs(1, level='step'),
                                                  motion_sequence[['vx', 'vy', 'vz']])
    steps = len(energy_sequence.index.unique('step'))
    energies = []

    for s in range(1, steps + 1):
        energy = pd.DataFrame([np.sum(energy_sequence.xs(s, level='step').to_numpy())],
                              columns=['kinetic'])
        energy['step'] = s
        energies.extend([energy])

    energies = pd.concat(energies)
    energies.set_index(['step'], inplace=True)

    return energies


def get_potential_energy(motion_sequence):
    """
    Computes the gravitational potential energy of the system as a function of manager step
    :param motion_sequence: a DataFrame containing planet data for many steps
    :return: the potential energy at each step (dataframe)
    """
    energy_sequence = engine.get_potential_sequence(motion_sequence['mass'].xs(1, level='step'),
                                                    motion_sequence[['x', 'y', 'z']])
    steps = len(energy_sequence.index.unique('step'))
    energies = []

    for s in range(1, steps + 1):
        energy = pd.DataFrame([np.sum(energy_sequence.xs(s, level='step').to_numpy())],
                              columns=['potential'])
        energy['step'] = s
        energies.extend([energy])

    energies = pd.concat(energies)
    energies.set_index(['step'], inplace=True)

    return energies


def get_total_energy(motion_sequence):
    """
    Computes the total energy of the system as a function of manager step
    :param motion_sequence: a DataFrame containing planet data for many steps
    :return: the total energy at each step (dataframe)
    """
    kinetic = get_kinetic_energy(motion_sequence)
    potential = get_potential_energy(motion_sequence)
    total = pd.DataFrame(kinetic.to_numpy() + potential.to_numpy())

    total.index = kinetic.index
    total.columns = ['energy']

    return total


def get_sun_distances(locations_sequence, planet_name):
    """
    Computes the distance of a planet from the sun as a function of time
    :param locations_sequence: object data returned from simulate
    :param planet_name: the name of the planet to find distances from the sun (a string)
    :return: list/array of distances of the planet from the sun (dataframe)
    """
    distances = engine.get_distance_sequence(locations_sequence.loc[['sun', planet_name]])

    return distances.loc['sun'][planet_name]


def get_orbit_periods(sun_distances):
    """
    Estimates the orbital period (in manager steps) of all planets using their distances from the sun
    :param sun_distances: a list of distances from the sun (as a function of time)
    :return: the estimated period of the motion
    """
    periods = []
    planets = sun_distances.columns

    for planet in planets:
        freq = engine.get_frequency(sun_distances[planet])

        T = pd.DataFrame({'frequency': [freq],
                          'period': [1 / freq],
                          'name': [planet]})

        periods.append(T)

    periods = pd.concat(periods)
    periods.set_index(['name'], inplace=True)

    return periods


if __name__ == "__main__":
    pass
