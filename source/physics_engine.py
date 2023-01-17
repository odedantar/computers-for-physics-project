import numpy as np
import pandas as pd
import scipy as sp
import scipy.fftpack
from scipy.constants import G as GRAVITY
from scipy.constants import g as G_EARTH


def get_distance_matrix(location):
    """
    Computes the mash of distances between planets
    :param location: dataframe containing stellar objects' locations at a certain time
    :return: dataframe containing distances between every combination of two planets.
    """
    coordinates = location[['x', 'y', 'z']].to_numpy()
    size = len(location.index)
    distances = np.zeros([size, size])

    for i, planet_crd in enumerate(coordinates):
        crd = coordinates - planet_crd  # Subtracts planet_crd from each row of coordinates
        distances[:, i] = np.sqrt(np.sum(crd ** 2, axis=1))

    distances = pd.DataFrame(distances)
    distances.index = location.index
    distances.columns = location.index

    return distances


def get_distance_sequence(location_sequence):
    """
    Computes the mash of distances between planets over a sequence.
    :param location_sequence: dataframe containing a sequence of planetary locations
    :return: dataframe containing sequence of mash distances between planet
    """
    steps = len(location_sequence.index.unique('step'))
    distances = []

    for s in range(1, steps + 1):
        distance_matrix = get_distance_matrix(location_sequence.xs(s, level='step'))
        distance_matrix['step'] = s
        distances.extend([distance_matrix])

    distances = pd.concat(distances)
    distances.set_index(['step'], append=True, inplace=True)

    return distances


def get_forces(motion_data):
    """
    Computes the gravitational forces applied on each planet
    :param motion_data: dataframe containing stellar objects' details at a certain time
    :return: dataframe containing total force vector for each planet
    """
    coordinates = motion_data.loc[:, ['x', 'y', 'z']].to_numpy()
    mass = np.array([motion_data.loc[:, 'mass'].to_numpy()])

    # mass_mat ∈ Mat(9X9)
    mass_mat = mass.T @ mass
    forces = np.zeros([9, 9, 3])
    distances = get_distance_matrix(motion_data.loc[:, ['x', 'y', 'z']]).to_numpy()

    for i, planet_crd in enumerate(coordinates):
        crd = coordinates - planet_crd  # Subtracts planet_crd from each row of coordinates
        forces[i, :, :] = crd

    # magnitudes ∈ Mat(9X9)
    magnitudes = np.divide(mass_mat, distances ** 3, out=np.zeros_like(mass_mat), where=distances != 0) * GRAVITY

    for i in range(3):
        forces[:, :, i] = np.multiply(forces[:, :, i], magnitudes)

    # total_forces ∈ Mat(9X3)
    total_forces = np.sum(forces, axis=1)

    # total_forces_df ∈ Mat(9X4)
    # Columns: [name, fx, fy, fz]
    total_forces_df = pd.DataFrame(total_forces, columns=['fx', 'fy', 'fz'])
    total_forces_df.index = motion_data.index

    return total_forces_df


def get_accelerations(force, mass):
    """
    Computes the acceleration of each planet
    :param force: dataframe containing force vector for each planet
    :param mass: dataframe containing mass of each planet
    :return: dataframe containing acceleration vector for each planet
    """
    forces = force.loc[:, ['fx', 'fy', 'fz']].to_numpy()
    masses = mass.to_numpy()
    accelerations = np.zeros([9, 3])

    for i in range(3):
        accelerations[:, i] = np.divide(forces[:, i], masses)

    # accelerations_df ∈ Mat(9X4)
    # Columns: [name, ax, ay, az]
    accelerations_df = pd.DataFrame(accelerations, columns=['ax', 'ay', 'az'])
    accelerations_df.index = force.index

    return accelerations_df


def get_velocities(initial_velocity, acceleration, time_delta):
    """
    Computes the velocity of each planet after time_delta seconds.
    :param initial_velocity: dataframe containing initial velocities vector for each planet
    :param acceleration: dataframe containing acceleration vector for each planet
    :param time_delta: time in seconds to advance positions
    :return: dataframe containing velocity vector for each planet
    """
    initial_velocities = initial_velocity.loc[:, ['vx', 'vy', 'vz']].to_numpy()
    accelerations = acceleration.loc[:, ['ax', 'ay', 'az']].to_numpy()

    velocities = initial_velocities + accelerations * time_delta

    # velocities_df ∈ Mat(9X4)
    # Columns: [name, vx, vy, vz]
    velocities_df = pd.DataFrame(velocities, columns=['vx', 'vy', 'vz'])
    velocities_df.index = initial_velocity.index

    return velocities_df


def get_locations(initial_location, velocity, time_delta):
    """
    Computes the location of each planet after time_delta seconds.
    :param initial_location: dataframe containing initial locations vector for each planet
    :param velocity: dataframe containing velocity vector for each planet
    :param time_delta: time in seconds to advance positions
    :return: dataframe containing location vector for each planet
    """
    initial_locations = initial_location.loc[:, ['x', 'y', 'z']].to_numpy()
    velocities = velocity.loc[:, ['vx', 'vy', 'vz']].to_numpy()

    locations = initial_locations + velocities * time_delta

    # locations_df ∈ Mat(9X4)
    # Columns: [name, x, y, z]
    locations_df = pd.DataFrame(locations, columns=['x', 'y', 'z'])
    locations_df.index = initial_location.index

    return locations_df


def get_kinetic_energy(mass, velocity):
    """
    Computes the kinetic energy of each planet at a given moment.
    :param mass: dataframe containing mass of each planet
    :param velocity: dataframe containing velocity vector for each planet
    :return: dataframe containing kinetic energy for each planet
    """
    coordinates = velocity[['vx', 'vy', 'vz']].to_numpy()
    masses = mass.to_numpy()
    energy = (1/2) * np.multiply(masses, np.sum(coordinates ** 2, axis=1))

    energy = pd.DataFrame(energy)
    energy.index = velocity.index
    energy.columns = ['kinetic']

    return energy


def get_kinetic_sequence(mass, velocity_sequence):
    """
    Computes the kinetic energy of each planet over a sequence.
    :param mass: dataframe containing mass of each planet
    :param velocity_sequence: dataframe containing a sequence of planetary velocities
    :return: dataframe containing a sequence of kinetic energies for each planet
    """
    steps = len(velocity_sequence.index.unique('step'))
    energies = []

    for s in range(1, steps + 1):
        energy = get_kinetic_energy(mass, velocity_sequence.xs(s, level='step'))
        energy['step'] = s
        energies.extend([energy])

    energies = pd.concat(energies)
    energies.set_index(['step'], append=True, inplace=True)

    return energies


def get_potential_energy(mass, location):
    """
    Computes the potential energy of each planet at a given moment.
    :param mass: dataframe containing mass of each planet
    :param location: dataframe containing velocity vector for each planet
    :return: dataframe containing potential energy for each planet
    """
    coordinates = location[['x', 'y', 'z']].to_numpy()
    coordinates = coordinates - location.loc['sun'][['x', 'y', 'z']].to_numpy()
    masses = mass.to_numpy()
    distances = np.sqrt(np.sum(coordinates ** 2, axis=1))
    energy = G_EARTH * np.multiply(masses, distances)

    energy = pd.DataFrame(energy)
    energy.index = location.index
    energy.columns = ['potential']

    return energy


def get_potential_sequence(mass, location_sequence):
    """
    Computes the potential energy of each planet over a sequence.
    :param mass: dataframe containing mass of each planet
    :param location_sequence: dataframe containing a sequence of planetary locations
    :return: dataframe containing a sequence of potential energies for each planet
    """
    steps = len(location_sequence.index.unique('step'))
    energies = []

    for s in range(1, steps + 1):
        energy = get_potential_energy(mass, location_sequence.xs(s, level='step'))
        energy['step'] = s
        energies.extend([energy])

    energies = pd.concat(energies)
    energies.set_index(['step'], append=True, inplace=True)

    return energies


def get_frequency(distance_sequence, scale=365):
    """
    Estimates the frequency of a planet given his distance from the sun over time.
    The estimation is using Fast Fourier Transform to calculate the coefficients of the
    discrete Fourier Sum representing the planet's distance from the sun function over time.
    The get_frequency function then takes the coefficient which represents the frequency with
    the highest amplitude to be the "base frequency" of the motion and returns it as the
    frequency of the planet.
                            ^
             .......- - - - 1
          ...   |   ...     |
        ..             ..   |      f(x) = sin(x)
       .        |        .  |
      .                   . |
     .          |          .|
    ----------pi/2----------.---------pi/2---------->
                            |.          |          .
                            | .                   .
                            |  .        |        .
                            |   ..             ..
                            |     ...   |   ...
                           -1- - - - .......

    :param distance_sequence: dataframe containing a sequence of the distance of the planet from the sun
    :param scale: integer representing the scale of the resulting frequencies in days.
    :return: float containing the base frequency as extracted using fft and as described above
    """
    dist_fft = sp.fftpack.fft(distance_sequence.to_numpy())
    dist_psd = np.abs(dist_fft) ** 2
    fft_freq = sp.fftpack.fftfreq(len(dist_psd), 1. / scale)

    pstv_freq_idx = fft_freq > 0

    base_freq_idx = (10 * np.log10(dist_psd[pstv_freq_idx])).argmax()
    base_freq = fft_freq[base_freq_idx]

    if base_freq == 0:
        dist_psd = np.delete(dist_psd, base_freq_idx)
        pstv_freq_idx = np.delete(pstv_freq_idx, base_freq_idx)
        fft_freq = np.delete(fft_freq, base_freq_idx)

        base_freq_idx = (10 * np.log10(dist_psd[pstv_freq_idx])).argmax()
        base_freq = fft_freq[base_freq_idx]

    return base_freq


