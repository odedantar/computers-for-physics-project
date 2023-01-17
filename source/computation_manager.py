import math
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join

from source import simulation as sim


MAX_STEPS = 365


def record_simulation(steps, time_delta, initial_data, dir_path):
    chunks = [MAX_STEPS for i in range(steps // MAX_STEPS)]
    if steps % MAX_STEPS > 0:
        chunks.append(steps % MAX_STEPS)
    chunks_count = len(chunks)
    id_len = int(math.ceil(math.log10(chunks_count)))

    print('Chunk Progress: ' + str(0) + ' / ' + str(chunks_count), end='\r')

    for i, chunk in enumerate(chunks):
        motion_sequence = sim.simulate(chunk, time_delta, initial_data)
        initial_data = motion_sequence.xs(chunk, level='step')
        initial_data.insert(loc=6, column='step', value=[0]*9)
        motion_sequence.to_csv(dir_path + 'chunk_' + format(i, '0' + str(id_len) + 'd'))
        print('Chunk Progress: ' + str(i+1) + ' / ' + str(chunks_count), end='\r')


def load_simulation(dir_path):
    files = [f for f in listdir(dir_path) if isfile(join(dir_path, f))]
    id_len = int(math.ceil(math.log10(len(files))))
    files.sort(key=lambda f: int(f[-id_len:]))
    motion_list = []

    step = 0
    for f in files:
        motion_sequence = pd.read_csv(join(dir_path, f))
        motion_sequence['step'] = motion_sequence['step'] + step
        step = max(motion_sequence['step'].unique())
        motion_list.extend([motion_sequence])

    motion_sequence = pd.concat(motion_list)
    motion_sequence.set_index(['name', 'step'], inplace=True)

    return motion_sequence


def record_distances(motion_sequence, dir_path):
    planets = motion_sequence.index.unique('name').drop('sun')
    distances = []

    for i, planet in enumerate(planets):
        distance = sim.get_sun_distances(motion_sequence[['x', 'y', 'z']], planet)
        distances.extend([distance])

    distances = pd.concat(distances, axis=1)
    distances.to_csv(join(dir_path, 'distances'))


def load_distances(dir_path):
    files = [f for f in listdir(dir_path) if isfile(join(dir_path, f))]

    distances = pd.read_csv(join(dir_path, files[0]))
    distances.set_index(['step'], inplace=True)

    return distances


if __name__ == '__main__':
    pass
