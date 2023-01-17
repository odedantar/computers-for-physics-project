import pandas as pd


def open_data_csv(path):
    """
    Opens a .csv file containing planetary data and returns its content as a pandas dataframe
    :param path: the path of the file to open
    :return: a Pandas dataframe with the planetary data from the file
    """
    data = pd.read_csv(path)
    data.set_index(['name'], inplace=True)
    data.columns = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'step', 'mass']

    return data

