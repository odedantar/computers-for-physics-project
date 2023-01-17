from source import simulation
from source import parse_data
from source import plot_solar_system as plotter
from source import computation_manager as manager


# FYI - Project's Files Tree:
# ---------------------------
# project/
# ├─ data/
# │  ├─ manager/
# │  │  ├─ distance/
# │  │  ├─ motion/
# │  ├─ data.csv
# ├─ source/
# │  ├─ __init__.py
# │  ├─ computation_manager.py
# │  ├─ parse_data.py
# │  ├─ physics_engine.py
# │  ├─ plot_solar_system.py
# │  ├─ simulation.py
# ├─ main.py


if __name__ == '__main__':

    ################################################################################
    # --------------------------- Run a new simulation --------------------------- #
    ################################################################################

    INITIAL_DATA_FILE = './data/data.csv'

    steps = 5*365  # days in 5 year
    time_delta = 24*60*60  # seconds in a day
    initial_data = parse_data.open_data_csv(INITIAL_DATA_FILE)

    motion_sequence = simulation.simulate(steps, time_delta, initial_data)

    ################################################################################
    # -------------------------------- Plot charts ------------------------------- #
    ################################################################################

    # Notes Regarding MatPlotLib:
    # ---------------------------
    # Make sure to run each plot separately, meaning - all the function calls are
    # commented except for the call of the plot you want to run. If you call them
    # one after the other on the same run, you might run into problems with
    # MatPlotLib. I'm sure there are better ways to make it all work together,
    # but it's the end of the semester and I must work on other assignments.

    # plotter.plot_trajectory(motion_sequence)
    # plotter.plot_sun_distances(motion_sequence)
    # plotter.plot_energy(motion_sequence)

    # Notes Regarding Kepler's Third Law:
    # -----------------------------------
    # Note 1 - After experimenting with Kepler's third law and period evaluation,
    #          I must say that a planet requires at least 1.5 or 2 full periods of
    #          simulation so the function would approximate a close enough period.
    #
    # Note 2 - If we want to plot Kepler's third law for only a subset of planets,
    #          we must not forget to also include the Sun in our subset, since we
    #          calculate all the distances of the planets from the Sun.

    # plotter.plot_kepler_third_law(motion_sequence.loc[['sun', 'mercury', 'venus', 'earth', 'mars']])

    # plotter.show()

    ################################################################################
    # ---------------------- Record & Load a new simulation ---------------------- #
    ################################################################################
    # I've had problems running the simulation on my computer, so I wrote          #
    # this function which runs the simulation in chunks of 365 steps and saves     #
    # the chunks into files.                                                       #
    # The purpose was twofold:                                                     #
    # 1. If my computer crashed in the middle, the simulation that ran up until    #
    #    the crash was saved +- one chunk of simulation data.                      #
    # 2. Working with the data is much faster when all I need to do is load the    #
    #    data from the files instead of simulating it all over again.              #
    #                                                                              #
    # Note - If you want to record a new simulation you might have to manually     #
    #        delete all the previous files in the directory you chose, or just     #
    #        create a new directory instead.                                       #
    ################################################################################

    # MOTION_DIR_PATH = './data/manager/motion/'

    # manager.record_simulation(365*100, 24*60*60, initial_data, MOTION_DIR_PATH)
    # motion_sequence = manager.load_simulation(MOTION_DIR_PATH)

    ################################################################################
    # ----------------------- Record & Load new distances ------------------------ #
    ################################################################################
    # Same reason here as in the Record & Load a new simulation, except the        #
    # record_distances function doesn't work in chunks, but it did save time       #
    # later when I didn't have to recalculate the distances everytime.             #
    ################################################################################

    # DISTANCE_DIR_PATH = './data/manager/distance/'

    # motion_sequence = manager.load_simulation(MOTION_DIR_PATH)
    # manager.record_distances(motion_sequence, DISTANCE_DIR_PATH)
    # distances = manager.load_distances(DISTANCE_DIR_PATH)

    ################################################################################
    # --------------------------- Plot preloaded data ---------------------------- #
    ################################################################################
    # I've also written an option for function plot_sun_distances and              #
    # plot_kepler_third_law to be called with preloaded data.                      #
    # This is the way to do it:                                                    #
    ################################################################################

    # distances = manager.load_distances(DISTANCE_DIR_PATH)

    # plotter.plot_sun_distances(distances, preload=True)

    # Notes Regarding Kepler's Third Law:
    # -----------------------------------
    # Notice that unlike the non-preloaded version, had we wanted to plot for only
    # a subset of planets we would have to include the location of the Sun, here the
    # distances are precalculated, so not only we can exclude the Sun - we must do so.
    # The reason being that the Sun's distance from itself is 0, therefor not useful,
    # and is not included in the precalculated distances.

    # plotter.plot_kepler_third_law(distances[['mercury', 'venus', 'earth', 'mars']], preload=True)

    # plotter.show()

    ################################################################################
