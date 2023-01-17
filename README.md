# Discrete simulation of the Solar System

This is the final assignment of the Computers for Physics course in Tel-Aviv University - 2023, Sem A

## Requirements

1. Simulate all the planetary motion in our solar system using Newtonian mechanics given number of steps, time delta between steps, 
and initial location, velocity and mass of each planet.
2. Plot the trajectories of all the planets in 3 dimensions.
3. Plot the distance of each planet from the Sun over time.
4. Plot the total energy (kinetic + potential) of the Solar System over time.
5. Show that your simulation adheres to Kepler's third law:
    * Calculate the time period of each planet's orbital period, denoted by $T$.
    * Calculate the mean distance of each planet from the Sun, denoted by $D$.
    * Show that $T^{2} \propto D^{3}$

*<ins>Note</ins>: using NumPy may be rewarded by up to 30 bonus pts.*

## Plots

### Trajectories
![Planetary trajectory over a 100-year period using 1 day time delta](/assets/trajectories-plot.png "Trajectories over 100 years").

### Distance from the Sun
![Planetary distance from the Sun over a 100-year period using 1 day time delta](/assets/distance-from-sun-plot.png "Distance from the Sun over 100 years").

### Total energy
![Total energy in the Solar System over a 100-year period using 1 day time delta](/assets/total-energy-plot.png "Total energy over 100 years").

### Kepler's third law (Mercury, Venus, Earth and Mars)
![Kepler's third law as calculated from the simulation](/assets/keplers-third-law-plot.png "Kepler's third law - linear correlation").
