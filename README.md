# Project A Solar System
Computer Modelling Project A Solar System

Describes an N-body system interacting through Newtonian gravity.

Units are kilograms, astronomical units, and days.
The derived units for force and energy are kg⋅au⋅day^-2 and kg⋅au^2⋅day^-2.

## Inputs
### Command Line
The names of:
1. The file containing simulation parameter
2. The file containing the particle data
3. The file the trajectory data are to be written to
4. The file the observables (like apsides, orbital periods) are to be written to.

Should look like: `python simulator.py parameters.txt particles.txt trajectory.xyz observables.txt`
### Simulation Parameters
A file with the following structure
```
    no_of_iterations  : The number of iterations in a simulation.
    timestep          : The time that passes over an iteration (d).
    no_of_bodies      : The number of bodies to be read from the body file.
    iter_for_output   : n such that every n-th iteration's data, like the positions for
                        the trajectory and total energy of the system, will be saved to file.
    body_orbits       : A list where the i-th item in the list (not including the sun)
                        is the body the i-th body orbits around. Depends on the ordering
                        in particle input. See below for examples.
```


### Particle Information
Each row represents a particle to be generated, and is laid out as such:
```
label, mass, x_position, y_position, z_position, x_velocity, y_velocity, z_velocity
```
### Examples
###### Example 1
A 3 body system, such as the Sun, Mercury, and Venus, over 300 days with a time interval of a day, saving the data each day.
The order of the particles in the input is:

<ol start="0">
  <li>The Sun</li>
  <li>Mercury</li>
  <li>Venus</li>
</ol>

As both Mercury and Venus orbit the Sun, the list for `body_orbits` will look like this: `0, 0`. The first item represents the body Mercury orbits, the second is the body Venus orbits. Since the Sun is the 0th body in the input, the value in the list is 0.

The simulation parameters file:
```
    no_of_iterations  : 300
    timestep          : 1
    no_of_bodies      : 3
    iter_for_output   : 1
    body_orbits       : 0, 0
```
###### Example 2
A 3 body system, such as the Sun, Earth, and Earth's Moon, over 1000 days with a time interval of a 2 days, saving the data every 10 days.
The order of the particles in the input is:

<ol start="0">
  <li>The Sun</li>
  <li>Earth</li>
  <li>Earth's Moon</li>
</ol>

As the Earth orbits the Sun and the Moon orbits the Earth, the list for `body_orbits` will look like this: `0, 1`. The first item represents the body the Earth orbits, the second is the body the Moon orbits. Since the Sun is the 0th body in the input and the Earth is the 1st, the values in the list are 0 and 1.

The simulation parameters file:
```
    no_of_iterations  : 500
    timestep          : 2
    no_of_bodies      : 3
    iter_for_output   : 5
    body_orbits       : 0, 1
```

## Outputs
### Trajectory File

### Observables

### Energy
