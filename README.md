# Project A Solar System
Computer Modelling Project A Solar System

Describes an N-body system interacting through Newtonian gravity.

Units are kilograms, astronomical units, and days.
The derived units for force and energy are kg⋅au⋅day^-2 and kg⋅au^2⋅day^-2.

## Structure

The velocity Verlet simulation and calculation of observables happen in the file `simulator.py`.

This uses a custom SolarSystem class, defined in `solar_system.py`. This is where the bulk of the methods are. Essentially it is a list of Particle3D objects from the first semester. This class is defined in the file `particle3d.py`.

## Running the code

To run the code there are sample parameter and particle information files. Using these will run a simulation of all 12 bodies for 100 years, with a time-step of 1 day. It will save the energies and trajectory every day. The planet data was mainly taken from [JPL Horizons](https://ssd.jpl.nasa.gov/horizons.cgi).

To run this simulation execute the following command from terminal/command line/preferred software (should take < 1min):

```
python simulator.py params.txt our_system.txt traj.xyz obs.txt
```

See below for more information about what each argument means. The energy.txt, obs.txt and traj.xyz files should match to the provided sample_energy.txt, sample_obs.txt and sample_traj.xyz files. See the outputs section for more information about the format of these outputs. It is expected that not all observables could be calculated as some planets have orbits over 100 years. The file 1000years_obs.txt is representative of what the observables look like when calculated over a long enough time, with a time-step of 1 day.

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

0. The Sun
1. Mercury
2. Venus

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

0. The Sun
1. The Earth
2. The Earth's Moon

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

The trajectory file is a `.xyz` file, so is made up of repeating blocks:

```
<No. of particles>
Comment line (Normally shows time)
<Particle 1 Label> <Particle 1 x> <Particle 1 y> <Particle 1 z>
<Particle 2 Label> <Particle 2 x> <Particle 2 y> <Particle 2 z>
<Particle 3 Label> <Particle 3 x> <Particle 3 y> <Particle 3 z>
...
<Particle N Label> <Particle N x> <Particle N y> <Particle N z>
```

For each time-step where the positions have been saved (determined by the simulation parameter `iter_for_output`) there is one of these blocks.

### Observables

The observables measured are the orbital period, the periapsis, and the apoapsis. The observables file is comma-delimited, where each row represents a planet or moon. The period is in days, and the peri- and apoapsides and respective errors are measured in au.

Each row looks like:

```
label, period, periapsis, periapsis error, apoapsis, apoapsis error
```

#### Averaging and errors

All of the observables are averaged. For the period, it is averaged by counting how many orbits there have been and dividing by the time. The error on the period is simple given by $$\sqrt{2}dt$$ for all bodies.

For the apsides, the peri- and apoapsis if each orbit is found, and these are averaged. The errors are the standard error on the mean, which is found by calculating the standard deviation of the set and dividing by the square root of the number of values found.

This means that if only one value is found, the standard deviation and therefore error is calculated as 0. It is very unlikely that an error would actually be 0. **An error of 0 probably just means only one value was found.**

#### Not enough data

If the simulation isn't run for long enough then planets won't make full orbits around the Sun, and might not even reach the closest/furthest points of their orbits. When this happens the program returns `nan`. This is to avoid returning an incorrect value for the observables.

### Energy

The total energy of the system is saved in a .txt file called `energy.txt`. It is a comma-delimited file of 4 columns. The first column is time in days, the second third and fourth are the kinetic energy, potential energy, and total energy of the system. All energies are in kg⋅au^2^⋅day^-2^.

Each row looks like:

```
time, kinetic energy, potential energy, total energy
```

