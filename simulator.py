# -*- coding: utf-8 -*-
"""
Author: Evan Jones - s1857441.
Computer Modelling Project A: Solar System

Uses velocity verlet integration to model a system of planets
acting under Newtonian gravity.

See README for information about I/O.
"""
import sys
import numpy as np
from solar_system import SolarSystem as SolSys
from scipy.signal import find_peaks
from tqdm import trange

# --- Constants ---
DAY = 86400  # seconds

# from NASA - https://cneos.jpl.nasa.gov/glossary/au.html
AU = 149597870700  # metres
# from CODATA - https://physics.nist.gov/cuu/Constants/index.html
G = 6.6743E-11 * DAY**2 / (AU**3)  # au^3 kg^-1 days^-2


def calculate_observables(body_distances, dt):
    """
    Calculate the observables for one body in the system.

    Calculates the period, mean periapsis (of the periapsis of each orbit),
    error of the periapsis, mean apoapsis, and error of the apoapsis.

    When the values cannot be calculated due to a lack of orbital data, nan is
    returned.

    Parameters
    ----------
    body_distances : An array of distance from the body to the orbiting body
    dt : The time step

    Returns
    -------
    A list of the observables and errors in the following order:
        period, periapsis, periapsis error, apoapsis, apoapsis error

    """
    prominence = np.max(body_distances) / 100
    peak_times = find_peaks(body_distances, prominence=prominence)[0]
    trough_times = find_peaks(-body_distances, prominence=prominence)[0]

    # Number of peaks and troughs
    Npeak = len(peak_times)
    Ntrough = len(trough_times)

    # Peri- and Apo- apsides -- If there isn't a peak or trough then return nan
    if Ntrough < 1:
        pmean = np.nan
        perr = np.nan
    else:
        peris = body_distances[trough_times]
        pmean = np.mean(peris)
        perr = np.std(peris) / (len(peris)**0.5)

    if Npeak < 1:
        amean = np.nan
        aerr = np.nan
    else:
        apos = body_distances[peak_times]
        amean = np.mean(apos)
        aerr = np.std(apos) / (len(apos)**0.5)

    # Period -- there need to be at least two peaks or troughs to calculate it
    if Npeak > 1:
        t_0, t_N = peak_times[[0, -1]]
        period = dt * (t_N - t_0) / (Npeak - 1)
    elif Ntrough > 1:
        t_0, t_N = trough_times[[0, -1]]
        period = dt * (t_N - t_0) / (Ntrough - 1)
    else:
        period = np.nan

    return [period, pmean, perr, amean, aerr]


def save_observables(fname, orbital_distances, dt, system, delimiter=','):
    """
    Create a file of the apsides and orbital periods of the solar system.

    Each row is a body, with the columns:
        name, orbital_period, periapsis, apoapsis

    Parameters
    ----------
    fname : The name of the file the data are to be saved to.
    orbital_distances : An array where each row is a timestep and each column
    the distance from the orbiting body to the body being orbited.
    dt : The time that passes over one iteration in the simulation.
    system : The SolarSystem instance being simulated.
    delimiter : Optional. The delimiter for the file, by default is a comma.

    """
    data_array = [['name', 'orbital_period (days)', 'periapsis (au)', 'error',
                   'apoapsis (au)', 'error']]
    names = []
    for body in system.bodies:
        if body.label.lower() != "sun":
            names.append(body.label)

    N = len(names)

    body_distances = orbital_distances.T

    for i in range(N):
        data_array.append([names[i]] + calculate_observables(body_distances[i], dt))
    np.savetxt(fname, data_array, delimiter=delimiter, fmt='%s')


def main():
    # --- Check and read data from command line arguments ---
    if len(sys.argv) != 5:
        print(("Not enough command line arguments given, there should be 4:\n"
               "a parameter file, a file with the celestial bodies' "
               "properties, a filename for the trajectory file and a filename "
               "for the observable file.\nFor example:\n\t"
               "python simulation.py parameters.txt particles.txt "
               "trajectory.xyz observables.txt"))
        sys.exit()
    else:
        param_fname, bodies_fname, traj_fname, obs_fname = sys.argv[1:]
        # Open and read in parameters
        with open(param_fname, "r") as filein:
            params = [line.split(":")[-1] for line in filein.readlines()]
        num_of_steps = int(params[0])
        dt = float(params[1])
        no_of_bodies = int(params[2])
        n = int(params[3])
        orbits_measured_from = [int(x) for x in params[4].split(",")]

    # --- Initialise Solar System ---
    solar_system = SolSys(bodies_fname, no_of_bodies)
    solar_system.com_correct()

    # --- Initialise data arrays ---
    # Open the trajectory file and create arrays to store energies and
    # orbital distances (the distance from the primary body to the secondary)
    # Energy array has columns: time, kinetic, potential, total energy.
    trajectory_file = open(traj_fname, "w")
    energy = np.zeros((num_of_steps + 1, 4), dtype=float)
    orbital_distance = np.zeros((num_of_steps + 1, no_of_bodies - 1),
                                dtype=float)

    # --- Time Integrator ---
    # Calculate forces at time=0, save energy, orbital distances, and positions
    time = 0
    trajectory_file.write(solar_system.xyz(f"time = {time}"))
    solar_system.get_pair_separations()
    forces = solar_system.gravitational_force(G)
    energy[0] = time, *solar_system.total_energy(G)

    orbital_distance[0] = solar_system.sep_mags[orbits_measured_from,
                                                range(1, solar_system.N)]

    # Start integrating
    for i in trange(num_of_steps):
        solar_system.update_position(dt, forces)
        solar_system.get_pair_separations()

        new_forces = solar_system.gravitational_force(G)

        avg_forces = (forces + new_forces) / 2

        solar_system.update_velocity(dt, avg_forces)

        # Updating the force variables for the next timestep
        forces = new_forces

        time += dt

        # --- Writing Data ---
        # Update distances data array
        orbital_distance[i+1] = solar_system.sep_mags[orbits_measured_from,
                                                      range(1, solar_system.N)]
        # Update trajctory file and energies every nth timestep
        # Use i+1 so that counting starts from t=0 not t=dt
        if (i + 1) % n == 0:
            # Output particle positions
            trajectory_file.write(solar_system.xyz(f"time = {time}"))
            # Calculate and store energy
            energy[i + 1] = time, *solar_system.total_energy(G)

    # --- Saving Data ---
    # Close trajectory file, save energies, calculate and save observables
    trajectory_file.close()
    np.savetxt("energy.txt", energy[::n], delimiter=",")
    save_observables(obs_fname, orbital_distance, dt, solar_system)


if __name__ == "__main__":
    main()
