# -*- coding: utf-8 -*-
"""
Author: Evan Jones - s1857441

A module for having a particle of mass m, with a position, velocity, etc.

The functions/information it has:
    particle’s mass,
    a label,
    position and velocity
    a __str__ method that prints in xyz format
    instance methods that find:
         kinetic energy
         linear momentum
         the velocity at the next timestep
         the first and second order positions at the next timestep
    static methods that:
        find total KE of a list of particles
        find total mass of a list, and CoM velocity of a list of particles
        creates a particle using data from a file

"""
import numpy as np


class Particle3D(object):
    """
    Stores information about a particle in 3D space.

    Properties: label, mass, position in 3D (x, y, z), velocity (vx, vy, vz).
    Methods: Print xyz format; calculate kinetic energy; linear momentum;
             velocity after a time dt; first order (symplectic Euler) and
             second order (velocity Verlet) positions after time dt;
             total kinetic energy in a list; total mass of a list AND velocity
             of the centre of mass of a list; and a method to create a particle
             instance using data from a file.
    """

    def __init__(self, label, mass, position, velocity):
        """
        Construct instance.

        Parameters
        ----------
        mass : Mass, a float
        label : The label, a strong
        position : 3D position, as a list (or numpy array)
        velocity : 3D velocity, as a list (or numpy array)

        """
        self.label = label
        self.mass = float(mass)
        self.pos = np.array(position, float)
        self.vel = np.array(velocity, float)

    def __str__(self):
        """
        Print label and position in xyz format.

        Returns
        -------
        <label> <x> <y> <z>

        """
        x, y, z = self.pos
        return self.label + f" {x} {y} {z}"

    def kinetic_energy(self):
        """
        Calculate kinetic energy of the particle.

        Returns
        -------
        Kinetic energy, a float.

        """
        return self.mass * np.dot(self.vel, self.vel) / 2

    def linear_momentum(self):
        """
        Calculate linear momentum of a particle.

        Returns
        -------
        The linear momentum as a numpy array.

        """
        return self.mass * self.vel

    def velocity_step(self, dt, force):
        """
        Update the velocity after a time interval dt and force f.

        v(t + dt) = v(t) + dt*f/m

        Parameters
        ----------
        dt : The time interval, a float or int
        force : Force vector, either as a numpy array

        """
        self.vel += dt * force / self.mass

    def position_step(self, dt, force):
        """
        Update position by 2nd order over time interval dt, under force f.

        Uses Verlet integration.

        r(t + dt) = r(t) + dt*v(t) + dt^2 * f / (2m)

        Parameters
        ----------
        dt : The time interval, a float or int
        force : Force vector, either as a numpy array

        """
        self.pos += dt * self.vel + dt * dt * force / (2 * self.mass)

    # Static Methods
    @staticmethod
    def total_KE(particles):
        """
        Calculate total kinetic energy of all the particles in the list.

        Parameters
        ----------
        particles : A list of particle3D objects

        Returns
        -------
        Total KE, a float

        """
        return sum([particle.kinetic_energy() for particle in particles])

    @staticmethod
    def totalmass_comvelocity(particle_list):
        """
        Calculate total mass and centre of mass velocity for particle list.

        v_com = total_momentum/total_mass

        Parameters
        ----------
        particle_list : A list of particle3D objects

        Returns
        -------
        Total mass, a float
        Centre of mass velocity, a numpy array

        """
        total_momentum = sum([particle.linear_momentum()
                              for particle in particle_list])
        total_mass = sum([particle.mass for particle in particle_list])

        return total_mass, total_momentum / total_mass

    @staticmethod
    def from_file(file_handle, delimiter=","):
        """
        Convert information in a file to particle3D objects.

        file layout:
            title, mass, x, y, z, vx, vy, vz
        Although it is possible to use any delimiter, the default is a comma.

        Parameters
        ----------
        file_handle : The handle of the file the particle information is in.
        delimiter : The seperator of data fields, by default is a comma.

        Returns
        -------
        A Particle3D object.

        """
        # Reading in the data
        line = file_handle.readline()
        line = line.rstrip("\n")
        line = line.split(delimiter)

        # Assigning the data
        label = line[0]
        mass = float(line[1])
        position = np.array(line[2:5], float)
        velocity = np.array(line[5:], float)

        return Particle3D(label, mass, position, velocity)