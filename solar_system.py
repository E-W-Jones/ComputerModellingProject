# -*- coding: utf-8 -*-
"""
Author: Evan Jones - s1857441
Computer Modelling Project A: Solar System

A Module for having a solar system, consisting of an arbitrary number of bodies,
that can include stars, planets, moons, etc, acting under Newtonian gravity,
with a velocity Verlet integration scheme.

The functions/information it has:
    The mass, posistion, velocity, and a label for every body in the system
    the number of bodies in the system
    the product of the masses of each pair of bodies in the system
    the vector and scalar separations between each pair of bodies
    the gravitational potential- and total- energy of the system
    the total gravitational force on each body
    Instance methods to:
        print the locations of each body in xyz format
        update the position of every body
        update the velocity of every body
        correct for motion of the centre of mass
"""
import numpy as np
from particle3D import Particle3D


class SolarSystem():
    """
    A class for simulating a solar system, as a list of particle3d instances.

    Propeties: The number of bodies, the list of all the particle3D objects in
               the system, an of the product of every pair of masses,
               an array of the vector separations between every pair of bodies,
               and an array of the magnitude of those separations.
    Methods:   Print an xyz compliant string for each body; create the mass
               matrix above; create the vector separation and magnitude arrays,
               calculate the total gravitational potential energy of the system;
               calculate the total gravitational force on each body; calculate
               the total energy of the system; update the positions of each body
               in the system (velocity Verlet); update the velocities of each
               body in the system; and correct for centre of mass drift of
               the system.
    """

    def __init__(self, filename, no_of_bodies, delimiter=','):
        """
        Construct instance based on information in a file.

        Parameters
        ----------
        filename : name of the file to read the particles from.
        no_of_bodies : The number of bodies in the solar system; to be read
        from the input file.
        delimiter : The delimiter used in the data file.

        """
        self.N = no_of_bodies

        with open(filename, 'r') as filein:
            self.bodies = [Particle3D.from_file(filein, delimiter)
                           for i in range(self.N)]

        self.get_mass_products()

    def xyz(self, comment=None):
        """
        Generate a .xyz string for all the bodies in the solar system.

        Uses the particle3d str() magic method to output an xyz compatible
        file of the locations of all the bodies in the solar system:

            <N, the no. of bodies>
            <optional comment>
            <body1 label> <body1 x> <body1 y> <body1 z>
            <body2 label> <body2 x> <body2 y> <body2 z>
            <body3 label> <body3 x> <body3 y> <body3 z>
            ...
            <bodyN label> <bodyN x> <bodyN y> <bodyN z>

        Parameters
        ----------
        comment : An optional comment to include, as a string

        Returns
        -------
        .xyz compatible string of body locations

        """
        # Create lists of the strings to write in xyz format, a header and body
        if comment is None:
            header = [str(self.N)]
        else:
            header = [str(self.N), comment]

        body = [str(body) for body in self.bodies]

        return "\n".join(header + body + [''])  # empty list adds final newline

    def get_mass_products(self):
        """
        Create an NxN array where each element is the product mi*mj.

        Creates a property called mass_products, where each element is the
        product of the masses of two different particles:
            self.mass_products[i, j] = mi * mj

        This is simply the outer product of a mass vector with itself.

        As these values will not change throught the simulation, it makes
        sense to only calculate them once, as they will be used at every
        timestep in Newton's law of gravitation:
            F_ij = G * mi * mj * (rj - ri) / abs|rj - ri|^3

        """
        masses = np.array([body.mass for body in self.bodies])
        self.mass_products = np.outer(masses, masses)

    def get_pair_separations(self):
        """
        Create numpy arrays of vector and normal separations between bodies.

        An NxN matrix where the element (i,j) is the separation of
        particle i and particle j: r_j - r_i, as a 3d vector.
        This vector points from particle i to particle j, designed with
        calculating F_ij, the force on particle i from particle j.

        """
        positions = np.array([body.pos for body in self.bodies], ndmin=3)
        separations = positions - positions.transpose((1, 0, 2))

        self.sep_vecs = separations
        self.sep_mags = np.linalg.norm(separations, axis=2)

    def gravitational_potential(self, G):
        """
        Calculate the total potential energy due to gravity of each body.

        Uses:
            U_ij = - G * mi * mj / abs|rj - ri|
        where:
            U_ij is the potential energy body i has due to body j.
            mi is the mass of body i and mj is the mass of body j.
            ri and rj are the position vectors of the bodies.
            G is the universal gravitation constant.

        Parameters
        ----------
        G : The Universal Constant of Gravitation

        Returns
        -------
        U-total : The total gravitational potential energy of the system.

        """
        mimj = self.mass_products
        # An array where each element is the reciprocal of self.sep_mags
        recip_r = np.divide(1, self.sep_mags, where=(self.sep_mags != 0),
                            out=np.zeros_like(self.sep_mags))
        # This is done to avoid division by 0, which would happen where ri=rj.
        # This is when the ith and jth body represent the same body.
        # There is no physical potential or force felt on a body by itself,
        # so set it to 0, as this means we can sum over the array.

        U_ij = -1 * G * mimj * recip_r
        # This array has the potential each body feels due to each other body,
        # but this actually overcounts the potential as two bodies do not each
        # feel a potential, they share a potential caused by both.

        # The total gravitational potential of the system
        U_total = np.sum(U_ij) / 2  # U_ij overcounts by a factor of 2
        return U_total

    def gravitational_force(self, G):
        """
        Calculate the force due to gravity of each body.

        Uses:
            F_ij = G * mi * mj * (rj - ri) / abs|rj - ri|^3
        where:
            F_ij is the force felt by body i due to body j.
            mi is the mass of body i and mj is the mass of body j.
            ri and rj are the position vectors of the bodies.
            G is the universal gravitation constant.


        Parameters
        ----------
        G : The Universal Constant of Gravitation

        Returns
        -------
        An array of the forces acting on the i-th particle; [F_1, F_2, ...]

        """
        rij = self.sep_vecs
        # Add another dimension to scalars so can be broadcast with vectors
        mimj = self.mass_products[:, :, np.newaxis]
        # An array where each element is the reciprocal of self.sep_mags
        recip_r = np.divide(1, self.sep_mags, where=(self.sep_mags != 0),
                            out=np.zeros_like(self.sep_mags))[:, :, np.newaxis]
        # This is done to avoid division by 0, which would happen where ri=rj.
        # This is when the ith and jth body represent the same body.
        # There is no physical potential or force felt on a body by itself,
        # so set it to 0, as this means we can sum over the array.
        F_ij = G * mimj * recip_r**3 * rij
        # An array where each element is the total force acting on a body
        F_i = np.sum(F_ij, axis=1)
        return F_i

    def update_velocity(self, dt, force_list):
        """
        Update the velocity of each body in the system.

        Updates each body in the system over a time interval dt and
        force f, using velocity verlet integration in the Particle3D
        velocity_step method.

        Parameters
        ----------
        dt : The time step, as a float or int.
        force_list : A numpy array of the forces on each body.

        """
        for i in range(self.N):
            self.bodies[i].velocity_step(dt, force_list[i])

    def update_position(self, dt, force_list):
        """
        Update the position of each body in the system.

        Updates each body in the system over a time interval dt and
        force f, using velocity verlet integration in the Particle3D
        position_step method.

        Parameters
        ----------
        dt : The time step, as a float or int.
        force_list : A numpy array of the forces on each body.

        """
        for i in range(self.N):
            self.bodies[i].position_step(dt, force_list[i])

    def com_correct(self):
        """Subtract motion of centre of mass from each body."""
        com_velocity = Particle3D.totalmass_comvelocity(self.bodies)[1]
        for body in self.bodies:
            body.vel -= com_velocity

    def total_energy(self, G):
        """
        Calculate the total energy of the system.

        Parameters
        ----------
        G : The Universal Constant of Gravitation.

        Returns
        -------
        The total energy of the system.

        """
        kinetic = Particle3D.total_KE(self.bodies)
        potential = self.gravitational_potential(G)
        return kinetic, potential, kinetic + potential
