#ifndef COULOMB_HPP
#define COULOMB_HPP

#include "../box/box.hpp"

void reset_forces(Box *domain);
void reset_forces_n(Box *domain);

void set_forces(Sphere *atom, double fx, double fy, double fz);
void set_forces_n(Sphere *atom, double fx, double fy, double fz);

void coulomb_force(Sphere *atom_1, Sphere *atom_2, bool n);
void lennard_jones_force(Sphere *atom_1, Sphere *atom_2, bool n);

void bond_force(Mol_2 *molecule, bool n);

void wall_force(Sphere *atom, bool n, double d);

void gravity(Sphere *atom, bool n);

#endif
