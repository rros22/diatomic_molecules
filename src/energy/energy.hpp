#ifndef ENERGY_HPP
#define ENERGY_HPP

#include "../diatomic_molecule/diatomic_molecule.hpp"

void ke_trans(Mol_2 *molecules, int i);

void ke_rot(Mol_2 *molecules, int i);

void ke_vibr(Mol_2 *molecules, int i);

void pe(Mol_2 *molecules, int i);

void e_total(Mol_2 *molecules, int i);

//run all functions
void calculate_energies(Mol_2 *molecules, int i);

#endif
