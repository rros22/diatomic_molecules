#ifndef POST_PROCESSING_HPP
#define POST_PROCESSING_HPP

#include <fstream>
#include <iostream>

#include "../box/box.hpp"

void labels(std::string path);
void append_sphere_csv(Sphere *sphere, std::string path);
void labels_mol(std::string path);
void append_molecule_csv(Mol_2 *molecule, std::string path);
void box_csv(Box *box, std::string path);
void delete_frame(std::string path);

std::string to_angstroms(double meters);

std::string white_spaces(int spaces);

void box_pdb(Box *box, std::string path);

#endif
