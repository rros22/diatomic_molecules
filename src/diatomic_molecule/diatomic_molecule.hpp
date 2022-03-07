#ifndef DIATOMIC_MOLECULE_HPP
#define DIATOMIC_MOLECULE_HPP

#include "../sphere/sphere.hpp"

typedef struct Mol_2{

    Sphere atom_1;
    Sphere atom_2;

    //bond parameters
    double k;
    double l_0;

    //positions
    double x;
    double y;
    double z;

    double y_rot;
    double z_rot;

    //velocities
    double u;
    double v;
    double w;

    double yy_rot;
    double zz_rot;

    //bond length
    double l;

    //energies
    double ke_trans;
    double ke_rot;
    double ke_vibr;

    double pe;

    double e_total;



} Mol_2;



#endif
