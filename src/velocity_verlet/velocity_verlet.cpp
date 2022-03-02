#include "velocity_verlet.hpp"
#include "../forces/forces.hpp"
#include <iostream>
#include <cmath>

#define k_e 8.988E9

void compute_forces(Box *domain, bool n){

    //reset force values to 0
    if (!n) {

        reset_forces(domain);

    }

    else {

        reset_forces_n(domain);

    }


    //parameters for reading ease
    int molecule_no = domain->molecule_no;
    Mol_2 *molecules = domain->molecules;


    //compute bonded interactions
    for (int i = 0; i < molecule_no; i++){

        bond_force(&(molecules[i]), n);

    }


    //compute non bonded interactions between non-co-molecular atoms
    for (int i = 0; i < molecule_no; i++){

        for (int j = i + 1; j < molecule_no; j++){

            lennard_jones_force(&(molecules[i].atom_1), &(molecules[j].atom_1), n);
            lennard_jones_force(&(molecules[i].atom_1), &(molecules[j].atom_2), n);
            lennard_jones_force(&(molecules[i].atom_2), &(molecules[j].atom_1), n);
            lennard_jones_force(&(molecules[i].atom_2), &(molecules[j].atom_2), n);

        }
    }

    double d = 2*7*4.33E-10; //domain side

    //compute reflective boundary conditions
    for (int i = 0; i < molecule_no; i++){

        wall_force(&(molecules[i].atom_1), n, d);
        wall_force(&(molecules[i].atom_2), n, d);

    }

    //compute gravity force
    for (int i = 0; i < molecule_no; i++){

        gravity(&(molecules[i].atom_1), n);
        gravity(&(molecules[i].atom_2), n);

    }

}

void next_position(Box *domain, double dt){

    int molecule_no = domain->molecule_no;
    Mol_2 *molecules = domain->molecules;

    for (int i = 0; i < molecule_no; i++){

        molecules[i].atom_1.x = molecules[i].atom_1.x + dt*(molecules[i].atom_1.u) + pow(dt, 2)/2/(molecules[i].atom_1.m)*(molecules[i].atom_1.f_x);
        molecules[i].atom_1.y = molecules[i].atom_1.y + dt*(molecules[i].atom_1.v) + pow(dt, 2)/2/(molecules[i].atom_1.m)*(molecules[i].atom_1.f_y);
        molecules[i].atom_1.z = molecules[i].atom_1.z + dt*(molecules[i].atom_1.w) + pow(dt, 2)/2/(molecules[i].atom_1.m)*(molecules[i].atom_1.f_z);

        molecules[i].atom_2.x = molecules[i].atom_2.x + dt*(molecules[i].atom_2.u) + pow(dt, 2)/2/(molecules[i].atom_2.m)*(molecules[i].atom_2.f_x);
        molecules[i].atom_2.y = molecules[i].atom_2.y + dt*(molecules[i].atom_2.v) + pow(dt, 2)/2/(molecules[i].atom_2.m)*(molecules[i].atom_2.f_y);
        molecules[i].atom_2.z = molecules[i].atom_2.z + dt*(molecules[i].atom_2.w) + pow(dt, 2)/2/(molecules[i].atom_2.m)*(molecules[i].atom_2.f_z);
    }

}

void next_velocity(Box *domain, double dt){

    int molecule_no = domain->molecule_no;
    Mol_2 *molecules = domain->molecules;

    for (int i = 0; i < molecule_no; i++){

        molecules[i].atom_1.u = molecules[i].atom_1.u + dt/2/(molecules[i].atom_1.m)*(molecules[i].atom_1.f_x + molecules[i].atom_1.f_x_n);
        molecules[i].atom_1.v = molecules[i].atom_1.v + dt/2/(molecules[i].atom_1.m)*(molecules[i].atom_1.f_y + molecules[i].atom_1.f_y_n);
        molecules[i].atom_1.w = molecules[i].atom_1.w + dt/2/(molecules[i].atom_1.m)*(molecules[i].atom_1.f_z + molecules[i].atom_1.f_z_n);

        molecules[i].atom_2.u = molecules[i].atom_2.u + dt/2/(molecules[i].atom_2.m)*(molecules[i].atom_2.f_x + molecules[i].atom_2.f_x_n);
        molecules[i].atom_2.v = molecules[i].atom_2.v + dt/2/(molecules[i].atom_2.m)*(molecules[i].atom_2.f_y + molecules[i].atom_2.f_y_n);
        molecules[i].atom_2.w = molecules[i].atom_2.w + dt/2/(molecules[i].atom_2.m)*(molecules[i].atom_2.f_z + molecules[i].atom_2.f_z_n);
    }

}

void verlet_integrate(Box *domain, double dt, double t){

    compute_forces(domain, 0);
    next_position(domain, dt);
    compute_forces(domain, 1);
    next_velocity(domain, dt);


}
