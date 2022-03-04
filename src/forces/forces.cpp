#include "forces.hpp"
#include <cmath>

void reset_forces(Box *domain){

    Mol_2 *molecules = domain->molecules;
    int molecule_no = domain->molecule_no;

    for (int i = 0; i < molecule_no; i++){

        molecules[i].atom_1.f_x = 0;
        molecules[i].atom_1.f_y = 0;
        molecules[i].atom_1.f_z = 0;

        molecules[i].atom_2.f_x = 0;
        molecules[i].atom_2.f_y = 0;
        molecules[i].atom_2.f_z = 0;

    }
}

void reset_forces_n(Box *domain){

    Mol_2 *molecules = domain->molecules;
    int molecule_no = domain->molecule_no;

    for (int i = 0; i < molecule_no; i++){

        molecules[i].atom_1.f_x_n = 0;
        molecules[i].atom_1.f_y_n = 0;
        molecules[i].atom_1.f_z_n = 0;

        molecules[i].atom_2.f_x_n = 0;
        molecules[i].atom_2.f_y_n = 0;
        molecules[i].atom_2.f_z_n = 0;

    }
}

void set_forces(Sphere *atom, double fx, double fy, double fz){

    atom->f_x += fx;
    atom->f_y += fy;
    atom->f_z += fz;

}

void set_forces_n(Sphere *atom, double fx, double fy, double fz){

    atom->f_x_n += fx;
    atom->f_y_n += fy;
    atom->f_z_n += fz;

}


void coulomb_force(Sphere *atom_1, Sphere *atom_2, bool n){

    double k_e = 8.988E9;

    double q_1 = atom_1->q;
    double q_2 = atom_2->q;

    double x_1 = atom_1->x;
    double y_1 = atom_1->y;

    double x_2 = atom_2->x;
    double y_2 = atom_2->y;

    double f_x = - k_e*q_1*q_2*(x_2 - x_1)/pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2), 3/2);
    double f_y = - k_e*q_1*q_2*(y_2 - y_1)/pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2), 3/2);
    double f_z = 0;

    if (!n){

        set_forces(atom_1, f_x, f_y, f_z);
        set_forces(atom_2, -f_x, -f_y, -f_z);

    }

    else{

        set_forces_n(atom_1, f_x, f_y, f_z);
        set_forces_n(atom_2, -f_x, -f_y, -f_z);

    }
}

//4*epsilon*((sigma/sqrt((x_2-x_1)^2+(y_2-y_1)^2))^12-(sigma/sqrt((x_2-x_1)^2+(y_2-y_1)^2))^6)

void lennard_jones_force(Sphere *atom_1, Sphere *atom_2, bool n){

    double epsilon =  6.57E-22; // 4.184 * 0.2385;  kJ/mol
    double sigma = 3.17E-10;//3.4E-10;

    double x_1 = atom_1->x;
    double y_1 = atom_1->y;
    double z_1 = atom_1->z;

    double x_2 = atom_2->x;
    double y_2 = atom_2->y;
    double z_2 = atom_2->z;

    double f_x = - 4*epsilon*((12*pow(sigma, 12)*(x_2 - x_1)/pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2), 7)) - (6*pow(sigma, 6)*(x_2 - x_1)/pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2), 4)));
    double f_y = - 4*epsilon*((12*pow(sigma, 12)*(y_2 - y_1)/pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2), 7)) - (6*pow(sigma, 6)*(y_2 - y_1)/pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2), 4)));
    double f_z = - 4*epsilon*((12*pow(sigma, 12)*(z_2 - z_1)/pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2), 7)) - (6*pow(sigma, 6)*(z_2 - z_1)/pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2), 4)));;

    if (!n){

        set_forces(atom_1, f_x, f_y, f_z);
        set_forces(atom_2, -f_x, -f_y, -f_z);

    }

    else{

        set_forces_n(atom_1, f_x, f_y, f_z);
        set_forces_n(atom_2, -f_x, -f_y, -f_z);

    }
}

void bond_force(Mol_2 *molecule, bool n){

    double k = molecule->k;

    double x_1 = molecule->atom_1.x;
    double y_1 = molecule->atom_1.y;
    double z_1 = molecule->atom_1.z;

    double x_2 = molecule->atom_2.x;
    double y_2 = molecule->atom_2.y;
    double z_2 = molecule->atom_2.z;

    double l_0 = molecule->l_0;

    double f_x = - k*(x_1 - x_2)*(sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2)) - l_0)/sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));
    double f_y = - k*(y_1 - y_2)*(sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2)) - l_0)/sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));
    double f_z = - k*(z_1 - z_2)*(sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2)) - l_0)/sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));

    if (!n){

        set_forces(&(molecule->atom_1), f_x, f_y, f_z);
        set_forces(&(molecule->atom_2), -f_x, -f_y, -f_z);

    }

    else{

        set_forces_n(&(molecule->atom_1), f_x, f_y, f_z);
        set_forces_n(&(molecule->atom_2), -f_x, -f_y, -f_z);

    }
}

void wall_force(Sphere *atom, bool n, double d){

    double epsilon =  6.57E-22; // 4.184 * 0.2385;  kJ/mol
    double sigma = 3.17E-10;//3.4E-10;

    double x_min = -0.1*d;
    double x_max = 1.1*d;

    double y_min = -0.1*d;
    double y_max = 1.1*d;

    double z_min = -0.1*d;
    double z_max = 1.1*d;

    double x = atom->x;
    double y = atom->y;
    double z = atom->z;

    double f_x = 0;
    double f_y = 0;
    double f_z = 0;

    //left/right wall
    f_x += 48*epsilon*pow(sigma, 12)/pow(std::abs(x - x_min) , 13);
    f_x -= 48*epsilon*pow(sigma, 12)/pow(std::abs(x_max - x), 13);

    //top/bottom wall
    f_y -= 48*epsilon*pow(sigma, 12)/pow(std::abs(y_max - y), 13);
    f_y += 48*epsilon*pow(sigma, 12)/pow(std::abs(y - y_min), 13);

    //front/back wall;
    f_z -= 48*epsilon*pow(sigma, 12)/pow(std::abs(z_max - z), 13);
    f_z += 48*epsilon*pow(sigma, 12)/pow(std::abs(z - z_min), 13);

    if (!n){

        set_forces(atom, f_x, f_y, f_z);

    }

    else{

        set_forces_n(atom, f_x, f_y, f_z);

    }


}

void gravity(Sphere *atom, bool n){

    double g = 9.81;

    double x = atom->x;
    double y = atom->y;
    double z = atom->z;

    double f_x = 0;
    double f_y = 0;
    double f_z = 0;

    //top/bottom wall
    f_y = -(atom->m)*g;

    if (!n){

        set_forces(atom, f_x, f_y, f_z);

    }

    else{

        set_forces_n(atom, f_x, f_y, f_z);

    }


}
