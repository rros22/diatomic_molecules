#include "box.hpp"
#include "../rand/rand.hpp"
#include <iostream>
#include <cmath>
#include "../energy/energy.hpp"

#define PI 3.1415926535

Box* create_box(double o, double x, double y, double z, int molecule_no){

    //create box
    Box *box = (Box *) calloc(1, sizeof(Box));

    //memory allocation for molecules
    Mol_2 *molecules = (Mol_2*) calloc(molecule_no, sizeof(Mol_2));

    //initialise wrapper box struct members
    box->o = o;

    box->x = x;
    box->y = y;
    box->z = z;

    box->molecule_no = molecule_no;
    box->molecules = molecules;

    int side = (double) cbrt(molecule_no);



    //initialise molecule positions randomly
    for (int i = 0; i < side; i++){

            for (int j = 0; j < side; j++){

                for (int k = 0; k < side; k++){

                    (box->molecules)[i + side*(j + side*k)].x = x/side + i*x/side;  //hypercube sampling
                    (box->molecules)[i + side*(j + side*k)].y = y/side + j*y/side;
                    (box->molecules)[i + side*(j + side*k)].z = z/side + k*z/side;

                    (box->molecules)[i + side*(j + side*k)].y_rot = rand_range(0, 2*PI);
                    (box->molecules)[i + side*(j + side*k)].z_rot = rand_range(0, 2*PI);

                }

            }
    }



    //initialise molecule velocities randomly
    for (int i = 0; i < molecule_no; i++){

        (box->molecules)[i].u = rand_range(-1000, 1000);
        (box->molecules)[i].v = rand_range(-1000, 1000);
        (box->molecules)[i].w = rand_range(-1000, 1000);

        (box->molecules)[i].yy_rot = rand_range(-6E12, 6E12);
        (box->molecules)[i].zz_rot = rand_range(-6E12, 6E12);
    }

    //initialise bond lengths randomly
    for (int i = 0; i < molecule_no; i++){

        (box->molecules)[i].l_0 = 1.09E-10;
        (box->molecules)[i].l = 1.15E-10;
    }

    //initialise stiffness
    for (int i = 0; i < molecule_no; i++){

        (box->molecules)[i].k = 2287;
    }

    //calculate atom parameters from molecules
    for (int i = 0; i < molecule_no; i++){

        //parameters for simplification
        double x_m = (box->molecules)[i].x;
        double y_m = (box->molecules)[i].y;
        double z_m = (box->molecules)[i].z;

        double y_rot = (box->molecules)[i].y_rot;
        double z_rot = (box->molecules)[i].z_rot;

        double u_m = (box->molecules)[i].u;
        double v_m = (box->molecules)[i].v;
        double w_m = (box->molecules)[i].w;

        double yy_rot = (box->molecules)[i].yy_rot;
        double zz_rot = (box->molecules)[i].zz_rot;

        double l = (box->molecules)[i].l;

        //positions
        (box->molecules)[i].atom_1.x = x_m - 0.5*l*cos(y_rot)*cos(z_rot);
        (box->molecules)[i].atom_1.y = y_m - 0.5*l*sin(z_rot);
        (box->molecules)[i].atom_1.z = z_m + 0.5*l*sin(y_rot)*cos(z_rot);

        (box->molecules)[i].atom_2.x = x_m + 0.5*l*cos(y_rot)*cos(z_rot);
        (box->molecules)[i].atom_2.y = y_m + 0.5*l*sin(z_rot);
        (box->molecules)[i].atom_2.z = z_m - 0.5*l*sin(y_rot)*cos(z_rot);

        //velocities
        (box->molecules)[i].atom_1.u = u_m + l/2*(cos(y_rot)*sin(z_rot)*zz_rot + cos(z_rot)*sin(y_rot)*yy_rot);
        (box->molecules)[i].atom_1.v = v_m - l/2*cos(z_rot)*zz_rot;
        (box->molecules)[i].atom_1.w = w_m + l/2*(cos(z_rot)*cos(y_rot)*yy_rot - sin(z_rot)*sin(y_rot)*zz_rot);

        (box->molecules)[i].atom_2.u = u_m - l/2*(cos(y_rot)*sin(z_rot)*zz_rot + cos(z_rot)*sin(y_rot)*yy_rot);
        (box->molecules)[i].atom_2.v = v_m + l/2*cos(z_rot)*zz_rot;
        (box->molecules)[i].atom_2.w = w_m - l/2*(cos(z_rot)*cos(y_rot)*yy_rot - sin(z_rot)*sin(y_rot)*zz_rot);

        //other parameters
        (box->molecules)[i].atom_1.m = 14.0067*1.66054E-27;
        (box->molecules)[i].atom_2.m = 14.0067*1.66054E-27;


    }


    //intitalise sphere forces to 0
    for (int i = 0; i < molecule_no; i++){

        (box->molecules)[i].atom_1.f_x = 0;
        (box->molecules)[i].atom_1.f_y = 0;
        (box->molecules)[i].atom_1.f_z = 0;

        (box->molecules)[i].atom_1.f_x_n = 0;
        (box->molecules)[i].atom_1.f_y_n = 0;
        (box->molecules)[i].atom_1.f_z_n = 0;

        (box->molecules)[i].atom_2.f_x = 0;
        (box->molecules)[i].atom_2.f_y = 0;
        (box->molecules)[i].atom_2.f_z = 0;

        (box->molecules)[i].atom_2.f_x_n = 0;
        (box->molecules)[i].atom_2.f_y_n = 0;
        (box->molecules)[i].atom_2.f_z_n = 0;

    }

    //initialise energies

    for (int i = 0; i < molecule_no; i++){

      calculate_energies(box->molecules, i);

    }

    return box;
}


void destroy_box(Box *box){

    free(box->molecules);
    free(box);
}

//Sphere functions

/*

void update_forces(Box *domain){

    Sphere *spheres = domain->spheres;
    int sphere_no = domain->sphere_no;

    for (int i = 0; i < sphere_no; i++){

        spheres[i].f_x = spheres[i].f_x_n;
        spheres[i].f_y = spheres[i].f_y_n;
        spheres[i].f_z = spheres[i].f_z_n;
    }

}

*/
