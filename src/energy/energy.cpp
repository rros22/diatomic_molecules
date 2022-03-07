#include "energy.hpp"
#include <cmath>

void ke_trans(Mol_2 *molecules, int i){

  double u_1 = molecules[i].atom_1.u;
  double u_2 = molecules[i].atom_2.u;

  double v_1 = molecules[i].atom_1.v;
  double v_2 = molecules[i].atom_2.v;

  double w_1 = molecules[i].atom_1.w;
  double w_2 = molecules[i].atom_2.w;

  double m_1 = molecules[i].atom_1.m;
  double m_2 = molecules[i].atom_2.m;

  double m = m_1;

  molecules[i].ke_trans = m/4*(pow(u_1 + u_2, 2) + pow(v_1 + v_2, 2) + pow(w_1 + w_2, 2));

}

void ke_rot(Mol_2 *molecules, int i){

  double x_1 = molecules[i].atom_1.x;
  double x_2 = molecules[i].atom_2.x;

  double y_1 = molecules[i].atom_1.y;
  double y_2 = molecules[i].atom_2.y;

  double z_1 = molecules[i].atom_1.z;
  double z_2 = molecules[i].atom_2.z;

  double u_1 = molecules[i].atom_1.u;
  double u_2 = molecules[i].atom_2.u;

  double v_1 = molecules[i].atom_1.v;
  double v_2 = molecules[i].atom_2.v;

  double w_1 = molecules[i].atom_1.w;
  double w_2 = molecules[i].atom_2.w;

  double m_1 = molecules[i].atom_1.m;
  double m_2 = molecules[i].atom_2.m;

  double m = m_1;

  //moment of intertia calculation
  double l = pow(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2), 0.5);

  double j = 0.5*m*pow(l, 2);

  double r_2 = 1/4*(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));

  //angular velocities

  double w_x = 1/pow(l/2, 2)/4*((w_2 - w_1)*(y_2 - y_1) - (v_2 - v_1)*(z_2 - z_1));
  double w_y = 1/pow(l/2, 2)/4*((u_2 - u_1)*(z_2 - z_1) - (w_2 - w_1)*(x_2 - x_1));
  double w_z = 1/pow(l/2, 2)/4*((v_2 - v_1)*(x_2 - x_1) - (u_2 - u_1)*(y_2 - y_1));

  double ke_rot = 0.5*j*(pow(w_x, 2) + pow(w_y, 2) + pow(w_z, 2));

  molecules[i].ke_rot = ke_rot;


}

void ke_vibr(Mol_2 *molecules, int i){

  double x_1 = molecules[i].atom_1.x;
  double x_2 = molecules[i].atom_2.x;

  double y_1 = molecules[i].atom_1.y;
  double y_2 = molecules[i].atom_2.y;

  double z_1 = molecules[i].atom_1.z;
  double z_2 = molecules[i].atom_2.z;

  double u_1 = molecules[i].atom_1.u;
  double u_2 = molecules[i].atom_2.u;

  double v_1 = molecules[i].atom_1.v;
  double v_2 = molecules[i].atom_2.v;

  double w_1 = molecules[i].atom_1.w;
  double w_2 = molecules[i].atom_2.w;

  double m_1 = molecules[i].atom_1.m;
  double m_2 = molecules[i].atom_2.m;

  double m = m_1;

  //fraction simplification
  double num = (u_2 - u_1)/2*(x_2 - x_1) + (v_2 - v_1)/2*(y_2 - y_1) + (w_2 - w_1)/2*(z_2 - z_1);
  double den = sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));

  molecules[i].ke_vibr = m*pow(num/den, 2);

}

void pe(Mol_2 *molecules, int i){

  double x_1 = molecules[i].atom_1.x;
  double x_2 = molecules[i].atom_2.x;

  double y_1 = molecules[i].atom_1.y;
  double y_2 = molecules[i].atom_2.y;

  double z_1 = molecules[i].atom_1.z;
  double z_2 = molecules[i].atom_2.z;

  double k = molecules[i].k;
  double l_0 = molecules[i].l_0;

  molecules[i].pe = 0.5*k*pow(sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2)) - l_0, 2);

}

void e_total(Mol_2 *molecules, int i){

  molecules[i].e_total = molecules[i].ke_trans + molecules[i].ke_vibr + molecules[i].ke_rot + molecules[i].pe;

}

void calculate_energies(Mol_2 *molecules, int i){

  ke_trans(molecules, i);
  ke_rot(molecules, i);
  ke_vibr(molecules, i);
  pe(molecules, i);
  e_total(molecules, i);

}
