#ifndef VELOCITY_VERLET_HPP
#define VELOCITY_VERLET_HPP

#include "../box/box.hpp"

void compute_forces(Box *domain, bool n);
void next_position(Box *domain, double dt);
void next_velocity(Box *domain, double dt);

void verlet_integrate(Box *domain, double dt, double t);


#endif
