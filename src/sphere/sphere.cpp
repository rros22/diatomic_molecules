#include "sphere.hpp"
#include <iostream>
#include <cmath>

//test functions
Sphere* create_sphere(double x, double y, double z, double u, double v, double w){

    Sphere *sphere = (Sphere*) calloc(1, sizeof(Sphere));

    sphere->x = x;
    sphere->y = y;
    sphere->z = z;

    sphere->u = u;
    sphere->v = v;
    sphere->w = w;

    return sphere;

}

void destroy_sphere(Sphere *sphere){

    free(sphere);

}

void set_coordinates(Sphere *sphere, double x, double y, double z){

    sphere->x = x;
    sphere->y = y;
    sphere->z = z;

}

void set_velocity(Sphere *sphere, double u, double v, double w){

    sphere->u = u;
    sphere->v = v;
    sphere->w = w;

}

double distance_spheres(Sphere *sphere_a, Sphere *sphere_b){

    double r[3];

    r[0] = sphere_b->x - sphere_a->x;
    r[1] = sphere_b->y - sphere_a->y;
    r[2] = sphere_b->z - sphere_a->z;

    return sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));

}
