#ifndef SPHERE_HPP
#define SPHERE_HPP

typedef struct Sphere{

    double x;
    double y;
    double z;

    double u;
    double v;
    double w;

    double q;

    double m;

    double f_x;
    double f_y;
    double f_z;

    double f_x_n;
    double f_y_n;
    double f_z_n;

} Sphere;


//test functions
Sphere* create_sphere(double x, double y, double z, double u, double v, double w);
void destroy_sphere(Sphere *sphere);

void set_coordinates(Sphere *sphere, double x, double y, double z);
void set_velocity(Sphere *sphere, double u, double v, double w);

double distance_spheres(Sphere *sphere_a, Sphere *sphere_b);


#endif
