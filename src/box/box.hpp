#ifndef BOX_HPP
#define BOX_HPP

#include "../sphere/sphere.hpp"
#include "../diatomic_molecule/diatomic_molecule.hpp"


typedef struct Box {

    //position
    double o;

    //dimensions
    double x;
    double y;
    double z;

    //molecule array pointer
    Mol_2 *molecules;

    int molecule_no;

} Box;

Box* create_box(double o, double x, double y, double z, int molecule_no);
void destroy_box(Box *box);


//Sphere functions



//change forces from old to new
void update_forces(Box *domain);



#endif
