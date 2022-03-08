#include <iostream>

#include "../box/box.hpp"
#include "../post_processing/post_processing.hpp"
#include "../rand/rand.hpp"
#include "../velocity_verlet/velocity_verlet.hpp"



int main(int argc, char *argv[]){

    std::string path = "results/frame_0.csv";
    std::string path_2 = "results/frame_0.pdb";

    std::string path_new;

    delete_frame(path);

    double dt = 1E-14;
    double t = 0;
    double t_end = 1E-8;

    int it = 1;

    double d = 2.06E-8;

    d = 2*7*4.32E-10;

    Box* domain = create_box(d/2, d, d, d, 3*3*3);












    box_csv(domain, path);
    box_pdb(domain, path_2);



    while (t < t_end){

        path = "results/frame_" + std::to_string(it) + ".csv";

        verlet_integrate(domain, dt, t);

        if (!(it % 40)){

            //box_pdb(domain, path_2);
            box_csv(domain, path);
            std::cout << "Frame: " << it << std::endl;

        }


        t += dt;
        it += 1;

    }




    destroy_box(domain);

    return 0;
}
