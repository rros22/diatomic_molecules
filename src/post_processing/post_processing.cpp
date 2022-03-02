#include "post_processing.hpp"

void labels(std::string path){

    std::ofstream file;
    file.open(path, std::fstream::app);

    file << "X";
    file << ",";

    file << "Y";
    file << ",";

    file << "Z";
    file << ",";

    file << "U";
    file << ",";

    file << "V";
    file << ",";

    file << "W";
    file << ",";

    file << "Fx";
    file << ",";

    file << "Fy";
    file << ",";

    file << "Fz";
    file << ",";

    file << "m";


    file << std::endl;

}



void append_sphere_csv(Sphere *sphere, std::string path){

    std::ofstream file;
    file.open(path, std::fstream::app);

    file << sphere->x;
    file << ",";

    file << sphere->y;
    file << ",";

    file << sphere->z;
    file << ",";

    file << sphere->u;
    file << ",";

    file << sphere->v;
    file << ",";

    file << sphere->w;
    file << ",";

    file << sphere->f_x;
    file << ",";

    file << sphere->f_y;
    file << ",";

    file << sphere->f_z;
    file << ",";

    file << sphere->m;

    file << std::endl;
}

void labels_mol(std::string path){

    std::ofstream file;
    file.open(path, std::fstream::app);

    file << "X_1" << "," << "X_2" << ",";

    file << "Y_1" << "," << "Y_2" << ",";

    file << "Z_1" << "," << "Z_2" << ",";

    file << "U_1" << "," << "U_2" << ",";

    file << "V_1" << "," << "V_2" << ",";

    file << "W_1" << "," << "W_2" << ",";

    file << "m_1" << "," << "m_2" << ",";

    file << "k" << "," << "l_0";

    file << std::endl;

}

void append_molecule_csv(Mol_2 *molecule, std::string path){

    std::ofstream file;
    file.open(path, std::fstream::app);

    //positions
    file << molecule->atom_1.x;
    file << ",";

    file << molecule->atom_2.x;
    file << ",";

    file << molecule->atom_1.y;
    file << ",";

    file << molecule->atom_2.y;
    file << ",";

    file << molecule->atom_1.z;
    file << ",";

    file << molecule->atom_2.z;
    file << ",";

    //velocities
    file << molecule->atom_1.u;
    file << ",";

    file << molecule->atom_2.u;
    file << ",";

    file << molecule->atom_1.v;
    file << ",";

    file << molecule->atom_2.v;
    file << ",";

    file << molecule->atom_1.w;
    file << ",";

    file << molecule->atom_2.w;
    file << ",";

    //atomic mass
    file << molecule->atom_1.m;
    file << ",";

    file << molecule->atom_2.m;
    file << ",";

    //bond properties
    file << molecule->k;
    file << ",";

    file << molecule->l_0;

    file << std::endl;


}

void box_csv(Box *box, std::string path){

    labels_mol(path);

    Mol_2 *molecules = box->molecules;
    int molecule_no = box->molecule_no;

    for (int i = 0; i < molecule_no; i++){

        append_molecule_csv(&(molecules[i]), path);
    }
}

void delete_frame(std::string path){

    const char *cstr = path.c_str();

    if (std::remove(cstr) != 0){

        std::cout << "No frame to be deleted" << std::endl;
    }

    else {

        std::cout << "Frame deleted succesfully" << std::endl;

    }
}

std::string to_angstroms(double meters){

    double temp = meters*1E10;

    std::string str(std::to_string(temp));

    int i = str.find('.');

    std::string str_1(str.substr(0, i));
    std::string str_2(str.substr(i, 4));

    return str_1 + str_2;
}

std::string white_spaces(int spaces){

    std::string str(spaces, ' ');

    return str;

}


void box_pdb(Box *box, std::string path){

    std::ofstream file;
    file.open(path, std::fstream::app);

    file << std::fixed;
    file << std::setprecision(3);

    Mol_2 *molecules = box->molecules;
    int molecule_no = box->molecule_no;

    int j = 1;

    int no_len;
    int name_len;
    int x_len;
    int y_len;
    int z_len;
    int sym_len;

    for (int i = 0; i < molecule_no; i++){

        file << "ATOM";

        //calculate no_length
        no_len = std::to_string(j).size();

        file << white_spaces(11 - 4 - no_len);
        file << std::to_string(j);

        //calculate name_length
        name_len = ("N" + std::to_string(j)).size();



        file << white_spaces(2);
        file << "N" + std::to_string(j);

        //cartesian coordinates
        x_len = to_angstroms(molecules[i].atom_1.x).size();
        y_len = to_angstroms(molecules[i].atom_1.y).size();
        z_len = to_angstroms(molecules[i].atom_1.z).size();

        file << white_spaces(38 - 13 - x_len - name_len);
        file << to_angstroms(molecules[i].atom_1.x);

        file << white_spaces(46 - 38 - y_len);
        file << to_angstroms(molecules[i].atom_1.y);

        file << white_spaces(54 - 46 - z_len);
        file << to_angstroms(molecules[i].atom_1.z);

        //element symbol
        sym_len = 1;

        file << white_spaces(78 - 54 - sym_len);
        file << "N";

        j += 1;

        file << std::endl;

        file << "ATOM";

        //calculate no_length
        no_len = std::to_string(j).size();

        file << white_spaces(11 - 4 - no_len);
        file << std::to_string(j);

        //calculate name_length
        name_len = ("N" + std::to_string(j)).size();



        file << white_spaces(2);
        file << "N" + std::to_string(j);

        //cartesian coordinates
        x_len = to_angstroms(molecules[i].atom_2.x).size();
        y_len = to_angstroms(molecules[i].atom_2.y).size();
        z_len = to_angstroms(molecules[i].atom_2.z).size();

        file << white_spaces(38 - 13 - x_len - name_len);
        file << to_angstroms(molecules[i].atom_2.x);

        file << white_spaces(46 - 38 - y_len);
        file << to_angstroms(molecules[i].atom_2.y);

        file << white_spaces(54 - 46 - z_len);
        file << to_angstroms(molecules[i].atom_2.z);

        //element symbol
        sym_len = 1;

        file << white_spaces(78 - 54 - sym_len);
        file << "N";

        file << std::endl;

        j += 1;
    }

    //new frame
    file << "ENDMDL";

    file << std::endl;
}
