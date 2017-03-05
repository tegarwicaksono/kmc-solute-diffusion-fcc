/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#ifndef KMC_LOGFILE_H_INCLUDED
#define KMC_LOGFILE_H_INCLUDED

#include "kmc_simulationbox.h"
#include <vector>
#include <string>
#include <fstream>

using namespace std;

class Logfile {
public:
    SimulationBox* kmc_box;

    vector<vector<vector<int> > > neighbour_list_nn;
	ofstream 			logfile;
	string folder_name = "log";

    Logfile() {}
    void assign_box(SimulationBox* const &sb);
//    void create_folder_for_logfile();
    void create_logfile_header();
    void start_logfile(SimulationBox* const &sb);
    void update(const unsigned long long int &step, const double &current_time);

    void build_neighbour_list_nn(const size_t &index);
    void build_neighbour_list();
    void print_neighbour_list();
    double compute_total_energy();
};


#endif // KMC_LOGFILE_H_INCLUDED
