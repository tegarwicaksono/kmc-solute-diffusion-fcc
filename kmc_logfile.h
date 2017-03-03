/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
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
    vector<vector<int>> neighbour_list_1nn;
    vector<vector<int>> neighbour_list_2nn;
    vector<vector<int>> neighbour_list_3nn;
	ofstream 			log;
	string folder_name = "log";

    Logfile() {}
    void assign_box(SimulationBox* const &sb);
    void create_folder_for_logfile();
    void create_logfile_header();
    void start_logfile(SimulationBox* const &sb);
    void update(const unsigned long long int &step, const double &current_time);
    void build_neighbour_list_1nn();
    void build_neighbour_list_2nn();
    void build_neighbour_list_3nn();
    void build_neighbour_list();
    void print_neighbour_list();
    double compute_total_energy();
};


#endif // KMC_LOGFILE_H_INCLUDED
