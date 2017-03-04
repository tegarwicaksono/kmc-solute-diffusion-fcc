/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#ifndef KMC_RESTART_H_INCLUDED
#define KMC_RESTART_H_INCLUDED

#include "kmc_simulationbox.h"
#include <fstream>
#include <string>

using namespace std;

class Restart {
public:
	SimulationBox* 			box;
	ofstream 				restart;
	string folder_name = "dump_restart";

	Restart();

	void initialize(SimulationBox* const &kmc_box);
//	void create_folder_for_restart();
	void update(const unsigned long long int &step, const double &current_time);
	void produce_restart(const unsigned long long int &step, const double &current_time);
	string pad_zeros(const unsigned long long int &step);
	void open_restart(const unsigned long long int &step);
	void write_restart(const double &current_time);
	void close_restart();
};

#endif // KMC_RESTART_H_INCLUDED
