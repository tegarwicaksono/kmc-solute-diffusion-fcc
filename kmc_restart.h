/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#ifndef KMC_RESTART_H_INCLUDED
#define KMC_RESTART_H_INCLUDED

#include "kmc_simulationbox.h"
#include <fstream>
#include <string>
#include <utility>

using namespace std;

class Restart {
public:
	SimulationBox* 			box;
	ofstream 				restart;
	string folder_name;

	Restart();
	Restart(const Restart& other);
	Restart(Restart&& other);
    Restart& operator= (Restart other);

	virtual ~Restart() = default;
    friend void swap(Restart &a, Restart &b);

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
