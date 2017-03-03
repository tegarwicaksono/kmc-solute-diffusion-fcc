/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#ifndef KMC_SNAPSHOT_H_INCLUDED
#define KMC_SNAPSHOT_H_INCLUDED

#include "kmc_simulationbox.h"
#include "kmc_movingspecies.h"
#include "kmc_latticesite.h"
#include <string>
#include <vector>
#include <fstream>

class Snapshot {
public:
	int  number_of_species_in_snapshot = 0;
	bool matrix_included = false;
	string folder_name = "dump_snapshot";

	SimulationBox* 			box;
	vector<MovingSpecies*> 	species_to_dump;
	vector<LatticeSite*> 	matrix_to_dump;
	vector<string> 			solute_name;
	vector<double>			solute_mass;
	vector<double>			box_dimension;
	ofstream 				dump;

	Snapshot();
	void initialize(SimulationBox* const &kmc_box);
	void update(const unsigned long long int &step);
	void create_folder_for_snapshot();
	void assign_solute_snapshot_properties();
	void tabulate_species_to_dump();
	void produce_snapshot(const unsigned long long int &step);
	string pad_zeros(const unsigned long long int &step);
	void print_header(const unsigned long long int &step);
	void print_species();
	void close_file();
};


#endif // KMC_SNAPSHOT_H_INCLUDED
