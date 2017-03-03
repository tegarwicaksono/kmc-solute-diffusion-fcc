/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#ifndef KMC_SIMULATIONBOX_H_INCLUDED
#define KMC_SIMULATIONBOX_H_INCLUDED

#include "kmc_latticesite.h"
#include "kmc_movingspecies.h"
#include "kmc_inputdata.h"

using namespace std;

class SimulationBox {
public:
	vector<LatticeSite>   lattice_sites;
	vector<MovingSpecies> solutes;
	unordered_set<int> 	  occupied_sites;
	int flag_for_identical_swap = 1;

	InputData* input;

	SimulationBox();
	void generate_box(InputData* const &id);
	LatticeSite* find_latt_id_fcc(const vector<double> &xyz);
	vector<vector<double>> atoms_in_fcc_unit_cell();
	vector<vector<double>> first_nearest_neighbours_fcc();
	vector<vector<double>> second_nearest_neighbours_fcc();
	vector<vector<double>> third_nearest_neighbours_fcc();

	void generate_sites_fcc();
	vector<vector<double>> find_nn_per_site(const LatticeSite &site, const vector<vector<double>> &neighbours);

	void assign_1nn_fcc_per_site(LatticeSite &site, const vector<vector<double>> &neighbours);
	void assign_2nn_fcc_per_site(LatticeSite &site, const vector<vector<double>> &neighbours);
	void assign_3nn_fcc_per_site(LatticeSite &site, const vector<vector<double>> &neighbours);

	void generate_first_nearest_neighbours();
	void generate_second_nearest_neighbours();
	void generate_third_nearest_neighbours();

	void generate_neighbours();
	void print_sites_fcc();
	void print_test_1nn_fcc();
	void print_1nn_fcc();
	void print_test_2nn_fcc();
	void generate_solute_location_from_scratch();
	void generate_solute_location_from_restart();
	double calculate_energy(LatticeSite* const &site, const int &ref_type);
	double calculate_energy(LatticeSite* const &target_site, LatticeSite* const &current_site, const int &ref_type, bool &flag);
	void generate_solute_energy_and_rate();
	void update_solute_energy_and_rate();
	void generate_solute_from_scratch();
	void generate_solute_from_restart();
	void print_solute();
	void print_solute_future_energy();

};

#endif // KMC_SIMULATIONBOX_H_INCLUDED
