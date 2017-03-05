/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#ifndef KMC_KMCSIMULATION_H_INCLUDED
#define KMC_KMCSIMULATION_H_INCLUDED

#include "kmc_simulationbox.h"
#include "kmc_rate.h"
#include "kmc_chosenevent.h"
#include "kmc_snapshot.h"
#include "kmc_restart.h"
#include "kmc_logfile.h"
#include <vector>

using namespace std;

class KMCSimulation {
public:
	double curr_time_clock;
	double prev_time_clock;
	double total_rate = 0.0;
	double chosen_cumulative_rate;
	double chosen_time_increment;
	vector<Rate> rates;
	SimulationBox* 	kmc_box;
	ChosenEvent 	chosen_event;
	Snapshot 		dump;
	Restart 		restart;
	Logfile         logfile;

	KMCSimulation();

	void print_rates();
	void assign_simulation_box(SimulationBox* const &sb);
	void prepare_for_snapshot();
	void prepare_for_restart();
	void calculate_total_rate();
	Rate find_rate(const double &rate);
	double draw_random_number();
	void choose_event();
	void execute_event();
	void update_simulation();
	void increment_time();
	void run_one_kmc_step();
	void check_to_produce_restart(const unsigned long long int &step);
	void check_to_produce_snapshot(const unsigned long long int &step);
	void check_to_produce_logfile(const unsigned long long int &step);
	void initialize_logfile();
	void run_kmc_simulation(SimulationBox* const &sb);

    void print_interaction_energy();
	void print_total_rate();
	void print_chosen_event();
	void print_execute_event();
	void print_current_time();
	void print_simulation_step_status(const unsigned long long int &step);
	void print_simulation_initial_status();
};

#endif // KMC_KMCSIMULATION_H_INCLUDED
