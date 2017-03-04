/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#include "kmc_restart.h"
//#include <windows.h>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

	Restart::Restart() {}

	void Restart::initialize(SimulationBox* const &kmc_box) {
		box  = kmc_box;
//		create_folder_for_restart();
	}
/*
	void Restart::create_folder_for_restart() {
		if (CreateDirectory(folder_name.c_str(), NULL)) {
		} else if (ERROR_ALREADY_EXISTS == GetLastError()) {
		} else {
		    cout << "can not create dump_restart folder for some reason\n";
		}
	}
*/
	void Restart::update(const unsigned long long int &step, const double &current_time) {
		if (step % box->input->period_restart == 0) {
			produce_restart(step, current_time);
		}
	}

	void Restart::produce_restart(const unsigned long long int &step, const double &current_time) {
		open_restart(step);
		write_restart(current_time);
		close_restart();
	}

	string Restart::pad_zeros(const unsigned long long int &step) {
		ostringstream step_index;
		step_index << setw(9) << setfill('0') << to_string(step);
		return step_index.str();
	}

	void Restart::open_restart(const unsigned long long int &step) {
		string output_name = "./" + folder_name + "/" + "restart.kmc_simulation." + pad_zeros(step);
  		restart.open(output_name);
	}

	void Restart::write_restart(const double &current_time) {
		for (size_t i = 0; i < box->solutes.size(); ++i) {
			restart << box->solutes[i].type << "\t" << box->solutes[i].id << "\t";
			restart << box->solutes[i].curr_location->id << endl;
		}

		restart << current_time << endl;
	}

	void Restart::close_restart() {
		restart.close();
	}
