/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#include "kmc_kmcsimulation.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <functional>
#include <chrono>
#include <random>

using namespace std;

	KMCSimulation::KMCSimulation() {}

	void KMCSimulation::print_rates() {
		for (size_t i = 0; i < rates.size(); ++i) {
			rates[i].print();
		}
	}

	void KMCSimulation::assign_simulation_box(SimulationBox* const &sb) {
		kmc_box = sb;

		//cout << "There are " << sb->solutes.size() << " number of solutes\n";
		for (size_t i = 0; i < kmc_box->solutes.size(); ++i) {
			for (size_t j = 0; j < kmc_box->solutes[i].rates.size(); ++j) {
				rates.push_back(Rate(&kmc_box->solutes[i], j));
			}
		}

		if (kmc_box->input->dump_snapshot) {
			prepare_for_snapshot();
		}

		if (kmc_box->input->dump_restart) {
			prepare_for_restart();
		}

		curr_time_clock = kmc_box->input->start_time;
		prev_time_clock = curr_time_clock;

        //print_interaction_energy();

	}

    void KMCSimulation::print_interaction_energy() {
        for (size_t i = 0; i < kmc_box->input->e_species.size(); ++i) {
            for (size_t j = 0; j < kmc_box->input->e_species[i].size(); ++j) {
                cout << "Interaction energy for solute " << i << "-" << j << " is : " << endl;
                cout << "  at 1st nearest ngb distance = " << kmc_box->input->e_species[i][j][1] << "\n";
                cout << "  at 2nd nearest ngb distance = " << kmc_box->input->e_species[i][j][2] << "\n";
                cout << "  at 3rd nearest ngb distance = " << kmc_box->input->e_species[i][j][3] << "\n";

            }
        }
    }

	void KMCSimulation::prepare_for_snapshot() {
		dump.initialize(kmc_box);
	}

	void KMCSimulation::prepare_for_restart() {
		restart.initialize(kmc_box);
	}

	void KMCSimulation::calculate_total_rate() {
		total_rate = 0.0;
		int count = 0;
		for (size_t i = 0; i < kmc_box->solutes.size(); ++i) {
			for (size_t j = 0; j < kmc_box->solutes[i].rates.size(); ++j) {
				total_rate += kmc_box->solutes[i].rates[j];
				rates[count++].cumulative = total_rate;
			}
		}

		if (total_rate == 0.0) {
            cout << "==============================" << endl;
            cout << "!!!!!!!!!!!!WARNING!!!!!!!!!!!" << endl;
            cout << "!!EACH EVENT HAS A ZERO RATE!!" << endl;
            cout << "!KMC DOES NOT KNOW WHAT TO DO!" << endl;
            cout << "SUGGESTION:  MODIFY INPUT FILE" << endl;
            cout << "==============================" << endl;
		}
	}

	Rate KMCSimulation::find_rate(const double &rate) {
		auto itr  = rates.begin();
		while (itr != rates.end()) {
			if (rate < itr->cumulative) {
				return *itr;
			} else {
				++itr;
			}
		}
		return *itr;
	}

	double KMCSimulation::draw_random_number() {
		auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		auto real_rand = std::bind(std::uniform_real_distribution<double>(0,1),
		                           mt19937(seed));
		return (1.0 - real_rand());
	}

	void KMCSimulation::choose_event() {
		chosen_cumulative_rate = draw_random_number()*total_rate;
		chosen_event.specify(find_rate(chosen_cumulative_rate));
	}

	void KMCSimulation::execute_event() {
		int id_of_current_site = chosen_event.species->curr_location->id;
		int id_of_future_site = chosen_event.species->next_locations[chosen_event.direction]->id;

		int current_occupant = chosen_event.species->id;
		int future_occupant = kmc_box->lattice_sites[id_of_future_site].occupant;
		//cout << "before swap, lattice id = " << id_of_future_site << " has occupant = " << future_occupant << "\n";

		if (future_occupant >= 0) {
			MovingSpecies* affected_solute = &kmc_box->solutes[future_occupant];
			chosen_event.species->update_location(chosen_event.direction);
			chosen_event.species->update_next_locations();
			affected_solute->update_location(chosen_event.species->prev_location);
			affected_solute->update_next_locations();
		} else {
			chosen_event.species->update_location(chosen_event.direction);
			chosen_event.species->update_next_locations();

		}

		kmc_box->lattice_sites[id_of_current_site].occupant = future_occupant;
		kmc_box->lattice_sites[id_of_future_site].occupant  = current_occupant;

		//cout << "after swap, lattice id = " << id_of_future_site << " has occupant = " << kmc_box->lattice_sites[id_of_future_site].occupant << endl;
	}

	void KMCSimulation::update_simulation() {
		kmc_box->update_solute_energy_and_rate();
	}

	void KMCSimulation::increment_time() {
		chosen_time_increment = (-log(draw_random_number()) / total_rate);
		prev_time_clock  = curr_time_clock;
		curr_time_clock += chosen_time_increment;
	}

	void KMCSimulation::run_one_kmc_step() {
		calculate_total_rate();
							//print_total_rate();
		choose_event();
							//print_chosen_event();
		execute_event();
							//print_execute_event();
		update_simulation();
		increment_time();
							//print_current_time();
	}

	void KMCSimulation::print_total_rate() {
		cout << " total_rate = " << total_rate << "\n";
	}

	void KMCSimulation::print_chosen_event() {
		cout << " chosen cumul rate = " << chosen_cumulative_rate << "\n";
		cout << "  this belongs to species #" << chosen_event.species->id << ", and direction of " << chosen_event.direction << endl;
	}

	void KMCSimulation::print_execute_event() {
		cout << " after executed, species #" << chosen_event.species->id << " is now at site " << chosen_event.species->curr_location-> id << "\n";
		cout << " its future sites are: \n";
		for (size_t i = 0; i < chosen_event.species->next_locations.size(); ++i) {
			cout << "   site id #" << chosen_event.species->next_locations[i]->id << "\n";
		}
	}

	void KMCSimulation::print_current_time() {
		cout << "increment time = " << chosen_time_increment << ", current simulation time = " << curr_time_clock << " s\n";
	}

	void KMCSimulation::initialize_logfile() {
        logfile.start_logfile(kmc_box);
		logfile.update(kmc_box->input->initial_timestep, curr_time_clock);
	}

	void KMCSimulation::check_to_produce_restart(const unsigned long long int &step) {
		if (kmc_box->input->dump_restart) {
			restart.update(step, curr_time_clock);
		}
	}

	void KMCSimulation::check_to_produce_logfile(const unsigned long long int &step) {
        logfile.update(step, curr_time_clock);
	}

	void KMCSimulation::check_to_produce_snapshot(const unsigned long long int &step) {
		if (kmc_box->input->dump_snapshot) {
			dump.update(step);
		}
	}

	void KMCSimulation::print_simulation_step_status(const unsigned long long int &step) {
		if (step % kmc_box->input->period_snapshot == 0) {
			cout << "kmc_step" << setw(9) << setfill(' ') << step;
			cout << ", chosen species is solute type " << chosen_event.species->type;
			cout << " with id = " << chosen_event.species->id;
			cout << " current time is = " << curr_time_clock;
			cout << endl;
		}
	}

	void KMCSimulation::print_simulation_initial_status() {
		if (kmc_box->input->initial_timestep == 0ull) {
			cout << "\nkmc simulation is initialized from scratch, with clock starts at " << curr_time_clock << endl << endl;
		} else {
			cout << "\nkmc simulation is initialized from a restart at step " << kmc_box->input->initial_timestep << ", with clock starts at " << curr_time_clock << endl << endl;
		}
	}

	void KMCSimulation::run_kmc_simulation(SimulationBox* const &sb) {
		assign_simulation_box(sb);
		check_to_produce_restart(kmc_box->input->initial_timestep);
		check_to_produce_snapshot(kmc_box->input->initial_timestep);
		initialize_logfile();
		print_simulation_initial_status();
		for (unsigned long long int step = kmc_box->input->initial_timestep + 1; step <= kmc_box->input->final_timestep; ++step) {
			run_one_kmc_step();
			print_simulation_step_status(step);
			if (step < kmc_box->input->final_timestep) {
                check_to_produce_logfile(step);
				check_to_produce_snapshot(step);
				check_to_produce_restart(step);
			}
		}
	}
