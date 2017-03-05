/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#include "kmc_snapshot.h"
//#include <windows.h>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

	Snapshot::Snapshot() {}

	void Snapshot::initialize(SimulationBox* const &kmc_box) {
		box  = kmc_box;
		//create_folder_for_snapshot();
		assign_solute_snapshot_properties();
		tabulate_species_to_dump();
	}

	void Snapshot::update(const unsigned long long int &step) {
		if (step % box->input->period_snapshot == 0) {
			produce_snapshot(step);
		}
	}

	/*
	void Snapshot::create_folder_for_snapshot() {
		if (CreateDirectory(folder_name.c_str(), NULL)) {
		} else if (ERROR_ALREADY_EXISTS == GetLastError()) {
		} else {
		    cout << "can not create dump_snapshot folder for some reason\n";
		}
	}
	*/

	void Snapshot::assign_solute_snapshot_properties() {
		vector<string> element_list{"Al","Ga","C","N","O","B","Fe","Mg","F","Ar"};

		for (int i = 0; i < box->input->number_of_solute_type + 1; ++i) {
			solute_name.push_back(element_list[i % element_list.size()]);
			solute_mass.push_back(1.0);
		}

	}

	void Snapshot::tabulate_species_to_dump() {
		if (box->input->include_species_in_snapshot[0]) {
			matrix_included = true;
			for (size_t i = 0; i < box->lattice_sites.size(); ++i) {
				matrix_to_dump.push_back(&box->lattice_sites[i]);
			}
			number_of_species_in_snapshot = box->lattice_sites.size();
		} else {
			for (size_t i = 0; i < box->solutes.size(); ++i) {
				if (box->input->include_species_in_snapshot[box->solutes[i].type]) {
					species_to_dump.push_back(&box->solutes[i]);
					++number_of_species_in_snapshot;
				}
			}
		}

		for (size_t i = 0; i < box->input->box_length.size(); ++i) {
			box_dimension.push_back(box->input->box_length[i]*box->input->unit_length[i]);
		}
	}

	void Snapshot::produce_snapshot(const unsigned long long int &step) {
		print_header(step);
		print_species();
		close_file();
	}

	string Snapshot::pad_zeros(const unsigned long long int &step) {
		ostringstream step_index;
		step_index << setw(9) << setfill('0') << to_string(step);
		return step_index.str();
	}

	void Snapshot::print_header(const unsigned long long int &step) {
		string output_name = "./" + folder_name + "/" + "snapshot.kmc_simulation." + pad_zeros(step);
  		dump.open(output_name);

  		dump << "Number of particles = " << number_of_species_in_snapshot << endl;
  		dump << "A = 1 Angstrom" << endl;
  		dump << "H0(1,1) = " << box_dimension[0] << " A" << endl;
  		dump << "H0(1,2) = 0 A" << endl;
  		dump << "H0(1,3) = 0 A" << endl;
  		dump << "H0(2,1) = 0 A" << endl;
  		dump << "H0(2,2) = " << box_dimension[1] << " A" << endl;
  		dump << "H0(2,3) = 0 A" << endl;
  		dump << "H0(3,1) = 0 A" << endl;
  		dump << "H0(3,2) = 0 A" << endl;
  		dump << "H0(3,3) = " << box_dimension[2] << " A" << endl;
  		dump << ".NO_VELOCITY." << endl;
        if (matrix_included) {

            dump << "entry_count = 3" << endl;
        } else {
            dump << "entry_count = 4" << endl;
        }
	}

	void Snapshot::print_species() {
		double kB = 1.38064852e-23;	//Boltzmann constant
		double eV = 1.60217662e-19;	//elementary charge
		if (matrix_included) {
			//cout << "printing species, matrix included\n";
			for (size_t i = 0; i < matrix_to_dump.size(); ++i) {
				int species_type = (matrix_to_dump[i]->occupant < 0) ? 0 : box->solutes[matrix_to_dump[i]->occupant].type;
				dump << solute_mass[species_type] << endl << solute_name[species_type] << endl;
				for (size_t j = 0; j < 3; ++j) {
					dump << matrix_to_dump[i]->xyz[j] / box->input->box_length[j] << " ";
				}
				dump << endl;
			}
		} else {
			for (size_t i = 0; i < species_to_dump.size(); ++i) {
				int species_type = species_to_dump[i]->type;
				dump << solute_mass[species_type] << endl << solute_name[species_type] << endl;

				for (size_t j = 0; j < species_to_dump[i]->curr_location->xyz.size(); ++j) {
					dump << species_to_dump[i]->curr_location->xyz[j] / static_cast<double>(box->input->box_length[j]) << "\t";
				}

				dump << species_to_dump[i]->current_energy *(kB * box->input->abs_temperature) / eV;
				dump << endl;
			}
		}
	}

	void Snapshot::close_file() {
		dump.close();
	}
