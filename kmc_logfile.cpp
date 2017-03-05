/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#include "kmc_logfile.h"
#include "kmc_simulationbox.h"
//#include <windows.h>
#include <sstream>
#include <iostream>

using namespace std;

    void Logfile::assign_box(SimulationBox* const &sb) {
        kmc_box = sb;
    }
/*
    void Logfile::create_folder_for_logfile() {
 		if (CreateDirectory(folder_name.c_str(), NULL)) {
		} else if (ERROR_ALREADY_EXISTS == GetLastError()) {
		} else {
		    cout << "can not create log folder for some reason\n";
		}
    }
*/
    void Logfile::create_logfile_header() {
		string output_name = "./" + folder_name + "/" + "log.kmc_simulation.from.step" + to_string(kmc_box->input->initial_timestep) + ".to.step" + to_string(kmc_box->input->final_timestep);
  		logfile.open(output_name);

  		logfile << "KMC_timestep\tTimeClock[second]\tTotalEnergy[eV]\n";
    }

    void Logfile::start_logfile(SimulationBox* const &sb) {
        assign_box(sb);
//        create_folder_for_logfile();
        create_logfile_header();
        build_neighbour_list();
    }

    void Logfile::update(const unsigned long long int &step, const double &current_time) {
        if (step % kmc_box->input->period_snapshot == 0) {
            logfile << step << "\t" << current_time << "\t" << compute_total_energy() << "\n";
        }
    }

    void Logfile::build_neighbour_list_nn(const size_t &index) {
        for (size_t i = 0; i < kmc_box->lattice_sites.size(); ++i) {
            vector<int> neighbour_per_site;
            int main_id = kmc_box->lattice_sites[i].id;
            for (size_t j = 0; j < kmc_box->lattice_sites[i].nth_neighbours[index].size(); ++j) {
                if (main_id < kmc_box->lattice_sites[i].nth_neighbours[index][j]->id) {
                    neighbour_per_site.push_back(kmc_box->lattice_sites[i].nth_neighbours[index][j]->id);
                }
            }
            neighbour_list_nn[index].push_back(neighbour_per_site);
        }
    }

    void Logfile::build_neighbour_list() {
        neighbour_list_nn.assign(kmc_box->input->eff_ngb_distance + 1, vector<vector<int> >());

        for (size_t i = 1; i < kmc_box->input->include_ngb.size(); ++i) {
            if (kmc_box->input->include_ngb[i]) {
                build_neighbour_list_nn(i);
            }
        }
    }

    void Logfile::print_neighbour_list() {
        for (size_t k = 1; k < kmc_box->input->include_ngb.size(); ++k) {
            int bond_count = 0;
            if (kmc_box->input->include_ngb[k]) {
                cout << "Pairs of " << k << "-th nearest ngb are: \n";
                for (size_t i = 0; i < neighbour_list_nn[k].size(); ++i) {
                    for (size_t j = 0; j < neighbour_list_nn[k][i].size(); ++j) {
                        cout << "bond #" << bond_count + 1 << " = " << i << "-" << neighbour_list_nn[k][i][j] << endl;
                        ++bond_count;
                    }
                }
            }
        }
    }

    double Logfile::compute_total_energy() {
        double total_energy = 0.0;

        for (size_t k = 1; k < kmc_box->input->include_ngb.size(); ++k) {
            int bond_count = 0;
            if (kmc_box->input->include_ngb[k]) {
                for (size_t i = 0; i < neighbour_list_nn[k].size(); ++i) {
                    for (size_t j = 0; j < neighbour_list_nn[k][i].size(); ++j) {
                        int sol1 = (kmc_box->lattice_sites[i].occupant < 0) ? 0 : kmc_box->solutes[kmc_box->lattice_sites[i].occupant].type;
                        int sol2 = (kmc_box->lattice_sites[neighbour_list_nn[k][i][j]].occupant < 0) ? 0 : kmc_box->solutes[kmc_box->lattice_sites[neighbour_list_nn[k][i][j]].occupant].type;

                        total_energy += kmc_box->input->e_species[sol1][sol2][k];
                        ++bond_count;
                    }
                }
            }
        }

        return total_energy* (kmc_box->input->kB * kmc_box->input->abs_temperature) / kmc_box->input->eV;
    }
