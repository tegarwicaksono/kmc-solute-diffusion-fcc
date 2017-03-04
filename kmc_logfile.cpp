/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
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
  		log.open(output_name);

  		log << "KMC_timestep\tTimeClock[second]\tTotalEnergy[eV]\n";
    }

    void Logfile::start_logfile(SimulationBox* const &sb) {
        assign_box(sb);
//        create_folder_for_logfile();
        create_logfile_header();
        build_neighbour_list();
    }

    void Logfile::update(const unsigned long long int &step, const double &current_time) {
        if (step % kmc_box->input->period_snapshot == 0) {
            log << step << "\t" << current_time << "\t" << compute_total_energy() << "\n";
        }
    }

    void Logfile::build_neighbour_list_1nn() {
        for (size_t i = 0; i < kmc_box->lattice_sites.size(); ++i) {
            vector<int> neighbour_per_site;
            int main_id = kmc_box->lattice_sites[i].id;
            for (size_t j = 0; j < kmc_box->lattice_sites[i].first_nn.size(); ++j) {
                if (main_id < kmc_box->lattice_sites[i].first_nn[j]->id) {
                    neighbour_per_site.push_back(kmc_box->lattice_sites[i].first_nn[j]->id);
                }
            }
            neighbour_list_1nn.push_back(neighbour_per_site);
        }
    }

    void Logfile::build_neighbour_list_2nn() {
        for (size_t i = 0; i < kmc_box->lattice_sites.size(); ++i) {
            vector<int> neighbour_per_site;
            int main_id = kmc_box->lattice_sites[i].id;
            for (size_t j = 0; j < kmc_box->lattice_sites[i].second_nn.size(); ++j) {
                if (main_id < kmc_box->lattice_sites[i].second_nn[j]->id) {
                    neighbour_per_site.push_back(kmc_box->lattice_sites[i].second_nn[j]->id);
                }
            }
            neighbour_list_2nn.push_back(neighbour_per_site);
        }
    }

    void Logfile::build_neighbour_list_3nn() {
        for (size_t i = 0; i < kmc_box->lattice_sites.size(); ++i) {
            vector<int> neighbour_per_site;
            int main_id = kmc_box->lattice_sites[i].id;
            for (size_t j = 0; j < kmc_box->lattice_sites[i].third_nn.size(); ++j) {
                if (main_id < kmc_box->lattice_sites[i].third_nn[j]->id) {
                    neighbour_per_site.push_back(kmc_box->lattice_sites[i].third_nn[j]->id);
                }
            }
            neighbour_list_3nn.push_back(neighbour_per_site);
        }
    }

    void Logfile::build_neighbour_list() {
        if (kmc_box->input->include_ngb[1]) {
            build_neighbour_list_1nn();
        }

        if (kmc_box->input->include_ngb[2]) {
            build_neighbour_list_2nn();
        }

        if (kmc_box->input->include_ngb[3]) {
            build_neighbour_list_3nn();
        }
    }

    void Logfile::print_neighbour_list() {
        int bond_count = 0;

        if (kmc_box->input->include_ngb[1]) {
            cout << "Pairs of 1st nearest ngb are: \n";
            for (size_t i = 0; i < neighbour_list_1nn.size(); ++i) {
                for (size_t j = 0; j < neighbour_list_1nn[i].size(); ++j) {
                    cout << "bond #" << bond_count + 1<< " = " << i << "-" << neighbour_list_1nn[i][j] << endl;
                    ++bond_count;
                }
            }
        }

        if (kmc_box->input->include_ngb[2]) {
            bond_count = 0;
            cout << "Pairs of 2nd nearest ngb are: \n";
            for (size_t i = 0; i < neighbour_list_2nn.size(); ++i) {
                for (size_t j = 0; j < neighbour_list_2nn[i].size(); ++j) {
                    cout << "bond #" << bond_count + 1<< " = " << i << "-" << neighbour_list_2nn[i][j] << endl;
                    ++bond_count;
                }
            }
        }

        if (kmc_box->input->include_ngb[3]) {
            bond_count = 0;
            cout << "Pairs of 3rd nearest ngb are: \n";
            for (size_t i = 0; i < neighbour_list_3nn.size(); ++i) {
                for (size_t j = 0; j < neighbour_list_3nn[i].size(); ++j) {
                    cout << "bond #" << bond_count + 1<< " = " << i << "-" << neighbour_list_3nn[i][j] << endl;
                    ++bond_count;
                }
            }
        }
    }

    double Logfile::compute_total_energy() {
        double total_energy = 0.0;

        int bond_count = 0;
        for (size_t i = 0; i < neighbour_list_1nn.size(); ++i) {
            for (size_t j = 0; j < neighbour_list_1nn[i].size(); ++j) {
                int sol1 = (kmc_box->lattice_sites[i].occupant < 0) ? 0 : kmc_box->solutes[kmc_box->lattice_sites[i].occupant].type;
                int sol2 = (kmc_box->lattice_sites[neighbour_list_1nn[i][j]].occupant < 0) ? 0 : kmc_box->solutes[kmc_box->lattice_sites[neighbour_list_1nn[i][j]].occupant].type;
                //if (kmc_box->input->e_species[sol1][sol2][1] > 0.0) {
                //    cout << "bond number #" << bond_count << ", between site " << i << " and " << neighbour_list[i][j] << endl;
                //}
                total_energy += kmc_box->input->e_species[sol1][sol2][1];
                ++bond_count;
            }
        }

        if (kmc_box->input->include_ngb[2]) {
            for (size_t i = 0; i < neighbour_list_2nn.size(); ++i) {
                for (size_t j = 0; j < neighbour_list_2nn[i].size(); ++j) {
                    int sol1 = (kmc_box->lattice_sites[i].occupant < 0) ? 0 : kmc_box->solutes[kmc_box->lattice_sites[i].occupant].type;
                    int sol2 = (kmc_box->lattice_sites[neighbour_list_2nn[i][j]].occupant < 0) ? 0 : kmc_box->solutes[kmc_box->lattice_sites[neighbour_list_2nn[i][j]].occupant].type;
                    //if (kmc_box->input->e_species[sol1][sol2][1] > 0.0) {
                    //    cout << "bond number #" << bond_count << ", between site " << i << " and " << neighbour_list[i][j] << endl;
                    //}
                    total_energy += kmc_box->input->e_species[sol1][sol2][2];
                    ++bond_count;
                }
            }
        }

        if (kmc_box->input->include_ngb[3]) {
            for (size_t i = 0; i < neighbour_list_3nn.size(); ++i) {
                for (size_t j = 0; j < neighbour_list_3nn[i].size(); ++j) {
                    int sol1 = (kmc_box->lattice_sites[i].occupant < 0) ? 0 : kmc_box->solutes[kmc_box->lattice_sites[i].occupant].type;
                    int sol2 = (kmc_box->lattice_sites[neighbour_list_3nn[i][j]].occupant < 0) ? 0 : kmc_box->solutes[kmc_box->lattice_sites[neighbour_list_3nn[i][j]].occupant].type;
                    //if (kmc_box->input->e_species[sol1][sol2][1] > 0.0) {
                    //    cout << "bond number #" << bond_count << ", between site " << i << " and " << neighbour_list[i][j] << endl;
                    //}
                    total_energy += kmc_box->input->e_species[sol1][sol2][3];
                    ++bond_count;
                }
            }
        }

        return total_energy;
    }
