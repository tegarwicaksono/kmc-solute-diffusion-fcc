/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#include "kmc_restart.h"
#include "kmc_simulationbox.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <functional>
#include <chrono>
#include <random>

using namespace std;

	SimulationBox::SimulationBox() {}

	void SimulationBox::generate_box(InputData* const &id) {
		input = id;
		generate_sites_fcc();
		generate_neighbours();
		if (input->start_from_restart) {
			generate_solute_from_restart();
		} else {
			generate_solute_from_scratch();
		}
	}

	LatticeSite* SimulationBox::find_latt_id_fcc(const vector<double> &xyz) {
		int base = 0;
		vector<int> xyz_int(3);

		for (int i = 0; i < 3; ++i) {
			xyz_int[i] = static_cast<int>(xyz[i]);
		}

		base += xyz_int[2]*input->box_length[1]*input->box_length[0];
		base += xyz_int[1]*input->box_length[0];
		base += xyz_int[0];
		base *= 4;

		for (int i = 0; i < 3; ++i) {
			if (xyz[i] == static_cast<double>(xyz_int[i])) {
				base += i;
			}
		}

		return &lattice_sites[base];
	}

	vector<vector<double> > SimulationBox::atoms_in_fcc_unit_cell() {
		vector<vector<double> > list;
		list.push_back(std::initializer_list<double>{0.0,0.5,0.5});
		list.push_back(std::initializer_list<double>{0.5,0.0,0.5});
		list.push_back(std::initializer_list<double>{0.5,0.5,0.0});
		list.push_back(std::initializer_list<double>{0.0,0.0,0.0});
		return list;
	}

	vector<vector<double> > SimulationBox::first_nearest_neighbours_fcc() {
		vector<vector<double> > list;
		list.push_back(std::initializer_list<double>{0.0,-0.5,-0.5});	//1
		list.push_back(std::initializer_list<double>{0.0,-0.5,+0.5});	//2
		list.push_back(std::initializer_list<double>{0.0,0.5,-0.5});		//3
		list.push_back(std::initializer_list<double>{0.0,0.5,0.5});		//4
		list.push_back(std::initializer_list<double>{-0.5,0.0,-0.5});	//5
		list.push_back(std::initializer_list<double>{-0.5,0.0,0.5});		//6
		list.push_back(std::initializer_list<double>{0.5,0.0,-0.5});		//7
		list.push_back(std::initializer_list<double>{0.5,0.0,0.5});		//8
		list.push_back(std::initializer_list<double>{-0.5,-0.5,0.0});	//9
		list.push_back(std::initializer_list<double>{-0.5,0.5,0.0});		//10
		list.push_back(std::initializer_list<double>{0.5,-0.5,0.0});		//11
		list.push_back(std::initializer_list<double>{0.5,0.5,0.0});		//12
		return list;
	}

	vector<vector<double> > SimulationBox::second_nearest_neighbours_fcc() {
		vector<vector<double> > list;
		list.push_back(std::initializer_list<double>{-1.0,0.0,0.0});	//1
		list.push_back(std::initializer_list<double>{1.0,0.0,0.0});	//2
		list.push_back(std::initializer_list<double>{0.0,-1.0,0.0});	//3
		list.push_back(std::initializer_list<double>{0.0,1.0,0.0});	//4
		list.push_back(std::initializer_list<double>{0.0,0.0,-1.0});	//5
		list.push_back(std::initializer_list<double>{0.0,0.0,1.0});	//6
		return list;
	}

	vector<vector<double> > SimulationBox::third_nearest_neighbours_fcc() {
        vector<vector<double> > list;
		list.push_back(std::initializer_list<double>{-0.5,-0.5,-1.0});   //1
		list.push_back(std::initializer_list<double>{-0.5,-0.5,1.0});    //2
		list.push_back(std::initializer_list<double>{-0.5,0.5,-1.0});    //3
		list.push_back(std::initializer_list<double>{-0.5,0.5,1.0});     //4
		list.push_back(std::initializer_list<double>{0.5,-0.5,-1.0});    //5
		list.push_back(std::initializer_list<double>{0.5,-0.5,1.0});     //6
		list.push_back(std::initializer_list<double>{0.5,0.5,-1.0});     //7
		list.push_back(std::initializer_list<double>{0.5,0.5,1.0});      //8

		list.push_back(std::initializer_list<double>{-0.5,-1.0,-0.5});   //9
		list.push_back(std::initializer_list<double>{-0.5,-1.0,0.5});    //10
		list.push_back(std::initializer_list<double>{-0.5,1.0,-0.5});    //11
		list.push_back(std::initializer_list<double>{-0.5,1.0,0.5});     //12
		list.push_back(std::initializer_list<double>{0.5,-1.0,-0.5});    //13
		list.push_back(std::initializer_list<double>{0.5,-1.0,0.5});     //14
		list.push_back(std::initializer_list<double>{0.5,1.0,-0.5});     //15
		list.push_back(std::initializer_list<double>{0.5,1.0,0.5});      //16

		list.push_back(std::initializer_list<double>{-1.0,-0.5,-0.5});   //17
		list.push_back(std::initializer_list<double>{-1.0,-0.5,0.5});    //18
		list.push_back(std::initializer_list<double>{-1.0,0.5,-0.5});    //19
		list.push_back(std::initializer_list<double>{-1.0,0.5,0.5});     //20
		list.push_back(std::initializer_list<double>{1.0,-0.5,-0.5});    //21
		list.push_back(std::initializer_list<double>{1.0,-0.5,0.5});     //22
		list.push_back(std::initializer_list<double>{1.0,0.5,-0.5});     //23
		list.push_back(std::initializer_list<double>{1.0,0.5,0.5});      //24

        return list;
	}

	vector<vector<double> > SimulationBox::fourth_nearest_neighbours_fcc() {
        vector<vector<double> > list;
		list.push_back(std::initializer_list<double>{0.0,-1.0,-1.0});	//1
		list.push_back(std::initializer_list<double>{0.0,-1.0,+1.0});	//2
		list.push_back(std::initializer_list<double>{0.0,1.0,-1.0});	//3
		list.push_back(std::initializer_list<double>{0.0,1.0,1.0});		//4

		list.push_back(std::initializer_list<double>{-1.0,0.0,-1.0});	//5
		list.push_back(std::initializer_list<double>{-1.0,0.0,1.0});	//6
		list.push_back(std::initializer_list<double>{1.0,0.0,-1.0});	//7
		list.push_back(std::initializer_list<double>{1.0,0.0,1.0});		//8

		list.push_back(std::initializer_list<double>{-1.0,-1.0,0.0});	//9
		list.push_back(std::initializer_list<double>{-1.0,1.0,0.0});	//10
		list.push_back(std::initializer_list<double>{1.0,-1.0,0.0});	//11
		list.push_back(std::initializer_list<double>{1.0,1.0,0.0});		//12

		return list;
	}

	vector<vector<double > > SimulationBox::fifth_nearest_neighbours_fcc() {
		list.push_back(std::initializer_list<double>{0.0,-0.5,-1.5});	//1
		list.push_back(std::initializer_list<double>{0.0,-0.5,+1.5});	//2
		list.push_back(std::initializer_list<double>{0.0,0.5,-1.5});	//3
		list.push_back(std::initializer_list<double>{0.0,0.5,1.5});		//4

		list.push_back(std::initializer_list<double>{-0.5,0.0,-1.5});	//5
		list.push_back(std::initializer_list<double>{-0.5,0.0,1.5});	//6
		list.push_back(std::initializer_list<double>{0.5,0.0,-1.5});	//7
		list.push_back(std::initializer_list<double>{0.5,0.0,1.5});		//8

		list.push_back(std::initializer_list<double>{-0.5,-1.5,0.0});	//9
		list.push_back(std::initializer_list<double>{-0.5,1.5,0.0});	//10
		list.push_back(std::initializer_list<double>{0.5,-1.5,0.0});	//11
		list.push_back(std::initializer_list<double>{0.5,1.5,0.0});		//12

		list.push_back(std::initializer_list<double>{0.0,-1.5,-0.5});	//13
		list.push_back(std::initializer_list<double>{0.0,-1.5,+0.5});	//14
		list.push_back(std::initializer_list<double>{0.0,1.5,-0.5});	//15
		list.push_back(std::initializer_list<double>{0.0,1.5,0.5});		//16

		list.push_back(std::initializer_list<double>{-1.5,0.0,-0.5});	//17
		list.push_back(std::initializer_list<double>{-1.5,0.0,0.5});	//18
		list.push_back(std::initializer_list<double>{1.5,0.0,-0.5});	//19
		list.push_back(std::initializer_list<double>{1.5,0.0,0.5});		//20

		list.push_back(std::initializer_list<double>{-1.5,-0.5,0.0});	//21
		list.push_back(std::initializer_list<double>{-1.5,0.5,0.0});	//22
		list.push_back(std::initializer_list<double>{1.5,-0.5,0.0});	//23
		list.push_back(std::initializer_list<double>{1.5,0.5,0.0});		//24
	}

	void SimulationBox::generate_sites_fcc() {
		int latt_id = 0;
		vector<vector<double>> fcc = atoms_in_fcc_unit_cell();
		for (int k = 0; k < input->box_length[2]; ++k) {
			for (int j = 0; j < input->box_length[1]; ++j) {
				for (int i = 0; i < input->box_length[0]; ++i) {
					for (auto&& coor : fcc) {
						vector<double> xyz({coor[0] + static_cast<double>(i),
											coor[1] + static_cast<double>(j),
											coor[2] + static_cast<double>(k)});
						lattice_sites.push_back(LatticeSite(latt_id++, xyz));
					}
				}
			}
		}
	}

	vector<vector<double>> SimulationBox::find_nn_per_site(const LatticeSite &site,
		const vector<vector<double>> &neighbours) {

		vector<vector<double>> list;
		for (size_t i = 0; i < neighbours.size(); ++i) {
			vector<double> neighbour_site(site.xyz);

			for (size_t j = 0; j < neighbours[i].size(); ++j) {
				neighbour_site[j] += neighbours[i][j];
				double periodic_image = 0.0;

				if (neighbour_site[j] >= input->box_length[j]) {
					periodic_image = -1.0*input->box_length[j];
				} else if (neighbour_site[j] < 0.0) {
					periodic_image = input->box_length[j];
				}

				neighbour_site[j] += periodic_image;
			}
			list.push_back(neighbour_site);
		}

		return list;
	}

	void SimulationBox::assign_1nn_fcc_per_site(LatticeSite &site, const vector<vector<double>> &neighbours) {
		vector<vector<double>> nn1s = find_nn_per_site(site, neighbours);
		for (const vector<double> &nn1 : nn1s) {
			site.first_nn.push_back(find_latt_id_fcc(nn1));
		}
	}

	void SimulationBox::assign_2nn_fcc_per_site(LatticeSite &site, const vector<vector<double>> &neighbours) {
		vector<vector<double>> nn2s = find_nn_per_site(site, neighbours);
		for (const vector<double> &nn2 : nn2s) {
			site.second_nn.push_back(find_latt_id_fcc(nn2));
		}
	}

    void SimulationBox::assign_3nn_fcc_per_site(LatticeSite &site, const vector<vector<double>> &neighbours) {
		vector<vector<double>> nn3s = find_nn_per_site(site, neighbours);
		for (const vector<double> &nn3 : nn3s) {
			site.third_nn.push_back(find_latt_id_fcc(nn3));
		}
	}

    void SimulationBox::assign_4nn_fcc_per_site(LatticeSite &site, const vector<vector<double>> &neighbours) {
		vector<vector<double>> nn4s = find_nn_per_site(site, neighbours);
		for (const vector<double> &nn4 : nn4s) {
			site.fourth_nn.push_back(find_latt_id_fcc(nn4));
		}
	}

    void SimulationBox::assign_5nn_fcc_per_site(LatticeSite &site, const vector<vector<double>> &neighbours) {
		vector<vector<double>> nn5s = find_nn_per_site(site, neighbours);
		for (const vector<double> &nn5 : nn5s) {
			site.fifth_nn.push_back(find_latt_id_fcc(nn5));
		}
	}

	void SimulationBox::generate_first_nearest_neighbours() {
		vector<vector<double>> nn1s = first_nearest_neighbours_fcc();

		for (size_t i = 0; i < lattice_sites.size(); ++i) {
			assign_1nn_fcc_per_site(lattice_sites[i], nn1s);
		}
	}

	void SimulationBox::generate_second_nearest_neighbours() {
		vector<vector<double>> nn2s = second_nearest_neighbours_fcc();

		for (size_t i = 0; i < lattice_sites.size(); ++i) {
			assign_2nn_fcc_per_site(lattice_sites[i], nn2s);
		}
	}

	void SimulationBox::generate_third_nearest_neighbours() {
		vector<vector<double>> nn3s = third_nearest_neighbours_fcc();

		for (size_t i = 0; i < lattice_sites.size(); ++i) {
			assign_3nn_fcc_per_site(lattice_sites[i], nn3s);
		}
	}

	void SimulationBox::generate_fourth_nearest_neighbours() {
		vector<vector<double> > nn4s = fourth_nearest_neighbours_fcc();

		for (size_t i = 0; i < lattice_sites.size(); ++i) {
			assign_4nn_fcc_per_site(lattice_sites[i], nn4s);
		}
	}

	void SimulationBox::generate_fifth_nearest_neighbours() {
		vector<vector<double>> nn5s = fifth_nearest_neighbours_fcc();

		for (size_t i = 0; i < lattice_sites.size(); ++i) {
			assign_5nn_fcc_per_site(lattice_sites[i], nn5s);
		}
	}

	void SimulationBox::generate_neighbours() {
		generate_first_nearest_neighbours();
		generate_second_nearest_neighbours();
		generate_third_nearest_neighbours();
		generate_fourth_nearest_neighbours();
		generate_fifth_nearest_neighbours();
	}

	void SimulationBox::print_sites_fcc() {
		cout << "There are " << lattice_sites.size() << " fcc sites in the box\n";
		for (auto &&site : lattice_sites) {
			site.print();
		}
	}

	void SimulationBox::print_test_1nn_fcc() {
		int sample_id;

		sample_id = 8;
		cout << "=========== SAMPLE TEST 1NN ===========\n";
		cout << "Reference site: ";
		lattice_sites[sample_id].print();
		cout << "In object sample, the first_nn has " << lattice_sites[sample_id].first_nn.size() << " elements\n";
		for (LatticeSite* const &nn1 : lattice_sites[sample_id].first_nn) {
			nn1->print();
		}
	}

	void SimulationBox::print_1nn_fcc() {
		for (size_t i = 0; i < lattice_sites.size(); ++i) {
			cout << "site id = " << lattice_sites[i].id << ", occup = " << lattice_sites[i].occupant << ", 1nn ngbs = ";
			for (size_t j = 0; j < lattice_sites[i].first_nn.size(); ++j) {
				cout << lattice_sites[i].first_nn[j]->id << ", ";
			}
			cout << "\n";

		}
	}

	void SimulationBox::print_test_2nn_fcc() {
		int sample_id;

		sample_id = 3;
		cout << "=========== SAMPLE TEST 2NN ===========\n";
		cout << "Reference site: ";
		lattice_sites[sample_id].print();
		cout << "In object sample, the second_nn has " << lattice_sites[sample_id].second_nn.size() << " elements\n";
		for (LatticeSite* const &nn2 : lattice_sites[sample_id].second_nn) {
			nn2->print();
		}
	}

	void SimulationBox::generate_solute_location_from_scratch() {
		auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		auto dice_rand = std::bind(std::uniform_int_distribution<int>(0,lattice_sites.size() - 1),
		                           mt19937(seed));

		int count = 0;
		for (int i = 0; i < input->number_of_solute_type; ++i) {
			for (int j = 0; j < input->number_of_solute_per_type[i]; ++j) {
				int site_for_solute = dice_rand();

				while (occupied_sites.find(site_for_solute) != occupied_sites.end()) {
					site_for_solute = dice_rand();
				}

				occupied_sites.insert(site_for_solute);
				lattice_sites[site_for_solute].occupant = count;
				solutes.push_back(MovingSpecies(i + 1, count, &lattice_sites[site_for_solute]));
				solutes[count].update_next_locations();
				++count;
			}
		}

		input->start_time = 0.0;
	}

	void SimulationBox::generate_solute_location_from_restart() {
		ostringstream step_index;
		step_index << setw(9) << setfill('0') << to_string(input->initial_timestep);

		Restart sample;
		string restart_filename = "./" + sample.folder_name + "/" + "restart.kmc_simulation." + step_index.str();
		ifstream inputfile;
		inputfile.open(restart_filename);

		int sol_type, sol_id, curr_location_id;
		for (size_t i = 0; i < input->number_of_solute_per_type.size(); ++i) {
			for (int j = 0; j < input->number_of_solute_per_type[i]; ++j) {
				inputfile >> sol_type >> sol_id >> curr_location_id;
				solutes.push_back(MovingSpecies(sol_type, sol_id, &lattice_sites[curr_location_id]));
				solutes[sol_id].update_next_locations();
			}
		}

		inputfile >> input->start_time;
	}

	double SimulationBox::calculate_energy(LatticeSite* const &site, const int &ref_type) {
		double total_energy = 0.0;

		if (input->include_ngb[1]) {
            for (size_t i = 0; i < site->first_nn.size(); ++i) {
                int ngb_type = (site->first_nn[i]->occupant < 0) ? 0 : solutes[site->first_nn[i]->occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][1];
            }
		}

		if (input->include_ngb[2]) {
            for (size_t i = 0; i < site->second_nn.size(); ++i) {
                int ngb_type = (site->second_nn[i]->occupant < 0) ? 0 : solutes[site->second_nn[i]->occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][2];
            }
		}

		if (input->include_ngb[3]) {
            for (size_t i = 0; i < site->third_nn.size(); ++i) {
                int ngb_type = (site->third_nn[i]->occupant < 0) ? 0 : solutes[site->third_nn[i]->occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][3];
            }
		}

		if (input->include_ngb[4]) {
            for (size_t i = 0; i < site->fourth_nn.size(); ++i) {
                int ngb_type = (site->fourth_nn[i]->occupant < 0) ? 0 : solutes[site->fourth_nn[i]->occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][4];
            }
		}

		if (input->include_ngb[5]) {
            for (size_t i = 0; i < site->fifth_nn.size(); ++i) {
                int ngb_type = (site->fifth_nn[i]->occupant < 0) ? 0 : solutes[site->fifth_nn[i]->occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][5];
            }
		}

		return total_energy;
	}

	double SimulationBox::calculate_energy(LatticeSite* const &target_site, LatticeSite* const &current_site, const int &ref_type, bool &flag) {
		double total_energy = 0.0;

		if (input->include_ngb[1]) {
            for (size_t i = 0; i < target_site->first_nn.size(); ++i) {
                int occupant = (target_site->first_nn[i]->id == current_site->id) ? target_site->occupant : target_site->first_nn[i]->occupant;
                int ngb_type = (occupant < 0) ? 0 : solutes[occupant].type;
                flag = (ref_type == ngb_type);
                total_energy += input->e_species[ref_type][ngb_type][1];
            }
		}

		if (input->include_ngb[2]) {
            for (size_t i = 0; i < target_site->second_nn.size(); ++i) {
                int occupant = (target_site->second_nn[i]->id == current_site->id) ? target_site->occupant : target_site->second_nn[i]->occupant;
                int ngb_type = (occupant < 0) ? 0 : solutes[occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][2];
            }
		}

		if (input->include_ngb[3]) {
            for (size_t i = 0; i < target_site->third_nn.size(); ++i) {
                int occupant = (target_site->third_nn[i]->id == current_site->id) ? target_site->occupant : target_site->third_nn[i]->occupant;
                int ngb_type = (occupant < 0) ? 0 : solutes[occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][3];
            }
		}

		if (input->include_ngb[4]) {
            for (size_t i = 0; i < target_site->fourth_nn.size(); ++i) {
                int occupant = (target_site->fourth_nn[i]->id == current_site->id) ? target_site->occupant : target_site->fourth_nn[i]->occupant;
                int ngb_type = (occupant < 0) ? 0 : solutes[occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][4];
            }
		}

		if (input->include_ngb[5]) {
            for (size_t i = 0; i < target_site->fifth_nn.size(); ++i) {
                int occupant = (target_site->fifth_nn[i]->id == current_site->id) ? target_site->occupant : target_site->fifth_nn[i]->occupant;
                int ngb_type = (occupant < 0) ? 0 : solutes[occupant].type;
                total_energy += input->e_species[ref_type][ngb_type][5];
            }
		}

		return total_energy;
	}

	void SimulationBox::generate_solute_energy_and_rate() {
		for (size_t i = 0; i < solutes.size(); ++i) {
			solutes[i].current_energy = calculate_energy(solutes[i].curr_location, solutes[i].type);
			for (size_t j = 0; j < solutes[i].next_locations.size(); ++j) {
				bool is_identical;
				solutes[i].next_energies.push_back(calculate_energy(solutes[i].next_locations[j], solutes[i].curr_location, solutes[i].type, is_identical));

				//calculate migration energy
				double migration_energy = input->e_migrate[solutes[i].type];
				migration_energy += (solutes[i].next_energies[j] - solutes[i].current_energy) / 2.0;

				//calculate event rate
				double event_rate = input->rate_factor*input->solute_rate[solutes[i].type]*exp(-1.0*migration_energy);
				if (is_identical) event_rate *= flag_for_identical_swap;

				solutes[i].rates.push_back(event_rate);
			}
		}
	}

	void SimulationBox::update_solute_energy_and_rate() {
		for (size_t i = 0; i < solutes.size(); ++i) {
			solutes[i].current_energy = calculate_energy(solutes[i].curr_location, solutes[i].type);
			for (size_t j = 0; j < solutes[i].next_locations.size(); ++j) {
				bool is_identical;
				solutes[i].next_energies[j] = calculate_energy(solutes[i].next_locations[j], solutes[i].curr_location, solutes[i].type, is_identical);

				//calculate migration energy
				double migration_energy = input->e_migrate[solutes[i].type];
				migration_energy += (solutes[i].next_energies[j] - solutes[i].current_energy) / 2.0;

				//calculate event rate
				double event_rate = input->rate_factor*input->solute_rate[solutes[i].type]*exp(-1.0*migration_energy);
				if (is_identical) event_rate *= flag_for_identical_swap;
				solutes[i].rates[j] = event_rate;
			}
		}
	}

	void SimulationBox::generate_solute_from_scratch() {
		generate_solute_location_from_scratch();
		generate_solute_energy_and_rate();
	}

	void SimulationBox::generate_solute_from_restart() {
		generate_solute_location_from_restart();
		generate_solute_energy_and_rate();
	}

	void SimulationBox::print_solute() {
		for (auto &&solute : solutes) {
			solute.print();
			solute.print_future_energy_and_rate();
		}
	}

	void SimulationBox::print_solute_future_energy() {
		for (auto &&solute : solutes) {
			solute.print_future_energy_and_rate();
		}
	}
