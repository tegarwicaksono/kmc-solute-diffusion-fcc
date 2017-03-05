/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
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

		//print_site_neighbours();
		//print_nn_fcc(2);
		if (input->start_from_restart) {
			generate_solute_from_restart();
		} else {
			generate_solute_from_scratch();
		}
	}

	double SimulationBox::calculate_distance(LatticeSite* const &site1, LatticeSite* const &site2) {
        double delta, distance = 0.0;

        for (int i = 0; i < 3; ++i) {
            delta = fabs(site1->xyz[i] - site2->xyz[i]);
            if (delta > input->box_length[i] / 2.0) {
                delta = fabs(delta - input->box_length[i]);
            }
            distance += delta*delta;
        }
        return sqrt(distance);
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

	vector<vector<double> > SimulationBox::zeroth_nearest_neighbours_fcc() {
        vector<vector<double> > list;
        list.push_back(std::initializer_list<double>{0.0,0.0,0.0});
        return list;
	}

	vector<vector<double> > SimulationBox::first_nearest_neighbours_fcc() {
		vector<vector<double> > list;
		list.push_back(std::initializer_list<double>{0.0,-0.5,-0.5});	//1
		list.push_back(std::initializer_list<double>{0.0,-0.5,+0.5});	//2
		list.push_back(std::initializer_list<double>{0.0,0.5,-0.5});	//3
		list.push_back(std::initializer_list<double>{0.0,0.5,0.5});		//4
		list.push_back(std::initializer_list<double>{-0.5,0.0,-0.5});	//5
		list.push_back(std::initializer_list<double>{-0.5,0.0,0.5});	//6
		list.push_back(std::initializer_list<double>{0.5,0.0,-0.5});	//7
		list.push_back(std::initializer_list<double>{0.5,0.0,0.5});		//8
		list.push_back(std::initializer_list<double>{-0.5,-0.5,0.0});	//9
		list.push_back(std::initializer_list<double>{-0.5,0.5,0.0});	//10
		list.push_back(std::initializer_list<double>{0.5,-0.5,0.0});	//11
		list.push_back(std::initializer_list<double>{0.5,0.5,0.0});		//12
		return list;
	}

	vector<vector<double> > SimulationBox::second_nearest_neighbours_fcc() {
		vector<vector<double> > list;
		list.push_back(std::initializer_list<double>{-1.0,0.0,0.0});	//1
		list.push_back(std::initializer_list<double>{1.0,0.0,0.0});	    //2
		list.push_back(std::initializer_list<double>{0.0,-1.0,0.0});	//3
		list.push_back(std::initializer_list<double>{0.0,1.0,0.0});	    //4
		list.push_back(std::initializer_list<double>{0.0,0.0,-1.0});	//5
		list.push_back(std::initializer_list<double>{0.0,0.0,1.0});	    //6
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
        vector<vector<double> > list;
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

		return list;
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
						lattice_sites.push_back(LatticeSite(latt_id, xyz));
						lattice_sites[latt_id].allocate_neighbours(input->max_ngb_distance + 1);
						++latt_id;
					}
				}
			}
		}
	}

	vector<vector<double>> SimulationBox::find_nn_per_site(const LatticeSite &site,
		const vector<vector<double>> &neighbours) {

		vector<vector<double> > list;
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

    void SimulationBox::assign_nth_ngbs_fcc_per_site(const size_t &index, LatticeSite &site, const vector<vector<double>> &nth_ngbs) {
        vector<vector<double> > nns = find_nn_per_site(site, nth_ngbs);
        for (const vector<double> &nn : nns) {
            site.nth_neighbours[index].push_back(find_latt_id_fcc(nn));
        }
	}

    void SimulationBox::generate_nth_nearest_neighbours() {
        vector<vector<vector<double> > > neighbours;

        neighbours.push_back(zeroth_nearest_neighbours_fcc());
        neighbours.push_back(first_nearest_neighbours_fcc());
        neighbours.push_back(second_nearest_neighbours_fcc());
        neighbours.push_back(third_nearest_neighbours_fcc());
        neighbours.push_back(fourth_nearest_neighbours_fcc());
        neighbours.push_back(fifth_nearest_neighbours_fcc());

/*
        cout << "for neighbours vector, there are " << endl;
        cout << neighbours[0].size() << " zeroth NNs" << endl;
        cout << neighbours[1].size() << " first NNs" << endl;
        cout << neighbours[2].size() << " second NNs" << endl;
        cout << neighbours[3].size() << " third NNs" << endl;
        cout << neighbours[4].size() << " fourth NNs" << endl;
        cout << neighbours[5].size() << " fifth NNs" << endl;
*/
        for (size_t i = 0; i < lattice_sites.size(); ++i) {
            for (size_t j = 0; j < neighbours.size(); ++j) {
                assign_nth_ngbs_fcc_per_site(j, lattice_sites[i], neighbours[j]);
            }
        }
    }

	void SimulationBox::generate_neighbours() {
		generate_nth_nearest_neighbours();
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
		cout << "In object sample, the first_nn has " << lattice_sites[sample_id].nth_neighbours[1].size() << " elements\n";
		for (LatticeSite* const &nn1 : lattice_sites[sample_id].nth_neighbours[1]) {
			nn1->print();
		}
	}

    void SimulationBox::print_nn_fcc(const int &index) {
		for (size_t i = 0; i < lattice_sites.size(); ++i) {
			cout << "site id = " << lattice_sites[i].id << ", occup = " << lattice_sites[i].occupant << ", " <<  index << "-th ngbs = ";
			for (size_t j = 0; j < lattice_sites[i].nth_neighbours[index].size(); ++j) {
				cout << lattice_sites[i].nth_neighbours[index][j]->id << ", ";
			}
			cout << "  distance = " << endl << "\t";
			for (size_t j = 0; j < lattice_sites[i].nth_neighbours[index].size(); ++j) {
				cout << calculate_distance(&lattice_sites[i], lattice_sites[i].nth_neighbours[index][j]) << ", ";
			}
			cout << endl;
		}
	}

	void SimulationBox::print_test_2nn_fcc() {
		int sample_id;

		sample_id = 3;
		cout << "=========== SAMPLE TEST 2NN ===========\n";
		cout << "Reference site: ";
		lattice_sites[sample_id].print();
		cout << "In object sample, the second_nn has " << lattice_sites[sample_id].nth_neighbours[2].size() << " elements\n";
		for (LatticeSite* const &nn2 : lattice_sites[sample_id].nth_neighbours[2]) {
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
		for (size_t i = 1; i < input->include_ngb.size(); ++i) {
		    double temp = 0.0;
            if (input->include_ngb[i]) {
                for (size_t j = 0; j < site->nth_neighbours[i].size(); ++j) {
                    int ngb_type = (site->nth_neighbours[i][j]->occupant < 0) ? 0 : solutes[site->nth_neighbours[i][j]->occupant].type;
                    temp += input->e_species[ref_type][ngb_type][i];
                }
            }
            total_energy += temp;
		}

		return 0.5*total_energy;
	}

	double SimulationBox::calculate_energy(LatticeSite* const &target_site, LatticeSite* const &current_site, const int &ref_type) {
		double total_energy = 0.0;

		for (size_t i = 1; i < input->include_ngb.size(); ++i) {
            if (input->include_ngb[i]) {
                for (size_t j = 0; j < target_site->nth_neighbours[i].size(); ++j) {
                    int occupant = (target_site->nth_neighbours[i][j]->id == current_site->id) ? target_site->occupant : target_site->nth_neighbours[i][j]->occupant;
                    int ngb_type = (occupant < 0) ? 0 : solutes[occupant].type;

                    total_energy += input->e_species[ref_type][ngb_type][i];
                }
            }
		}

		return 0.5*total_energy;
	}

	void SimulationBox::generate_solute_energy_and_rate() {
		for (size_t i = 0; i < solutes.size(); ++i) {
			solutes[i].current_energy = calculate_energy(solutes[i].curr_location, solutes[i].type);
			//cout << "solute #" << i << " current energy = " << solutes[i].current_energy << endl;
			for (size_t j = 0; j < solutes[i].next_locations.size(); ++j) {
				//bool is_identical;
				//solutes[i].next_energies.push_back(calculate_energy(solutes[i].next_locations[j], solutes[i].curr_location, solutes[i].type, is_identical));
				solutes[i].next_energies.push_back(calculate_energy(solutes[i].next_locations[j], solutes[i].curr_location, solutes[i].type));

				//calculate migration energy
				double migration_energy = input->e_migrate[solutes[i].type];
				migration_energy += (solutes[i].next_energies[j] - solutes[i].current_energy) / 2.0;

				//calculate event rate
				double event_rate = input->rate_factor*input->solute_rate[solutes[i].type]*exp(-1.0*migration_energy);
				//if (is_identical) event_rate *= flag_for_identical_swap;

				solutes[i].rates.push_back(event_rate);
			}
		}
	}

	void SimulationBox::update_solute_energy_and_rate() {
		for (size_t i = 0; i < solutes.size(); ++i) {
			solutes[i].current_energy = calculate_energy(solutes[i].curr_location, solutes[i].type);
			for (size_t j = 0; j < solutes[i].next_locations.size(); ++j) {
				//bool is_identical;
				//solutes[i].next_energies[j] = calculate_energy(solutes[i].next_locations[j], solutes[i].curr_location, solutes[i].type, is_identical);
				solutes[i].next_energies[j] = calculate_energy(solutes[i].next_locations[j], solutes[i].curr_location, solutes[i].type);

				//calculate migration energy
				double migration_energy = input->e_migrate[solutes[i].type];
				migration_energy += (solutes[i].next_energies[j] - solutes[i].current_energy) / 2.0;

				//calculate event rate
				double event_rate = input->rate_factor*input->solute_rate[solutes[i].type]*exp(-1.0*migration_energy);
//				if (is_identical) event_rate *= flag_for_identical_swap;
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

	void SimulationBox::print_site_neighbours() {
        for (size_t i = 0; i < lattice_sites.size(); ++i) {
            cout << "For lattice site = " << lattice_sites[i].id << ", there are : ";
            for (size_t j = 0; j < lattice_sites[i].nth_neighbours.size(); ++j) {
                cout << lattice_sites[i].nth_neighbours[j].size() << " ngbs, ";
            }
            cout << endl;
        }
	}
