/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#include "kmc_inputdata.h"
#include <iostream>
#include <cmath>
#include <regex>
#include <string>
#include <fstream>
#include <iomanip>
#include <limits>

using namespace std;


	InputData::InputData() : box_length(3), unit_length(3) {}

	void InputData::assign_nearest_neighbour_distance() {
        nn_distance.assign(4,0.0);
        nn_distance[0] = 0.5*sqrt(0.0);
        nn_distance[1] = 0.5*sqrt(2.0);
        nn_distance[2] = 0.5*sqrt(4.0);
        nn_distance[3] = 0.5*sqrt(6.0);

        include_ngb.assign(4,false);
	}

	void InputData::assign_parameter_name() {
		parameter_name.push_back("box_length_x");	//in unit cell //0
		parameter_name.push_back("box_length_y");	//in unit cell //1
		parameter_name.push_back("box_length_z");	//in unit cell //2
		parameter_name.push_back("unit_length_x");	//in angstrom  //3
		parameter_name.push_back("unit_length_y");	//in angstrom  //4
		parameter_name.push_back("unit_length_z");	//in angstrom  //5
		parameter_name.push_back("absolute_temp");	//in kelvin    //6
		parameter_name.push_back("number_of_solute_type");		   //7
		parameter_name.push_back("number_of_solute_of_type"); //8	//format: number_of_solute_of_type 1 10
		parameter_name.push_back("pre_exp_rate_term"); //9 //in second^-1
		parameter_name.push_back("solute_migration_energy"); //10 //in eV
		parameter_name.push_back("solute_solute_interaction_energy"); //11 //in eV, format: solute_solute_interaction_energy 1 2 0.5
		parameter_name.push_back("matrix_matrix_interaction_energy"); //12 //in eV
		parameter_name.push_back("solute_matrix_interaction_energy"); //13 //in eV
		parameter_name.push_back("number_of_KMC_steps"); //14
		parameter_name.push_back("dump_snapshot");	//15
		parameter_name.push_back("dump_restart");	//16
		parameter_name.push_back("read_from_restart"); //17
		parameter_name.push_back("effective_interaction_distance"); //18
	}

	vector<string> InputData::split_token(const string &input) {
		std::regex rgx("[\\s]+");
		std::sregex_token_iterator it(input.begin(), input.end(), rgx, -1);
		std::sregex_token_iterator reg_end;
		vector<string> tokens;
		for (; it != reg_end; ++it) {
			if (it->str().size() > 0) {
				tokens.push_back(it->str());
			}
		}
		return tokens;
	}

	void InputData::assign_parameter(const int &it, const vector<string> &tokens, unordered_set<int> &spec_param) {
		if (it < 3) {
			box_length[it] = std::stoi(tokens[1]);
			spec_param.insert(it);
		} else if (it >= 3 && it < 6) {
			unit_length[it - 3] = std::stod(tokens[1]);
			spec_param.insert(it);
		} else if (it == 6) {
			abs_temperature = std::stod(tokens[1]);
			spec_param.insert(it);
		} else if (it == 7) {
			number_of_solute_type = std::stoi(tokens[1]);
			spec_param.insert(it);

			//assigning space for the vectors of solute/matrix properties
			number_of_solute_per_type.assign(number_of_solute_type, int());

			species_interaction_energy.assign(number_of_solute_type + 1, vector<vector<double> >(number_of_solute_type + 1, vector<double>(4, 0.0)));
			e_species.assign(number_of_solute_type + 1, vector<vector<double> >(number_of_solute_type + 1, vector<double>(4, 0.0)));
			rate_pre_exponential.assign(number_of_solute_type + 1, double());
			solute_rate.assign(number_of_solute_type + 1, double());
			solute_migration_energy.assign(number_of_solute_type + 1, double());
			e_migrate.assign(number_of_solute_type + 1, double());
			include_species_in_snapshot.assign(number_of_solute_type + 1, false);

			//assign element at index zero
			rate_pre_exponential[0] = 0.0;
			solute_rate[0] = 0.0;
			solute_migration_energy[0] = 0.0;
			e_migrate[0] = 0.0;


		} else if (it == 8) {
			int sol = std::stoi(tokens[1]) - 1;
			int n_sol = std::stoi(tokens[2]);
			number_of_solute_per_type[sol] = n_sol;
			int id = std::stoi(std::to_string(it) + std::to_string(sol));
			spec_param.insert(id);
			//cout << "number_of_solute_per_type[" << sol << "] = " << n_sol << "\n";
			//cout << "inserted_id = " << id << "\n";
		} else if (it == 9) {
			int sol = std::stoi(tokens[1]);
			double rate_sol = std::stod(tokens[2]);
			rate_pre_exponential[sol] = rate_sol;
			int id = std::stoi(std::to_string(it) + std::to_string(sol - 1));
			spec_param.insert(id);
			//cout << "rate_pre_exponential[" << sol << "] = " << rate_sol << "\n";
			//cout << "inserted_id = " << id << "\n";
		} else if (it == 10) {
			int sol = std::stoi(tokens[1]);
			double diff_sol = std::stod(tokens[2]);
			solute_migration_energy[sol] = diff_sol;
			int id = std::stoi(std::to_string(it) + std::to_string(sol - 1));
			spec_param.insert(id);
			//cout << "solute_migration_energy[" << sol << "] = " << diff_sol << "\n";
			//cout << "inserted_id = " << id << "\n";
		} else if (it == 11) {
			int sol1 = std::stoi(tokens[1]);
			int sol2 = std::stoi(tokens[2]);
			double energy = std::stod(tokens[3]);

			int nn_distance = 0;
			if (tokens[4].compare("first") == 0) {
                nn_distance = 1;
			} else if (tokens[4].compare("second") == 0) {
                nn_distance = 2;
			} else if (tokens[4].compare("third") == 0) {
                nn_distance = 3;
			}
			species_interaction_energy[sol1][sol2][nn_distance] = energy;
			species_interaction_energy[sol2][sol1][nn_distance] = energy;

			int id1 = std::stoi(std::to_string(it) + std::to_string(sol1) + std::to_string(sol2));
			int id2 = std::stoi(std::to_string(it) + std::to_string(sol2) + std::to_string(sol1));
			spec_param.insert(id1);
			spec_param.insert(id2);
			//cout << "species_interaction_energy[" << sol1 << "][" << sol2 << " = " << energy << "\n";
			//cout << "inserted_id = " << id1 << " and " << id2 << "\n";
		} else if (it == 12) {
			double energy = std::stod(tokens[1]);

			int nn_distance = 0;
			if (tokens[2].compare("first") == 0) {
                nn_distance = 1;
			} else if (tokens[2].compare("second") == 0) {
                nn_distance = 2;
			} else if (tokens[2].compare("third") == 0) {
                nn_distance = 3;
			}

			species_interaction_energy[0][0][nn_distance] = energy;
			spec_param.insert(1100);
			//cout << "species_interaction_energy[0][0] = " << energy << "\n";
			//cout << "inserted_id = 1100\n";
		} else if (it == 13) {
			int sol = std::stoi(tokens[1]);
			double energy = std::stod(tokens[2]);

			int nn_distance = 0;
			if (tokens[3].compare("first") == 0) {
                nn_distance = 1;
			} else if (tokens[3].compare("second") == 0) {
                nn_distance = 2;
			} else if (tokens[3].compare("third") == 0) {
                nn_distance = 3;
			}

			species_interaction_energy[0][sol][nn_distance] = energy;
			species_interaction_energy[sol][0][nn_distance] = energy;

			int id1 = std::stoi(std::to_string(110) + std::to_string(sol));
			int id2 = std::stoi(std::to_string(11) + std::to_string(sol) + "0");
			spec_param.insert(id1);
			spec_param.insert(id2);
			//cout << "species_interaction_energy[0][" << sol << "] = " << energy << "\n";
			//cout << "inserted_id = " << id1 << " and " << id2 << "\n";
		} else if (it == 14) {
			total_KMC_steps = std::stoull(tokens[1]);
			spec_param.insert(it);
		} else if (it == 15) {
			dump_snapshot = true;
			period_snapshot = std::stoull(tokens[1]);

			for (size_t i = 2; i < tokens.size(); ++i) {
				int species_to_dump = std::stoi(tokens[i]);
				if (species_to_dump <= number_of_solute_type) {
					include_species_in_snapshot[species_to_dump] = true;
				}
			}
			spec_param.insert(it);

		} else if (it == 16) {
			dump_restart = true;
			period_restart = std::stoull(tokens[1]);
			spec_param.insert(it);
		} else if (it == 17) {
			start_from_restart = true;
			restart_timestep = std::stoull(tokens[1]);
			spec_param.insert(it);
		} else if (it == 18) {
            if (tokens[1].compare("first") == 0) {
                eff_ngb_distance = 1;
            } else if (tokens[1].compare("second") == 0) {
                eff_ngb_distance = 2;
            } else if (tokens[1].compare("third") == 0) {
                eff_ngb_distance = 3;
            } else {
                eff_ngb_distance = 1;
            }

            for (int i = 0; i <= eff_ngb_distance; ++i) {
                include_ngb[i] = true;
            }
		}
	}

	bool InputData::check_if_parameters_assigned(const unordered_set<int> &spec_param) {
		int count = 0;
		while (count <= 7) {
			if (spec_param.find(count) == spec_param.end()) {
				cout << "Parameter " << parameter_name[count] << " has not been specified\n";
				return false;
			} else {
				++count;
			}
		}

		while (count <= 10) {
			int linemax = count*10 + (number_of_solute_type - 1);
			int lineini = count*10;
			for (int i = lineini; i <= linemax; ++i) {
				if (spec_param.find(i) == spec_param.end()) {
					cout << "Parameter " << parameter_name[count] << " for solute type " << i%lineini + 1 << " has not been specified\n";
					return false;
				}
			}
			++count;
		}

		int lineini = count*100;

		if (spec_param.find(lineini) == spec_param.end()) {
			cout << "Parameter " << parameter_name[12] << " has not been defined\n";
			return false;
		}

		for (int i = 0; i < number_of_solute_type; ++i) {
			if (spec_param.find(lineini + i + 1) == spec_param.end()) {
				cout << "Parameter " << parameter_name[13] << " for solute type " << i + 1 << " has not been defined\n";
				return false;
			}
		}

		for (int i = 0; i < number_of_solute_type; ++i) {
			for (int j = i; j < number_of_solute_type; ++j) {
				int id = std::stoi(std::to_string(i + 1) + std::to_string(j + 1));
				if (spec_param.find(lineini + id) == spec_param.end()) {
					cout << "Parameter " << parameter_name[11] << " for solute types " << i + 1 << " and " << j + 1 << " has not been defined\n";
					return false;
				}
			}
		}

		count = 14;
		if (spec_param.find(count) == spec_param.end()) {
			cout << "Parameter " << parameter_name[14] << " has not been defined\n";
		}

		return true;
	}

	bool InputData::process_input_line(const string &input, unordered_set<int> &spec_param) {
		if (input.at(0) == '#' || input.size() == 1) {
			//cout << "comment line ignored\n";
			return true;
		}

		vector<string> tokens = split_token(input);

		bool found_match = false;
		int  it = 0;

		while(!found_match) {
			if (tokens[0].compare(parameter_name[it]) == 0) {
				found_match = true;
			} else {
				++it;
			}
		}

		if (!found_match) {
			return false;
		}

		if (it > 7 && spec_param.find(7) == spec_param.end()) {
			cout << "You need to define number of solute type before defining any of solute properties\n";
			return false;
		}

		assign_parameter(it, tokens, spec_param);
		return true;
	}

	bool InputData::read_from_file(const string &filename) {
		unordered_set<int> specified_param;
		ifstream inputfile;
		inputfile.open(filename);

		string line;
		for (int i = 0; i < 3; ++i) {
			getline(inputfile, line);
		}

		int count = 4;
		while (getline(inputfile, line)) {
			if (!process_input_line(line, specified_param)) {
				break;
			} else {
				++count;
			}
		}
		return check_if_parameters_assigned(specified_param);
	}

	void InputData::convert_input_data() {
		double kB = 1.38064852e-23;	//Boltzmann constant
		double eV = 1.60217662e-19;	//elementary charge

		e_species = species_interaction_energy;
		for (size_t i = 0; i < e_species.size(); ++i) {
			for (size_t j = 0; j < e_species[i].size(); ++j) {
                for (size_t k = 0; k < e_species[i][j].size(); ++k) {
                    e_species[i][j][k] *= eV / (kB * abs_temperature);
                }
			}
		}

		e_migrate = solute_migration_energy;
		for (size_t i = 0; i < e_migrate.size(); ++i) {
			e_migrate[i] *= eV / (kB * abs_temperature);
		}

		solute_rate = rate_pre_exponential;

		rate_factor = std::numeric_limits<double>::max();

		for (size_t i = 1; i < rate_pre_exponential.size(); ++i) {
            if (rate_pre_exponential[i] > 0.0 && rate_factor > rate_pre_exponential[i]) {
                rate_factor = rate_pre_exponential[i];
            }
		}

		for (size_t i = 0; i < solute_rate.size(); ++i) {
			solute_rate[i] /= rate_factor;
		}

		if (start_from_restart) {
			initial_timestep = restart_timestep;
			final_timestep = total_KMC_steps + initial_timestep;
		} else {
			initial_timestep = 0ull;
			final_timestep = total_KMC_steps;
		}
	}

	void InputData::initialize(const string &filename) {
		assign_parameter_name();
        assign_nearest_neighbour_distance();

		if (!read_from_file(filename)) {
			return ;
		}
		convert_input_data();
		//print_scaled_properties();
	}

	void InputData::print_scaled_properties() {
		cout << "Printing scaled rate and energy terms\n\n";
		cout << ">> Rate terms\n";
		cout << ">> Rate factor = " << rate_factor << "\n";
		for (size_t i = 1; i < solute_rate.size(); ++i) {
			cout << ">> Scaled rate for solute " << i << " = " << solute_rate[i] << "\n";
		}

		cout << "\n>> Diffusion energy terms\n";
		for (size_t i = 1; i < e_migrate.size(); ++i) {
			cout << ">> Scaled diffusion energy for solute " << i << " = " << e_migrate[i] << "\n";
		}

		cout << "\n>> Solute interaction with matrix terms\n";
		for (size_t i = 1; i < e_species[0].size(); ++i) {
			cout << ">> Scaled interaction between matrix and solute " << i << " : " << endl;
            cout << "     1st nearest ngb distance = " << e_species[0][i][1] << "\n";
            cout << "     2nd nearest ngb distance = " << e_species[0][i][2] << "\n";
            cout << "     3rd nearest ngb distance = " << e_species[0][i][3] << "\n";

		}

		cout << "\n>> Solute-solute interaction terms\n";
		for (int i = 0; i < number_of_solute_type; ++i) {
			for (int j = i; j < number_of_solute_type; ++j) {
				cout << ">> Scaled interaction between solutes " << i + 1 << " and " << j + 1 << " : " << endl;
				cout << "   1st nearest ngb distance = " << e_species[i+1][j+1][1] << "\n";
				cout << "   2nd nearest ngb distance = " << e_species[i+1][j+1][2] << "\n";
				cout << "   3rd nearest ngb distance = " << e_species[i+1][j+1][3] << "\n";
			}
		}
	}
