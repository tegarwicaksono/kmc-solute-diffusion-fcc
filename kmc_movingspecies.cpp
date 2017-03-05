/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#include "kmc_movingspecies.h"
#include <iostream>

using namespace std;

	MovingSpecies::MovingSpecies() {}
	MovingSpecies::MovingSpecies(const int &type, const int &id) : type(type), id(id) {}
	MovingSpecies::MovingSpecies(const int &type, const int &id, LatticeSite* const &curr_location)
		: type(type), id(id), curr_location(curr_location) {}

	void MovingSpecies::print() {
		cout << "sol. type = " << type << ", id = " << id << ", site = " << curr_location->id;
		cout << ", current_energy: " << current_energy;
		cout << "\n";
	}

	void MovingSpecies::print_future_energy_and_rate() {
		for (size_t i = 0; i < next_energies.size(); ++i) {
			cout << " future site = " << next_locations[i]->id;
			cout << ", future energy = " << next_energies[i] ;
			cout << ", future rate = " << rates[i];
			cout << "\n";
		}
	}

	void MovingSpecies::update_location(const int &direction) {
		prev_location = curr_location;
		curr_location = next_locations[direction];
	}

	void MovingSpecies::update_location(LatticeSite* const &assigned_site) {
		prev_location = curr_location;
		curr_location = assigned_site;
	}

	void MovingSpecies::update_next_locations() {
		next_locations = curr_location->nth_neighbours[1];
	}
