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

	MovingSpecies::MovingSpecies(const MovingSpecies& other)
        : type{other.type}
        , id{other.id}
        , prev_location{other.prev_location}
        , curr_location{other.curr_location}
        , current_energy{other.current_energy}
        , next_locations{other.next_locations}
        , rates{other.rates}
        , next_energies{other.next_energies}
	{ }

	MovingSpecies::MovingSpecies(MovingSpecies&& other)
        : type{std::move(other.type)}
        , id{std::move(other.id)}
        , prev_location{std::move(other.prev_location)}
        , curr_location{std::move(other.curr_location)}
        , current_energy{std::move(other.current_energy)}
        , next_locations{std::move(other.next_locations)}
        , rates{std::move(other.rates)}
        , next_energies{std::move(other.next_energies)}
    {
        other.prev_location = nullptr;
        other.curr_location = nullptr;
    }

    MovingSpecies& MovingSpecies::operator= (MovingSpecies other) {
        swap(*this, other);
        return *this;
    }

	void swap(MovingSpecies& a, MovingSpecies& b) {
        using std::swap;
        swap(a.type, b.type);	//1 = solute of type 1, 2 = solute of type 2
        swap(a.id, b.id);		//id ranges from 0 to total number of solutes
        swap(a.prev_location, b.prev_location);
        swap(a.curr_location, b.curr_location);
        swap(a.current_energy, b.current_energy);

        swap(a.next_locations, b.next_locations);
        swap(a.rates, b.rates);
        swap(a.next_energies, b.next_energies);
	}


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
