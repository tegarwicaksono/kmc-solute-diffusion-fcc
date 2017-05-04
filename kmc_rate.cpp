/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#include "kmc_rate.h"
#include <iostream>

using namespace std;

	Rate::Rate() : cumulative{0.0}
	{}

	Rate::Rate(MovingSpecies* const &species, const int &direction)
		: species(species)
		, direction(direction)
		, cumulative(0.0)
    {}

    Rate::Rate(const Rate& other)
        : species{other.species}
        , direction{other.direction}
        , cumulative{other.cumulative}
    { }

    Rate::Rate(Rate&& other)
        : species{std::move(other.species)}
        , direction{std::move(other.direction)}
        , cumulative{std::move(other.cumulative)}
    {
        other.species = nullptr;
    }

    Rate& Rate::operator= (Rate other) {
        swap(*this, other);
        return *this;
    }

	void swap(Rate& a, Rate& b) {
        using std::swap;
        swap(a.species, b.species);
        swap(a.direction, b.direction);
        swap(a.cumulative, b.cumulative);
	}

	void Rate::print() {
		cout << "species = " << species->id << " is type " << species->type << ", direction = " << direction << ", cumul rate = " << cumulative << endl;
	}
