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

	Rate::Rate() {}
	Rate::Rate(MovingSpecies* const &species, const int &direction) :
		species(species), direction(direction),
		cumulative(0.0) {}

	void Rate::print() {
		cout << "species = " << species->id << " is type " << species->type << ", direction = " << direction << ", cumul rate = " << cumulative << endl;
	}
