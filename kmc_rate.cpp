/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
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
