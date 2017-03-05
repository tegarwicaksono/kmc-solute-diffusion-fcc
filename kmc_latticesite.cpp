/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#include "kmc_latticesite.h"
#include <iostream>

using namespace std;

	LatticeSite::LatticeSite() {}

	LatticeSite::LatticeSite(const int &id, const vector<double> &xyz) : id(id), xyz(xyz){}
	void LatticeSite::print() {
		cout << "site id = " << id;
		cout << ", coor = " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "\n";
	}

	void LatticeSite::allocate_neighbours(const int& nth_order) {
        nth_neighbours.assign(nth_order, vector<LatticeSite*>());
	}
