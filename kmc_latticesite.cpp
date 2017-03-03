/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
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
