/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#ifndef KMC_LATTICESITE_H_INCLUDED
#define KMC_LATTICESITE_H_INCLUDED

#include<vector>

using namespace std;

class LatticeSite {
public:
	int id;
	int occupant = -1;	//by default, it's -1 (matrix). else, species = id of the solute
	vector<double> xyz;
	vector<vector<LatticeSite*> > nth_neighbours;
	// i.e. nth_neighbours[1] = 1st nearest neighbours, nth_neighbours[2] = 2nd nearest neighbours, up to 6

	LatticeSite();
	LatticeSite(const int &id, const vector<double> &xyz);

    void allocate_neighbours(const int &);
	void print();
};


#endif // KMC_LATTICESITE_H_INCLUDED
