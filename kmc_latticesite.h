/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
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
	vector<LatticeSite*> first_nn;
	vector<LatticeSite*> second_nn;
	vector<LatticeSite*> third_nn;
	vector<LatticeSite*> fourth_nn;
	vector<LatticeSite*> fifth_nn;

	LatticeSite();
	LatticeSite(const int &id, const vector<double> &xyz);

	void print();
};


#endif // KMC_LATTICESITE_H_INCLUDED
