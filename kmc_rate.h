/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#ifndef KMC_RATE_H_INCLUDED
#define KMC_RATE_H_INCLUDED
#include "kmc_movingspecies.h"

class Rate {
public:
	MovingSpecies* species;
	int 	direction;
	double 	cumulative;

	Rate();
	Rate(MovingSpecies* const &species, const int &direction);

	void print();
};

#endif // KMC_RATE_H_INCLUDED
