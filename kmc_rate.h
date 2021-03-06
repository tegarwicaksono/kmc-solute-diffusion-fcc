/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#ifndef KMC_RATE_H_INCLUDED
#define KMC_RATE_H_INCLUDED
#include "kmc_movingspecies.h"
#include <utility>


class Rate {
public:
	MovingSpecies* species;
	int 	direction;
	double 	cumulative;

	Rate();
	Rate(MovingSpecies* const &species, const int &direction);
	Rate(const Rate& other);
	Rate(Rate&& other);

    Rate& operator= (Rate other);

	virtual ~Rate() = default;
    friend void swap(Rate &a, Rate &b);

	void print();
};

#endif // KMC_RATE_H_INCLUDED
