/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#ifndef KMC_CHOSENEVENT_H_INCLUDED
#define KMC_CHOSENEVENT_H_INCLUDED

#include "kmc_movingspecies.h"
#include "kmc_rate.h"

class ChosenEvent {
public:
	MovingSpecies* species;
	int direction;

	ChosenEvent();
	ChosenEvent(MovingSpecies* const &species, const int &direction);

	void specify(const Rate &rate);
};

#endif // KMC_CHOSENEVENT_H_INCLUDED
