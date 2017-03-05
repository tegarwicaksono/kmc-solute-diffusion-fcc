/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#include "kmc_chosenevent.h"

	ChosenEvent::ChosenEvent() {}
	ChosenEvent::ChosenEvent(MovingSpecies* const &species, const int &direction) :
		species(species), direction(direction) {}

	void ChosenEvent::specify(const Rate &rate) {
		species = rate.species;
		direction = rate.direction;
	}
