/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#include "kmc_chosenevent.h"

	ChosenEvent::ChosenEvent() {}
	ChosenEvent::ChosenEvent(MovingSpecies* const &species, const int &direction) :
		species(species), direction(direction) {}

	void ChosenEvent::specify(const Rate &rate) {
		species = rate.species;
		direction = rate.direction;
	}
