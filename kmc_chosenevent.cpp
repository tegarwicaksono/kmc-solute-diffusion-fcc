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

    ChosenEvent::ChosenEvent(const ChosenEvent& other)
        : species(other.species)
        , direction(other.direction)
    { }

    ChosenEvent::ChosenEvent(ChosenEvent&& other)
        : species{std::move(other.species)}
        , direction{std::move(other.direction)}
    {
        other.species = nullptr;
    }

    ChosenEvent& ChosenEvent::operator= (ChosenEvent other) {
        swap(*this, other);
        return *this;
    }

	void ChosenEvent::specify(const Rate &rate) {
		species = rate.species;
		direction = rate.direction;
	}

    void swap(ChosenEvent &a, ChosenEvent &b) {
        using std::swap;
        swap(a.species, b.species);
        swap(a.direction, b.direction);
    }

