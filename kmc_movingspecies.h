/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#ifndef KMC_MOVINGSPECIES_H_INCLUDED
#define KMC_MOVINGSPECIES_H_INCLUDED

#include "kmc_latticesite.h"
#include <utility>


using namespace std;

class MovingSpecies {
public:
	int type;	//1 = solute of type 1, 2 = solute of type 2
	int id;		//id ranges from 0 to total number of solutes
	LatticeSite* prev_location;
	LatticeSite* curr_location;
	double    current_energy;

	vector<LatticeSite*> next_locations;
	vector<double> rates;
	vector<double> next_energies;

	MovingSpecies();
	MovingSpecies(const int &type, const int &id);
	MovingSpecies(const int &type, const int &id, LatticeSite* const &curr_location);
    MovingSpecies(const MovingSpecies& other);
    MovingSpecies(MovingSpecies&& other);
    MovingSpecies& operator= (MovingSpecies other);

	virtual ~MovingSpecies() = default;
    friend void swap(MovingSpecies &a, MovingSpecies &b);

	void print();
	void print_future_energy_and_rate();
	void update_location(const int &direction);
	void update_location(LatticeSite* const &assigned_site);
	void update_next_locations();
};

#endif // KMC_MOVINGSPECIES_H_INCLUDED
