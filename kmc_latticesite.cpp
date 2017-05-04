/*  KMC Simulation for FCC lattice with diffusion
    by species swap and/or vacancy exchange
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017

    Check repository below for the most updated version:
    https://github.com/tegarwicaksono/kmc-solute-diffusion-fcc
*/

#include "kmc_latticesite.h"
#include <iostream>

using namespace std;

	LatticeSite::LatticeSite() : occupant(-1)
	{}

	LatticeSite::LatticeSite(const LatticeSite& other)
        : id{other.id}
        , occupant{other.occupant}	//by default, it's -1 (matrix). else, species = id of the solute
        , xyz{other.xyz}
        , nth_neighbours{other.nth_neighbours}
    { }

	LatticeSite::LatticeSite(LatticeSite&& other)
        : id{std::move(other.id)}
        , occupant{std::move(other.occupant)}	//by default, it's -1 (matrix). else, species = id of the solute
        , xyz{std::move(other.xyz)}
        , nth_neighbours{std::move(other.nth_neighbours)}
	{ }

    LatticeSite& LatticeSite::operator= (LatticeSite other) {
        swap(*this, other);
        return *this;
    }

	LatticeSite::LatticeSite(const int &id, const vector<double> &xyz)
        : id(id)
        , occupant(-1)
        , xyz(xyz)
    {}

	void LatticeSite::print() {
		cout << "site id = " << id;
		cout << ", coor = " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "\n";
	}

	void LatticeSite::allocate_neighbours(const int& nth_order) {
        nth_neighbours.assign(nth_order, vector<LatticeSite*>());
	}

	void swap(LatticeSite& a, LatticeSite& b) {
        using std::swap;
        swap(a.id, b.id);
        swap(a.occupant, b.occupant);
        swap(a.xyz, b.xyz);
        swap(a.nth_neighbours, b.nth_neighbours);
	}
