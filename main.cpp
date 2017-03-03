/*  KMC Simulation for FCC lattice with diffusion by swapping
    Author: Tegar Wicaksono (tegar@alumni.ubc.ca)
    Written: March 2017
*/

#include <iostream>
#include <string>
#include "kmc_latticesite.h"
#include "kmc_movingspecies.h"
#include "kmc_rate.h"
#include "kmc_chosenevent.h"
#include "kmc_inputdata.h"
#include "kmc_restart.h"
#include "kmc_snapshot.h"
#include "kmc_simulationbox.h"
#include "kmc_kmcsimulation.h"

using namespace std;

int main()
{
    cout << "Starting KMC simulation for solute atoms in FCC\n";
    cout << "\nPlease enter the name of your input file, e.g. test_input.txt (include .txt part too)\n";

    string input_filename;
    cin >> input_filename;

	InputData  id;
	cout << "\nreading input file ...\n";
	id.initialize(input_filename);
	cout << "finish reading input file\n\n";
	cout << "generating simulation box ...\n";
	cout << "(this may take a while... be patient...)\n";

	SimulationBox sb;
	sb.generate_box(&id);
	cout << "finish generating simulation box\n";

	KMCSimulation kmc;
	kmc.run_kmc_simulation(&sb);

    cout << endl << "simulation completed" << endl;
    return 0;
}
