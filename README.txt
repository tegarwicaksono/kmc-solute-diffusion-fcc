/*****************  KMC simulation for solute diffusion in FCC crystal *******************/

Purpose: 
The codes are designed to simulate the diffusion of atoms in an FCC crystal lattice using a kinetic Monte Carlo method. The atoms are assumed to be in a fix set of lattice sites. User needs to provide an input file to run the simulation. A sample of an input file is stored in /input/test_input.txt

Compilation:
The codes are written in C++ language and requires C++11 compilers to build. User needs to create .o files for each kmc_**.cpp files to which the main.cpp file can be linked to create an executable binary.

Pre-requisite for simulations:
To make the simulation running headache-free, user is advised to do the following before running the simulation:
1. Make sure an input file is created that contains all the necessary parameters. Please see /input/test_input.txt to view an example of an input file.
2. Create a dedicated folder where user can store the executable binary.
3. In that folder, create these 3 subfolders (the name has to be matched): (1) dump_snapshot, (2) dump_restart, (3) log
4. If interested in running simulations with the same input file, but different spatial configuration of solute atoms, create another folder to put the input file and subfolders from Point #3 above.

Contact:
If user finds error (e.g. segmented fault, etc), please contact Tegar at tegar@alumni.ubc.ca 


Thanks
