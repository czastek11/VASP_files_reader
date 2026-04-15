#ifndef VASP_READ_H 
#define VASP_READ_H
#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <algorithm>
#include <math.h>

std::vector<double> moving_average(const std::vector<double> data, int window_size);

class VASP_data
{
public:
	VASP_data(); //empty constructor
	VASP_data(std::string file_path, bool read_POSCAR, bool read_CHGCAR, bool read_LOCPOT, bool read_DOS, bool read_EIGENVAL); //
	// file_path the directory where the files are, bool read_XXX true or false to read specific file
	// CHGCAR and LOCPOT already contain all POSCAR informations so it's redundant to enbale reading POSCAr while reading those
	// constructor strarting with loading any selected files with default options
	~VASP_data(); //destructor
	static arma::mat sorting_positions(arma::mat positions, std::string method);
	//sort positions of ions inside each atom type group, available methods:
	//z_rising - sorts by the trhird reduced coordinate ascending, usefull for TMDS
	std::vector<int> get_mesh_indices(arma::vec pos);
	// gives position in mesh grid correspodning to given cartesian psotion or closest fit
	void read_POSCAR(std::string filename);
	//reads all informations from POSCAR file
	void read_bestsqs(std::string filename);
	//converts best sqs output into POSCAR data
	void read_CHGCAR(std::string filename);
	//reads the data from CHGCAR file, automaticlly fills POSCAR file informations
	void read_LOCPOT(std::string filename);
	//reads the data from LOCPOT file, automaticlly fills POSCAR file informations
	void write_POSCAR(std::string filename);
	//writes out POSCAR file with informations stored in object
	double count_total_electrons_double();
	//count total electrons - exact result
	int count_total_electrons();
	//count total electrons - rounded to nearest intiger


	std::vector<double> sum_potential_averaged_xy_z(std::string period_type, int period);
	//when setting otpion "manual", user need to provide number of mesh points to average over
	std::vector<double> sum_potential_averaged_xy_z(std::string period_type);
	// options are "primitive" - window is the size of the cell (only in 3 direction now)
	// "layered" - finds third coordinate of first and third type of first ion (it was used for layered TMDS, this why it's so specific for now)
	// and then take distance between them as window
	// "none" - skip calculating moving average
	std::vector<double> average_potential_over(int direction);

	std::vector<double> moving_average_potential_over(std::vector<double>, int direction, std::string period_type);
	std::vector<double> moving_average_potential_over(std::vector<double>, int direction, std::string period_type, int period);
	std::vector<double> moving_average_potential_over(std::vector<double> av_pot, int direction, std::string period_type, int ion1, int ion2);

	void write_potential_over(std::string filename, std::vector<double> potential_av, int direction);

	//methods to average potential in xy (first and second lattice vector directions) direction and then do moving average over z direction. The range of moving average is defined
	// by period type and period (number of points in mesh in third lattice vector direction)
	void write_potential_z(std::string filename, std::vector<double> potential_z);
	//writes out averaged potential that you got from previous method
	arma::vec calc_dipole_moment(arma::vec center, std::vector<int> start, std::vector<int> end);
	//calculate the diople moment relative to center position, start and end are starting and ending range for mesh points to calculate over
	void write_potential(std::string filename);
	//write out potential in CHGCAR in format that is easy to plot : x y z value
	void read_DOS(std::string filename, bool spin_orbit);
	//read DOS to object. SO version not yet implemented
	arma::mat sum_DOS_types(int atoms_sep_type,int orbitals_sep_type);
	//sum DOS contributions in specific way for atoms and orbitals:
	// 0 - sum all ions/orbitals(sum over l and m)
	// 1 - sum all atoms(same atom ions)/orbital types (sum over m)
	// 2 - don't sum at all, get seperate contribution for each ion/orbital 
	void read_EIGENVAL(std::string filename);
	//read information from EIGENVAL: kpoints, their indexes, their wieght, all the band energies for each and all their occupations
	void read_BS(std::string filename, bool header, bool verbose_kpts);
	void write_BS(std::string filename, bool verbose_kpts, bool only_path);
	//write out band structure in more compact way for graphing, verbose_kpt is to enable writing out kpoints values
	// in reduced form (not cartesian), only_path skips writing points that have integration weight, so that only band structure path is written
	void write_BS(std::string filename);
	// write out band structure with default settings: false false 
	void write_DOS_sum_types(std::string id, const arma::mat& dos_summed, int atoms_sep_type,int orbitals_sep_type, bool header);
	// write down the summed DOS from sum_DOS_types. You need to provide the same settings for atoms_sep_type
	// and orbitals_sep_type as before. header controls if the header with the name of the columns is written into the file
	arma::mat get_cell_matrix();
	arma::mat get_BS();
	arma::mat get_occupations();
	// get methods
	arma::rowvec find_kpoint_energy(arma::rowvec kpt, bool weight, int& index);
	// get rowvector of all bands energy at specified kpoint, index return what index that kpoints is. Weight controls if the moethod
	// looks trough weighted kpoints
	double find_band_extremum(int band_index, bool weight, int& kpt_index, bool maxormin);
	// get the value and kpoint index corresponding to global extremum of given band. maxormin controls if it's maximum(true)
	//or minimum(false)
	int find_valence_band();
	//gives the index of highest occupied band
	VASP_data supercell_grid(int rep_x, int rep_y, int rep_z,std::vector<double> add_vacuum); 
	//generate supercell as new VASP_data object. rep_x,rep_y,rep_z controls how many times in given direction the cell is multiplied
	// vacuum is array of double values to set how many multiplies of vectors are to be added as vacuum

    void alloy_geometry(VASP_data data2, double percentage, std::vector<std::string> mixed_atom_names1, std::vector<std::string> mixed_atom_names2, std::string filename);
	//generate alloy geometry as new VASP_data object. this and data2 are the two parent structures, percentage is the percentage of data2 in the alloy
	// then write it in mcsqs format
	// make sure that the order of postions that are to be mixed is the same in both structures, also make sure that all atoms are either in both compounds
	// or listed in with its to be mixed counterpart as mixed_atom_names1 and mixed_atom_names2
	// the pair of mixed atoms should in the same positons in the vectors mixed_atom_names1 and mixed_atom_names2,
	// so that the first element of mixed_atom_names1 is mixed with the first element of mixed_atom_names2 and so on


	static double calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R);
	// calculate potential between two dipoles
	// dip1 and dip2 are two dipoles and R is the vector of their separation R=r_1-r_2
	static arma::vec calc_dip_dip_force(arma::vec dip_1, arma::vec dip_2, arma::vec R);
	// calculate force from the potential betwee  two dipoles
	// dip1 and dip2 are two dipoles and R is the vector of their separation R=r_1-r_2

private:
	arma::mat cell_matrix; //cell base vecotrs as written in POSCAR file, transpose before using for transformation from fractional to cartesian
	std::vector<int> NGiF; //dimension of point mesh, it can be different between CHGCAR and LOCPOT so proceed with caution when reading both
	std::vector<int> atoms_per_type; // how many ions per atom there is
	std::vector<std::string> atom_names; // names of atoms in system
	std::vector<arma::mat> types_atom_positions; // postion of ions separated into groups of the same atoms
	arma::mat atom_positions; // all atom postions as written in POSCAR, by default it is in fractional coordinates
	arma::cube charge_density_raw; // charge density from CHGCAR before normalisation
	arma::cube charge_density; // normalised charge density from CHGCAR
	arma::cube potential; // potential from LOCPOT
	std::vector<arma::mat> dos_data; //energy and DOS read from DOSCAR
	arma::mat KPOINTS; //list of kpoints and their weight
	arma::mat BS; // energies of all bands and kpoints
	arma::mat occupations; // occupations for each band and kpoint
	void read_POSCAR_like(std::string file_name, std::fstream& file); //reading POSCAR like header in POSCAR,CHGCAR,LOCPOT
	bool checkgeo();
	bool checkmesh();
	bool checkcharge();
	bool checkpot();
	bool checkdos();
	bool checkKPOINTS();
	bool checkBS();
	// safety check to be sure if the data we want to use is loaded already
	
};


#endif // VASP_READ_H




