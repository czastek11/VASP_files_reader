#ifndef VASP_READ_H 
#define VASP_READ_H

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>

std::vector<double> moving_average(const std::vector<double> data, int window_size);

class VASP_data
{
public:
	VASP_data();
	VASP_data(std::string file_path, int ions, std::string format, bool read_CHGCAR, bool read_LOCPOT, bool read_DOS, bool read_EIGENVAL);
	~VASP_data();
	static arma::mat sorting_positions(arma::mat positions, std::string method);
	std::vector<int> get_mesh_indices(arma::vec pos);
	void read_POSCAR(std::string filename);
	void read_CHGCAR(std::string filename);
	void read_LOCPOT(std::string filename);
	void write_POSCAR(std::string filename);
	double count_total_electrons_double();
	int count_total_electrons();
	void write_potential_averaged_xy_z(std::string filename, std::string period_type, int period);
	void write_potential_averaged_xy_z(std::string filename, std::string period_type);
	arma::vec calc_dipole_moment(arma::vec center, std::vector<int> start, std::vector<int> end);
	void write_potential(std::string filename);
	void read_DOS(std::string filename, int ions, bool spin_orbit, int m);
	arma::mat sum_DOS_types(int atoms_sep_type,int orbitals_sep_type);
	void read_EIGENVAL(std::string filename);
	void write_BS(std::string filename, bool verbose_kpts, bool only_path);
	void write_BS(std::string filename);
	void write_DOS_sum_types(std::string id, const arma::mat& dos_summed, int atoms_sep_type,int orbitals_sep_type, bool header);
	arma::mat get_cell_matrix();
	VASP_data supercell_grid(int rep_x, int rep_y, int rep_z,std::vector<bool> add_vacuum); 
	static double calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R);
	static arma::vec calc_dip_dip_force(arma::vec dip_1, arma::vec dip_2, arma::vec R);


private:
	arma::mat cell_matrix;
	std::vector<int> NGiF;
	std::vector<int> atoms_per_type;
	std::vector<std::string> atom_names;
	std::vector<arma::mat> types_atom_positions;
	arma::mat atom_positions;
	arma::cube charge_density_raw;
	arma::cube charge_density;
	arma::cube potential;
	std::vector<arma::mat> dos_data;
	arma::mat KPOINTS;
	arma::mat BS;
	void read_POSCAR_like(std::string file_name, std::fstream& file);
	bool checkgeo();
	bool checkcharge();
	bool checkpot();
	bool checkdos();
	bool checkKPOINTS();
	bool checkBS();
	
};


#endif // VASP_READ_H




