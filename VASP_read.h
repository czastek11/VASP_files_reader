#ifndef VASP_READ_H 
#define VASP_READ_H

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <fstream>
#include <sstream>
#include <iomanip>

std::vector<double> moving_average(const std::vector<double> data, int window_size);

class VASP_data
{
public:
	VASP_data();
	VASP_data(std::string file_path, int ions, std::string format, bool read_CHGCAR, bool read_LOCPOT, bool read_DOS);
	~VASP_data();
	static arma::mat sorting_positions(arma::mat positions, std::string method);
	std::vector<int> get_mesh_indices(arma::vec pos);
	void read_CHGCAR(std::string filename);
	void read_LOCPOT(std::string filename);
	double count_total_electrons_double();
	int count_total_electrons();
	void write_potential_averaged_xy_z(std::string filename, std::string period_type, double period);
	void write_potential_averaged_xy_z(std::string filename, std::string period_type);
	arma::vec calc_dipole_moment(arma::vec center, std::vector<int> start, std::vector<int> end);
	void write_potential(std::string filename);
	void read_DOS(std::string filename, int ions, std::string format);
	std::vector<std::vector<std::vector<double>>> sum_DOS_types(std::vector<int>& sets);
	void read_EIGENVAL(std::string filename);
	void write_BS(std::string filename);
	arma::mat get_cell_matrix(); 
	static double calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R);
	static arma::vec calc_dip_dip_force(arma::vec dip_1, arma::vec dip_2, arma::vec R);
	static void write_DOS_sum_types(std::string id, const std::vector<std::vector<std::vector<double>>>& dos_summed, const std::vector<std::string>& names);

private:
	int kpoints;
	int NBANDS;
	arma::mat cell_matrix;
	std::vector<int> NGiF;
	std::vector<int> atoms_per_type;
	std::vector<arma::mat> types_atom_positions;
	arma::mat atom_positions;
	double*** charge_density_raw;
	double*** charge_density;
	double*** potential;
	std::vector<std::vector<std::vector<double>>> dos_data;
	double** BS;
	void read_POSCAR_like(std::string file_name, std::fstream& file);
	bool checkgeo();
	bool checkcharge();
	bool checkpot();
	bool checkdos();
	bool checkBS();
	
};


#endif // VASP_READ_H




