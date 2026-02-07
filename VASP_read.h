#ifndef VASP_READ_H 
#define VASP_READ_H

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <fstream>
#include <sstream>

std::vector<double> moving_average(const std::vector<double> data, int window_size); 
double calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R);
arma::vec calc_dip_dip_force(arma::vec dip_1, arma::vec dip_2, arma::vec R);
void write_DOS_sum_types(std::string id, const std::vector<std::vector<std::vector<double>>>& dos_summed, const std::vector<std::string>& names);

class VASP_data
{
public:
	VASP_data();
	~VASP_data();
	static arma::mat sorting_positions(arma::mat positions, std::string method);
	std::vector<int> get_mesh_indices(arma::vec pos);
	void read_POSCAR_like(std::string file_name, std::fstream& file);
	void read_CHGCAR(std::string path, std::string body, std::string id);
	void read_LOCPOT(std::string path, std::string body, std::string id);
	double count_total_electrons_double();
	int count_total_electrons();
	void write_potential_averaged_xy_z(std::string id, bool primitive);
	arma::vec calc_dipole_moment(arma::vec center, std::vector<int> start, std::vector<int> end);
	void write_potential(std::string id);
	void read_DOS(std::string id, int ions, std::string format);
	std::vector<std::vector<std::vector<double>>> sum_DOS_types(std::vector<int>& sets);
	arma::mat get_cell_matrix();

private:
	arma::mat cell_matrix;
	std::vector<int> NGiF;
	std::vector<int> atoms_per_type;
	std::vector<arma::mat> types_atom_positions;
	arma::mat atom_positions;
	double*** charge_density_raw;
	double*** charge_density;
	double*** potential;
	std::vector<std::vector<std::vector<double>>> dos_data;
	bool checkgeo();
	bool checkcharge();
	bool checkpot();
	bool checkdos();
	
};


#endif // VASP_READ_H




