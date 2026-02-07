#include "VASP_read.h"

std::vector<double> moving_average(std:: vector<double> data, int window_size)
{
	std::vector<double> averaged(data.size(), 0.0);
	int half_window = window_size / 2;
	for (int i = 0; i < data.size(); i++)
	{
		double sum = 0.0;
		int count = 0;
		for (int j = -half_window; j <= half_window; j++)
		{
			int idx = i + j;
			if (idx >= 0 && idx < data.size())
			{
				sum += data[idx];
				count++;
			}
		}
		averaged[i] = sum / count;
	}
	return averaged;
}

double calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R)
{
	//double epsilon_0 = 8.854187817e-12; // vacuum permittivity in F/m
	// Using atomic units where 1/(4pieps0) = 1
	double R_len = arma::norm(R);
	arma::vec r_hat = R / R_len;
	double potential = (arma::dot(dip_1, dip_2) - 3 * arma::dot(dip_1, r_hat) * arma::dot(dip_2, r_hat)) / (R_len * R_len * R_len);
	return potential;
}

arma::vec calc_dip_dip_force(arma::vec dip_1, arma::vec dip_2, arma::vec R) //check the derivation : -grad_r calc_dip_dip_potential
{
	//double epsilon_0 = 8.854187817e-12; // vacuum permittivity in F/m
	// Using atomic units where 1/(4pieps0) = 1
	double R_len = arma::norm(R);
	arma::vec r_hat = R / R_len;
	arma::vec grad_r_3 = -3 * r_hat / (R_len * R_len * R_len * R_len);
	arma::vec term1 = grad_r_3 * arma::dot(dip_1, dip_2);
	double term2_1 = arma::dot(dip_1, r_hat);
	double term2_2 = arma::dot(dip_2, r_hat);
	arma::vec term2_1_grad = (dip_1 - r_hat * term2_1) / R_len;
	arma::vec term2_2_grad = (dip_2 - r_hat * term2_2) / R_len;
	arma::vec term2 = -3 * (term2_1_grad * term2_2 + term2_1 * term2_2_grad) / (R_len * R_len * R_len) - 3 * term2_1 * term2_2 * grad_r_3;
	arma::vec force = -(term1 + term2);
	return force;
}

void write_DOS_sum_types(std::string id, const std::vector<std::vector<std::vector<double>>>& dos_summed, const std::vector<std::string>& names)
{
	std::fstream file;

	int num_types = dos_summed.size();
	int NDOS = dos_summed[0].size();
	for (int type_id = 0; type_id < num_types; type_id++)
	{
		std::string file_name = "workspace\\" + id + "_" + names[type_id] + "_DOS_sum.txt";
		file.open(file_name, std::ios::out);
		file.precision(12);
		file << "# DOS for type: " << names[type_id] << "\n";
		file << "# Energy (eV)   DOS\n";
		for (int i = 0; i < NDOS; i++)
		{
			file << dos_summed[type_id][i][0] << " " << dos_summed[type_id][i][1] << "\n";
		}
		file << "\n\n";
		file.close();
	}

}

VASP_data::VASP_data() :  NGiF(), atoms_per_type(), types_atom_positions(), atom_positions(), charge_density_raw(nullptr), charge_density(nullptr), potential(nullptr), dos_data()
{
	cell_matrix = arma::mat(3, 3, arma::fill::zeros);
}

VASP_data::~VASP_data()
{
	if (charge_density_raw != nullptr)
	{
		for (int i = 0; i < NGiF[0]; i++)
		{
			for (int j = 0; j < NGiF[1]; j++)
			{
				delete[] charge_density_raw[i][j];
			}
			delete[] charge_density_raw[i];
		}
		delete[] charge_density_raw;
	}
	if (charge_density != nullptr)
	{
		for (int i = 0; i < NGiF[0]; i++)
		{
			for (int j = 0; j < NGiF[1]; j++)
			{
				delete[] charge_density[i][j];
			}
			delete[] charge_density[i];
		}
		delete[] charge_density;
	}
	if (potential != nullptr)
	{
		for (int i = 0; i < NGiF[0]; i++)
		{
			for (int j = 0; j < NGiF[1]; j++)
			{
				delete[] potential[i][j];
			}
			delete[] potential[i];
		}
		delete[] potential;
	}
}

bool VASP_data::checkgeo()
{
	if (NGiF.size() != 3)
	{
		std::cerr << "Error: Either you did not load data first or did something wrong that NGiF does not have 3 elements." << std::endl;
		throw std::runtime_error("Error: NGiF does not have 3 elements");
		return false;
	}
	else return true;
}

bool VASP_data::checkcharge()
{
	if (charge_density_raw == nullptr)
	{
		std::cerr << "Error: charge density data not loaded. Please load data before using charge density." << std::endl;
		throw std::runtime_error("Error: charge density data not loaded");
		return false;
	}
	else return true;
}

bool VASP_data::checkpot()
{
	if (potential == nullptr)
	{
		std::cerr << "Error: potential data not loaded. Please load data before using potential." << std::endl;
		throw std::runtime_error("Error: potential data not loaded");
		return false;
	}
	else return true;
}

bool VASP_data::checkdos()
{
	if (dos_data.size() == 0)
	{
		std::cerr << "Error: DOS data not loaded. Please load data before using DOS data." << std::endl;
		throw std::runtime_error("Error: DOS data not loaded");
		return false;
	}
	else return true;
}

arma::mat VASP_data::sorting_positions(arma::mat positions, std::string method)
{
	if (method == "z_rising")
	{
		arma::uvec sorted_indices = arma::sort_index(positions.row(2).t(), "ascend");
		return positions.cols(sorted_indices);
	}
	else
	{
		std::cerr << "Unknown sorting method: " << method << std::endl;
		throw std::runtime_error("Unknown sorting method");
	}

}

std::vector<int> VASP_data::get_mesh_indices(arma::vec pos)
{

	if (checkgeo())
	{
		arma::vec a = cell_matrix.row(0).t();
		arma::vec b = cell_matrix.row(1).t();
		arma::vec c = cell_matrix.row(2).t();
		std::vector<int> indices(3);
		arma::mat transform_matrix = arma::mat(3, 3);
		transform_matrix.col(0) = a / NGiF.at(0);
		transform_matrix.col(1) = b / NGiF.at(1);
		transform_matrix.col(2) = c / NGiF.at(2);
		//transform_matrix.print("Transform matrix:");
		arma::vec fractional_coords = arma::inv(transform_matrix) * pos;
		for (int i = 0; i < 3; i++)
		{
			indices[i] = static_cast<int>(floor(fractional_coords(i) + 0.5));
			if (indices[i] < 0) indices[i] = 0;
			if (indices[i] >= NGiF[i]) indices[i] = NGiF[i] - 1;
		}
		return indices;
	}
	else return std::vector<int>();
	
}

void VASP_data::read_POSCAR_like(std::string file_name, std::fstream& file)
{
	//fstream file;
	file.open(file_name, std::ios::in);
	if (!file.is_open())
	{
		std::cerr << "Error opening file: " << file_name << std::endl;
		throw std::runtime_error("Error opening file");
	}
	std::string line, word;
	double number_read;
	int int_read;
	getline(file, line); // Skip the first line (comment)
	getline(file, line);

	for (int i = 0; i < 9; i++)
	{
		file >> number_read;
		cell_matrix(i / 3, i % 3) = number_read;
	}
	getline(file, line);
	getline(file, line); // Skip the line after cell matrix
	getline(file, line);
	int num_atom = 0;
	std::stringstream ss(line);
	//vector<int> atoms_per_type;

	atoms_per_type.clear();
	while (ss >> int_read)
	{
		atoms_per_type.push_back(int_read);
		num_atom += int_read;
	}
	//arma::mat atom_positions = arma::mat(3, 0);
	//vector<arma::mat> types_atom_positions;
	getline(file, line); //skip direct
	types_atom_positions.clear();
	for (int i = 0; i < atoms_per_type.size(); i++)
	{
		types_atom_positions.push_back(arma::mat(3, atoms_per_type[i]));
		for (int j = 0; j < atoms_per_type[i]; j++)
		{
			getline(file, line);
			std::stringstream ss2(line);
			double x, y, z;
			ss2 >> x >> y >> z;
			arma::vec pos = { x, y, z };
			types_atom_positions[i].col(j) = pos;
		}
		types_atom_positions[i] = sorting_positions(types_atom_positions[i], "z_rising");
		atom_positions = arma::join_rows(atom_positions, types_atom_positions[i]);
	}
	getline(file, line); // Skip the line before charge density data

	getline(file, line);

	std::stringstream ss2(line);
	NGiF.clear();
	for (int i = 0; i < 3; i++)
	{
		ss2 >> int_read;
		NGiF.push_back(int_read);
	}
}

void VASP_data::read_CHGCAR(std::string path, std::string body, std::string id)
{
	std::string file_name = path + body + id;
	std::fstream file;
	if (charge_density_raw != nullptr)
	{
		for (int i = 0; i < NGiF[0]; i++)
		{
			for (int j = 0; j < NGiF[1]; j++)
			{
				delete[] charge_density_raw[i][j];
				delete[] charge_density[i][j];
			}
			delete[] charge_density_raw[i];
			delete[] charge_density[i];
		}
		delete[] charge_density_raw;
		delete[] charge_density;
	}
	read_POSCAR_like(file_name, file);

	int total_grid_points = NGiF[0] * NGiF[1] * NGiF[2];

	
	charge_density_raw = new double** [NGiF[0]];
	charge_density = new double** [NGiF[0]];
	for (int i = 0; i < NGiF[0]; i++)
	{
		charge_density_raw[i] = new double* [NGiF[1]];
		charge_density[i] = new double* [NGiF[1]];
		for (int j = 0; j < NGiF[1]; j++)
		{
			charge_density_raw[i][j] = new double[NGiF[2]];
			charge_density[i][j] = new double[NGiF[2]];
		}
	}
	for (int i = 0; i < total_grid_points; i++)
	{
		int ix = i % NGiF[0];
		int iy = (i / NGiF[0]) % NGiF[1];
		int iz = i / (NGiF[0] * NGiF[1]);
		file >> charge_density_raw[ix][iy][iz];
		//charge_density[ix][iy][iz] = charge_density_raw[ix][iy][iz] / total_grid_points / cell_volume;
		charge_density[ix][iy][iz] = charge_density_raw[ix][iy][iz] / total_grid_points;
	}
	file.close();
}

void VASP_data::read_LOCPOT(std::string path, std::string body, std::string id)
{
	std::string file_name = path + body + id;
	std::fstream file;
	if (potential != nullptr)
	{
		for (int i = 0; i < NGiF[0]; i++)
		{
			potential[i] = new double* [NGiF[1]];
			for (int j = 0; j < NGiF[1]; j++)
			{
				delete[] potential[i][j];
			}
			delete[] potential[i];
		}
		delete[] potential;
	}
	read_POSCAR_like(file_name, file);


	int total_grid_points = NGiF[0] * NGiF[1] * NGiF[2];
	
	potential = new double** [NGiF[0]];
	for (int i = 0; i < NGiF[0]; i++)
	{
		potential[i] = new double* [NGiF[1]];
		for (int j = 0; j < NGiF[1]; j++)
		{
			potential[i][j] = new double[NGiF[2]];
		}
	}
	for (int i = 0; i < total_grid_points; i++)
	{
		int ix = i % NGiF[0];
		int iy = (i / NGiF[0]) % NGiF[1];
		int iz = i / (NGiF[0] * NGiF[1]);
		file >> potential[ix][iy][iz];
		//potentai lcan have up to 4 datasets for spin polarized calculations, but we only read the first one here
		//current example has 4 but magmom is set to 0 for all atoms, so the other 3 datasets are mostly zeros
	}
	file.close();
}

double VASP_data::count_total_electrons_double()
{
	if (checkcharge())
	{
		double total_electrons = 0;
		for (int i = 0; i < atoms_per_type.size(); i++)
		{
			// sanity check
			double nel = 0.0;
			for (int i = 0; i < NGiF[0]; i++)
			{
				for (int j = 0; j < NGiF[1]; j++)
				{
					for (int k = 0; k < NGiF[2]; k++)
					{
						nel += charge_density_raw[i][j][k];
					}
				}
			}
			nel /= (NGiF[0] * NGiF[1] * NGiF[2]);
			return nel; // rounding to nearest integer
		}
		return total_electrons;
	}
	else return 0.0;
}

int VASP_data::count_total_electrons()
{
	int total_electrons = 0;
	return static_cast<int>(count_total_electrons_double() + 0.5); // rounding to nearest integer
}

void VASP_data::write_potential_averaged_xy_z(std::string id, bool primitive)
{
	if (checkpot())
	{
		std::vector<double> potential_z(NGiF[2], 0.0);
		int total_xy_points = NGiF[0] * NGiF[1];
		for (int k = 0; k < NGiF[2]; k++) { // z index
			double sum = 0.0;
			for (int i = 0; i < NGiF[0]; i++) { // x index
				for (int j = 0; j < NGiF[1]; j++) { // y index
					sum += potential[i][j][k];
				}
			}
			potential_z[k] = sum / total_xy_points;
		}

		// average potential in z over moving window. Window size is one periodic layer thickness.
		int window = 1;
		if (primitive)
		{
			window = NGiF[2];
		}
		else
		{
			arma::vec distance = cell_matrix.t() * (atom_positions.col(2) - atom_positions.col(0));
			window = get_mesh_indices(distance)[2];
		}

		potential_z = moving_average(potential_z, window);


		std::string file_name = "workspace\\" + id + "_potential_z.txt";
		std::fstream file;
		file.open(file_name, std::ios::out);
		double z_real;
		for (int i = 0; i < NGiF[2]; i++)
		{
			z_real = i * cell_matrix(2, 2) / NGiF[2];
			file << z_real << " " << potential_z[i] << "\n";

		}
		file.close();
	}
}

arma::vec VASP_data::calc_dipole_moment(arma::vec center, std::vector<int> start, std::vector<int> end)
{
	if (checkcharge())
	{
		arma::vec dipole_mom(3, arma::fill::zeros);
		arma::vec pos(3, arma::fill::zeros);
		arma::vec a = cell_matrix.row(0).t();
		arma::vec b = cell_matrix.row(1).t();
		arma::vec c = cell_matrix.row(2).t();
		for (int i = start[0]; i < end[0]; i++)
		{
			for (int j = start[1]; j < end[1]; j++)
			{
				for (int k = start[2]; k < end[2]; k++)
				{
					pos = i * a / NGiF[0] + j * b / NGiF[1] + k * c / NGiF[2];
					dipole_mom += -charge_density[i][j][k] * (pos - center);
				}
			}
		}
		return dipole_mom;
	}
	else return arma::vec(3, arma::fill::zeros);
}

void VASP_data::write_potential(std::string id)
{
	if (checkpot())
	{
		arma::vec a = cell_matrix.col(0), b = cell_matrix.col(1), c = cell_matrix.col(2), pos;
		std::string file_name = "workspace\\" + id + "_potential.txt";
		std::fstream file;
		file.open(file_name, std::ios::out);
		for (int i = 0; i < NGiF[0]; i++)
		{
			for (int j = 0; j < NGiF[1]; j++)
			{
				for (int k = 0; k < NGiF[2]; k++)
				{
					pos = i * 1.0 / NGiF[0] * a + j * 1.0 / NGiF[1] * b + k * 1.0 / NGiF[2] * c;
					file << pos(0) << " " << pos(1) << " " << pos(2) << " " << potential[i][j][k] << "\n";
				}
			}
		}
		file.close();
	}
}

void VASP_data::read_DOS(std::string id, int ions, std::string format)
{
	if (format == "LORBIT=11,no_SO")
	{
		std::fstream file;
		file.open(id, std::ios::in);
		if (!file.is_open())
		{
			std::cerr << "Error opening file: " << id << std::endl;
			throw std::runtime_error("Error opening file");
		}
		else
		{
			std::string line;
			int NDOS;
			double pom;
			for (int i = 0; i < 5; i++)// skip first 5 lines
			{
				getline(file, line);
			}
			getline(file, line); // sixth line contains number of energy points
			std::stringstream ss(line);
			for (int i = 0; i < 3; i++)// third number is NDOS
			{
				ss >> pom;
			}
			NDOS = static_cast<int>(pom);

			// allocate memory for DOS data
			dos_data.clear();
			std::vector<std::vector<std::vector<double>>> dos_data(ions, std::vector<std::vector<double>>(NDOS, std::vector<double>(10, 0.0)));
			//skip first set with two columns (total DOS)
			for (int i = 0; i < NDOS; i++) getline(file, line);

			// read DOS data
			for (int ion = 0; ion < ions; ion++)
			{
				getline(file, line); // skip header line for each ion
				for (int i = 0; i < NDOS; i++)
				{
					getline(file, line);
					std::stringstream ss2(line);
					for (int j = 0; j < 10; j++)
					{
						ss2 >> dos_data[ion][i][j];
					}
				}
			}
			file.close();
		}
	}
	else
	{
		std::cerr << "Unknown DOS format: " << format << std::endl;
		throw std::runtime_error("Unknown DOS format");
	}
}

std::vector<std::vector<std::vector<double>>> VASP_data::sum_DOS_types(std::vector<int>& sets)
{
	if (checkdos())
	{
		int ions = dos_data.size(), NDOS = dos_data[0].size(), num_types = sets.size();
		// sanity check
		int total_ions = 0;
		if (ions == 0) return { {} };
		for (int set_size : sets) total_ions += set_size;
		if (total_ions != ions)
		{
			std::cerr << "Error: Sum of sets doesn't match total ions" << std::endl;
			return { {} };
		}

		std::vector<std::vector<std::vector<double>>> results(num_types, std::vector<std::vector<double>>(NDOS, std::vector<double>(2, 0.0)));
		int ion_index = 0;
		for (int type_id = 0; type_id < num_types; type_id++)
		{
			int set_size = sets[type_id];

			for (int i = 0; i < set_size; i++)
			{
				for (int j = 0; j < NDOS; j++)
				{
					results[type_id][j][0] = dos_data[ion_index][j][0]; // energy value
					// sum over all orbitals (index 1 to 9)
					for (int k = 1; k < 10; k++)
					{
						results[type_id][j][1] += dos_data[ion_index][j][k];
					}
				}
				ion_index++;
			}
		}
		return results;
	}
	else return {};
}

arma::mat VASP_data::get_cell_matrix()
{
	return cell_matrix;
}