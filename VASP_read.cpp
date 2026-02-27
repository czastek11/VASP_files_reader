#include "VASP_read.h"

std::vector<double> moving_average(std::vector<double> data, int window_size)
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

void multiply_cell_in_direction(std::vector<arma::mat>& cart_types, arma::vec base_vec, int rep, bool add_vacuum_below, bool add_vacuum_above)
{
	arma::mat cart_new;
	std::vector<arma::mat> starting_cart_types = cart_types;
	std::vector<int> starting_num_types;
	for(int i =0 ;i<cart_types.size();i++) starting_num_types.push_back(cart_types.at(i).n_cols);
	for (int i = 0; i < rep; i++)
	{
		if(add_vacuum_below && i==0)
		{
			for(int t = 0; t < cart_types.size(); t++)
			{
				
				for(int pos = 0 ; pos < starting_num_types.at(t); pos++)
				{
					starting_cart_types[t].col(pos) += base_vec;
				}
				
			}
			cart_types = starting_cart_types;
		}
		else if(i>0)
		{
			for(int t = 0; t < cart_types.size(); t++)
			{
				cart_new = starting_cart_types[t];
				for(int pos = 0 ; pos < starting_num_types.at(t); pos++)
				{
					cart_new.col(pos) += base_vec * i;
				}
				cart_types[t] = arma::join_rows(cart_types[t], cart_new);
			}
		}
	}
}

VASP_data::VASP_data() :
 	NGiF(), atoms_per_type(), types_atom_positions(), atom_positions(),
	charge_density_raw(), charge_density(), potential(), dos_data(), 
	KPOINTS(), BS(), atom_names(), occupations()
{
	cell_matrix = arma::mat(3, 3, arma::fill::zeros);
}

VASP_data::VASP_data(std::string file_path, int ions, std::string format, bool read_CHGCAR, bool read_LOCPOT, bool read_DOS, bool read_EIGENVAL) : VASP_data()
{
	if (read_CHGCAR) this->read_CHGCAR(file_path + "CHGCAR");
	if (read_LOCPOT) this->read_LOCPOT(file_path + "LOCPOT");
	if (read_DOS) this->read_DOS(file_path + "DOSCAR", ions, false, 2); // currently only support non-spin-orbit DOSCAR with s,p,d,f orbitals separated. Will add more options in the future if needed.
	if (read_EIGENVAL) this->read_EIGENVAL(file_path + "EIGENVAL");
}

VASP_data::~VASP_data()
{
	
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
	if (charge_density_raw.is_empty())
	{
		std::cerr << "Error: charge density data not loaded. Please load data before using charge density." << std::endl;
		throw std::runtime_error("Error: charge density data not loaded");
		return false;
	}
	else return true;
}

bool VASP_data::checkpot()
{
	if (potential.is_empty())
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

bool VASP_data::checkKPOINTS()
{
	if (KPOINTS.is_empty())
	{
		std::cerr << "Error: KPOINTS data not loaded. Please load data before using KPOINTS data." << std::endl;
		throw std::runtime_error("Error: KPOINTS data not loaded");
		return false;
	}
	else return true;
}

bool VASP_data::checkBS()
{
	if (BS.is_empty())
	{
		std::cerr << "Error: Band structure data not loaded. Please load data before using band structure data." << std::endl;
		throw std::runtime_error("Error: Band structure data not loaded");
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

void VASP_data::read_POSCAR_like(std::string filename, std::fstream& file)
{
	//fstream file;
	double number_read;
	int int_read, num_atom = 0;
	std::string line, word;

	file.open(filename, std::ios::in);
	if (!file.is_open())
	{
		std::cerr << "Error opening file: " << filename << std::endl;
		throw std::runtime_error("Error opening file");
	}
	//clear previous data after successfully opening the file, to avoid partial clearing if file opening fails
	atoms_per_type.clear();
	types_atom_positions.clear();
	atom_positions.clear();
	NGiF.clear();
	
	getline(file, line); // Skip the first line (comment)
	getline(file, line);

	for (int i = 0; i < 9; i++)
	{
		file >> number_read;
		cell_matrix(i / 3, i % 3) = number_read;
	}

	getline(file, line);

	getline(file, line); // read the line with atom names, which is optional in VASP POSCAR format. If not present, we will just have an empty atom_names vector and rely on atoms_per_type for counting.
	std::stringstream ss4(line);
	atom_names.clear();
	while (ss4 >> word)
	{
		atom_names.push_back(word);
	}

	getline(file, line);
	std::stringstream ss(line);
	while (ss >> int_read)
	{
		atoms_per_type.push_back(int_read);
		num_atom += int_read;
	}

	getline(file, line); //skip direct

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
	for (int i = 0; i < 3; i++)
	{
		ss2 >> int_read;
		NGiF.push_back(int_read);
	}
}

void VASP_data::read_POSCAR(std::string filename)
{
	std::fstream file;
	read_POSCAR_like(filename,file);
	if(file.is_open()) file.close();
	else 
	{
		std::cout << "POSCAR not probably not read!\n";
		std::cerr << "POSCAR not probably not read!\n";
	}
}

void VASP_data::read_CHGCAR(std::string filename)
{
	std::fstream file;
	read_POSCAR_like(filename, file);

	int total_grid_points = NGiF[0] * NGiF[1] * NGiF[2];

	
	charge_density_raw = arma::cube(NGiF[0], NGiF[1], NGiF[2], arma::fill::zeros);
	charge_density = arma::cube(NGiF[0], NGiF[1], NGiF[2], arma::fill::zeros);

	for (int i = 0; i < total_grid_points; i++)
	{
		int ix = i % NGiF[0];
		int iy = (i / NGiF[0]) % NGiF[1];
		int iz = i / (NGiF[0] * NGiF[1]);
		file >> charge_density_raw(ix, iy, iz);
		//charge_density[ix][iy][iz] = charge_density_raw[ix][iy][iz] / total_grid_points / cell_volume;
		charge_density(ix, iy, iz) = charge_density_raw(ix, iy, iz) / total_grid_points;
	}
	file.close();
}

void VASP_data::read_LOCPOT(std::string filename)
{
	std::fstream file;
	read_POSCAR_like(filename, file);


	int total_grid_points = NGiF[0] * NGiF[1] * NGiF[2];
	int ix = 0, iy = 0, iz = 0;
	
	potential = arma::cube(NGiF[0], NGiF[1], NGiF[2], arma::fill::zeros);
	for (int i = 0; i < total_grid_points; i++)
	{
		ix = i % NGiF[0];
		iy = (i / NGiF[0]) % NGiF[1];
		iz = i / (NGiF[0] * NGiF[1]);
		file >> potential(ix,iy,iz);
		//potentai lcan have up to 4 datasets for spin polarized calculations, but we only read the first one here
		//current example has 4 but magmom is set to 0 for all atoms, so the other 3 datasets are mostly zeros
	}
	file.close();
}

void VASP_data::write_POSCAR(std::string filename)
{
	std::filesystem::path fullPath = std::filesystem::path("workspace") / (filename + "_POSCAR");
	std::fstream file;
	file.open(fullPath, std::ios::out);
	if (!file.is_open())
	{
		std::cerr << "Error opening file for writing: " << fullPath << std::endl;
		throw std::runtime_error("Error opening file for writing");
	}
	for(int i=0 ; i< atoms_per_type.size(); i++)
	{
		file << atom_names[i] << atoms_per_type[i] << " ";
	}
	file<<std::fixed << std::setprecision(16);
	file<< "\n";
	file << "   "<<1.0<<"\n";
	for (int i = 0; i < 3; i++)
	{
		file<<" ";
		for (int j = 0; j < 3; j++)
		{
			if(cell_matrix(i,j)>=0) file<<" ";
			file << "   " << cell_matrix(i, j);
		}
		file << "\n";
	}

	for(int i =0; i<atom_names.size(); i++)
	{
		file << "   "<< atom_names[i] ;
	}
	file << "\n";

	for(int i = 0; i<atoms_per_type.size();i++)
	{
		file<<"     "<<atoms_per_type.at(i);
	}
	file<<"\n";

	file<<"Direct\n";
	for (int i = 0; i < types_atom_positions.size(); i++)
	{
		for (int j = 0; j < types_atom_positions[i].n_cols; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if(types_atom_positions[i](k,j)>=0) file<<" ";
				file << " " << std::fixed << std::setprecision(16) << types_atom_positions[i](k, j);
			}
			file << "\n";
		}
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
						nel += charge_density_raw(i,j,k);
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

std::vector<double> VASP_data::sum_potential_averaged_xy_z(std::string period_type, int period)
{
	if (checkpot())
	{
		std::vector<double> potential_z(NGiF[2], 0.0);
		int total_xy_points = NGiF[0] * NGiF[1];
		double sum;
		for (int k = 0; k < NGiF[2]; k++)  // z index
		{ 
			sum = 0.0;
			for (int i = 0; i < NGiF[0]; i++)// x index
			{ 
				for (int j = 0; j < NGiF[1]; j++) // y index
				{ 
					sum += potential(i,j,k);
				}
			}
			potential_z[k] = sum / total_xy_points;
		}

		// average potential in z over moving window. Window size is one periodic layer thickness.
		int window = 1;
		if (period_type=="primitive")
		{
			window = NGiF[2];
		}
		else if (period_type == "manual")
		{
			window = period;
		}
		else if (period_type == "layered")
		{
			arma::vec distance = cell_matrix.t() * (atom_positions.col(2) - atom_positions.col(0));
			window = get_mesh_indices(distance)[2];
		}

		if (period_type!="none") potential_z = moving_average(potential_z, window);

		return potential_z;
	}
	else return {};
}

std::vector<double> VASP_data::sum_potential_averaged_xy_z(std::string period_type)
{
	if(period_type!="manual") return sum_potential_averaged_xy_z(period_type, 0.0);
	else
	{
		std::cerr << "Error: For manual period type, you must provide a period value." << std::endl;
		throw std::runtime_error("Error: Missing period value for manual period type");
		return {};
	}
}

void VASP_data::write_potential_z(std::string filename, std::vector<double> potential_z)
{
	std::filesystem::path fullPath = std::filesystem::path("workspace") / (filename + "_potential_z.txt");
	std::fstream file;
	file.open(fullPath, std::ios::out);
	double z_real;
	for (int i = 0; i < NGiF[2]; i++)
	{
		z_real = i * cell_matrix(2, 2) / NGiF[2];
		file << z_real << " " << potential_z[i] << "\n";
	}
	file.close();
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
					dipole_mom += -charge_density(i,j,k) * (pos - center);
				}
			}
		}
		return dipole_mom;
	}
	else return arma::vec(3, arma::fill::zeros);
}

void VASP_data::write_potential(std::string filename)
{
	if (checkpot())
	{
		arma::vec a = cell_matrix.col(0), b = cell_matrix.col(1), c = cell_matrix.col(2), pos;
		//std::string full_filename = "workspace\\" + filename + "_potential.txt";
		std::filesystem::path fullPath = std::filesystem::path("workspace") / (filename + "_potential.txt");
		std::fstream file;
		file.open(fullPath, std::ios::out);
		for (int i = 0; i < NGiF[0]; i++)
		{
			for (int j = 0; j < NGiF[1]; j++)
			{
				for (int k = 0; k < NGiF[2]; k++)
				{
					pos = i * 1.0 / NGiF[0] * a + j * 1.0 / NGiF[1] * b + k * 1.0 / NGiF[2] * c;
					file << pos(0) << " " << pos(1) << " " << pos(2) << " " << potential(i,j,k) << "\n";
				}
			}
		}
		file.close();
	}
}

void VASP_data::read_DOS(std::string filename, int ions, bool spin_orbit, int m)
{
	if (!spin_orbit)
	{
		std::fstream file;
		file.open(filename, std::ios::in);
		if (!file.is_open())
		{
			std::cerr << "Error opening file: " << filename << std::endl;
			throw std::runtime_error("Error opening file");
		}
		else
		{
			std::string line;
			int l_tot = (m+1)*(m+1);
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
			arma::mat ion_dos; 
			//skip first set with two columns (total DOS)
			for (int i = 0; i < NDOS; i++) getline(file, line);

			// read DOS data
			for (int ion = 0; ion < ions; ion++)
			{
				ion_dos = arma::mat(NDOS, l_tot+1, arma::fill::zeros); // energy + l_tot orbitals 1s,3p,5d,7f,...
				getline(file, line); // skip header line for each ion
				for (int i = 0; i < NDOS; i++)
				{
					
					getline(file, line);
					std::stringstream ss2(line);
					for (int j = 0; j < l_tot+1; j++)
					{
						ss2 >> ion_dos(i,j);
					}
				}
				dos_data.push_back(ion_dos);
			}
			file.close();
		}
	}
	//else
	//{
	//	std::cerr << "Unknown DOS format: " << format << std::endl;
	//	throw std::runtime_error("Unknown DOS format");
	//}
}

arma::mat VASP_data::sum_DOS_types(int atoms_sep_type,int orbitals_sep_type)
{
	if (checkdos())
	{
		int ions = dos_data.size(), NDOS = dos_data[0].n_rows, atom_types= atoms_per_type.size();
		int l_tot = dos_data[0].n_cols -1 , m = sqrt(l_tot)-1;
		arma::mat results;
		std::vector<int> atom_sets;
		std::vector<int> orb_sets;

		if(atoms_sep_type == 0) atom_sets.push_back(ions);
		else if (atoms_sep_type == 1) atom_sets = atoms_per_type;
		else if(atoms_sep_type == 2) for(int i=0 ; i<ions ; i++) atom_sets.push_back(1);

		if(orbitals_sep_type == 0) orb_sets.push_back(l_tot);
		else if (orbitals_sep_type == 1) for(int i=0 ; i<=m ; i++) orb_sets.push_back(i*2+1);
		else if(orbitals_sep_type == 2) for(int i=0 ; i<l_tot ; i++) orb_sets.push_back(1);

		int at_col = atom_sets.size(), orb_col= orb_sets.size();
		int tot = at_col * orb_col;
		results = arma::mat(NDOS, 1 + tot, arma::fill::zeros); // first column for energy, next columns for all sets combinations
		results.col(0) = dos_data[0].col(0);
		int ion_index = 0, orb_index = 0, col = 0;
		for(int ion_set = 0 ; ion_set < atom_sets.size(); ion_set ++)
		{
			for(int ion = 0 ; ion < atom_sets.at(ion_set); ion++)
			{
				orb_index = 0 ;
				for(int orb_set = 0 ; orb_set < orb_sets.size(); orb_set++)
				{
					for(int orb = 0 ; orb < orb_sets.at(orb_set); orb++)
					{
						col = orb_col * ion_set + orb_set;
						results.col(col) += dos_data.at(ion_index).col(1 + orb_index);
						//for(int en =0 ; en<NDOS; en++)
						//{
						//	results(en, 1 + col) += dos_data.at(ion_index)(en,1 + orb_index);
						//}
						orb_index ++;
					}
				}
				ion_index ++;
			}
		}






		//for (int type_id = 0; type_id < num_types; type_id++)
		//{
		//	int set_size = sets[type_id];
//
		//	for (int i = 0; i < set_size; i++)
		//	{
		//		for (int j = 0; j < NDOS; j++)
		//		{
		//			results(j, 0) = dos_data.at(ion_index)(j, 0); //energy column, should be the same for all ion; could probably just skip this step after the first ion
		//			// sum over all orbitals (index 1 to 9)
		//			for (int k = 1; k < 10; k++)
		//			{
		//				results(j, type_id) += dos_data[ion_index](j, k);
		//			}
		//		}
		//		ion_index++;
		//	}
		//}
		return results;
	}
	else return {};
}

void VASP_data::read_EIGENVAL(std::string filename)
{
	std::fstream file;
	file.open(filename, std::ios::in);
	if (file.is_open())
	{
		// skip first 5 lines
		std::string line;
		double energy;
		int occ,kpoints,NBANDS;
		for (int i = 0; i < 5; i++)
		{
			getline(file, line);
		}
		getline(file, line);
		std::stringstream ss(line);
		ss >> kpoints >> kpoints >> NBANDS; // second number is kpoints, third number is NBANDS
		BS = arma::mat(kpoints, NBANDS, arma::fill::zeros);
		occupations = arma::mat(kpoints, NBANDS, arma::fill::zeros);
		KPOINTS = arma::mat(kpoints, 4, arma::fill::zeros); // 3 for k-point coordinates and 1 for weight (if needed in the future)
		for (int k = 0; k < kpoints; k++)
		{
			getline(file, line); // skip empty line before each k-point block
			getline(file, line); // read k-point coordinates
			std::stringstream ss(line);
			ss >> KPOINTS(k,0) >> KPOINTS(k,1) >> KPOINTS(k,2) >> KPOINTS(k,3); // k-point coordinates in reciprocal space, fourth number is weight
			for (int b = 0; b < NBANDS; b++)
			{
				getline(file, line);
				std::stringstream ss2(line);

				ss2 >> occ >> energy >> occ; // first number band number, second number energy, third number occupation. We only care about energy here.
				BS(k,b) = energy;
				occupations(k, b) = occ;
			}
		}
		file.close();
	}
	else
	{
		std::cerr << "Error opening file: " << filename << std::endl;
		throw std::runtime_error("Error opening file");
	}
}

void VASP_data::write_BS(std::string filename, bool verbose_kpts, bool only_path)
{
	if (checkBS() && checkKPOINTS()) //band structure and k-points must be loaded before writing band structure data
	{
		//std::string full_filename = "workspace/" + filename + "_BS.txt";
		std::filesystem::path fullPath = std::filesystem::path("workspace") / (filename + "_BS.txt");
		std::ofstream file;

		file.open(fullPath, std::ios::out);
		file << std::fixed << std::setprecision(8);
		int kpoints = this->KPOINTS.n_rows, NBANDS = this->BS.n_cols;

		int max_k_width = std::to_string(kpoints).length();
		int max_band_width = 12;  // Enough for -XX.XXXXXXXX format (-12.34567800)
		int print_k=0;
		for (int k = 0; k < kpoints; k++)
		{
			if (only_path && (KPOINTS(k,3) != 0.0)) continue; // skip k-points that are not on the path (weight not zero)
			file <<std::setw(max_k_width)<< print_k + 1;
			if (verbose_kpts) file << " " << std::setw(10) << KPOINTS(k,0) << " " << std::setw(10) << KPOINTS(k,1) << " " << std::setw(10) << KPOINTS(k,2);
			for (int b = 0; b < NBANDS; b++)
			{
				file  << " " << std::setw(max_band_width) << BS(k,b);
			}
			file << "\n";
			print_k++;
		}
		file.close();
	}
}

void VASP_data::write_BS(std::string filename)
{
	write_BS(filename, false, false);
}

arma::mat VASP_data::get_cell_matrix()
{
	return cell_matrix;
}

arma::mat VASP_data::get_BS()
{
	return BS;
}

arma::mat VASP_data::get_occupations()
{
	return occupations;
}

arma::rowvec VASP_data::find_kpoint_energy(arma::rowvec kpt, bool weight, int& index)
{
	if (checkBS() && checkKPOINTS())
	{
		for (int i = 0; i < BS.n_rows; i++)
		{
			if (arma::approx_equal(KPOINTS.row(i).subvec(0, 2), kpt, "absdiff", 1e-6))
			{
				if (!weight && KPOINTS(i, 3) != 0.0) continue; // if weight is not considered, skip k-points that are not on the path (weight not zero)
				index = i;
				return BS.row(i);
			}
		}
		return arma::rowvec();
	}
	else return arma::rowvec();
}

double VASP_data::find_band_extremum(int band_index, bool weight, int& kpt_index, bool maxormin)
{
	if (checkBS() && checkKPOINTS())
	{
		double extremum_energy = maxormin ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
		for (int i = 0; i < BS.n_rows; i++)
		{
			if (!weight && KPOINTS(i, 3) != 0.0) continue; // if weight is not considered, skip k-points that are not on the path (weight not zero)
			if (maxormin)
			{
				if (BS(i, band_index) > extremum_energy)
				{
					extremum_energy = BS(i, band_index);
					kpt_index = i;
				}
			}
			else
			{
				if (BS(i, band_index) < extremum_energy)
				{
					extremum_energy = BS(i, band_index);
					kpt_index = i;
				}
			}

		}
		return extremum_energy;
	}
	return -1e14;
}

int VASP_data::find_valence_band()
{
	if (checkBS() && checkKPOINTS())
	{
		int valence_band_index = -1;
		for (int b = 0; b < BS.n_cols; b++)
		{

			if (occupations(0, b) > 0.5) // assuming occupation is either 0 or 1, with possible small numerical noise
			{
				valence_band_index = b;
			}
			else
			{
				break; // since bands are usually ordered by energy, we can break once we find an unoccupied state
			}

		}
		return valence_band_index;
	}
	else return -1;
}

VASP_data VASP_data::supercell_grid(int rep_x, int rep_y, int rep_z,std::vector<bool> add_vacuum)
{
	if(rep_x==0 || rep_y==0 || rep_z==0)
	{
		std::cerr << "One of dimensions is reduced to 0 !. Can't find proper transformation"<<std::endl;
		throw std::runtime_error("Improper cell multiplication");
	}
	arma::mat new_cell_matrix = cell_matrix;
	arma::vec a = cell_matrix.row(0).t();
	arma::vec b = cell_matrix.row(1).t();
	arma::vec c = cell_matrix.row(2).t();
	std::vector<arma::mat> new_types_atom_positions; // start with original positions, then we will add new positions for the supercell
	std::vector<int> new_atoms_per_type = atoms_per_type;
	std::vector<arma::mat> cart_types;
	for (int i = 0; i < types_atom_positions.size(); i++)
	{
		cart_types.push_back(cell_matrix.t() * types_atom_positions[i]); // convert to cartesian coordinates for easier manipulation
	}
	arma::mat cart_new;
	new_cell_matrix.row(0) *= rep_x + add_vacuum[0] + add_vacuum[1]; // additonal lattice vector for vacuum below and above the original cell in x direction
	new_cell_matrix.row(1) *= rep_y + add_vacuum[2] + add_vacuum[3]; // additonal lattice vector for vacuum below and above the original cell in y direction
	new_cell_matrix.row(2) *= rep_z + add_vacuum[4] + add_vacuum[5]; // additonal lattice vector for vacuum below and above the original cell in z direction
	arma::mat frac = arma::inv(new_cell_matrix);
	for ( int i=0; i< atoms_per_type.size(); i++)
	{
		new_atoms_per_type[i] *= (rep_x * rep_y * rep_z);
	}
	if(rep_x >1 || (add_vacuum[0]||add_vacuum[1])) multiply_cell_in_direction(cart_types, a, rep_x, add_vacuum[0], add_vacuum[1]);
	if(rep_y >1 || (add_vacuum[2]||add_vacuum[3])) multiply_cell_in_direction(cart_types, b, rep_y, add_vacuum[2], add_vacuum[3]);
	if(rep_z >1 || (add_vacuum[4]||add_vacuum[5])) multiply_cell_in_direction(cart_types, c, rep_z, add_vacuum[4], add_vacuum[5]);
	for (int i = 0; i < cart_types.size(); i++)
	{
		new_types_atom_positions.push_back(frac.t() * cart_types[i]); // convert back to fractional coordinates
	}
	arma::mat new_atom_positions;
	for (int i = 0; i < new_types_atom_positions.size(); i++)
	{
		new_atom_positions = arma::join_rows(new_atom_positions, new_types_atom_positions[i]);
	}
	VASP_data supercell;
	supercell.cell_matrix = new_cell_matrix;
	supercell.atoms_per_type = new_atoms_per_type;
	supercell.types_atom_positions = new_types_atom_positions;
	supercell.atom_positions = new_atom_positions;
	supercell.atom_names = atom_names;
	return supercell;

}

double VASP_data::calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R)
{
	//double epsilon_0 = 8.854187817e-12; // vacuum permittivity in F/m
	// Using atomic units where 1/(4pieps0) = 1
	double R_len = arma::norm(R);
	arma::vec r_hat = R / R_len;
	double potential = (arma::dot(dip_1, dip_2) - 3 * arma::dot(dip_1, r_hat) * arma::dot(dip_2, r_hat)) / (R_len * R_len * R_len);
	return potential;
}

arma::vec VASP_data::calc_dip_dip_force(arma::vec dip_1, arma::vec dip_2, arma::vec R) //check the derivation : -grad_r calc_dip_dip_potential
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

void VASP_data::write_DOS_sum_types(std::string id, const arma::mat& dos_summed, int atoms_sep_type,int orbitals_sep_type, bool header)
{
	std::vector<char> orb = {'s', 'p' , 'd' , 'f' , 'g'};
	std::vector<std::string> orb_spec = {
		"s", 
		"px","py","pz", 
		"dx2_xy2", "dx2_xy2", "dx2_xy2", "dx2_xy2", "dx2_xy2",
		"fxxx","fxxx","fxxx","fxxx","fxxx","fxxx","fxxx",
		"gxxx","gxxx","gxxx","gxxx","gxxx","gxxx","gxxx","gxxx","gxxx"
	}; //TO DO : check how those are called and sorted on VASP wiki

	std::vector<int> atom_sets;
	std::vector<int> orb_sets;
	int ions = dos_data.size(), NDOS = dos_data[0].n_rows, atom_types= atoms_per_type.size();
	int l_tot = dos_data[0].n_cols -1 , m = sqrt(l_tot) -1;
	std::vector<std::string> atom_names;
	std::vector<std::string> orb_names;
	bool poscar_names;
	if(atom_names.size() == 0) // if atom names are not provided in POSCAR, we will just use generic names like Atom1, Atom2, etc.
	{
		poscar_names = false;
		for(int i=0 ; i<atom_types; i++)
		{
			atom_names.push_back("Atom" + std::to_string(i+1));
		}
	}
	else
	{
		if(atom_names.size() != atom_types)
		{
			std::cerr << "Error: The number of atom names does not match the number of atom types. Please check your POSCAR file." << std::endl;
			throw std::runtime_error("Error: Mismatch between atom names and atom types");
		}
	}


	if(atoms_sep_type == 0) 
	{
		atom_sets.push_back(ions);
		atom_names.push_back("Total");
	}
	else if (atoms_sep_type == 1) 
	{
		atom_sets = atoms_per_type;
		if(poscar_names) atom_names = this->atom_names; // use atom names from POSCAR if available, otherwise use generic names
	}
	else if(atoms_sep_type == 2) 
	{
		for(int i=0 ; i<ions ; i++) atom_sets.push_back(1);
		if(poscar_names)
		{
			for(int i=0 ; i<atom_types; i++) 
			{
				for(int j=0 ; j<atoms_per_type[i]; j++)
				{
					atom_names.push_back(this->atom_names[i] + "_" + std::to_string(j+1));
				}
			}
		} 
	}

	if(orbitals_sep_type == 0) 
	{
		orb_sets.push_back(l_tot);
		orb_names.push_back("Total");
	}
	else if (orbitals_sep_type == 1) 
	{
		for(int i=0 ; i<=m ; i++) orb_sets.push_back(i*2+1);
		for(int i=0 ; i<=m ; i++) orb_names.push_back(std::string(1, orb.at(i)) + "-orbitals");

	}
	else if(orbitals_sep_type == 2) 
	{
		for(int i=0 ; i<l_tot ; i++) orb_sets.push_back(1);
		orb_names = orb_spec; // use specific orbital names for each column if orbitals are separated
	}

	std::fstream file;
	//std::string file_name = "workspace\\" + id + "_DOS_sum.txt";
	std::filesystem::path fullPath = std::filesystem::path("workspace") / (id + "_DOS_sum.txt");
	file.open(fullPath, std::ios::out);
	file.precision(12);
	int type_num = dos_summed.n_cols - 1;
	int orb_col = orb_sets.size();
	if (header)
	{
		std::string name;
		file << "# DOS summed " << "\n";
		file << "# Energy (eV)  ";
		for (int type_id = 0; type_id < type_num; type_id++)
		{
			name = atom_names.at(type_id / orb_col) + "_" + orb_names.at(type_id % orb_col);
			file << name << "  ";
		}
		file << "\n";
	}
	for (int i = 0; i < NDOS; i++)
		{
			for (int type_id = 0; type_id < type_num; type_id++)
			{
				file << dos_summed(i, type_id) << " "; 	
			}
			file << "\n";
		}
	file.close();
}
