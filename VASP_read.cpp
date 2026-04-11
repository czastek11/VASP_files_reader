#include "VASP_read.h"

std::vector<double> moving_average(std::vector<double> data, int window_size)
{
	std::vector<double> averaged(data.size(), 0.0);
	int half_window = window_size / 2; //half it as it is from both sides, this way the lenght of the window is not effectively doubled
	for (int i = 0; i < data.size(); i++) //for each point
	{
		double sum = 0.0;
		int count = 0;
		for (int j = -half_window; j <= half_window; j++) //average all the points in the window
		{
			int idx = i + j;
			if (idx >= 0 && idx < data.size()) //safety check to not go out of the array range at the start and end due to window size
			{
				sum += data[idx];
				count++;
			}
		}
		averaged[i] = sum / count;
	}
	return averaged;
}

void multiply_cell_in_direction(std::vector<arma::mat>& cart_types, arma::vec base_vec, int rep, double add_vacuum_below, double add_vacuum_above)
{
	arma::mat cart_new;
	std::vector<arma::mat> starting_cart_types = cart_types;
	std::vector<int> starting_num_types;
	for(int i =0 ;i<cart_types.size();i++) starting_num_types.push_back(cart_types.at(i).n_cols);
	for (int i = 0; i < rep; i++) //repeat all of this for each repetition 
	{
		if(add_vacuum_below > 0 && i==0) // we need to move up all ions up without adding new, if we want to create vacuum below
		{
			for(int t = 0; t < cart_types.size(); t++) //separate the atom types to not mix them and keep correct order in supercell
			{
				
				for(int pos = 0 ; pos < starting_num_types.at(t); pos++)
				{
					starting_cart_types[t].col(pos) += base_vec * add_vacuum_below; //move up all ions in base_vec direction
				}
				
			}
			cart_types = starting_cart_types; //repleace base positon instead of adding
		}
		else if(i>0) //case for multiplying cell in give direction. First repetition is avoided because it's the base cell.
		{	//In other words don't do anything if repetition is set to 1
			for(int t = 0; t < cart_types.size(); t++)
			{
				cart_new = starting_cart_types[t];
				for(int pos = 0 ; pos < starting_num_types.at(t); pos++)
				{
					cart_new.col(pos) += base_vec * i; //each repetition is n times futher in base_vec direction
				}
				cart_types[t] = arma::join_rows(cart_types[t], cart_new); // add new potions to results
			}
		}
	}
}

VASP_data::VASP_data() :
 	NGiF(), atoms_per_type(), types_atom_positions(), atom_positions(),
	charge_density_raw(), charge_density(), potential(), dos_data(), 
	KPOINTS(), BS(), atom_names(), occupations() //construct using empty constructor of all members
{
	cell_matrix = arma::mat(3, 3, arma::fill::zeros); //except cell matrix which always is going to be 3x3. So just make 0 3x3 matrix.
}

VASP_data::VASP_data(std::string file_path,bool read_POSCAR, bool read_CHGCAR, bool read_LOCPOT, bool read_DOS, bool read_EIGENVAL) : VASP_data()
{
	if (read_POSCAR) this->read_POSCAR(file_path + "POSCAR");
	if (read_CHGCAR) this->read_CHGCAR(file_path + "CHGCAR");
	if (read_LOCPOT) this->read_LOCPOT(file_path + "LOCPOT");
	if (read_DOS) this->read_DOS(file_path + "DOSCAR", false); // currently only support non-spin-orbit DOSCAR with s,p,d,f orbitals separated. Will add more options in the future if needed.
	if (read_EIGENVAL) this->read_EIGENVAL(file_path + "EIGENVAL");
	//read all selected files when constructing the object with default settings
}

VASP_data::~VASP_data() //destructor. All attributes are objects of classes with build in destructors, so we don't need to delete anyting here
{
	
}

bool VASP_data::checkgeo()
{
	if (cell_matrix.is_empty() || atoms_per_type.size() == 0 || types_atom_positions.size() == 0 || atom_positions.is_empty())
	{
		std::cerr << "Error: POSCAR data not loaded. Please load POSCAR data before using geometry data." << std::endl;
		throw std::runtime_error("Error: POSCAR data not loaded");
		return false;
		//return cerr and runtime error messages and give false
	}
	else return true;
}

bool VASP_data::checkmesh()
{
	if (NGiF.size() != 3) // ??? why geo checks NGiF which is mesh for example of CHGCAR and cell matrix ? Should make separate checkmesh
	{
		std::cerr << "Error: Either you did not load data first or did something wrong that NGiF does not have 3 elements." << std::endl;
		throw std::runtime_error("Error: NGiF does not have 3 elements");
		return false;
		//return cerr and runtime error messages and give false
	}
	else return true;
}

bool VASP_data::checkcharge()
{
	if (charge_density_raw.is_empty()) //by default this attribute is empty so after loading the data it will satisfy this condition
	{
		std::cerr << "Error: charge density data not loaded. Please load data before using charge density." << std::endl;
		throw std::runtime_error("Error: charge density data not loaded");
		return false;
		//same as other checkxxx()
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
//same as before

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
//saem as before

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
//same as before

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
//same as before

arma::mat VASP_data::sorting_positions(arma::mat positions, std::string method)
{
	if (method == "z_rising")// one method implemented so far - sort by third lattice vector direction
	{
		arma::uvec sorted_indices = arma::sort_index(positions.row(2).t(), "ascend"); //arma sorts the rows by the specified column and mode
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
		// get the cell lattise base vectors
		std::vector<int> indices(3);
		arma::mat transform_matrix = arma::mat(3, 3);
		transform_matrix.col(0) = a / NGiF.at(0);
		transform_matrix.col(1) = b / NGiF.at(1);
		transform_matrix.col(2) = c / NGiF.at(2);
		// construct the matrix of transformation R = x_1 * a / n1 + x_2 * b / n2 + x_3 * c / n3
		//transform_matrix.print("Transform matrix:");
		arma::vec fractional_coords = arma::inv(transform_matrix) * pos; //perform the transformation into mesh postions
		for (int i = 0; i < 3; i++)
		{
			indices[i] = static_cast<int>(floor(fractional_coords(i) + 0.5)); //round to nearest indice - intiger
			if (indices[i] < 0) indices[i] = 0; //all negatives results are set as 0
			if (indices[i] >= NGiF[i]) indices[i] = NGiF[i] - 1; //all indices above the cell are set to highest indice
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
	//clear previous data only after successfully opening the file, to avoid partial clearing if file opening fails
	atoms_per_type.clear();
	types_atom_positions.clear();
	atom_positions.clear();
	//clear all the data to be sure, a lot of them are vectors or arma::mat added by joing col/row 
	//so we need to avoid adding to the old data - a lot of bugs were caused by this in repeated usage of VASP_data object
	
	getline(file, line); // Skip the first line (comment)
	getline(file, line);

	for (int i = 0; i < 9; i++) //read up t
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
	while (ss >> int_read) // readout the number of ions for each atom
	{
		atoms_per_type.push_back(int_read);
		num_atom += int_read;
	}

	getline(file, line); //skip direct

	for (int i = 0; i < atoms_per_type.size(); i++)
	{
		types_atom_positions.push_back(arma::mat(3, atoms_per_type[i])); //prepare the matrix to read positions of each type of atoms, we will sort them later by the method we want, but for now we just read them in the order they are in the file
		for (int j = 0; j < atoms_per_type[i]; j++)
		{
			getline(file, line);
			std::stringstream ss2(line);
			double x, y, z;
			ss2 >> x >> y >> z;
			arma::vec pos = { x, y, z };
			types_atom_positions[i].col(j) = pos;
			// read out the position of each ion
			// and store them in the matrix of the corresponding atom type, each column is one ion
		}
		types_atom_positions[i] = sorting_positions(types_atom_positions[i], "z_rising"); //sort the atoms of the same type by their z coordinate, this is useful for layered materials like TMDS, but can be changed to other methods in the future if needed
		// TO DO: add option to choose other mathod after we implement them
		atom_positions = arma::join_rows(atom_positions, types_atom_positions[i]);
	}
	
}

void VASP_data::read_POSCAR(std::string filename)
{
	std::fstream file;
	read_POSCAR_like(filename,file); //the POSCAR header is the same as POSCCAR file itself, so just run this function
	if(file.is_open()) file.close(); // and then close the file after
	else 
	{
		std::cout << "POSCAR not probably not read!\n";
		std::cerr << "POSCAR not probably not read!\n";
	}
}

void VASP_data::read_bestsqs(std::string filename)
{
	std::fstream file;
	std::string line,word;
	double number_read;
	int int_read, num_atom = 0;

	file.open(filename, std::ios::in);
	if (!file.is_open())
	{
		std::cerr << "Error opening file: " << filename << std::endl;
		throw std::runtime_error("Error opening file");
	}
	//clear previous data only after successfully opening the file, to avoid partial clearing if file opening fails
	atoms_per_type.clear();
	types_atom_positions.clear();
	atom_positions.clear();
	atom_names.clear();
	//clear all the data to be sure, a lot of them are vectors or arma::mat added by joing col/row 
	//so we need to avoid adding to the old data - a lot of bugs were caused by this in repeated usage of VASP_data object

	arma::mat old_cell_matrix = arma::mat(3, 3, arma::fill::zeros);
	arma::mat trans_matrix = arma::mat(3, 3, arma::fill::zeros);
	
	for (int i = 0; i < 9; i++) //read up base cell matrix
	{
		file >> number_read;
		old_cell_matrix(i / 3, i % 3) = number_read;
	}
	for (int i = 0; i < 9; i++) //read up trans matrix
	{
		file >> number_read;
		trans_matrix(i / 3, i % 3) = number_read;
	}
	arma::mat inv = arma::inv(trans_matrix);
	arma::vec pos = arma::vec(3,arma::fill::zeros);
	std::map <std::string, int> id_name;
	std::vector<std::vector<arma::vec>> list;
	list.clear();
	int ii = 0, point=0;
	getline(file, line);
	while (getline(file, line))
	{
		std::stringstream ss(line);
		ss >> pos(0) >> pos(1) >> pos(2); //read up coordiantes
		ss >> word; //atom name 
		if (std::find(atom_names.begin(), atom_names.end(), word) != atom_names.end()) //check if already encountered this name
		{
			point = id_name[word];
			atoms_per_type.at(point)++;
			list.at(point).push_back(pos);
		}
		else //create new set for given name if one does not exist
		{
			atom_names.push_back(word);
			id_name[word] = ii;
			atoms_per_type.push_back(1);
			std::vector<arma::vec> dummy;
			list.push_back(dummy);
			list.at(ii).push_back(pos);
			ii++;
		}
	}
	//now we need to insert this data into class members
	cell_matrix = trans_matrix * old_cell_matrix;
	for (int i = 0; i < atoms_per_type.size(); i++)
	{
		int n = atoms_per_type.at(i);
		arma::mat loc_positions = arma::mat(3, n, arma::fill::zeros);
		for (int j = 0; j < n; j++) loc_positions.col(j) = list.at(i).at(j);
		loc_positions = inv * loc_positions;
		loc_positions = sorting_positions(loc_positions, "z_rising");
		types_atom_positions.push_back(loc_positions);
		atom_positions = arma::join_rows(atom_positions, loc_positions);
	}
	file.close();
}

void VASP_data::read_CHGCAR(std::string filename)
{
	std::fstream file;
	std::string line;
	int int_read;
	read_POSCAR_like(filename, file); //read POSCAR like header
	NGiF.clear();
	getline(file, line); // Skip the line before charge density data
	getline(file, line);

	std::stringstream ss2(line);
	for (int i = 0; i < 3; i++)
	{
		ss2 >> int_read;
		NGiF.push_back(int_read);
	}

	int total_grid_points = NGiF[0] * NGiF[1] * NGiF[2];

	
	charge_density_raw = arma::cube(NGiF[0], NGiF[1], NGiF[2], arma::fill::zeros);
	charge_density = arma::cube(NGiF[0], NGiF[1], NGiF[2], arma::fill::zeros);
	//by defining new cube (3d array) old data is cleared automaticly 

	for (int i = 0; i < total_grid_points; i++) //reading out strem of values written by VASP using fortran. x is fastest changing (inside loop)
	{
		int ix = i % NGiF[0]; // index x is least significant digit in NGiF base: i = n2 * n1 * iz + n1 * iy + ix ==> i mod n1 = ix
		int iy = (i / NGiF[0]) % NGiF[1]; // index y is middle significant digit in NGiF base: i = n2 * n1 * iz + n1 * iy + ix ==> i / n1 mod n2 = iy, ix<n1
		int iz = i / (NGiF[0] * NGiF[1]); // index z is most significant digit in NGiF base: i = n2 * n1 * iz + n1 * iy + ix ==> i / n2 /n1  = iz, ix<n1 * n2, iy*n1 < n1 * n2
		file >> charge_density_raw(ix, iy, iz); //write in correct mesh postion read data value
		//charge_density[ix][iy][iz] = charge_density_raw[ix][iy][iz] / total_grid_points / cell_volume; //wrong normalization
		charge_density(ix, iy, iz) = charge_density_raw(ix, iy, iz) / total_grid_points; // normalisation. Sum of all those values gives number of electrons
	}
	file.close();
}



void VASP_data::read_LOCPOT(std::string filename)
{
	std::fstream file;
	std::string line;
	int int_read;
	read_POSCAR_like(filename, file); //read POSCAR like header
	NGiF.clear();
	getline(file, line); // Skip the line before charge density data
	getline(file, line);

	std::stringstream ss2(line);
	for (int i = 0; i < 3; i++)
	{
		ss2 >> int_read;
		NGiF.push_back(int_read);
	}

	int total_grid_points = NGiF[0] * NGiF[1] * NGiF[2];
	int ix = 0, iy = 0, iz = 0;
	
	potential = arma::cube(NGiF[0], NGiF[1], NGiF[2], arma::fill::zeros);
	//by defining new cube (3d array) old data is cleared automaticly 


	for (int i = 0; i < total_grid_points; i++)
	{
		ix = i % NGiF[0]; // index x is least significant digit in NGiF base: i = n2 * n1 * iz + n1 * iy + ix ==> i mod n1 = ix
		iy = (i / NGiF[0]) % NGiF[1]; // index y is middle significant digit in NGiF base: i = n2 * n1 * iz + n1 * iy + ix ==> i / n1 mod n2 = iy, ix<n1
		iz = i / (NGiF[0] * NGiF[1]); // index z is most significant digit in NGiF base: i = n2 * n1 * iz + n1 * iy + ix ==> i / n2 /n1  = iz, ix<n1 * n2, iy*n1 < n1 * n2
		file >> potential(ix,iy,iz); //write in correct mesh postion read data value
		//potentai lcan have up to 4 datasets for spin polarized calculations, but we only read the first one here
		//current example has 4 but magmom is set to 0 for all atoms, so the other 3 datasets are mostly zeros
	}
	file.close();
}

void VASP_data::write_POSCAR(std::string filename)
{
	std::filesystem::path fullPath = std::filesystem::path("workspace") / (filename + "_POSCAR"); //used filesystem::path to handle windows/unix path formatting automaticly
	std::fstream file;
	file.open(fullPath, std::ios::out);
	if (!file.is_open()) //error handling
	{
		std::cerr << "Error opening file for writing: " << fullPath << std::endl;
		throw std::runtime_error("Error opening file for writing");
	}
	for(int i=0 ; i< atoms_per_type.size(); i++) //read out atoms names and their quantity into header
	{
		file << atom_names[i] << atoms_per_type[i] << " "; // M2 S4 ...
	}
	file<<std::fixed << std::setprecision(16); //write all entries in VASP precisions
	file<< "\n";
	file << "   "<<1.0<<"\n"; //scaling factor. To this day I've always got 1.0 and never put any other myself.
	for (int i = 0; i < 3; i++) // write out cell matrix in three lines with columns separated by two spaces
	{
		file<<" ";
		for (int j = 0; j < 3; j++)
		{
			if(cell_matrix(i,j)>=0) file<<" "; // leave space for minus so that positive and negative entries share the same length
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
		for (int j = 0; j < types_atom_positions[i].n_cols; j++) //write out ion postion in each row
		{
			for (int k = 0; k < 3; k++)
			{
				if(types_atom_positions[i](k,j)>=0) file<<" "; // space for negative sign
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
		// sanity check
		double nel = 0.0;
		for (int i = 0; i < NGiF[0]; i++) //just summing all the values in charge density
		{
			for (int j = 0; j < NGiF[1]; j++)
			{
				for (int k = 0; k < NGiF[2]; k++)
				{
					nel += charge_density_raw(i,j,k);
				}
			}
		}
		nel /= (NGiF[0] * NGiF[1] * NGiF[2]); //normalisation
		return nel; 
	}
	else return 0.0;
}

int VASP_data::count_total_electrons()
{
	return static_cast<int>(count_total_electrons_double() + 0.5); // rounding to nearest integer
}

std::vector<double> VASP_data::sum_potential_averaged_xy_z(std::string period_type, int period)
{
	if (checkpot())
	{
		std::vector<double> potential_z(NGiF[2], 0.0);
		int total_xy_points = NGiF[0] * NGiF[1];
		double sum;
		// for each point in third direction
		for (int k = 0; k < NGiF[2]; k++)  // z index 
		{ 
			sum = 0.0; //sum all over xy plane for given z index
			for (int i = 0; i < NGiF[0]; i++)// x index
			{ 
				for (int j = 0; j < NGiF[1]; j++) // y index
				{ 
					sum += potential(i,j,k);
				}
			}
			potential_z[k] = sum / total_xy_points; // and then average it
		}

		// average potential in z over moving window. Window size is to be set
		int window = 1;
		if (period_type=="primitive") //for normal bulk cell, averaging over whole length
		{
			window = NGiF[2];
		}
		else if (period_type == "manual") // providing windows size in number of mesh points in z direction manually
		{
			window = period;
		}
		else if (period_type == "layered") // for many layers (at least 3). It checks the difference betwen first and third ion of the first type
		{
			arma::vec distance = cell_matrix.t() * (atom_positions.col(2) - atom_positions.col(0));
			window = get_mesh_indices(distance)[2];
		}
		//this method is specific for TMDS and should be later generalized for at least other directions and number of layers in bulk

		if (period_type!="none") potential_z = moving_average(potential_z, window); //don't perform moving average if period type is "none"

		return potential_z;
	}
	else return {};
}

std::vector<double> VASP_data::sum_potential_averaged_xy_z(std::string period_type) //overload for when user doesn't provide period value, but sets period type to something other than "manual". If "manual" is set without providing period value, it raises an error
{
	if(period_type!="manual") return sum_potential_averaged_xy_z(period_type, 0.0);
	else
	{
		std::cerr << "Error: For manual period type, you must provide a period value." << std::endl;
		throw std::runtime_error("Error: Missing period value for manual period type");
		return {};
	}
}


std::vector<double> VASP_data::average_potential_over(int direction)
{
	if (checkpot())
	{
		std::vector<int> dirs = { 0,1,2 };
		int dir = dirs.at(direction - 1);
		dirs.erase(dirs.begin() + dir);
		std::vector<double> potential_av(NGiF[dir], 0.0);
		int total_par_point = NGiF[dirs.at(0)] * NGiF[dirs.at(1)];
		double sum;
		// Create an array to map loop indices to potential arguments
		int indices[3];  // Will hold {x, y, z} indices
		// for each point in third direction
		for (int k = 0; k < NGiF[dir]; k++)  // perpendicular index 
		{
			sum = 0.0; //sum all over plane for given perpendicular index

			indices[dir] = k;
			for (int i = 0; i < NGiF[dirs.at(0)]; i++)// first pararell index
			{
				indices[dirs.at(0)] = i;
				for (int j = 0; j < NGiF[dirs.at(1)]; j++) // second pararell index
				{
					indices[dirs.at(1)] = j;
					sum += potential(indices[0], indices[1], indices[2]);
				}
			}
			potential_av[k] = sum / total_par_point; // and then average it
		}
		return potential_av;
	}
}

std::vector<double> VASP_data::moving_average_potential_over(std::vector<double> av_pot, int direction, std::string period_type)
{
	// average potential in z over moving window. Window size is to be set
	int window = 1;
	std::vector<int> dirs = { 0,1,2 };
	int dir = dirs.at(direction - 1);
	dirs.erase(dirs.begin() + dir);
	std::vector<double> pot_mov_av;
	if (period_type == "primitive") //for normal bulk cell, averaging over whole length
	{
		window = NGiF[dir];
	}
	else
	{
		std::cerr << "Wrong avering type" << std::endl;
		throw std::runtime_error("No matching type in moving_average_potential_over");
	}
	//this method is specific for TMDS and should be later generalized for at least other directions and number of layers in bulk
	pot_mov_av = moving_average(av_pot, window); //don't perform moving average if period type is "none"
	return pot_mov_av;
}

std::vector<double> VASP_data::moving_average_potential_over(std::vector<double> av_pot, int direction, std::string period_type, int period)
{
	// average potential in z over moving window. Window size is to be set
	int window = 1;
	std::vector<int> dirs = { 0,1,2 };
	int dir = dirs.at(direction - 1);
	dirs.erase(dirs.begin() + dir);
	std::vector<double> pot_mov_av;
	if (period_type == "manual") // providing windows size in number of mesh points in z direction manually
	{
		window = period;
	}
	else
	{
		std::cerr << "Wrong avering type" << std::endl;
		throw std::runtime_error("No matching type in moving_average_potential_over");
	}
	pot_mov_av = moving_average(av_pot, window);
	return pot_mov_av;
}

std::vector<double> VASP_data::moving_average_potential_over(std::vector<double> av_pot, int direction, std::string period_type, int ion1 , int ion2)
{
	// average potential in z over moving window. Window size is to be set
	int window = 1;
	std::vector<int> dirs = { 0,1,2 };
	int dir = dirs.at(direction - 1);
	dirs.erase(dirs.begin() + dir);
	std::vector<double> pot_mov_av;
	if (period_type == "layered") // for many layers (at least 3). It checks the difference betwen first and third ion of the first type
	{
		// Get the real-space vector between atoms
		arma::vec distance_vec = cell_matrix.t() * (atom_positions.col(ion2) - atom_positions.col(ion1));

		// Get the lattice basis vector for the direction we care about
		arma::rowvec basis_vec = cell_matrix.row(dir);

		// Project distance onto this basis direction
		// This gives how many times the basis vector fits into the distance component
		double projection = arma::dot(distance_vec, basis_vec) / arma::dot(basis_vec, basis_vec);

		// Now convert to mesh indices
		// Create a vector that only has component in this direction
		arma::vec dir_only_vec = projection * basis_vec.t();  // or projection * unit_vec.t()

		// Get mesh indices for this direction-only vector
		window = get_mesh_indices(dir_only_vec)[dir];

		// Make sure window is at least 1
		window = std::max(1, window);
	}
	else
	{
		std::cerr << "Wrong avering type" << std::endl;
		throw std::runtime_error("No matching type in moving_average_potential_over");
	}
	pot_mov_av = moving_average(av_pot, window); 
	return pot_mov_av;
}

void VASP_data::write_potential_z(std::string filename, std::vector<double> potential_z)
{
	std::filesystem::path fullPath = std::filesystem::path("workspace") / (filename + "_potential_z.txt"); //save in workspace folder with filename + _potential_z.txt
	std::fstream file;
	file.open(fullPath, std::ios::out);
	double z_real;
	for (int i = 0; i < NGiF[2]; i++)
	{
		z_real = i * cell_matrix(2, 2) / NGiF[2]; //transform from mesh index to real space coordinate in direction
		//since for know TMDS have 0 0 1 vector this is just z coordinate,
		// but later it should be generalised to the legth of the third lattice vector direction
		// then it would go from 0 to |c| and not necessary from 0 to cell_matrix(2,2) if the third lattice vector is not along z direction
		file << z_real << " " << potential_z[i] << "\n";
	}
	file.close();
}

void VASP_data::write_potential_over(std::string filename, std::vector<double> potential_av, int direction)
{
	std::vector<int> dirs = { 0,1,2 };
	int dir = dirs.at(direction - 1);
	dirs.erase(dirs.begin() + dir);
	std::filesystem::path fullPath = std::filesystem::path("workspace") / (filename + "_potential_" + std::to_string(direction) + ".txt"); //save in workspace folder with filename + _potential_z.txt
	std::fstream file;
	file.open(fullPath, std::ios::out);
	double dir_real;
	arma::mat base = cell_matrix.t();
	arma::vec base_vec = base.col(dir);
	double norm = arma::norm(base_vec);
	for (int i = 0; i < NGiF[dir]; i++)
	{
		dir_real = i * norm / NGiF[dir];
		file << dir_real << " " << potential_av[i] << "\n";
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
		//lattice base vectors
		for (int i = start[0]; i < end[0]; i++) //loop over mesh points in the specified range
		{
			for (int j = start[1]; j < end[1]; j++)
			{
				for (int k = start[2]; k < end[2]; k++)
				{
					pos = i * a / NGiF[0] + j * b / NGiF[1] + k * c / NGiF[2];  //transform from mesh indices to real space coordinates
					dipole_mom += -charge_density(i,j,k) * (pos - center); //calculate contribution to dipole moment from this mesh point and add it to the total dipole moment. Charge density is negative for electrons, so we add minus sign to get correct direction of dipole moment vector
					// the contribution is negetaive since charge density describes electrons, which have negative charge
					// this contribution will be postive for ions, but those will summed outside, cause it's easy
					// in future we can add option to include ions contribution to dipole moment as well,
					// since we just need to sum ion position from POSCAR
					// but we need to know the amount of valence electrons for each ion type,
					// this is not obvious how since there are different pseudopotentials with different number of valence electrons for the same element.
					// either user needs to provide this or OUTCAR/POTCAR can be read?
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
		arma::vec pos;
		arma::vec a = cell_matrix.row(0).t();
		arma::vec b = cell_matrix.row(1).t();
		arma::vec c = cell_matrix.row(2).t();
		//base vectors and position vector to write out potential in real space coordinates
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
					pos = i * 1.0 / NGiF[0] * a + j * 1.0 / NGiF[1] * b + k * 1.0 / NGiF[2] * c; //transform to cartesian coordinates
					file << pos(0) << " " << pos(1) << " " << pos(2) << " " << potential(i,j,k) << "\n";
					//write out x,y,z coordinates and potential value in each line, separated by space
				}
			}
		}
		file.close();
	}
}

void VASP_data::read_DOS(std::string filename, bool spin_orbit)
{
	// only no spin orbit interaction case 
	// for read it should be ieither bool or maybe can be read from DOSCAR itself
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
			int l_tot = 0;
			int NDOS;
			double pom;
			int ions;
			getline(file, line);
			std::stringstream ss4(line);
			ss4 >> ions >> ions; // to get number of ions, which is the second number in the first line of DOSCAR, first number includes empty spheres
			for (int i = 0; i < 4; i++)// skip next 4 lines
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
			// skip first set with two columns (total DOS)
			// since total DOS is written sa two columns there is no need to handle this by program,
			// besides the header it is already sorted and ready for plotting
			for (int i = 0; i < NDOS; i++) getline(file, line);

			// read DOS data
			for (int ion = 0; ion < ions; ion++)
			{
				
				getline(file, line); // skip header line for each ion
				getline(file, line); // first line
				std::stringstream ss2(line);
				if (ion == 0) // check max l value from first line once, initialise the matrix as well, since all ions have the same output we don't need to clear the matrix
				{
					
					while (ss2 >> pom) l_tot++;
					l_tot--; // first column is energy, so total number of orbitals is total columns - 1
					ion_dos = arma::mat(NDOS, l_tot + 1, arma::fill::zeros); // energy + l_tot orbitals 1s,3p,5d,7f,...
					ss2.clear();
					ss2.seekg(0); // reset stringstream to read the line again for the first ion after determining l_tot
				}
				for (int j = 0; j < l_tot + 1; j++) //first line
				{
					ss2 >> ion_dos(0, j);
				}
				for (int i = 1; i < NDOS; i++) //rest
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
	// for spin there will be another option:
	// -1 no spin, 0 sum the spins, 1 keep the spins separetly
	if (checkdos())
	{
		int ions = dos_data.size(), NDOS = dos_data[0].n_rows, atom_types= atoms_per_type.size();
		int l_tot = dos_data[0].n_cols -1 , m = sqrt(l_tot)-1; // biggest m number in the set
		arma::mat results;
		std::vector<int> atom_sets;
		std::vector<int> orb_sets;

		// setting up the sets for ions and orbitals
		// 0 is all in one group
		// 1 is groups of same atom/same orbital type
		// 2 is all separate
		if(atoms_sep_type == 0) atom_sets.push_back(ions);
		else if (atoms_sep_type == 1) 
		{
			if (checkgeo()) //becuse we need to know the atom groups for this separation, having POSCAR data is required
			{
				atom_sets = atoms_per_type;
			}
			else
			{
				std::cerr << "Error: Cannot separate atoms by type because geometric data not loaded. Please load POSCAR or CHGCAR data before using this option." << std::endl;
				throw std::runtime_error("Error: Geometric data not loaded for atom separation by type");
			}
		}
		else if(atoms_sep_type == 2) for(int i=0 ; i<ions ; i++) atom_sets.push_back(1);

		if(orbitals_sep_type == 0) orb_sets.push_back(l_tot);
		else if (orbitals_sep_type == 1) for(int i=0 ; i<=m ; i++) orb_sets.push_back(i*2+1);
		else if(orbitals_sep_type == 2) for(int i=0 ; i<l_tot ; i++) orb_sets.push_back(1);

		int at_col = atom_sets.size(), orb_col= orb_sets.size();
		int tot = at_col * orb_col;
		results = arma::mat(NDOS, 1 + tot, arma::fill::zeros); // first column for energy, next columns for all sets combinations
		results.col(0) = dos_data[0].col(0); // fill in energy already
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
						col = orb_col * ion_set + orb_set; //there are orb_col orbitals type for each ion
						// so next ions is numerating not from 0 but from number of ions times orb_col
						results.col(1  + col) += dos_data.at(ion_index).col(1 + orb_index); // we can add columns
						// because energy and it's mesh is the same for all
						orb_index ++;
					}
				}
				ion_index ++;
			}
		}
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

				ss2 >> occ >> energy >> occ; // first number band number, second number energy, third number occupation.
				BS(k,b) = energy;
				occupations(k,b) = occ;
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

void VASP_data::read_BS(std::string filename, bool header, bool verbose_kpts)
{
	std::fstream file;
	file.open(filename, std::ios::in);
	if (file.is_open())
	{
		
		std::string line;
		double energy;
		int occ, kpoints, NBANDS;
		std::vector<std::vector<double>> read_BS;
		std::vector<arma::vec> read_kpoints;
		arma::vec kpoint = arma::vec(3, arma::fill::zeros);
		if (header)
		{
			// skip first 5 lines
			for (int i = 0; i < 4; i++)
			{
				getline(file, line);
			}
		}
		while (getline(file, line))
		{
			std::stringstream ss(line);
			
			if (verbose_kpts)
			{
				ss >> kpoints >> kpoint(0) >> kpoint(1) >> kpoint(2);
				read_kpoints.push_back(kpoint);
			}
			else ss >> kpoints;
			std::vector<double> dummy;
			while (ss >> energy)
			{
				dummy.push_back(energy);
			}
			read_BS.push_back(dummy);
		}
		file.close();
		NBANDS = read_BS.at(0).size();
		KPOINTS = arma::mat(kpoints, 4, arma::fill::zeros);
		BS = arma::mat(kpoints, NBANDS,arma::fill::zeros);
		for (int i = 0; i < kpoints; i++)
		{
			if(verbose_kpts)
			{
				KPOINTS(i, 0) = read_kpoints.at(i)(0);
				KPOINTS(i, 1) = read_kpoints.at(i)(1);
				KPOINTS(i, 2) = read_kpoints.at(i)(2);
			}
			for (int j = 0; j < NBANDS; j++)
			{
				BS(i, j) = read_BS.at(i).at(j);
			}
		}
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
			//write out kpoints themselves if verbose
			for (int b = 0; b < NBANDS; b++) // write out energies themselves
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
	write_BS(filename, false, false); //default settings n verbose kpoints and write out all kpoints
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
			if (arma::approx_equal(KPOINTS.row(i).subvec(0, 2), kpt, "absdiff", 1e-6)) //check if this row is match with provide kpoint to the tolerance
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
		// we set starting value for comparison as -infity for maximum and plus infinity for mininmum
		// could also set the first values as starting but copilot already suggested this nice option
		for (int i = 0; i < BS.n_rows; i++) //go trough all kpoints
		{
			if (!weight && KPOINTS(i, 3) != 0.0) continue; // if weight is not considered, skip k-points that are not on the path (weight not zero)
			if (maxormin) //strandar comparison logic tree to find max or min value
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

VASP_data VASP_data::supercell_grid(int rep_x, int rep_y, int rep_z,std::vector<double> add_vacuum)
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
	//base lattice vectors
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
	for ( int i=0; i< atoms_per_type.size(); i++) // multiply number of each atom by the volume of grid
	{
		new_atoms_per_type[i] *= (rep_x * rep_y * rep_z);
	}

	if(rep_x >1 || (add_vacuum[0]||add_vacuum[1])) multiply_cell_in_direction(cart_types, a, rep_x, add_vacuum[0], add_vacuum[1]);
	if(rep_y >1 || (add_vacuum[2]||add_vacuum[3])) multiply_cell_in_direction(cart_types, b, rep_y, add_vacuum[2], add_vacuum[3]);
	if(rep_z >1 || (add_vacuum[4]||add_vacuum[5])) multiply_cell_in_direction(cart_types, c, rep_z, add_vacuum[4], add_vacuum[5]);
	// multiply in each direction if there is muliplication or vacuum from either side
	// it's important that the result of multiplication is fed again for the next direction
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
	// return the supercell as new VASP_data object with initialised POSCAR informations

}


void VASP_data::alloy_geometry(VASP_data data2, double percentage, std::vector<std::string> mixed_atom_names1, std::vector<std::string> mixed_atom_names2, std::string filename)
{
	if(checkgeo() && data2.checkgeo())
	{
		arma::mat new_cell_matrix = cell_matrix * percentage + data2.cell_matrix * (1-percentage);
		arma::vec new_atom_positions;
		std::vector<std::vector<std::string>> atom_groups;
		std::string name,name2;
		bool present;
		std::vector<std::string>::iterator position1 , position2 ;
		int index, id1,id2;
		for(int i=0; i<atom_names.size(); i++)
		{
			std::vector<std::string> group;
			name = atom_names.at(i);
			if(std::find(data2.atom_names.begin(), data2.atom_names.end(), atom_names.at(i)) != data2.atom_names.end()) // if this atom type is present in both structures, we can mix them
			{
				group.push_back(name);
			}
			else 
			{
				position1 = std::find(mixed_atom_names1.begin(), mixed_atom_names1.end(), name);
				position2 = std::find(mixed_atom_names2.begin(), mixed_atom_names2.end(), name);
				if(position1 != mixed_atom_names1.end())
				{
					index = std::distance(mixed_atom_names1.begin(), position1);
					group.push_back(mixed_atom_names1.at(index));
					group.push_back(mixed_atom_names2.at(index));
				}
				else if(position2 != mixed_atom_names2.end())
				{
					index = std::distance(mixed_atom_names2.begin(), position2);
					group.push_back(mixed_atom_names1.at(index));
					group.push_back(mixed_atom_names2.at(index));
				}
				else
				{
					std::cerr << "Error: Atom type " << name << " is not present in both structures and not listed in mixed_atom_names1 or mixed_atom_names2. Cannot determine how to mix this atom type." << std::endl;
					throw std::runtime_error("Error: Atom type not found for mixing");
				}
			}
			atom_groups.push_back(group);
			// we create groups of atom types that should be mixed together
		}




		double a = arma::norm(new_cell_matrix.row(0));
		double b = arma::norm(new_cell_matrix.row(1));
		double c = arma::norm(new_cell_matrix.row(2));
		double alpha = std::acos(arma::dot(new_cell_matrix.row(1), new_cell_matrix.row(2)) / (b * c)) * 180.0 / M_PI;
		double beta = std::acos(arma::dot(new_cell_matrix.row(0), new_cell_matrix.row(2)) / (a * c)) * 180.0 / M_PI;
		double gamma = std::acos(arma::dot(new_cell_matrix.row(0), new_cell_matrix.row(1)) / (a * b)) * 180.0 / M_PI;
		std::filesystem::path fullPath = std::filesystem::path("workspace") / (filename + "_lat.in");
		std::ofstream file;
		file.open(fullPath, std::ios::out);
		file <<std::fixed << std::setprecision(16)<< a << " " << b << " " << c << " " <<std::setprecision(2)<< alpha << " " << beta << " " << gamma << "\n";
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				if (i==j) file << 1 << " ";
				else file << 0 << " ";
			}
			file << "\n";
		}
		for(int i =0 ; i<atom_groups.size(); i++)
		{
			//std::fixed << std::setprecision(16)
			//find index
			if(atom_groups.at(i).size() == 1) // if this group has only one atom type, we mix only the positions
			{
				name = atom_groups.at(i).at(0);
				id1 = std::distance(atom_names.begin(), std::find(atom_names.begin(), atom_names.end(), name));
				id2 = std::distance(data2.atom_names.begin(), std::find(data2.atom_names.begin(), data2.atom_names.end(), name));
				for(int j=0; j<atoms_per_type.at(id1); j++)
				{
					new_atom_positions = cell_matrix.t() *types_atom_positions.at(id1).col(j) *  percentage + data2.cell_matrix.t() * types_atom_positions.at(id2).col(j) * (1-percentage); // mix the positions of the this atom type
					arma::mat inverted = arma::inv(new_cell_matrix.t());
					new_atom_positions = inverted * new_atom_positions; // convert back to fractional coordinates
					file <<std::fixed << std::setprecision(16)<< new_atom_positions[0]<< " " << new_atom_positions[1] << " " << new_atom_positions[2] ;
					file << " " << name << "="<<std::setprecision(6)<< 1.0000 << "\n"; // write out the mixed number of atoms for this type
				}
			}
			else if (atom_groups.at(i).size() == 2) // if this group has two atom types, we need to mix them
			{
				name = atom_groups.at(i).at(0);
				id1 = std::distance(atom_names.begin(), std::find(atom_names.begin(), atom_names.end(), name));
				name2 = atom_groups.at(i).at(1);
				id2 = std::distance(data2.atom_names.begin(), std::find(data2.atom_names.begin(), data2.atom_names.end(), name2));
				for(int j=0; j<atoms_per_type.at(id1); j++)
				{
					new_atom_positions = cell_matrix.t() *types_atom_positions.at(id1).col(j) *  percentage + data2.cell_matrix.t() * types_atom_positions.at(id2).col(j) * (1-percentage); // mix the positions of the this atom type
					arma::mat inverted = arma::inv(new_cell_matrix.t());
					new_atom_positions = inverted * new_atom_positions; // convert back to fractional coordinates
					file <<std::fixed << std::setprecision(16)<< new_atom_positions[0]<< " " << new_atom_positions[1] << " " << new_atom_positions[2];
					file << " " << name << "="<<std::setprecision(6)<< percentage << " , " << name2 << "=" <<std::setprecision(6)<< (1-percentage) << "\n"; // write out the mixed number of atoms
				}
			}
		}
	}
}




double VASP_data::calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R)
{
	//double epsilon_0 = 8.854187817e-12; // vacuum permittivity in F/m
	// Using atomic units where 1/(4pieps0) = 1
	double R_len = arma::norm(R);
	arma::vec r_hat = R / R_len;
	double potential = (arma::dot(dip_1, dip_2) - 3 * arma::dot(dip_1, r_hat) * arma::dot(dip_2, r_hat)) / (R_len * R_len * R_len);
	return potential;
	// formula for dipole potential (d1 . d2 - 3 * d1 . r^ * d2 . r^) / |r|^3
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
	//this derived -grad of potential from function above
}

void VASP_data::write_DOS_sum_types(std::string id, const arma::mat& dos_summed, int atoms_sep_type,int orbitals_sep_type, bool header)
{
	std::vector<char> orb = {'s', 'p' , 'd' , 'f' , 'g'};
	std::vector<std::string> orb_spec = {
		"s", 
		"py","pz","px", 
		"dxy", "dyz", "dz2_r2", "dx2_xy2", "dx2_xy2",
		"fy(3x2-y2)","fxyz","fyz2","fz3","fxz2","z(x2-y2)","fx(x2-3y2)",
		"gy2(5x2-y2)","gxxxx","gxxxx","gxxxx","gz4","gx4-y4","gxxxx","gxxxx","gxxxx"
	}; //g is only for very heavy elements, and the specific order of g orbitals can vary, so we just use generic names here. The order of orbitals in DOSCAR is determined by VASP and is usually s, p, d, f, g in order of increasing l, with m values ordered as described in VASP manual. We assume this order when assigning names to orbitals.
	// g names are not verified

	std::vector<int> atom_sets;
	std::vector<int> orb_sets;
	int ions = dos_data.size(), NDOS = dos_data[0].n_rows, atom_types= atoms_per_type.size();
	int l_tot = dos_data[0].n_cols -1 , m = sqrt(l_tot) -1;
	std::vector<std::string> ions_names;
	std::vector<std::string> orb_names;
	bool poscar_names = true;
	if(atom_names.size() == 0) // if atom names are not provided in POSCAR, we will just use generic names like Atom1, Atom2, etc.
	{
		poscar_names = false;
		for(int i=0 ; i<ions; i++)
		{
			ions_names.push_back("Atom" + std::to_string(i+1));
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
		ions_names.push_back("Total");
	}
	else if (atoms_sep_type == 1) 
	{
		if (checkgeo())
		{
			atom_sets = atoms_per_type;
			ions_names = atom_names; // use atom names from POSCAR if available, otherwise use generic names
		}
		else
		{
			std::cerr << "Error: Cannot separate atoms by type because geometric data not loaded. Please load POSCAR or CHGCAR data before using this option." << std::endl;
			throw std::runtime_error("Error: Geometric data not loaded for atom separation by type");
		}
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
					ions_names.push_back(atom_names[i] + "_" + std::to_string(j+1));
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
	if (header) //writng headers with names for each column
	{
		std::string name;
		file << "# DOS summed " << "\n";
		file << "# Energy (eV)  ";
		for (int type_id = 0; type_id < type_num; type_id++)
		{
			name = ions_names.at(type_id / orb_col) + "_" + orb_names.at(type_id % orb_col);
			file << name << "  ";
		}
		file << "\n";
	}
	file << std::fixed << std::setprecision(8);
	for (int i = 0; i < NDOS; i++) //writing out values
		{
			for (int type_id = 0; type_id < type_num+1; type_id++)
			{
				file << std::setw(10)<< dos_summed(i, type_id) << " ";
			}
			file << "\n";
		}
	file.close();
}
