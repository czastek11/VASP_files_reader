# VASP_files_reader
This program contains class and it's method that allow to read and calculate certain properties out of VASP output files. Most of functionality of implemented methods can be done by other program or manually, however this usually takes a lot of time. The purpose of this program is the automisation of such tasks with higher control. Some methods at new functionalities that usually need to be done manually.
Currently one can do following tasks using methods implemented:
-read POSCAR file and then make multiplied grid of given cell with options to add vacuum on either side and direction.
-read band structure from EIGENVAL and write it out in more compact form for plotting or futher operation.
-read and sort data from CHGCAR and LOCPOT and then write them out or sum them in specific ways.
-read DOSCAR and sum up contributions according to provided settings
-automaticly find valence band and any band global maxima
-calculate dipoles, their potnetial and forces in specified range (this option is working correctly but it's physical interpretation is not clear to me if it's correct)



Documentation of all methods and functions.

VASP_data: There are two constructors right now.
==================================================================================================================================
- () no arguments. Defauly empty constructor.
- (std::string file_path, bool read_POSCAR, bool read_CHGCAR, bool read_LOCPOT, bool read_DOS, bool read_EIGENVAL)

    arguments are string of file path, and logical flags for reading given file type
    file_path the directory where the files are, bool read_XXX true or false to read specific file
	CHGCAR and LOCPOT already contain all POSCAR informations so it's redundant to enbale reading POSCAr while reading those
	constructor strarting with loading any selected files with default options


sorting_positions(arma::mat positions, std::string method)
==================================================================================================================================
-arguments are matrix of positions and a string naming the method of sorting

	sort positions of ions inside each atom type group, available methods:
	"z_rising" - sorts by the third reduced coordinate ascending, usefull for TMDS


std::vector<int> get_mesh_indices(arma::vec pos)
==================================================================================================================================
-argument is vector of position in cartesian coordinates
-output is array of intigers pointing to position in the mesh space

    gives position in mesh grid correspodning to given cartesian psotion or closest fit


read_POSCAR(std::string filename)
==================================================================================================================================
-argument is the name and path of the file

	reads all informations from POSCAR file


read_CHGCAR(std::string filename)
==================================================================================================================================
-argument is the name and path of the file

	reads the data from CHGCAR file, automaticlly fills POSCAR file informations


read_LOCPOT(std::string filename)
==================================================================================================================================
-argument is the name and path of the file

	reads the data from LOCPOT file, automaticlly fills POSCAR file informations


write_POSCAR(std::string filename)
==================================================================================================================================
-argument is the string name of the written file and it's path

    writes out POSCAR file with informations stored in object, with the same formating as VASP POSCAR file


int count_total_electrons()
==================================================================================================================================
-no arguments
-output is intiger number of electrons

    counts number of total electrons from loaded CHGCAR data


std::vector<double> sum_potential_averaged_xy_z(std::string period_type, int period)
==================================================================================================================================
-arguments are the string naming the scheme of moving average and if selected with "manual" the int length of moving average window in mesh space
-output is an array of averaged potential along a,b with whole array going trough c direction

    This calculated the potential averaged over first and second lattice vector and averaged of a moving window in third direction
    The options for controlling the window size are:
    "manual" - window is specified by the user
    "primitive" - window is the size of the cell (only in 3 direction now)
	"layered" - finds third coordinate of first and third type of first ion (it was used for layered TMDS, this why it's so specific for now) and then take distance between them as window
	"none" - skip calculating moving average


write_potential_z(std::string filename, std::vector<double> potential_z)
==================================================================================================================================
-arguments  are string the name of the file where it is written into and previosly calculated double array (vector<double>) of averaged potential

    writes out averaged potential that was calculated from method sum_potential_averaged_xy_z


arma::vec calc_dipole_moment(arma::vec center, std::vector<int> start, std::vector<int> end)
==================================================================================================================================
-arguments are the vector where the dipole moment is calculated, start and end are the arrays (vector<int>) of starting and ending point of 'prostopadłościan' from which the charge denisty are summed from
-output is the dipole vector in cartesian coordinates

	calculates the diople moment relative to center position, start and end are starting and ending range for mesh points to calculate over
	

write_potential(std::string filename)
==================================================================================================================================
-argument is the string file name and path where potnetial is to be written

    write out potential in CHGCAR in format that is easy to plot : x y z value


read_DOS(std::string filename, bool spin_orbit)
==================================================================================================================================
-arguments are the string path and filename of DOSCAR file and bool flag if the calculations were done with spin orbit
    
    read DOS to object. SO version not yet implemented


arma::mat sum_DOS_types(int atoms_sep_type,int orbitals_sep_type)
==================================================================================================================================
-arguments are int type of summation over ions and int type of sumation over orbitals
-output is matrix(arma::mat) of summed DOS values, where n_rows is NDOS and first column is energy while rest is specified combiantion of ion and orbital

	sum DOS contributions in specific way for atoms and orbitals:
	0 - sum all ions/orbitals(sum over l and m)
	1 - sum all ions of atoms(same atom ions)/orbital types (sum over m)
	2 - don't sum at all, get seperate contribution for each ion/orbital 


read_EIGENVAL(std::string filename)
==================================================================================================================================
-argument is string name of the file and it's path

	read information from EIGENVAL: kpoints, their indexes, their wieght, all the band energies for each and all their occupations


write_BS(std::string filename, bool verbose_kpts, bool only_path)
==================================================================================================================================
-arguments are string filename in /workspace/ where the bands structure is supposed to be written, bool flag whetever to write out kpoints themself and bool flag wheteever to skip weighted kpoints
	
    write out band structure in more compact way for graphing, verbose_kpt is to enable writing out kpoints values
	in reduced form (not cartesian), only_path skips writing points that have integration weight, so that only band structure path is written. User can skip writing the falgs to use default setting of false false
	

write_DOS_sum_types(std::string id, const arma::mat& dos_summed, int atoms_sep_type,int orbitals_sep_type, bool header)
==================================================================================================================================
-arguments are string how should the output file in /workspace/ be called, the summed DOS as the matrix (arma::mat) obtained from sum_DOS_types, types of summation for ions and for orbitals used in summation and bool flag whatever to write out header with column names in the file

	writes down the summed DOS from sum_DOS_types. You need to provide the same settings for atoms_sep_type
	and orbitals_sep_type as before. header controls if the header with the name of the columns is written into the file   


arma::mat get_cell_matrix()
==================================================================================================================================
-no arguments
-output is the matrix containing all lattice base vectors of the system as written in POSCAR

    get function for the cell base vectors 


arma::mat get_BS()
==================================================================================================================================
-no arguments
-output is the matrix containing all of bandstructure energies

    get function for the bandstructure of the system





arma::mat get_occupations()
==================================================================================================================================
-no arguments
-output is the matrix containing occupations for all bands and kpoints

    get function for the occupations in the system


arma::rowvec find_kpoint_energy(arma::rowvec kpt, bool weight, int& index)
==================================================================================================================================
-arguments are row vector with kpoint reduced coordinates, logic flag if to search within weighted kpoints and passed by reference index of kpoint
-output is row vector of band structure for specified kpoint

	get rowvector of all bands energy at specified kpoint, index return what index that kpoints is. Weight controls if the method
	looks trough weighted kpoints


double find_band_extremum(int band_index, bool weight, int& kpt_index, bool maxormin)
==================================================================================================================================
-arguments are int index of the band where extremum is searched, logic flag if to look trough weighted kpoints, int by refrence index of found kpoint and logic flag if to look for maximum or minimum
-output is double value valalue of found extremum energy

	get the value and kpoint index corresponding to global extremum of given band. maxormin controls if it's maximum(true) or minimum(false)


int find_valence_band()
==================================================================================================================================
-no arguments
-output is int band index of found band

	gives the index of highest occupied band


VASP_data supercell_grid(int rep_x, int rep_y, int rep_z,std::vector<bool> add_vacuum)
==================================================================================================================================
-arguments are int multiplication of cell in first, second and third direction provided by cell lattice base vectors, a 6 element long array of flags (vectro<bool>) if to add vacuum of the length of given vector to below and/or above cell in given direction {belowx,abovex,belowy,abovey,belowz,abovez}
-output is new object of VASP_data class with generated supercell geometry sotred as new POSCAR

	generate supercell as new VASP_data object. rep_x,rep_y,rep_z controls how many times in given direction the cell is multiplied vacuum is array of logical values to set if vacum is to be generated. two for each direction for below and above cell.


static double calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R)
==================================================================================================================================
-arguments are  three vectors: first dipole, second dipole and postion of one from point of view fo second (separtion between them)
-output is double value of calculated dipole ponetial

	calculate potential between two dipoles, dip1 and dip2 are two dipoles and R is the vector of their separation R=r_1-r_2


static arma::vec calc_dip_dip_force(arma::vec dip_1, arma::vec dip_2, arma::vec R)
==================================================================================================================================
-arguments are  three vectors: first dipole, second dipole and postion of one from point of view fo second (separtion between them)
-output is vector of force acting between those dipoles
	calculate force from the potential betwee  two dipoles
	dip1 and dip2 are two dipoles and R is the vector of their separation R=r_1-r_2


==================================================================================================================================
==================================================================================================================================
Examples of usage in code:


Generation of stacked MoS2 layers:
	string id;
	string body = "POSCAR";
	vector<string> layers = { "2layer","4layer","6layer", "8layer" };
	vector<int> layer_values = { 2,4,6,8 };
	VASP_data data_og, data_mod;
	data_og.read_POSCAR("workspace/" + body + "_MoS2");
	for (int i = 0; i < layers.size(); i++)
	{
		data_mod = data_og.supercell_grid(1, 1, layer_values.at(i) / 2, { 0,0,0,0,1,1 });
		id = body + "_MoS2_" + layers.at(i);
		data_mod.write_POSCAR(id);
		cout << "Generated supercell POSCAR for " << id << endl;
	}


Calculating averaged potentian along z direction in MoS2 with calculation of ionisation energy:
	string body1 = "LOCPOT_", body2 = "EIGENVAL_";
	string id = "MoS2_8layer";
	int pom;
	VASP_data data = VASP_data();
	data.read_LOCPOT("workspace/" + body1 + id);
	data.read_EIGENVAL("workspace/" + body2 + id);
	vector<double> potential_z = data.sum_potential_averaged_xy_z("layered");
	double valence_band_max = data.find_band_extremum(data.find_valence_band(), true, pom, true);
	double vacuum_level = potential_z.at(0);
	double ionisation_energy = vacuum_level - valence_band_max;
	data.write_potential_z("potential_z_" + id, potential_z);
	cout << "Processed " << id << ": Valence band maximum energy = " << valence_band_max << " eV, Vacuum level = " << vacuum_level << " eV";


Writing out band structure in format usable for plotting:
	string body = "EIGENVAL_";
	string id = "MoS2";
	VASP_data data = VASP_data();
	data.read_EIGENVAL("workspace/" + body + id);
	data.write_BS(id, true, true);


Reading andthen writting out summed DOS contribution. In this example all possible combinations modes for summing are performed:
	VASP_data data = VASP_data();
	int pom, dummy;
	double res;
	data.read_POSCAR("workspace/POSCAR_MoSSe2");
	data.read_DOS("workspace/DOSCAR", false);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			arma::mat result = data.sum_DOS_types(i, j);
			data.write_DOS_sum_types("test_" + to_string(i) + "_" + to_string(j), result, i, j, true);
		}
	}