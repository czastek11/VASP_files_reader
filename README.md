# VASP_files_reader

This program contains a class and its methods that allow reading and calculating certain properties from VASP output files. Most functionality of implemented methods can be done by other programs or manually, however this usually takes a lot of time. The purpose of this program is the automation of such tasks with higher control. Some methods offer new functionalities that without this program usually need to be done manually.

Currently one can do the following tasks using implemented methods:
- Read POSCAR file and then create a multiplied grid of a given cell with options to add vacuum on either side and direction
- Read band structure from EIGENVAL and write it out in a more compact form for plotting or further operations
- Read and sort data from CHGCAR and LOCPOT, then write them out or sum them in specific ways
- Read DOSCAR and sum up contributions according to provided settings
- Automatically find the valence band and any band global maxima
- Calculate dipoles, their potential, and forces in a specified range (this option is working correctly but its physical interpretation is not clear to me if it's correct)

## Documentation of all methods and functions

### VASP_data: Constructors

There are two constructors currently:

- **`()`** - No arguments. Default empty constructor.
- **`(std::string file_path, bool read_POSCAR, bool read_CHGCAR, bool read_LOCPOT, bool read_DOS, bool read_EIGENVAL)`**
  
  Arguments:
  - `file_path`: the directory where the files are located
  - `read_XXX`: boolean flags (true/false) to read specific file types
  
  Note: CHGCAR and LOCPOT already contain all POSCAR information, so it's redundant to enable reading POSCAR while reading those files. This constructor starts by loading any selected files with default options.

---

### `sorting_positions(arma::mat positions, std::string method)`

- **Arguments**: matrix of positions and a string naming the method of sorting
- **Description**: Sorts positions of ions inside each atom type group. Available methods:
  - `"z_rising"` - sorts by the third reduced coordinate ascending, useful for TMDs

---

### `std::vector<int> get_mesh_indices(arma::vec pos)`

- **Arguments**: vector of position in cartesian coordinates
- **Output**: array of integers pointing to position in the mesh space
- **Description**: Gives position in mesh grid corresponding to given cartesian position or closest fit

---

### `read_POSCAR(std::string filename)`

- **Arguments**: name and path of the file
- **Description**: Reads all information from a POSCAR file

---

### `read_bestsqs(std::string filename)`

- **Arguments**: name and path of the file
- **Description**: Reads all information from best sqs output and then converts into POSCAR data

---

### `read_CHGCAR(std::string filename)`

- **Arguments**: name and path of the file
- **Description**: Reads the data from a CHGCAR file, automatically fills POSCAR file information

---

### `read_LOCPOT(std::string filename)`

- **Arguments**: name and path of the file
- **Description**: Reads the data from a LOCPOT file, automatically fills POSCAR file information

---

### `write_POSCAR(std::string filename)`

- **Arguments**: string name of the written file and its path
- **Description**: Writes out a POSCAR file with information stored in the object, with the same formatting as a VASP POSCAR file

---

### `int count_total_electrons()`

- **Arguments**: none
- **Output**: integer number of electrons
- **Description**: Counts the total number of electrons from loaded CHGCAR data

---

### `std::vector<double> sum_potential_averaged_xy_z(std::string period_type, int period)`

- **Arguments**: 
  - `period_type`: string naming the scheme of moving average
  - `period`: if selected with "manual", the integer length of moving average window in mesh space
- **Output**: array of averaged potential along a,b with whole array going through c direction
- **Description**: Calculates the potential averaged over first and second lattice vectors and averaged over a moving window in the third direction. Options for controlling the window size:
  - `"manual"` - window is specified by the user
  - `"primitive"` - window is the size of the cell (only in 3 direction now)
  - `"layered"` - finds third coordinate of first and third type of first ion (used for layered TMDs, specific for now) and then takes distance between them as window
  - `"none"` - skip calculating moving average

---

### `write_potential_z(std::string filename, std::vector<double> potential_z)`

- **Arguments**: 
  - `filename`: name of the file where it is written
  - `potential_z`: previously calculated double array (vector<double>) of averaged potential
- **Description**: Writes out averaged potential that was calculated from method `sum_potential_averaged_xy_z`

---

### `arma::vec calc_dipole_moment(arma::vec center, std::vector<int> start, std::vector<int> end)`

- **Arguments**: 
  - `center`: vector where the dipole moment is calculated
  - `start`, `end`: arrays (vector<int>) of starting and ending points of the cuboid from which charge densities are summed
- **Output**: dipole vector in cartesian coordinates
- **Description**: Calculates the dipole moment relative to center position; start and end define the range of mesh points to calculate over

---

### `write_potential(std::string filename)`

- **Arguments**: string file name and path where potential is to be written
- **Description**: Writes out potential in CHGCAR in a format that is easy to plot: x y z value

---

### `read_DOS(std::string filename, bool spin_orbit)`

- **Arguments**: 
  - `filename`: string path and filename of DOSCAR file
  - `spin_orbit`: boolean flag if calculations were done with spin-orbit coupling
- **Description**: Reads DOS into object. SO version not yet implemented

---

### `arma::mat sum_DOS_types(int atoms_sep_type, int orbitals_sep_type)`

- **Arguments**: 
  - `atoms_sep_type`: integer type of summation over ions
  - `orbitals_sep_type`: integer type of summation over orbitals
- **Output**: matrix (arma::mat) of summed DOS values, where n_rows is NDOS and first column is energy while the rest are specified combinations of ion and orbital
- **Description**: Sums DOS contributions in specific ways for atoms and orbitals:
  - `0` - sum all ions/orbitals (sum over l and m)
  - `1` - sum all ions of the same atom type/orbital types (sum over m)
  - `2` - don't sum at all, get separate contribution for each ion/orbital

---

### `read_EIGENVAL(std::string filename)`

- **Arguments**: string name of the file and its path
- **Description**: Reads information from EIGENVAL: k-points, their indices, their weight, all band energies for each, and all their occupations

---

### `write_BS(std::string filename, bool verbose_kpts, bool only_path)`

- **Arguments**: 
  - `filename`: string filename in `/workspace/` where the band structure is to be written
  - `verbose_kpts`: boolean flag whether to write out k-points themselves
  - `only_path`: boolean flag whether to skip weighted k-points
- **Description**: Writes out band structure in a more compact way for graphing. `verbose_kpts` enables writing out k-point values in reduced form (not cartesian); `only_path` skips writing points that have integration weight, so that only the band structure path is written. User can skip writing flags to use default setting of false, false

---

### `write_DOS_sum_types(std::string id, const arma::mat& dos_summed, int atoms_sep_type, int orbitals_sep_type, bool header)`

- **Arguments**: 
  - `id`: string name for the output file in `/workspace/`
  - `dos_summed`: summed DOS matrix (arma::mat) obtained from `sum_DOS_types`
  - `atoms_sep_type`, `orbitals_sep_type`: types of summation for ions and orbitals used in summation
  - `header`: boolean flag whether to write a header with column names in the file
- **Description**: Writes down the summed DOS from `sum_DOS_types`. You need to provide the same settings for `atoms_sep_type` and `orbitals_sep_type` as before. `header` controls if the header with column names is written into the file

---

### `arma::mat get_cell_matrix()`

- **Arguments**: none
- **Output**: matrix containing all lattice base vectors of the system as written in POSCAR
- **Description**: Get function for the cell base vectors

---

### `arma::mat get_BS()`

- **Arguments**: none
- **Output**: matrix containing all band structure energies
- **Description**: Get function for the band structure of the system

---

### `arma::mat get_occupations()`

- **Arguments**: none
- **Output**: matrix containing occupations for all bands and k-points
- **Description**: Get function for the occupations in the system

---

### `arma::rowvec find_kpoint_energy(arma::rowvec kpt, bool weight, int& index)`

- **Arguments**: 
  - `kpt`: row vector with k-point reduced coordinates
  - `weight`: logical flag whether to search within weighted k-points
  - `index`: passed by reference, returns the index of the k-point
- **Output**: row vector of band structure for specified k-point
- **Description**: Gets row vector of all band energies at specified k-point. `weight` controls if the method looks through weighted k-points

---

### `double find_band_extremum(int band_index, bool weight, int& kpt_index, bool maxormin)`

- **Arguments**: 
  - `band_index`: index of the band where extremum is searched
  - `weight`: logical flag whether to look through weighted k-points
  - `kpt_index`: passed by reference, returns the index of found k-point
  - `maxormin`: logical flag whether to look for maximum (true) or minimum (false)
- **Output**: double value of found extremum energy
- **Description**: Gets the value and k-point index corresponding to global extremum of given band

---

### `int find_valence_band()`

- **Arguments**: none
- **Output**: integer band index of found band
- **Description**: Returns the index of the highest occupied band

---

### `VASP_data supercell_grid(int rep_x, int rep_y, int rep_z, std::vector<double> add_vacuum)`

- **Arguments**: 
  - `rep_x`, `rep_y`, `rep_z`: integer multiplication of cell in first, second, and third direction
  - `add_vacuum`: 6-element array of mulitplies (vector<double>) whether to add vacuum of mulitple of the length of the given vector to below and/or above cell in given direction {below_x, above_x, below_y, above_y, below_z, above_z}
- **Output**: new object of VASP_data class with generated supercell geometry stored as new POSCAR
- **Description**: Generates supercell as a new VASP_data object. `rep_x`, `rep_y`, `rep_z` control how many times in each direction the cell is multiplied. `add_vacuum` is array of double values to set how many multiplies of vectors are to be added as vacuum

---

### `static double calc_dip_dip_potential(arma::vec dip_1, arma::vec dip_2, arma::vec R)`

- **Arguments**: three vectors: first dipole, second dipole, and position of one from the point of view of the second (separation between them)
- **Output**: double value of calculated dipole potential
- **Description**: Calculates potential between two dipoles. `dip_1` and `dip_2` are two dipoles, and `R` is the vector of their separation R = r₁ - r₂

---

### `static arma::vec calc_dip_dip_force(arma::vec dip_1, arma::vec dip_2, arma::vec R)`

- **Arguments**: three vectors: first dipole, second dipole, and position of one from the point of view of the second (separation between them)
- **Output**: vector of force acting between the dipoles
- **Description**: Calculates force from the potential between two dipoles. `dip_1` and `dip_2` are two dipoles, and `R` is the vector of their separation R = r₁ - r₂

---
## Examples of usage in code

### Generation of stacked MoS2 layers:

```cpp
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
```
### Calculating averaged potentian along z direction in MoS2 with calculation of ionisation energy:
```cpp
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
```

### Writing out band structure in format usable for plotting:
```cpp
string body = "EIGENVAL_";
string id = "MoS2";
VASP_data data = VASP_data();
data.read_EIGENVAL("workspace/" + body + id);
data.write_BS(id, true, true);
```

### Reading andthen writting out summed DOS contribution. In this example all possible combinations modes for summing are performed:
```cpp
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
```
