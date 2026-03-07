#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <armadillo>
#include "VASP_read.h"

using namespace std;


int main()
{
	map<int, string> jobs;
	jobs[1] = "Generate supercells and corresponding POSCAR files";
	jobs[2] = "Charge density analysis and its derviatives (potential, dipole moment, etc.)";
	jobs[3] = "Potential and ionisation energy";
	jobs[4] = "Band structure";
	jobs[5] = "Density of states";
	jobs[6] = "Custom / debugging";
	string job_name;
	cout << "Choose job type:\n";
	for (const auto& job : jobs)
	{
		cout << job.first << ": " << job.second << endl;
	}
	int job_type;
	cin >> job_type;
	if( jobs.find(job_type) != jobs.end())
	{
		job_name = jobs[job_type];
		cout << "You chose: " << job_name << endl;
		try
		{
			switch (job_type)
			{
				case 1:
				{
					/* //supercels for TMDS all from 2 to 8
					string id;
					string body = "POSCAR_";
					vector<string> metals = { "Mo","W" };
					vector<string> chalcogenides = { "S2", "Se2" };
					vector<string> compounds = { "MoSe2", "WS2", "WSe2" };
					vector<string> layers = { "2layer","4layer","6layer", "8layer" };//{ "8layer" } "bulk",;
					vector<int> layer_values = { 2,4,6,8 }; //2,
					VASP_data data_og, data_mod;
					for (const auto& metal : metals)
					{
						for (const auto& chalcogenide : chalcogenides)
						{
							data_og.read_POSCAR("workspace\\" + body + metal + chalcogenide);
							cout << "Processed POSCAR for " << metal + chalcogenide << endl;
							for (int i = 0; i < layers.size(); i++)
							{
								data_mod = data_og.supercell_grid(1, 1, layer_values.at(i) / 2, { 0,0,0,0,1,1 });
								id = metal + chalcogenide + "_" + layers.at(i);
								data_mod.write_POSCAR(id);
							}
						}
					}*/
					//presentation
					/*
					string id;
					string body = "POSCAR";
					vector<string> layers = { "2layer","4layer","6layer", "8layer" };//{ "8layer" } "bulk",;
					vector<int> layer_values = { 2,4,6,8 }; //2,
					VASP_data data_og, data_mod;
					data_og.read_POSCAR("workspace/" + body + "_MoS2");
					for (int i = 0; i < layers.size(); i++)
					{
						data_mod = data_og.supercell_grid(1, 1, layer_values.at(i) / 2, { 0,0,0,0,1,1 });
						id = body + "_MoS2_" + layers.at(i);
						data_mod.write_POSCAR(id);
						cout << "Generated supercell POSCAR for " << id << endl;
					}
					//*/
					VASP_data data_in = VASP_data(), data_out;
					data_in.read_POSCAR("workspace/POSCAR_MoSSe2_4_250");
					data_out = data_in.supercell_grid(1, 1, 1, { 0,0,0,0,1.0/4.0,1.0 / 4.0 });
					data_out.write_POSCAR("MoSSe2_25_vac");

					data_in.read_POSCAR("workspace/POSCAR_MoSSe2_4_500");
					data_out = data_in.supercell_grid(1, 1, 1, { 0,0,0,0,1.0 / 4.0,1.0 / 4.0 });
					data_out.write_POSCAR("MoSSe2_50_vac");

					data_in.read_POSCAR("workspace/POSCAR_MoSSe2_4_750");
					data_out = data_in.supercell_grid(1, 1, 1, { 0,0,0,0,1.0 / 4.0,1.0 / 4.0 });
					data_out.write_POSCAR("MoSSe2_75_vac");

					break;
				}
				case 2:
			{
				//*/
		/*
	try
	{

		int columns = metals.size() * chalcogenides.size();
		int rows = 0, between_rows = 0;
		for (const int& num : layer_values)
		{
			rows += num;
			between_rows += num - 1;
		}
		arma::vec** dipole_moments = new arma::vec * [rows];
		vector<arma::vec> positions, loc_dip;
		double** potenials = new double* [between_rows];
		arma::vec** forces = new arma::vec * [between_rows];
		for (int i = 0; i < rows; i++)
		{
			dipole_moments[i] = new arma::vec[columns];
			if (i < between_rows)
			{
				potenials[i] = new double[columns];
				forces[i] = new arma::vec[columns];
			}
		}

		bool real_calculation = true;
		if (real_calculation)
		{
			// real calculation commented out for faster testing
			int row_id, between_row_id, col_id = 0;
			for (const auto& chalcogenide : chalcogenides)
			{
				for (const auto& metal : metals)
				{
					row_id = 0;
					between_row_id = 0;
					for (const auto& lay : layer)
					{
						string id = metal + chalcogenide + lay;
						try {

							readout data = read_CHGCAR(path, body, id);
							arma::vec Mo, S1, S2, ion_dipol, el_dipol, dipol;
							vector<int> start = { 0,0,0 };
							vector<int> end = data.NGiF;
							arma::mat base = data.cell_matrix.t();
							double height = (base * data.types_atom_positions[0].col(1) - base * data.types_atom_positions[0].col(0))[2];
							int height_in_mesh = static_cast<int>(floor(height / base(2, 2) * data.NGiF[2] + 0.5));
							positions.clear();
							loc_dip.clear();
							for (int i = 0; i < data.atoms_per_type[0]; i++)
							{
								Mo = base * data.types_atom_positions[0].col(i);
								S1 = base * data.types_atom_positions[1].col(2 * i);
								S2 = base * data.types_atom_positions[1].col(2 * i + 1);
								int middle = get_mesh_indices(data.cell_matrix, data.NGiF, Mo)[2];
								start.at(2) = middle - height_in_mesh / 2;
								end.at(2) = middle + height_in_mesh / 2;
								positions.push_back(Mo);
								ion_dipol = (S1 - Mo) * 6 + (S2 - Mo) * 6;
								el_dipol = calc_dipole_mom(data, Mo, start, end);
								dipol = ion_dipol + el_dipol;

								dipole_moments[row_id][col_id] = dipol;
								loc_dip.push_back(dipol);
								row_id++;
								cout << "Calculated dipole moments for " << id << ", layer " << i + 1 << endl;
							}
							for (int i = 0; i < positions.size()-1; i++)
							{
								arma::vec R = positions.at(i + 1) - positions.at(i);
								arma::vec mu_1 = loc_dip.at(i);
								arma::vec mu_2 = loc_dip.at(i + 1);
								potenials[between_row_id][col_id] = calc_dip_dip_potential(mu_1, mu_2, R);
								forces[between_row_id][col_id] = calc_dip_dip_force(mu_1, mu_2, R);
								between_row_id++;

							}
							cleanup_memory(data);
						}
						catch (const std::exception& ex) {
							cerr << "Error at " << id << ex.what() << endl;
						}
					}
					col_id++;
				}
			}
			cout << "\n\n";
		}
		else
		{
			// Mock data for testing
			for (int i = 0; i < columns; i++)
			{
				for (int j = 0; j < rows; j++)
				{
					dipole_moments[j][i] = arma::vec(3, arma::fill::randu);
					if (j < between_rows)
					{
						potenials[j][i] = static_cast<double>(rand()) / RAND_MAX;
						forces[j][i] = arma::vec(3, arma::fill::randu);
					}
				}
			}
		}

		// Writing dipole moments to file
		fstream file;
		arma::vec dip,force;
		file.open("dipole_moments_and_potentials.txt", ios::out);
		bool writing_header = true;
		if (writing_header)
		{
			file << "\tcompund\tWS2\t\t\t\t\t\t\tMoS2\t\t\t\t\t\t\tWSe2\t\t\t\t\t\t\tMoSe2\t\t\t\t\t\t\t\n";
			file << "layers\tlayer    quantity\td_x\td_y\td_z\tU\tF_x\tF_y\tF_z\td_x\td_y\td_z\tU\tF_x\tF_y\tF_z";
			file << "\td_x\td_y\td_z\tU\tF_x\tF_y\tF_z\td_x\td_y\td_z\tU\tF_x\tF_y\tF_z\n";
		}
		int count = 0,count_off = 0;
		for (int pom = 0; pom < layer_values.size(); pom++)
		{
			int j = layer_values[pom];
			for (int k = 0; k < j; k++)
			{

				if (writing_header)
				{
					if (k == 0)
					{

						file << layer.at(pom) << "\t" << k + 1;
					}
					else
					{
						file << "\t" << k + 1;
					}
					file << "\t";
				}
				for (int i = 0; i < columns; i++)
				{
					dip = dipole_moments[count][i];
					file << dip(0) << "\t" << dip(1) << "\t" << dip(2) << "\t";
					file << "\t\t\t\t";
					//cout << "Written dipole moments for " << layer.at(pom) << " layer " << k + 1 << ", type " << dip_type.at(l) << ", material:" <<
					//	metals.at(i % 2) + chalcogenides.at(i / 2 % 2) << endl;
				}
				file << endl;
				count++;
				if (writing_header) file << "\t\t";
				if (k< j - 1)
				{
					for (int i = 0; i < columns; i++)
					{
						file << "\t\t\t";
						file << potenials[count_off][i] << "\t";
						force = forces[count_off][i];
						file << force(0) << "\t" << force(1) << "\t" << force(2) << "\t";
						//cout << "Written dipole moments for " << layer.at(pom) << " layer " << k + 1 << ", type " << dip_type.at(l) << ", material:" <<
						//	metals.at(i % 2) + chalcogenides.at(i / 2 % 2) << endl;
					}
					file << endl;
					count_off++;
				}
				else
				{
					file << "\t\t\t\t\t\t\t";
					file << "\t\t\t\t\t\t\t";
					file << "\t\t\t\t\t\t\t";
					file << "\t\t\t\t\t\t\t";
					file << endl;
				}



			}

		}
		file.close();

		for (int i = 0; i < rows; i++) {
			delete[] dipole_moments[i];
			if (i < between_rows)
			{
				delete[] potenials[i];
				delete[] forces[i];
			}
		}
		delete[] dipole_moments;
		delete[] potenials;
		delete[] forces;

	}
	catch (const std::exception& ex) {
		cerr << "Error at " << "MoS2_bulk:" << ex.what() << endl;
	}
	//*/	
				break;
			}			
				case 3:
				{
					/*
	string body1 = "LOCPOT_", body2 = "EIGENVAL_";
	vector<string> metals = { "Mo","W" };
	vector<string> chalcogenides = { "S2", "Se2" };
	vector<string> compounds = { "MoSe2", "WS2", "WSe2" };
	vector<string> layers = { "2layer","4layer","6layer", "8layer" };//{ "8layer" } "bulk",;
	vector<int> layer_values = { 2,4,6,8 }; //2,
	string  id,filename1,filename2;
	VASP_data data = VASP_data();
	vector<double> potential_z, val, vac;
	double pom;
	int zzz;
	fstream file;
	for (int i = 0; i < compounds.size(); i++)
	{
		val.clear();
		vac.clear();
		for (int j = 0; j < layers.size(); j++)
		{
			id = compounds.at(i) + "_" + layers.at(j);
			filename1 = "workspace\\" + body1  + id;
			filename2 = "workspace\\" + body2  + id;
			data.read_LOCPOT(filename1);
			data.read_EIGENVAL(filename2);
			if(layers.at(j) == "2layer") potential_z = data.sum_potential_averaged_xy_z("primitive");
			else potential_z = data.sum_potential_averaged_xy_z("layered");
			pom = data.find_band_extremum(data.find_valence_band(), true, zzz, true);
			val.push_back(pom);
			pom = potential_z[0];
			vac.push_back(pom);
			data.write_potential_z("potential_z_" + id, potential_z);
			cout << "Processed " << id << ": Valence band maximum energy = " << val.back() << " eV, Vacuum level = " << vac.back() << " eV" << endl;
		}
		filename1 = "workspace\\potential_and_ionisation_energies_" + compounds.at(i) + ".txt";
		file.open(filename1, ios::out);
		file << "#Layer\tValence Band Maximum (eV)\tVacuum Level (eV)\tIonisation Energy (eV)\n";
		for (int j = 0; j < layers.size(); j++)
		{
			double ionisation_energy = vac.at(j) - val.at(j);
			file << layers.at(j) << "\t" << val.at(j) << "\t" << vac.at(j) << "\t" << ionisation_energy << "\n";
		}
		file.close();
	}


	//*/

					/*
job_name = "MoS2_vacum_pot_layers";
VASP_data data = VASP_data();
arma::mat cell_matrix;
int pom;
for(const auto& layer: layers)
{
	data.read_LOCPOT("workspace\\MoS2_POT\\" + body  + layer);
	cell_matrix = data.get_cell_matrix().t();
	if(layer == "bulk")
	{
		data.write_potential_averaged_xy_z("MoS2_avg_pot_z_" + layer, "primitive");
	}
	else if(layer == "2layer")
	{
		pom = floor(0.5 + (data.get_mesh_indices(cell_matrix.col(2))[2] / 3.0));
		data.write_potential_averaged_xy_z("MoS2_avg_pot_z_" + layer, "manual", pom);
	}
	else
	{
		data.write_potential_averaged_xy_z("MoS2_avg_pot_z_" + layer, "layered");
	}
}
//*/

					/*
string body1 = "LOCPOT_", body2 = "EIGENVAL_";
string id1 = "MoS2_2_no_so", id2 = "MoS2_2_so";
vector<double> en, pot;
data.read_LOCPOT("workspace\\" + body1 + id1);
data.read_EIGENVAL("workspace\\" + body2 + id1);
arma::mat cell_matrix = data.get_cell_matrix().t();
pom = floor(0.5 + (data.get_mesh_indices(cell_matrix.col(2))[2] / 3.0));
vector<double> pot_z = data.sum_potential_averaged_xy_z( "manual", pom);
res = data.find_band_extremum(data.find_valence_band(), true, dummy, true);
en.push_back(res);
pot.push_back(pot_z.at(0));
data.write_potential_z(id1, pot_z);

data.read_LOCPOT("workspace\\" + body1 + id2);
data.read_EIGENVAL("workspace\\" + body2 + id2);
pot_z = data.sum_potential_averaged_xy_z("manual", pom);
res = data.find_band_extremum(data.find_valence_band(), true, dummy, true);
en.push_back(res);
pot.push_back(pot_z.at(0));
data.write_potential_z(id2, pot_z);

cout << "Processed " << id1 << ": Valence band maximum energy = " << en.at(0) << " eV, Vacuum level = " << pot.at(0) << " eV";
cout << " ionisation energy = " << pot.at(0) - en.at(0) << " eV" << endl;
cout << "Processed " << id2 << ": Valence band maximum energy = " << en.at(1) << " eV, Vacuum level = " << pot.at(1) << " eV";
cout << " ionisation energy = " << pot.at(1) - en.at(1) << " eV" << endl;
*/

					//presentation
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
					break;
				}
				case 4:
				{
					/* // band structure for TMDS with and without spin orbit coupling bulk and monolayer
				string id1 = "MoSe2", id2 = "WS2", layer = "mono", so1 = "no_so", so2 = "so";
				vector<string> id_list = { id1, id2 };
				vector<string> so_list = { so1, so2 };
				vector<string> layer_list = { layer, ""};
				string body = "EIGENVAL";
				VASP_data data = VASP_data();

				for (int i = 0; i < id_list.size(); i++)
				{
					for(int j = 0; j < so_list.size(); j++)
					{
						//for(int k = 0; k < layer_list.size(); k++)
						//{
						int k = 0;
							if(layer_list.at(k) == "") job_name = body + "_" + id_list.at(i) + "_" + so_list.at(j);
							else job_name = body + "_" + layer_list.at(k) + "_" + id_list.at(i) + "_" + so_list.at(j);
							cout<< "Processing " << job_name << endl;
							data.read_EIGENVAL("workspace/"+job_name);
							data.write_BS(job_name, true, true);
						//}
					}
				}
				//*/

					//presentation
					string body = "EIGENVAL_";
					string id = "MoS2";
					VASP_data data = VASP_data();
					data.read_EIGENVAL("workspace/" + body + id);
					data.write_BS(id, true, true);
					break;
				}
				case 5:
				{
					//test / presentation
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
					break;
				}
				case 6:
				{
					///*
					job_name = "debug";
					VASP_data data = VASP_data();
					data.read_bestsqs("workspace/bestsqs.out");
					data.write_POSCAR("bestsqs");

					/*
					int pom, dummy;
					double res;
					data.read_POSCAR("workspace/POSCAR");
					data.read_DOS("workspace/DOSCAR", false);

					for (int i = 0; i < 3; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							arma::mat result = data.sum_DOS_types(i, j);

							data.write_DOS_sum_types("test_" + to_string(i) + "_" + to_string(j), result, i, j, true);
						}
					}
					*/
					//POSCAR is needed to write it properly
					int stop = 0;
					//Add check to see if POSCAr is read before hand, this makes sense without POSCAR or set of atoms it's impossible to separate the ions how intended

					//data.read_POSCAR("workspace/POSCAR");
					//VASP_data data2 = data.supercell_grid(3,3,2,{1,1,1,1,1,1});
					//data2.write_POSCAR("test");


					//*/
					break;
				}
			}
		}
		catch (const std::exception& ex) {
			cerr << "Error at " << job_name << ": " << ex.what() << endl;
		}
	}
	else
	{
		cerr << "Invalid job type. Exiting." << endl;
		return 1;
	}
	return 0;
}