#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <armadillo>
#include "VASP_read.h"

using namespace std;

string format_percentage(double value, int digits)
{
	string p = to_string(value);
	size_t dot_pos = p.find('.');
	p.at(dot_pos) = '_'; // replace dot with underscore for filename compatibility
	int length_after_dot = p.length() - dot_pos - 1;
	if (length_after_dot > digits) p = p.substr(0, dot_pos + 1 + digits); // cut off extra digits if there are more than desired
	else if (length_after_dot < digits)
	{
		for (int i = 0; i < digits - length_after_dot; i++)
		{
			p += '0'; // add zeros if needed to reach the desired number of digits after the decimal point
		}
	}
	return p;
}



void job1()
{
	VASP_data data_in = VASP_data(), data_out = VASP_data();
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
	/*
	data_in.read_POSCAR("workspace/POSCAR_MoSSe2_4_250");
	data_out = data_in.supercell_grid(1, 1, 1, { 0,0,0,0,1.0/4.0,1.0 / 4.0 });
	data_out.write_POSCAR("MoSSe2_25_vac");

	data_in.read_POSCAR("workspace/POSCAR_MoSSe2_4_500");
	data_out = data_in.supercell_grid(1, 1, 1, { 0,0,0,0,1.0 / 4.0,1.0 / 4.0 });
	data_out.write_POSCAR("MoSSe2_50_vac");

	data_in.read_POSCAR("workspace/POSCAR_MoSSe2_4_750");
	data_out = data_in.supercell_grid(1, 1, 1, { 0,0,0,0,1.0 / 4.0,1.0 / 4.0 });
	data_out.write_POSCAR("MoSSe2_75_vac");
	//*/
	/*
	//CONTCAR_WS_Se2_0_500_man
	std::vector<double> percentages = { 0.25, 0.50, 0.75 };
	std::vector<std::string> alloys = { "WS_Se2", "W_MoSe2", "W_MoS2" };
	std::vector<std::string> method = { "rand","man" };
	string perc,filename;
	VASP_data datain = VASP_data();
	VASP_data dataout = VASP_data();

	for (int i = 0; i<alloys.size(); i++)
	{
		for (int j = 0; j<percentages.size(); j++)
		{
			perc = format_percentage(percentages.at(j), 3);
			for (int k = 0; k<method.size(); k++)
			{
				filename = "workspace/CONTCAR_" + alloys.at(i) + "_" + perc + "_" + method.at(k);
				datain.read_POSCAR(filename);
				dataout = datain.supercell_grid(1, 1, 1, { 0,0,0,0,1.0,1.0});
				filename = alloys.at(i) + "_" + perc + "_" + method.at(k);
				dataout.write_POSCAR(filename);

			}
		}
	}
	*/

	vector<string> compounds = { "MoSe2","MoS2 ","WS2", "WSe2" };
	for (const auto& compound : compounds)
	{
		data_in.read_POSCAR("POSCAR_" + compound);
		data_out = data_in.supercell_grid(2, 2, 1, { 0,0,0,0,0.0,0.0 });
		data_out.write_POSCAR(compound + "_2x2");
	}

	


}

void job2()
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
}

void job3()
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
					/*
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
					*/
	/*
	std::vector<double> av_pot, percentages;
	std::vector<std::string> names;
	for (int i = 0; i < 4; i++)
	{
		percentages.push_back(i * 0.25 + 0.125);
		names.push_back("");
	}
	std::vector<arma::vec> results;
	string filename1, filename2, filename3, perc;
	int val, pom, win;
	VASP_data data = VASP_data();
	double val_g, vac, ionis;
	arma::vec result = arma::vec(3, arma::fill::zeros);
	for (int i = 0; i < 4; i++)
	{
		pom = static_cast<int>(percentages.at(i) * 1000);
		std::ostringstream ss;
		ss << std::setfill('0') << std::setw(4) << pom;
		perc = ss.str();

		//if (percentages.at(i) < 1.0) perc = "0" + perc;
		filename1 = "workspace/EIGENVAL_MoSSe2_" + perc + names.at(i);
		filename2 = "workspace/energies_MoSSe2_" + perc + names.at(i);
		filename3 = "workspace/LOCPOT_MoSSe2_" + perc + names.at(i);

		cout << "Processing: " << perc + names.at(i) << "\n";

		data.read_EIGENVAL(filename1);
		val = data.find_valence_band();

		data.read_BS(filename2, false, true);
		val_g = data.find_band_extremum(val, true, pom, true);

		data.read_LOCPOT(filename3);
		av_pot = data.average_potential_over(3);
		win = data.get_mesh_indices(data.get_cell_matrix().row(2).t() / 3.0).at(2);
		av_pot = data.moving_average_potential_over(av_pot, 3, "manual", win);
		vac = av_pot.front();//max(av_pot.front(), av_pot.back());

		ionis = vac - val_g;

		cout << "Valence max: " << val_g << " vacuum pot.: " << vac << " ionisation potential: " << ionis << "\n";
		result(0) = val_g; result(1) = vac; result(2) = ionis;
		results.push_back(result);


		//data.write_potential_over("MoSSe2_" + perc + "_2x2" + names.at(i), av_pot, 3);
	}
	fstream file;
	file.open("workspace/ionisation_energy_2x2.txt", ios::out);
	for (int i = 0; i < 4; i++)
	{
		file << percentages.at(i) << "\t" << results.at(i)(0) << "\t" << results.at(i)(1) << "\t" << results.at(i)(2) << "\n";
	}
	//*/
	/*
	std::vector<std::string> percentages = { "0250", "0500", "0750" };
	std::vector<std::string> supercells = { "2x2", "4" };
	std::map<std::string, std::vector<std::vector<double>>> results;
	for (const auto& supercell : supercells)
	{
		std::string filename_output = "workspace/ionisation_energies_" + supercell + ".txt";
		std::ofstream outfile(filename_output);

		if (!outfile.is_open()) {
			std::cerr << "Cannot open output file: " << filename_output << std::endl;
			continue;
		}
		outfile << "# x_percent\tval_G\tvacuum_potential\tionisation_energy\n";
		std::cout << "\nResults for " << supercell << " supercell:\n";

		std::vector<std::vector<double>> supercell_results;

		for (const auto& perc : percentages)
		{
			// Convert percentage to decimal for output
			double x_percent;
			if (perc == "0250") x_percent = 0.25;
			else if (perc == "0500") x_percent = 0.50;
			else if (perc == "0750") x_percent = 0.75;

			std::cout << "\nProcessing " << perc << " " << supercell << std::endl;

			// Construct filenames
			std::string eigenval_file = "workspace/EIGENVAL_MoSSe2_" + perc + "_" + supercell;
			std::string bs_file = "workspace/BS_MoSSe2_" + perc + "_" + supercell;
			std::string locpot_file = "workspace/LOCPOT_MoSSe2_" + perc + "_" + supercell;

			VASP_data data = VASP_data();

			// Read EIGENVAL
			data.read_EIGENVAL(eigenval_file);
			int val = data.find_valence_band();
			int pom;
			data.read_BS(bs_file, false, true);
			double val_G = data.find_band_extremum(val, true, pom, true);


			data.read_LOCPOT(locpot_file);
			std::vector<double> av_pot = data.average_potential_over(3);
			//double ionis_en = *std::max_element(av_pot.begin(), av_pot.end());

			if (supercell == "2x2")
			{
				int win = data.get_mesh_indices(data.get_cell_matrix().row(2).t() / 3.0).at(2);
				av_pot = data.moving_average_potential_over(av_pot, 3, "manual", win);
			}
			else if (supercell == "4") av_pot = data.moving_average_potential_over(av_pot, 3, "layered", 16, 18);

			double ionis_en = av_pot.front();// max(av_pot.front(), av_pot.back());
			double ionisation_energy = ionis_en - val_G;

			// Output to console
			std::cout << "  x = " << x_percent
				<< "  val_G = " << val_G
				<< "  vacuum pot = " << ionis_en
				<< "  ionisation energy = " << ionisation_energy << std::endl;

			outfile << x_percent << "\t"
				<< val_G << "\t"
				<< ionis_en << "\t"
				<< ionisation_energy << "\n";

			// Write potential to file
			//if(x_percent == 0.25 && supercell == "2x2") data.write_potential_over("MoSSe2_" + perc + "_" + supercell, av_pot, 3);

			// Store for summary
			supercell_results.push_back({ x_percent, val_G, ionis_en, ionisation_energy });
		}
		outfile.close();
		results[supercell] = supercell_results;
		std::cout << "\nResults saved to " << filename_output << std::endl;
	}
	*/
	std::vector<double> percentages = { 0,1.0 };// { 0.25, 0.50, 0.75 }, ;
	std::vector<double> av_pot;
	std::vector<std::string> alloys = { "WS_Se2" };// , "W_MoSe2", "W_MoS2"};
	std::vector<std::string> method = { "" };// { "rand", "man" };
	string perc, filename1,filename2,filename3,filename_alloy;
	VASP_data datain = VASP_data();
	VASP_data dataout = VASP_data();
	fstream file;
	arma::mat ionisation_results;
	int val, pom, win;
	double val_G, ionis_en, ionisation_energy;

	for (int k = 0; k < method.size(); k++)
	{
		for (int i = 0; i < alloys.size(); i++)
		{
			ionisation_results = arma::mat(percentages.size(), 4, arma::fill::zeros);
			for (int j = 0; j < percentages.size(); j++)
			{
				perc = format_percentage(percentages.at(j), 3);

				{
					filename1 = "EIGENVAL_" + alloys.at(i) + "_" + perc + "_" + method.at(k);
					filename2 = "energies_" + alloys.at(i) + "_" + perc + "_" + method.at(k);
					filename3 = "LOCPOT_" + alloys.at(i) + "_" + perc + "_" + method.at(k);
					datain.read_EIGENVAL(filename1);
					val = datain.find_valence_band();
					datain.read_BS(filename2, false, true);
					val_G = datain.find_band_extremum(val, true, pom, true);
					datain.read_LOCPOT(filename3);
					av_pot = datain.average_potential_over(3);

					win = datain.get_mesh_indices(datain.get_cell_matrix().row(2).t() / 3.0).at(2);
					av_pot = datain.moving_average_potential_over(av_pot, 3, "manual", win);

					ionis_en = max(av_pot.front(), av_pot.back()); // (av_pot.front() + av_pot.back()) / 2.0;
					ionisation_energy = ionis_en - val_G;
					cout << "Processed " << alloys.at(i) << " " << perc << " " << method.at(k) << ": Valence band max = " << val_G << " eV, Vacuum level = " << ionis_en << " eV, Ionisation energy = " << ionisation_energy << " eV\n";

					ionisation_results(j, 0) = percentages.at(j);
					ionisation_results(j, 1) = val_G;
					ionisation_results(j, 2) = ionis_en;
					ionisation_results(j, 3) = ionisation_energy;

					datain.write_potential_over(alloys.at(i) + "_" + perc + "_" + method.at(k), av_pot, 3);
				}
			}

			filename_alloy = "output/ionisation_energies_" + alloys.at(i) + "_" + method.at(k);
			file.open(filename_alloy, ios::out);
			if (!file.is_open()) {
				std::cerr << "Cannot open output file: " << filename_alloy << std::endl;
				continue;
			}
			else
			{
				file << "# x_percent\tval_G\tvacuum_potential\tionisation_energy\n";
				for (int j = 0; j < percentages.size(); j++)
				{
					file << ionisation_results(j, 0) << "\t" << ionisation_results(j, 1) << "\t" << ionisation_results(j, 2) << "\t" << ionisation_results(j, 3) << "\n";
				}
			}
			file.close();
		}
	}
}

void job4()
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
	/*
	std::vector<double> percentages = { 0.25, 0.50, 0.75 };
	std::vector<double> av_pot;
	std::vector<std::string> alloys = { "MoS_Se2","WS_Se2", "W_MoSe2", "W_MoS2"};
	std::vector<std::string> method = { "man" };// { "rand", "man" };
	string perc, filename1, filename2, filename_alloy;
	VASP_data datain = VASP_data();
	fstream file;
	int val, pom;
	arma::mat result;
	arma::rowvec en_k;
	for (int k = 0; k < method.size(); k++)
	{
		for (int i = 0; i < alloys.size(); i++)
		{
			result = arma::mat(6, 3, arma::fill::zeros);
			for (int j = 0; j < percentages.size(); j++)
			{
				perc = format_percentage(percentages.at(j), 3);
				
				filename1 = "EIGENVAL_" + alloys.at(i) + "_" + perc + "_" + method.at(k);
				filename2 = "energies_" + alloys.at(i) + "_" + perc + "_" + method.at(k);
				datain.read_EIGENVAL(filename1);
				val = datain.find_valence_band();
				datain.read_BS(filename2, false, true);
				result(j*2, 0) = percentages.at(j);
				result(j*2+1, 0) = percentages.at(j);
				en_k = datain.find_kpoint_energy({ 0,0,0 }, true, pom);
				result(j*2, 1) = en_k(val);
				result(j*2, 2) = en_k(val + 1);
				en_k = datain.find_kpoint_energy({ 0.333333,0.333333,0 }, true, pom);
				result(j*2+1, 1) = en_k(val);
				result(j*2+1, 2) = en_k(val + 1);
			}

			filename_alloy = "output/band" + alloys.at(i) + ".txt";
			file.open(filename_alloy, ios::out);
			if (!file.is_open()) {
				std::cerr << "Cannot open output file: " << filename_alloy << std::endl;
				continue;
			}
			else
			{
				file << "# x_percent\tVB\tCB\n";
				for (int j = 0; j < percentages.size(); j++)
				{
					file << std::fixed << std::setprecision(8);
					file << result(j*2, 0) << "\t" << result(j*2, 1) << "\t" << result(j*2, 2) << " #Gamma" << "\n";
				}
				for (int j = 0; j < percentages.size(); j++)
				{
					file << std::fixed << std::setprecision(8);
					file << result(j * 2 + 1, 0) << "\t" << result(j * 2 + 1, 1) << "\t" << result(j * 2 + 1, 2) << " #K" << "\n";
				}
			}
			file.close();
		}
	}
	*/
	vector<string> compunds = { "MoS2", "MoSe2", "WS2", "WSe2" };
	vector<string> meth = { "no_so", "so" };
	VASP_data data = VASP_data();
	for (const auto& comp : compunds)
	{
		for (const auto& method : meth)
		{
			string job_name = method + "_" + comp;
			string filename = "EIGENVAL_" + comp + "__BS_" + method;
			cout << "Processing " << job_name << endl;
			data.read_EIGENVAL(filename);
			data.write_BS("EIGENVAL_R2SCAN_" + comp + "_" + method, true, true);
		}
	}
}

void job5()
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
}

void job6()
{
	/*
					VASP_data data = VASP_data();
					data.read_EIGENVAL("workspace/EIGENVAL_MoSSe2_0250_2x2");
					int val = data.find_valence_band(),pom;
					data.read_BS("workspace/BS_MoSSe2_0250_2x2", false, true);
					double val_G = data.find_band_extremum(val, true, pom, true);
					data.read_LOCPOT("workspace/LOCPOT_MoSSe2_0250_2x2");
					vector<double> av_pot = data.average_potential_over(3);
					double ionis_en = *max_element(av_pot.begin(), av_pot.end());


					int win = data.get_mesh_indices(data.get_cell_matrix().row(2).t()/3.0).at(2);


					av_pot = data.moving_average_potential_over(av_pot, 3, "manual", win); //for 2x2
					av_pot = data.moving_average_potential_over(av_pot, 3, "layered", 16,18);
					//max(av_pot.at(0), av_pot.at(av_pot.size() - 1));
					cout << "Valence band max: " << val_G << " vacuum potential: " << ionis_en << " Ionisation energy: " << ionis_en - val_G << "\n";
					data.write_potential_over("MoSSe2_0250_2x2", av_pot, 3);
					//data.read_bestsqs("workspace/bestsqs.out");
					//data.write_POSCAR("bestsqs");
					//*/
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
					//*/
					//POSCAR is needed to write it properly
					//Add check to see if POSCAr is read before hand, this makes sense without POSCAR or set of atoms it's impossible to separate the ions how intended

					//data.read_POSCAR("workspace/POSCAR");
					//VASP_data data2 = data.supercell_grid(3,3,2,{1,1,1,1,1,1});
					//data2.write_POSCAR("test");
					
					//*/
					VASP_data POSCAR1 = VASP_data();
					VASP_data POSCAR2 = VASP_data();
					vector<string> mixinga, mixingb, mixingcombinations, mixing_at1, mixing_at2;
					string filename1, filename2, filename3;
					mixinga = { "WS2", "WSe2", "WS2" };
					mixingb = { "WSe2", "MoSe2", "MoS2" };
					mixing_at1 = { "S", "W", "W" };
					mixing_at2 = { "Se", "Mo", "Mo" };
					mixingcombinations = { "WSSe2", "WMoSe2", "WMoS2" };
					string job_name_;
					for (int i = 0; i < mixinga.size(); i++)
					{
						filename1 = "workspace/POSCAR_" + mixinga.at(i);
						filename2 = "workspace/POSCAR_" + mixingb.at(i);
						POSCAR1.read_POSCAR(filename1);
						POSCAR2.read_POSCAR(filename2);
						for (double x = 0.25; x < 1.0; x += 0.25)
						{
							job_name_ = "Alloy geometry for " + mixinga.at(i) + " and " + mixingb.at(i) + " with x = " + to_string(x);
							cout << "Processing: " << job_name_ << endl;
							POSCAR1.alloy_geometry(POSCAR2, x, { mixing_at1.at(i) }, { mixing_at2.at(i) }, to_string(x) + "_" + mixingcombinations.at(i));
						}
					}

					//*/

	//POSCAR1.read_POSCAR("workspace/POSCAR_WS2");
	//POSCAR2.read_POSCAR("workspace/POSCAR_WSe2");
	//POSCAR1.alloy_geometry(POSCAR2,0.875,{"S"},{"Se"},"WSSe2");

}

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
					job1();
					break;
				}
				case 2:
				{
					job2();
					break;
				}			
				case 3:
				{
					job3();
					break;
				}
				case 4:
				{
					job4();
					break;
				}
				case 5:
				{
					job5();
					break;
				}
				case 6:
				{
					job6();
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