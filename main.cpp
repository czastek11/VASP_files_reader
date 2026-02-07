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
	string path = "C:\\nauka\\pwr\\DFT\\TDMS\\sprezynki_cd\\charge densities\\";
	//string body = "CHGCAR_";
	string body = "LOCPOT_";
	vector<string> metals = { "W", "Mo" };
	vector<string> chalcogenides = { "S2_", "Se2_" };
	vector<string> layer = { "bulk","2layer","6layer" };
	vector<int> layer_values = { 2,2,6 };

	/*
	for (const auto& metal : metals)
	{
		for (const auto& chalcogenide : chalcogenides)
		{
			for (const auto& lay : layer)
			{
				try {
					string id = metal + chalcogenide + lay;
					readout data = read_CHGCAR(path, body, id);
					write_charge_density_integrated_xyz(id, data);
				}
				catch (...) {
					cerr << "Error processing file for " << metal << " " << chalcogenide << " " << lay << endl;
				}
			}
		}
	}
	*/
	/*
	for (const auto& lay : layer)
	{
		try {
			string id = "MoS2_" + lay;
			readout data = read_CHGCAR(path, body, id);
			cout << "Total electrons for " << id << ": " << count_total_electrons(data) << endl;
			cout << "Volume for " << id << ": " << arma::det(data.cell_matrix) << endl;
			cleanup_memory(data);
		}
		catch (...) {
			cerr << "Error processing file for " << lay << endl;
		}
	}
	*/
	/*
	try
	{
		readout data1 = read_CHGCAR(path, body, "MoS2_bulk");
		readout data2 = read_CHGCAR("", body, "MoS2_bulk_new");
		write_charge_density_integrated_xyz("MoS2_bulk_old", data1);
		write_charge_density_integrated_xyz("MoS2_bulk_new", data2);
		cout << "Total electrons for MoS2_bulk old: " << count_total_electrons_flout(data1) << endl;
		cout << "Total electrons for MoS2_bulk new: " << count_total_electrons_flout(data2) << endl;
	}
	catch (const std::exception& ex) {
		cerr << "Error at " << "MoS2_bulk:" <<ex.what()<< endl;
	}
	//*/
	/*try
	{
		readout data = read_CHGCAR(path, body, "WSe2_6layer");
		arma::vec Mo = { 1.0/3.0,2.0/3.0,0.312506 };
		arma::vec S1 = { 2.0 / 3.0,1.0 / 3.0,0.280910 };
		arma::vec S2 = { 2.0 / 3.0,1.0 / 3.0,0.344103 };
		Mo = data.cell_matrix.t() * Mo;
		S1 = data.cell_matrix.t() * S1;
		S2 = data.cell_matrix.t() * S2;
		vector<int> start = { 0,0,0 };
		vector<int> end = data.NGiF;
		start.at(2) = data.NGiF[2] / 4;
		end.at(2) = data.NGiF[2] * 3 / 8;
		arma::vec dip_mom_ion, dip_mom_elec, dip_mom_total;
		dip_mom_ion = (S1 - Mo) * 6 + (S2 - Mo) * 6;
		dip_mom_elec = calc_dipole_mom(data, Mo, start, end);
		dip_mom_total = dip_mom_ion + dip_mom_elec;
		cout << "Calculating dipole moment around Mo:\n";
		cout << "Ionic dipole moment:\n";
		dip_mom_ion.print();
		cout << "Electronic dipole moment:\n";
		dip_mom_elec.print();
		cout << "Total dipole moment:\n";
		dip_mom_total.print();
	}
	catch (const std::exception& ex) {
		cerr << "Error at dipole moment calculation: " << ex.what() << endl;
	}
	//*/
	/*
	try
	{
		vector<string> dip_type = { "ion","elec","total" };
		int columns = metals.size() * chalcogenides.size();
		int rows = 0;
		for (const int& num : layer_values) rows += num;
		rows *= 3; // three types of dipole moments: ionic, electronic, total
		arma::vec** dipole_moments = new arma::vec * [rows];
		for (int i = 0; i < rows; i++)
		{
			dipole_moments[i] = new arma::vec[columns];
		}

		bool real_calculation = true;
		if (real_calculation)
		{
			// real calculation commented out for faster testing
			int row_id, col_id = 0;
			for (const auto& chalcogenide : chalcogenides)
			{
				for (const auto& metal : metals)
				{
					row_id = 0;
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
							for (int i = 0; i < data.atoms_per_type[0]; i++)
							{
								Mo = base * data.types_atom_positions[0].col(i);
								S1 = base * data.types_atom_positions[1].col(2 * i);
								S2 = base * data.types_atom_positions[1].col(2 * i + 1);
								int middle = get_mesh_indices(data.cell_matrix, data.NGiF, Mo)[2];
								start.at(2) = middle - height_in_mesh / 2;
								end.at(2) = middle + height_in_mesh / 2;

								ion_dipol = (S1 - Mo) * 6 + (S2 - Mo) * 6;
								el_dipol = calc_dipole_mom(data, Mo, start, end);
								dipol = ion_dipol + el_dipol;

								dipole_moments[row_id][col_id] = ion_dipol;
								dipole_moments[row_id + 1][col_id] = el_dipol;
								dipole_moments[row_id + 2][col_id] = dipol;
								row_id += 3;
								cout << "Calculated dipole moments for " << id << ", layer " << i + 1 << endl;
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
				}
			}
		}

		// Writing dipole moments to file
		fstream file;
		arma::vec dip;
		file.open("dipole_moments.txt", ios::out);
		bool writing_header = false;
		if (writing_header)
		{
			file << "\t\tcompund\tWS2\t\t\tMoS2\t\t\tWSe2\t\t\tMoSe2\t\t\n";
			file << "layers\tlayer\ttype/direction\tx\ty\tz\tx\ty\tz\tx\ty\tz\tx\ty\tz\t\n";
		}
		int count = 0;
		for (int pom = 0; pom < layer_values.size(); pom++)
		{
			int j = layer_values[pom];
			for (int k = 0; k < j; k++)
			{

				for (int l = 0; l < 3; l++)
				{
					if (writing_header)
					{
						if (k == 0 && l == 0)
						{

							file << layer.at(pom) << "\t" << k + 1;
						}
						else if (l == 0)
						{
							file << "\t" << k + 1;
						}
						else
						{
							file << "\t";
						}
						file << "\t" << dip_type.at(l) << "\t";
					}
					for (int i = 0; i < columns; i++)
					{
						dip = dipole_moments[count][i];
						file << dip(0) << "\t" << dip(1) << "\t" << dip(2) << "\t";
						//cout << "Written dipole moments for " << layer.at(pom) << " layer " << k + 1 << ", type " << dip_type.at(l) << ", material:" <<
						//	metals.at(i % 2) + chalcogenides.at(i / 2 % 2) << endl;
					}
					file << endl;
					count++;


				}


			}

		}
		file.close();

		for (int i = 0; i < rows; i++) {
			delete[] dipole_moments[i];
		}
		delete[] dipole_moments;

	}
	catch (const std::exception& ex) {
		cerr << "Error at " << "MoS2_bulk:" << ex.what() << endl;
	}
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
	///*
	try
	{
		string id = "MoS2_slab";
		VASP_data data = VASP_data();
		data.read_LOCPOT("", body, id);
		cout << "Volume for "+id+": " << arma::det(data.get_cell_matrix()) << endl;
		data.write_potential_averaged_xy_z(id + "_avg", false);
	}
	catch (const std::exception& ex) {
		cerr << "Error at " << "MoS2_bulk:" << ex.what() << endl;
	}
	//*/
	/*
	try
	{
		string id = "MoS2_slab";
		readout_pot data = read_LOCPOT("", body, id);
		cout << "Volume for " + id + ": " << arma::det(data.cell_matrix) << endl;
		write_potential(id, data);
		cleanup_memory(data);
	}
	catch (const std::exception& ex) {
		cerr << "Error at::"<<ex.what() << endl;
	}
	//*/
	/*
	try
	{
		vector<string> type_names = { "Se","S","Mo" };
		vector<vector<vector<double>>> dos_data = load_DOS("DOSCAR", 24, "LORBIT=11,no_SO");
		vector<int> sets = { 8,8,8 };
		vector<vector<vector<double>>> dos_summed = sum_DOS_types(dos_data, sets);
		write_DOS_sum_types("MoSeS2_0500", dos_summed, type_names);

	}
	catch (const std::exception& ex)
	{
		cerr << "Error at::" << ex.what() << endl;
	}
	//*/
	return 0;
}