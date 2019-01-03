#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <memory>
#include <algorithm>
#include "Molecule.h"
#include "FinderAtom.h"
#include "FinderBond.h"

constexpr int BLANK = 2;
constexpr int DATAWIDTH = 15;
constexpr int DATAPRECISION = 6;

using namespace std;

ofstream debug;

int main(int argc, char **argv)
{
#ifdef DEBUG_MOLECULE
	debug.open("debug.txt", ofstream::out);
#endif // DEBUG_MOLECULE

	{
		ifstream cfg;
		string sys_name = getenv("MOLECULE");
		string cfg_file = getenv("MOLECULE_DIR");
		cfg_file += "/." + sys_name;

		cfg.open(cfg_file.c_str(), ifstream::in);

		Molecule::usingBond();
		Molecule::InputInfo(cfg);
		cfg.close();
	}

	// -r: find bond length according to rule file
	bool opt_r = false;
	// -h: do not print header line
	bool opt_h = false;
	// -e: do not print energy
	bool opt_e = false;
	// -f: find bond length according to one rule line
	bool opt_f = false;	

	{
		char cmd;
		while ((cmd = getopt(argc, argv, "rfhe")) != -1) {
			switch (cmd)
			{
			case 'r':
				opt_r = true;
				break;
			case 'f':
				opt_f = true;
				break;
			case 'h':
				opt_h = true;
				break;
			case 'e':
				opt_e = true;
				break;
			}
		}
	}

	vector<shared_ptr<FinderBase>> rule;
	if (opt_f) {

		Molecule::usingMatrixR();

		istringstream sin;
		sin.str(argv[optind]);

		string var;
		sin >> var;
		if (var == "bond") {
			string iE, jE, sort_type;
			int num;

			sin >> iE >> jE >> sort_type >> num;
			if (sin)
				rule.emplace_back(new FinderBond(iE, jE, sort_type, num));
			else {
				cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
				cerr << "invalid rule" << endl;
				exit(1);
			}
		}
		else if (var == "atom") {
			string aE, aiE, ajE, aSortType;
			string bE, biE, bjE, bSortType;
			int aNum, bNum;

			sin >> aE >> aiE >> ajE >> aSortType >> aNum;
			sin >> bE >> biE >> bjE >> bSortType >> bNum;

			if (aiE == ajE) {
				cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
				cerr << "aiE == ajE" << endl;
				exit(1);
			}
			else if (aE != aiE && aE != ajE) {
				cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
				cerr << "aE != aiE && aE != ajE" << endl;
				exit(1);
			}

			if (biE == bjE) {
				cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
				cerr << "biE == bjE" << endl;
				exit(1);
			}
			else if (bE != biE && bE != bjE) {
				cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
				cerr << "bE != biE && bE != bjE" << endl;
				exit(1);
			}

			if (sin)
				rule.emplace_back(new FinderAtom(aE, aiE, ajE, aSortType, aNum, bE, biE, bjE, bSortType, bNum));
			else {
				cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
				cerr << "invalid rule" << endl;
				exit(1);
			}
		}
		else {
			cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
			cerr << "invalid rule" << endl;
			exit(1);
		}


	}
	else if (opt_r) {

		if (argc - optind < 1) {
			cerr << "BondAnalyze: missing operand" << endl;
			exit(1);
		}

		Molecule::usingMatrixR();

		ifstream fin;
		fin.open(argv[optind], ifstream::in);
		if (!fin) {
			cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
			cerr << "rule file open failed!" << endl;
			exit(1);
		}

		while (fin) {
			int pk = fin.peek();
			if (pk == '#' || pk == '\n' || pk == ' ' || pk == '\r')
				fin.ignore(1024, '\n');

			string var;
			fin >> var;
			if (var == "bond") {
				string iE, jE, sort_type;
				int num;

				fin >> iE >> jE >> sort_type >> num;
				rule.emplace_back(new FinderBond(iE, jE, sort_type, num));
			}
			else if (var == "atom") {
				string aE, aiE, ajE, aSortType;
				string bE, biE, bjE, bSortType;
				int aNum, bNum;

				fin >> aE >> aiE >> ajE >> aSortType >> aNum;
				fin >> bE >> biE >> bjE >> bSortType >> bNum;

				if (aiE == ajE) {
					cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
					cerr << "aiE == ajE" << endl;
					exit(1);
				}
				else if (aE != aiE && aE != ajE) {
					cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
					cerr << "aE != aiE && aE != ajE" << endl;
					exit(1);
				}

				if (biE == bjE) {
					cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
					cerr << "biE == bjE" << endl;
					exit(1);
				}
				else if (bE != biE && bE != bjE) {
					cerr << "Error: " << __FILE__ << " : " << __LINE__ << endl;
					cerr << "bE != biE && bE != bjE" << endl;
					exit(1);
				}

				rule.emplace_back(new FinderAtom(aE, aiE, ajE, aSortType, aNum, bE, biE, bjE, bSortType, bNum));
			}
		}
		fin.close();
	}
	else {

	}

	Molecule molc;

	if ((opt_r || opt_f) && (argc - optind > 1) || !(opt_r || opt_f) && (argc - optind > 0)) 
	{
		string in_file = argv[argc - 1];
		string out_file = in_file;
		
		const auto pos = out_file.rfind('.');
		if (pos < out_file.size()) {
			out_file.replace(pos, out_file.size() - pos, ".anly");
		}
		else {
			out_file.append(".anly");
		}

		ifstream fin;
		ofstream fout;

		fin.open(in_file.c_str(), ifstream::in);		
		fout.open(out_file.c_str(), ofstream::out);

		fout << setprecision(DATAPRECISION);

		if (opt_f) {

			if (opt_h) {
				fout << setw(BLANK) << left << '#';
				for (int i = 0; i < rule.size(); ++i) {
					ostringstream sout;
					sout << "rule" << i + 1;
					fout << setw(DATAWIDTH) << left << sout.str();
				}
				if (!opt_e)
					fout << endl;
				else
					fout << "Energy" << endl;
			}
		}
		else if (opt_r) {

			if (!opt_h) {
				fout << setw(BLANK) << left << '#';
				for (int i = 0; i < rule.size(); ++i) {
					ostringstream sout;
					sout << "rule" << i + 1;
					fout << setw(DATAWIDTH) << left << sout.str();
				}
				if (opt_e)
					fout << endl;
				else
					fout << "Energy" << endl;
			}
		}
		else {

			if (!opt_h) {
				fout << setw(BLANK) << left << '#';
				for (int iBondType = 0; iBondType < Molecule::nBondtype; ++iBondType) {
					for (int iBond = 0; iBond < Molecule::nBond[iBondType]; ++iBond) {
						ostringstream sout;
						sout << Molecule::bondtype_list[iBondType] << '(' << iBond + 1 << ')';
						fout << setw(DATAWIDTH) << left << sout.str();
					}
				}
				if (opt_e)
					fout << endl;
				else
					fout << "Energy" << endl;
			}
		}

		int tmp;
		while (fin >> tmp) {
			molc.InputEnergy(fin);
			molc.InputX(fin);

			if (opt_f) {
				if (opt_h)
					fout << setw(BLANK) << left << ' ';

				fout << setw(DATAWIDTH) << left << rule[0]->GetBond(molc);
				if (opt_e)
					fout << molc.refEnergy() << endl;
				else
					fout << endl;
			}
			else if (opt_r) {
				if (!opt_h)
					fout << setw(BLANK) << left << ' ';

				for (size_t i = 0; i < rule.size(); ++i) {
					fout << setw(DATAWIDTH) << left << rule[i]->GetBond(molc);
				}
				if (opt_e)
					fout << endl;
				else
					fout << molc.refEnergy() << endl;
			}
			else {
				if (!opt_h)
					fout << setw(BLANK) << left << ' ';

				vector<vector<Molecule::Bond>> bond = molc.refBond();
				for (int iBondtype = 0; iBondtype < Molecule::nBondtype; ++iBondtype) {
					sort(bond[iBondtype].begin(), bond[iBondtype].end());

					for (int iBond = 0; iBond < Molecule::nBond[iBondtype]; ++iBond) {
						fout << setw(DATAWIDTH) << left << bond[iBondtype][iBond].getLen();
					}
				}
				if (opt_e)
					fout << endl;
				else
					fout << molc.refEnergy() << endl;
			}
		}
		fin.close();
		fout.close();
	}
	else {

	cout << setprecision(DATAPRECISION);

	if (opt_f) {

		if (opt_h) {
			cout << setw(BLANK) << left << '#';
			for (int i = 0; i < rule.size(); ++i) {
				ostringstream sout;
				sout << "rule" << i + 1;
				cout << setw(DATAWIDTH) << left << sout.str();
			}
			if (!opt_e)
				cout << endl;
			else
				cout << "Energy" << endl;
		}
	}
	else if (opt_r) {

		if (!opt_h) {
			cout << setw(BLANK) << left << '#';
			for (int i = 0; i < rule.size(); ++i) {
				ostringstream sout;
				sout << "rule" << i + 1;
				cout << setw(DATAWIDTH) << left << sout.str();
			}
			if (opt_e)
				cout << endl;
			else
				cout << "Energy" << endl;
		}
	}
	else {

		if (!opt_h) {
			cout << setw(BLANK) << left << '#';
			for (int iBondType = 0; iBondType < Molecule::nBondtype; ++iBondType) {
				for (int iBond = 0; iBond < Molecule::nBond[iBondType]; ++iBond) {
					ostringstream sout;
					sout << '(' << Molecule::bondtype_list[iBondType] << ')' << iBond + 1;
					cout << setw(DATAWIDTH) << left << sout.str();
				}
			}
			if (opt_e)
				cout << endl;
			else
				cout << "Energy" << endl;
		}
	}

		int tmp;
		while (cin >> tmp) {
			molc.InputEnergy(cin);
			molc.InputX(cin);

			if (opt_f) {
				if (opt_h)
					cout << setw(BLANK) << left << ' ';

				cout << setw(DATAWIDTH) << left << rule[0]->GetBond(molc);
				if (opt_e)
					cout << molc.refEnergy() << endl;
				else
					cout << endl;
			}
			else if (opt_r) {
				if (!opt_h)
					cout << setw(BLANK) << left << ' ';

				for (size_t i = 0; i < rule.size(); ++i) {
					cout << setw(DATAWIDTH) << left << rule[i]->GetBond(molc);
				}
				if (opt_e)
					cout << endl;
				else
					cout << molc.refEnergy() << endl;
			}
			else {
				if (!opt_h)
					cout << setw(BLANK) << left << ' ';

				vector<vector<Molecule::Bond>> bond = molc.refBond();
				for (int iBondtype = 0; iBondtype < Molecule::nBondtype; ++iBondtype) {
					sort(bond[iBondtype].begin(), bond[iBondtype].end());

					for (int iBond = 0; iBond < Molecule::nBond[iBondtype]; ++iBond) {
						cout << setw(DATAWIDTH) << left << bond[iBondtype][iBond].getLen();
					}
				}
				if (opt_e)
					cout << endl;
				else
					cout << molc.refEnergy() << endl;
			}
		}
	}

#ifdef DEBUG_MOLECULE
	debug.close();
#endif // DEBUG_MOLECULE

	return 0;
}