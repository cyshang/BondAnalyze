#include "Molecule.h"
#include <sstream>
#include <stdexcept>

using namespace Eigen;
using std::vector;
using std::string;
using std::istream;
using std::ostream;
using std::istringstream;
using std::ostringstream;
using std::endl;
using std::cerr;

// ---------- define static member variable ----------
string Molecule::name;

int Molecule::nElem;
vector<string> Molecule::elem_list;
Molecule::str2int Molecule::elem2num;
Molecule::int2str Molecule::num2elem;

int Molecule::totAtom;
vector<int> Molecule::nAtom;
vector<int> Molecule::atom_list;
vector<vector<int>> Molecule::atomTravlist;

int Molecule::totBond;
int Molecule::nBondtype;
vector<Molecule::BondType> Molecule::bondtype_list;
Molecule::bdtype2int Molecule::bdtype2num;
Molecule::int2bdtype Molecule::num2bdtype;
vector<int> Molecule::nBond;
vector<vector<Molecule::Array2>> Molecule::bondTravlist;

bool Molecule::ifString = false;
bool Molecule::ifVectorR = false;
bool Molecule::ifMatrixR = false;
bool Molecule::ifBond = false;

// =============== construct =============== 

Molecule::Molecule()
	:energy_str(), atom_str(),
	Energy(0.0),
	X(3, totAtom),
	vectorR(),
	bond()
{
	if (ifString)
		atom_str.resize(totAtom);

	if (ifVectorR)
		vectorR.resize(totAtom * (totAtom - 1) / 2);

	if (ifMatrixR)
		matrixR.resize(totAtom, totAtom);

	if (ifBond) {
		bond.resize(nBondtype);
		for (int iBondtype = 0; iBondtype < nBondtype; ++iBondtype) {
			bond[iBondtype].resize(nBond[iBondtype]);
		}
	}
}
// =========================================

// ======================================================
// =============== static member function ===============
// ======================================================

void Molecule::InputInfo(istream & fin)
{
	string line;
	while (std::getline(fin, line)) {
		if (line[0] == '#' || line[0] == '\n' || line[0] == ' ')
			continue;
		
		for (size_t i = 0; i < line.size(); ++i) {
			if (line[i] == '=' || line[i] == '(' || line[i] == ')' || line[i]=='\'' || line[i]=='"')
				line[i] = ' ';
		}
		
		istringstream sin(line);
		string key_word;

		sin >> key_word;
		if (key_word == "molecule") {
			sin >> name;
		}
		else if (key_word == "elem_num") {
			sin >> nElem;
		}
		else if (key_word == "elem_list") {
			elem_list.resize(nElem);

			for (size_t i = 0; i < elem_list.size(); ++i) {
				sin >> elem_list[i];
				elem2num[elem_list[i]] = i;
				num2elem[i] = elem_list[i];
			}
		}
		else if (key_word == "atom_num") {
			sin >> totAtom;
		}
		else if (key_word == "atom_list") {
			atom_list.resize(totAtom);

			nAtom.resize(nElem);
			for (size_t i = 0; i < nAtom.size(); ++i) {
				nAtom[i] = 0;
			}

			string elem;
			for (size_t i = 0; i < atom_list.size(); ++i) {
				sin >> elem;
				atom_list[i] = elem2num[elem];
				nAtom[atom_list[i]]++;
			}
		}
	}

	atomTravlist.resize(nElem, vector<int>());
	for (int iAtom = 0; iAtom < totAtom; ++iAtom) {
		atomTravlist[atom_list[iAtom]].push_back(iAtom);
	}

#ifdef DEBUG_MOLECULE

	debug << "name: " << name << endl;
	debug << "nElem: " << nElem << endl;
	debug << "elem_list: ";
	for (auto & i : elem_list) debug << i << ' ';
	debug << endl;
	for (auto & i : elem_list) {
		debug << elem2num[i] << "<->" << num2elem[elem2num[i]] << endl;
	}
	debug << "totAtom: " << totAtom << endl;
	debug << "nAtom: ";
	for (auto & i : nAtom) debug << i << ' ';
	debug << endl;
	debug << "atom_list: ";
	for (auto & i : atom_list) debug << i << ' ';
	debug << endl;
	debug << "atomTravlist:" << endl;
	for (size_t i = 0; i < atomTravlist.size(); ++i) {
		debug << num2elem[i] << ": ";
		for (const auto & j : atomTravlist[i])
			debug << j << ' ';
		debug << endl;
	}
#endif // DEBUG_MOLECULE

	if (ifBond)
		BondInfo();
}

void Molecule::BondInfo()
{
	totBond = (totAtom * (totAtom - 1)) / 2;

	nBondtype = 0;
	bondtype_list.clear();
	for (int iE = 0; iE < nElem; ++iE) {
		if (nAtom[iE] > 1) {
			const BondType tmp_type(iE, iE);
			bdtype2num[tmp_type] = nBondtype;
			num2bdtype[nBondtype] = tmp_type;
			bondtype_list.push_back(BondType(iE, iE));
			nBondtype++;
		}
		for (int jE = iE + 1; jE < nElem; ++jE) {
			const BondType tmp_type(iE, jE);
			bdtype2num[tmp_type] = nBondtype;
			num2bdtype[nBondtype] = tmp_type;
			bondtype_list.push_back(BondType(iE, jE));
			nBondtype++;
		}
	}

	bondTravlist.clear();
	nBond.clear();
	int iBondtype = 0;
	for (int iE = 0; iE < nElem; ++iE) {
		if (nAtom[iE] > 1) {
			nBond.push_back((nAtom[iE] * (nAtom[iE] - 1)) / 2);
			bondTravlist.push_back(vector<Array2>());
			for (int i = 0; i < nAtom[iE] - 1; ++i) {
				for (int j = i + 1; j < nAtom[iE]; ++j) {
					bondTravlist[iBondtype].push_back(Array2(atomTravlist[iE][i], atomTravlist[iE][j]));
				}
			}
			iBondtype++;
		}
		for (int jE = iE + 1; jE < nElem; ++jE) {
			nBond.push_back(nAtom[iE] * nAtom[jE]);
			bondTravlist.push_back(vector<Array2>());
			for (int i = 0; i < nAtom[iE]; ++i) {
				for (int j = 0; j < nAtom[jE]; ++j) {
					bondTravlist[iBondtype].push_back(Array2(atomTravlist[iE][i], atomTravlist[jE][j]));
				}
			}
			iBondtype++;
		}
	}

#ifdef DEBUG_MOLECULE

	debug << "totBond: " << totBond << endl;
	debug << "nBondtype: " << nBondtype << endl;
	debug << "num to bondtype:" << endl;
	for (int i = 0; i < nBondtype; ++i) {
		debug << i << " ==> " << Num2BondType(i) << endl;
	}
	debug << "bondtype to num: " << endl;
	for (int i = 0; i < nBondtype; ++i) {
		debug << Num2Elem(bondtype_list[i].iElem) << '-' << Num2Elem(bondtype_list[i].jElem);
		debug << " ==> " << BondType2Num(bondtype_list[i].iElem, bondtype_list[i].jElem) << endl;
	}
	debug << "nBond: ";
	for (const auto & i : nBond) debug << i << ' ';
	debug << endl;
	debug << "bondTravlist:" << endl;
	for (const auto & i : bondTravlist) {
		for (const auto & j : i) {
			debug << '(' << j.iAtom << ", " << j.jAtom << ") ";
		}
		debug << endl;
	}

#endif // DEBUG_MOLECULE

}

// ======================================================
// =================== input function ===================
// ======================================================

void Molecule::InputX(std::istream & fin)
{
	string elem;
	for (int i = 0; i < totAtom; ++i) {
		fin >> elem;

		for (int j = 0; j < 3; ++j) {
			fin >> X(j, i);
		}
	}

	if (ifVectorR)
		CalcVectorR();
	if (ifMatrixR)
		CalcMatrixR();
	if (ifBond)
		CalcBond();
}

std::istream & operator >> (std::istream & fin, Molecule & m)
{
	using std::getline;
	
	string line;
	string elem;
	if (getline(fin, line)) {

		getline(fin, m.energy_str);
		for (int i = 0; i < m.totAtom; ++i) {
			getline(fin, m.atom_str[i]);

			istringstream sin(m.atom_str[i]);

			sin >> elem;
			for (int j = 0; j < 3; ++j) {
				sin >> m.X(j, i);
			}
		}
	}

	if (m.ifVectorR)
		m.CalcVectorR();
	if (m.ifMatrixR)
		m.CalcMatrixR();
	if (m.ifBond)
		m.CalcBond();

	return fin;
}

void Molecule::CalcBond()
{
	for (int iBondtype = 0; iBondtype < nBondtype; ++iBondtype) {
		for (int iBond = 0; iBond < nBond[iBondtype]; ++iBond) {
			const auto & ij = bondTravlist[iBondtype][iBond];
			bond[iBondtype][iBond].assign((X.col(ij.iAtom) - X.col(ij.jAtom)).norm(), ij.iAtom, ij.jAtom);
		}
	}

#ifdef DEBUG_MOLECULE
	debug << "bond:" << endl;
	for (int iBondtype = 0; iBondtype < nBondtype; ++iBondtype) {
		debug << num2bdtype[iBondtype] << " : ";
		for (int iBond = 0; iBond < nBond[iBondtype]; ++iBond) {
			debug << bond[iBondtype][iBond].len << ", ";
		}
		debug << endl;
	}
	debug << endl;
#endif // DEBUG_MOLECULE

}

void Molecule::CalcVectorR()
{
	int pos = 0;
	for (int i = 0; i < totAtom - 1; ++i) {
		for (int j = i + 1; j < totAtom; ++j) {
			vectorR(pos++) = (X.col(i) - X.col(j)).norm();
		}
	}

#ifdef DEBUG_MOLECULE
	debug << "vectorR:" << endl;
	debug << vectorR << endl;
	debug << endl;
#endif // DEBUG_MOLECULE
}

void Molecule::CalcMatrixR()
{
	matrixR.setZero();
	for (int i = 0; i < totAtom - 1; ++i) {
		for (int j = i + 1; j < totAtom; ++j) {
			matrixR(i, j) = (X.col(i) - X.col(j)).norm();
			matrixR(j, i) = matrixR(i, j);
		}
	}

#ifdef DEBUG_MOLECULE
	debug << "matrixR:" << endl;
	debug << matrixR << endl;
	debug << endl;
#endif // DEBUG_MOLECULE

}

ostream & operator << (ostream & fout, const Molecule & m)
{
	fout << m.totAtom << endl;
	fout << m.energy_str << endl;
	for (size_t i = 0; i < m.atom_str.size(); ++i) {
		fout << m.atom_str[i] << endl;
	}

	return fout;
}

int Molecule::Elem2Num(const std::string & elem)
{
	try {
		return elem2num.at(elem);
	}
	catch (std::out_of_range & e) {
		cerr << "Error: " << __FILE__ << ": " << __LINE__ << endl;
		cerr << "Element: " << elem << " doesn't exist!" << endl;
		exit(1);
	}
}

std::string Molecule::Num2Elem(const int & iE)
{
	try {
		return num2elem.at(iE);
	}
	catch(std::out_of_range & e){
		cerr << "Error: " << __FILE__ << ": " << __LINE__ << endl;
		cerr << "Element num: " << iE << " out of range!" << endl;
		exit(1);
	}
	
}

int Molecule::BondType2Num(const int & iE, const int & jE)
{
	try {
		return bdtype2num.at(BondType(iE, jE));
	}
	catch (std::out_of_range & e) {
		cerr << "Error: " << __FILE__ << ": " << __LINE__ << endl;
		cerr << "The BondType: " << BondType(iE, jE) << " doesn't exist!" << endl;
		exit(1);
	}
}

Molecule::BondType Molecule::Num2BondType(const int & iBondtype)
{
	try {
		return num2bdtype.at(iBondtype);
	}
	catch (std::out_of_range & e) {
		cerr << "Error: " << __FILE__ << ": " << __LINE__ << endl;
		cerr << "iBondType: " << iBondtype << " out of range!" << endl;
		exit(1);
	}	
}