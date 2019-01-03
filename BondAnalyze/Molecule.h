#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Core>

extern std::ofstream debug;

#define DEBUG_MOLECULE

class Molecule
{
	// ===========================================================
	// ========================= public ==========================
	// ===========================================================
public:
	// ---------- typedef ----------
	class Bond;
	class BondType;
	struct Array2;
	typedef std::map<std::string, int> str2int;
	typedef std::map<int, std::string> int2str;
	typedef std::map<BondType, int> bdtype2int;
	typedef std::map<int, BondType> int2bdtype;

	// ---------- construct ----------
	Molecule();

	// =============== variable usage =============== 

	// use string data
	inline static void usingString() { ifString = true; }
	// use vectorR
	inline static void usingVectorR() { ifVectorR = true; }
	// use matrixR
	inline static void usingMatrixR() { ifMatrixR = true; }
	// use bond
	inline static void usingBond() { ifBond = true; }

	// =============== input function ===============

	// input static info
	static void InputInfo(std::istream &);
	// input energy
	inline bool InputEnergy(std::istream & fin) { return (fin >> Energy) ? true : false; }
	// input X
	void InputX(std::istream &);
	// input string data
	friend std::istream & operator >> (std::istream &, Molecule &);

	// =============== calculate function ===============

	// calculate vectorR
	void CalcVectorR();
	// calculate vectorR, result stored in argument
	template<typename Derived>
	void CalcVectorR(const Eigen::MatrixBase<Derived> & R);
	// calculate matrixR
	void CalcMatrixR();
	// calculate Euclidean distance between two vectorR
	inline double operator - (const Molecule & m) const { return (vectorR - m.vectorR).norm(); }
	// calculate bond
	void CalcBond();

	// =============== data pointer ===============

	// return X.data()
	inline double* X_ptr() { return X.data(); }
	// return vectorR.data()
	inline double* vectorR_ptr() { return vectorR.data(); }
	// return matrixR.data()
	inline double* matrixR_prt() { return matrixR.data(); }
	
	// =============== data reference ===============

	// return reference of energy
	inline double & refEnergy() { return Energy; }
	// return reference of X
	inline Eigen::MatrixXd & refX() { return X; }
	// return reference of vectorR
	inline Eigen::VectorXd & refVectorR() { return vectorR; }
	// return reference of matrixR
	inline Eigen::MatrixXd & refMatrixR() { return matrixR; }
	// return reference of bond
	inline std::vector<std::vector<Bond>> & refBond() { return bond; }
	// convert BondType to iBondtype;
	static int BondType2Num(const int & iE, const int & jE);
	static BondType Num2BondType(const int & iBondtype);
	static int Elem2Num(const std::string &);
	static std::string Num2Elem(const int &);

	// =============== output ===============

	// print string data of Molecule
	friend std::ostream & operator << (std::ostream &, const Molecule &);

	// =============== molecule description ===============

	// name of molecule
	static std::string name;
	// number of elements
	static int nElem;
	// a list of elements
	static std::vector<std::string> elem_list;
	// convert element string to element id
private:
	static str2int elem2num;
	// convert element id to element string
	static int2str num2elem;
	// total number of atoms
public:
	static int totAtom;
	// number of atoms for each element
	static std::vector<int> nAtom;
	// a list of atom <-> element id
	static std::vector<int> atom_list;
	// a list of atom id, sort by element
	static std::vector<std::vector<int>> atomTravlist;

	// =============== bond description ===============

	// total number of bond
	static int totBond;
	// total number of bond types
	static int nBondtype;
	static std::vector<BondType> bondtype_list;
private:
	// convert BondType to iBondType
	static bdtype2int bdtype2num;
	// convert iBondType to BondType
	static int2bdtype num2bdtype;
public:
	// number of bonds for each BondType
	static std::vector<int> nBond;
	// bond traversal list, a list of atom ids for each BondType, used to calculate bond data
	static std::vector<std::vector<Array2>> bondTravlist;

	// ===========================================================
	// ========================= private =========================
	// ===========================================================
private:
	// =============== molecule string data ===============

	// the line of energy
	std::string energy_str;
	// lines of atoms
	std::vector<std::string> atom_str;

	// =============== molcule data ===============
	double Energy;
	Eigen::MatrixXd X;
	Eigen::VectorXd vectorR;
	std::vector<std::vector<Bond>> bond;

	Eigen::MatrixXd matrixR;
	//Eigen::MatrixXd matrixR2;
	//std::vector<Eigen::MatrixXd> cos0;

	// =============== variable usage ===============
	static bool ifString;
	static bool ifVectorR;
	static bool ifMatrixR;
	static bool ifBond;
	static void BondInfo();
};

// ========== template functions ==========

template<typename Derived>
void Molecule::CalcVectorR(const Eigen::MatrixBase<Derived> & R)
{
	int pos = 0;
	for (int i = 0; i < totAtom - 1; ++i) {
		for (int j = i + 1; j < totAtom; ++j) {
			const_cast<Eigen::MatrixBase<Derived>&>(R)(pos++) = (X.col(i) - X.col(j)).norm();
		}
	}
}

// ========================================

// ============================================================
// ========================= BondType =========================
// ============================================================

class Molecule::BondType
{
public:

	int iElem;
	int jElem;

	BondType() {}
	BondType(const int & iE, const int & jE) {
		iElem = (iE < jE) ? iE : jE;
		jElem = (iE < jE) ? jE : iE;
	}
	BondType(const BondType & t) {
		iElem = t.iElem;
		jElem = t.jElem;
	}

	inline BondType & operator = (const BondType & t) {
		iElem = t.iElem;
		jElem = t.jElem;
		return *this;
	}

	inline bool operator == (const BondType type) { return iElem == type.iElem && jElem == type.jElem; }
	inline bool operator != (const BondType type) { return iElem != type.iElem || jElem != type.jElem; }
	friend inline bool operator < (const BondType a, const BondType b) {
		return ((a.iElem < b.iElem) ? true : (a.iElem == b.iElem && a.jElem < b.jElem));
	}
	friend inline std::ostream & operator << (std::ostream & os, const BondType & t) { 
		os << Molecule::Num2Elem(t.iElem) << '-' << Molecule::Num2Elem(t.jElem);
		return os;
	}
};
// ============================================================

struct Molecule::Array2
{
	int iAtom;
	int jAtom;

	Array2(const int & i, const int & j) :iAtom(i), jAtom(j) {}
};


// ============================================================
// ========================= Bond =============================
// ============================================================
class Molecule::Bond
{
	friend Molecule;

public:
	Bond() :len(0.0), iAtom(0), jAtom(0),iElem(0),jElem(0) {}
	Bond(const double & _len, const int & i, const int & j) :len(_len), iAtom(i), jAtom(j) 
	{
		iElem = Molecule::atom_list[iAtom];
		jElem = Molecule::atom_list[jAtom];
	}

	inline void assign(const double & _len, const int & i, const int & j) {
		len = _len;
		iAtom = i;
		jAtom = j;
	}

	inline double getLen() { return len; }

	inline int getAtom(const int & elem) {
		return ((elem == iElem) ? iAtom : jAtom);
	}

	inline int getOtherAtom(const int & elem) {
		return ((elem == iElem) ? jAtom : iAtom);
	}

	friend inline bool operator < (const Bond & a, const Bond & b) {
		return a.len < b.len;
	}

	friend inline bool operator > (const Bond & a, const Bond & b) {
		return a.len > b.len;
	}

	friend inline std::ostream & operator << (std::ostream & os, const Bond & a) {
		os << a.iAtom << ", " << a.jAtom << ", " << a.len;
		return os;
	}

private:
	double len;
	int iAtom, jAtom;
	int iElem, jElem;
};
// ============================================================

#endif // !MOLECULE_H_