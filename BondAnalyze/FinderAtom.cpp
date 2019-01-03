#include <vector>
#include <algorithm>
#include "FinderAtom.h"
#include "Molecule.h"

using std::vector;
using std::sort;

FinderAtom::FinderAtom(
	const string & aE, const string & aiE, const string & ajE, const string & aSort, const int & aNum,
	const string & bE, const string & biE, const string & bjE, const string & bSort, const int & bNum
) :aBondnum(aNum - 1), bBondnum(bNum - 1)
{
	aElem = Molecule::Elem2Num(aE);
	bElem = Molecule::Elem2Num(bE);
	aBondType = Molecule::BondType2Num(Molecule::Elem2Num(aiE), Molecule::Elem2Num(ajE));
	bBondType = Molecule::BondType2Num(Molecule::Elem2Num(biE), Molecule::Elem2Num(bjE));
	aSortGreat = (aSort == "max") ? true : false;
	bSortGreat = (bSort == "max") ? true : false;
}

double FinderAtom::GetBond(Molecule & molc)
{
	vector<Molecule::Bond> aBond = molc.refBond()[aBondType];
	vector<Molecule::Bond> bBond = molc.refBond()[aBondType];

	if (aSortGreat)
		sort(aBond.begin(), aBond.end(), std::greater<Molecule::Bond>());
	else
		sort(aBond.begin(), aBond.end());
	
	if (bSortGreat)
		sort(bBond.begin(), bBond.end(), std::greater<Molecule::Bond>());
	else
		sort(bBond.begin(), bBond.end());

	int aAtom = aBond[aBondnum].getAtom(aElem);
	int bAtom = bBond[bBondnum].getAtom(bElem);

	return (molc.refMatrixR())(aAtom, bAtom);
}