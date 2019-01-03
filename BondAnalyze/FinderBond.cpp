#include <vector>
#include <algorithm>
#include "FinderBond.h"
#include "Molecule.h"

using std::vector;
using std::sort;

FinderBond::FinderBond(const string & iE, const string & jE, const string & sort_type, const int & num)
	:bondnum(num - 1)
{
	bondType = Molecule::BondType2Num(Molecule::Elem2Num(iE), Molecule::Elem2Num(jE));
	sortGreat = (sort_type == "max") ? true : false;
}

double FinderBond::GetBond(Molecule & molc)
{
	vector<Molecule::Bond> bond = molc.refBond()[bondType];
	
	if (sortGreat)
		sort(bond.begin(), bond.end(), std::greater<Molecule::Bond>());
	else
		sort(bond.begin(), bond.end());

	return bond[bondnum].getLen();
}