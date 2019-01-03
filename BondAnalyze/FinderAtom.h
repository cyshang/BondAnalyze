#ifndef FINDERATOM_H_
#define FINDERATOM_H_

#include <string>
#include "FinderBase.h"

using std::string;

class FinderAtom : public FinderBase
{
	int aElem;
	int bElem;
	// a atom's iBondType
	int aBondType;
	// b atom's iBondType
	int bBondType;
	int aBondnum;
	int bBondnum;
	bool aSortGreat;
	bool bSortGreat;

public:
	FinderAtom(
		const string & aE, const string & aiE, const string & ajE, const string & aSort, const int & aNum,
		const string & bE, const string & biE, const string & bjE, const string & bSort, const int & bNum
		);
	virtual double GetBond(Molecule & molc);
};

#endif // !FINDERATOM_H_