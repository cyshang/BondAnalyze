#ifndef FINDERBOND_H_
#define FINDERBOND_H_

#include <string>
#include "FinderBase.h"

using std::string;

class FinderBond : public FinderBase
{
	int bondType;
	int bondnum;
	bool sortGreat;

public:
	FinderBond(
		const string & iE, const string & jE, const string & sort_type, const int & num
	);
	virtual double GetBond(Molecule & molc);
};

#endif // !FINDERBOND_H_
