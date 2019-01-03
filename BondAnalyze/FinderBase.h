#ifndef FINDERBASE_H_
#define FINDERBASE_H_

class Molecule;

class FinderBase
{
public:
	virtual double GetBond(Molecule & molc) = 0;
	virtual ~FinderBase() {}
};

#endif // !FINDERBASE_H_

