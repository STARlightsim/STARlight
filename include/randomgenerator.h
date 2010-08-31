#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

class Randomgenerator
{
	public:
	void SetSeed(unsigned int seed);
	double Rndom(int i=0);
	
	private:
	unsigned int fMt[624];
	int fCount624;
};

#endif //RANDOMGENERATOR_H
	
