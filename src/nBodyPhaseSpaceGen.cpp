//
// calculates n-body phase space (constant matrix element) using various algorithms
//
// the n-body decay is split up into (n - 2) successive 2-body decays
// each 2-body decay is considered in its own center-of-mass frame thereby
// seperating the mass from the (trivial) angular dependence
//
// the event is boosted into the same frame in which the n-body system is
// given
// 
// based on:
// GENBOD (CERNLIB W515), see F. James, "Monte Carlo Phase Space", CERN 68-15 (1968)
// NUPHAZ, see M. M. Block, "Monte Carlo phase space evaluation", Comp. Phys. Commun. 69, 459 (1992)
// S. U. Chung, "Spin Formalism", CERN Yellow Report
// S. U. Chung et. al., "Diffractive Dissociation for COMPASS"
//
// index convention:
// - all vectors have the same size (= number of decay daughters)
// - index i corresponds to the respective value in the (i + 1)-body system: effective mass M, break-up momentum, angles
// - thus some vector elements are not used like breakupMom[0], theta[0], phi[0], ...
//   this overhead is negligible compared to the ease of notation
//
// the following graph illustrates how the n-body decay is decomposed into a sequence of two-body decays
//
// n-body       ...                   3-body                 2-body                 single daughter
//
// m[n - 1]                           m[2]                   m[1]
//  ^                                  ^                      ^
//  |                                  |                      |
//  |                                  |                      |
// M[n - 1] --> ... -->               M[2] -->               M[1] -->               M    [0] = m[0]
// theta[n - 1] ...                   theta[2]               theta[1]               theta[0] = 0 (not used)
// phi  [n - 1] ...                   phi  [2]               phi  [1]               phi  [0] = 0 (not used)
// mSum [n - 1] ...                   mSum [2]               mSum [1]               mSum [0] = m[0]
// = sum_0^(n - 1) m[i]               = m[2] + m[1] + m[0]   = m[1] + m[0]
// breakUpMom[n - 1] ...              breakUpMom[2]          breakUpMom[1]          breakUpMom[0] = 0 (not used)
// = q(M[n - 1], m[n - 1], M[n - 2])  = q(M[2], m[2], M[1])  = q(M[1], m[1], m[0])
// 


#include <algorithm>

#include "nBodyPhaseSpaceGen.h"


using namespace std;
using namespace StarlightConstants;


nBodyPhaseSpaceGen::nBodyPhaseSpaceGen()
	: _n                (0),
	  _norm             (0),
	  _weight           (0),
	  _maxWeightObserved(0),
	  _maxWeight        (0)
{ }


nBodyPhaseSpaceGen::~nBodyPhaseSpaceGen()
{ }


// sets decay constants and prepares internal variables
bool
nBodyPhaseSpaceGen::setDecay(const vector<double>& daughterMasses)  // array of daughter particle masses
{
	_n = daughterMasses.size();
	if (_n < 2) {
		cout << "error: number of daughters = " << _n << " does not make sense." << endl;
		return false;
	}
	// copy daughter masses
	_m.clear();
	_m = daughterMasses;
	// prepare effective mass vector
	_M.clear();
	_M.resize(_n, 0);
	_M[0] = _m[0];
	// prepare angle vectors
	_cosTheta.clear();
	_cosTheta.resize(_n, 0);
	_phi.clear();
	_phi.resize(_n, 0);
	// calculate daughter mass sums
	_mSum.clear();
	_mSum.resize(_n, 0);
	_mSum[0] = _m[0];
	for (unsigned int i = 1; i < _n; ++i)
		_mSum[i] = _mSum[i - 1] + _m[i];
	// prepare breakup momentum vector
	_breakupMom.clear();
	_breakupMom.resize(_n, 0);
	// prepare vector for daughter Lorentz vectors
	_daughters.clear();
	_daughters.resize(_n, LorentzVector(0, 0, 0, 0));
	// calculate normalization
	_norm = 1 / (2 * pow(twoPi, 2 * (int)_n - 3) * factorial(_n - 2));
	resetMaxWeightObserved();
	return true;
}


// set decay constants and prepare internal variables
bool
nBodyPhaseSpaceGen::setDecay(const unsigned int nmbOfDaughters,  // number of daughter particles
                             const double*      daughterMasses)  // array of daughter particle masses
{
	vector <double> m;
	m.resize(nmbOfDaughters, 0);
	for (unsigned int i = 0; i < nmbOfDaughters; ++i)
		m[i] = daughterMasses[i];
	return setDecay(m);
}


// generates event with certain n-body mass and momentum and returns event weigth
// general purpose function
double
nBodyPhaseSpaceGen::generateDecay(const LorentzVector& nBody)  // Lorentz vector of n-body system in lab frame
{
	const double nBodyMass = nBody.M();
	if (_n < 2) {
		cout << "warning: number of daughter particles = " << _n << " is smaller than 2. weight is set to 0." << endl;
		_weight = 0;
	} else if (nBodyMass < _mSum[_n - 1]) {
		cout << "warning: n-body mass = " << nBodyMass << " is smaller than sum of daughter masses = "
		     << _mSum[_n - 1] << ". weight is set to 0." << endl;
		_weight = 0;
	} else {
		pickMasses(nBodyMass);
		calcWeight();
		pickAngles();
		calcEventKinematics(nBody);
	}
	return _weight;
}


// generates full event with certain n-body mass and momentum only, when event is accepted (return value = true)
// this function is more efficient, if only weighted evens are needed
bool
nBodyPhaseSpaceGen::generateDecayAccepted(const LorentzVector& nBody,      // Lorentz vector of n-body system in lab frame
                                          const double         maxWeight)  // if positive, given value is used as maximum weight, otherwise _maxWeight
{
	const double nBodyMass = nBody.M();
	if (_n < 2) {
		cout << "warning: number of daughter particles = " << _n << " is smaller than 2. no event generated." << endl;
		return false;
	} else if (nBodyMass < _mSum[_n - 1]) {
		cout << "warning: n-body mass = " << nBodyMass << " is smaller than sum of daughter masses = "
		     << _mSum[_n - 1] << ". no event generated." << endl;
		return false;
	}
	pickMasses(nBodyMass);
	calcWeight();
	if (!eventAccepted(maxWeight))
		return false;
	pickAngles();
	calcEventKinematics(nBody);
	return true;
}


// randomly choses the (n - 2) effective masses of the respective (i + 1)-body systems
void
nBodyPhaseSpaceGen::pickMasses(const double nBodyMass)  // total energy of the system in its RF
{
	_M[_n - 1] = nBodyMass;
	// create vector of sorted random values
	vector<double> r(_n - 2, 0);  // (n - 2) values needed for 2- through (n - 1)-body systems
	for (unsigned int i = 0; i < (_n - 2); ++i)
		r[i] = random();
	sort(r.begin(), r.end());
	// set effective masses of (intermediate) two-body decays
	const double massInterval = nBodyMass - _mSum[_n - 1];  // kinematically allowed mass interval
	for (unsigned int i = 1; i < (_n - 1); ++i)             // loop over intermediate 2- to (n - 1)-bodies
		_M[i] = _mSum[i] + r[i - 1] * massInterval;           // _mSum[i] is minimum effective mass
}


// computes event weight (= integrand value) and breakup momenta
// uses vector of intermediate two-body masses prepared by pickMasses()
double
nBodyPhaseSpaceGen::calcWeight()
{
	for (unsigned int i = 1; i < _n; ++i)  // loop over 2- to n-bodies
		_breakupMom[i] = breakupMomentum(_M[i], _M[i - 1], _m[i]);
	double momProd = 1;                    // product of breakup momenta
	for (unsigned int i = 1; i < _n; ++i)  // loop over 2- to n-bodies
		momProd *= _breakupMom[i];
	const double massInterval = _M[_n - 1] - _mSum[_n - 1];  // kinematically allowed mass interval
	_weight = _norm * pow(massInterval, (int)_n - 2) * momProd / _M[_n - 1];
	if (_weight > _maxWeightObserved)
		_maxWeightObserved = _weight;
	if (isnan(_weight))
		cout << "warning: weight = " << _weight << endl;
	return _weight;
}


// calculates complete event from the effective masses of the (i + 1)-body
// systems, the Lorentz vector of the decaying system, and the decay angles
// uses the break-up momenta calculated by calcWeight()
void
nBodyPhaseSpaceGen::calcEventKinematics(const LorentzVector& nBody)  // Lorentz vector of n-body system in lab frame
{
	// build event starting in n-body RF going down to 2-body RF
	// is more efficicient than Raubold-Lynch method, since it requitres only one rotation and boost per daughter
	LorentzVector P = nBody;  // Lorentz of (i + 1)-body system in lab frame
	for (unsigned int i = _n - 1; i >= 1; --i) {  // loop from n-body down to 2-body
		// construct Lorentz vector of daughter _m[i] in (i + 1)-body RF
		const double   sinTheta = sqrt(1 - _cosTheta[i] * _cosTheta[i]);
		const double   pT       = _breakupMom[i] * sinTheta;
		LorentzVector& daughter = _daughters[i];
		daughter.SetPxPyPzE(pT * cos(_phi[i]),
		                    pT * sin(_phi[i]),
		                    _breakupMom[i] * _cosTheta[i],
		                    sqrt(_m[i] * _m[i] + _breakupMom[i] * _breakupMom[i]));
		// boost daughter into lab frame
		daughter.Boost(P.BoostVector());
		// calculate Lorentz vector of i-body system in lab frame
		P -= daughter;
	}
	// set last daughter
	_daughters[0] = P;
}


// calculates maximum weight for given n-body mass
double
nBodyPhaseSpaceGen::estimateMaxWeight(const double       nBodyMass,        // sic!
                                      const unsigned int nmbOfIterations)  // number of generated events
{
	double maxWeight = 0;
	for (unsigned int i = 0; i < nmbOfIterations; ++i) {
		pickMasses(nBodyMass);
		calcWeight();
		maxWeight = max(_weight, maxWeight);
	}
	return maxWeight;
}


ostream&
nBodyPhaseSpaceGen::print(ostream& out) const
{
	out << "nBodyPhaseSpaceGen parameters:" << endl
	    << "    number of daughter particles ............... " << _n                 << endl
	    << "    masses of the daughter particles ........... " << _m                 << endl
	    << "    sums of daughter particle masses ........... " << _mSum              << endl
	    << "    effective masses of (i + 1)-body systems ... " << _M                 << endl
	    << "    cos(polar angle) in (i + 1)-body systems ... " << _cosTheta          << endl
	    << "    azimuth in (i + 1)-body systems ............ " << _phi               << endl
	    << "    breakup momenta in (i + 1)-body systems .... " << _breakupMom        << endl
	    << "    normalization value ........................ " << _norm              << endl
	    << "    weight of generated event .................. " << _weight            << endl
	    << "    maximum weight used in hit-miss MC ......... " << _maxWeight         << endl
	    << "    maximum weight since instantiation ......... " << _maxWeightObserved << endl
	    << "    daughter four-momenta:" << endl;
	for (unsigned int i = 0; i < _n; ++i)
		out << "        daughter " << i << ": " << _daughters[i] << endl;
	return out;
}
