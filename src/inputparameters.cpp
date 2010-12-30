///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "reportingUtils.h"
#include "starlightconstants.h"
#include "inputparameters.h"


using namespace std;


//______________________________________________________________________________
inputParameters::inputParameters()
	: _configFileName        ("slight.in"),
	  _beam1Z                (0),
	  _beam1A                (0),
	  _beam2Z                (0),
	  _beam2A                (0),
	  _beamLorentzGamma      (0),
	  _maxW                  (0),
	  _minW                  (0),
	  _nmbWBins              (0),
	  _maxRapidity           (0),
	  _nmbRapidityBins       (0),
	  _productionMode        (0),
	  _nmbEventsTot          (0),
	  _prodParticleId        (0),
	  _randomSeed            (0),
	  _outputFormat          (0),
	  _beamBreakupMode       (0),
	  _interferenceEnabled   (0),
	  _interferenceStrength  (0),
	  _coherentProduction    (0),
	  _incoherentFactor      (0),
	  _bford                 (0),
	  _maxPtInterference     (0),
	  _nmbPtBinsInterference (0),
	  _ptBinWidthInterference(0),
	  _protonEnergy          (0)
{ }


//______________________________________________________________________________
inputParameters::~inputParameters()
{ }


//______________________________________________________________________________
bool
inputParameters::init(const string& configFileName)
{
	// open config file
	_configFileName = configFileName;
	ifstream configFile(_configFileName.c_str());
	if ((!configFile) || (!configFile.good())) {
		printWarn << "cannot read parameters from file '" << _configFileName << "'" << endl;
		return false;
	}

	// read config file
	double Wmaxinput = 0;
	double Wmininput = 0;
	if (configFile
	    >> _beam1Z >> _beam1A
	    >> _beam2Z>> _beam2A
	    >> _beamLorentzGamma
	    >> Wmaxinput>> Wmininput >> _nmbWBins
	    >> _maxRapidity >> _nmbRapidityBins
	    >> _productionMode
	    >> _nmbEventsTot
	    >> _prodParticleId
	    >> _randomSeed
	    >> _outputFormat
	    >> _beamBreakupMode
	    >> _interferenceEnabled >> _interferenceStrength
	    >> _coherentProduction >> _incoherentFactor
	    >> _bford
	    >> _maxPtInterference
	    >> _nmbPtBinsInterference)
		printInfo << "successfully read input parameters from '" << _configFileName << "'" << endl;
	else {
		printWarn << "problems reading input parameters from '" << _configFileName << "'" << endl
		          << *this;
		return false;
	}
	configFile.close();

	_ptBinWidthInterference = _maxPtInterference / _nmbPtBinsInterference;
	_protonEnergy           = _beamLorentzGamma * starlightConstants::protonMass;

	// define interaction type
	switch (_productionMode) {
	case 1:
		_interactionType = starlightConstants::PHOTONPHOTON;
		break;
	case 2:
		_interactionType = starlightConstants::PHOTONPOMERONNARROW;
		break;
	case 3:
		_interactionType = starlightConstants::PHOTONPOMERONWIDE;
		break;
	default:
		printWarn << "unknown production mode " << _productionMode << endl;
		return false;
	}
  
	//Trying to define the proper Wmins and Wmaxs. a TEMPORARY fix....Better solution=??
	double mass         = 0;
	double Wmin_default = 0; //This is the dault for _minW, unless defined later
	double width        = 0;
	if(_prodParticleId==11){
		//electron
		Wmin_default = 0.01; //_minW is settable in constants--default is 0.01; _minW up to 0.15 is safe for Summer 2000 triggering for e+e- pairs
		_maxW=Wmaxinput;
		_particleType=starlightConstants::ELECTRON;
		_decayType=starlightConstants::LEPTONPAIR;
	}
	else if(_prodParticleId==13){
		//muon
		Wmin_default   =  2.*starlightConstants::mmu;
		_maxW=Wmaxinput;
		_particleType=starlightConstants::MUON;
		_decayType=starlightConstants::LEPTONPAIR;
	}
	else if(_prodParticleId==15){
		//tauon
		Wmin_default   =  2*starlightConstants::mtau;
		_maxW=Wmaxinput;
		_particleType=starlightConstants::TAUON;
		_decayType=starlightConstants::LEPTONPAIR;
	}
	else if(_prodParticleId==115){
		//a2(1320)
		_maxW=Wmaxinput;
		_particleType=starlightConstants::A2;
		_decayType=starlightConstants::SINGLEMESON;
	}
	else if(_prodParticleId==221){
		//eta
		_maxW=Wmaxinput;
		_particleType=starlightConstants::ETA;
		_decayType=starlightConstants::SINGLEMESON;
	}
	else if(_prodParticleId==225){
		//f2(1270)
		_maxW=Wmaxinput;
		_particleType=starlightConstants::F2;
		_decayType=starlightConstants::SINGLEMESON;
	}
	else if(_prodParticleId==331){
		//eta'(958)
		_maxW=Wmaxinput;
		_particleType=starlightConstants::ETAPRIME;
		_decayType=starlightConstants::SINGLEMESON;
	}
	else if(_prodParticleId==335){
		//f2'(1525)
		_maxW=Wmaxinput;
		_particleType=starlightConstants::F2PRIME;
		_decayType=starlightConstants::SINGLEMESON;
	}
	else if(_prodParticleId==441){
		//eta-c
		_maxW=Wmaxinput;
		_particleType=starlightConstants::ETAC;
		_decayType=starlightConstants::SINGLEMESON;
	}
	else if(_prodParticleId==9010221){
		//f0(980)--was orginally called 10221? updated to standard number
		_maxW=Wmaxinput;
		_particleType=starlightConstants::F0;
		_decayType=starlightConstants::SINGLEMESON;
	}
	else if(_prodParticleId==33){
		//Z"/Z03
		_maxW=Wmaxinput;
		_particleType=starlightConstants::ZOVERZ03;
		_decayType=starlightConstants::SINGLEMESON;
	}
	else if(_prodParticleId==113)   {
		mass        =  0.7685;                 //rho(770)0
		width       =  0.1507;
		Wmin_default=  2.*starlightConstants::mpi;
		_maxW        =  mass+5.*width;
		_particleType=starlightConstants::RHO;
		_decayType=starlightConstants::WIDEVMDEFAULT;
	}
	else if(_prodParticleId==913)   {
		mass        =  0.7685;                      //rho(770)0 direct pi+pi- decay, interference given by ZEUS
		width       =  0.1507;
		Wmin_default=  2.*starlightConstants::mpi;
		_maxW        =  mass+5.*width;    //Use the same (ZEUS) 1.5GeV max mass
		_particleType=starlightConstants::RHOZEUS;
		_decayType=starlightConstants::WIDEVMDEFAULT;
	}
	else if (_prodParticleId == 999) {  // pi+pi-pi+pi- decay
		mass         = 1.350;
		width        = 0.360;
		Wmin_default = 4 * starlightConstants::mpi;
		_maxW         = 3;
		_particleType      = starlightConstants::FOURPRONG;
		_decayType    = starlightConstants::WIDEVMDEFAULT;
	}
	else if(_prodParticleId==223){
		mass          =  0.78194;          //omega(782)
		width         =  0.00843;
		_maxW          =  mass+5.*width;
		Wmin_default  =  mass - 5.*width;
		_particleType=starlightConstants::OMEGA;
		_decayType=starlightConstants::NARROWVMDEFAULT;//will probably be shifted to 3body decay
	}
	else if(_prodParticleId==333) {
		mass          =  1.019413;     //phi(1020)
		width         =  0.00443;
		Wmin_default  =  2*starlightConstants::mK;
		_maxW          =  mass+5.*width;
		_particleType=starlightConstants::PHI;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==443){
		mass          =  3.09692;  // JN  3.09688;      //J/psi
		width         =  0.000091; //JN  0.000087;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::JPSI;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==444){
		mass          =  3.686093;      //J/psi
		width         =  0.000337;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::JPSI2S;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==444011){
		mass          =  3.686093;      //J/psi
		width         =  0.000337;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::JPSI2S_ee;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==444013){
		mass          =  3.686093;      //J/psi
		width         =  0.000337;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::JPSI2S_mumu;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==443011){
		mass          =  3.09692;  // JN  3.09688;      //J/psi
		width         =  0.000091; //JN  0.000087;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::JPSI_ee;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==443013){
		mass          =  3.09692;  // JN  3.09688;      //J/psi
		width         =  0.000091; //JN  0.000087;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::JPSI_mumu;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==553){
		mass          =  9.46030;      //Upsilon
		width         =  0.00005402;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==553011){
		mass          =  9.46030;      //Upsilon
		width         =  0.00005402;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON_ee;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==553013){
		mass          =  9.46030;      //Upsilon
		width         =  0.00005402;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON_mumu;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==554){
		mass          =  10.02326;      //Upsilon(2S)
		width         =  0.00003198;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON2S;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==554011){
		mass          =  10.02326;      //Upsilon(2S)
		width         =  0.00003198;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON2S_ee;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==554013){
		mass          =  10.02326;      //Upsilon(2S)
		width         =  0.00003198;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON2S_mumu;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==555){
		mass          =  10.3552;      //Upsilon(3S)
		width         =  0.00002032;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON3S;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==555011){
		mass          =  10.3552;      //Upsilon(3S)
		width         =  0.00002032;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON3S_ee;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else if(_prodParticleId==555013){
		mass          =  10.3552;      //Upsilon(3S)
		width         =  0.00002032;
		Wmin_default  =  mass - 5.*width;
		_maxW          =  mass + 5.*width;
		_particleType=starlightConstants::UPSILON3S_mumu;
		_decayType=starlightConstants::NARROWVMDEFAULT;
	}
	else{
		_maxW=Wmaxinput;
	}
	if (Wmininput==-1){ 
		_minW = Wmin_default;
	}
	else{  
		_minW = Wmininput;
	}

	printInfo << "using the following " << *this;
	return true;
}


//______________________________________________________________________________
starlightConstants::particle inputParameters::getPidTest()
{
	return _particleType;
}


//______________________________________________________________________________
starlightConstants::decaytype inputParameters::getDecayTest()
{
	return _decayType;
}


//______________________________________________________________________________
starlightConstants::interactiontype inputParameters::getInteractionTest()
{
	return _interactionType;
}


//______________________________________________________________________________
double inputParameters::getProtonEnergy()
{
	return _protonEnergy;
}


//______________________________________________________________________________
ostream&
inputParameters::print(ostream& out) const
{
	out << "starlight parameters:" << endl
	    << "    config file name  ...................... '" << _configFileName << "'" << endl
	    << "    beam 1 atomic number ................... " << _beam1Z << endl
	    << "    beam 1 atomic mass number .............. " << _beam1A << endl
	    << "    beam 2 atomic number ................... " << _beam2Z << endl
	    << "    beam 2 atomic mass number .............. " << _beam2A << endl
	    << "    Lorentz gamma of beams in CM frame ..... " << _beamLorentzGamma << endl
	    << "    mass W of produced hadronic system ..... " << _minW << " < W < " << _maxW << " GeV/c^2" << endl
	    << "    # of W bins ............................ " << _nmbWBins << endl
	    << "    maximum absolute value for rapidity .... " << _maxRapidity << endl
	    << "    # of rapidity bins ..................... " << _nmbRapidityBins << endl
	    << "    meson production mode .................. " << _productionMode << endl
	    << "    number of events to generate ........... " << _nmbEventsTot << endl
	    << "    PDG ID of produced particle ............ " << _prodParticleId << endl
	    << "    seed for random generator .............. " << _randomSeed << endl
	    << "    output format .......................... " << _outputFormat << endl
	    << "    breakup mode for beam particles ........ " << _beamBreakupMode << endl
	    << "    interference enabled ................... " << yesNo(_interferenceEnabled) << endl
	    << "    interference strength .................. " << _interferenceStrength << endl
	    << "    coherent scattering off nucleus ........ " << yesNo(_coherentProduction) << endl
	    << "    scaling factor for incoh. VM prod. ..... " << _incoherentFactor << endl
	    << "    bford .................................. " << _bford << endl
	    << "    maximum p_T for interference calc. ..... " << _maxPtInterference << " GeV/c" << endl
	    << "    # of p_T bins for interference calc. ... " << _nmbPtBinsInterference << endl;
	return out;
}


//______________________________________________________________________________
ostream&
inputParameters::write(ostream& out) const
{
	out << beam1Z               () <<endl
	    << beam1A               () <<endl
	    << beam2Z               () <<endl
	    << beam2A               () <<endl
	    << beamLorentzGamma     () <<endl
	    << maxW                 () <<endl
	    << minW                 () <<endl
	    << numWBins             () <<endl
	    << maxRapidity          () <<endl
	    << nmbRapidityBins      () <<endl
	    << productionMode       () <<endl
	    << nmbEvents            () <<endl
	    << prodParticleId       () <<endl
	    << randomSeed           () <<endl
	    << outputFormat         () <<endl
	    << beamBreakupMode      () <<endl
	    << interferenceEnabled  () <<endl
	    << interferenceStrength () <<endl
	    << coherentProduction   () <<endl
	    << incoherentFactor     () <<endl
	    << getbford             () <<endl
	    << maxPtInterference    () <<endl
	    << nmbPtBinsInterference() <<endl;
	return out;
}
