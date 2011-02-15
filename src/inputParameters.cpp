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
#include "inputParameters.h"
#include "inputParser.h"


using namespace std;
using namespace starlightConstants;


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
	  _accCutPt              (0),
	  _minPt                 (0),
	  _maxPt                 (0),
	  _accCutEta             (0),
	  _minEta                (0),
	  _maxEta                (0),
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
	  _deuteronSlopePar      (0),
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

	double minWConfigFile = 0;
	double maxWConfigFile = 0;
	
	inputParser ip;
	
	ip.addUintParameter(string("BEAM_1_Z"), &_beam1Z);
	ip.addUintParameter(string("BEAM_2_Z"), &_beam2Z);
	ip.addUintParameter(string("BEAM_1_A"), &_beam1A);
	ip.addUintParameter(string("BEAM_2_A"), &_beam2A);

	ip.addDoubleParameter(string("BEAM_GAMMA"), &_beamLorentzGamma);
	
	ip.addDoubleParameter(string("W_MAX"), &maxWConfigFile);
	ip.addDoubleParameter(string("W_MIN"), &minWConfigFile);
	ip.addUintParameter(string("W_N_BINS"), &_nmbWBins);;
	
	ip.addDoubleParameter(string("RAP_MAX"), &_maxRapidity);
	ip.addUintParameter(string("RAP_N_BINS"), &_nmbRapidityBins);
	
	ip.addBoolParameter(string("CUT_PT"), &_accCutPt);
	ip.addDoubleParameter(string("PT_MIN"), &_minPt);
	ip.addDoubleParameter(string("PT_MAX"), &_maxPt);
	
	ip.addBoolParameter(string("CUT_ETA"), &_accCutEta);
	ip.addDoubleParameter(string("ETA_MIN"), &_minEta);
	ip.addDoubleParameter(string("ETA_MAX"), &_maxEta);
	
	ip.addIntParameter(string("PROD_MODE"), &_productionMode);
	
	ip.addUintParameter(string("N_EVENTS"), &_nmbEventsTot);
	
	ip.addIntParameter(string("PROD_PID"), &_prodParticleId);
	
	ip.addIntParameter(string("RND_SEED"), &_randomSeed);
	
	ip.addIntParameter(string("OUTPUT_FORMAT"), &_outputFormat);
	
	ip.addIntParameter(string("BREAKUP_MODE"), &_beamBreakupMode);
	
	ip.addBoolParameter(string("INTERFERENCE"), &_interferenceEnabled);
	ip.addDoubleParameter(string("IF_STRENGTH"), &_interferenceStrength);
	
	ip.addBoolParameter(string("COHERENT"), &_coherentProduction);
	ip.addDoubleParameter(string("INCO_FACTOR"), &_incoherentFactor);
	
	ip.addDoubleParameter(string("BFORD"), &_deuteronSlopePar);

	ip.addDoubleParameter(string("INT_PT_MAX"), &_maxPtInterference);
	ip.addIntParameter(string("INT_PT_N_BINS"), &_nmbPtBinsInterference);
	
	int nParameters = ip.parseFile(_configFileName);
	if(nParameters == -1) 
	{
		printWarn << "could not open file '" << _configFileName << "'" << endl;
	  return false;
	}
	
	//ip.printParameterInfo(printInfo);
	
	if(ip.validateParameters(printErr))
	{
	  printInfo << "successfully read input parameters from '" << _configFileName << "'" << endl;
	}
 	else {
 		printWarn << "problems reading input parameters from '" << _configFileName << "'" << endl
 		          << *this;
 		return false;
 	}

	_ptBinWidthInterference = _maxPtInterference / _nmbPtBinsInterference;
	_protonEnergy           = _beamLorentzGamma * protonMass;

	// define interaction type
	switch (_productionMode) {
	case 1:
		_interactionType = PHOTONPHOTON;
		break;
	case 2:
		_interactionType = PHOTONPOMERONNARROW;
		break;
	case 3:
		_interactionType = PHOTONPOMERONWIDE;
		break;
	default:
		printWarn << "unknown production mode '" << _productionMode << "'" << endl;
		return false;
	}
  
	//Trying to define the proper Wmins and Wmaxs. a TEMPORARY fix....Better solution=??
	double mass        = 0;
	double width       = 0;
	double defaultMinW = 0;  // default for _minW, unless it is defined later [GeV/c^2]
	switch (_prodParticleId) {
	case 11:  // e+e- pair
		_particleType = ELECTRON;
		_decayType    = LEPTONPAIR;
		defaultMinW   = 0.01; // default is 0.01; up to 0.15 is safe for Summer 2000 triggering for e+e- pairs
		_maxW         = maxWConfigFile;
		break;
	case 13:  // mu+mu- pair
		_particleType = MUON;
		_decayType    = LEPTONPAIR;
		defaultMinW   = 2 * muonMass;
		_maxW         = maxWConfigFile;
		break;
	case 15:  // tau+tau- pair
		_particleType = TAUON;
		_decayType    = LEPTONPAIR;
		defaultMinW   = 2 * tauMass;
		_maxW         = maxWConfigFile;
		break;
	case 115:  // a_2(1320)
		_particleType = A2;
		_decayType    = SINGLEMESON;
		_maxW         = maxWConfigFile;
		break;
	case 221:  // eta
		_particleType = ETA;
		_decayType    = SINGLEMESON;
		_maxW         = maxWConfigFile;
		break;
	case 225:  // f_2(1270)
		_particleType = F2;
		_decayType    = SINGLEMESON;
		_maxW         = maxWConfigFile;
		break;
	case 331:  // eta'(958)
		_particleType = ETAPRIME;
		_decayType    = SINGLEMESON;
		_maxW         = maxWConfigFile;
		break;
	case 335:  // f_2'(1525)
		_particleType = F2PRIME;
		_decayType    = SINGLEMESON;
		_maxW         = maxWConfigFile;
		break;
	case 441:  // eta_c(1s)
		_particleType = ETAC;
		_decayType    = SINGLEMESON;
		_maxW         = maxWConfigFile;
		break;
	case 9010221:  // f_0(980), was orginally called 10221? updated to standard number
		_particleType = F0;
		_decayType    = SINGLEMESON;
		_maxW         = maxWConfigFile;
		break;
	case 33:  // Z"/Z03
		_particleType = ZOVERZ03;
		_decayType    = SINGLEMESON;
		_maxW         = maxWConfigFile;
		break;
	case 113:  // rho(770)
		_particleType = RHO;
		_decayType    = WIDEVMDEFAULT;
		mass          = 0.7685;
		width         = 0.1507;
		defaultMinW   = 2 * pionChargedMass;
		_maxW         = mass + 5 * width;
		break;
	case 913:  // rho(770) with direct pi+pi- decay, interference given by ZEUS data
		_particleType = RHOZEUS;
		_decayType    = WIDEVMDEFAULT;
		mass          = 0.7685;                    
		width         = 0.1507;
		defaultMinW   = 2 * pionChargedMass;
		_maxW         = mass + 5 * width;  // use the same 1.5GeV max mass as ZEUS
		break;
	case 999:  // pi+pi-pi+pi- phase space decay
		_particleType = FOURPRONG;
		_decayType    = WIDEVMDEFAULT;
		mass          = 1.350;
		width         = 0.360;
		defaultMinW   = 4 * pionChargedMass;
		_maxW         = 3;
		break;
	case 223:  // omega(782)
		_particleType = OMEGA;
		_decayType    = NARROWVMDEFAULT;  // will probably be moved to 3-body decay
		mass          = 0.78194;
		width         = 0.00843;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 333:  // phi(1020)
		_particleType = PHI;
		_decayType    = NARROWVMDEFAULT;
		mass          = 1.019413;
		width         = 0.00443;
		defaultMinW   = 2 * kaonChargedMass;
		_maxW         = mass + 5 * width;
		break;
	case 443:  // J/psi
		_particleType = JPSI;
		_decayType    = NARROWVMDEFAULT;
		mass          = 3.09692;   // JN  3.09688;
		width         = 0.000091;  // JN  0.000087;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 443011:  // J/psi
		_particleType = JPSI_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = 3.09692;   // JN  3.09688;
		width         = 0.000091;  // JN  0.000087;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 443013:  // J/psi
		_particleType = JPSI_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = 3.09692;   // JN  3.09688;
		width         = 0.000091;  // JN  0.000087;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 444:  // J/psi
		_particleType = JPSI2S;
		_decayType    = NARROWVMDEFAULT;
		mass          = 3.686093;
		width         = 0.000337;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 444011: // J/psi
		_particleType = JPSI2S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = 3.686093;     
		width         = 0.000337;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 444013:  // J/psi
		_particleType = JPSI2S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = 3.686093;
		width         = 0.000337;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 553:  // Upsilon
		_particleType = UPSILON;
		_decayType    = NARROWVMDEFAULT;
		mass          = 9.46030;
		width         = 0.00005402;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 553011:  // Upsilon
		_particleType = UPSILON_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = 9.46030;
		width         = 0.00005402;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 553013:  // Upsilon
		_particleType = UPSILON_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = 9.46030;
		width         = 0.00005402;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 554:  // Upsilon(2S)
		_particleType = UPSILON2S;
		_decayType    = NARROWVMDEFAULT;
		mass          = 10.02326;
		width         = 0.00003198;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 554011:  // Upsilon(2S)
		_particleType = UPSILON2S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = 10.02326;
		width         = 0.00003198;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 554013:  // Upsilon(2S)
		_particleType = UPSILON2S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = 10.02326;
		width         = 0.00003198;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 555:  // Upsilon(3S)
		mass          = 10.3552;
		width         = 0.00002032;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		_particleType = UPSILON3S;
		_decayType    = NARROWVMDEFAULT;
		break;
	case 555011:  // Upsilon(3S)
		_particleType = UPSILON3S_ee;
		_decayType    = NARROWVMDEFAULT;
		mass          = 10.3552;
		width         = 0.00002032;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	case 555013:  // Upsilon(3S)
		_particleType = UPSILON3S_mumu;
		_decayType    = NARROWVMDEFAULT;
		mass          = 10.3552;
		width         = 0.00002032;
		defaultMinW   = mass - 5 * width;
		_maxW         = mass + 5 * width;
		break;
	default:
		printWarn << "unknown particle ID " << _prodParticleId << endl;
		return false;
	}  // _prodParticleId

	if (minWConfigFile == -1)
		_minW = defaultMinW;
	else
		_minW = minWConfigFile;

	printInfo << "using the following " << *this;
	return true;
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
	    << "    cut in pT............................... " << yesNo(_accCutPt) << endl
	    << "    minumum pT.............................. " << _minPt << endl
            << "    maximum pT.............................. " << _maxPt << endl
	    << "    cut in eta.............................. " << yesNo(_accCutEta) << endl
	    << "    minumum eta............................. " << _minEta << endl
            << "    maximum eta............................. " << _maxEta << endl
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
	    << "    deuteron slope parameter ............... " << _deuteronSlopePar << " (GeV/c)^{-2}" << endl
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
	    << nmbWBins             () <<endl
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
	    << deuteronSlopePar     () <<endl
	    << maxPtInterference    () <<endl
	    << nmbPtBinsInterference() <<endl;
	return out;
}
