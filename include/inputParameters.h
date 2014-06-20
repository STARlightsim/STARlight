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


#ifndef INPUTPARAMETERS_H
#define INPUTPARAMETERS_H


#include "starlightconstants.h"
#include "inputParser.h"
#include "singleton.h"
#include <string>
#include <ostream>
#include <vector>
#include <sstream>

class parameterbase;


class parameterlist
{
public:

    parameterlist() : _parameters(0) {}

    void add(parameterbase* p) {
        _parameters.push_back(p);
    }

    // Returns a string with a key of the current state of the parameter list
    // only
    inline std::string validationKey();
    

private:

    std::vector<parameterbase*> _parameters;

};

// Base class for parameters, needed to keep a list of parameters
class parameterbase
{
public:

    // Add this to parameter list
    parameterbase()
    {
        _parameters.add(this);
    }
    virtual std::string validationkey() = 0;

    template<typename T>
    std::string toString(T v)
    {
        std::stringstream s;
        s << v;
        return s.str();
    }
    inline friend std::ostream& operator<<(std::ostream& os, const parameterbase& par);
 
    // List of all parameters
    static parameterlist _parameters;


   
};
// Need to init the static variable
// parameterlist parameterbase::_parameters;


// The actual parameter class
// validate parameter specifies if the parameter should be a part of the validity check of the current parameters
template<typename T, bool validate>
class parameter : public parameterbase
{
public:

    // Constructor
    parameter(const std::string &name, T value, bool required = true) :parameterbase(),_name(name), _value(value), _validate(validate), _required(required) {}

//     T operator()() const {
//         return _value;
//     }

    parameter &operator=(T v) { _value = v; return *this;}
    T* ptr() const {
        return const_cast<T*>(&_value);
    }
    
    T value() const { return _value; }
    
    std::string name() const { return _name;}
    
    bool required() const { return _required; }
    
    void setValue(T v) { _value = v; }
    
    void setName(std::string name) { _name = name; }
    
    void setRequired(bool r) { _required = r; }
    
    // Validation key for this parameter
    std::string validationkey()
    {
        return (_validate ? _name + ":" + toString(_value) + "-" : std::string(""));
    }

    template<typename S, bool v>
    inline friend std::ostream& operator<<(std::ostream& os, const parameter<S,v>& par);



private:
    std::string _name;

    T _value; // Value
    bool _validate; // true if a change in the parameter invalidates x-sec tables
    bool _required; // true if this is required option.

    parameter();
};

template<typename S, bool v>
std::ostream& operator<<(std::ostream& os, const parameter<S,v>& par)
{
    os << par._value;
    return os;
}

std::ostream& operator<<(std::ostream& os, const parameterbase& par)
{
    os << par._parameters.validationKey(); 
    return os;
}
std::string parameterlist::validationKey()
{
    std::stringstream s;
    for(unsigned int i = 0; i < _parameters.size(); ++i)
    {
        s << _parameters[i]->validationkey(); // Will print names and values of validation parameters
    }
    return s.str();
}

class inputParameters {

private:
	// inputParameters is now a singleton
	friend class Singleton<inputParameters>;
	inputParameters();
public:

	~inputParameters();

	bool init();
	bool configureFromFile(const std::string &configFileName = "./config/slight.in");

	unsigned int beam1Z                () const { return _beam1Z.value();                 }  ///< returns atomic number of beam particle 1
	unsigned int beam1A                () const { return _beam1A.value();                 }  ///< returns atomic mass number of beam particle 1
	unsigned int beam2Z                () const { return _beam2Z.value();                 }  ///< returns atomic number of beam particle 2
	unsigned int beam2A                () const { return _beam2A.value();                 }  ///< returns atomic mass number of beam particle 2
	double       beamLorentzGamma      () const { return _beamLorentzGamma;       }  ///< returns Lorentz gamma factor of both beams in beam CMS frame
	double       beam1LorentzGamma     () const { return _beam1LorentzGamma.value();      }  ///< returns Lorentz gamma factor of beam 1 in collider frame
	double       beam2LorentzGamma     () const { return _beam2LorentzGamma.value();      }  ///< returns Lorentz gamma factor of beam 2 in collider frame
	double       maxW                  () const { return _maxW.value();                   }  ///< returns maximum mass W of produced hadronic system [GeV/c^2]
	double       minW                  () const { return _minW.value();                   }  ///< returns minimum mass W of produced hadronic system [GeV/c^2]
	unsigned int nmbWBins              () const { return _nmbWBins.value();               }  ///< returns number of W bins in lookup table
	double       maxRapidity           () const { return _maxRapidity.value();            }  ///< returns maximum absolute value of rapidity
	unsigned int nmbRapidityBins       () const { return _nmbRapidityBins.value();        }  ///< returns number of rapidity bins in lookup table
	bool         ptCutEnabled          () const { return _ptCutEnabled.value();           }  ///< returns cut in pt
	double       ptCutMin              () const { return _ptCutMin.value();               }  ///< returns minimum pt
	double       ptCutMax              () const { return _ptCutMax.value();               }  ///< returns maximum pt
	bool         etaCutEnabled         () const { return _etaCutEnabled.value();          }  ///< returns cut in eta
	double       etaCutMin             () const { return _etaCutMin.value();              }  ///< returns minimum eta
	double       etaCutMax             () const { return _etaCutMax.value();              }  ///< returns maximum eta
	int          productionMode        () const { return _productionMode.value();         }  ///< returns production mode
	unsigned int nmbEvents             () const { return _nmbEventsTot.value();           }  ///< returns total number of events to generate
	int          prodParticleId        () const { return _prodParticleId.value();         }  ///< returns PDG particle ID of produced particle
	int          randomSeed            () const { return _randomSeed.value();             }  ///< returns seed for random number generator
	int          beamBreakupMode       () const { return _beamBreakupMode.value();        }  ///< returns breakup mode for beam particles
	bool         interferenceEnabled   () const { return _interferenceEnabled.value();    }  ///< returns whether interference is taken into account
	double       interferenceStrength  () const { return _interferenceStrength.value();   }  ///< returns percentage of interference
	double       maxPtInterference     () const { return _maxPtInterference.value();      }  ///< returns maximum p_T for interference calculation [GeV/c]
	int          nmbPtBinsInterference () const { return _nmbPtBinsInterference.value();  }  ///< returns number of p_T bins for interference calculation
	double       ptBinWidthInterference() const { return _ptBinWidthInterference.value(); }  ///< returns width of p_T bins for interference calculation [GeV/c]
	bool         coherentProduction    () const { return _coherentProduction.value();     }  ///< returns whether production is coherent or incoherent
	double       incoherentFactor      () const { return _incoherentFactor.value();       }  ///< returns incoherent contribution in vector meson production
	double 	     minGammaEnergy        () const { return _minGammaEnergy.value();         }  ///< returns minimum gamma energy in case of photo nuclear processes [GeV]
	double       maxGammaEnergy        () const { return _maxGammaEnergy.value();         }  ///< returns maximum gamma energy in case of photo nuclear processes [GeV]
	std::string  pythiaParams          () const { return _pythiaParams.value();           }  ///< returns parameters to be passed to pythia
	bool         pythiaFullEventRecord () const { return _pythiaFullEventRecord.value();  }  ///< returns if the full pythia event record should be printed
	int	     xsecCalcMethod        () const { return _xsecCalcMethod.value();         }  ///< returns the method used for the x-sec calculation
	int          nThreads              () const { return _nThreads.value();               }  ///< returns the number of threads in case method 1 is used for the x-sec calc
	unsigned int nBinsQKniehl          () const { return _nBinsQKniehl.value();           }  ///< Number of bins in Q used for the transformation to the impact paramter space of the Kniehl function
        unsigned int nBinsEKniehl          () const { return _nBinsEKniehl.value();           }  ///< Number of bins in photon energy used for the Kniehl function
        unsigned int nBinsBKniehl          () const { return _nBinsBKniehl.value();           }  ///< Number of bins in impact parameter used for the Kniehl function
        double       qMaxKniehl            () const { return _qMaxKniehl.value();             }  ///< Max value of Q used for the Kniehl funcion
        double       eGammaMinKniehl       () const { return _eGammaMinKniehl.value();        }  ///< Min value of gamma energy used for the Kniehl funcion
        double       eGammaMaxKniehl       () const { return _eGammaMaxKniehl.value();        }  ///< Max value of gamma energy used for the Kniehl funcion
        double       bMinKniehl            () const { return _bMinKniehl.value();             }  ///< Min value of impact parameter used for the Kniehl funcion
        double       bMaxKniehl            () const { return _bMaxKniehl.value();             }  ///< Max value of impact parameter used for the Kniehl funcion
        
	starlightConstants::particleTypeEnum    prodParticleType     () const { return _particleType;    }  ///< returns type of produced particle
	starlightConstants::decayTypeEnum       prodParticleDecayType() const { return _decayType;       }  ///< returns decay type of produced particle
	starlightConstants::interactionTypeEnum interactionType      () const { return _interactionType; }  ///< returns interaction type
	// double vmPhotonCoupling();
	// double slopeParameter();
	double protonEnergy                () const { return _protonEnergy.value(); }

	void setBeam1Z                (unsigned int v)  {  _beam1Z = v;                 }  ///< returns atomic number of beam particle 1
	void setBeam1A                (unsigned int v)  {  _beam1A = v;                 }  ///< returns atomic mass number of beam particle 1
	void setBeam2Z                (unsigned int v)  {  _beam2Z = v;                 }  ///< returns atomic number of beam particle 2
	void setBeam2A                (unsigned int v)  {  _beam2A = v;                 }  ///< returns atomic mass number of beam particle 2
	void setBeamLorentzGamma      (double v)  {  _beamLorentzGamma = v;       }  ///< returns Lorentz gamma factor of both beams in beam CMS frame
	void setBeam1LorentzGamma     (double v)  {  _beam1LorentzGamma = v;      }  ///< returns Lorentz gamma factor of beam 1 in collider frame
	void setBeam2LorentzGamma     (double v)  {  _beam2LorentzGamma = v;      }  ///< returns Lorentz gamma factor of beam 2 in collider frame
	void setMaxW                  (double v)  {  _maxW = v;                   }  ///< returns maximum mass W of produced hadronic system [GeV/c^2]
	void setMinW                  (double v)  {  _minW = v;                   }  ///< returns minimum mass W of produced hadronic system [GeV/c^2]
	void setNmbWBins              (unsigned int v)  {  _nmbWBins = v;               }  ///< returns number of W bins in lookup table
	void setMaxRapidity           (double v)  {  _maxRapidity = v;            }  ///< returns maximum absolute value of rapidity
	void setNmbRapidityBins       (unsigned int v)  {  _nmbRapidityBins = v;        }  ///< returns number of rapidity bins in lookup table
	void setPtCutEnabled          (bool v)  {  _ptCutEnabled = v;           }  ///< returns cut in pt
	void setPtCutMin              (double v)  {  _ptCutMin = v;               }  ///< returns minimum pt
	void setPtCutMax              (double v)  {  _ptCutMax = v;               }  ///< returns maximum pt
	void setEtaCutEnabled         (bool v)  {  _etaCutEnabled = v;          }  ///< returns cut in eta
	void setEtaCutMin             (double v)  {  _etaCutMin = v;              }  ///< returns minimum eta
	void setEtaCutMax             (double v)  {  _etaCutMax = v;              }  ///< returns maximum eta
	void setProductionMode        (int v)  {  _productionMode = v;         }  ///< returns production mode
	void setNmbEvents             (unsigned int v)  {  _nmbEventsTot = v;           }  ///< returns total number of events to generate
	void setProdParticleId        (int v)  {  _prodParticleId = v;         }  ///< returns PDG particle ID of produced particle
	void setRandomSeed            (int v)  {  _randomSeed = v;             }  ///< returns seed for random number generator
	void setBeamBreakupMode       (int v)  {  _beamBreakupMode = v;        }  ///< returns breakup mode for beam particles
	void setInterferenceEnabled   (bool v)  {  _interferenceEnabled = v;    }  ///< returns whether interference is taken into account
	void setInterferenceStrength  (double v)  {  _interferenceStrength = v;   }  ///< returns percentage of interference
	void setMaxPtInterference     (double v)  {  _maxPtInterference = v;      }  ///< returns maximum p_T for voiderference calculation [GeV/c]
	void setNmbPtBinsInterference (int v)  {  _nmbPtBinsInterference = v;  }  ///< returns number of p_T bins for interference calculation
	void setPtBinWidthInterference(double v)  {  _ptBinWidthInterference = v; }  ///< returns width of p_T bins for voiderference calculation [GeV/c]
	void setCoherentProduction    (bool v)  {  _coherentProduction = v;     }  ///< returns whether production is coherent or incoherent
	void setIncoherentFactor      (double v)  {  _incoherentFactor = v;       }  ///< returns incoherent contribution in vector meson production
	void setMinGammaEnergy        (double v)  {  _minGammaEnergy = v;         }  ///< returns minimum gamma energy in case of photo nuclear processes [GeV]
	void setMaxGammaEnergy        (double v)  {  _maxGammaEnergy = v;         }  ///< returns maximum gamma energy in case of photo nuclear processes [GeV]
	void setPythiaParams          (std::string v)  {  _pythiaParams = v;           }  ///< returns parameters to be passed to pythia
	void setPythiaFullEventRecord (bool v)  {  _pythiaFullEventRecord = v;  }  ///< returns if the full pythia event record should be prvoided
	void setXsecCalcMethod        (int v)  {  _xsecCalcMethod = v;         }  ///< returns the method used for the x-sec calculation
	void setNThreads              (int v)  {  _nThreads = v;               }  ///< returns the number of threads in case method 1 is used for the x-sec calc
	void setNBinsQKniehl          (unsigned int v)  {  _nBinsQKniehl = v;           }  ///< Number of bins in Q used for the transformation to the impact paramter space of the Kniehl function
        void setNBinsEKniehl          (unsigned int v)  {  _nBinsEKniehl = v;           }  ///< Number of bins in photon energy used for the Kniehl function
        void setNBinsBKniehl          (unsigned int v)  {  _nBinsBKniehl = v;           }  ///< Number of bins in impact parameter used for the Kniehl function
        void setQMaxKniehl            (double v)  {  _qMaxKniehl = v;             }  ///< Max value of Q used for the Kniehl funcion
        void setEGammaMinKniehl       (double v)  {  _eGammaMinKniehl = v;        }  ///< Min value of gamma energy used for the Kniehl funcion
        void setEGammaMaxKniehl       (double v)  {  _eGammaMaxKniehl = v;        }  ///< Max value of gamma energy used for the Kniehl funcion
        void setBMinKniehl            (double v)  {  _bMinKniehl = v;             }  ///< Min value of impact parameter used for the Kniehl funcion
        void setBMaxKniehl            (double v)  {  _bMaxKniehl = v;             }  ///< Max value of impact parameter used for the Kniehl funcion
        
	void setProdParticleType      (starlightConstants::particleTypeEnum v)   { _particleType = v;    }  ///< returns type of produced particle
	void setProdParticleDecayType (starlightConstants::decayTypeEnum v)        { _decayType = v;       }  ///< returns decay type of produced particle
	void setInteractionType       (starlightConstants::interactionTypeEnum v)  { _interactionType = v; }  ///< returns interaction type
	 
	// double vmPhotonCoupling();
	// double slopeParameter();
	void setProtonEnergy        (double v)  { _protonEnergy = v; }
	
	template<typename T>
	inline bool setParameter(std::string expression);
	
	std::ostream& print(std::ostream& out) const;  ///< prints parameter summary
	std::ostream& write(std::ostream& out) const;  ///< writes parameters back to an ostream
	
	std::string parameterValueKey() const; ///< Generates key for the current parameters

  
private:

    
// To indicate if the crossection table should be re-calculated if parameter changes
#define VALIDITY_CHECK true
#define NO_VALIDITY_CHECK false
	
	std::string _configFileName;  ///< path to configuration file (default = ./config/slight.in)

	// config file parameters
	parameter<unsigned int,VALIDITY_CHECK>     _beam1Z;                  ///< atomic number of beam particle 1
	parameter<unsigned int,VALIDITY_CHECK>     _beam1A;                  ///< atomic mass number of beam particle 1
	parameter<unsigned int,VALIDITY_CHECK>     _beam2Z;                  ///< atomic number of beam particle 2
	parameter<unsigned int,VALIDITY_CHECK>     _beam2A;                  ///< atomic mass number of beam particle 2
	parameter<double, VALIDITY_CHECK>          _beam1LorentzGamma;       ///< Lorentz gamma factor of beam 1 in collider frame
	parameter<double, VALIDITY_CHECK>          _beam2LorentzGamma;       ///< Lorentz gamma factor of beam 2 in collider frame
	parameter<double, VALIDITY_CHECK>          _maxW;                    ///< maximum mass W of produced hadronic system [GeV/c^2]
	parameter<double, VALIDITY_CHECK>          _minW;                    ///< minimum mass W of produced hadronic system; if set to -1 default value is taken [GeV/c^2]
	parameter<unsigned int, VALIDITY_CHECK>    _nmbWBins;                ///< number of W bins in lookup table
	parameter<double, VALIDITY_CHECK>          _maxRapidity;             ///< maximum absolute value of rapidity
	parameter<unsigned int, VALIDITY_CHECK>    _nmbRapidityBins;         ///< number of rapidity bins in lookup table
	parameter<bool, VALIDITY_CHECK>            _ptCutEnabled;            ///< en/disables cut in pt
	parameter<double, VALIDITY_CHECK>          _ptCutMin;                ///< minimum pt, if cut is enabled
	parameter<double, VALIDITY_CHECK>          _ptCutMax;                ///< maximum pt, if cut is enabled
	parameter<bool, VALIDITY_CHECK>            _etaCutEnabled;           ///< en/disables cut in eta
	parameter<double, VALIDITY_CHECK>          _etaCutMin;               ///< minimum eta, if cut is enabled
	parameter<double, VALIDITY_CHECK>          _etaCutMax;               ///< maximum eta, if cut is enabled
	parameter<unsigned int, VALIDITY_CHECK>    _productionMode;          ///< \brief production mode
	                                                                     ///<
	                                                                     ///< 1 = photon-photon fusion,
	                                                                     ///< 2 = narrow vector meson resonance in photon-Pomeron fusion,
	                                                                     ///< 3 = Breit-Wigner vector meson resonance in photon-Pomeron fusion
	parameter<unsigned int, VALIDITY_CHECK>    _nmbEventsTot;            ///< total number of events to generate
	parameter<unsigned int, VALIDITY_CHECK>    _prodParticleId;          ///< PDG particle ID of produced particle
	parameter<unsigned int, VALIDITY_CHECK>    _randomSeed;              ///< seed for random number generator
	                                                                     ///<
	                                                                     ///< 1 = ASCII
	                                                                     ///< 2 = GSTARtext,
	                                                                     ///< 3 = PAW ntuple (not working)
	parameter<unsigned int, VALIDITY_CHECK>    _beamBreakupMode;         ///< \brief breakup mode for beam particles
	                                                                     ///<
	                                                                     ///< 1 = hard sphere nuclei (b > 2R),
	                                                                     ///< 2 = both nuclei break up (XnXn),
	                                                                     ///< 3 = a single neutron from each nucleus (1n1n),
	                                                                     ///< 4 = neither nucleon breaks up (with b > 2R),
	                                                                     ///< 5 = no hadronic break up (similar to option 1, but with the actual hadronic interaction)
	parameter<bool, VALIDITY_CHECK>            _interferenceEnabled;     ///< if VALIDITY_CHECK, interference is taken into account
	parameter<double, VALIDITY_CHECK>          _interferenceStrength;    ///< percentage of interference: from 0 = none to 1 = full
	parameter<double, VALIDITY_CHECK>          _maxPtInterference;       ///< maximum p_T for interference calculation [GeV/c]
	parameter<unsigned int, VALIDITY_CHECK>    _nmbPtBinsInterference;   ///< number of p_T bins for interference calculation
	parameter<double, VALIDITY_CHECK>          _ptBinWidthInterference;  ///< width of p_T bins for interference calculation [GeV/c]
	parameter<bool, VALIDITY_CHECK>            _coherentProduction;      ///< if VALIDITY_CHECK, production is coherent, else incoherent
	parameter<double, VALIDITY_CHECK>          _incoherentFactor;        ///< allows to scale the incoherent contribution in vector meson production
	parameter<double, VALIDITY_CHECK>          _protonEnergy;
	parameter<double, VALIDITY_CHECK>          _minGammaEnergy;          ///< minimum gamma energy in case of photo nuclear processes [GeV]
	parameter<double, VALIDITY_CHECK>          _maxGammaEnergy;          ///< maximum gamma energy in case of photo nuclear processes [GeV]
	parameter<std::string,NO_VALIDITY_CHECK>   _pythiaParams;            ///< semi-colon separated parameters to pass to pythia, e.g. "mstj(1)=0;paru(13)=0.1" 
	parameter<bool, NO_VALIDITY_CHECK>         _pythiaFullEventRecord;   ///< if the full pythia event record should be in the output
	parameter<unsigned int, VALIDITY_CHECK>    _xsecCalcMethod;	     ///< Select x-sec calc method. (0 is standard starlight method, 1 must be used for assym. collisions (e.g. p-A), but is slow)
	parameter<unsigned int, NO_VALIDITY_CHECK> _nThreads;	             ///< Number of threads used in the case of using method 1 for calculating the x-sections
        parameter<unsigned int, VALIDITY_CHECK>    _nBinsQKniehl;            ///< Number of bins in Q used for the transformation to the impact paramter space of the Kniehl function
        parameter<unsigned int, VALIDITY_CHECK>    _nBinsEKniehl;            ///< Number of bins in photon energy used for the Kniehl function
        parameter<unsigned int, VALIDITY_CHECK>    _nBinsBKniehl;            ///< Number of bins in impact parameter used for the Kniehl function
        parameter<double, VALIDITY_CHECK>          _qMaxKniehl;              ///< Max value of Q used for the Kniehl funcion
        parameter<double, VALIDITY_CHECK>          _eGammaMinKniehl;         ///< Min value of gamma energy used for the Kniehl funcion
        parameter<double, VALIDITY_CHECK>          _eGammaMaxKniehl;         ///< Max value of gamma energy used for the Kniehl funcion
        parameter<double, VALIDITY_CHECK>          _bMinKniehl;              ///< Min value of impact parameter used for the Kniehl funcion
        parameter<double, VALIDITY_CHECK>          _bMaxKniehl;              ///< Max value of impact parameter used for the Kniehl funcion
        
	
	starlightConstants::particleTypeEnum       _particleType;
	starlightConstants::decayTypeEnum          _decayType;
	starlightConstants::interactionTypeEnum    _interactionType;

	double                         _beamLorentzGamma;                    ///< Lorentz gamma factor of the beams in CMS frame, not an input parameter
	
	inputParser _ip;
	
};

#define inputParametersInstance Singleton<inputParameters>::instance()

template<typename T>
inline 
bool inputParameters::setParameter(std::string expression)
{
   
    return _ip.parseString(expression);
   
   
}

inline
std::ostream&
operator <<(std::ostream&          out,
            const inputParameters& par)
{
	return par.print(out);
}
 

#endif  // INPUTPARAMETERS_H
