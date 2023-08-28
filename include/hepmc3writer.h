
#ifndef HEPMC3WRITER_H
#define HEPMC3WRITER_H

#include "inputParameters.h"
#include "upcXevent.h"
#include  "HepMC3/WriterAscii.h"
#include "HepMC3/FourVector.h"

using HepMC3::WriterAscii;
using HepMC3::FourVector;
using namespace starlightConstants;

class hepMC3Writer
{

 public:
  /** Default Constructor **/
  hepMC3Writer();
  ~hepMC3Writer();
  
  /** Init HepMC3 Writer with input parameters **/
  int initWriter(const inputParameters & param);
  
  /** Write an eX event to file  **/
  int writeEvent(const upcXEvent &event, int eventnumber);

  /* closes file */
  void close(){ _hepmc3_output->close(); };
  
 private:

  /** Init HepMC3 Beam Four Vectors **/
  int initBeamHepMC3(const inputParameters &param);
  
  WriterAscii * _hepmc3_output;
  FourVector hepmc3_beam1_four_vector;
  FourVector hepmc3_beam2_four_vector;
  int beam1_pdg_id;
  int beam2_pdg_id;
  particleTypeEnum VM_pdg_id;//Particle type or pdg code for the produced vector meson
  particleTypeEnum PID;//channel id or prod_pid for the channel of interest
  decayTypeEnum _decay;
  
};

#endif // HEPMC3WRITER_H
