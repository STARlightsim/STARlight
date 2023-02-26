
#ifndef HEPMC3WRITER_H
#define HEPMC3WRITER_H

#include "inputParameters.h"
#include "upcXevent.h"
#include  "HepMC3/WriterAscii.h"
#include "HepMC3/FourVector.h"

using HepMC3::WriterAscii;
using HepMC3::FourVector;

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
  FourVector hepmc3_electronBeam_four_vector_;
  FourVector hepmc3_targetBeam_four_vector_;
  int electronBeam_pdg_id_;
  int targetBeam_pdg_id_;
  
};

#endif // HEPMC3WRITER_H
