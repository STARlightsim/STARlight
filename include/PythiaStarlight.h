
#ifndef PYTHIASTARLIGHT_H
#define PYTHIASTARLIGHT_H

#include "Pythia.h"
#include <string>

class PythiaStarlight
{
   public:

      PythiaStarlight();
      int Init(std::string xmldocpath);

      Pythia8::Pythia* GetPythia() const { return fPythia; }
      
   private:
      
      Pythia8::Pythia* fPythia;
      
};


      // Generator; shorthand for event.
        //Pythia pythia("/home/butter/pythia/pythia8120/xmldoc");
        //Event& event = pythia.event;

#endif//PYTHIASTARLIGHT_H
