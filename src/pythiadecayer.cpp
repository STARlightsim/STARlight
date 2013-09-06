#include "pythiadecayer.h"
#include "reportingUtils.h"
#include "starlightconfig.h"
using namespace Pythia8;



pythiaDecayer::pythiaDecayer() :
    _pythia(PYTHIA8_SETTINGS_DIR)
{}
/*
pythiaDecayer::pythiaDecayer(const pythiaDecayer &obj) :
    _pythia(obj._pythia)
{}
*/
pythiaDecayer::~pythiaDecayer()
{}
/*
pythiaDecayer& pythiaDecayer::operator=(const pythiaDecayer &other)
{
    if (this != &other)
    {
        _pythia = other._pythia;
    }
    return *this;
}
*/
void pythiaDecayer::init()
{
  _pythia.readString("ProcessLevel:all = off"); 
  _pythia.readString("Standalone:allowResDec = on");
  _pythia.readString("Next:numberShowEvent = 0");
  _pythia.init();                               
  _pythia.event.reset();
}

void pythiaDecayer::addParticle(const starlightParticle &p)
{
  
  Event &pyEvent = _pythia.event;
  int status = 23; // Outgoing particle from the hardest sub-process
  int col = 0;
  int acol = 0;
  int code = p.getPdgCode();

  pyEvent.append(code, status, col, acol, p.GetPx(), p.GetPy(), p.GetPz(), p.GetE(), p.M());
  
}

upcEvent pythiaDecayer::execute()
{
  upcEvent slEvent;
  
  Event &pyEvent = _pythia.event;
  _pythia.forceTimeShower(1, 2, 100000.0);
//  pyEvent.list();
  if(!_pythia.next())
  {
    printWarn << "Pythia::next() failed" << std::endl;
    return upcEvent();
  }
  
  for(int i = 0; i < pyEvent.size(); ++i)
  {
    
      Particle p = pyEvent[i];
      starlightParticle slPart(p.px(), p.py(), p.pz(), p.e(), p.mass(), p.idAbs()*(p.charge()<0?-1:1), p.charge(),
			       p.xProd(), p.yProd(), p.zProd(), p.tProd(),
			       p.mother1(), p.mother2(), p.daughter1(), p.daughter2(), p.status());
      slEvent.addParticle(slPart);
  }
  pyEvent.clear();
  pyEvent.reset(); 
  return slEvent;
  
}