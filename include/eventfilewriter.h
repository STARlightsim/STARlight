
#ifndef EVENTFILEWRITER_H
#define EVENTFILEWRITER_H

#include <string>
#include "filewriter.h"

class EventFileWriter : public FileWriter
{
   public:
      
      /** Default constructor */
      EventFileWriter();
      
      /** Constructor with name */
      EventFileWriter(std::string filename);

      /** Write an UPC event to file */
      int WriteEvent(UPCEvent &event, int eventnumber);
      
      
};

#endif // EVENTFILEWRITER_H
