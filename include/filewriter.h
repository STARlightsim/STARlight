
#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <string>
#include "upcevent.h"
#include <fstream>

class FileWriter
{
public:

    /** Default constructor */
    FileWriter();

    /** Constructor with filename */
    FileWriter(std::string filename);

    /** Destructor */
    virtual ~FileWriter();

    /** Open the file */
    int Open();

    /** Open file with given filename */
    int Open(std::string filename);
    
    /** Close the file */
    int Close();

    /** Set the filename we're writing to */
    void SetFileName(std::string filename) {fFilename = filename; }

protected:

   /** The file name */
    std::string fFilename;

    /** The file stream */ 
    std::ofstream fFileStream;

};

#endif // FILEWRITER_H
