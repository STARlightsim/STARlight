// beam.cpp
/*
 * $Id: beam.cpp,v 1.0 2010/07/04   $
 *
 * /author you name or location where you obtained this code 
 *
 * $Log: $
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *

 */

#include "filewriter.h"
#include <iostream>
#include <exception>
#include <cstdlib>

FileWriter::FileWriter() :
        fFilename("")
        ,fFileStream()
{
}

FileWriter::FileWriter(std::string filename) :
        fFilename(filename)
        ,fFileStream(filename.c_str())
{

}

FileWriter::~FileWriter()
{
}

int FileWriter::Open()
{
    try
    {
        fFileStream.open(fFilename.c_str());
    }
    catch (const std::ios::failure & error)
    {
        std::cerr << "I/O exception: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}

int FileWriter::Open(std::string filename)
{
    fFilename = filename;
    return Open();
}

int FileWriter::Close()
{
    try
    {
        fFileStream.close();
    }
    catch (const std::ios::failure & error)
    {
        std::cerr << "I/O exception: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }
    return 0;

}
