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
//      some simple streams for reporting plus some stream operators
//      for common STL classes
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef REPORTINGUTILS_H
#define REPORTINGUTILS_H


#include <iostream>
#include <string>
#include <vector>


//////////////////////////////////////////////////////////////////////////////
// macros for printing errors, warnings, and infos

// cuts out block "className::methodName" from __PRETTY_FUNCTION__ output
inline
std::string
getClassMethod__(std::string prettyFunction)
{
	size_t pos = prettyFunction.find("(");
	if (pos == std::string::npos)
		return prettyFunction;           // something is not right
	prettyFunction.erase(pos);         // cut away signature
	pos = prettyFunction.rfind(" ");
	if (pos == std::string::npos)
		return prettyFunction;           // something is not right
	prettyFunction.erase(0, pos + 1);  // cut away return type
	return prettyFunction;
}

#define printErr  std::cerr << "!!! " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: error: "   << std::flush
#define printWarn std::cerr << "??? " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: warning: " << std::flush
#define printInfo std::cout << ">>> " << getClassMethod__(__PRETTY_FUNCTION__) << "(): info: "  << std::flush


//////////////////////////////////////////////////////////////////////////////
// functions to print version and compilation info

#ifndef SVN_VERSION  // SVN_VERSION set by Makefile
#define SVN_VERSION "undefined"
#endif
inline std::string svnVersion() { return SVN_VERSION; }

inline
void
printSvnVersion()
{
	const std::string ver = svnVersion();
	if (ver == "")
		printInfo << "subversion repository revision is unknown." << std::endl;
	else
		printInfo << "subversion repository revision is '" << ver << "'" << std::endl;
}


#ifndef CMAKE_SOURCE_DIR  // CMAKE_SOURCE_DIR set by Makefile
#define CMAKE_SOURCE_DIR "undefined"
#endif
inline std::string compileDir() { return CMAKE_SOURCE_DIR; }

inline
void
printCompilerInfo()
{
	const std::string date = __DATE__;
	const std::string time = __TIME__;
	const std::string ver  = __VERSION__;
	const std::string dir  = compileDir();
	printInfo << "this executable was compiled in ";
	if (dir != "")
		std::cout << "'" << dir << "'";
	else
		std::cout << "unknown directory";
	std::cout << " on " << date << " " << time << " by compiler " << ver << std::endl;
}


//////////////////////////////////////////////////////////////////////////////
// simple stream operators for some STL classes

template<typename T1, typename T2>
inline
std::ostream&
operator << (std::ostream&            out,
             const std::pair<T1, T2>& pair)
{
	return out << "(" << pair.first << ", " << pair.second << ")";
}
	
	
template<typename T>
inline
std::ostream&
operator << (std::ostream&         out,
             const std::vector<T>& vec)
{
	out << "{";
	for (unsigned int i = 0; i < (vec.size() - 1); ++i)
		out << "[" << i << "] = " << vec[i] << ", ";
	return out << "[" << vec.size() - 1 << "] = " << vec[vec.size() - 1] << "}";
}


#endif  // REPORTINGUTILS_H
