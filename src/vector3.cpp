// vector3.cpp
/*
 * $Id: vector3.cpp,v 1.0 2010/07/04  $
 *
 *
 * /author 
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
 * 
 * $Log: $
 */

#include "vector3.h"

Vector3::Vector3() 
{
   fVec[0] = 0;
   fVec[1] = 0;
   fVec[2] = 0;
}

Vector3::Vector3(double x, double y, double z)
{
   fVec[0] = x;
   fVec[1] = y;
   fVec[2] = z;
}

Vector3::~Vector3()
{
}

void Vector3::SetVector(double x, double y, double z)
{
   fVec[0] = x;
   fVec[1] = y;
   fVec[2] = z;
}

void Vector3::SetVector(double *vec)
{
   fVec[0] = vec[0];
   fVec[1] = vec[1];
   fVec[2] = vec[2];
}
