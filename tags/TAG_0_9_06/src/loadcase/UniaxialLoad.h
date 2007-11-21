/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 2, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License along   *
*   with this program; if not, write to the Free Software Foundation, Inc.,   *
*   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#ifndef _UNIAXIALLOAD_H
#define _UNIAXIALLOAD_H

#include <ElementalLoad.h>

/******************************************************************************
* Uniaxial Load
******************************************************************************/
class UniaxialLoad: public ElementalLoad
{
protected:
	Vector myDirection;
	Vector projections;
	double L;
public:
	UniaxialLoad();
	UniaxialLoad(int elemID,const char* dir);
	~UniaxialLoad();
};

/******************************************************************************
* Beam load point
******************************************************************************/
class BeamLoadPoint: public UniaxialLoad
{
private:
	double a0;
	double p0;
public:
	BeamLoadPoint();
	BeamLoadPoint(int elemID,const char* dir,double a0,double p0);
 	~BeamLoadPoint();
	const Vector& getP();
};

/******************************************************************************
* Beam load uniform
******************************************************************************/
class BeamLoadUniform: public UniaxialLoad
{
private:
	double p0;
public:
	BeamLoadUniform();
	BeamLoadUniform(int elemID,const char* dir,double p0);
 	~BeamLoadUniform();
    const Vector& getP();
};
#endif