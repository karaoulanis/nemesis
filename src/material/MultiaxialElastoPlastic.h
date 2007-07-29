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

#ifndef _MULTIAXIALELASTOPLASTIC_H
#define _MULTIAXIALELASTOPLASTIC_H

#include <MultiaxialMaterial.h>
#include <Surface.h>
/**
 * The MultiaxialElastoPlastic Class.
 */
class MultiaxialElastoPlastic: public MultiaxialMaterial
{
protected:
	MultiaxialMaterial* myElastic;

	static Matrix C;
	Vector eTrial;
	Vector ePTrial,ePConvg;
	Vector qTrial,qConvg;
	double aTrial,aConvg;
	bool plastic;
	int nHardeningVariables;

	std::vector<Surface*> fSurfaces;
	std::vector<Surface*> gSurfaces;
	inline std::vector<Surface*> getfSurfaces()		{return fSurfaces;}
	inline std::vector<Surface*> getgSurfaces()		{return gSurfaces;}
	void returnMapSYS(const Vector& De);
	void returnMapMYS(const Vector& De);
	void returnMapMYS2(const Vector& De);
	void returnMapTest(const Vector& De);
public:
	MultiaxialElastoPlastic();
	MultiaxialElastoPlastic(int ID,int elasticID);
	~MultiaxialElastoPlastic();
	
	void setStrain(const Vector& De);
	void commit();
	const Matrix& getC();
	bool isPlastic()							{return plastic;}
	
	// Tracker member functions
	void track();
};
#endif
