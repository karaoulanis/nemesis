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

#ifndef _MODELOBJECT_H
#define _MODELOBJECT_H

#include <Element.h>

class ModelObject
{
protected:
	IDContainer theFTable;
	static Matrix** theStaticMatrices;
	static Vector** theStaticVectors;
	Matrix* myMatrix;
	Vector* myVector;
public:
	// Constructors
	ModelObject();
	ModelObject(const IDContainer& FTable);
	virtual ~ModelObject();
	
	// Access to data members
	const IDContainer& getFTable() const;
	void setFTable(const IDContainer& FTable);
	void setFTable(int index,int val);

	inline void zeroMatrix()				{myMatrix->clear();}
	inline void zeroVector()				{myVector->clear();}
	inline const Matrix& getMatrix() const	{return *myMatrix; }
	inline const Vector& getVector() const	{return *myVector; }
};

#endif
