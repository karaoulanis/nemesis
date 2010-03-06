/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 3, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#ifndef _MOHRCOULOMB_H
#define _MOHRCOULOMB_H

#include <MultiaxialMaterial.h>

/**
 * The Mohr-Coulomb Class.
 */
class MohrCoulomb: public MultiaxialMaterial
{
private:
	MultiaxialMaterial* myElastic;
	static Matrix C;
	static Matrix C3;
	bool plastic;
	int inaccurate;
public:
	MohrCoulomb();
	MohrCoulomb(int ID,int elasticID,double c,double phi,double alpha);
	~MohrCoulomb();

	MultiaxialMaterial* getClone();
	void setStrain(const Vector& De);
	void commit();
	const Matrix& getC();
	bool isPlastic();

	// Tracker member functions
	void track();
};
#endif
