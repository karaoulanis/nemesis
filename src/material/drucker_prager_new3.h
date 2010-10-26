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

#ifndef _DRUCKERPRAGERNEW3_H
#define _DRUCKERPRAGERNEW3_H

#include "material/MultiaxialMaterial.h"
#include "material/ys.h"
#include "material/Hardening.h"

/**
 * The Drucker-Prager Class.
 */
class DruckerPragerNew3: public MultiaxialMaterial
{
private:
	MultiaxialMaterial* myElastic;
	static Matrix C;
	static Matrix C3;
	bool plastic;
	int inaccurate;
	double aTrial,aConvg;

	vector<YS*> fSurfaces;
	vector<YS*> gSurfaces;
	Hardening EL;
public:
	DruckerPragerNew3();
	DruckerPragerNew3(int ID,int elasticID,double c,double phi,double psi,double Kci,double Kphi,double T);
	~DruckerPragerNew3();

	MultiaxialMaterial* getClone();
	void setStrain(const Vector& De);
	void commit();
	const Matrix& getC();
	bool isPlastic();

	// Tracker member functions
	void track();
};
#endif