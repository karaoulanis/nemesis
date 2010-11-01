/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#ifndef NEMESIS_MATERIAL_DRUCKER_PRAGER_NEW_H_
#define NEMESIS_MATERIAL_DRUCKER_PRAGER_NEW_H_

#include "material/multiaxial_material.h"
#include "material/ys.h"
#include "material/hardening.h"

/**
 * The Drucker-Prager Class.
 */
class DruckerPragerNew: public MultiaxialMaterial
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
	DruckerPragerNew();
	DruckerPragerNew(int ID,int elasticID,double c,double phi,double psi,double Kci,double Kphi,double T);
	~DruckerPragerNew();

	MultiaxialMaterial* getClone();
	void setStrain(const Vector& De);
	void commit();
	const Matrix& getC();
	bool isPlastic();

	// Tracker member functions
	void track();
};
#endif //NEMESIS_MATERIAL_DRUCKER_PRAGER_NEW_H_
