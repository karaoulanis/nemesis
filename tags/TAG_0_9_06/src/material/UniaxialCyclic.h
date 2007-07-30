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

#ifndef _UNIAXIALCYCLIC_H
#define _UNIAXIALCYCLIC_H

#include <UniaxialMaterial.h>

class UniaxialCyclic: public UniaxialMaterial
{
private:
	double sr;
	double er;
	double eTrial;
	double eConvg;
public:
	UniaxialCyclic();
	UniaxialCyclic(int ID,double E,double nu,double rho,double aT,double tmax,double Gmax);
	~UniaxialCyclic();
	UniaxialMaterial* getClone();
	void setStrain(const double De);
	const double getC();
	void commit();
};
#endif
