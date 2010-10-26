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

#ifndef _GROUNDMOTIONSIN_H
#define _GROUNDMOTIONSIN_H

#include <iostream>
#include <map>
#include <vector>
#include "elements/Element.h"
#include "loadcase/Load.h"

class GroundMotionSin: public Load
{
protected:
	int dof;
	double a;
	double omega;
	double phi;
public:
	GroundMotionSin();
	GroundMotionSin(int dof_,double a_,double omega_,double phi_=0.);

	void apply(double fact,double time);
};

#endif
