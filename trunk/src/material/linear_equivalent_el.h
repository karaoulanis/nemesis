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

#ifndef NEMESIS_MATERIAL_LINEAR_EQUIVALENT_EL_H_
#define NEMESIS_MATERIAL_LINEAR_EQUIVALENT_EL_H_

#include "material/evolution_law.h"

/**
 * The Linear Equivalent Evolution Law Class.
 */
class LinearEquivalentEL: public EvolutionLaw
{
private:
public:
	LinearEquivalentEL();
	~LinearEquivalentEL();
	
	double get_h(const Vector& v);
	double get_dhds(const Vector& sTrial,const Vector& ePTrial);
	double get_dhda(const Vector& sTrial,const Vector& ePTrial);
};
#endif //NEMESIS_MATERIAL_LINEAR_EQUIVALENT_EL_H_
