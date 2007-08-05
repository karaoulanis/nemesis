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

#ifndef _MULTIAXIALELASTIC_H
#define _MULTIAXIALELASTIC_H

#include <MultiaxialMaterial.h>

/**
 * The Elastic Class.
 */
class MultiaxialElastic: public MultiaxialMaterial
{
public:
	MultiaxialElastic();
	MultiaxialElastic(int ID,double E,double nu,double rho,double aT);
	MultiaxialMaterial* getClone();
	void setStrain(const Vector& De);
	const Matrix& getC();
	void commit();
	// Tracker member functions
	void track();
};
#endif
