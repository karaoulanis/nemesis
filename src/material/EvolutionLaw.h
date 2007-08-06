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

#ifndef _EVOLUTIONLAW_H
#define _EVOLUTIONLAW_H

#include <Vector.h>

/**
 * The Evolution Law Class.
 */
class EvolutionLaw
{
private:
public:
	EvolutionLaw();
	virtual ~EvolutionLaw();
	
	virtual const double get_h(const Vector& v)=0;
	virtual const double get_dhds(const Vector& sTrial,const Vector& ePTrial)=0;
	virtual const double get_dhda(const Vector& sTrial,const Vector& ePTrial)=0;
};
#endif
