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

#ifndef NEMESIS_MATERIAL_LADE_DUNCAN_H_
#define NEMESIS_MATERIAL_LADE_DUNCAN_H_

#include "material/multiaxial_elastic_plastic.h"

/**
 * The Lade-Duncan Material Class.
 */
class LadeDuncan: public MultiaxialElastoPlastic
{
private:
public:
	LadeDuncan();
	LadeDuncan(int ID,int elasticID,double K);
	MultiaxialMaterial* getClone();
	~LadeDuncan();
};
#endif //NEMESIS_MATERIAL_LADE_DUNCAN_H_
