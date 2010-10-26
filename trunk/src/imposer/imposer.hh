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

#ifndef _IMPOSER_H
#define _IMPOSER_H

#include "analysis/analysis.h"
#include "analysis/analysis_object.h"
#include "model/Model.h"
#include "model/StandardModelElement.h"
#include "model/StandardModelNode.h"

// Forward declarations
class Model;
class Domain;

class Imposer: public AnalysisObject
{
protected:
	Model* theModel;
	Domain* theDomain;
	ConstraintContainer* theConstraints;
	IDContainer myNodalIDs;
	IDContainer theNodalGlobalDofs;
public:
	Imposer();
	virtual ~Imposer();
	int createGlobalDofNumbering();
	int getGlobalDof(int NodeID,int localDof);
	const IDContainer getGlobalDofs(int NodeID);
	virtual int impose()=0;

};
#endif
