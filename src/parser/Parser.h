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

#ifndef _PARSER_H
#define _PARSER_H

#include "algorithm/bfgs.h"
#include "algorithm/linear_algorithm.h"
#include "algorithm/newton_raphson_full.h"
#include "algorithm/newton_raphson_initial.h"
#include "algorithm/newton_raphson_modified.h"
#include "analysis/Analysis.h"
#include "analysis/EigenAnalysis.h"
#include "analysis/SensitivityStaticAnalysis.h"
#include "analysis/StaticAnalysis.h"
#include "analysis/TransientAnalysis.h"
#include "analysis/XFemAnalysis.h"
#include "constraints/Constraint.h"
#include "control/ArcLengthSpherical.h"
#include "control/ArcLengthUNP.h"
#include "control/DisplacementControl.h"
#include "control/LoadControl.h"
#include "control/Newmark.h"
#include "crack/Crack.h"
#include "domain/Domain.h"
#include "elements/Bar2s.h"
#include "elements/Bar2t.h"
#include "elements/Beam2e.h"
#include "elements/Brick8b.h"
#include "elements/Brick8d.h"
#include "elements/Brick8i.h"
#include "elements/Quad4b.h"
#include "elements/Quad4d.h"
#include "elements/Quad4DispPlain.h"
#include "elements/Quad4DispAxisymmetric.h"
#include "elements/Quad4e.h"
#include "elements/Quad4i.h"
#include "elements/SDofElement.h"
#include "elements/Spring.h"
#include "elements/Tetrahedron4Disp.h"
#include "elements/Timoshenko2d.h"
#include "elements/Triangle3.h"
#include "elements/Triangle3XFem.h"
#include "elements/Triangle6.h"
#include "group/Group.h"
#include "imposer/EliminationImposer.h"
#include "imposer/PenaltyImposer.h"
#include "imposer/LagrangeImposer.h"
#include "loadcase/ElementSensitivityParameter.h"
#include "loadcase/GroundMotionFile.h"
#include "loadcase/GroundMotionSin.h"
#include "loadcase/InitialDisplacement.h"
#include "loadcase/InitialVelocity.h"
#include "loadcase/InitialStresses.h"
#include "loadcase/NodalLoadConstant.h"
#include "loadcase/NodalLoadLinear.h"
#include "loadcase/NodalLoadSin.h"
#include "loadcase/UniaxialLoad.h"
#include "material/Creep.h"
#include "material/DruckerPrager.h"
#include "material/DruckerPragerNew.h"
#include "material/DruckerPragerNew2.h"
#include "material/DruckerPragerNew3.h"
#include "material/DuncanChang.h"
#include "material/HoekBrown.h"
#include "material/LadeDuncan.h"
#include "material/ModifiedCamClay.h"
#include "material/MohrCoulomb.h"
#include "material/MultiaxialElastic.h"
#include "material/PlaneStress.h"
#include "material/SDofMaterial.h"
#include "material/SpringElastic.h"
#include "material/SpringContact.h"
#include "material/SpringMaterial.h"
#include "material/Tresca.h"
#include "material/UniaxialCyclic.h"
#include "material/UniaxialElastic.h"
#include "material/UniaxialElastoPlastic.h"
#include "material/UniaxialGap.h"
#include "material/VonMises.h"
#include "model/Model.h"
#include "node/Node.h"
#include "reorderer/ForwardCuthillMckee.h"
#include "reorderer/ForwardSloan.h"
//#include "reorderer/King.h"
//#include "reorderer/MinimumDegreeOrdering.h"
#include "reorderer/ReverseCuthillMckee.h"
#include "reorderer/ReverseSloan.h"
#include "reorderer/Reorderer.h"
#include "soe/FullLinearSOE.h"
#include "soe/SymmLinearSOE.h"
#include "soe/BandLinearSOE.h"

class Parser
{
protected:
	Domain D;
	Analysis A;    
public:
	// Constructors and destructor
	Parser();
	virtual ~Parser();

	// Parse the problem
	virtual int parse()=0;
	virtual int parse(char* filename)=0;
};

#endif
