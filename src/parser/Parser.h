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
#include "analysis/analysis.h"
#include "analysis/eigen_analysis.h"
#include "analysis/sensitivity_static_analysis.h"
#include "analysis/static_analysis.h"
#include "analysis/transient_analysis.h"
#include "analysis/xfem_analysis.h"
#include "constraints/constraint.h"
#include "control/arc_length_spherical.h"
#include "control/arc_length_unp.h"
#include "control/displacement_control.h"
#include "control/load_control.h"
#include "control/newmark.h"
#include "crack/crack.h"
#include "domain/domain.h"
#include "elements/bar2s.h"
#include "elements/bar2t.h"
#include "elements/beam2e.h"
#include "elements/brick8b.h"
#include "elements/brick8d.h"
#include "elements/brick8i.h"
#include "elements/quad4b.h"
#include "elements/quad4d.h"
#include "elements/quad4_disp_plain.h"
#include "elements/quad4_disp_axisymmetric.h"
#include "elements/quad4e.h"
#include "elements/quad4i.h"
#include "elements/sdof_element.h"
#include "elements/spring.h"
#include "elements/tetrahedron4_disp.h"
#include "elements/timoshenko2d.h"
#include "elements/triangle3.h"
#include "elements/triangle3_xfem.h"
#include "elements/triangle6.h"
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
