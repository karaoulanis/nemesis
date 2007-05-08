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

#ifndef _PARSER_H
#define _PARSER_H

#include <Domain.h>
#include <Analysis.h>
#include <Model.h>
#include <Node.h>
#include <Group.h>
#include <UniaxialElastic.h>
#include <UniaxialElastoPlastic.h>
#include <UniaxialCyclic.h>
#include <MultiaxialElastic.h>
#include <VonMises.h>
#include <MohrCoulomb.h>
#include <Bar2s.h>
#include <Bar2t.h>
#include <Beam2e.h>
//#include <Beam2t.h>
//#include <Beam2t3d.h>
#include <SDofElement.h>
#include <Brick8Disp.h>
#include <Quad4DispPlain.h>
#include <Quad4DispAxisymmetric.h>
#include <Triangle3.h>
#include <Tetrahedron4Disp.h>
#include <Constraint.h>
#include <NodalLoadConstant.h>
#include <NodalLoadLinear.h>
#include <NodalLoadSin.h>
#include <UniaxialLoad.h>
#include <ElementSensitivityParameter.h>
#include <StaticAnalysis.h>
#include <TransientAnalysis.h>
#include <EigenAnalysis.h>
#include <SensitivityStaticAnalysis.h>
#include <EliminationImposer.h>
#include <PenaltyImposer.h>
#include <LagrangeImposer.h>
#include <ArcLengthSpherical.h>
#include <ArcLengthUNP.h>
#include <DisplacementControl.h>
#include <LoadControl.h>
#include <Newmark.h>
#include <LinearAlgorithm.h>
#include <BFGS.h>
#include <NewtonRaphsonFull.h>
#include <NewtonRaphsonModified.h>
#include <NewtonRaphsonInitial.h>
#include <FullLinearSOE.h>
#include <SymmLinearSOE.h>
#include <BandLinearSOE.h>
#include <Tracker.h>
#include <InitialDisplacement.h>
#include <InitialVelocity.h>
#include <InitialStresses.h>
#include <GroundMotionFile.h>
#include <GroundMotionSin.h>
#include <Reorderer.h>
#include <ReverseCuthillMckee.h>
#include <ForwardCuthillMckee.h>
#include <ReverseSloan.h>
#include <ForwardSloan.h>
//#include <King.h>
//#include <MinimumDegreeOrdering.h>

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
