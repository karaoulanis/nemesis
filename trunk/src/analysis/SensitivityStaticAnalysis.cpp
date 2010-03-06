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
//            C.G. Panagiotopoulos (pchr@civil.auth.gr)
//*****************************************************************************

#include <SensitivityStaticAnalysis.h>

SensitivityStaticAnalysis::SensitivityStaticAnalysis()
	:AnalysisType()
{
	myTag=TAG_ANALYSIS_STATIC;
	// defaults
	pA->setImposer(new EliminationImposer());
	pA->setControl(new LoadControl(1.,1.,1.,1,0.5,0.));
	pA->setAlgorithm(new LinearAlgorithm());
	pA->setSOE(new FullLinearSOE());
	theSensitivityControl=new SensitivityControl;
}
SensitivityStaticAnalysis::~SensitivityStaticAnalysis()
{
	delete theSensitivityControl;
}
bool SensitivityStaticAnalysis::checkIfAllows(FEObject* f)
{
	return false;
}
int SensitivityStaticAnalysis::run(int nLC,int nLoadSteps)
{

	// Create model by applying the constraints
	pA->getImposer()->impose();

	// Now that model is complete, reorder the model
	if(pA->getReorderer()!=0) pA->getReorderer()->reorder();

	// Now that model is complete, the SOE can be initialized
	pA->getSOE()->setTheSize();
	
	// Apply the loads carried by the given loadcase
//	pA->getControl()->
///		setTheLoadCase(pA->getDomain()->get<LoadCase>(pA->getDomain()->getLoadCases(),nLC));
//	theSensitivityControl->setTheLoadCase(pA->getDomain()->get<LoadCase>(pA->getDomain()->getLoadCases(),nLC));

	// Initialize control
	pA->getControl()->init();
	theSensitivityControl->init();

	// Initialize the convergence check
	pA->getConvergenceNorm()->init(nLC,nLoadSteps);

	int ret=0;
	for(int i=0;i<nLoadSteps;i++)
	{
		// Call algorithm to solve step
		int check=pA->getAlgorithm()->solveStep(i);
		// Algorithm failed
		if(check<0)
		{
			if(check==-1)		cout<<"Warning  : Solution is diverging."<<endl;
			else if(check==-2)	cout<<"Warning  : Maximum number of iteration was exceeded."<<endl;
			ret=-1;
		}
		// Algorithm succeeded
		pA->getControl()->commit();
		
		// Sensitivity
		int nParams=pA->getDomain()->get<LoadCase>(pA->getDomain()->getLoadCases(),nLC)->getnSensitivityParameters();
		for(int j=0;j<nParams;j++)
		{
			theSensitivityControl->formTangent();
			theSensitivityControl->formResidual(0.);
			//pA->getSOE()->print();
			pA->getSOE()->solve();
			//cout<<pA->getSOE()->getX();
			theSensitivityControl->commit();
		}
	}

	// Finalize loadcase
	pA->getModel()->setNodalStress();
	return ret;
}
