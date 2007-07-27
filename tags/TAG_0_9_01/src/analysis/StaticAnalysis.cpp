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

#include <StaticAnalysis.h>

StaticAnalysis::StaticAnalysis()
	:AnalysisType()
{
	myTag=TAG_ANALYSIS_STATIC;
	// defaults
	pA->setImposer(new EliminationImposer());
	pA->setControl(new LoadControl(1.,1.,1.,1,0.5));
	pA->setAlgorithm(new LinearAlgorithm());
	pA->setSOE(new FullLinearSOE());
}
bool StaticAnalysis::checkIfAllows(FEObject* f)
{
	if(f->getTag()==TAG_CONTROL_LOAD						||
					TAG_CONTROL_DISPLACEMENT				||
					TAG_CONTROL_ARC_LENGTH_SPHERICAL		||
					TAG_CONTROL_ARC_LENGTH_UNP				||
					TAG_ALGORITHM_LINEAR					||
					TAG_ALGORITHM_NEWTON_RAPHSON_FULL		||
					TAG_ALGORITHM_NEWTON_RAPHSON_MODIFED	||
					TAG_ALGORITHM_NEWTON_RAPHSON_INITIAL	||
					TAG_ALGORITHM_NEWTON_RAPHSON_PERIODIC	||
					TAG_IMPOSER_PENALTY						||
					TAG_IMPOSER_LAGRANGE					||
					TAG_SOE_FULL_GENERIC_POSITIVE_DEFINE	||
					TAG_SOE_LINEAR_FULL						||
					TAG_SOE_LINEAR_SYMM						||
					TAG_SOE_LINEAR_BAND						||
					TAG_CONVERGENCE_NORM_EUCLIDEAN			||
					TAG_CONVERGENCE_NORM_MAXIMUM			||
					TAG_REORDERER_NONE						||
					TAG_REORDERER_FORWARD_CUTHILL_MCKEE		||
					TAG_REORDERER_REVERSE_CUTHILL_MCKEE		||
					TAG_REORDERER_FORWARD_SLOAN				||
					TAG_REORDERER_REVERSE_SLOAN				 )

			return true;
	return false;
}
int StaticAnalysis::run(int nLC,int nLoadSteps)
{
	// Check the imposer
	if(pA->getImposer()==0) 
		throw SolverException(9999,"No imposer has been set.");
	if(!this->checkIfAllows(pA->getImposer()))
		throw SolverException(9999,"The imposer type is incorrect.");
	// Check the control
	if(pA->getControl()==0) 
		throw SolverException(9999,"No control has been set.");
	if(!this->checkIfAllows(pA->getControl()))
		throw SolverException(9999,"The control type is incorrect.");
	// Check the algorithm
	if(pA->getAlgorithm()==0) 
		throw SolverException(9999,"No algorithm has been set.");
	if(!this->checkIfAllows(pA->getAlgorithm()))
		throw SolverException(9999,"The algorithm type is incorrect.");
	// Check the SOE
	if(pA->getSOE()==0) 
		throw SolverException(9999,"No soe has been set.");
	if(!this->checkIfAllows(pA->getSOE()))
		throw SolverException(9999,"The soe type is incorrect.");
	// Check the Reorderer
	if(pA->getReorderer()!=0&&(!this->checkIfAllows(pA->getReorderer())))
		throw SolverException(9999,"The reorderer type is incorrect.");

	// Create model by applying the constraints
	pA->getImposer()->impose();

	// Now that model is complete, reorder the model
	if(pA->getReorderer()!=0) pA->getReorderer()->reorder();

	// Now that model is complete, the SOE can be initialized
	pA->getSOE()->setTheSize();
	
	// Initialize
	pA->getDomain()->get<LoadCase>(pA->getDomain()->getLoadCases(),nLC)->init();
	pA->getControl()->init();
	pA->getConvergenceNorm()->init(nLC,nLoadSteps);
	pA->getDomain()->keepTrack(pA->getControl()->getLambda(),pA->getControl()->getLambda());

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
//			pA->getControl()->returnToConverged();
			pA->getControl()->commit();
//			pA->getDomain()->keepTrack(pA->getControl()->getLambda(),pA->getControl()->getLambda());
//			break;
			ret=-1;
		}
		// Algorithm succeeded
		pA->getControl()->commit();
		pA->getDomain()->keepTrack(pA->getControl()->getLambda(),pA->getControl()->getLambda());
	}
	// Finalize
	pA->getDomain()->get<LoadCase>(pA->getDomain()->getLoadCases(),nLC)->commit();
	pA->getModel()->setNodalStress();
	return ret;
}