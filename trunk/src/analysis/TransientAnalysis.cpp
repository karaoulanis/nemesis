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

#include <TransientAnalysis.h>

TransientAnalysis::TransientAnalysis()
	:AnalysisType()
{
	myTag=TAG_ANALYSIS_TRANSIENT;
}
bool TransientAnalysis::checkIfAllows(FEObject* f)
{
/*	if(	f->getTag()==TAG_CONTROL_GENERALIZED_A				||
		f->getTag()==TAG_ALGORITHM_LINEAR					||
		f->getTag()==TAG_ALGORITHM_NEWTON_RAPHSON_FULL		||
		f->getTag()==TAG_ALGORITHM_NEWTON_RAPHSON_MODIFED	||
		f->getTag()==TAG_ALGORITHM_NEWTON_RAPHSON_INITIAL	||
		f->getTag()==TAG_ALGORITHM_NEWTON_RAPHSON_PERIODIC	||
		f->getTag()==TAG_IMPOSER_ELIMINATION				||
		f->getTag()==TAG_IMPOSER_PENALTY					||
		f->getTag()==TAG_IMPOSER_LAGRANGE					||
		f->getTag()==TAG_SOE_FULL_GENERIC_POSITIVE_DEFINE	||
		f->getTag()==TAG_SOE_LINEAR_FULL					||
		f->getTag()==TAG_SOE_LINEAR_SYMM					||
		f->getTag()==TAG_SOE_LINEAR_BAND					||
		f->getTag()==TAG_CONVERGENCE_NORM_EUCLIDEAN			||
		f->getTag()==TAG_CONVERGENCE_NORM_MAXIMUM			||
		f->getTag()==TAG_REORDERER_NONE						||
		f->getTag()==TAG_REORDERER_FORWARD_CUTHILL_MCKEE	||
		f->getTag()==TAG_REORDERER_REVERSE_CUTHILL_MCKEE	||
		f->getTag()==TAG_REORDERER_FORWARD_SLOAN			||
		f->getTag()==TAG_REORDERER_REVERSE_SLOAN			 )
			return true;
	return false;*/
	return true;
}
int TransientAnalysis::run(int nLC,int nLoadSteps)
{
	// Check the imposer
	if(pA->getImposer()==0) 
		throw SException("[nemesis:%d] %s",9999,"No imposer has been set.");
	if(!this->checkIfAllows(pA->getImposer()))
		throw SException("[nemesis:%d] %s",9999,"The imposer type is incorrect.");
	// Check the control
	if(pA->getControl()==0) 
		throw SException("[nemesis:%d] %s",9999,"No control has been set.");
	if(!this->checkIfAllows(pA->getControl()))
		throw SException("[nemesis:%d] %s",9999,"The control type is incorrect.");
	// Check the algorithm
	if(pA->getAlgorithm()==0) 
		throw SException("[nemesis:%d] %s",9999,"No algorithm has been set.");
	if(!this->checkIfAllows(pA->getAlgorithm()))
		throw SException("[nemesis:%d] %s",9999,"The algorithm type is incorrect.");
	// Check the SOE
	if(pA->getSOE()==0) 
		throw SException("[nemesis:%d] %s",9999,"No soe has been set.");
	if(!this->checkIfAllows(pA->getSOE()))
		throw SException("[nemesis:%d] %s",9999,"The soe type is incorrect.");
	// Check the Norm
	if(pA->getConvergenceNorm()==0) 
		throw SException("[nemesis:%d] %s",9999,"No norm has been set.");
	if(!this->checkIfAllows(pA->getConvergenceNorm()))
		throw SException("[nemesis:%d] %s",9999,"The norm type is incorrect.");
	// Check the Reorderer
	if(pA->getReorderer()!=0&&(!this->checkIfAllows(pA->getReorderer())))
		throw SException("[nemesis:%d] %s",9999,"The reorderer type is incorrect.");

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
	pA->getDomain()->keepTrack(pA->getControl()->getLambda(),pA->getControl()->getTime());

	for(int i=0;i<nLoadSteps;i++)
	{
		// Call algorithm to solve step
		int check=pA->getAlgorithm()->solveStep(i);
		
		// Algorithm failed
		if(check<0)
		{
			if(check==-1)		cout<<"Warning  : Solution is diverging."<<endl;
			else if(check==-2)	cout<<"Warning  : Maximum number of iteration was exceeded."<<endl;
			pA->getControl()->rollback();
			break;
		}
		// Algorithm succeeded
		pA->getControl()->commit();
		pA->getDomain()->keepTrack(pA->getControl()->getLambda(),pA->getControl()->getTime());
	}
	pA->getDomain()->get<LoadCase>(pA->getDomain()->getLoadCases(),nLC)->commit();
	return 0;
}
