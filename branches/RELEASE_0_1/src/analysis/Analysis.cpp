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

#include <Analysis.h>

Analysis::Analysis(Domain* pDomain)
:M(pDomain),theDomain(pDomain)
{
	M.setAnalysis(this);
	theNorm=new ConvergenceNorm();
	theAnalysisType=0;
	theImposer=0;
	theControl=0;
	theAlgorithm=0;
	theSOE=0;
	theReorderer=0;
}
Analysis::~Analysis()
{
	delete theNorm;
	this->clear();
}
int Analysis::analyze(int lcID,int nLoadSteps)
{
	if(theAnalysisType==0) throw SolverException(9999,"No analysis type set.");
	// Run the analysis
	return theAnalysisType->run(lcID,nLoadSteps);
}
void Analysis::clear()
{
	M.clear();
	if(theAnalysisType!=0)	{delete theAnalysisType;	theAnalysisType=0;}
	if(theImposer!=0)		{delete theImposer;			theImposer=0;}
	if(theControl!=0)		{delete theControl;			theControl=0;}
	if(theAlgorithm!=0)		{delete theAlgorithm;		theAlgorithm=0;}
	if(theSOE!=0)			{delete theSOE;				theSOE=0;}
	if(theReorderer!=0)		{delete theReorderer;		theReorderer=0;}
}
void Analysis::setAnalysisType(AnalysisType* p)
{
	if(theAnalysisType!=0) delete theAnalysisType;
	theAnalysisType=p;
}
void Analysis::setAlgorithm(Algorithm* p)
{
	if(theAlgorithm!=0) delete theAlgorithm;
	theAlgorithm=p;
}
void Analysis::setControl(Control* p)
{
	if(theControl!=0) delete theControl;
	theControl=p;
}
void Analysis::setImposer(Imposer* p)
{
	if(theImposer!=0) delete theImposer;
	theImposer=p;
}
void Analysis::setReorderer(Reorderer* p)
{
	if(theReorderer!=0) delete theReorderer;
	theReorderer=p;
}
void Analysis::setSOE(SOE* p)
{
	if(theSOE!=0) delete theSOE;
	theSOE=p;
}
