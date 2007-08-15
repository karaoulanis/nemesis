/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
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

#ifndef _ANALYSIS_H
#define _ANALYSIS_H

#include <Domain.h>
#include <Model.h>
#include <AnalysisType.h>
#include <Model.h>
#include <Imposer.h>
#include <Control.h>
#include <Algorithm.h>
#include <ConvergenceNorm.h>
#include <Reorderer.h>
#include <SOE.h>

class Model;
class Control;
class Imposer;
class SOE;
class Reorderer;
class AnalysisType;
class Algorithm;
class ConvergenceNorm;

class Analysis
{
private:
	Model M;
	Domain* theDomain;
	AnalysisType* theAnalysisType;
	Algorithm* theAlgorithm;
	Control* theControl;
	SOE* theSOE;
	Imposer* theImposer;
	ConvergenceNorm* theNorm;
	Reorderer* theReorderer;
public:
	Analysis(Domain* pDomain);
	~Analysis();
	
	inline Model* getModel()						{return &M;}
	inline AnalysisType* getAnalysisType()			{return theAnalysisType;}
	inline Algorithm* getAlgorithm()				{return theAlgorithm;}
	inline Control* getControl()					{return theControl;}
	inline Imposer* getImposer()					{return theImposer;}
	inline ConvergenceNorm* getConvergenceNorm()	{return theNorm;}
	inline Reorderer* getReorderer()				{return theReorderer;}
	inline SOE* getSOE()							{return theSOE;}
	inline Domain* getDomain()						{return theDomain;}

	void setAnalysisType(AnalysisType* p);
	void setAlgorithm(Algorithm* p);
	void setControl(Control* p);
	void setImposer(Imposer* p);
	void setReorderer(Reorderer* p);
	void setSOE(SOE* p);

	int analyze(int lcID,int nLoadSteps);
	void clear();
};
#endif
