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

#ifndef _CONTROL_H
#define _CONTROL_H

#include "analysis/analysis.h"
#include "analysis/analysis_object.h"
#include "model/ModelElement.h"
#include "model/ModelNode.h"

class Control: public AnalysisObject
{
protected:
	double lambdaTrial;
	double lambdaConvg;
	double DLambda;
	double dLambda;
	Vector qRef;
public:
	Control();
	virtual ~Control();

	virtual void formTangent();    
	virtual void formResidual(double factor)=0;

	// Functions to build Element by Element or Node by Node it's contribution
	virtual void formElementalTangent(ModelElement* pModelElement)=0;
	virtual void formElementalResidual(ModelElement* pModelElement,double time=0.)=0;
	virtual void formNodalResidual(ModelNode* pModelNode)=0;

	virtual double getLambda();
	virtual double getTime() {return 0;} ///@todo: implement this better
	
	virtual void init()=0;
	virtual void predict()=0;
	virtual void correct()=0;
	virtual void commit()=0;
	virtual void rollback()=0;
};
#endif
