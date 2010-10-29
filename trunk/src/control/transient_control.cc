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

#include "control/transient_control.h"

/**
 * Constructor.
 */
TransientControl::TransientControl()
{
}
/**
 * Destructor.
 */
TransientControl::~TransientControl()
{
}
/**
 * @return 0 if everything is ok.
 */
void TransientControl::init()
{
	int size=pA->getModel()->getnEquations();
	u.resize(size);			u.clear();
	v.resize(size);			v.clear();
	a.resize(size);			a.clear();
	ut.resize(size);		ut.clear();
	vt.resize(size);		vt.clear();
	at.resize(size);		at.clear();
	lambdaTrial=1.0;
	lambdaConvg=1.0;
/*
	//-------------------------------------------------------------------------
	// Apply initial conditions and gather displacements and accelerations
	//-------------------------------------------------------------------------
	theLoadCase->initialize();
	int i,n;
	n=pA->getModel()->getModelNodes().size();
	for(i=0;i<n;i++)
	{
		ModelNode* pModelNode=pA->getModel()->getModelNodes()[i];
		IDContainer FTable=pModelNode->getFTable();
		// Copy displacements
		pModelNode->zeroVector();
		pModelNode->add_uTrial();
		for(unsigned j=0;j<FTable.size();j++)
			if((FTable[j])>=0)
				u[FTable[j]]=pModelNode->getVector()[j];
		// Copy velocities
		pModelNode->zeroVector();
		pModelNode->add_vTrial();
		for(unsigned j=0;j<FTable.size();j++)
			if((FTable[j])>=0) 
				v[FTable[j]]=pModelNode->getVector()[j];
	}
	//-------------------------------------------------------------------------
	// Apply loads
	//-------------------------------------------------------------------------
	pA->getDomain()->zeroLoads();
	theLoadCase->setFactor(1.0);
	theLoadCase->applyLoads(0.);

	//-------------------------------------------------------------------------
	// Build the system of equations M.a=R-K.u-C.v
	//-------------------------------------------------------------------------
	pA->getSOE()->zeroA();
	pA->getSOE()->zeroB();
	// Take contribution from Elements
	n=pA->getModel()->getModelElements().size();
	for(i=0;i<n;i++)
	{
		ModelElement* pModelElement=pA->getModel()->getModelElements()[i];
		pModelElement->zeroMatrix();
		pModelElement->add_M();
		pA->getSOE()->insertMatrixIntoA(pModelElement->getMatrix(),
			pModelElement->getFTable(),1.0);

		pModelElement->zeroVector();
//		pModelElement->add_KuTrial(1.0);
//		pModelElement->add_CvTrial(1.0);
		pModelElement->add_Reff(1.0);
		pA->getSOE()->insertVectorIntoB(pModelElement->getVector(),
			pModelElement->getFTable(),-1.0);
	}
	// Take contribution from Nodes
	n=pA->getModel()->getModelNodes().size();
	for(i=0;i<n;i++)
	{
		ModelNode* p=pA->getModel()->getModelNodes()[i];
		this->formNodalResidual(p);
		pA->getSOE()->insertVectorIntoB(p->getVector(),p->getFTable(),-1.0);
	}
	//-------------------------------------------------------------------------
	// Solve the system, update accelerations and commit
	//-------------------------------------------------------------------------
	pA->getSOE()->solve();
	a=pA->getSOE()->getX();
	// Only accelerations are updated here, u and v are already updated from db
	pA->getModel()->setTrialVecs(u,v,a);
	this->commit();
*/
}
/**
 * Forms the elemental tangent for a static analysis, element by element.
 * @param pModelElement	A pointer to the ModelElement that is treated.
 * @return 0 if everything is ok.
 */
void TransientControl::formElementalTangent(ModelElement* pModelElement) 
{
	pModelElement->zeroMatrix();
	pModelElement->add_K(c[0]);
	pModelElement->add_M(c[1]);
	pModelElement->add_C(c[2]);
}
void TransientControl::formElementalResidual(ModelElement* pModelElement,double /*time*/)
{
	pModelElement->zeroVector();
	pModelElement->add_Reff(1.0);
}
void TransientControl::formNodalResidual(ModelNode* pModelNode)
{
	pModelNode->zeroVector();
	pModelNode->add_R(1.0);
}
void TransientControl::commit()
{
	pA->getModel()->commit();
}
/**
 * Aborts step.
 */
void TransientControl::rollback()
{
	lambdaTrial=lambdaConvg;
}
void TransientControl::formResidual(double fac)	
{
	pA->getSOE()->zeroB();
	pA->getDomain()->zeroLoads();
	pA->getDomain()->applyLoads(fac,pA->getDomain()->getTimeCurr());

	// Take contribution from Nodes
	for(unsigned i=0;i<pA->getModel()->getModelNodes().size();i++)		
	{
		ModelNode* p=pA->getModel()->getModelNodes()[i];
		this->formNodalResidual(p);
		pA->getSOE()->insertVectorIntoB(p->getVector(),p->getFTable(),-1.0);
	}
	// Take contribution from Elements
	for(unsigned i=0;i<pA->getModel()->getModelElements().size();i++)		
	{
		ModelElement* p=pA->getModel()->getModelElements()[i];
		this->formElementalResidual(p);
		pA->getSOE()->insertVectorIntoB(p->getVector(),p->getFTable(),-1.0);
	}
}
