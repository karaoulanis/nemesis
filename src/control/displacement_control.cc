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

#include "control/displacement_control.h"
#include <cmath>

/**
 * Constructor.
 * \param nodeID		Reference node.
 * \param dofID			Reference dof.
 * \param Du0			Initial Delta u.
 * \param minDu			Lower bound for Delta u.
 * \param maxDu			Upper bound for Delta u.
 * \param IterDesired	Desired number of iterations.
 * \param n				Exponent parameter.
 * \param DeltaTime     Timestep for viscoplastic solutions
*/
DisplacementControl::DisplacementControl(int nodeID,int dofID,
		double Du0,double minDu,double maxDu,int IterDesired,double n,double DeltaTime)
:StaticControl(Du0,minDu,maxDu,IterDesired,n,DeltaTime),DeltaU(Du0)
{
	theNodeID=nodeID;
	theDofID=dofID;
	myTag=TAG_CONTROL_DISPLACEMENT;
}
/**
 * Destructor.
 */
DisplacementControl::~DisplacementControl()
{
	// Does nothing
}
/**
 * Creates new step for a Displacement control based static analysis.
 * It does two things:
 * \li Auto-incrementation \n
 * \f$\Delta u(k)\f$ at each new step is given as
 * \f[\Delta u_{n}(k)=\Delta u_{0}(k)
 *						\left(
 *						\frac{I_d}{I_o}
 *                      \right)^n
 * \f]
 * \f$I_d\f$ is the desired number of iterations within each step (Crisfield 
 * suggests ~3), \f$I_o\f$ is the number of iterations in the last step and 
 * \f$n\f$ is an exponent, usually set to 0.5 as suggested by Ramm.\n
 * \f$\Delta u(k)\f$ is also limited within min\f$\Delta l\f$ and
 * max\f$\Delta u(k)\f$. By setting those equal to \f$\Delta u(k)_{0}\f$,
 * then no auto-incrementation takes place and \f$\Delta l=\Delta u(k)_{0}\f$. 
 * It should be noted that min\f$\Delta u(k)\f$ and min\f$\Delta u(k)\f$ are 
 * given as absolute values.
 * \li Predictor Step\n
 * The predictor step is based on a forward Euler scheme. It is based on 
 * p.276 (eq.9.34) of Cridfield's book. As soon this is computed the domain is 
 * updated.
 */
void DisplacementControl::predict()
{
	// Set the reference dof
	///@todo Since changes in the dof schemes have occurred, check if this is ok.
	theRefDof=pA->getModel()->getSOEDof(theNodeID,theDofID);
	
	// Find DeltaU increment
	///@todo Auto-incrementation involves abs() and this might be a problem...
//	DeltaU*=pow(((double)Id/(double)Io),nExp);
//	if(DeltaU<minDelta) DeltaU=minDelta;
//	else if (DeltaU>maxDelta) DeltaU=maxDelta;
	
	// Find duT (tangent) at coefficient k
//	theSOE->print();
	pA->getSOE()->setB(qRef);
	pA->getSOE()->solve();
	duT=pA->getSOE()->getX();
	duT_k=duT[theRefDof];
	
	// Now DLambda can be found, see Crisfield, (9.34). 
	DLambda=DeltaU/duT_k;
	lambdaTrial+=DLambda;

	// Find displacement vectors and update the model
	du=DLambda*duT;
	Du_k=du[theRefDof];
	pA->getModel()->incTrialDisp(du);
	pA->getModel()->update();

	// Set num of achieved iterations to one
	Io=1;

	// Set du to the SOE so the norm can find it there
	pA->getSOE()->setX(du);
}
/**
 * Updates that occur within each iterative step.
 * \f$\delta\lambda\f$ is given in (9.32) p.276, Crisfield, Vol.1.  
  */
void DisplacementControl::correct()
{
	// Find du_bar
	du=pA->getSOE()->getX();
	double duBar_k=du[theRefDof];

	// Find du_t_k
	pA->getSOE()->setB(qRef);
	pA->getSOE()->solve();
	duT=pA->getSOE()->getX();
	duT_k=duT[theRefDof];

	// Find dLambda and update quantities in the control
	dLambda=(DeltaU-Du_k-duBar_k)/duT_k;
	lambdaTrial+=dLambda;
	Du_k+=duBar_k+dLambda*duT_k;

	// Update displacements in the model
	du+=dLambda*duT;
	pA->getModel()->incTrialDisp(du);
	pA->getModel()->update();

	// Increase number of iterations
	Io++;
//	pA->getSOE()->setX(du);
}
