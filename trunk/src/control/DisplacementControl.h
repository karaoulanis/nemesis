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

#ifndef _DISPLACEMENTCONTROL_H
#define _DISPLACEMENTCONTROL_H

#include "control/StaticControl.h"

/**
 * The DisplacementControl Class.
 * DisplacementControl is a static control that constraints the displacement
 * increment at a particular variable \f$k\f$ to a fixed quantity 
 * \f$\Delta k\f$. For more details see p.275-276 of "Non-linear Finite Element 
 * Analysis of Solids and Strucures", Vol.1.\n
 * Auto incrementation is done here similarly to that used in LoadControl, i.e.
 * \f[\Delta u_{n}=\Delta u_{0}
 *						\left(
 *						\frac{I_d}{I_o}
 *                      \right)^n
 * \f]
 * \f$I_d\f$ is the desired number of
 * iterations within each step (Crisfield suggests ~3), \f$I_o\f$ is the number
 * of iterations in the last step and \f$n\f$ is an exponent, usually set to 
 * 0.5 as suggested by Ramm.\n
 * \f$\Delta u\f$ is also limited within min\f$\Delta u\f$ and
 * max\f$\Delta u\f$. By setting those equal to \f$\Delta u_0\f$, then no 
 * auto-incrementation takes place and \f$\Delta u=\Delta u_0\f$. 
 */
class DisplacementControl :public StaticControl
{
private:
	double DeltaU;				///< Current Delta u for this step
	double Du_k;				///< Accumulative Delta u at k
	double duT_k;				///< Tangent delta u at k
	int theNodeID;				///< The Node that determines the Du
	int theDofID;				///< The local Dof that determines the Du
	int theRefDof;				///< From Node and Dof comes the global dof
public:
	DisplacementControl(int nodeID,int dofID,
		double Du0,double minDu,double maxDu,int IterDesired,double n,double DeltaTime);
	~DisplacementControl();

	// Methods for incremental/iterative algorithms
	void predict();
	void correct();
};

#endif
