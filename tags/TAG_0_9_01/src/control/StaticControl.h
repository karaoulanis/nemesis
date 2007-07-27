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

#ifndef _STATICCONTROL_H
#define _STATICCONTROL_H

#include <Control.h>

/**
 * The StaticControl Class.
 * Static Control is a class derived from the Control class. It should
 * control the model through a static analysis, i.e.
 * \li 1. It should form for each ModelElement and for each ModelNode its
 * tangent matrix and its residual vector and it should insert those into the 
 * SOE. I.e. It must create the SOE for a static analysis.
 * \li 2. It should create a predictor step. (new step)
 * \li 3. It should be able to create corrective steps (update).
 * \li 4. It should update the model with the iterative displacements.
 * \li 5. It should commit the converged displacements to the model after 
 * each successive increment.
 * \li 6. It should return to the previous converged step, if the current 
 * step fails.
 */
class StaticControl :public Control
{
protected:
	double Delta0;				///< Initial Delta
	double minDelta;			///< Lower bound for Delta (absolute);
	double maxDelta;			///< Upper bound for Delta (absolute);
	Vector du;					///< Iterative displacements
	Vector duT;					///< SOE solution for q        
	Vector duBar;				///< SOE solution for lambda*q
	Vector Du;					///< Accumulative solution within each step
	int Io;						///< Number of iterations in the last step
	int Id;						///< Desired number of iterations in this step
	double nExp;				///< Exponent for the auto incrementation				
public:

	// Constructors and destructor
	StaticControl();
	StaticControl(double D0,double minD,double maxD,
		int IterDesired,double n);
	virtual ~StaticControl();

	void formResidual(double factor);

	// Form tangent and residual element by element 
	virtual void formElementalTangent(ModelElement* pModelElement);
	virtual void formElementalResidual(ModelElement* pModelElement,double time=0.);
	
	// Form residual node by node
	void formNodalResidual(ModelNode* pModelNode);

	// Methods that are used through analysis
	virtual void init();
	virtual void commit();
	virtual void rollback();
};

#endif