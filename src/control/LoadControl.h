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

#ifndef _LOADCONTROL_H
#define _LOADCONTROL_H

#include <StaticControl.h>

/**
 * The LoadControl Class.
 * LoadControl is a static control that increases \f$\lambda_n\f$ according 
 * to the relation \f$\lambda_n=\lambda_o+\Delta\lambda\f$, where 'o' stands 
 * for old (previous) step and 'n' stands for new (current step). 
 * \f$\Delta\lambda\f$ may vary between different steps (auto-incrementation)
 * but is constant within each step, i.e. \f$\delta\lambda\f$=0. 
 * \f$\lambda_n\f$ is in fact a trial \f$\lambda_n\f$, and becomes 
 * \f$\lambda_o\f$ for the next step if and only if this step converges.\n
 */
class LoadControl :public StaticControl
{
private:
public:
	// Constructor and destructor
	LoadControl(double DL0,double minDL,double maxDL,int IterDesired,double n,double DeltaTime);
	~LoadControl();

	// Methods for incremental/iterative algorithms
	void predict();
	void correct();
};

#endif
