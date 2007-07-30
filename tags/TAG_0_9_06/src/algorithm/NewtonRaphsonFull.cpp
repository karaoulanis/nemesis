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

#include <NewtonRaphsonFull.h>

NewtonRaphsonFull::NewtonRaphsonFull()
{
	myTag=TAG_ALGORITHM_NEWTON_RAPHSON_FULL;
}
NewtonRaphsonFull::~NewtonRaphsonFull()
{
}
int NewtonRaphsonFull::solveStep(int n)
{
	// Predictor phase
	pA->getControl()->formTangent();
	pA->getControl()->predict();
	pA->getConvergenceNorm()->newStep();
	pA->getControl()->formResidual(pA->getControl()->getLambda());

	// Corrector phase
	int check;
	while((check=pA->getConvergenceNorm()->update())>0)
	{
		pA->getControl()->formTangent();
		pA->getSOE()->solve();
		pA->getControl()->correct();
 		pA->getControl()->formResidual(pA->getControl()->getLambda());
	}
	return check;
}
