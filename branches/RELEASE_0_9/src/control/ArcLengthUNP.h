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
//            C.G. Panagiotopoulos (pchr@civil.auth.gr)
//            E.A. Paraskevopoulos (eapcivil@yahoo.gr)
//*****************************************************************************

#ifndef _ARCLENGTHUNP_H
#define _ARCLENGTHUNP_H

#include <StaticControl.h>

/**
 * The ArcLengthUNP Class.
 * ArcLengthUNP (Updated Normal Plane)is a static control that is based on 
 * Crisfield's book "Non-linear Finite Element Analysis of Solids and 
 * Strucures", Vol.1.
 * For more details see:
 * \li Theory p.274	\n
 * \li Implementation p.276-278 (Similar to AecLengthSpherical)\n
 * \li Predictor p.285-286	\n
 */
class ArcLengthUNP :public StaticControl
{
private:
public:
	ArcLengthUNP(double DL0,double minDL,double maxDL,
								int IterDesired,double n);
	~ArcLengthUNP();

	// Methods for incremental/iterative algorithms
	void predict();
	void correct();
};

#endif
