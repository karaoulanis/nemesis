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

#include "crosssection/CrossSection.h"
#include "numeric/Numeric.h"

RectangularCrossSection::RectangularCrossSection(int ID,double width,double height)
:CrossSection(ID)
{
	w=width;
	h=height;
	A=w*h;
	As2=5./6.*A;
	As3=5./6.*A;
	J1=h*w*w*w*(num::d13-0.21*w/h*(1-w*w*w*w/(12*h*h*h*h)));
	J2=h*w*w*w/12.;
	J3=w*h*h*h/12.;
	myTag=TAG_CROSSSECTION_RECTANGLE;
}

UserDefinedCrossSection::UserDefinedCrossSection(int ID,
													double A_in,
													double As2_in,
													double As3_in,
													double J1_in,
													double J2_in,
													double J3_in,
													double h2_in,
													double h3_in)
:CrossSection(ID)
{
	A=A_in;
	As2=As2_in;
	As3=As3_in;
	J1=J1_in;
	J2=J2_in;
	J3=J3_in;
	h2=h2_in;
	h3=h3_in;
	myTag=TAG_CROSSSECTION_USER_DEFINED;
}
