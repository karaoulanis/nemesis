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

#ifndef _LADEDUNCANSURFACE_H
#define _LADEDUNCANSURFACE_H

#include <Surface.h>

class LadeDuncanSurface: public Surface
{
private:
	double K;
	void find_C(const Vector& s,const double a);
public:
	LadeDuncanSurface();
	LadeDuncanSurface(double K_);
	~LadeDuncanSurface();
	
	double get_f(const Vector& s,const double q);
	const double  get_dfda(const Vector& s,const double a);
	const Vector& get_df2dsa(const Vector& s,const double a);
	const double  get_df2daa(const Vector& s,const double a);
};
#endif
