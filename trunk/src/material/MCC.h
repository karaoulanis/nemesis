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

#ifndef _MCC_H
#define _MCC_H

#include <Surface.h>
#include <Vector.h>

class MCC: public Surface
{
protected:
	double M;
	double po;
	double kappa;
	double lambda;
	void find_C(const Vector& s,const Vector& e);
public:
	MCC();
	MCC(double M_,double po_,double kappa_,double lambda_);
	~MCC();
	double get_f(const Vector& s,const Vector& e);
	const Vector& get_dfdq(const Vector& s,const Vector& e);
	const Matrix& get_df2dsq(const Vector& s,const Vector& e);
	const Matrix& get_df2dqq(const Vector& s,const Vector& e);
};
#endif
