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

#include <MCC.h>

MCC::MCC()
{
}
MCC::MCC(double M_,double pc_)
{
	M=M_;
	pc=pc_;
}
MCC::~MCC()
{
}
double MCC::get_f(const Vector& s)
{
	double q=s.q();
	double p=s.p();
	return q*q+M*M*p*(p-pc);
}
const Vector& MCC::get_dfds(const Vector& s)
{
	static Vector ret(6,0.);
	return ret;
}
const Vector& MCC::get_dfdq(const Vector& s)
{
	static Vector ret(6,0.);
	return ret;
}
const Matrix& MCC::get_df2dss(const Vector& s)
{
	static Matrix ret(6,6,0.);
	return ret;
}
const Matrix& MCC::get_df2dsq(const Vector& s)
{
	static Matrix ret(6,6,0.);
	return ret;
}
const Matrix& MCC::get_df2dqq(const Vector& s)
{
	static Matrix ret(6,6,0.);
	return ret;
}
