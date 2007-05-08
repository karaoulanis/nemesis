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

#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

#include <Domain.h>
#include <vector>

struct cDof
{
	Node* pNode;
	int dof;
	double coeff;
};
class Constraint: public DomainObject
{
protected:
	static int nConstraints;
	std::vector<cDof> theCDofs;
	double cVal;
	double fTrial;
	double fConvg;
public:
	// Constructors
	Constraint();

	void setcDof(int NodeID,int dof,double coeff);
	void setcVal(double val);

	int getncDofs();
	const cDof& getcDof(int i);
	virtual double getcVal(double time=0.);

	double getDisp(int i);
	double getVelc(int i);
	double getDispConvg(int i);
	double getDispTrial(int i);
	double getVelcConvg(int i);
	double getAcclConvg(int i);
	double getF();

	void update(double f);
	void commit();
};
#endif
