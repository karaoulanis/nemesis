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

#ifndef _SOE_H
#define _SOE_H

// C++ system files
#include <cmath>
#include <iostream>
#include <valarray>

// Project files (alphabetically)
#include "analysis/analysis.h"
#include "analysis/analysis_object.h"
#include "numeric/Lapack.h"
#include "numeric/Matrix.h"
#include "numeric/Vector.h"
#include "containers/Containers.h"

// Forward declarations
class Analysis;
class Model;

class SOE: public AnalysisObject
{
protected:
	int theSize;
	valarray<double> A;
	Vector X;
	Vector B;
	bool isLUFactored;
	valarray<int> IPIV;
public:
	SOE();
	virtual ~SOE();

	virtual int insertMatrixIntoA(const Matrix& Ke,const IDContainer& EFTable,
								  double factor=1.0)=0;
	virtual int insertVectorIntoB(const Vector& Ve,const IDContainer& EFTable,
								  double factor=1.0);
	virtual int insertMatrixIntoA(const Matrix& Be,const IDContainer& EFTable,
								  const IDContainer& SFTable,double factor=1.0);

	virtual void zeroA();
	virtual void zeroB();	
	virtual void zeroX();
	virtual void zero();

	const Vector& getX();
	const Vector& getB();
	void addB(const Vector& v);
	void setB(const Vector& v);
	void setX(const Vector& v);

	virtual void setTheSize()=0;
	virtual int getEigenSign()=0;
	virtual void print()=0;
	virtual void printSolution();

	virtual int solve()=0;
	int plotGraph(const char* s);
};
#endif
