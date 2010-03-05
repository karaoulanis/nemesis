/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
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

#ifndef _HOEKBROWN_H
#define _HOEKBROWN_H

#include <MultiaxialElastoPlastic.h>

/**
 * The Hoek-Brown Class.
 */
class HoekBrown: public MultiaxialMaterial
{
private:
	MultiaxialMaterial* myElastic;
	static Matrix C;
	static Matrix C3;
	bool plastic;
	int inaccurate;

	Vector f;
	vector<Vector> dfds;
	vector<Vector> dgds;
	vector<Matrix> d2gdsds;

	double aTrial,aConvg;
public:
	HoekBrown();
	HoekBrown(int ID,int elasticID,double si,double sp,double mb,double mbb,double alpha);
	~HoekBrown();

	MultiaxialMaterial* getClone();
	void setStrain(const Vector& De);
	void commit();
	const Matrix& getC();
	bool isPlastic();

	void find_f(const Vector& s,double q);
	void find_dfds(const Vector& s,double q);
	void find_dgds(const Vector& s,double q);
	void find_d2gdsds(const Vector& s,double q);

	// Tracker member functions
	void track();
};
#endif