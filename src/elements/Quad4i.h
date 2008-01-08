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

#ifndef _QUAD4I_H
#define _QUAD4I_H

#include <Quad4.h>

class Quad4i: public Quad4
{
private:
	static double shpStd[4][3][4];
	static double shpInc[2][3][4];

	static double detJ[4];
	static std::vector<int> perm;
	
	Vector aTrial;
	Vector aConvg;
	void shapeFunctions();
	void getBStd(Matrix& B,int node,int gPoint);
	void getBInc(Matrix& B,int node,int gPoint);
	void getKdd(Matrix& K);
	void getKda(Matrix& K);
	void getKaa(Matrix& K);
public:
	// Constructors and Destructor
	Quad4i();
	Quad4i(int ID,int Node_1,int Node_2,int Node_3,int Node_4,int MatID);
	~Quad4i();
	
	const Matrix& getK();
    const Matrix& getM();
	const Vector& getR();
	void update();
	void commit();
};

#endif
