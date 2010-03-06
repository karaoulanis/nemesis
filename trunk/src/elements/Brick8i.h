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

#ifndef _BRICK8I_H
#define _BRICK8I_H

#include <Brick8.h>

class Brick8i: public Brick8
{
private:
	static double shpStd[8][4][8];
	static double shpInc[3][4][8];
	static double detJ[8];
	
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
	Brick8i();
	Brick8i(int ID,
				int Node_1,int Node_2,int Node_3,int Node_4,	
				int Node_5,int Node_6,int Node_7,int Node_8,
				int matID);
	~Brick8i();
	
	const Matrix& getK();
    const Matrix& getM();
	const Vector& getR();
	void update();
	void commit();
	void getB(Matrix& B,int node,int gPoint) {}

};

#endif
