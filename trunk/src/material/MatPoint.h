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

#ifndef _MATPOINT_H
#define _MATPOINT_H

#include <DomainObject.h>
#include <MultiaxialMaterial.h>
#include <InitialStresses.h>

class MatPoint: public DomainObject
{
private:
	double x1,x2,x3;
	double xi1,xi2,xi3;
	double weight;
	static const double GaussCoords[7][7];
	static const double GaussWeights[7][7];
	MultiaxialMaterial* myMaterial;
	static int IDCounter;
public:
	MatPoint();
	MatPoint(MultiaxialMaterial* mat,int index,int p1);
	MatPoint(MultiaxialMaterial* mat,int index1,int index2,int p1,int p2);
	MatPoint(MultiaxialMaterial* mat,int index1,int index2,int index3,int p1,int p2,int p3);
	~MatPoint();

	void setX(double x1_,double x2_=0,double x3_=0);
	inline MultiaxialMaterial* getMaterial()	{return myMaterial;}
	bool isPlastic()							{return myMaterial->isPlastic();}
	void setInitialStresses(InitialStresses* pInitialStresses);

	inline double getx1()			{return x1;			}
	inline double getx2()			{return x2;			}
	inline double getx3()			{return x3;			}
	inline double getxi1()			{return xi1;		}
	inline double getxi2()			{return xi2;		}
	inline double getxi3()			{return xi3;		}
	inline double getWeight()		{return weight;		}
	const Packet& getPacket()		{return thePacket;	}
	void setPacket(const Packet& p) {/*does nothing*/	}
};

#endif
