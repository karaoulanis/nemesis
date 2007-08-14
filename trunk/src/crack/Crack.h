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

#ifndef _CRACK_H
#define _CRACK_H

#include <DomainObject.h>
#include <vector>
#include <Containers.h>

using namespace std;

struct CrackPoint
{
	double x,y,z;
	CrackPoint(double x_,double y_,double z_=0.)
	{
		x=x_;
		y=y_;
		z=z_;
	}
};
/**
 * The Crack class.
 */
class Crack: public DomainObject
{
protected:
	vector<CrackPoint*> myCrackPoints;
public:
	Crack();
	Crack(int ID,double xS_,double yS_,double xT_,double yT_);
	~Crack();
	const vector<CrackPoint*>& getCrack();
	CrackPoint* getCrackPoint(int n);
};
#endif
