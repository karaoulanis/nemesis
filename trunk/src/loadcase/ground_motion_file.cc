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

#include "loadcase/GroundMotionFile.h"

GroundMotionFile::GroundMotionFile()
{
}
GroundMotionFile::GroundMotionFile(int dof_,std::istream& s,double dt_,double scale_)
	:Load()
{
	dof=dof_-1;
	dt=dt_;
	scale=scale_;
	while(!s.eof())
	{
		double d;
		s>>d;
		data.push_back(d);
	}
}
void GroundMotionFile::apply(double fact,double time)
{
	unsigned n=(int)floor(time/dt);
	if(n<0||n>=data.size()) throw SException("[nemesis:%d] %s",9999,"Time exceeds given values.");
	double t1=n*dt;
	double t2=(n+1)*dt;
	double d1=data[n];
	double d2=data[n+1];
    double d=d1+(time-t1)*(d2-d1)/(t2-t1);
	const std::map<int,Element*>& c=pD->getElements();
	std::map<int,Element*>::const_iterator i;
	for(i=c.begin();i!=c.end();i++) i->second->addGroundMotion(dof,fact*scale*d);
}
