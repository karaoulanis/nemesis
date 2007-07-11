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

#include <Tracker.h>
#include <iostream>

Tracker::Tracker()
{
}
Tracker::~Tracker()
{
}
const int Tracker::getSteps()
{
	return myRecords.size();
}
void Tracker::track(double lambda_,double time_,string data_)
{
	TrackerRecord record;
	record.lambda=lambda_;
	record.time=time_;
	record.data=data_;
	myRecords.push_back(record);
}
void Tracker::save(std::ostream& s)
{
	s<<"TRACKER "	<<' ';
	s<<"steps "	<<1000<<' '<<myRecords.size()<<' ';
	for(unsigned i=0;i<myRecords.size();i++)
	{
		s<<"lambda "	<<1010<<' '<<myRecords[i].lambda	<<' ';
		s<<"time "		<<1010<<' '<<myRecords[i].time		<<' ';
		s<<"data "		<<1020<<' '<<myRecords[i].data		<<' ';
	}
	s<<"END "<<' ';
}
