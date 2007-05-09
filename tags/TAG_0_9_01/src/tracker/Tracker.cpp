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

Tracker::Tracker()
{
}
Tracker::Tracker(int ID,int nodeID)
:DomainObject(ID)
{
	///@todo: check for errors
	myNode=pD->get<Node>(pD->getNodes(),nodeID);
}
Tracker::~Tracker()
{
}
Node* Tracker::getNode()
{
	return myNode;
}
const int Tracker::getSteps()
{
	return lambda.size();
}
const double Tracker::getLambda(int step)
{
	return lambda.at(step);
}
const double Tracker::getTime(int step)
{
	return time.at(step);
}
const Packet& Tracker::getPacket(int step)
{
	return data.at(step);
}
void Tracker::keepTrack(double lambda_,double time_)
{
	lambda.push_back(lambda_);
	time.push_back(time_);
	data.push_back(myNode->getPacket());
}
