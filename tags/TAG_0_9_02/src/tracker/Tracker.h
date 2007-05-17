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

#ifndef _TRACKER_H
#define _TRACKER_H

#include <Domain.h>
#include <DomainObject.h>
#include <vector>

// Forward Declerations
class Domain;
class Node;

class Tracker: public DomainObject
{
private:
	Node* myNode;
	int dof;
	std::vector<double> lambda;
	std::vector<double> time;
	std::vector<Packet> data;
public:
	Tracker();
	Tracker(int ID,int nodeID);
	virtual ~Tracker();
	Node* getNode();
	const int getSteps();
	const double getLambda(int step);
	const double getTime(int step);
	const Packet& getPacket(int step);

	virtual void keepTrack(double lambda_,double time_);
};

#endif
