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

#include <NodeTracker.h>

NodeTracker::NodeTracker()
{
}
NodeTracker::NodeTracker(int ID,int nodeID)
{
	myNode=pD->get<Node>(pD->getNodes(),nodeID);
}
NodeTracker::~NodeTracker()
{
}
Node* NodeTracker::getNode()
{
	return myNode;
}
void NodeTracker::keepTrack(double lambda_,double time_)
{
	lambda.push_back(lambda_);
	time.push_back(time_);
	data.push_back(myNode->getPacket());
}
void NodeTracker::save(std::ostream& s)
{
/*	static Vector crds(3);
	crds[0]=x1;
	crds[1]=x2;
	crds[2]=x3;
	s<<"NODE "	<<' ';
	s<<"tag	"	<<1000<<' '<<myTag<<' ';
	s<<"id "	<<1000<<' '<<myID<<' ';
	s<<"crds "	<<' '<<crds;
	s<<"disp "	<<' '<<dispConvg;
	s<<"velc "	<<' '<<velcConvg;
	s<<"accl "	<<' '<<acclConvg;
	s<<"stress "<<' '<<stress;
	s<<"strain "<<' '<<strain;
	s<<"dsens  "<<' '<<dispSensi;
	s<<"eigen  "<<' '<<eigenVecs;
	s<<"END "<<' ';*/
}