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

#ifndef _QUAD4DISPPLAIN_H
#define _QUAD4DISPPLAIN_H

#include <Quad4.h>

class Quad4DispPlain: public Quad4
{
public:
	// Constructors and Destructor
	Quad4DispPlain();
	Quad4DispPlain(int ID,int Node_1,int Node_2,int Node_3,int Node_4,int MatID,
				 int integrationRuleXi,int integrationRuleEta);
	~Quad4DispPlain();
	
	const Matrix& getK();
    const Matrix& getM();
	const Vector& getR();

	void update();
};

#endif