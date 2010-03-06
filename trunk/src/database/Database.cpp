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

#include <Database.h>
#include <cstring>
#include <fstream>
using namespace std;

Packet Database::myPacket;

Database::Database()
	:isConnected(false)
{
}
Database::~Database()
{
}
bool Database::existsFile(const char* filename)
{
	bool ret=true;
	fstream tmpFile;
	tmpFile.open(filename,ios::binary|ios::in);
	if(tmpFile.fail())
	{
		tmpFile.clear();
		ret=false;
	}
	tmpFile.close();
	return ret;
}
int Database::beginTransaction()
{
	return 0;
}
int Database::commitTransaction()
{
	return 0;
}
