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

#ifndef _DATABASE_H
#define _DATABASE_H

#include <Packet.h>

class Database
{
protected:
	char dbName[512];
	char tableInUse[512];
	bool isConnected;
	static Packet myPacket;
public:
	Database();
	virtual ~Database();
	virtual int closeDB()=0;

	virtual int createTable(const char* tableName)=0;
	virtual int deleteTable(const char* tableName)=0;
	virtual int useTable(const char* tableName)=0;
	virtual int storeData(const Packet& p)=0;
	virtual bool existsTable(const char* tableName)=0;
	virtual const Packet& retrieveData(int tag,int id)=0;

	virtual void exportToVtk(const char* tableName)=0;

	virtual int beginTransaction();
	virtual int commitTransaction();

	virtual bool existsFile(const char* filename);
};

#endif
