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

#include <iostream>
#include <stdio.h>
#include "database/sqlite_database.h"
#include "exception/sexception.h"

using namespace std;

/**
 * Default constructor.
 */
SQLiteDatabase::SQLiteDatabase()
:Database()
{
	///@ todo Use synch/not synch to increase performance.
	///@ todo Use PRIMARY KEY index.
	///@ todo Check if db locked.
}
/**
 * Constructor.
 */
SQLiteDatabase::SQLiteDatabase(const char* workname)
{
	// Check if connected
	if(isConnected) this->closeDB();
	// Check if file exists and overwrite previous file
	char tmpName[512];
	strcpy(tmpName,workname);
	strcat(tmpName,".sdb");
	if(existsFile(tmpName)) 
		if(remove(tmpName)) 
			throw SException("[nemesis:%d] %s",9999,"Cannot overwrite database.");
	// Now it is ok to create the database
	if(sqlite3_open(tmpName,&db)) 
		throw SException("[nemesis:%d] %s",9999,"Cannot open database.");
	// Success
	strcpy(dbName,workname);
	isConnected=true;
}
/**
 * Destructor.
 */
SQLiteDatabase::~SQLiteDatabase()
{
}
/**
 * Callback function
 */
int SQLiteDatabase::
callback(void* /*NotUsed*/,int argc,char** argv,char** azColName)
{
	for(int i=0; i<argc; i++)
		printf("%s = %s\n",azColName[i],argv[i]?argv[i]:"NULL");
	printf("\n");
	return 0;
}
/**
 * Close database.
 */
int SQLiteDatabase::closeDB()
{
	if(isConnected)
	{
		int errCode=sqlite3_close(db);
		if(errCode) 
			throw SException("[nemesis:%d] %s",9999,"Error while closing database.");
		dbName[0]='\0';
		isConnected=false;
	}
	return 0;
}
/**
 * Create new table.
 */
int SQLiteDatabase::createTable(const char* tableName)
{
	if(!isConnected) 
		throw SException("[nemesis:%d] %s",9999,"Not connected to a database.");

	char query[2048];
	sprintf(query,
			"CREATE TABLE %s (tag INTEGER,id INTEGER,",tableName);

	char temp[64];
	for(int i=0;i<myPacket.dblSize;i++) 
	{
		if(i<10)		sprintf(temp,"dbl_0%d  FLOAT, ",i);
		else			sprintf(temp,"dbl_%d   FLOAT, ",i);
		strcat(query,temp);
	}
	for(int i=0;i<myPacket.intSize;i++) 
	{
		if(i<10)		sprintf(temp,"int_0%d  INTEGER, ",i);
		else			sprintf(temp,"int_%d   INTEGER, ",i);
		strcat(query,temp);
	}
	sprintf(temp,"chr CHAR(%d) );",myPacket.chrSize);
	strcat(query,temp);
	executeQuery(query);
	return 0;
}
/**
 * Delete given table.
 */
int SQLiteDatabase::deleteTable(const char* tableName)
{
	if(this->existsTable(tableName))
	{
		char query[64];
		sprintf(query,"DROP TABLE %s;",tableName);
		this->executeQuery(query);
	}
	return 0;
}
/**
 * Check if table exists.
 */
bool SQLiteDatabase::existsTable(const char* tableName)
{
	char query[128];
	char** result;
	char *errMsg;
	int nrow,ncol;
	bool ret=false;
	sprintf(query,"SELECT name FROM sqlite_master WHERE (type='table');");
	sqlite3_get_table(db,query,&result,&nrow,&ncol,&errMsg);
	for(int i=1;i<nrow+1;i++)
		if(!strcmp(result[i],tableName)) {ret=true; break;}
	sqlite3_free_table(result);
	return ret;
}
/**
 * Use database table.
 */
int SQLiteDatabase::useTable(const char* tableName)
{
	if(this->existsTable(tableName)) strcpy(tableInUse,tableName);	
	return 0;
}
/**
 * Begin transaction
 */
int SQLiteDatabase::beginTransaction()
{
	executeQuery("begin;");
	return 0;
}
/**
 * Commit transaction
 */
int SQLiteDatabase::commitTransaction()
{
	executeQuery("commit;");
	return 0;
}
/**
 * Store data.
 */
int SQLiteDatabase::storeData(const Packet& p)
{
	static char query[66000];
	sprintf(query,"INSERT INTO %s VALUES (%d,%d,",tableInUse,p.tag,p.id);  
	char *pos=query+strlen(query);
	int i;
	for(i=0;i<p.dblSize;i++) 
	{
		if(p.dblArray[i]!=p.dblDefault)	pos+=sprintf(pos, "%f, ",p.dblArray[i]);
		else							pos+=sprintf(pos, "NULL, ");
	}
	for(i=0;i<p.intSize;i++)
	{
		if(p.intArray[i]!=p.intDefault)	pos+=sprintf(pos, "%d, ",p.intArray[i]);
		else							pos+=sprintf(pos, "NULL, ");
	}
	pos+=sprintf(pos, "\"%s \"",p.chrArray);
	strcpy(pos,");");
	executeQuery(query);
	return 0;
}
/**
 * Retrieve data.
 */
const Packet& SQLiteDatabase::retrieveData(int tag,int id)
{
	char query[128];
	sprintf(query,"SELECT * FROM %s WHERE tag=%d AND id=%d",tableInUse,tag,id);
	char** result;
	char *errMsg;
	int i,nrow,ncol;
	sqlite3_get_table(db,query,&result,&nrow,&ncol,&errMsg);
	if(nrow!=1) exit(-8796);
	int count=ncol;
	myPacket.tag =atoi(result[count++]);
	myPacket.id  =atoi(result[count++]);
	for(i=0;i<myPacket.dblSize;i++)
		if(result[count++]!=0) myPacket.dblArray[i]=atof(result[count]);
	for(i=0;i<myPacket.intSize;i++) 
		if(result[count++]!=0) myPacket.intArray[i]=atoi(result[count]);
	if(result[count++]!=0) strcpy(myPacket.chrArray,result[count]);
	sqlite3_free_table(result);
	return myPacket;
}
/**
 * Execute sql statement.
 */
int SQLiteDatabase::executeQuery(const char* query)
{
	char *errMsg;
	int rc = sqlite3_exec(db,query,callback,0,&errMsg);
	if(rc!=SQLITE_OK)
	{
		cout<<"SQL error : "<<errMsg<<endl;
		cout<<"in query  : "<<query <<endl;
	}
	return rc;
}

void SQLiteDatabase::exportToVtk(const char* tableName)
{
	if(!(this->existsTable(tableName))) return;

	char vtkFileName[512];
	sprintf(vtkFileName,"%s.vtk",tableName);  
	ofstream vtkFile(vtkFileName);
	vtkFile.setf(ios_base::fixed,ios_base::floatfield);
	vtkFile.precision(6);

	// Print preamble
	vtkFile<<"# vtk DataFile Version 3.0"<<endl;
	vtkFile<<"nemesis vtk output"<<endl;
	vtkFile<<"ASCII"<<endl;
	vtkFile<<"DATASET UNSTRUCTURED_GRID"<<endl;
	vtkFile<<endl;

	char query[128];
	char *errMsg;
	int nrow,ncol;

	// Print nodes
	char** resNodes;
	sprintf(query,"SELECT * FROM %s WHERE tag=%d",tableName,1000);
	sqlite3_get_table(db,query,&resNodes,&nrow,&ncol,&errMsg);
	int nNodes=nrow;
	IDContainer nodalIds(nrow);
	vtkFile<<"POINTS "<<nrow<<" float"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		nodalIds[i-1]=atoi(resNodes[i*ncol+1]);
		vtkFile	<<atof(resNodes[i*ncol+2])<<' '
				<<atof(resNodes[i*ncol+3])<<' '
				<<atof(resNodes[i*ncol+4])<<endl;
	}
	vtkFile<<endl;

	//*************************************************************************
	// Print Elements
	//*************************************************************************
	char** resElems;
	sprintf(query,"SELECT * FROM %s WHERE tag BETWEEN 2000 AND 2999;",tableName);
	sqlite3_get_table(db,query,&resElems,&nrow,&ncol,&errMsg);
	int nValues=0;
	int nElems=nrow;
	IDContainer elemType(nElems);
	for(int i=1;i<nElems+1;i++) 
		nValues+=atoi(resElems[i*ncol+2+myPacket.dblSize])+1;
	vtkFile<<"CELLS "<<nElems<<" "<<nValues<<endl;
	for(int i=1;i<nElems+1;i++)
	{
		elemType[i-1]=atoi(resElems[i*ncol]);
		int nVertices=atoi(resElems[i*ncol+2+myPacket.dblSize]);
		vtkFile<<nVertices<<' ';
		for(int j=0;j<nVertices;j++)
			vtkFile<<Containers::index_find(nodalIds,atoi(resElems[i*ncol+2+myPacket.dblSize+1+j]))<<' ';
		vtkFile<<endl;
	}
	vtkFile<<endl;

	vtkFile<<"CELL_TYPES "<<elemType.size()<<endl;
	for(unsigned ii=0;ii<elemType.size();ii++)
	{
		switch(elemType[ii])
		{
		case 2001:
		case 2002:
		case 2003:
		case 2004:
		case 2005:
		case 2006:
			vtkFile<< 3<<endl; break;
		case 2007:
			vtkFile<< 9<<endl; break;
		case 2008:
			vtkFile<< 5<<endl; break;
		case 2009:
			vtkFile<<22<<endl; break;
		case 2010:
			vtkFile<<12<<endl; break;
		case 2011:
			vtkFile<<10<<endl; break;
		case 2012:
			vtkFile<<3<<endl; break;
		case 2013:
			vtkFile<<21<<endl; break;
		default:
			vtkFile<<0<<endl;
		}
	}
	vtkFile<<endl;

	//*************************************************************************
	// Cell data
	//*************************************************************************
	// Materials 
	vtkFile<<"CELL_DATA     "<<nElems<<endl;
	vtkFile<<"SCALARS materials int"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nElems+1;i++)
			vtkFile<<atoi(resElems[i*ncol+2+myPacket.dblSize+29])<<endl;
	vtkFile<<endl;

	// Plastic points 
	vtkFile<<"SCALARS plastic_points int"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nElems+1;i++)
			vtkFile<<atoi(resElems[i*ncol+2+myPacket.dblSize+30])<<endl;
	vtkFile<<endl;

	//*************************************************************************
	// Point data
	//*************************************************************************

	// Print ux
	vtkFile<<"POINT_DATA     "<<nNodes<<endl;
	vtkFile<<endl;

	vtkFile<<"VECTORS ux float"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3]==0) vtkFile<<0.<<" ";
		else vtkFile<<atof(resNodes[i*ncol+2+3])<<" ";
		vtkFile<<0.<<" ";
		vtkFile<<0.<<" ";
		vtkFile<<endl;
	}
	vtkFile<<endl;

	// Print uy
	vtkFile<<"VECTORS uy float"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		vtkFile<<0.<<" ";
		if(resNodes[i*ncol+2+4]==0) vtkFile<<0.<<" ";
		else vtkFile<<atof(resNodes[i*ncol+2+4])<<" ";
		vtkFile<<0.<<" ";
		vtkFile<<endl;
	}
	vtkFile<<endl;

	// Print uz
	vtkFile<<"VECTORS uz float"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		vtkFile<<0.<<" ";
		vtkFile<<0.<<" ";
		if(resNodes[i*ncol+2+5]==0) vtkFile<<0.<<" ";
		else vtkFile<<atof(resNodes[i*ncol+2+5])<<" ";
		vtkFile<<endl;
	}
	vtkFile<<endl;

	// Print (ux,uy,uz)
	vtkFile<<"VECTORS u float"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3]==0) vtkFile<<0.<<' ';
		else vtkFile<<atof(resNodes[i*ncol+2+3])<<' ';
		if(resNodes[i*ncol+2+4]==0) vtkFile<<0.<<' ';
		else vtkFile<<atof(resNodes[i*ncol+2+4])<<' ';
		if(resNodes[i*ncol+2+5]==0) vtkFile<<0.<<' ';
		else vtkFile<<atof(resNodes[i*ncol+2+5])<<' ';
		vtkFile<<endl;
	}
	vtkFile<<endl;

	// Print ux
	vtkFile<<"SCALARS uxp float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+3])<<endl;
	}
	vtkFile<<endl;

	// Print uy
	vtkFile<<"SCALARS uyp float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+4]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+4])<<endl;
	}
	vtkFile<<endl;

	// Print uz
	vtkFile<<"SCALARS uzp float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+5]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+5])<<endl;
	}
	vtkFile<<endl;

	// Print sx
	vtkFile<<"SCALARS sx float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3+3*16]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+3+3*16])<<endl;
	}
	vtkFile<<endl;

	// Print sy
	vtkFile<<"SCALARS sy float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3+3*16+1]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+3+3*16+1])<<endl;
	}
	vtkFile<<endl;

	// Print sz
	vtkFile<<"SCALARS sz float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3+3*16+2]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+3+3*16+2])<<endl;
	}
	vtkFile<<endl;

	// Print sxy
	vtkFile<<"SCALARS sxy float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3+3*16+3]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+3+3*16+3])<<endl;
	}
	vtkFile<<endl;

	// Print syz
	vtkFile<<"SCALARS syz float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3+3*16+4]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+3+3*16+4])<<endl;
	}
	vtkFile<<endl;

	// Print szx
	vtkFile<<"SCALARS szy float"<<endl;
	vtkFile<<"LOOKUP_TABLE DEFAULT"<<endl;
	for(int i=1;i<nNodes+1;i++) 
	{
		if(resNodes[i*ncol+2+3+3*16+5]==0) vtkFile<<0.<<endl;
		else vtkFile<<atof(resNodes[i*ncol+2+3+3*16+5])<<endl;
	}
	vtkFile<<endl;
	vtkFile.close();
	
	sqlite3_free_table(resNodes);
	sqlite3_free_table(resElems);
}
