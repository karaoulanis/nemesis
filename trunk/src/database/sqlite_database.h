/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#ifndef SRC_DATABASE_SQLITE_DATABASE_H_
#define SRC_DATABASE_SQLITE_DATABASE_H_

#include <fstream>
#include <sqlite3.h>
#include "containers/containers.h"
#include "database/database.h"

class SQLiteDatabase: public Database {
  private:
  sqlite3* db;
  int executeQuery(const char* query);
  static int callback(void *NotUsed, int argc, char **argv, char **azColName);
  public:
  SQLiteDatabase();
  SQLiteDatabase(const char* workname);
  ~SQLiteDatabase();

  int closeDB();
  int createTable(const char* tableName);
  int deleteTable(const char* tableName);
  bool existsTable(const char* tableName);
  int useTable(const char* tableName);
  int storeData(const Packet& p);
  const Packet& retrieveData(int tag, int id);

  int beginTransaction();
  int commitTransaction();

  void exportToVtk(const char* tableName);

};

#endif  // SRC_DATABASE_SQLITE_DATABASE_H_
