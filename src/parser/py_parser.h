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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#ifndef NEMESIS_PARSER_PY_PARSER_H_
#define NEMESIS_PARSER_PY_PARSER_H_

#ifdef _DEBUG
  #include < Python.h>
  #undef Py_DEBUG
#else
  #include < Python.h>
#endif

#include "parser/parser.h"

class PyParser: public Parser {
  private:
  int initModules();
  public:
  // Constructors and destructor
  PyParser();
  virtual ~PyParser();
  int parse();
  int parse(char* filename);
};
#endif  // NEMESIS_PARSER_PY_PARSER_H_
