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

#include "main/nemesis_config.h"
#include "parser/py_parser.h"

int main(int argc, char* argv[]) {
  /*****************************************************************************
  * Info
  *****************************************************************************/
  cout  << "nemesis-"   << NEMESIS_VERSION
        << "[nemesis-project.org] "
        << "Copyright (C) 2004-2010 F.E. Karaoulanis." << endl;
  cout  << "Built on "  << NEMESIS_PLATFORM
        << " with "     << NEMESIS_COMPILER_NAME
        << " ver."      << NEMESIS_COMPILER_VERSION
        << " on "       << __DATE__
        << " at "       << __TIME__ << " ." <<endl;
  cout  << "Licensed under GPL v3.0; "
        << "provided \"AS IS\"; "
        << "comes with ABSOLUTELY NO WARRANTY." <<endl;

  if (argc == 1) {
  /*****************************************************************************
  * Python parser interactive
  *****************************************************************************/
    cout << "Type help(nemesis) for additional info." << endl;
    cout << "Use Ctrl-Z plus Return to exit." << endl;
    Parser* theParser = new PyParser();
    theParser->parse();
    delete theParser;
  } else if (argc == 3 && !strcmp(argv[1], "-p")) {
  /*****************************************************************************
  // Python parser with file
  *****************************************************************************/
    Parser* theParser = new PyParser();
    theParser->parse(argv[2]);
    delete theParser;
  } else {
  /*****************************************************************************
  // Wrong arguments
  *****************************************************************************/
    cout  << "\nUsage: nemesis" <<endl
          << "nemesis -p file"  <<endl << endl;
  }
  // Exit
  return 0;
}
