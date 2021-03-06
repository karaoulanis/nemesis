/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2011 F.E.Karaoulanis [http://www.nemesis-project.org]     *
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

#include <stdio.h>
#include "main/nemesis_config.h"
#include "parser/py_parser.h"

int main(int argc, char* argv[]) {
  /*****************************************************************************
  * Info
  *****************************************************************************/
  printf("nemesis %s. ", NEMESIS_VERSION);
  printf("Copyright (C) 2004-2011 Fotios E. Karaoulanis.\n");
  printf("Built on %s ", NEMESIS_PLATFORM);
  printf("with %s ", NEMESIS_COMPILER_NAME);
  printf("ver.%s ", NEMESIS_COMPILER_VERSION);
  printf("on %s at %s.\n", __DATE__, __TIME__);
  printf("Licensed under GPL v3.0; ");
  printf("provided \"AS IS\"; ");
  printf("comes with ABSOLUTELY NO WARRANTY.\n");
  printf("Please visit http://www.nemesis-project.org for more information.\n");

  if (argc == 1) {
  /*****************************************************************************
  * Python parser interactive
  *****************************************************************************/
    printf("Type help(nemesis) for additional info.\n");
    printf("Use Ctrl-Z plus Return to exit.\n");
    Parser* theParser = new PyParser();
    theParser->Parse();
    delete theParser;
  } else if (argc == 3 && !strcmp(argv[1], "-p")) {
  /*****************************************************************************
  // Python parser with file
  *****************************************************************************/
    Parser* theParser = new PyParser();
    theParser->Parse(argv[2]);
    delete theParser;
  } else {
  /*****************************************************************************
  // Wrong arguments
  *****************************************************************************/
    printf("\nUsage: nemesis\nnemesis -p file\n\n");
  }
  // Exit
  return 0;
}
