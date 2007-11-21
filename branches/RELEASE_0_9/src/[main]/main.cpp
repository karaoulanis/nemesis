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

#include <NemesisConfig.h>
#include <PyParser.h>

int main(int argc,char* argv[])
{
	/*****************************************************************************
	* Info
	*****************************************************************************/
	cout<<"nemesis "<<_NEMESIS_VERSION
	    <<" [www.nemesis-project.org] (C)2004-2007 F.E. Karaoulanis."<<endl;
	cout<<"Built on "	<<_NEMESIS_PLATFORM
		<<" with "		<<_NEMESIS_COMPILER_NAME
		<<" ver."		<<_NEMESIS_COMPILER_VERSION
		<<" on "		<<__DATE__
		<<" at "		<<__TIME__<<"."<<endl;
	cout<<"Licensed under GPL v2.0; provided \"AS IS\"; comes with ABSOLUTELY NO WARRANTY."<<endl;
	/*****************************************************************************
	* Python parser interactive
	*****************************************************************************/
	if(argc==1)
	{
		cout<<"Type help(nemesis) for additional info."<<endl;
		cout<<"Use Ctrl-Z plus Return to exit."<<endl;
		Parser* theParser=new PyParser();
		theParser->parse();
		delete theParser;
	}
	/*****************************************************************************
	// Python parser with file
	*****************************************************************************/
	else if(argc==3 && !strcmp(argv[1],"-p"))
	{
		Parser* theParser=new PyParser();
		theParser->parse(argv[2]);
		delete theParser;
	}
	/*****************************************************************************
	// Wrong arguments
	*****************************************************************************/
	else
	{
		cout<<"\nUsage: nemesis"<<endl
			<<"         nemesis -p file"<<endl<<endl;
	}
	// Exit
	return 0;
}