/*
	File:			Basics.cc

	Function:		Implements Basics.h

	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, Andrew Willmott

	Notes:			

*/

#include "Basics.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;


// --- Error functions for range and routine checking -------------------------


ostream &operator << (ostream &s, Bool &b)
{
	if (b)
		s << "true";
	else
		s << "false";

	return(s);
}

Void _Assert(Int condition, char *errorMessage, Char *file, Int line)
{
	if (!condition)
	{
		char reply;
		
		cerr << "\n*** Assert failed (line " << line << " in " << 
			file << "): " << errorMessage << endl;
		cerr << "    Continue? [y/n] ";
		cin >> reply;
		
		if (reply != 'y')
		{
			*((long *) 0) = 0; // Force a core dump/debugger break
			exit(1);
		}
	}
}

Void _Expect(Int condition, char *warningMessage, Char *file, Int line)
{
	if (!condition)
		cerr << "\n*** Warning (line " << line << " in " << file << "): " <<
			warningMessage << endl;
}

Void _CheckRange(Int i, Int lowerBound, Int upperBound, 
				 Char *rangeMessage, Char *file, Int line)
{
	if (i < lowerBound || i >= upperBound)
	{
		char reply;
		
		cerr << "\n*** Range Error (line " << line << " in " << file <<
			"): " << rangeMessage << endl;	
		cerr << "    Continue? [y/n] ";
		cin >> reply;
		
		if (reply != 'y')
		{
			*((long *) 0) = 0;
			exit(1);
		}
	}
}

