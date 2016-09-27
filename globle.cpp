#include <iostream>
#include <stdlib.h>	// standard C library function
#include "Space.h"
#include "globle.h"
#include <vector>


#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include <cmath>
#include <iostream>


using namespace std;

// this file is only for globle variables

bool bUniformDiel;
bool bVerbose;
float fRMid;
float fScale;
float fRadiusMax;
float fIonStrenth;
float fExternRadius;
float fRadPrb[2];

int iTestGloble;
int iBoundNum;
int ibmx;
int iGrid;
int iNObject;
int iGaussian;
int iNatom;
int iNMedia;

bool bDebug;
bool bDebMap[65][65][65];
bool bOnlyMol;

//float * fCoord;
//float * fRadius;
//float * fCharge;



SExtrema <float> * sLimGridUnit;
SExtrema <float> * sLimObject;

delphipdb_struc * sDelPhiPDB;

SGrid <int> iEpsMap[65][65][65];

SGrid <float> cOldMid;
SGrid <float> cMin;
SGrid <float> cMax;

SGrid <float> * xn1;
SGrid <float> * xn2;

int * iAtomMed;

//semi_global:

int iBNumSurf;
float radpmax;



int * ast;

float * r0;
float * r02;
float * rs2;



SExtrema <int> LimEps;




