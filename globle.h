#ifndef GLOBAL_H // header guards
#define GLOBAL_H


#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include <cmath>
#include <iostream>



// extern tells the compiler this variable is declared elsewhere
extern bool bUniformDiel,bVerbose,bOnlyMol,bDebug;

extern int iTestGloble,iGrid,ibmx,iNObject,iBoundNum,iGaussian,iNatom,iNMedia;

extern float fRMid,fScale,fRadiusMax,fIonStrenth,fExternRadius,fRadPrb[2];


//extern struct coord
/*
struct coord{
    float x,y,z;
    };
struct int_coord{
    int i,j,k;
    };

struct object_min_max{
    coord cMin, cMax;
};

struct int_coord_min_max{
    int_coord iMin,iMax;
};
*/

struct delphipdb_struc{
    SGrid <float> xyz;
    float charge;
    float radius;
};



extern SGrid <float> cOldMid,cMin,cMax;


extern SExtrema <float> * sLimGridUnit;
extern SExtrema <float> * sLimObject;
extern SGrid <int> iEpsMap[65][65][65];
extern SGrid <float> * xn1;
extern SGrid <float> * xn2;
extern int * iAtomMed;



extern bool bDebMap[65][65][65];
extern delphipdb_struc * sDelPhiPDB;


// semi global

extern int iBNumSurf;

extern float radpmax;

extern SExtrema <int> LimEps;

extern int * ast;

extern float * r0;
extern float * r02;
extern float * rs2;



#endif

