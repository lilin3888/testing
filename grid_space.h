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


//int_coord(): covert 3 integer numbers to be a SGrid <integer> structure
SGrid <int> int_coord( int a,  int b,  int c)
{
    SGrid <int>  d;
    d.nX=a;
    d.nY=b;
    d.nZ=c;

    return(d);
}

//coord(): covert 3 float numbers to be a SGrid <float> structure
SGrid <float> coord( float a,  float b,  float c)
{
    SGrid <float>  d;
    d.nX=a;
    d.nY=b;
    d.nZ=c;

    return(d);
}

//Float2Int(): 1. convert SGrid <float> structure to SGrid <int> structure; 2. convert float to int
SGrid <int> Float2Int( SGrid <float> a )
{
    SGrid <int> b;

    b.nX=int(a.nX);
    b.nY=int(a.nY);
    b.nZ=int(a.nZ);

    return(b);
}

int Float2Int(float a){

    int b;
    b=int(a);
    return(b);

}

//Int2Float(): 1. convert SGrid <int> structure to SGrid <float> structure; 2. convert int to float
SGrid <float> Int2Float( SGrid <int> a )
{
    SGrid <float> b;

    b.nX=float(a.nX);
    b.nY=float(a.nY);
    b.nZ=float(a.nZ);

    return(b);
}

float Int2Float(int a){

    float b;
    b=float(a);
    return(b);

}


