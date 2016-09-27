#include <iostream>
#include <stdlib.h>	// standard C library function
#include "Space.h"
#include "globle.h"
#include "grid_space.h"
#include <vector>
#include <cmath>
#include <iostream>

#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"



using namespace std;

void epsmak()
{
    cout << "This is in epsmak:" << endl;

    SGrid <float> amin,amax;


    float fRmid,fStart,fFinish,fRMaxTemp,fCputime;
    int i,j,k;
    int test[3][3], test2[6];

    test[1][1]=6;
    test2[0]=0;
    cout << "######test" << test[test2[0]+1][1] << endl << endl;

    //SExtrema <float> * sLimGridUnit;
    // assign values for all global variables:
    amin.nX=1.0;
    amin.nY=1.0;


    fScale=2.0;
    iTestGloble = 0;
    ibmx=50000000;
    iGrid=65;
    iNObject=2;
    iBoundNum=0;
    bUniformDiel=false;
    bVerbose=true;
    bOnlyMol=true;
    bDebug=false;
    fRadiusMax=1.908;
    iGaussian=0;
    // iNatom=10;
    iNMedia=1;
    fRadPrb[0]=1.4;
    fRadPrb[1]=-1.0;
    iBNumSurf=0;


    sLimGridUnit = new SExtrema <float> [iNObject];
    sLimObject = new SExtrema <float> [iNObject];
    iAtomMed = new int [iNatom];
    //bDebMap = new bool [iGrid,iGrid,iGrid];
    //iEpsMap = new SGrid <int>[iGrid,iGrid,iGrid];
    SGrid <float> fs;
    SGrid <int> is;
    fs=coord(1.3,2.5,6.2);
    is=Float2Int(fs);

    cout << "fs:" << fs.nX << " " << fs.nY << " " << fs.nZ << " " << endl;
    cout << "is:" << is.nX << " " << is.nY << " " << is.nZ << " " << endl;
    cout << "float2int(3.5) " << Float2Int(3.5) << endl;

    for(i=0; i<=iNatom; i++)
    {
        iAtomMed[i] = 1;
    }



    cOldMid.nX=0.1;
    cOldMid.nY=0.2;
    cOldMid.nZ=0.3;
    //end of global variables

    //sLimGridUnit=new SExtrema <float> [iNObjects];
    //sLimGridUnit[0].cMax.x=1.0;
    //sLimGridUnit[0].cMin.x=0.0;
    //sLimGridUnit[1]=2;
    //here limobject is expressed in grid units
    fRMid=(iGrid+1)/2.0;


    for(i=1; i<=iNObject-1; i++)
    {
        //     sLimGridUnit[i].cMin=
        // to be finished using SGrid <float> operations:
        //in module operators_on_coordinates
        //limgunit(ii)%min=(limobject(ii)%min-oldmid)*scale+rmid
        //limgunit(ii)%max=(limobject(ii)%max-oldmid)*scale+rmid

        cout << "to be finished using SGrid <float> operations:" << endl;

        //sLimGridUnit[i].cMin.x=(sLimObject[i].cMin.x-cOldMid.x)*fScale+fRMid;
        //sLimGridUnit[i].cMin.y=(sLimObject[i].cMin.y-cOldMid.y)*fScale+fRMid;
        //sLimGridUnit[i].cMin.z=(sLimObject[i].cMin.z-cOldMid.z)*fScale+fRMid;

        sLimGridUnit[i].nMin=(sLimObject[i].nMin-cOldMid)*fScale+fRMid;

        //sLimGridUnit[i].cMax.x=(sLimObject[i].cMax.x-cOldMid.x)*fScale+fRMid;
        //sLimGridUnit[i].cMax.y=(sLimObject[i].cMax.y-cOldMid.y)*fScale+fRMid;
        //sLimGridUnit[i].cMax.z=(sLimObject[i].cMax.z-cOldMid.z)*fScale+fRMid;

        sLimGridUnit[i].nMax=(sLimObject[i].nMax-cOldMid)*fScale+fRMid;
        //cout << "i is:" << i << endl;

    }

    if(bUniformDiel)
    {
        cout << "not going to calculate boundary elements since" << endl;
        cout << "uniform dielectric" << endl;
        iBoundNum=0;
        return;
    }
    //lepsx.y.z and uepsx.y.z should be the upper and lower limits of
    //the expanded box. if the molecule is smaller than this then
    //reduce leps and upeps accordingly
    //note leps/ueps not yet defined..

    //2011-05-10 Converted to SGrid <float> derived type
    amin=sLimGridUnit[0].nMin;
    amax=sLimGridUnit[0].nMin;

    //find global limits IN GRID UNITS, both, molecule and objects,
    //are considered

    if(iNObject > 1)
    {
        for(i=1; i<=iNObject-1; i++)
        {
            cout << "to be finished using SGrid <float> operations:" << endl;
        }
    }
    fRMaxTemp=fRadiusMax;


    if(fIonStrenth !=0 )
    {
        fRMaxTemp=max(fRMaxTemp,fExternRadius);
    }

    fRMaxTemp=fRMaxTemp*fScale;

    //Using operations on SGrid <float> type variables defined
    //in module operators_on_coordinates

    cout << "to be finished using SGrid <float> operations:" << endl;
    //amin=amin-fRMaxTemp;
    //amax=amax+fRMaxTemp;


    cout << "to be finished using SGrid <float> operations:" << endl;
    //Using operations on SGrid <float> and int_coord
    //type variables defined in module operators_on_coordinates
    //limeps%min=max(int(amin)-2,1)
    //limeps%max=min(int(amax)+2,igrid)

//    LimEps.iMin=max(int(amin)-2,1)
//    LimEps.iMax=min(int(amax)+2,iGrid)


    //hanged to array operations
    // initialize bDebmap and iEpsMap:
    for(i=0; i<=iGrid-1; i++)
    {
        for(j=0; j<=iGrid-1; j++)
        {
            for(k=0; k<=iGrid-1; k++)
            {
                //cout << "i" << i << "j" << j << "k" << k <<endl;
                bDebMap[i][j][k]=true;
            }
        }
    }

    cout << "to be finished using SGrid <float> operations in epsmak:" << endl;
    for(i=0; i<=iGrid-1; i++)
    {
        for(j=0; j<=iGrid-1; j++)
        {
            for(k=0; k<=iGrid-1; k++)
            {
                iEpsMap[i][j][k].nX=0;
                iEpsMap[i][j][k].nY=0;
                iEpsMap[i][j][k].nZ=0;

            }
        }
    }

    //point is out of any kind of object (probably in solution)

    //if radprb is less than half of grid spacing, then use old
    //algorithm (sri 29March 93)
    //The new algorithm should be able to handle all scales; still
    //remains to be tested (Sri Apr 95)


    cout << "print time 1:" << endl;

    if(iGaussian==0)
    {
        cout << "go to setout..." << endl;
        //setout();
    }
    else
    {
        cout << "go to setgaussian..." << endl;
        //setgaussian();
    }

    cout << "print time 2:" << endl;

    cout << "deallocate iLimGridUnit:" << endl;

    delete [] sLimGridUnit;


    if(iGaussian==0)
    {
        cout <<"go to VdwToMs..." << endl;

        // VdwToMs();

    }


    cout << "print time 3:" << endl;


    cout << "amin.nX:" << amin.nX << endl;
    cout << "amin.nY:" << amin.nY << endl;
    cout << "iTestGloble:" << iTestGloble << endl;
    cout << "ibmx:" << ibmx << endl;
    cout << "iGrid:" << iGrid << endl;
    cout << "fRMid:" << fRMid << endl;
    cout << "sLimGridUnit[0].nMax.nX:" << sLimGridUnit[0].nMax.nX << endl;
    cout << "sLimGridUnit[0].nMin.nY:" << sLimGridUnit[0].nMin.nY << endl;
    cout << "cOldMid.nX:" << cOldMid.nX << endl;
    cout << "bUniformDiel:" << bUniformDiel << endl;


}





