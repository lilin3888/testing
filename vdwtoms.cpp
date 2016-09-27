#include <iostream>
#include <stdlib.h>	// standard C library function
#include "Space.h"
#include "globle.h"
#include <vector>
#include <math.h>

#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include <cmath>
#include <iostream>



using namespace std;

void VdwToMs(){


    bool outcb[5][5][5];
    string line;
    int nbra[1000];



    int eps[6],nt,itmp[6],dim1,cont,iaprec,objecttype,imap[4][6],dim,sign,epsdim,isign;

    int kind,eps2[6],j3;
    bool remov;
    string strtmp;


    SGrid <float> itemp, rtemp,xq,dxyz,dr123,dx123,u123;
    SGrid <float> goff[6],xg,x1xyz,s123,xxyyzz;
    SGrid <int>ixyz,it123,ixyz2,jxyz;
    SGrid <int>* ibgrd_temp;
    bool out,nbe[7],intb,exists,flag;

    float offf,cbln,cba,del,dis,dmn,dist,ds2,off,dsr,r0a;
    float s1,s2,s3,x1,radp,prbrd12,prbrd22,rsm,rsm1,prbrd2;
    int ia,i,iac,ibgp,ii,iext,iarv,iacv,iiii,imedia;
    int iiord,ios,it1,it2,it3,iv,ix2,iy2,iz2,ix,iy,iz;
    int iii,iord,j,jj,jjj,jx,jy,jz,limu,liml,kk,m,k,n;
    int mr,mpr,n2,nacc,natm,nmmt,ndv,ncms,ncav,nnn,nnbr;
    int nn,ntt,nmt,n1,Nar,iqqq;

    SGrid <int>* iBnd;


    cout << "this is in the VdwToMs..." << endl;

    iBnd = new SGrid <int>[ibmx];
    r0 = new float [iNatom];
    r02 = new float [iNatom];
    rs2 = new float [iNatom];

    ast = new int [iNatom];

 // should be dynamically allocated later
    int bndeps[65][65][6][2];

    offf=(iGrid+1.)/2.;
    epsdim=iNatom+iNObject+2;



    //imap maps from midpoint position to iepsmap entry positions
    //Changed to array operations
    for(i=0;i<=3;i++){
        for(j=0;j<=5;j++){
            imap[i][j]=0;
        }
    }

    imap[1][4]=-1;
    imap[2][5]=-1;
    imap[3][6]=-1;
    imap[4][1]=1;
    imap[4][2]=2;
    imap[4][3]=3;
    imap[4][4]=1;
    imap[4][5]=2;
    imap[4][6]=3;




    //Changed to array operations
    cout << "need to be re done ... Lin Li " << endl;
     // outcb=.true.; outcb(-1:1,-1:1,-1:1)=.false.
     // nbe=.false.; nbe(1:5)=.true.
    //goff=coord(0.,0.,0.)
    //off=0.5/scale

    //hanged to SGrid <float> derived variable type
    //goff(1)%x=off; goff(2)%y=off; goff(3)%z=off
    //goff(4)%x=-off; goff(5)%y=-off; goff(6)%z=-off
    radpmax=max(fRadPrb[0],fRadPrb[1]);

    //convertion from grid to real coordinates(can also use routine gtoc)
    //Using operations on SGrid <float> and SGrid <int>type
    //variables defined in module operators_on_coordinates

    x1=1.0/fScale;
    x1xyz.x=cOldMid.x-(0.5*x1*(iGrid+1));
    x1xyz.y=cOldMid.y-(0.5*x1*(iGrid+1));
    x1xyz.z=cOldMid.z-(0.5*x1*(iGrid+1));

    //find extrema

    //find global extrema

    cMin.x=6000.0;
    cMin.y=6000.0;
    cMin.z=6000.0;

    cMax.x=-6000.0;
    cMax.y=-6000.0;
    cMax.z=-6000.0;

    //here might be problems: Lin Li:
    for(ii=1;ii<=iNObject-1;ii++){
        cMin.x=min(cMin.x,sLimObject[ii].cMin.x);
        cMin.y=min(cMin.y,sLimObject[ii].cMin.y);
        cMin.z=min(cMin.z,sLimObject[ii].cMin.z);

        cMax.x=min(cMax.x,sLimObject[ii].cMax.x);
        cMax.y=min(cMax.y,sLimObject[ii].cMax.y);
        cMax.z=min(cMax.z,sLimObject[ii].cMax.z);


    }
    //find vanderwaals boundary
    n=0; nn=0; nmt=0; nmmt=0;
    //NB change limits to those of the molecule.
    //set for iepsmp NOT equal to unity



    for(k=LimEps.iMin.k+1;k <= LimEps.iMax.k-1; k++){
        for(j=LimEps.iMin.j+1;j <= LimEps.iMax.j-1; j++){
            for(i=LimEps.iMin.i+1;i <= LimEps.iMax.i-1; i++){

                //one distinguishes between internal,external,
                //internal bgp and external bgp
                iext=0;
                ibgp=0;

                //2011-05-17 Changed to iepsmp to SGrid <int>derived type
                itmp[0]=abs(iEpsMap[i][j][k].i)/epsdim;
                itmp[1]=abs(iEpsMap[i][j][k].j)/epsdim;
                itmp[2]=abs(iEpsMap[i][j][k].k)/epsdim;
                itmp[3]=abs(iEpsMap[i-1][j][k].i)/epsdim;
                itmp[4]=abs(iEpsMap[i][j-1][k].j)/epsdim;
                itmp[5]=abs(iEpsMap[i][j][k-1].k)/epsdim;

                if(itmp[0] == 0) {iext=1;}
                if(itmp[0] != itmp[5]) {ibgp=1;}

                for (cont=1;cont<=5;cont++){
                    if(itmp[cont-1] == 0) iext=1;
                    if(itmp[cont-1] != itmp[cont-2]) ibgp=1;
                }
                //assignement of right values to bndeps according to the
                //point nature
                //from now ibnum is the total number of internal and
                //external boundary grid points

                //here index is messy
                if (ibgp > 0){
                    n=n+1;
                    bndeps[i][j][k][1]=n;
                    bndeps[i][j][k][2]=iext;
                    if (iext > 0) nn=nn+1;
                    iBnd[n].i=i;
                    iBnd[n].j=j;
                    iBnd[n].k=k;

                }
                else{
                    bndeps[i][j][k][0]=0;
                    bndeps[i][j][k][1]=0;

                }
                //if (debug) {
                    //for debug

                //}


            }
        }
    }



    iBoundNum=n;
    iBNumSurf=nn;
    nn=0;

    if(bVerbose){

        cout << "boundary points facing continuum solvent= " << iBNumSurf << endl;
        cout << "total number of boundary points before elab.= " << iBoundNum << endl;
    }


    if(iBoundNum > ibmx){
        cout << "ibnum=  " << iBoundNum << " is greater than ibmx = " << ibmx << endl;
        cout << "increase ibmx in vwtms.f" << endl;
        exit(1);
    }


    if (bDebug) {
        //debug
    }














    cout << "going to quit VdwToMs..." << endl;

}

