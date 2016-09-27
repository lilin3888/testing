#include <iostream>
#include <stdlib.h>	// standard C library function
#include "Protein.h"
#include "Space.h"
#include "globle.h"
#include <vector>

#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include <cmath>
#include <iostream>


using namespace std;

//int add(int a, int b);


int main()
{

    int i;

    iTestGloble=0;

    /*
     cout << "begin!" << endl;

     size_t size = 10;
     int isize =20;
     cout << "test start:" << endl;
     //std::vector <int> v(size);
     //int v[size];

     int v[isize];
     cout << "test start1:" << endl;

     for (i = 0; i < 20; i++) {
     v[i] = i;

     cout << "test " << i << ": *v is: " << *v << " and v[] is: " << v[i] << endl;

     }
     //cout << "size of v is: " << v.size() << endl;
     */

    Protein p;			// construct a new Protein object
    cout << "LinLi 1" << endl;
    p.readPDB("1ABA_2.pdb");		// read PDB file
    cout << "LinLi 2" << endl;
    for (i = 0; i < 10; i++)
    {
        cout << p.atom_p[i]->name << p.atom_p[i]->xyz[0] << endl;
    }

    iNatom = p.atom_p.size();

    cout << iNatom << endl;

//	float fCoord[3][iNatom], fRadius[iNatom], fCharge[iNatom];

    sDelPhiPDB = new delphipdb_struc[iNatom];
    xn1 = new SGrid <float> [iNatom];
    xn2 = new SGrid <float> [iNatom];



    for (i = 0; i <= iNatom-1; i++)
    {

        sDelPhiPDB[i].xyz.nX = p.atom_p[i]->xyz[0];
        sDelPhiPDB[i].xyz.nY = p.atom_p[i]->xyz[1];
        sDelPhiPDB[i].xyz.nZ = p.atom_p[i]->xyz[2];
        sDelPhiPDB[i].charge = p.atom_p[i]->fCharge;
        sDelPhiPDB[i].radius = p.atom_p[i]->fRadius;
    }
    /*
    for (i = 0; i < iNatom; i++) {
    	fCoord[0][i] = p.atom_p[i]->xyz[0];
    	fCoord[1][i] = p.atom_p[i]->xyz[1];
    	fCoord[2][i] = p.atom_p[i]->xyz[2];
    	fCharge[i] = p.atom_p[i]->fCharge;
    	fRadius[i] = p.atom_p[i]->fRadius;
    }
    */
    cout << showpoint << fixed;

    //xn2(i)=(delphipdb(i)%xyz-oldmid)*scale+rmid
    //assign values for xn1,xn2
    for (i = 0; i <= iNatom-1; i++)
    {
        //xn1[i].x=sDelPhiPDB[i].xyz.x;
        //xn1[i].y=sDelPhiPDB[i].xyz.y;
        //xn1[i].z=sDelPhiPDB[i].xyz.z;
        xn1[i]=sDelPhiPDB[i].xyz;

        //xn2[i].x=(sDelPhiPDB[i].xyz.x-cOldMid.x)*fScale+fRMid;
        //xn2[i].y=(sDelPhiPDB[i].xyz.y-cOldMid.y)*fScale+fRMid;
        //xn2[i].z=(sDelPhiPDB[i].xyz.z-cOldMid.z)*fScale+fRMid;
        xn2[i]=(sDelPhiPDB[i].xyz-cOldMid)*fScale+fRMid;

    }


    /*
    for (i = 0; i < 20; i++) {
    	cout << setprecision(3) << setw(8) << fCoord[0][i] << setw(8)
    			<< fCoord[1][i] << setw(8) << fCoord[2][i] << setw(8)
    			<< setprecision(4) << fCharge[i] << setw(8) << fRadius[i]
    			<< endl;
    }
    */
    for (i = 0; i < 20; i++)
    {
        cout << setprecision(3) << setw(8) << sDelPhiPDB[i].xyz.nX << setw(8)
             << sDelPhiPDB[i].xyz.nY << setw(8) << sDelPhiPDB[i].xyz.nZ << setw(8)
             << setprecision(4) << sDelPhiPDB[i].charge << setw(8) << sDelPhiPDB[i].radius
             << endl;
    }


    //###############################  ##############################

    int iNgrid, iNobject;


    iNobject = 1;

    cout << setw(8) << iNobject << setw(8) << iNatom
         << endl;

    epsmak();

    cout << "out of epsmak:" << endl;
    //cout << 'out of epsmak:' << endl;

    cout << "iTestGloble:" << iTestGloble << endl;
    cout << "cOldMid.nX:" << cOldMid.nX << endl;
    if(true) cout << "1==1" << endl;

    //cout << "sLimGridUnit[0].nMax.nX:" << sLimGridUnit[0].nMax.nX << endl;
    //cout << "sLimGridUnit[0].cMin.x:" << sLimGridUnit[0].cMin.x << endl;




}

