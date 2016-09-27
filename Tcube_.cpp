#include <iostream>
#include <stdlib.h>	// standard C library function
#include "Space.h"
#include "globle.h"
#include <vector>
#include <math.h>


using namespace std;

void cube_temp(){

    SGrid <float>  sqtemp[31],rad2aavtemp[31];
    bool itobig,itest2,ipore,ionlymoldebug;
    SGrid <float> * sq=sqtemp+15;
    SGrid <float> * rad2aav=rad2aavtemp+15;
    SGrid <int>* ioff;
    string strtmp,strtmp1;
    //Non-standard integer variable, thus necessary
    int epsdim, objecttype;
    //Non-standard real variable, thus necessary
    float modul,modul2, mod2,modx,mody;
    //Non-standard type of variable, thus necessary
    SGrid <float> xa,xb,xc,xd,xp,ddist,xyz2,dxyz,ddxyz,xn,vecdist;
    SGrid <float> tmpvect1,tmpvect2;
    SGrid <float> ectx,vecty,vectz,rad2av,fxn,vtemp;
    SGrid <int>ismin,ismax,idist,idist1,ixyz,itest,ixn,i123;

    //here radprb is not zero only if one wants to map the extended
    //surf. that has to be done with radprb(1) if solute-solvent
    //interface is concerned
    //imedia = medium number for a object
    //a non-zero entry in iepsmp indicates an atom # plus 1 (to
    //properly treat re-entrant mid-points later (15 Aug 93)

    int iboxt,iac,ibox,ii,igrdc,i,imedia,iv,ix,iy,iz,kind,j,k;
    int limmax,lim;
    float alpha,dis2min2,dis2min1,dot,distsq,exrd,dist,rad,rad2;
    float rad2a,rad4,radius,radmax2,rmid,rad5,radp2,radtest,radprobe;
    float radp,tan2,temp,tmp,tmp1,tmp2;




    radprobe=fRadPrb[0];

    cout << "this is in the setout..." << endl;

    if(bVerbose){
        cout << "Starting creating Van der Waals  Epsilon Map" << endl;
    }

    epsdim=iNatom+iNObject+2;
    iboxt=0;
    radmax2=0.0;
    rmid=float((iGrid+1)/2.0);
    itest2=false;


    if (iNatom > 0){
        for(ix=0;ix<=iNatom-1;ix++){
            //Chaged to derive-type array delphipdb (pointers module)
            //radmax2=max(radmax2,delphipdb(ix)%rad3)
            radmax2=max(radmax2,sDelPhiPDB[ix].radius);
        }

        //this is probably the best way to do it,depending on which surf. is desired
        temp=max(radprobe,fExternRadius);

        radmax2=fScale*(radmax2+temp);
        lim=1+radmax2;
        limmax = 12;
        itobig = false;


        if(lim > limmax) {
            itobig=true;
        }

        igrdc = pow((2*lim+1),3.0);
        ioff = new SGrid <int>[igrdc];

        if (!itobig) {

            //radtest= (radmax2 + 0.5* sqrt (3.0))**2
            radtest= pow ((radmax2 + 0.5 * sqrt (3.0)),2.0);
            ibox=0;

            //Strange statement. May allocate or may not
            //allocate array that used later in the program
            //irrespectively of itobig value, thus moved array
            //allocation before if condition



            for(ix = -lim;ix <= lim; ix++){
                for(iy = -lim;iy <= lim; iy++){
                    for(iz = -lim;iz <= lim; iz++){
                        //Replaced by faster operation
                        //Using operations on SGrid <float> and
                        //SGrid <int>type variables defined in module
                        //operators_on_coordinates
                        idist.i=ix;
                        idist.j=iy;
                        idist.k=iz;


                        //cout << "to be finished using SGrid <float> operations in setout:" << endl;
                        //dist=real(idist.dot.idist);
                        dist=idist.i*idist.i+idist.j*idist.j+idist.k*idist.k;

                        //ddist = dist + 0.25 + float(idist);
                        ddist.x = dist + 0.25 + float(idist.i);
                        ddist.y = dist + 0.25 + float(idist.j);
                        ddist.z = dist + 0.25 + float(idist.k);

                        //cout << "to be finished using SGrid <float> operations in setout(to be finished):" << endl;
                        //#################################
                        //if ((dist < radtest)||(ddist.vorlt.radtest))then{
                        if ((dist < radtest)){
                        //###################################

                                ibox=ibox+1;
                                ioff[ibox]=idist;
                        }


                        //cout << "ix,iy,iz:" << ix << iy << iz << endl;
                    }
                }
            }
        }

    }




    //set interiors in OBJECTS
    // OBJECTS are removed


    //set interiors in MOLECULES
    if(itest2||itobig) {
        //write(6,*)'setout method 1',itest2,itobig
        cout << "setout method 1  " << itest2 << "  " << itobig;
    }

    for(iv=0;i<=iNatom-1;i++){

        //restore values
        rad = sDelPhiPDB[iv].radius;

        //Using operations on SGrid <float> and SGrid <int>type variables defined
        //in module operators_on_coordinates
        xn=xn2[iv];

        //if (rad.lt.1.e-6) then
        //    cycle DoATOMS
        //end if

        if(rad < 0.000001){
            continue;
        }


        //scale radius to grid

        //rad=rad*scale; rad5=(rad+0.5)**2; radp=rad+exrad*scale
        //rad=rad+radprobe*scale
        //rad4=(rad+0.5)**2; rad2=rad*rad; radp2=radp*radp

        rad=rad*fScale;
        rad5=pow((rad+0.5),2.0);
        radp=rad+fExternRadius*fScale;
        rad=rad+radprobe*fScale;
        rad4=pow((rad+0.5),2.0);
        rad2=rad*rad;
        radp2=radp*radp;


        //set dielectric map
        //check if sphere sits within limits of box
        itest2=false;


        //Using operations on SGrid <float> and SGrid <int>type
        //variables defined in module operators_on_coordinates
        cout << "to be finished using SGrid <float> operations in setout:" << endl;

        //ismin=int(xn-radmax2-1.0);
        //ismax=int(xn+radmax2+1.0);
        //itest=ismin;
        //ismin=min(ismin,igrid);
        //ismin=max(ismin,1);

        ismin.i=int(xn.x-radmax2-1.0);
        ismin.j=int(xn.y-radmax2-1.0);
        ismin.k=int(xn.z-radmax2-1.0);

        ismax.i=int(xn.x+radmax2+1.0);
        ismax.j=int(xn.y+radmax2+1.0);
        ismax.k=int(xn.z+radmax2+1.0);
        itest=ismin;

        ismin.i=min(ismin.i,iGrid);
        ismin.j=min(ismin.j,iGrid);
        ismin.k=min(ismin.k,iGrid);

        ismin.i=max(ismin.i,1);
        ismin.j=max(ismin.j,1);
        ismin.k=max(ismin.k,1);

        cout << "to be finished using SGrid <float> operations in setout:" << endl;

        //if(itest.vorne.ismin) itest2=.true. //Lin Li: this is to be done;
        itest=ismax;
        ismax.i=min(ismax.i,iGrid) ;
        ismax.j=min(ismax.j,iGrid) ;
        ismax.k=min(ismax.k,iGrid) ;
        ismax.i=max(ismax.i,1);
        ismax.j=max(ismax.j,1);
        ismax.k=max(ismax.k,1);
        //if(itest.vorne.ismax) itest2=.true.  //Lin Li: this is to be done;

        //slow method
        if (itest2||itobig) {
            rad2a = rad2 - 0.25;

            for(iz=ismin.k;iz<=ismax.k;iz++){
                for(iz=ismin.j;iz<=ismax.j;iz++){
                    for(iz=ismin.i;iz<=ismax.i;iz++){
                        //Using operations on SGrid <float> and
                        //SGrid <int>type variables defined in module
                        //operators_on_coordinates
//                        ixyz=int_coord(ix,iy,iz); Lin Li:
//                        dxyz=float(ixyz)-xn;      Lin Li:
//                        distsq=dxyz.dot.dxyz;     Lin Li:
//                        dxyz=dxyz+distsq;         Lin Li:



                        if (dxyz.x < rad2a) {
                            iEpsMap[ix][iy][iz].i=iv+1+iAtomMed[iv]*epsdim;
                        }
                        if (dxyz.y < rad2a) {
                            iEpsMap[ix][iy][iz].j=iv+1+iAtomMed[iv]*epsdim;
                        }
                        if (dxyz.z < rad2a) {
                            iEpsMap[ix][iy][iz].k=iv+1+iAtomMed[iv]*epsdim;
                        }

                        if(distsq < radp2) {
                            bDebMap[ix][iy][iz] = false;
                        }
                    }
                }
            }


        }
        else{ //faster method
        //IT HAS PROBLEMS!!!! Walter (be careful before using
        //also in multidilectric case!!!.and..not.isitmd

            rad2a=rad2-0.25;

            //Using operations on SGrid <float> and SGrid <int>type
            //variables defined in module operators_on_coordinates
            cout << "to be finished using SGrid <float> operations in setout:" << endl;
            //ixn=nint(xn);
            //fxn=float(ixn)-xn;
            //rad2av=rad2a-fxn;

            for(ix=-lim;ix<=lim;ix++){
                //vtemp=real(ix)+fxn;
                //sq(ix)=vtemp*vtemp;
                //rad2aav(ix)=rad2a-vtemp;


            }



            //adjust inter-atom, different epsilon bgps+++04/2004 Walter

            if (iNMedia > 1 && bOnlyMol){ //multiple dielectic constant

                    /*
               do i=1,ibox
                  !2011-05-14  Using operations on SGrid <float> and int_coord
                  !type variables defined in module
                  !operators_on_coordinates
                  i123=ioff(i)
                  ixyz=ixn+i123
                  ix=ixyz%i; iy=ixyz%j; iz=ixyz%k
                  distsq=sq(i123%i)%x+sq(i123%j)%y+sq(i123%k)%z

                  if (distsq.lt.rad2aav(i123%i)%x)  then
                     iac=mod(iepsmp(ix,iy,iz)%i,epsdim)-1

                     if (iac.eq.-1.or.iac.gt.natom) then
                        !occhio! non so cosa fa con i pori!!
                        iepsmp(ix,iy,iz)%i=iv+1+iatmmed(iv)*epsdim
                     else
                        !2011-05-14  Using operations on SGrid <float> and
                        !SGrid <int>type variables defined in module
                        !operators_on_coordinates
                        ddxyz=float(ixyz)-xn; ddxyz%x= ddxyz%x+0.5
                        dis2min1=ddxyz.dot.ddxyz-rad2

                        ddxyz=float(ixyz)-xn2(iac); ddxyz%x= ddxyz%x+0.5
                        dis2min2=ddxyz.dot.ddxyz-&
                                       &(delphipdb(iac)%rad3*scale)**2

                        if (dis2min2.gt.dis2min1) iac=iv
                        iepsmp(ix,iy,iz)%i=iac+1+iatmmed(iac)*epsdim
                     end if
                  end if

                  if (distsq.lt.rad2aav(i123%j)%y) then
                     iac=mod(iepsmp(ix,iy,iz)%j,epsdim)-1
                     if (iac.eq.-1.or.iac.gt.natom) then
                        !occhio! non so cosa fa con i pori!!
                        iepsmp(ix,iy,iz)%j=iv+1+iatmmed(iv)*epsdim
                     else
                        !2011-05-14  Using operations on SGrid <float> and
                        !SGrid <int>type variables defined in module
                        !operators_on_coordinates
                        ddxyz=float(ixyz)-xn; ddxyz%y= ddxyz%y+0.5
                        dis2min1=ddxyz.dot.ddxyz-rad2

                        ddxyz=float(ixyz)-xn2(iac); ddxyz%y= ddxyz%y+0.5
                        dis2min2=ddxyz.dot.ddxyz&
                                      &-(delphipdb(iac)%rad3*scale)**2

                        if (dis2min2.gt.dis2min1) iac=iv
                        iepsmp(ix,iy,iz)%j=iac+1+iatmmed(iac)*epsdim
                     end if
                  end if

                  if (distsq.lt.rad2aav(i123%k)%z) then
                     iac=mod(iepsmp(ix,iy,iz)%k,epsdim)-1
                     if (iac.eq.-1.or.iac.gt.natom) then
                        !occhio! non so cosa fa con i pori!!
                        iepsmp(ix,iy,iz)%k=iv+1+iatmmed(iv)*epsdim
                     else
                        !2011-05-14  Using operations on SGrid <float> and
                        !SGrid <int>type variables defined in module
                        !operators_on_coordinates
                        ddxyz=float(ixyz)-xn; ddxyz%z= ddxyz%z+0.5
                        dis2min1=ddxyz.dot.ddxyz-rad2

                        ddxyz=float(ixyz)-xn2(iac); ddxyz%z=ddxyz%z+0.5
                        dis2min2=ddxyz.dot.ddxyz&
                                      &-(delphipdb(iac)%rad3*scale)**2

                        if (dis2min2.gt.dis2min1) iac=iv
                        iepsmp(ix,iy,iz)%k=iac+1+iatmmed(iac)*epsdim
                     end if
                  end if

                  if(distsq.lt.radp2) idebmap(ix,iy,iz)=.false.
               end do
            */

            }
            else{
                for(i=0;i<=ibox-1;i++){
                    //Using operations on SGrid <float> and int_coord
                    //type variables defined in module
                    //operators_on_coordinates

                    cout << "to be finished using SGrid <float> operations in setout:" << endl;
                    /*
                    i123=ioff(i);
                    ixyz=ixn+i123;
                    ix=ixyz%i;
                    iy=ixyz%j;
                    iz=ixyz%k;
                    distsq = sq(i123%i)%x +sq(i123%j)%y + sq(i123%k)%z;
                    */

                    if(distsq < rad2aav[i123.i].x){
                        iEpsMap[ix][iy][iz].i=iv+1+iAtomMed[iv]*epsdim;
                    }

                    if(distsq < rad2aav[i123.j].y){
                        iEpsMap[ix][iy][iz].j=iv+1+iAtomMed[iv]*epsdim;
                    }

                    if(distsq < rad2aav[i123.k].z){
                        iEpsMap[ix][iy][iz].k=iv+1+iAtomMed[iv]*epsdim;
                    }

                    if (distsq < radp2) {
                        bDebMap[ix][iy][iz]=false;
                    }



                }


            } //end if



        } //end if



    } //end for

    //Array deallocation is made by ordinary F95 statement
    delete ioff;

    if(bVerbose){
        cout << "Ending creating Van der Waals  Epsilon Map " << endl;
    }



    sqtemp[0].x=1.123;

    cout << "sqtemp[0].x:" << sqtemp[0].x << endl;
    cout << "sq[-15].x:" << sq[-15].x << endl;
    cout << "iNatom:" << iNatom << endl;

    cout << "going to quit setout..." << endl;

    return;
}

