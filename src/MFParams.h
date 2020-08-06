#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams{ 
public:
    /* Convention
 0==yz
 1==xz
 2==xy
 0==up
 1==dn
 index = orb + spin*3 + site*6;

yz_up(site=0),xz_up(site=0),xy_up(site=0), yz_dn(site=0),xz_dn(site=0),xy_dn(site=0)....site(n)...
*/
    // Define Fields
    Mat_3_Complex_doub OParams_;

    Matrix<double> Disorder;

    // Constructor
    MFParams(Parameters& Parameters__, Coordinates&  Coordinates__, mt19937_64& Generator1__ , mt19937_64& Generator2__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }


    double random1();
    double random2();
    void initialize();


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int lx_,ly_,ns_, no_dof_;

    uniform_real_distribution<double> dis1_;//for random fields
    uniform_real_distribution<double> dis2_;//for random disorder

    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};



double MFParams::random1(){

    return dis1_(Generator1_);

}

double MFParams::random2(){

    return dis2_(Generator2_);

}


void MFParams::initialize(){

    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;
    no_dof_=Coordinates_.no_dof_;
    ns_=Coordinates_.ns_;
    int site_i;

    int YZ_,XZ_,XY_;
    YZ_=0; XZ_=1; XY_=2;

    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    // srand(Parameters_.RandomSeed);

    Disorder.resize(lx_,ly_);

    OParams_.resize(ns_);
    for(int i=0;i<ns_;i++){
        OParams_[i].resize(6);
        for(int state=0;state<6;state++){
            OParams_[i][state].resize(6);
        }
    }


    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file<<"#seed="<<Parameters_.RandomDisorderSeed<<
                        " for mt19937_64 Generator is used"<<endl;
    Disorder_conf_file<<"#ix   iy    Dis[ix,iy]"<<endl;

    ofstream Initial_OrderParams_file("Initial_OrderParams_values_generated.txt");

    if(!Parameters_.Read_OPs){

        if(!Parameters_.Create_OPs){
            for(int i=0;i<ns_;i++){
                for(int state1=0;state1<6;state1++){
                    for(int state2=state1;state2<6;state2++){
                        if(state2 !=state1){
                            OParams_[i][state1][state2].real(random1());
                            OParams_[i][state1][state2].imag(random1());
                        }
                        else{
                            OParams_[i][state1][state2].real(random1());
                            OParams_[i][state1][state2].imag(0.0);
                        }
                    }
                }
            }

            Initial_OrderParams_file<<"#seed="<<Parameters_.RandomSeed<<
                                      " for mt19937_64 Generator is used"<<endl;
        }
        else{

            if(Parameters_.Create_OPs_type=="T1_SMixed"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);


                    if(ix%2==0){  //alternating holes in xz and xy

                        if(ix%4==0){
                            spin_direction = 1.0;
                            if(iy%2==0){ //yz,xy

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                            else{ //yz, xz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                        }
                        else{
                            spin_direction = -1.0;
                            if(iy%2==0){//yz,xz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                            else{ //yz,xy
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                        }

                    }
                    else{  //stripe along y for this x with hole in yz
                        spin_direction = -1.0*(pow(-1.0, 1.0*iy));
                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T1_SMixed_perp"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indy(i);
                    iy=Coordinates_.indx(i);


                    if(ix%2==0){  //alternating holes in yz and xy

                        if(ix%4==0){
                            spin_direction = 1.0;
                            if(iy%2==0){ //xz,xy

                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                            else{ //xz, yz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                        }
                        else{
                            spin_direction = -1.0;
                            if(iy%2==0){//yz,xz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                            else{ //xz,xy
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                        }

                    }
                    else{  //stripe along x for this y with hole in xz
                        spin_direction = -1.0*(pow(-1.0, 1.0*iy));
                        OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T1_SMixed_perp_OBC"){ //Not implemented

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indy(i);
                    iy=Coordinates_.indx(i);


                    if(ix%2==0){  //alternating holes in yz and xy

                        if(ix%4==0){
                            spin_direction = 1.0;
                            if(iy%2==0){ //xz,xy

                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                            else{ //xz, yz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                        }
                        else{
                            spin_direction = -1.0;
                            if(iy%2==0){//yz,xz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                            else{ //xz,xy
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                        }

                    }
                    else{  //stripe along x for this y with hole in xz
                        spin_direction = -1.0*(pow(-1.0, 1.0*iy));
                        OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T1p_AFM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0*(pow(-1.0, 1.0*(ix+iy)));

                    if(ix%2==0){  //alternating holes in xz and xy

                       // spin_direction = 1.0;
                        if(iy%2==0){ //yz,xy

                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                            OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                        }
                        else{ //yz, xz
                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                            OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                        }

                    }
                    else{  //stripe along y for this x with hole in yz
                       // spin_direction = -1.0*(pow(-1.0, 1.0*iy));
                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }


            if(Parameters_.Create_OPs_type=="T1p_Spipi_0pi"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);


                    if(ix%2==0){  //alternating holes in xz and xy

                        spin_direction = 1.0;
                        if(iy%2==0){ //yz,xy

                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                            OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                        }
                        else{ //yz, xz
                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                            OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                        }

                    }
                    else{  //stripe along y for this x with hole in yz
                        spin_direction = -1.0*(pow(-1.0, 1.0*iy));
                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T1_FM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0;

                    if(ix%2==0){  //alternating holes in xz and xy

                        if(ix%4==0){
                            if(iy%2==0){ //yz,xy
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                            else{ //yz, xz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                        }
                        else{
                            if(iy%2==0){//yz,xz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                            else{ //yz,xy
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                        }

                    }
                    else{  //stripe along y for this x with hole in yz
                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T1_AFM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.1;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0*(pow(-1.0, 1.0*(ix+iy)));

                    if(ix%2==0){  //alternating holes in xz and xy

                        if(ix%4==0){
                            if(iy%2==0){ //yz,xy
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                            else{ //yz, xz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                            }
                        }
                        else{
                            if(iy%2==0){//yz,xz
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                            else{ //yz,xy
                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                                OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                                OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                                OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                                OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            }
                        }

                    }
                    else{  //stripe along y for this x with hole in yz
                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T1p_AFM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0*(pow(-1.0, 1.0*(ix+iy)));

                    if(ix%2==0){  //alternating holes in xz and xy

                        if(iy%2==0){ //yz,xy
                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                            OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                        }
                        else{ //yz, xz
                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                            OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                            OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);

                            OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                            OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);

                        }

                    }
                    else{  //stripe along y for this x with hole in yz
                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T1_AFM_perp"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                int spin;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin = int((-1.0*(pow(-1.0, 1.0*(ix+iy))) + 1.0)/2.0);

                    if(iy%2==0){

                        if(iy%4==0){
                            if(ix%2==0){
                                OParams_[i][YZ_ + 3*spin][YZ_ + 3*spin]=complex<double>(0.5,0);
                                OParams_[i][XY_ + 3*spin][XY_ + 3*spin]=complex<double>(0.5,0);
                            }
                            else{
                                OParams_[i][YZ_ + 3*spin][YZ_ + 3*spin]=complex<double>(0.5,0);
                                OParams_[i][XZ_ + 3*spin][XZ_ + 3*spin]=complex<double>(0.5,0);
                            }
                        }
                        else{
                            if(ix%2==0){
                                OParams_[i][YZ_ + 3*spin][YZ_ + 3*spin]=complex<double>(0.5,0);
                                OParams_[i][XZ_ + 3*spin][XZ_ + 3*spin]=complex<double>(0.5,0);
                            }
                            else{
                                OParams_[i][YZ_ + 3*spin][YZ_ + 3*spin]=complex<double>(0.5,0);
                                OParams_[i][XY_ + 3*spin][XY_ + 3*spin]=complex<double>(0.5,0);
                            }
                        }

                    }
                    else{
                        OParams_[i][XZ_ + 3*spin][XZ_ + 3*spin]=complex<double>(0.5,0);
                        OParams_[i][XY_ + 3*spin][XY_ + 3*spin]=complex<double>(0.5,0);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T4_AFM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0*(pow(-1.0, 1.0*(ix+iy)));

                    if(iy%2==0){  //stripe along x with e in yz,xy
                        OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                    else{ //stripe along x with e in xz,xy
                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }

            if(Parameters_.Create_OPs_type=="T4_AFM_perp"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indy(i);
                    iy=Coordinates_.indx(i);
                    spin_direction = 1.0*(pow(-1.0, 1.0*(ix+iy)));

                    if(iy%2==0){  //stripe along x with e in yz,xy
                        OParams_[i][YZ_ + 3*UP_][YZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][YZ_ + 3*DOWN_][YZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][YZ_ + 3*UP_][YZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                    else{ //stripe along x with e in xz,xy
                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                        OParams_[i][XZ_ + 3*DOWN_][XZ_ + 3*DOWN_]=complex<double>(OP_temp,0);
                        OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                        OParams_[i][XZ_ + 3*UP_][XZ_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                        OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    }
                }

            }



            if(Parameters_.Create_OPs_type=="T5_AFM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                int ORB_TEMP;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0*(pow(-1.0, 1.0*(ix+iy)));
                    if(spin_direction>0.0){
                        ORB_TEMP=YZ_;
                    }
                    else{
                        ORB_TEMP=XZ_;
                    }

                    OParams_[i][ORB_TEMP + 3*UP_][ORB_TEMP + 3*UP_]=complex<double>(OP_temp,0);
                    OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                    OParams_[i][ORB_TEMP + 3*DOWN_][ORB_TEMP + 3*DOWN_]=complex<double>(OP_temp,0);
                    OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                    OParams_[i][ORB_TEMP + 3*UP_][ORB_TEMP + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);


                }

            }


            if(Parameters_.Create_OPs_type=="T5_FM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                int ORB_TEMP;
                double OP_temp;
                OP_temp=0.4;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0 ; //1.0*(pow(-1.0, 1.0*(ix+iy)));
                    if(pow(-1.0, 1.0*(ix+iy))>0.0){
                        ORB_TEMP=YZ_;
                    }
                    else{
                        ORB_TEMP=XZ_;
                    }

                    OParams_[i][ORB_TEMP + 3*UP_][ORB_TEMP + 3*UP_]=complex<double>(OP_temp,0);
                    OParams_[i][XY_ + 3*UP_][XY_ + 3*UP_]=complex<double>(OP_temp,0);
                    OParams_[i][ORB_TEMP + 3*DOWN_][ORB_TEMP + 3*DOWN_]=complex<double>(OP_temp,0);
                    OParams_[i][XY_ + 3*DOWN_][XY_ + 3*DOWN_]=complex<double>(OP_temp,0);

                    OParams_[i][ORB_TEMP + 3*UP_][ORB_TEMP + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    OParams_[i][XY_ + 3*UP_][XY_ + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);


                }

            }


            if(Parameters_.Create_OPs_type=="T3_FM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                int ORB_TEMP1, ORB_TEMP2;
                double OP_temp;
                OP_temp=0.5;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0; //1.0*(pow(-1.0, 1.0*(ix+iy)));


                    if(pow(-1.0, 1.0*(ix+iy))<0.0){
                        ORB_TEMP1=YZ_;
                        ORB_TEMP2=XZ_;
                    }
                    else{
                       if(iy%2==0){
                           ORB_TEMP1=YZ_;
                           ORB_TEMP2=XY_;
                       }
                       else{
                           ORB_TEMP1=XZ_;
                           ORB_TEMP2=XY_;
                       }

                    }

                    OParams_[i][ORB_TEMP1 + 3*UP_][ORB_TEMP1 + 3*UP_]=complex<double>(OP_temp,0);
                    OParams_[i][ORB_TEMP2 + 3*UP_][ORB_TEMP2 + 3*UP_]=complex<double>(OP_temp,0);
                    OParams_[i][ORB_TEMP1 + 3*DOWN_][ORB_TEMP1 + 3*DOWN_]=complex<double>(OP_temp,0);
                    OParams_[i][ORB_TEMP2 + 3*DOWN_][ORB_TEMP2 + 3*DOWN_]=complex<double>(OP_temp,0);

                    OParams_[i][ORB_TEMP1 + 3*UP_][ORB_TEMP1 + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    OParams_[i][ORB_TEMP2 + 3*UP_][ORB_TEMP2 + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);


                }

            }

            if(Parameters_.Create_OPs_type=="T3_AFM"){

                for(int i=0;i<ns_;i++){
                    for(int state1=0;state1<6;state1++){
                        for(int state2=0;state2<6;state2++){
                            OParams_[i][state1][state2]=zero_complex;
                        }
                    }
                }

                int ix, iy;
                double spin_direction;
                int ORB_TEMP1, ORB_TEMP2;
                double OP_temp;
                OP_temp=0.1;
                for(int i=0;i<ns_;i++){
                    ix=Coordinates_.indx(i);
                    iy=Coordinates_.indy(i);
                    spin_direction = 1.0*(pow(-1.0, 1.0*(ix+iy)));


                    if(pow(-1.0, 1.0*(ix+iy))<0.0){
                        ORB_TEMP1=YZ_;
                        ORB_TEMP2=XZ_;
                    }
                    else{
                       if(iy%2==0){
                           ORB_TEMP1=YZ_;
                           ORB_TEMP2=XY_;
                       }
                       else{
                           ORB_TEMP1=XZ_;
                           ORB_TEMP2=XY_;
                       }

                    }

                    OParams_[i][ORB_TEMP1 + 3*UP_][ORB_TEMP1 + 3*UP_]=complex<double>(OP_temp,0);
                    OParams_[i][ORB_TEMP2 + 3*UP_][ORB_TEMP2 + 3*UP_]=complex<double>(OP_temp,0);
                    OParams_[i][ORB_TEMP1 + 3*DOWN_][ORB_TEMP1 + 3*DOWN_]=complex<double>(OP_temp,0);
                    OParams_[i][ORB_TEMP2 + 3*DOWN_][ORB_TEMP2 + 3*DOWN_]=complex<double>(OP_temp,0);

                    OParams_[i][ORB_TEMP1 + 3*UP_][ORB_TEMP1 + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);
                    OParams_[i][ORB_TEMP2 + 3*UP_][ORB_TEMP2 + 3*DOWN_]=complex<double>(0.0,spin_direction*OP_temp);


                }

            }



        }


    }
    else{
        vector<string> OPstring;
        OPstring.clear();
        OPstring.push_back("All_OP");

        for(int op_no=0;op_no<OPstring.size();op_no++){
            string fl_initial_OP_in = Parameters_.File_OPs_in;
            ifstream file_initial_OP_in(fl_initial_OP_in.c_str());
            string temp1;//,temp2,temp3,temp4,temp5,temp6,temp7;
            int site_temp, x,y, state1_temp, state2_temp;
            double val_real,val_imag;

            for(int i=0;i<7;i++){
                file_initial_OP_in>>temp1;
                cout<<temp1<<"   ";
            }
            cout<<endl;


            // file_initial_OP_in>>temp1>>temp2>>temp3;
            // cout<<temp1<<"  "<<temp2<<" "<<temp3<<endl;

            for(int iy=0;iy<ly_;iy++){
                for(int ix=0;ix<lx_;ix++){
                    site_i=Coordinates_.Nc(ix,iy);
                    for(int state1=0;state1<6;state1++){
                        for(int state2=state1;state2<6;state2++){
                            file_initial_OP_in>>site_temp>>x>>y>>state1_temp>>state2_temp>>val_real>>val_imag;
                            //cout<<ix<<"\t"<<iy<<"\t"<<x<<"\t"<<y<<"\t"<<val<<endl;
                            assert(site_temp==site_i);
                            assert(x==ix);
                            assert(y==iy);
                            assert(state1==state1_temp);
                            assert(state2==state2_temp);
                            OParams_[site_i][state1][state2].real(val_real);
                            OParams_[site_i][state1][state2].imag(val_imag);
                        }
                    }
                }
            }
        }

        Initial_OrderParams_file<<"#OParams are read from "<<Parameters_.File_OPs_in<<" file"<<endl;

    }



    Initial_OrderParams_file<<"#site     lx      ly     state1     state2   OParams_[site][state1][state2].real()       OParams_[site][state1][state2].imag()"<<endl;


    for(int iy=0;iy<ly_;iy++){
        for(int ix=0;ix<lx_;ix++){
            site_i=Coordinates_.Nc(ix,iy);
            for(int state1=0;state1<6;state1++){
                for(int state2=state1;state2<6;state2++){

                    Initial_OrderParams_file<<site_i<<setw(15)<<ix<<setw(15)<<iy<<setw(15)<<state1<<setw(15)<<state2<<setw(15)<<
                                              OParams_[site_i][state1][state2].real()<<setw(15)<<OParams_[site_i][state1][state2].imag()<<endl;
                }
            }
        }
    }


    //Creating Lower Half of Single particle density matrix
    for(int iy=0;iy<ly_;iy++){
        for(int ix=0;ix<lx_;ix++){
            site_i=Coordinates_.Nc(ix,iy);
            for(int state1=0;state1<6;state1++){
                for(int state2=0;state2<state1;state2++){
                    OParams_[site_i][state1][state2] = conj(OParams_[site_i][state2][state1]);
                }
            }
        }
    }


    int temp_i, temp_j;
    if(Parameters_.ReadDisorder==false){
        cout <<"Disorder conf is initialized using random seed given in input file"<< endl;
        for(int j=0;j<ly_;j++){
            for(int i=0;i<lx_;i++){
                //RANDOM Disorder
                Disorder(i,j)=Parameters_.Disorder_Strength*((2.0*random2())-1.0);
                Disorder_conf_file<<i<<"  "<<j<<"  "<<Disorder(i,j)<<endl;
            }
            Disorder_conf_file<<endl;
        }
    }
    else{
        cout<<"Disorder conf is read from file path given in the input file"<<endl;
        ifstream Disorder_in_file(Parameters_.DisorderSeedFile);
        for(int j=0;j<ly_;j++){
            for(int i=0;i<lx_;i++){
                //RANDOM Disorder
                Disorder_in_file>>temp_i>>temp_j>>Disorder(i,j);
                assert(i==temp_i);
                assert(j==temp_j);
                Disorder(i,j)=Parameters_.Disorder_Strength*(Disorder(i,j));
                Disorder_conf_file<<i<<"  "<<j<<"  "<<Disorder(i,j)<<endl;
            }
            Disorder_conf_file<<endl;
        }
    }


} // ----------

#endif
