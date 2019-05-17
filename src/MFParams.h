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


    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            //RANDOM Disorder
            Disorder(i,j)=Parameters_.Disorder_Strength*((2.0*random2())-1.0);
            Disorder_conf_file<<i<<"  "<<j<<"  "<<Disorder(i,j)<<endl;
        }
        Disorder_conf_file<<endl;
    }


} // ----------

#endif
