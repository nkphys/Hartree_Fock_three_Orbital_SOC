#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"
#include "functions.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

//n, a, lda, ipvt, work, lwork, info
extern "C" void   zgetri_(int *,std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);

//zgetrf_ (&n, &n, &(_TEMP(0,0)), &n, &(ipvt[0]), &info);
extern "C" void   zgetrf_(int *,int *, std::complex<double> *, int *, int *, int *);


//zhetri (character UPLO, integer N, complex*16 dimension( lda, * ) A, integer LDA,
//integer IPIV, complex*16 dimension( * ) WORK, integer INFO)
extern "C" void   zhetri_(char *, int *, std::complex<double> *, int *, int *, std::complex<double> *, int *);


//zhetrf(uplo,n,a,lda,ipiv,work,lwork,info)
extern "C" void   zhetrf_(char *, int *, std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);

//dsptrf(character UPLO, integer N, double precision dimension( * ) AP
//       , integer dimension( * ) IPIV, integer INFO)
extern "C" void dsptrf_(char *, int *, double *, int *, int *);

//dsptri	(character 	UPLO, integer N, double precision dimension( * ) AP,
//            integer dimension( * ) IPIV, double precision dimension( * ) WORK,
//            integer INFO)
extern "C" void dsptri_(char *, int *, double *, int *, double *, int *);


class Observables{
public:

    Observables(Parameters& Parameters__, Coordinates& Coordinates__,
                MFParams& MFParams__, Hamiltonian& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns)
    {
        Initialize();
    }

    void Initialize();

    void OccDensity();
    void Calculate_Akw_t2g();
    void Calculate_Akw_jm();
    void Calculate_Nw_t2g();
    void Calculate_Nw_jm();
    void Calculate_SpinSpincorrelations();
    void Calculate_SpinSpincorrelations_Smartly();
    void Calculate_Orbitalcorrelations_Smartly();
    void Calculate_Excitoncorrelations_Smartly();
    void Calculate_IPR();
    void Calculate_Optical_Conductivity();
    void Get_Non_Interacting_dispersion();
    double Lorentzian(double x, double brd);
    void TotalOccDensity();
    void DensityOfStates();
    complex<double> DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right);
    double DOT_P(Mat_1_doub left, Mat_1_doub right);

    void Total_Energy_Average(double Curr_QuantE, double CurrE);

    void OccDensity(int tlabel);
    void Calculate_Local_Density();
    void Calculate_Order_Params();
    void Calculate_Exciton_Matrix();
    void Calculate_Local_n_orb_resolved();
    void Get_OrderParameters_diffs();
    void Update_OrderParameters(int iter);
    void Update_OrderParameters_Second_Broyden(int iter_count);
    void Invert_Beta();
    void Invert_Beta_double();

    double Omega(int i);

    double BandWidth;
    double nia_t,nib_t,nic_t,n_t;
    Matrix<complex<double>> SiSjQ_, SiSjQ_Mean_, SiSjQ_square_Mean_;
    Matrix<complex<double>> Transformation;
    Matrix<double> SiSj_Mean_, SiSj_square_Mean_;
    double Nematic_order_mean_, Nematic_order_square_mean_;
    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    int lx_,ly_,ns_;
    double dosincr_,tpi_;
    vector<double> nia_,nib_,nic_;
    Matrix<double> SiSj_,dos;
    vector<double> sx_,sy_,sz_;
    Mat_3_Complex_doub OParams;
    Mat_1_doub Local_n_orb_resolved;
    Mat_4_Complex_doub F_Exciton;

    // Declare Fields
    Matrix<double> Sz_obs, Sx_obs, Sy_obs;
    Matrix<double> Local_density_obs;
    double Error_OP_;
    double AVG_Total_Energy, AVG_Total_Energy_sqr;


    //Declare Broyden_Mixing vectors
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;


    //Declarations for Second Broyden Method
    Mat_2_doub _Delta_F;
    Mat_2_doub _u;
    Mat_2_doub _A;
    Mat_1_doub _Fm, _Fm_minus1;
    Mat_1_doub _Delta_OPm, _Delta_OPm_minus1;
    Matrix<double> _Beta;
    Mat_1_doub _cm, _gammam;
    double w_minus1;
    Mat_1_doub w;




};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/



void Observables::Invert_Beta(){
/*
    int n=_Beta.n_row();
    int lda=_Beta.n_col();
    int info;
    vector<int> ipvt;
    ipvt.resize(_Beta.n_col());
    vector<complex<double>> work(n);
    int lwork= n;


    char uplo='U';
    zhetrf_(&uplo, &n, &(_Beta(0,0)),&lda,&(ipvt[0]),&(work[0]),&lwork,&info);
    //cout<<"FACTORIZATION OF MATRIX:"<<endl;
    //_TEMP.print();

    zhetri_(&uplo, &n, &(_Beta(0,0)),&lda,&(ipvt[0]),&(work[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("Inverse: zgetri: failed with info!=0.\n");
    }

    //cout<<"INVERSE OF MATRIX:"<<endl;
    //_TEMP.print();
*/
}

void Observables::Invert_Beta_double(){

    int n=_Beta.n_row();
    int info;
    vector<int> ipvt;
    ipvt.resize(_Beta.n_col());
    vector<double> work(n);

    char uplo='U';
    dsptrf_(&uplo, &n, &(_Beta(0,0)),&(ipvt[0]),&info);
    //cout<<"FACTORIZATION OF MATRIX:"<<endl;
    //_TEMP.print();

    dsptri_(&uplo, &n, &(_Beta(0,0)),&(ipvt[0]),&(work[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("Inverse: dsptri: failed with info!=0.\n");
    }

    //cout<<"INVERSE OF MATRIX:"<<endl;
    //_TEMP.print();

}



complex<double> Observables::DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right){
    complex<double> temp_;
    temp_=zero_complex;

    assert(left.size()==right.size());

    for(int i=0;i<left.size();i++){
        temp_ += conj(left[i])*right[i];
    }
    return temp_;

}

double Observables::DOT_P(Mat_1_doub left, Mat_1_doub right){
    double temp_;
    temp_=0.0;

    assert(left.size()==right.size());

    for(int i=0;i<left.size();i++){
        temp_ += left[i]*right[i];
    }
    return temp_;
}

void Observables::Get_OrderParameters_diffs(){

    Error_OP_=0.0;

    for(int site=0;site<ns_;site++){
        for(int state1=0;state1<6;state1++){
            for(int state2=state1;state2<6;state2++){
                Error_OP_ += abs(OParams[site][state1][state2] - MFParams_.OParams_[site][state1][state2])*
                        abs(OParams[site][state1][state2] - MFParams_.OParams_[site][state1][state2]);
            }
        }
    }

    Error_OP_=sqrt(Error_OP_);

}

void Observables::Update_OrderParameters(int iter){


    //Simple mixing
    double alpha_OP=Parameters_.alpha_OP;

    if(Parameters_.Simple_Mixing==true){

        if(iter==0){
            cout<<"Using Simple Mixing to gain Self-Consistency"<<endl;
        }

        for(int site=0;site<ns_;site++){
            for(int state1=0;state1<6;state1++){
                for(int state2=state1;state2<6;state2++){
                    MFParams_.OParams_[site][state1][state2] = (1-alpha_OP)*MFParams_.OParams_[site][state1][state2]
                            + alpha_OP*OParams[site][state1][state2];
                }
            }
        }

        //      Parameters_.mu_old = (1-alpha_OP)*Parameters_.mu_old + alpha_OP*Parameters_.mus;

    }



    /*
    //Declared in initialize function();
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;
     */

    vector<double> vec_V, vec_U;
    double Denominator_;
    vector<double> vec_L;
    Mat_1_int Offsets_;
    Offsets_.resize(5);
    Offsets_[0]=5;Offsets_[1]=10;Offsets_[2]=14;Offsets_[3]=17;Offsets_[4]=19;


    if(Parameters_.Broyden_Mixing==true){
        //assert(false);

        if(iter==0){
            cout<<"Using Broyden Mixing to gain Self-Consistency"<<endl;


            //Get Jinv_np1
            for(int i=0;i<(36*ns_);i++){
                for(int j=0;j<(36*ns_);j++){
                    if(i==j){
                        Jinv_np1[i][j]=-1.0*alpha_OP;
                    }
                    else{
                        Jinv_np1[i][j]=0.0;
                    }
                }}

            //Get F_n
            for(int site=0;site<ns_;site++){
                for(int state1=0;state1<6;state1++){
                    for(int state2=state1;state2<6;state2++){
                        if(state1==state2){
                            F_n[36*(site) + state2]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
                        }
                        else{
                            F_n[36*(site) + Offsets_[state1] + (state2-state1)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
                            F_n[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).imag();
                        }
                    }
                }
            }


            for(int i=0;i<(36*ns_);i++){
                Delta_x_n[i] =0.0;
                for(int j=0;j<(36*ns_);j++){
                    Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                }
            }


            for(int site=0;site<ns_;site++){
                for(int state1=0;state1<6;state1++){
                    for(int state2=state1;state2<6;state2++){
                        if(state1==state2){
                            MFParams_.OParams_[site][state1][state2] += one_complex*(Delta_x_n[36*(site) + state2]);
                        }
                        else{
                            MFParams_.OParams_[site][state1][state2] += complex<double>(Delta_x_n[36*(site) + Offsets_[state1] + (state2-state1)],
                                    Delta_x_n[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]);
                        }
                    }
                }
            }

            //Copy Jinv_np1 to Jinv_n
            Jinv_n = Jinv_np1;

            //Copy F_n to F_nm1
            F_nm1=F_n;
        }

        else{
            //Get F_n
            for(int site=0;site<ns_;site++){
                for(int state1=0;state1<6;state1++){
                    for(int state2=state1;state2<6;state2++){
                        if(state1==state2){
                            F_n[36*(site) + state2]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
                        }
                        else{
                            F_n[36*(site) + Offsets_[state1] + (state2-state1)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
                            F_n[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).imag();
                        }
                    }
                }
            }

            //Get DeltaF_n
            for (int i=0;i<(36*ns_);i++){
                DeltaF_n[i] = F_n[i] - F_nm1[i];
            }

            //Get vec_V = Jinv_n*DeltaF_n
            vec_V.clear();
            vec_V.resize(36*ns_);

            for(int i=0;i<36*ns_;i++){
                vec_V[i] =0.0;
                for(int j=0;j<36*ns_;j++){
                    vec_V[i] += Jinv_n[i][j]*DeltaF_n[j];
                }
            }

            //Get vec_U = Delta_x_n^dagg*Jinv_n
            vec_U.clear();
            vec_U.resize(36*ns_);

            for(int i=0;i<36*ns_;i++){
                vec_U[i] =0.0;
                for(int j=0;j<36*ns_;j++){
                    vec_U[i] += Delta_x_n[j]*Jinv_n[j][i];
                }
            }

            // Get Denominator_=<Delta_x_n|vec_V>
            Denominator_=0.0;
            for(int i=0;i<36*ns_;i++){
                Denominator_ +=Delta_x_n[i]*vec_V[i];
            }


            //Get vec_L=  Delta_x_n - vec_V;
            vec_L.clear();
            vec_L.resize(36*ns_);
            for(int i=0;i<36*ns_;i++){
                vec_L[i] = Delta_x_n[i] - vec_V[i];
            }


            //Get Mat_Temp [Remember to clear later on];
            Mat_2_doub Mat_Temp;
            Mat_Temp.resize(36*ns_);
            for(int i=0;i<36*ns_;i++){
                Mat_Temp[i].resize(36*ns_);
                for(int j=0;j<36*ns_;j++){
                    Mat_Temp[i][j] = (vec_L[i]*vec_U[j])/(Denominator_);
                }
            }


            //Get Jinv_np1

            for(int i=0;i<36*ns_;i++){
                for(int j=0;j<36*ns_;j++){
                    Jinv_np1[i][j]  = Jinv_n[i][j]  + Mat_Temp[i][j];
                }
            }

            for(int i=0;i<36*ns_;i++){
                Delta_x_n[i] =0.0;
                for(int j=0;j<36*ns_;j++){
                    Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                }
            }


            for(int site=0;site<ns_;site++){
                for(int state1=0;state1<6;state1++){
                    for(int state2=state1;state2<6;state2++){
                        if(state1==state2){
                            MFParams_.OParams_[site][state1][state2] += one_complex*(Delta_x_n[36*(site) + state2]);
                        }
                        else{
                            MFParams_.OParams_[site][state1][state2] += complex<double>(Delta_x_n[36*(site) + Offsets_[state1] + (state2-state1)],
                                    Delta_x_n[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]);
                        }
                    }
                }
            }


            //Copy Jinv_np1 to Jinv_n
            Jinv_n = Jinv_np1;

            //Copy F_n to F_nm1
            F_nm1=F_n;


            //Clear Mat_Temp
            for(int i=0;i<36*ns_;i++){
                Mat_Temp[i].clear();
            }
            Mat_Temp.clear();
        }



        //Checking MFParams_.OParams symmetries

        //        for(int i=0;i<lx_;i++){
        //            for(int j=0;j<ly_;j++){
        //                site=Coordinates_.Nc(i,j);
        //                for(int state_i=0;state_i<6;state_i++){
        //                    for(int state_j=0;state_j<6;state_j++){
        //                    if(state_i==state_j){
        //                    if(MFParams_.OParams[site][state_i][state_j].imag()!=0.0){
        //                        cout<<"local densities order params have problem"<<endl;
        //                    }
        //                      assert(MFParams_.OParams[site][state_i][state_j].imag()==0.0);
        //                    }
        //                    else{
        //                    if(
        //                      abs(MFParams_.OParams[site][state_i][state_j] - conj(MFParams_.OParams[site][state_j][state_i]))>0.00001
        //                       ){
        //                        cout<<"Hermiticity of order params is an issue"<<endl;

        //                        cout<<site<<"   "<<state_i<<"    "<<state_j<<"    "<<MFParams_.OParams[site][state_i][state_j]<<endl;
        //                        cout<<site<<"   "<<state_j<<"    "<<state_i<<"    "<<MFParams_.OParams[site][state_j][state_i]<<endl;
        //                    }
        //                    assert(
        //                          abs(MFParams_.OParams[site][state_i][state_j] - conj(MFParams_.OParams[site][state_j][state_i]))<0.00001
        //                          );
        //                    }

        //                    }}}}



        //Maintain hermiticity of order params, because after many iterations they tend to loose it
        //        for(int i=0;i<lx_;i++){
        //            for(int j=0;j<ly_;j++){
        //                site=Coordinates_.Nc(i,j);
        //                for(int state_i=0;state_i<6;state_i++){
        //                    for(int state_j=state_i;state_j<6;state_j++){
        //                        if(state_j>state_i){
        //                            MFParams_.OParams[site][state_i][state_j]=conj(MFParams_.OParams[site][state_j][state_i]);
        //                        }
        //                        else{
        //                            assert(state_j==state_i);
        //                            MFParams_.OParams[site][state_i][state_j].imag(0.0);
        //                        }
        //                    }}}}


    }



    for(int site=0;site<ns_;site++){
        for(int state1=0;state1<6;state1++){
            for(int state2=state1;state2<6;state2++){
                if(state2!=state1){
                    MFParams_.OParams_[site][state2][state1] = conj(MFParams_.OParams_[site][state1][state2]);
                }
            }
        }
    }



}

void Observables::Update_OrderParameters_Second_Broyden(int iter_count){


    //For details see your own notes at
    //"https://github.com/nkphys/3_ORB_SOC_Hartree_Fock/tree/master/Notes/Modified_Broyden_Method"

    //XXXXXXXXXXLiteratureXXXXXXXXXXXXX
    //Second Broyden is used from "D. D. Johnson, Phys. Rev. B 38, 12807, 1988".
    //Look into this "https://arxiv.org/pdf/0805.4446.pdf" as well.
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    double alpha_OP=Parameters_.alpha_OP;
    double normalization_;
    int site;
    int iter;
    Mat_1_int Offsets_;
    Offsets_.resize(5);
    Offsets_[0]=5;Offsets_[1]=10;Offsets_[2]=14;Offsets_[3]=17;Offsets_[4]=19;


    iter=iter_count%Parameters_.BroydenSecondMethodCounter;


    if(iter==0){
        //****Getting ready for iters>0*********
        _Delta_F.clear();
        _u.clear();
        _A.clear();
        _cm.clear();
        _gammam.clear();
        w_minus1=Parameters_.w_minus1;
        w.clear();
        //************************************

        cout<<"Using Modified Broyden Mixing to gain Self-Consistency"<<endl;

        /*
        //Get Fm
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);

                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=state_i;state_j<6;state_j++){
                        _Fm[site*36 + state_i*6 + state_j]=OParams[site][state_i][state_j]
                                - MFParams_.OParams[site][state_i][state_j];
                    }
                }
            }
        }
        */

        //Get Fm
        for(int site=0;site<ns_;site++){
            for(int state1=0;state1<6;state1++){
                for(int state2=state1;state2<6;state2++){
                    if(state1==state2){
                        _Fm[36*(site) + state2]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
                    }
                    else{
                        _Fm[36*(site) + Offsets_[state1] + (state2-state1)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
                        _Fm[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).imag();
                    }
                }
            }
        }


        for(int j=0;j<36*ns_;j++){
            _Delta_OPm[j] = alpha_OP*_Fm[j];
        }


        /*
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        MFParams_.OParams[site][state_i][state_j] +=_Delta_OPm[site*36 + state_i*6 + state_j];

                    }
                }
            }
        }
        */
        for(int site=0;site<ns_;site++){
            for(int state1=0;state1<6;state1++){
                for(int state2=state1;state2<6;state2++){
                    if(state1==state2){
                        MFParams_.OParams_[site][state1][state2] += one_complex*(_Delta_OPm[36*(site) + state2]);
                    }
                    else{
                        MFParams_.OParams_[site][state1][state2] += complex<double>(_Delta_OPm[36*(site) + Offsets_[state1] + (state2-state1)],
                                _Delta_OPm[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]);
                    }
                }
            }
        }

        //Copy Jinv_np1 to Jinv_n
        _Delta_OPm_minus1 = _Delta_OPm;

        //Copy F_n to F_nm1
        _Fm_minus1=_Fm;

    }

    else{

        w.resize(iter);
        w[iter-1]=Parameters_.wn;

        //Get Fm******************
        /*
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);

                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        _Fm[site*36 + state_i*6 + state_j]=OParams_[site][state_i][state_j]
                                - MFParams_.OParams[site][state_i][state_j];
                    }
                }
            }
        }
        */


        for(int site=0;site<ns_;site++){
            for(int state1=0;state1<6;state1++){
                for(int state2=state1;state2<6;state2++){
                    if(state1==state2){
                        _Fm[36*(site) + state2]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
                    }
                    else{
                        _Fm[36*(site) + Offsets_[state1] + (state2-state1)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
                        _Fm[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).imag();
                    }
                }
            }
        }
        //******************************


        //Get DeltaFm/|DeltaFm|-------------------------//
        _Delta_F.resize(iter);
        _Delta_F[iter-1].resize(36*ns_);

        /*
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);

                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        _Delta_F[iter-1][site*36 + state_i*6 + state_j]=_Fm[site*36 + state_i*6 + state_j]
                                - _Fm_minus1[site*36 + state_i*6 + state_j];
                    }
                }
            }
        }*/
        for(int i_=0;i_<36*ns_;i_++){
            _Delta_F[iter-1][i_]=_Fm[i_] - _Fm_minus1[i_];
        }

        //Getting sqrt(<DeltaFm|DeltaFm>)
        /*
        normalization_=0.0;
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        normalization_ += (conj(_Delta_F[iter-1][site*36 + state_i*6 + state_j])*
                                _Delta_F[iter-1][site*36 + state_i*6 + state_j]).real();
                    }
                }
            }
        }*/

        normalization_=0.0;
        for(int i_=0;i_<36*ns_;i_++){
            normalization_ += ((_Delta_F[iter-1][i_])*_Delta_F[iter-1][i_]);
        }
        normalization_=sqrt(normalization_);



        for(int i_=0;i_<36*ns_;i_++){
            _Delta_F[iter-1][i_]=_Delta_F[iter-1][i_]*(1.0/normalization_);
        }
        //--------------------------------------------//


        //Getting Delta_n/|DeltaFm|-------------------//
        for(int i=0;i<36*ns_;i++){
            _Delta_OPm_minus1[i]=_Delta_OPm_minus1[i]*(1.0/normalization_);
        }
        //-----------------------------------------------//


        //Get u[iter-1]------------------------------//
        _u.resize(iter);
        _u[iter-1].resize(36*ns_);
        for(int i=0;i<36*ns_;i++){
            _u[iter-1][i] = (alpha_OP*_Delta_F[iter-1][i])  +  _Delta_OPm_minus1[i];
        }
        //-------------------------------------------------//


        //UPDATE _A----------------------------//

        Mat_1_doub temp_vec;
        temp_vec.resize(iter);
        _A.push_back(temp_vec);
        temp_vec.clear();

        for(int i=0;i<_A.size();i++){
            _A[i].resize(iter);
        }


        for(int i=0;i<iter;i++){
            if(i==(iter-1)){
                _A[i][i]=w[i]*w[i]*DOT_P(_Delta_F[i],_Delta_F[i]);
            }
            else{
                _A[iter-1][i]=w[iter-1]*w[i]*DOT_P(_Delta_F[i],_Delta_F[iter-1]);
                _A[i][iter-1]=w[iter-1]*w[i]*DOT_P(_Delta_F[iter-1],_Delta_F[i]);
            }
        }
        //---------------------------------------------//

        //Get Beta--------------------------------------//
        _Beta.resize(iter,iter);
        for(int i=0;i<iter;i++){
            for(int j=0;j<iter;j++){
                if(i==j){
                    _Beta(i,j) = (w_minus1*w_minus1) + _A[i][j];
                }
                else{
                    _Beta(i,j) = _A[i][j];
                }
            }

        }

        //  cout<<"Before Inversion:"<<endl;
        //Matrix<complex<double>> Beta_before;
        //Beta_before=_Beta;
        //Beta_before.print();

        Invert_Beta_double();
        for(int i=0;i<_Beta.n_col();i++){
            for(int j=0;j<i;j++){
                _Beta(i,j)=_Beta(j,i);
            }
        }


        //cout<<"After Inversion:"<<endl;
        //Matrix<complex<double>> Beta_after;
        //Beta_after=_Beta;
        //Beta_after.print();

        //cout<<"Identity_check:"<<endl;
        //Matrix<complex<double>> Identity_check;
        //Identity_check=product(Beta_before, Beta_after);
        //Identity_check.print();


        //-----------------------------------------------//

        //Get _cm-------------------------------------------//
        _cm.clear();
        _cm.resize(iter);
        for(int i=0;i<iter;i++){
            _cm[i]=w[i]*DOT_P(_Delta_F[i],_Fm);
        }
        //---------------------------------------------------//

        //Get _gammam------------------------------------------//
        _gammam.clear();
        _gammam.resize(iter);
        for(int l=0;l<iter;l++){
            _gammam[l]=0.0;
            for(int k=0;k<iter;k++){
                _gammam[l] += _cm[k]*_Beta(k,l);
            }
        }
        //--------------------------------------------------//


        //Get _Delta_OPm-----------------------------------------//
        for(int i=0;i<36*ns_;i++){
            _Delta_OPm[i]=0.0;
            for(int n=0;n<iter;n++){
                _Delta_OPm[i] += (-1.0)*w[n]*_gammam[n]*_u[iter-1][i];
            }
        }

        for(int i=0;i<36*ns_;i++){
            _Delta_OPm[i] += alpha_OP*_Fm[i];
        }
        //---------------------------------------------//



        /*
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        MFParams_.OParams[site][state_i][state_j] +=_Delta_OPm[site*36 + state_i*6 + state_j];

                    }
                }
            }
        }
        */

        for(int site=0;site<ns_;site++){
            for(int state1=0;state1<6;state1++){
                for(int state2=state1;state2<6;state2++){
                    if(state1==state2){
                        MFParams_.OParams_[site][state1][state2] += one_complex*(_Delta_OPm[36*(site) + state2]);
                    }
                    else{
                        MFParams_.OParams_[site][state1][state2] += complex<double>(_Delta_OPm[36*(site) + Offsets_[state1] + (state2-state1)],
                                _Delta_OPm[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]);
                    }
                }
            }
        }

        //Copy Jinv_np1 to Jinv_n
        _Delta_OPm_minus1 = _Delta_OPm;

        //Copy F_n to F_nm1
        _Fm_minus1=_Fm;

    }

    //Maintain hermiticity of order params, because after many iterations they tend to loose it

    /*for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            for(int state_i=0;state_i<6;state_i++){
                for(int state_j=state_i;state_j<6;state_j++){
                    if(state_j>state_i){
                        MFParams_.OParams[site][state_i][state_j]=conj(MFParams_.OParams[site][state_j][state_i]);
                    }
                    else{
                        assert(state_j==state_i);
                        MFParams_.OParams[site][state_i][state_j].imag(0.0);
                    }
                }}}}
                */
    for(int site=0;site<ns_;site++){
        for(int state1=0;state1<6;state1++){
            for(int state2=state1;state2<6;state2++){
                if(state2!=state1){
                    MFParams_.OParams_[site][state2][state1] = conj(MFParams_.OParams_[site][state1][state2]);
                }
            }
        }
    }



}

void Observables::Calculate_Local_Density(){

    int c1, c2;
    int site;
    complex<double> value;
    vector<double> avg_den_jm, avg_den_t2g;
    avg_den_jm.resize(6);avg_den_t2g.resize(6);

    ofstream local_density_out_file("local_density_out.txt");

    local_density_out_file<<"#ix   iy   3by2_m3by2 3by2_3by2   3by2_m1by2   3by2_1by2    1by2_m1by2    1by2_1by2   yz_up    xz_up      xy_up    yz_dn    xz_dn      xy_dn"<<endl;



    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            local_density_out_file<<i<<setw(10)<<j;
            site=Coordinates_.Nc(i,j);

            for(int state=0;state<6;state++){
                c1=Coordinates_.Nc_dof(site,state);
                value=zero_complex;
                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    value+=( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                             (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                             );

                }
                local_density_out_file<<setw(15)<<value.real();
                avg_den_jm[state] +=value.real();
            }


            for(int t2g_type=0;t2g_type<6;t2g_type++){

                value=zero_complex;
                for(int jm_type1=0;jm_type1<6;jm_type1++){
                    c1=Coordinates_.Nc_dof_(site,jm_type1);

                    for(int jm_type2=0;jm_type2<6;jm_type2++){
                        c2=Coordinates_.Nc_dof_(site,jm_type2);

                        if( (Transformation(t2g_type,jm_type1) != zero_complex)
                                &&
                                (Transformation(t2g_type,jm_type2) != zero_complex) ){

                            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                                value += conj(Hamiltonian_.Ham_(c2,n))*conj(Transformation(t2g_type,jm_type1))*
                                        Hamiltonian_.Ham_(c1,n)*Transformation(t2g_type,jm_type2)*
                                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                            }

                        }
                    }
                }
                local_density_out_file<<setw(15)<<value.real();
                avg_den_t2g[t2g_type] +=value.real();

            }


            local_density_out_file<<endl;

        }
        local_density_out_file<<endl;
    }

    local_density_out_file<<"#Average densities in same order"<<endl;
    local_density_out_file<<"#"<<setw(15);
    for(int state=0;state<6;state++){
        local_density_out_file<<avg_den_jm[state]/(1.0*ns_)<<setw(15);
    }
    for(int state=0;state<6;state++){
        local_density_out_file<<avg_den_t2g[state]/(1.0*ns_)<<setw(15);
    }
    local_density_out_file<<endl;


}



void Observables::Calculate_IPR(){
    /*
    double IPR;
    int c1, site;
    double eta = 0.001;
    int n_chosen=(Parameters_.ns*Parameters_.Fill*2.0) - 1;
    IPR=0.0;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;
                //  IPR += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                //         abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                //         Lorentzian( Parameters_.mus - Hamiltonian_.eigs_[n], eta);

                IPR += abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen))*
                        abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen));




            }
        }
    }

    cout<<"IPR for state no. "<<n_chosen<<" = "<<IPR<<endl;



    IPR=0.0;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){


                    IPR += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            Lorentzian( Parameters_.mus - Hamiltonian_.eigs_[n], eta);

                }
            }
        }
    }

    cout<<"IPR for near mu(with eta = "<<eta<<") = "<<IPR<<endl;




    string fileout_FermiState="Fermi_state_probability.txt";
    ofstream file_FermiState_out(fileout_FermiState.c_str());
    file_FermiState_out<<"#ix   iy   site   |Psi_{Fermi,up}(ix,iy)|^2 + |Psi_{Fermi,dn}(ix,iy)|^2"<<endl;

    double value;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            value = 0.0;

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;

                value += abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen));
            }

            file_FermiState_out<<i<<"\t"<<j<<"\t"<<site<<"\t"<<value<<endl;

        }
        file_FermiState_out<<endl;
    }





    string fileout_Near_mu="Near_mu_probability.txt";
    ofstream file_Near_mu_out(fileout_Near_mu.c_str());
    file_Near_mu_out<<"#ix   iy   site   sum_{n}Lorentz(near_mu)*|Psi_{n,up}(ix,iy)|^2 + |Psi_{n,dn}(ix,iy)|^2"<<endl;


    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            value = 0.0;

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;
                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    value += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            Lorentzian(Parameters_.mus - Hamiltonian_.eigs_[n], eta);;
                }
            }

            file_Near_mu_out<<i<<"\t"<<j<<"\t"<<site<<"\t"<<value<<endl;

        }
        file_Near_mu_out<<endl;
    }



*/
}


void Observables::Calculate_Local_n_orb_resolved(){
    Local_n_orb_resolved.resize(ns_*6);

    int c1;

    for(int site=0;site<ns_;site++){
        for(int orb=0;orb<3;orb++){
            for(int spin=0;spin<2;spin++){

                c1=Coordinates_.Nc_dof(site,orb+3*spin);
                Local_n_orb_resolved[c1]=0;

                for(int n=0;n<Hamiltonian_.eigs_.size();n++){
                    Local_n_orb_resolved[c1] += (
                                ( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n))*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                                ).real();
                }

            }
        }

    }

}

void Observables::Calculate_Order_Params(){


    int c1,c2;

    int YZ_,XZ_,XY_;
    YZ_=0;XZ_=1;XY_=2;


    OParams.resize(ns_);
    for(int i=0;i<ns_;i++){
        OParams[i].resize(6);
        for(int state=0;state<6;state++){
            OParams[i][state].resize(6);
        }
    }



    for(int site=0;site<ns_;site++){
        for(int state1=0;state1<6;state1++){
            for(int state2=state1;state2<6;state2++){
                OParams[site][state1][state2]=zero_complex;

                for(int n=0;n<Hamiltonian_.eigs_.size();n++){
                    c1=Coordinates_.Nc_dof(site,state1);
                    c2=Coordinates_.Nc_dof(site,state2);
                    OParams[site][state1][state2] += (
                                (( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)))*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                                );

                }
                if(state1==state2){
                    assert(OParams[site][state1][state2].imag()<10e-6);
                    OParams[site][state1][state2].imag(0.0);
                }
                else{
                    OParams[site][state2][state1]=conj(OParams[site][state1][state2]);}


            }
        }
    }


    //checking Hermiticity of order parameters


}

void Observables::Calculate_Exciton_Matrix(){

    int c1,c2;
    int jm_state1,jm_state2;
    double Coherence_Length_11;


    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            for(int type1=0;type1<2;type1++){
                if (type1==0){
                    jm_state1=4; // i.e. 1by2_m1by2
                }
                if(type1==1){
                    jm_state1=5; // i.e. 1by2_1by2
                }

                for(int type2=0;type2<2;type2++){
                    if (type2==0){
                        jm_state2=2; // i.e. 3by2_m1by2
                    }
                    if(type2==1){
                        jm_state2=3; // i.e. 3by2_1by2
                    }

                    F_Exciton[i][j][type1][type2]=zero_complex;


                    for(int orb1=0;orb1<3;orb1++){
                        for(int orb2=0;orb2<3;orb2++){
                            for(int spin1=0;spin1<2;spin1++){
                                for(int spin2=0;spin2<2;spin2++){

                                    if(Transformation(jm_state1,orb1+3*spin1) !=zero_complex
                                            &&
                                            Transformation(jm_state2,orb2+3*spin2) !=zero_complex
                                            ){
                                        c1=Coordinates_.Nc_dof(i,orb1+3*spin1); //up
                                        c2=Coordinates_.Nc_dof(j,orb2+3*spin2); //down

                                        for(int n=0;n<Hamiltonian_.eigs_.size();n++){
                                            F_Exciton[i][j][type1][type2] += (
                                                        (conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n))*
                                                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                                                        ).real();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    Coherence_Length_11=0.0;
    double num_, den_;

    int ix,iy, jx,jy;
    num_=0.0;
    den_=0.0;
    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            ix=Coordinates_.indx(i);iy=Coordinates_.indy(i);
            jx=Coordinates_.indx(j);jy=Coordinates_.indy(j);
            num_ += ( ((jx-ix)*(jx-ix)) + ((jy-iy)*(jy-iy)) )*abs(F_Exciton[i][j][1][1])*abs(F_Exciton[i][j][1][1]);
            den_ += abs(F_Exciton[i][j][1][1])*abs(F_Exciton[i][j][1][1]);
        }
    }

    Coherence_Length_11 = sqrt(num_/den_);


    cout<<"Coherence Length for excitons with m=1/2 : "<< Coherence_Length_11<<endl;

}

void Observables::Calculate_Akw_t2g(){


    //---------Read from input file-----------------------//
    string fileout="Akw_t2g.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.05;
    //omega_min=Hamiltonian_.eigs_[0]-0.5-Parameters_.mus;omega_max=Hamiltonian_.eigs_[6*ns_ -1]+0.5-Parameters_.mus;d_omega=0.0005;
    omega_min=-8;omega_max=8;d_omega=0.01;
	//---------------------------------------------------//


    int UP_=0;
    int DN_=1;
    int YZ_, XZ_, XY_;
    YZ_=0; XZ_=1; XY_=2;

    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Akw_out(fileout.c_str());

    int c1,c2;

    Mat_3_Complex_doub A_YZ, A_XZ,A_XY;
    A_YZ.resize(Parameters_.ns);
    A_XZ.resize(Parameters_.ns);
    A_XY.resize(Parameters_.ns);

    for (int i=0;i<Parameters_.ns;i++){
        A_YZ[i].resize(Parameters_.ns);
        A_XZ[i].resize(Parameters_.ns);
        A_XY[i].resize(Parameters_.ns);

        for(int j=0;j<Parameters_.ns;j++){
            A_YZ[i][j].resize(omega_index_max);
            A_XZ[i][j].resize(omega_index_max);
            A_XY[i][j].resize(omega_index_max);
        }
    }


    complex<double> Nup_check(0,0);
    complex<double> Ndn_check(0,0);

    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            cout<<"Akw for "<<l<<"  "<<j<<" done"<<endl;
            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
                A_YZ[j][l][omega_ind]=zero_complex;
                A_XZ[j][l][omega_ind]=zero_complex;
                A_XY[j][l][omega_ind]=zero_complex;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                    //c= l + or1*ns_ + ns_*orbs_*spin;
                    for(int spin=0;spin<2;spin++){
                        //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                        c1 = Coordinates_.Nc_dof(l,YZ_ + 3*spin);
                        c2 = Coordinates_.Nc_dof(j,YZ_ + 3*spin);
                        A_YZ[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                        c1 = Coordinates_.Nc_dof(l,XZ_ + 3*spin);
                        c2 = Coordinates_.Nc_dof(j,XZ_ + 3*spin);
                        A_XZ[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                        c1 = Coordinates_.Nc_dof(l,XY_ + 3*spin);
                        c2 = Coordinates_.Nc_dof(j,XY_ + 3*spin);
                        A_XY[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);
                    }
                }
            }
        }
    }

    cout << "Nup_check = "<<Nup_check<<endl;
    cout << "Ndn_check = "<<Ndn_check<<endl;

    complex<double> temp_YZ, temp_XZ, temp_XY;

    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------For 1D in x direction-----------------
    ky_i=0;
    for(kx_i=-1;kx_i<=(Parameters_.lx);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    /*
    //--------\Gamma to X-----------------
    ky_i=0;
    for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i=(Parameters_.lx/2);
    for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i=(Parameters_.lx/2) - 1;
    ky_i=(Parameters_.lx/2) - 1;
    for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------
    */


    double k22_offset=0;
    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
            temp_YZ=zero_complex;
            temp_XZ=zero_complex;
            temp_XY=zero_complex;


            for(int j=0;j<ns_;j++){
                for(int l=0;l<ns_;l++){
                    temp_YZ += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_YZ[j][l][omega_ind];

                    temp_XZ += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_XZ[j][l][omega_ind];

                    temp_XY += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_XY[j][l][omega_ind];

                }
            }
            //Use 1:6:7----for gnuplot
            file_Akw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    "
                        <<temp_YZ.real()<<"    "<<temp_XZ.real()<<"    "<<temp_XY.real()<<"    "
                       <<temp_YZ.imag()<<"    "<<temp_XZ.imag()<<"    "<<temp_XY.imag()<<"    "<<endl;

        }
        file_Akw_out<<endl;
    }



}


void Observables::Calculate_Akw_jm(){


    Matrix<complex<double>> _HAM = Hamiltonian_.Ham_;

    //Transformation(jm,orb+3*spin)


    //---------Read from input file-----------------------//
    string fileout="Akw_jm.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.08;
    //omega_min=Hamiltonian_.eigs_[0]-0.5-Parameters_.mus;omega_max=Hamiltonian_.eigs_[6*ns_ - 1]+0.5-Parameters_.mus;d_omega=0.0005;
    omega_min=-8.0;omega_max=8.0;d_omega=0.02;
	//---------------------------------------------------//

    int UP_=0;
    int DN_=1;
    int YZ_, XZ_, XY_;
    YZ_=0; XZ_=1; XY_=2;

    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Akw_out(fileout.c_str());

    int c1,c2;

    Mat_3_Complex_doub A_j3by2_pm3by2, A_j3by2_pm1by2, A_j1by2_pm1by2;
    A_j3by2_pm3by2.resize(Parameters_.ns);
    A_j3by2_pm1by2.resize(Parameters_.ns);
    A_j1by2_pm1by2.resize(Parameters_.ns);

    for (int i=0;i<Parameters_.ns;i++){
        A_j3by2_pm3by2[i].resize(Parameters_.ns);
        A_j3by2_pm1by2[i].resize(Parameters_.ns);
        A_j1by2_pm1by2[i].resize(Parameters_.ns);

        for(int j=0;j<Parameters_.ns;j++){
            A_j3by2_pm3by2[i][j].resize(omega_index_max);
            A_j3by2_pm1by2[i][j].resize(omega_index_max);
            A_j1by2_pm1by2[i][j].resize(omega_index_max);
        }
    }


    complex<double> Nup_check(0,0);
    complex<double> Ndn_check(0,0);

    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            cout<<"Akw for "<<l<<"  "<<j<<" done"<<endl;
            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
                A_j3by2_pm3by2[j][l][omega_ind]=zero_complex;
                A_j3by2_pm1by2[j][l][omega_ind]=zero_complex;
                A_j1by2_pm1by2[j][l][omega_ind]=zero_complex;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                    //c= l + or1*ns_ + ns_*orbs_*spin;
                    for(int spin1=0;spin1<2;spin1++){
                        for(int orb1=0;orb1<3;orb1++){

                            for(int spin2=0;spin2<2;spin2++){
                                for(int orb2=0;orb2<3;orb2++){

                                    //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                                    c1 = Coordinates_.Nc_dof(l,orb1 + 3*spin1);
                                    c2 = Coordinates_.Nc_dof(j,orb2 + 3*spin2);

                                    if(Transformation(1,orb1+3*spin1) !=zero_complex
                                            &&
                                            Transformation(1,orb2+3*spin2) !=zero_complex
                                            ){
                                        A_j3by2_pm3by2[j][l][omega_ind] +=  conj(Transformation(1,orb1+3*spin1))*
                                                conj(_HAM(c1,n))*
                                                Transformation(1,orb2+3*spin2)*
                                                _HAM(c2,n)*
                                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                                    }
                                    if(Transformation(0,orb1+3*spin1) !=zero_complex
                                            &&
                                            Transformation(0,orb2+3*spin2) !=zero_complex
                                            ){
                                        A_j3by2_pm3by2[j][l][omega_ind] +=  conj(Transformation(0,orb1+3*spin1))*
                                                conj(_HAM(c1,n))*
                                                Transformation(0,orb2+3*spin2)*
                                                _HAM(c2,n)*
                                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                                    }

                                    if(Transformation(3,orb1+3*spin1) !=zero_complex
                                            &&
                                            Transformation(3,orb2+3*spin2) !=zero_complex
                                            ){
                                        A_j3by2_pm1by2[j][l][omega_ind] +=  conj(Transformation(3,orb1+3*spin1))*
                                                conj(_HAM(c1,n))*
                                                Transformation(3,orb2+3*spin2)*
                                                _HAM(c2,n)*
                                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                                    }
                                    if(Transformation(2,orb1+3*spin1) !=zero_complex
                                            &&
                                            Transformation(2,orb2+3*spin2) !=zero_complex
                                            ){
                                        A_j3by2_pm1by2[j][l][omega_ind] +=  conj(Transformation(2,orb1+3*spin1))*
                                                conj(_HAM(c1,n))*
                                                Transformation(2,orb2+3*spin2)*
                                                _HAM(c2,n)*
                                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                                    }

                                    if(Transformation(5,orb1+3*spin1) !=zero_complex
                                            &&
                                            Transformation(5,orb2+3*spin2) !=zero_complex
                                            ){
                                        A_j1by2_pm1by2[j][l][omega_ind] +=  conj(Transformation(5,orb1+3*spin1))*
                                                conj(_HAM(c1,n))*
                                                Transformation(5,orb2+3*spin2)*
                                                _HAM(c2,n)*
                                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                                    }
                                    if(Transformation(4,orb1+3*spin1) !=zero_complex
                                            &&
                                            Transformation(4,orb2+3*spin2) !=zero_complex
                                            ){
                                        A_j1by2_pm1by2[j][l][omega_ind] +=  conj(Transformation(4,orb1+3*spin1))*
                                                conj(_HAM(c1,n))*
                                                Transformation(4,orb2+3*spin2)*
                                                _HAM(c2,n)*
                                                Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }



    complex<double> temp_j3by2_pm3by2, temp_j3by2_pm1by2, temp_j1by2_pm1by2;

    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------For 1D in x direction-----------------
    ky_i=0;
    for(kx_i=-1;kx_i<=(Parameters_.lx);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    /*
    //--------\Gamma to X-----------------
    ky_i=0;
    for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i=(Parameters_.lx/2);
    for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i=(Parameters_.lx/2) - 1;
    ky_i=(Parameters_.lx/2) - 1;
    for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------
    */


    double k22_offset=0;
    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
            temp_j3by2_pm3by2=zero_complex;
            temp_j3by2_pm1by2=zero_complex;
            temp_j1by2_pm1by2=zero_complex;


            for(int j=0;j<ns_;j++){
                for(int l=0;l<ns_;l++){
                    temp_j3by2_pm3by2 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_j3by2_pm3by2[j][l][omega_ind];

                    temp_j3by2_pm1by2 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_j3by2_pm1by2[j][l][omega_ind];

                    temp_j1by2_pm1by2 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_j1by2_pm1by2[j][l][omega_ind];

                }
            }
            //Use 1:6:7----for gnuplot
            file_Akw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    "
                        <<temp_j3by2_pm3by2.real()<<"    "<<temp_j3by2_pm1by2.real()<<"    "<<temp_j1by2_pm1by2.real()<<"    "
                       <<temp_j3by2_pm3by2.imag()<<"    "<<temp_j3by2_pm1by2.imag()<<"    "<<temp_j1by2_pm1by2.imag()<<"    "<<endl;

        }
        file_Akw_out<<endl;
    }



    _HAM.resize(1,1);
}


void Observables::Calculate_Nw_t2g(){

    //---------Read from input file-----------------------//
    double omega_min, omega_max, d_omega;
    double eta = 0.01;
    omega_min=Hamiltonian_.eigs_[0]-0.5-Parameters_.mus;
    omega_max=Hamiltonian_.eigs_[6*ns_ - 1]+0.5-Parameters_.mus;
    //omega_min=-1.0;
    //omega_max=1.0;
    d_omega=0.005;
    //---------------------------------------------------//

    int c1;
    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );
    double temp_val ;


    //---------------------------------------------------------------------------------//
    //************************************Nw_jm****************************************//
    //---------------------------------------------------------------------------------//
    /*
    ofstream file_Nw_out(fileout.c_str());


    file_Nw_out<<"#(w-mu)    jm_3by2_m3by2     jm_3by2_3by2     ";
    file_Nw_out<<"jm_3by2_m1by2       jm_3by2_1by2     jm_1by2_m1by2   jm_1by2_1by2"<<endl;

    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        file_Nw_out<<omega_min + (omega_ind*d_omega)<<"       ";


        for(int state_type=0;state_type<6;state_type++){
            temp_val=0.0;

            for(int site=0;site<Coordinates_.ns_;site++){
                c1=Coordinates_.Nc_dof_(site,state_type);

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    temp_val +=  (conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                                  Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta)).real();

                }
            }
            file_Nw_out<<temp_val<<"      ";
        }
        file_Nw_out<<endl;
    }

    file_Nw_out<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

    */
    //---------------------------------------------------------------------------------//
    //*********************************************************************************//



    //---------------------------------------------------------------------------------//
    //************************************Nw_jm****************************************//
    //---------------------------------------------------------------------------------//

    int c2;

    string fileout_t2g="Nw_t2g.txt";
    ofstream file_Nw_out_t2g(fileout_t2g.c_str());

    file_Nw_out_t2g<<"#(w-mu)    yz_up    xz_up      xy_up     ";
    file_Nw_out_t2g<<"yz_dn    xz_dn      xy_dn"<<endl;

    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        file_Nw_out_t2g<<omega_min + (omega_ind*d_omega)<<"       ";


        for(int t2g_type=0;t2g_type<6;t2g_type++){
            temp_val=0.0;

            for(int site=0;site<Coordinates_.ns_;site++){
                c1=Coordinates_.Nc_dof_(site,t2g_type);

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    temp_val += ( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                                  Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta)).real();

                }

            }
            file_Nw_out_t2g<<temp_val<<"      ";
        }
        file_Nw_out_t2g<<endl;
    }

    file_Nw_out_t2g<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

    //---------------------------------------------------------------------------------//
    //*********************************************************************************//






    int n_chosen=Parameters_.Total_Particles - 1;
    cout<<"Gap = "<<Hamiltonian_.eigs_[n_chosen+1] - Hamiltonian_.eigs_[n_chosen]<<endl;

    string fileout_Eigen="Eigen_spectrum.txt";
    ofstream file_Eigen_out(fileout_Eigen.c_str());

    for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
        file_Eigen_out<<n<<"\t"<<Hamiltonian_.eigs_[n]<<"\t"<<(1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)) <<endl;
    }


}


void Observables::Calculate_Nw_jm(){

    //---------Read from input file-----------------------//
    string fileout="Nw_jm.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.01;
    omega_min=Hamiltonian_.eigs_[0]-0.5-Parameters_.mus;
    omega_max=Hamiltonian_.eigs_[6*ns_ - 1]+0.5-Parameters_.mus;
    //omega_min=-1.0;
    //omega_max=1.0;
    d_omega=0.005;
    //---------------------------------------------------//

    int c1;
    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );
    double temp_val ;


    int c2;


    ofstream file_Nw_out(fileout.c_str());

    file_Nw_out<<"#(w-mu)    jm_3by2_m3by2     jm_3by2_3by2     ";
    file_Nw_out<<"jm_3by2_m1by2       jm_3by2_1by2     jm_1by2_m1by2   jm_1by2_1by2"<<endl;

    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        file_Nw_out<<omega_min + (omega_ind*d_omega)<<"       ";


        for(int state_type=0;state_type<6;state_type++){
            temp_val=0.0;


            for(int site=0;site<Coordinates_.ns_;site++){

                for(int orb1=0;orb1<3;orb1++){
                    for(int spin1=0;spin1<2;spin1++){
                        for(int orb2=0;orb2<3;orb2++){
                            for(int spin2=0;spin2<2;spin2++){

                                if(Transformation(state_type,orb1 + 3*spin1) != zero_complex
                                        &&
                                        Transformation(state_type,orb2 + 3*spin2) != zero_complex
                                        ){


                                    c1=Coordinates_.Nc_dof_(site,orb1 + 3*spin1);
                                    c2=Coordinates_.Nc_dof_(site,orb2 + 3*spin2);

                                    for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                                        temp_val += ( conj(Transformation(state_type,orb1 + 3*spin1)*Hamiltonian_.Ham_(c1,n))*Transformation(state_type,orb2 + 3*spin2)*Hamiltonian_.Ham_(c2,n)*
                                                      Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta)).real();

                                    }

                                }

                            }
                        }
                    }
                }

            }
            file_Nw_out<<temp_val<<"      ";
        }
        file_Nw_out<<endl;
    }

    file_Nw_out<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

}

void Observables::Calculate_Optical_Conductivity(){
    /*

    string fileout_sigma_w = "Optical_conductivity.txt";
    ofstream file_sigma_w_out(fileout_sigma_w.c_str());
    file_sigma_w_out<<"#omega   Sigma_xx  Sigma_yy"<<endl;

    //--------------------------------------------------//
    double omega_min, omega_max, d_omega;
    double eta = 0.05;
    omega_min=0.00001;omega_max=10.0;d_omega=0.001;
    //---------------------------------------------------//


    Mat_2_doub PSI_x, PSI_y;
    complex<double> value_x, value_y;
    int ipx, ipy;




    PSI_x.resize(2*ns_);
    PSI_y.resize(2*ns_);
    for(int n=0;n<2*ns_;n++){
        PSI_x[n].resize(2*ns_);
        PSI_y[n].resize(2*ns_);
    }



    for(int n=0;n<2*ns_;n++){
        for(int m=0;m<2*ns_;m++){

            value_x=zero_complex;
            value_y=zero_complex;

            for(int i=0;i<ns_;i++){
                ipx = Coordinates_.neigh(i,0);
                ipy = Coordinates_.neigh(i,2);



                for(int spin=0;spin<2;spin++){
                    value_x += ( conj(Hamiltonian_.Ham_(ipx + (ns_*spin),n))*Hamiltonian_.Ham_(i + (ns_*spin),m) )
                            -  ( conj(Hamiltonian_.Ham_(i + (ns_*spin),n))*Hamiltonian_.Ham_(ipx + (ns_*spin),m) );

                    value_y += ( conj(Hamiltonian_.Ham_(ipy + (ns_*spin),n))*Hamiltonian_.Ham_(i + (ns_*spin),m) )
                            -  ( conj(Hamiltonian_.Ham_(i + (ns_*spin),n))*Hamiltonian_.Ham_(ipy + (ns_*spin),m) );

                }


            }


            PSI_x[n][m] = abs(value_x)*abs(value_x);
            PSI_y[n][m] = abs(value_y)*abs(value_y);


        }
    }




    double sigma_x, sigma_y;
    double omega_val;


    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );


    //cout<<"omega_index_max = "<<omega_index_max<<endl;
    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        omega_val = omega_min + (omega_ind*d_omega);

        //cout<<omega_ind<<endl;

        sigma_x=0.0; sigma_y=0.0;
        for(int n=0;n<2*ns_;n++){
            for(int m=0;m<2*ns_;m++){

                if(n!=m){
                    sigma_x += (PSI_x[n][m])
                            *((1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)))
                            *((1.0/( exp((Parameters_.mus-Hamiltonian_.eigs_[m])*Parameters_.beta ) + 1.0)))
                            *Lorentzian( omega_min + (omega_ind*d_omega) + Hamiltonian_.eigs_[n] - Hamiltonian_.eigs_[m], eta);

                    sigma_y += (PSI_y[n][m])
                            *((1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)))
                            *((1.0/( exp((Parameters_.mus-Hamiltonian_.eigs_[m])*Parameters_.beta ) + 1.0)))
                            *Lorentzian( omega_min + (omega_ind*d_omega) + Hamiltonian_.eigs_[n] - Hamiltonian_.eigs_[m], eta);


                }
            }
        }

        sigma_x = sigma_x*PI*(1.0 - exp(-1.0*Parameters_.beta*(omega_val)))*(1.0/(omega_val*ns_));
        sigma_y = sigma_y*PI*(1.0 - exp(-1.0*Parameters_.beta*(omega_val)))*(1.0/(omega_val*ns_));

        file_sigma_w_out<<omega_val<<"     "<<sigma_x<<"     "<<sigma_y<<endl;

    }



*/
}


void Observables::Calculate_SpinSpincorrelations(){
    /*

    string S2_out = "Local_S2.txt";
    ofstream file_S2_out(S2_out.c_str());
    file_S2_out<<"#site_i    ix   iy   S^2[site_i]"<<endl;


    string Sq_out = "Sq.txt";
    ofstream file_Sq_out(Sq_out.c_str());
    file_Sq_out<<"#qx  qy   qx_index    qy_index   S(qx,qy)"<<endl;

    string SSr_out = "SSr.txt";
    ofstream file_SSr_out(SSr_out.c_str());
    file_SSr_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)  SS[site_i][site_j]"<<endl;

    int c_iup, c_idn, c_jup, c_jdn;

    Mat_2_Complex_doub SS_nm;
    SS_nm.resize(Hamiltonian_.Ham_.n_row());
    for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
        SS_nm[n].resize(Hamiltonian_.Ham_.n_row());
    }


    Mat_2_Complex_doub SS_ri_rj;
    SS_ri_rj.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ri_rj[site_i].resize(ns_);
    }


    for(int site_i=0;site_i<ns_;site_i++){
        for(int site_j=0;site_j<ns_;site_j++){

            c_iup=site_i;
            c_jup=site_j;
            c_idn=site_i + ns_;
            c_jdn=site_j + ns_;

            SS_ri_rj[site_i][site_j]=zero_complex;
            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                for(int m=0;m<Hamiltonian_.Ham_.n_row();m++){

                    if(n==m){

                        //SzSz*f(En)*..
                        SS_ri_rj[site_i][site_j] +=(0.25*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,n))*Hamiltonian_.Ham_(c_jup,n) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,n))*Hamiltonian_.Ham_(c_jdn,n) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,n))*Hamiltonian_.Ham_(c_jdn,n) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,n))*Hamiltonian_.Ham_(c_jup,n) )

                                    );

                        //0.5*(S+S- + S-S+)*f(En)*..
                        SS_ri_rj[site_i][site_j] +=(0.5*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,n))*Hamiltonian_.Ham_(c_jdn,n) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,n))*Hamiltonian_.Ham_(c_jup,n) )
                                    );



                    }

                    else{

                        //SzSz*f(En)*f(Em)..
                        SS_ri_rj[site_i][site_j] +=(0.25*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jup,m) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jdn,m) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jdn,m) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jup,m) )

                                    );


                        //SzSz*f(En)*(1.0 - f(Em))..
                        SS_ri_rj[site_i][site_j] +=(0.25*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0 - (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0)) )*
                                (
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,m)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jup,n) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_iup,m)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jdn,n) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,m)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jdn,n) )
                                    -
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_idn,m)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jup,n) )

                                    );


                        //0.5*(S+S- + S-S+)*f(En)*f(Em)..
                        SS_ri_rj[site_i][site_j] +=(0.5*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_iup,n)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jdn,m) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_idn,n)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jup,m) )

                                    );

                        //0.5*(S+S- + S-S+)*f(En)*(1-f(Em))..
                        SS_ri_rj[site_i][site_j] +=(0.5*one_complex)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0)) )*
                                (
                                    (conj(Hamiltonian_.Ham_(c_idn,n))*Hamiltonian_.Ham_(c_iup,m)*
                                     conj(Hamiltonian_.Ham_(c_jup,m))*Hamiltonian_.Ham_(c_jdn,n) )
                                    +
                                    (conj(Hamiltonian_.Ham_(c_iup,n))*Hamiltonian_.Ham_(c_idn,m)*
                                     conj(Hamiltonian_.Ham_(c_jdn,m))*Hamiltonian_.Ham_(c_jup,n) )

                                    );


                    }


                }



            }



            cout<<site_i<<"\t"<<site_j<<" done"<<endl;
        }

    }


    int site_i, site_j;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);

            for(int jx=0;jx<lx_;jx++){
                for(int jy=0;jy<ly_;jy++){
                    site_j=Coordinates_.Nc(jx,jy);
                    file_SSr_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<site_j<<"\t"<<
                                  jx<<"\t"<<jy<<"\t"<<
                                  real(SS_ri_rj[site_i][site_j])<<"\t"<<
                                  imag(SS_ri_rj[site_i][site_j])<<endl;
                }
            }

        }
    }


    double qx_, qy_;
    complex<double> value;



    for(int qx=0;qx<lx_;qx++){
        for(int qy=0;qy<ly_;qy++){
            value = zero_complex;
            qx_ = (2.0*qx*PI)*(1.0/(1.0*lx_));
            qy_ = (2.0*qy*PI)*(1.0/(1.0*ly_));

            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){
                    site_i=Coordinates_.Nc(ix,iy);

                    for(int jx=0;jx<lx_;jx++){
                        for(int jy=0;jy<ly_;jy++){
                            site_j=Coordinates_.Nc(jx,jy);

                            value += one_complex*(exp(iota_complex*( (qx_*(ix - jx)) + (qy_*(iy-jy)) )))*
                                    SS_ri_rj[site_i][site_j]*(1.0/(1.0*ns_));

                        }}
                }}

            file_Sq_out<<qx_<<"\t"<<qy_<<"\t"<<qx<<"\t"<<qy<<"\t"<<real(value)<<"\t"<<imag(value)<<endl;

        }
    }



    complex<double> Avg_S2=0.0;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);
            file_S2_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<real(SS_ri_rj[site_i][site_i])
                      <<"\t"<<imag(SS_ri_rj[site_i][site_i])<<endl;

            Avg_S2 += one_complex*(SS_ri_rj[site_i][site_i]);
        }}

    cout<<"Avg Local Moment (S^2) = "<<real(Avg_S2)/(1.0*ns_)<<"\t"<<imag(Avg_S2)/(1.0*ns_)<<endl;


*/
}


void Observables::Calculate_SpinSpincorrelations_Smartly(){


    string S2_out = "Local_S2.txt";
    ofstream file_S2_out(S2_out.c_str());
    file_S2_out<<"#site_i    ix   iy   S^2[site_i]"<<endl;


    string Sq_out = "Sq.txt";
    ofstream file_Sq_out(Sq_out.c_str());
    file_Sq_out<<"#qx  qy   qx_index    qy_index   S(qx,qy)"<<endl;

    string SSr_out = "SSr.txt";
    ofstream file_SSr_out(SSr_out.c_str());
    file_SSr_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)  SS[site_i][site_j]"<<endl;

    int c_iup, c_idn, c_jup, c_jdn;
    int ci,cj, spin_i,spin_j;

    Mat_4_Complex_doub UP_UP_Fermi, DOWN_DOWN_Fermi, UP_DOWN_Fermi, DOWN_UP_Fermi ;
    Mat_4_Complex_doub UP_UP_1mFermi, DOWN_DOWN_1mFermi, UP_DOWN_1mFermi, DOWN_UP_1mFermi ;

    UP_UP_Fermi.resize(ns_);DOWN_DOWN_Fermi.resize(ns_);UP_DOWN_Fermi.resize(ns_);DOWN_UP_Fermi.resize(ns_);
    UP_UP_1mFermi.resize(ns_);DOWN_DOWN_1mFermi.resize(ns_);UP_DOWN_1mFermi.resize(ns_);DOWN_UP_1mFermi.resize(ns_);

    for(int n=0;n<ns_;n++){
        UP_UP_Fermi[n].resize(ns_);DOWN_DOWN_Fermi[n].resize(ns_);
        UP_DOWN_Fermi[n].resize(ns_);DOWN_UP_Fermi[n].resize(ns_);
        UP_UP_1mFermi[n].resize(ns_);DOWN_DOWN_1mFermi[n].resize(ns_);
        UP_DOWN_1mFermi[n].resize(ns_);DOWN_UP_1mFermi[n].resize(ns_);

        for(int m=0;m<ns_;m++){
            UP_UP_Fermi[n][m].resize(3);DOWN_DOWN_Fermi[n][m].resize(3);
            UP_DOWN_Fermi[n][m].resize(3);DOWN_UP_Fermi[n][m].resize(3);
            UP_UP_1mFermi[n][m].resize(3);DOWN_DOWN_1mFermi[n][m].resize(3);
            UP_DOWN_1mFermi[n][m].resize(3);DOWN_UP_1mFermi[n][m].resize(3);
            for(int orb=0;orb<3;orb++){
                UP_UP_Fermi[n][m][orb].resize(3);DOWN_DOWN_Fermi[n][m][orb].resize(3);
                UP_DOWN_Fermi[n][m][orb].resize(3);DOWN_UP_Fermi[n][m][orb].resize(3);
                UP_UP_1mFermi[n][m][orb].resize(3);DOWN_DOWN_1mFermi[n][m][orb].resize(3);
                UP_DOWN_1mFermi[n][m][orb].resize(3);DOWN_UP_1mFermi[n][m][orb].resize(3);
            }
        }
    }


    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            for(int orbi=0;orbi<3;orbi++){
                for(int orbj=0;orbj<3;orbj++){

                    UP_UP_Fermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_DOWN_Fermi[i][j][orbi][orbj] = zero_complex;
                    UP_DOWN_Fermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_UP_Fermi[i][j][orbi][orbj] = zero_complex;

                    UP_UP_1mFermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_DOWN_1mFermi[i][j][orbi][orbj] = zero_complex;
                    UP_DOWN_1mFermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_UP_1mFermi[i][j][orbi][orbj] = zero_complex;


                    for(int n=0;n<Hamiltonian_.eigs_.size();n++){

                        spin_i=0;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 3*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 3*spin_j);
                        UP_UP_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=1;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 3*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 3*spin_j);
                        DOWN_DOWN_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=0;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 3*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 3*spin_j);
                        UP_DOWN_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=1;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 3*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 3*spin_j);
                        DOWN_UP_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=0;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 3*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 3*spin_j);
                        UP_UP_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=1;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 3*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 3*spin_j);
                        DOWN_DOWN_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=0;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 3*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 3*spin_j);
                        UP_DOWN_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=1;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 3*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 3*spin_j);
                        DOWN_UP_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                    }
                }
            }
        }
    }

    Mat_4_Complex_doub SS_ri_rj;
    SS_ri_rj.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ri_rj[site_i].resize(ns_);
        for(int site_j=0;site_j<ns_;site_j++){
            SS_ri_rj[site_i][site_j].resize(3);
            for(int orb=0;orb<3;orb++){
                SS_ri_rj[site_i][site_j][orb].resize(3);
            }
        }
    }


    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            for(int orbi=0;orbi<3;orbi++){
                for(int orbj=0;orbj<3;orbj++){

                    SS_ri_rj[i][j][orbi][orbj]=zero_complex;

                    //SzSz..
                    SS_ri_rj[i][j][orbi][orbj] +=(0.25*one_complex)*(
                                (UP_UP_Fermi[i][j][orbi][orbj]*UP_UP_1mFermi[j][i][orbj][orbi])
                                - (UP_DOWN_Fermi[i][j][orbi][orbj]*DOWN_UP_1mFermi[j][i][orbj][orbi])
                                + (DOWN_DOWN_Fermi[i][j][orbi][orbj]*DOWN_DOWN_1mFermi[j][i][orbj][orbi])
                                - (DOWN_UP_Fermi[i][j][orbi][orbj]*UP_DOWN_1mFermi[j][i][orbj][orbi])
                                + (UP_UP_Fermi[i][i][orbi][orbi]*UP_UP_Fermi[j][j][orbj][orbj])
                                - (UP_UP_Fermi[i][i][orbi][orbi]*DOWN_DOWN_Fermi[j][j][orbj][orbj])
                                + (DOWN_DOWN_Fermi[i][i][orbi][orbi]*DOWN_DOWN_Fermi[j][j][orbj][orbj])
                                - (DOWN_DOWN_Fermi[i][i][orbi][orbi]*UP_UP_Fermi[j][j][orbj][orbj])
                                );


                    //0.5*(S+S- + S-S+)
                    SS_ri_rj[i][j][orbi][orbj] +=(0.5*one_complex)*(

                                (UP_UP_Fermi[i][j][orbi][orbj]*DOWN_DOWN_1mFermi[j][i][orbj][orbi])
                                + (DOWN_DOWN_Fermi[i][j][orbi][orbj]*UP_UP_1mFermi[j][i][orbj][orbi])
                                + (UP_DOWN_Fermi[i][i][orbi][orbi]*DOWN_UP_Fermi[j][j][orbj][orbj])
                                + (DOWN_UP_Fermi[i][i][orbi][orbi]*UP_DOWN_Fermi[j][j][orbj][orbj])

                                );


                    //cout<<i<<"\t"<<j<<" done"<<endl;
                }
            }
        }

    }



    Mat_2_Complex_doub SS_ij;
    SS_ij.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ij[site_i].resize(ns_);
    }

    for(int site_i=0;site_i<ns_;site_i++){
        for(int site_j=0;site_j<ns_;site_j++){
            SS_ij[site_i][site_j]=zero_complex;
            for(int orbi=0;orbi<3;orbi++){
                for(int orbj=0;orbj<3;orbj++){
                    SS_ij[site_i][site_j] += SS_ri_rj[site_i][site_j][orbi][orbj];
                }
            }
        }
    }



    int site_i, site_j;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);

            for(int jx=0;jx<lx_;jx++){
                for(int jy=0;jy<ly_;jy++){
                    site_j=Coordinates_.Nc(jx,jy);
                    file_SSr_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<site_j<<"\t"<<
                                  jx<<"\t"<<jy<<"\t"<<
                                  real(SS_ij[site_i][site_j])<<"\t"<<
                                  imag(SS_ij[site_i][site_j])<<endl;


                }
            }

        }
    }


    double qx_, qy_;
    complex<double> value;


    for(int qy=0;qy<ly_;qy++){
        for(int qx=0;qx<lx_;qx++){
            value = zero_complex;
            qx_ = (2.0*qx*PI)*(1.0/(1.0*lx_));
            qy_ = (2.0*qy*PI)*(1.0/(1.0*ly_));

            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){
                    site_i=Coordinates_.Nc(ix,iy);

                    for(int jx=0;jx<lx_;jx++){
                        for(int jy=0;jy<ly_;jy++){
                            site_j=Coordinates_.Nc(jx,jy);

                            value += one_complex*(exp(iota_complex*( (qx_*(ix - jx)) + (qy_*(iy-jy)) )))*
                                    SS_ij[site_i][site_j]*(1.0/(1.0*ns_));

                        }}
                }}

            file_Sq_out<<qx_<<"\t"<<qy_<<"\t"<<qx<<"\t"<<qy<<"\t"<<real(value)<<"\t"<<imag(value)<<endl;

        }
        file_Sq_out<<endl;

    }



    complex<double> Avg_S2=0.0;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);
            file_S2_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<real(SS_ij[site_i][site_i])
                      <<"\t"<<imag(SS_ij[site_i][site_i])<<endl;

            Avg_S2 += one_complex*(SS_ij[site_i][site_i]);
        }
        file_S2_out<<endl;
    }

    cout<<"Avg Local Moment (S^2) = "<<real(Avg_S2)/(1.0*ns_)<<"\t"<<imag(Avg_S2)/(1.0*ns_)<<endl;



}



void Observables::Calculate_Orbitalcorrelations_Smartly(){


    //    0==yz
    //    1==xz
    //    2==xy
    //    0==up
    //    1==dn

    int YZ_,XZ_,XY_;
    YZ_=0; XZ_=1; XY_=2;


    string L2_out = "Local_L2.txt";
    ofstream file_L2_out(L2_out.c_str());
    file_L2_out<<"#site_i    ix   iy   L^2[site_i]"<<endl;


    string Lq_out = "Lq.txt";
    ofstream file_Lq_out(Lq_out.c_str());
    file_Lq_out<<"#qx  qy   qx_index    qy_index   L(qx,qy)"<<endl;

    string LLr_out = "LLr.txt";
    ofstream file_LLr_out(LLr_out.c_str());
    file_LLr_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)  LL[site_i][site_j]"<<endl;


    int ci,cj;

    Mat_4_Complex_doub Fermi_, Inv_Fermi_ ;

    Fermi_.resize(ns_); Inv_Fermi_.resize(ns_);
    for(int n=0;n<ns_;n++){
        Fermi_[n].resize(ns_);
        Inv_Fermi_[n].resize(ns_);
        for(int m=0;m<ns_;m++){
            Fermi_[n][m].resize(6);
            Inv_Fermi_[n][m].resize(6);
            for(int state=0;state<6;state++){
                Fermi_[n][m][state].resize(6);
                Inv_Fermi_[n][m][state].resize(6);
            }
        }
    }

    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            for(int alpha=0;alpha<3;alpha++){
                for(int sigma=0;sigma<2;sigma++){
                    for(int alpha_p=0;alpha_p<3;alpha_p++){
                        for(int sigma_p=0;sigma_p<2;sigma_p++){

                            ci=Coordinates_.Nc_dof(i,alpha + 3*sigma);
                            cj=Coordinates_.Nc_dof(j,alpha_p + 3*sigma_p);

                            Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p]=zero_complex;
                            Inv_Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p]=zero_complex;
                            for(int n=0;n<Hamiltonian_.eigs_.size();n++){
                                Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p] += conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                                Inv_Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p] += conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                        (1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                                //  Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p] += conj(Hamiltonian_.Ham_(n,ci))*Hamiltonian_.Ham_(n,cj)*
                                //            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                                //  Inv_Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p] += conj(Hamiltonian_.Ham_(n,ci))*Hamiltonian_.Ham_(n,cj)*
                                //         (1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                            }

                        }
                    }
                }
            }
        }
    }



    Mat_2_Complex_doub LL_ri_rj;
    LL_ri_rj.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        LL_ri_rj[site_i].resize(ns_);
    }





    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            LL_ri_rj[i][j]=zero_complex;

            for(int sigma=0;sigma<2;sigma++){
                for(int sigma_p=0;sigma_p<2;sigma_p++){

                    //<Lz[i]Lz[j]>
                    LL_ri_rj[i][j] += (Fermi_[i][i][YZ_ + 3*sigma][XZ_ + 3*sigma]*
                            Fermi_[j][j][XZ_ + 3*sigma_p][YZ_ + 3*sigma_p])
                            +
                            (Fermi_[i][j][YZ_ + 3*sigma][YZ_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][XZ_ + 3*sigma_p][XZ_ + 3*sigma])
                            -
                            (Fermi_[i][i][YZ_ + 3*sigma][XZ_ + 3*sigma]*
                            Fermi_[j][j][YZ_ + 3*sigma_p][XZ_ + 3*sigma_p])
                            -
                            (Fermi_[i][j][YZ_ + 3*sigma][XZ_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][YZ_ + 3*sigma_p][XZ_ + 3*sigma])
                            +
                            (Fermi_[i][i][XZ_ + 3*sigma][YZ_ + 3*sigma]*
                            Fermi_[j][j][YZ_ + 3*sigma_p][XZ_ + 3*sigma_p])
                            +
                            (Fermi_[i][j][XZ_ + 3*sigma][XZ_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][YZ_ + 3*sigma_p][YZ_ + 3*sigma])
                            -
                            (Fermi_[i][i][XZ_ + 3*sigma][YZ_ + 3*sigma]*
                            Fermi_[j][j][XZ_ + 3*sigma_p][YZ_ + 3*sigma_p])
                            -
                            (Fermi_[i][j][XZ_ + 3*sigma][YZ_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][XZ_ + 3*sigma_p][YZ_ + 3*sigma])
                            ;



                    //<Lx[i]Lx[j]>
                    LL_ri_rj[i][j] += (Fermi_[i][i][XZ_ + 3*sigma][XY_ + 3*sigma]*
                            Fermi_[j][j][XY_ + 3*sigma_p][XZ_ + 3*sigma_p])
                            +
                            (Fermi_[i][j][XZ_ + 3*sigma][XZ_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][XY_ + 3*sigma_p][XY_ + 3*sigma])
                            -
                            (Fermi_[i][i][XZ_ + 3*sigma][XY_ + 3*sigma]*
                            Fermi_[j][j][XZ_ + 3*sigma_p][XY_ + 3*sigma_p])
                            -
                            (Fermi_[i][j][XZ_ + 3*sigma][XY_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][XZ_ + 3*sigma_p][XY_ + 3*sigma])
                            +
                            (Fermi_[i][i][XY_ + 3*sigma][XZ_ + 3*sigma]*
                            Fermi_[j][j][XZ_ + 3*sigma_p][XY_ + 3*sigma_p])
                            +
                            (Fermi_[i][j][XY_ + 3*sigma][XY_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][XZ_ + 3*sigma_p][XZ_ + 3*sigma])
                            -
                            (Fermi_[i][i][XY_ + 3*sigma][XZ_ + 3*sigma]*
                            Fermi_[j][j][XY_ + 3*sigma_p][XZ_ + 3*sigma_p])
                            -
                            (Fermi_[i][j][XY_ + 3*sigma][XZ_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][XY_ + 3*sigma_p][XZ_ + 3*sigma])
                            ;




                    //<Ly[i]Ly[j]>
                    LL_ri_rj[i][j] += (Fermi_[i][i][XY_ + 3*sigma][YZ_ + 3*sigma]*
                            Fermi_[j][j][YZ_ + 3*sigma_p][XY_ + 3*sigma_p])
                            +
                            (Fermi_[i][j][XY_ + 3*sigma][XY_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][YZ_ + 3*sigma_p][YZ_ + 3*sigma])
                            -
                            (Fermi_[i][i][XY_ + 3*sigma][YZ_ + 3*sigma]*
                            Fermi_[j][j][XY_ + 3*sigma_p][YZ_ + 3*sigma_p])
                            -
                            (Fermi_[i][j][XY_ + 3*sigma][YZ_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][XY_ + 3*sigma_p][YZ_ + 3*sigma])
                            +
                            (Fermi_[i][i][YZ_ + 3*sigma][XY_ + 3*sigma]*
                            Fermi_[j][j][XY_ + 3*sigma_p][YZ_ + 3*sigma_p])
                            +
                            (Fermi_[i][j][YZ_ + 3*sigma][YZ_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][XY_ + 3*sigma_p][XY_ + 3*sigma])
                            -
                            (Fermi_[i][i][YZ_ + 3*sigma][XY_ + 3*sigma]*
                            Fermi_[j][j][YZ_ + 3*sigma_p][XY_ + 3*sigma_p])
                            -
                            (Fermi_[i][j][YZ_ + 3*sigma][XY_ + 3*sigma_p]*
                            Inv_Fermi_[j][i][YZ_ + 3*sigma_p][XY_ + 3*sigma])
                            ;



                }
            }
        }

    }


    int site_i, site_j;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);

            for(int jx=0;jx<lx_;jx++){
                for(int jy=0;jy<ly_;jy++){
                    site_j=Coordinates_.Nc(jx,jy);
                    file_LLr_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<site_j<<"\t"<<
                                  jx<<"\t"<<jy<<"\t"<<
                                  real(LL_ri_rj[site_i][site_j])<<"\t"<<
                                  imag(LL_ri_rj[site_i][site_j])<<endl;


                }
            }

        }
    }


    double qx_, qy_;
    complex<double> value;


    for(int qy=0;qy<ly_;qy++){
        for(int qx=0;qx<lx_;qx++){
            value = zero_complex;
            qx_ = (2.0*qx*PI)*(1.0/(1.0*lx_));
            qy_ = (2.0*qy*PI)*(1.0/(1.0*ly_));

            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){
                    site_i=Coordinates_.Nc(ix,iy);

                    for(int jx=0;jx<lx_;jx++){
                        for(int jy=0;jy<ly_;jy++){
                            site_j=Coordinates_.Nc(jx,jy);

                            value += one_complex*(exp(iota_complex*( (qx_*(ix - jx)) + (qy_*(iy-jy)) )))*
                                    LL_ri_rj[site_i][site_j]*(1.0/(1.0*ns_));

                        }}
                }}

            file_Lq_out<<qx_<<"\t"<<qy_<<"\t"<<qx<<"\t"<<qy<<"\t"<<real(value)<<"\t"<<imag(value)<<endl;

        }
        file_Lq_out<<endl;

    }



    complex<double> Avg_L2=0.0;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);
            file_L2_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<real(LL_ri_rj[site_i][site_i])
                      <<"\t"<<imag(LL_ri_rj[site_i][site_i])<<endl;

            Avg_L2 += one_complex*(LL_ri_rj[site_i][site_i]);
        }
        file_L2_out<<endl;
    }

    cout<<"Avg Local Moment (L^2) = "<<real(Avg_L2)/(1.0*ns_)<<"\t"<<imag(Avg_L2)/(1.0*ns_)<<endl;



    /*
    //-----------TEMP-----------------//
    complex<double> Local_Lz2_temp;
    int c1,c2,c3,c4;

    cout<<"XXXXXXXXX--Local Lz^2---XXXXXXXXXXXXXXXX"<<endl;
    int ORB1_=XY_;
    int ORB2_=YZ_;
    for(int site=0;site<ns_;site++){
        Local_Lz2_temp=zero_complex;

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                for(int n=0;n<Hamiltonian_.eigs_.size();n++){
                    for(int m=0;m<Hamiltonian_.eigs_.size();m++){

                        c1=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma);
                        c2=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma);
                        c3=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma_p);
                        c4=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma_p);

                        Local_Lz2_temp += (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    1.0*conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                    conj(Hamiltonian_.Ham_(c3,m))*Hamiltonian_.Ham_(c4,m)
                                    );

                        c1=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma);
                        c2=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma);
                        c3=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma_p);
                        c4=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma_p);

                        Local_Lz2_temp += (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    -1.0*conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                    conj(Hamiltonian_.Ham_(c3,m))*Hamiltonian_.Ham_(c4,m)
                                    );

                        c1=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma);
                        c2=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma);
                        c3=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma_p);
                        c4=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma_p);

                        Local_Lz2_temp += (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    +1.0*conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                    conj(Hamiltonian_.Ham_(c3,m))*Hamiltonian_.Ham_(c4,m)
                                    );

                        c1=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma);
                        c2=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma);
                        c3=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma_p);
                        c4=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma_p);

                        Local_Lz2_temp += (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                (
                                    -1.0*conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                    conj(Hamiltonian_.Ham_(c3,m))*Hamiltonian_.Ham_(c4,m)
                                    );


                        if(n!=m){
                            c1=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma);
                            c2=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma);
                            c3=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma_p);
                            c4=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma_p);

                            Local_Lz2_temp += (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                    (1.0 - (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0)))*
                                    (
                                        1.0*conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,m)*
                                        conj(Hamiltonian_.Ham_(c3,m))*Hamiltonian_.Ham_(c4,n)
                                        );

                            c1=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma);
                            c2=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma);
                            c3=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma_p);
                            c4=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma_p);

                            Local_Lz2_temp += (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                    (1.0 - (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0)))*
                                    (
                                        -1.0*conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,m)*
                                        conj(Hamiltonian_.Ham_(c3,m))*Hamiltonian_.Ham_(c4,n)
                                        );

                            c1=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma);
                            c2=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma);
                            c3=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma_p);
                            c4=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma_p);

                            Local_Lz2_temp += (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                    (1.0 - (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0)))*
                                    (
                                        1.0*conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,m)*
                                        conj(Hamiltonian_.Ham_(c3,m))*Hamiltonian_.Ham_(c4,n)
                                        );

                            c1=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma);
                            c2=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma);
                            c3=Coordinates_.Nc_dof(site, ORB2_ + 3*sigma_p);
                            c4=Coordinates_.Nc_dof(site, ORB1_ + 3*sigma_p);

                            Local_Lz2_temp += (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))*
                                    (1.0 - (1.0/( exp((Hamiltonian_.eigs_[m]-Parameters_.mus)*Parameters_.beta ) + 1.0)))*
                                    (
                                        -1.0*conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,m)*
                                        conj(Hamiltonian_.Ham_(c3,m))*Hamiltonian_.Ham_(c4,n)
                                        );

                        }


                    }

                }

            }

        }
        cout << site<< "    "<<Local_Lz2_temp<<endl;

    }
    */


}




void Observables::Calculate_Excitoncorrelations_Smartly(){


    //    0==yz
    //    1==xz
    //    2==xy
    //    0==up
    //    1==dn

    int YZ_,XZ_,XY_;
    YZ_=0; XZ_=1; XY_=2;

    string Delta_m1by2_q_out = "Delta_m1by2_q.txt";
    ofstream file_Delta_m1by2_q_out(Delta_m1by2_q_out.c_str());
    file_Delta_m1by2_q_out<<"#qx  qy   qx_index    qy_index   Delta_m1by2(qx,qy)"<<endl;

    string Delta_dag_m1by2_Delta_m1by2_r_out = "Delta_dag_m1by2_Delta_m1by2_r.txt";
    ofstream file_Delta_dag_m1by2_Delta_m1by2_r_out(Delta_dag_m1by2_Delta_m1by2_r_out.c_str());
    file_Delta_dag_m1by2_Delta_m1by2_r_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)  Delta_dag_m1by2_Delta_m1by2[site_i][site_j]"<<endl;


    int ci,cj;

    Mat_4_Complex_doub Fermi_, Inv_Fermi_ ;

    Fermi_.resize(ns_); Inv_Fermi_.resize(ns_);
    for(int n=0;n<ns_;n++){
        Fermi_[n].resize(ns_);
        Inv_Fermi_[n].resize(ns_);
        for(int m=0;m<ns_;m++){
            Fermi_[n][m].resize(6);
            Inv_Fermi_[n][m].resize(6);
            for(int state=0;state<6;state++){
                Fermi_[n][m][state].resize(6);
                Inv_Fermi_[n][m][state].resize(6);
            }
        }
    }

    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            for(int alpha=0;alpha<3;alpha++){
                for(int sigma=0;sigma<2;sigma++){
                    for(int alpha_p=0;alpha_p<3;alpha_p++){
                        for(int sigma_p=0;sigma_p<2;sigma_p++){

                            ci=Coordinates_.Nc_dof(i,alpha + 3*sigma);
                            cj=Coordinates_.Nc_dof(j,alpha_p + 3*sigma_p);

                            Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p]=zero_complex;
                            Inv_Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p]=zero_complex;
                            for(int n=0;n<Hamiltonian_.eigs_.size();n++){
                                Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p] += conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                                Inv_Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p] += conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                        (1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                                //  Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p] += conj(Hamiltonian_.Ham_(n,ci))*Hamiltonian_.Ham_(n,cj)*
                                //            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                                //  Inv_Fermi_[i][j][alpha + 3*sigma][alpha_p + 3*sigma_p] += conj(Hamiltonian_.Ham_(n,ci))*Hamiltonian_.Ham_(n,cj)*
                                //         (1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                            }

                        }
                    }
                }
            }
        }
    }



    Mat_2_Complex_doub Delta_dag_m1by2_Delta_m1by2_ri_rj;
    Delta_dag_m1by2_Delta_m1by2_ri_rj.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        Delta_dag_m1by2_Delta_m1by2_ri_rj[site_i].resize(ns_);
    }





    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            Delta_dag_m1by2_Delta_m1by2_ri_rj[i][j]=zero_complex;

            for(int alpha1=0;alpha1<3;alpha1++){
                for(int alpha2=0;alpha2<3;alpha2++){
                    for(int alpha3=0;alpha3<3;alpha3++){
                        for(int alpha4=0;alpha4<3;alpha4++){

                            for(int sigma1=0;sigma1<2;sigma1++){
                                for(int sigma2=0;sigma2<2;sigma2++){
                                    for(int sigma3=0;sigma3<2;sigma3++){
                                        for(int sigma4=0;sigma4<2;sigma4++){

                                            if(Transformation(5,alpha1+3*sigma1) != zero_complex
                                                    &&
                                                    Transformation(3,alpha2+3*sigma2) != zero_complex
                                                    &&
                                                    Transformation(3,alpha3+3*sigma3) != zero_complex
                                                    &&
                                                    Transformation(5,alpha4+3*sigma4) != zero_complex
                                                    ){
                                                Delta_dag_m1by2_Delta_m1by2_ri_rj[i][j] += conj(Transformation(5,alpha1+3*sigma1))*Transformation(3,alpha2+3*sigma2)*
                                                        conj(Transformation(3,alpha3+3*sigma3))*Transformation(5,alpha4+3*sigma4)*
                                                        ((Fermi_[i][i][alpha1+3*sigma1][alpha2+3*sigma2]*
                                                        Fermi_[j][j][alpha3+3*sigma3][alpha4+3*sigma4])
                                                        +
                                                        (Fermi_[i][j][alpha1+3*sigma1][alpha4+3*sigma4]*
                                                        Inv_Fermi_[j][i][alpha3+3*sigma3][alpha2+3*sigma2]));

                                            }




                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }

        }

    }


    int site_i, site_j;

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            site_i=Coordinates_.Nc(ix,iy);

            for(int jx=0;jx<lx_;jx++){
                for(int jy=0;jy<ly_;jy++){
                    site_j=Coordinates_.Nc(jx,jy);
                    file_Delta_dag_m1by2_Delta_m1by2_r_out<<site_i<<"\t"<<ix<<"\t"<<iy<<"\t"<<site_j<<"\t"<<
                                                            jx<<"\t"<<jy<<"\t"<<
                                                            real(Delta_dag_m1by2_Delta_m1by2_ri_rj[site_i][site_j])<<"\t"<<
                                                            imag(Delta_dag_m1by2_Delta_m1by2_ri_rj[site_i][site_j])<<endl;


                }
            }

        }
    }


    double qx_, qy_;
    complex<double> value;


    for(int qy=0;qy<ly_;qy++){
        for(int qx=0;qx<lx_;qx++){
            value = zero_complex;
            qx_ = (2.0*qx*PI)*(1.0/(1.0*lx_));
            qy_ = (2.0*qy*PI)*(1.0/(1.0*ly_));

            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){
                    site_i=Coordinates_.Nc(ix,iy);

                    for(int jx=0;jx<lx_;jx++){
                        for(int jy=0;jy<ly_;jy++){
                            site_j=Coordinates_.Nc(jx,jy);

                            value += one_complex*(exp(iota_complex*( (qx_*(ix - jx)) + (qy_*(iy-jy)) )))*
                                    Delta_dag_m1by2_Delta_m1by2_ri_rj[site_i][site_j]*(1.0/(1.0*ns_));

                        }}
                }}

            file_Delta_m1by2_q_out<<qx_<<"\t"<<qy_<<"\t"<<qx<<"\t"<<qy<<"\t"<<real(value)<<"\t"<<imag(value)<<endl;

        }
        file_Delta_m1by2_q_out<<endl;

    }


}



void Observables::Get_Non_Interacting_dispersion(){

}


double Observables::Lorentzian(double x, double brd){
    double temp;

    temp = (1.0/PI)*( (brd/2.0)/ ( (x*x) + ((brd*brd)/4.0) ) );

    return temp;

}

void Observables::DensityOfStates(){
    //-----------Calculate Bandwidth------//
    BandWidth=2.0;
    //-----------------------------------//

} // ----------


void Observables::OccDensity(){

} // ----------


void Observables::TotalOccDensity(){

} // ----------



void Observables::Total_Energy_Average(double Curr_QuantE, double CurrE){

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE)*(Curr_QuantE + CurrE);
}

void Observables::OccDensity(int tlabel){

} // ----------



void Observables::Initialize(){


    F_n.resize(ns_*36); //F_n=x_n_out - x_n_in []
    F_nm1.resize(ns_*36);
    DeltaF_n.resize(ns_*36); //DeltaF_n=F_n - F-nm1;
    Delta_x_n.resize(ns_*36); //Delta_x_n= x_n_in - x_nm1_in;
    Jinv_n.resize(ns_*36);
    Jinv_np1.resize(ns_*36);


    _Fm.resize((ns_*36));
    _Delta_OPm.resize((ns_*36));
    _Fm_minus1.resize((ns_*36));
    _Delta_OPm_minus1.resize((ns_*36));




    for(int i=0;i<ns_*36;i++){
        Jinv_n[i]. resize(ns_*36);
        Jinv_np1[i]. resize(ns_*36);
    }


    Transformation.resize(6,6);
    //Saves a_{jm} ---to---> c_{spin,orbital} [NOT dagger]
    /*
    Transformation(0,1)=iota_complex*(-1.0/sqrt(2));
    Transformation(0,2)=one_complex*(1.0/sqrt(6));
    Transformation(0,4)=one_complex*(-1.0/sqrt(3));

    Transformation(1,1)=one_complex*(1.0/sqrt(2));
    Transformation(1,2)=iota_complex*(-1.0/sqrt(6));
    Transformation(1,4)=iota_complex*(1.0/sqrt(3));

    Transformation(2,3)=one_complex*(2.0/sqrt(6));
    Transformation(2,5)=one_complex*(1.0/sqrt(3));

    Transformation(3,0)=iota_complex*(1.0/sqrt(2));
    Transformation(3,3)=one_complex*(-1.0/sqrt(6));
    Transformation(3,5)=one_complex*(1.0/sqrt(3));

    Transformation(4,0)=one_complex*(1.0/sqrt(2));
    Transformation(4,3)=iota_complex*(-1.0/sqrt(6));
    Transformation(4,5)=iota_complex*(1.0/sqrt(3));

    Transformation(5,2)=one_complex*(2.0/sqrt(6));
    Transformation(5,4)=one_complex*(1.0/sqrt(3));
    */

    Transformation(1,0)=iota_complex*(1.0/sqrt(2));
    Transformation(2,0)=one_complex*(1.0/sqrt(6));
    Transformation(4,0)=one_complex*(-1.0/sqrt(3));

    Transformation(1,1)=one_complex*(1.0/sqrt(2));
    Transformation(2,1)=iota_complex*(1.0/sqrt(6));
    Transformation(4,1)=iota_complex*(-1.0/sqrt(3));

    Transformation(3,2)=one_complex*(2.0/sqrt(6));
    Transformation(5,2)=one_complex*(1.0/sqrt(3));

    Transformation(0,3)=iota_complex*(-1.0/sqrt(2));
    Transformation(3,3)=one_complex*(-1.0/sqrt(6));
    Transformation(5,3)=one_complex*(1.0/sqrt(3));

    Transformation(0,4)=one_complex*(1.0/sqrt(2));
    Transformation(3,4)=iota_complex*(1.0/sqrt(6));
    Transformation(5,4)=iota_complex*(-1.0/sqrt(3));

    Transformation(2,5)=one_complex*(2.0/sqrt(6));
    Transformation(4,5)=one_complex*(1.0/sqrt(3));


    F_Exciton.resize(ns_);
    for(int i =0;i<ns_;i++){
        F_Exciton[i].resize(ns_);
        for(int j=0;j<ns_;j++){
            F_Exciton[i][j].resize(2);
            for(int type=0;type<2;type++){
                F_Exciton[i][j][type].resize(2);
            }
        }
    }



} // ----------


double Observables::Omega(int i){
    return -20.0+double(i)*dosincr_;
} // ----------









#endif // OBSERVABLES_H
