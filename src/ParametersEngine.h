#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{

public:
    int lx, ly, ns, IterMax, RandomSeed;
    double Convergence_Error;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus,Total_Particles,pi, mu_old;
    double MuValueFixed;
    bool FixingMu;
    double J_Hund;
    double U_onsite;
    double U_prime_onsite;
    double Lambda_SOC;
    double Disorder_Strength, RandomDisorderSeed;
    bool PBC_X, PBC_Y;
    bool PNICTIDES_HOPPING;

    double dw_dos, eta_dos;
    double w_min, w_max;

    bool Read_OPs;
    bool Create_OPs;
    bool Just_Hartree;
    string File_OPs_in, File_OPs_out;
    string Create_OPs_type;


    Matrix<double> t2g_hopping_NN_X;
    Matrix<double> t2g_hopping_NN_Y;
    Matrix<double> t2g_hopping_NNN_PXPY;
    Matrix<double> t2g_hopping_NNN_PXMY;

    Mat_1_doub Crystal_Field;

    bool Simple_Mixing;
    bool Broyden_Mixing;
    bool BroydenSecondMethodMixing;
    double w_minus1,wn;
    int BroydenSecondMethodCounter;
    double alpha_OP;

    double Temperature,beta,Eav,maxmoment;

    bool ReadDisorder;
    string ReadDisorderString, FixingMuString;
    string DisorderSeedFile;


    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){



    double Simple_Mixing_double, Broyden_Mixing_double, BroydenSecondMethodMixing_double;
    double Read_OPs_double;
    double Create_OPs_double;
    double Just_Hartree_double;
    string PBC_X_string, PBC_Y_string, Pnictides_Hopping_string;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;



    lx = int(matchstring(inputfile_,"Xsite"));
    ly = int(matchstring(inputfile_,"Ysite"));
    TBC_mx = int(matchstring(inputfile_,"TwistedBoundaryCond_mx"));
    TBC_my = int(matchstring(inputfile_,"TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_,"TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_,"TBC_cellsY"));
    BroydenSecondMethodCounter = int(matchstring(inputfile_,"BroydenSecondMethodCounter"));

    ns = lx*ly;
    cout << "TotalNumberOfSites = "<< ns << endl;

    Total_Particles = matchstring(inputfile_,"Total No. of particles");
    cout << "TotalNumberOfParticles = "<< Total_Particles << endl;

    IterMax = int(matchstring(inputfile_,"No_of_SelfConsistency_iters"));
    Convergence_Error=matchstring(inputfile_,"Convergence_Error");
    RandomSeed = matchstring(inputfile_,"RandomSeed");
    FixingMuString = matchstring2(inputfile_, "FixingMu");
    RandomDisorderSeed = matchstring(inputfile_,"RandomDisorderSeed");
    ReadDisorderString = matchstring2(inputfile_,"ReadDisorderConf");
    Disorder_Strength = matchstring(inputfile_,"Disorder_Strength");
    J_Hund = matchstring(inputfile_,"J_HUND");
    U_onsite = matchstring(inputfile_,"U_Onsite");
    U_prime_onsite = matchstring(inputfile_,"U_prime_Onsite");
    Lambda_SOC = matchstring(inputfile_, "Lambda_SOC");
    Temperature = matchstring(inputfile_,"Temperature");

    //dw_dos, eta_dos
    dw_dos = matchstring(inputfile_, "dw_dos");
    eta_dos = matchstring(inputfile_, "eta_dos");
    w_min = matchstring(inputfile_, "w_min");
    w_max = matchstring(inputfile_, "w_max");

    alpha_OP = matchstring(inputfile_,"alpha_OP");
    w_minus1 = matchstring(inputfile_,"w_minus1");
    wn = matchstring(inputfile_,"wn");



    Dflag = 'N';

    Simple_Mixing_double=double(matchstring(inputfile_,"Simple_Mixing"));
    Broyden_Mixing_double=double(matchstring(inputfile_,"Broyden_Mixing"));
    BroydenSecondMethodMixing_double=double(matchstring(inputfile_,"Broyden_Second_Method_Mixing"));



    if(FixingMuString=="true"){
        FixingMu=true;
        MuValueFixed = matchstring(inputfile_,"MuValueFixed");
    }
    else{
        FixingMu=false;
        MuValueFixed=-100000;
    }

    if(ReadDisorderString=="true"){
        ReadDisorder = true;
        DisorderSeedFile = matchstring2(inputfile_,"DisorderSeedFile");
    }
    else{
        ReadDisorder = false;
    }

    if(BroydenSecondMethodMixing_double==1.0){
         BroydenSecondMethodMixing=true;
        Broyden_Mixing=false;
        Simple_Mixing=false;
    }
    else{
         BroydenSecondMethodMixing=false;
        if(Broyden_Mixing_double==1.0){
            Broyden_Mixing=true;
            Simple_Mixing=false;

        }
        else if(Broyden_Mixing_double==0.0){
            Broyden_Mixing=false;
            Simple_Mixing=true;
            cout<<"Broyden_Mixing and  BroydenSecondMethodMixing, both are 0(false). So Simple mixing is used"<<endl;

        }

    }



    Read_OPs_double=double(matchstring(inputfile_,"Read_initial_OPvalues"));
    if(Read_OPs_double==1.0){
        Read_OPs=true;
    }
    else{
        Read_OPs=false;
    }




    Just_Hartree_double=double(matchstring(inputfile_,"Just_Hartree"));
    if(Just_Hartree_double==1.0){
        Just_Hartree=true;
    }
    else{
        Just_Hartree=false;
    }

    Create_OPs_double=double(matchstring(inputfile_,"Create_OPvalues"));
    if(Create_OPs_double==1.0){
        assert(!Read_OPs);
        Create_OPs=true;
    Create_OPs_type=matchstring2(inputfile_,"Create_OPType");
    }
    else{
        Create_OPs=false;
    }



    PBC_X_string=matchstring2(inputfile_,"PBC_X");
    if(PBC_X_string=="true"){
        PBC_X=true;
    }
    else{
        PBC_X=false;
    }
    PBC_Y_string=matchstring2(inputfile_,"PBC_Y");
    if(PBC_Y_string=="true"){
        PBC_Y=true;
    }
    else{
        PBC_Y=false;
    }

    Pnictides_Hopping_string=matchstring2(inputfile_,"Pnictides_Hopping");
    if(Pnictides_Hopping_string=="true"){
        PNICTIDES_HOPPING=true;
    }
    else{
        PNICTIDES_HOPPING=false;
    }

    File_OPs_in=matchstring2(inputfile_,"Read_initial_OPvalues_file");
    File_OPs_out=matchstring2(inputfile_,"Write_Final_OPvalues_file");

    string Nearest_Neigh_Hopping_t2g_basis_row0_X;
    string Nearest_Neigh_Hopping_t2g_basis_row1_X;
    string Nearest_Neigh_Hopping_t2g_basis_row2_X;
    string Nearest_Neigh_Hopping_t2g_basis_row0_Y;
    string Nearest_Neigh_Hopping_t2g_basis_row1_Y;
    string Nearest_Neigh_Hopping_t2g_basis_row2_Y;


    Nearest_Neigh_Hopping_t2g_basis_row0_X=matchstring2(inputfile_, "Nearest_Neigh_Hopping_X_t2g_basis_row0");
    Nearest_Neigh_Hopping_t2g_basis_row1_X=matchstring2(inputfile_, "Nearest_Neigh_Hopping_X_t2g_basis_row1");
    Nearest_Neigh_Hopping_t2g_basis_row2_X=matchstring2(inputfile_, "Nearest_Neigh_Hopping_X_t2g_basis_row2");
    Nearest_Neigh_Hopping_t2g_basis_row0_Y=matchstring2(inputfile_, "Nearest_Neigh_Hopping_Y_t2g_basis_row0");
    Nearest_Neigh_Hopping_t2g_basis_row1_Y=matchstring2(inputfile_, "Nearest_Neigh_Hopping_Y_t2g_basis_row1");
    Nearest_Neigh_Hopping_t2g_basis_row2_Y=matchstring2(inputfile_, "Nearest_Neigh_Hopping_Y_t2g_basis_row2");

    stringstream t2g_row0_stream_X(Nearest_Neigh_Hopping_t2g_basis_row0_X);
    stringstream t2g_row1_stream_X(Nearest_Neigh_Hopping_t2g_basis_row1_X);
    stringstream t2g_row2_stream_X(Nearest_Neigh_Hopping_t2g_basis_row2_X);
    stringstream t2g_row0_stream_Y(Nearest_Neigh_Hopping_t2g_basis_row0_Y);
    stringstream t2g_row1_stream_Y(Nearest_Neigh_Hopping_t2g_basis_row1_Y);
    stringstream t2g_row2_stream_Y(Nearest_Neigh_Hopping_t2g_basis_row2_Y);

    t2g_hopping_NN_X.resize(3,3);
    for(int n=0;n<3;n++){
        t2g_row0_stream_X >> t2g_hopping_NN_X(0,n);
        t2g_row1_stream_X >> t2g_hopping_NN_X(1,n);
        t2g_row2_stream_X >> t2g_hopping_NN_X(2,n);
    }

    t2g_hopping_NN_Y.resize(3,3);
    for(int n=0;n<3;n++){
        t2g_row0_stream_Y >> t2g_hopping_NN_Y(0,n);
        t2g_row1_stream_Y >> t2g_hopping_NN_Y(1,n);
        t2g_row2_stream_Y >> t2g_hopping_NN_Y(2,n);
    }



    //Next Nearest hopping------------

    string NNN_t2g_basis_row0_PXPY;
    string NNN_t2g_basis_row1_PXPY;
    string NNN_t2g_basis_row2_PXPY;
    string NNN_t2g_basis_row0_PXMY;
    string NNN_t2g_basis_row1_PXMY;
    string NNN_t2g_basis_row2_PXMY;


    NNN_t2g_basis_row0_PXPY=matchstring2(inputfile_, "NNN_Hopping_PXPY_t2g_basis_row0");
    NNN_t2g_basis_row1_PXPY=matchstring2(inputfile_, "NNN_Hopping_PXPY_t2g_basis_row1");
    NNN_t2g_basis_row2_PXPY=matchstring2(inputfile_, "NNN_Hopping_PXPY_t2g_basis_row2");
    NNN_t2g_basis_row0_PXMY=matchstring2(inputfile_, "NNN_Hopping_PXMY_t2g_basis_row0");
    NNN_t2g_basis_row1_PXMY=matchstring2(inputfile_, "NNN_Hopping_PXMY_t2g_basis_row1");
    NNN_t2g_basis_row2_PXMY=matchstring2(inputfile_, "NNN_Hopping_PXMY_t2g_basis_row2");

    stringstream t2g_row0_stream_PXPY(NNN_t2g_basis_row0_PXPY);
    stringstream t2g_row1_stream_PXPY(NNN_t2g_basis_row1_PXPY);
    stringstream t2g_row2_stream_PXPY(NNN_t2g_basis_row2_PXPY);
    stringstream t2g_row0_stream_PXMY(NNN_t2g_basis_row0_PXMY);
    stringstream t2g_row1_stream_PXMY(NNN_t2g_basis_row1_PXMY);
    stringstream t2g_row2_stream_PXMY(NNN_t2g_basis_row2_PXMY);

    t2g_hopping_NNN_PXPY.resize(3,3);
    for(int n=0;n<3;n++){
        t2g_row0_stream_PXPY >> t2g_hopping_NNN_PXPY(0,n);
        t2g_row1_stream_PXPY >> t2g_hopping_NNN_PXPY(1,n);
        t2g_row2_stream_PXPY >> t2g_hopping_NNN_PXPY(2,n);
    }

    t2g_hopping_NNN_PXMY.resize(3,3);
    for(int n=0;n<3;n++){
        t2g_row0_stream_PXMY >> t2g_hopping_NNN_PXMY(0,n);
        t2g_row1_stream_PXMY >> t2g_hopping_NNN_PXMY(1,n);
        t2g_row2_stream_PXMY >> t2g_hopping_NNN_PXMY(2,n);
    }


    //----------------------------------

    string Crystal_Field_t2g;
    Crystal_Field_t2g=matchstring2(inputfile_, "Crystal_Field_t2g");

    stringstream Crystal_Field_t2g_stream(Crystal_Field_t2g);


    Crystal_Field.resize(3);
    for(int n=0;n<3;n++){
        Crystal_Field_t2g_stream >> Crystal_Field[n];
    }

    pi=4.00*atan(double(1.0));
    Eav=0.0;


    //beta=(11605.0/Temperature);
     beta=(1.0/Temperature);

    mus=0.25;
    mu_old=0.25;
    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file,string match) {

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if(readFile.is_open())
    {
        while(!readFile.eof())
        {
            getline(readFile,line);

            if ((offset = line.find(match, 0)) != string::npos) {
                amount = line.substr (offset+match.length()+1);				}

        }
        readFile.close();
    }
    else
    {cout<<"Unable to open input file while in the Parameters class."<<endl;}




    cout << match << " = " << amount << endl;
    return amount;
}

#endif



