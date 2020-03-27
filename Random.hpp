//
//  Random.hpp
//
//  Created by Raffaele Marino on 2018-08-25.
//  Copyright © 2018 Raffaele Marino. All rights reserved.
//

/*Class Random Using WELL Algorithm, minimal modifiications for compatibility*/

/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

#ifndef Random_h
#define Random_h
#include "Header.h"
#define W 32
#define R 32
#define M1 3
#define M2 24
#define M3 10

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define Identity(v) (v)

#define FACT 2.32830643653869628906e-10

static unsigned int state_i = 0;
static unsigned int STATE[R];
static unsigned int z0, z1, z2;


#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000001fU]
#define VM2           STATE[(state_i+M2) & 0x0000001fU]
#define VM3           STATE[(state_i+M3) & 0x0000001fU]
#define VRm1          STATE[(state_i+31) & 0x0000001fU]
#define newV0         STATE[(state_i+31) & 0x0000001fU]
#define newV1         STATE[state_i                   ]


#define INITLEN 32


class Random{
    
    
public:
    
    /*declaration public members for the Random.hpp*/
    
    Random(){
        WellRandomInitialization();
    };
    
    ~Random(){};
    
    void print_seed();
    void print_seeds_on_file();
    double random_number();
    unsigned int seed();
    void GetNSITE(const long int NSite);
    
private:
    
    /*declaration private variables for the Random.hpp*/
    unsigned int _seed;
    long int _NSITES;
    vector<unsigned int> _vecinittoprint;
    
    void InitWELLRNG1024a (unsigned int *init);
    
    double WELLRNG1024a (void);
    
    unsigned int GetRandom();
    
    void WellRandomInitialization();
    
    
};

inline void Random::GetNSITE(const long int NSite){
    _NSITES=NSite;
}

/*random number generator*/
inline double Random::random_number(){return WELLRNG1024a();}

/*print seed*/
inline void Random::print_seed(){cout<<_seed<<endl;}

inline unsigned int Random::seed(){return _seed;}
/*initialization random number generator*/
inline void Random::InitWELLRNG1024a (unsigned int *init){
    int j;
    state_i = 0;
    for (j = 0; j < R; j++)
        STATE[j] = init[j];
}
/*generator random number*/
inline double Random::WELLRNG1024a (void){
    z0    = VRm1;
    z1    = Identity(V0)       ^ MAT0POS (8, VM1);
    z2    = MAT0NEG (-19, VM2) ^ MAT0NEG(-14,VM3);
    newV1 = z1                 ^ z2;
    newV0 = MAT0NEG (-11,z0)   ^ MAT0NEG(-7,z1)    ^ MAT0NEG(-13,z2) ;
    state_i = (state_i + 31) & 0x0000001fU;
    return ((double) STATE[state_i]  * FACT);
}

/*generator random seed*/
inline unsigned int Random::GetRandom(){/*inizialzation seed*/
    unsigned int ris;
    
    FILE *devran = fopen("/dev/urandom","r");
    
    fread(&ris, 4, 1, devran);
    
    fclose(devran);
    
    return ris;
}

/*initialization random generator*/
inline void Random::WellRandomInitialization(){
    _seed=GetRandom();
    if (_seed < 1) exit(-1);
    srandom(_seed);
    unsigned int init[INITLEN];
    double j=0;
    for (int i=0;i<INITLEN;++i) {
        init[i] = (unsigned int) random();
        _vecinittoprint.push_back(init[i]);
    }
    InitWELLRNG1024a(init);
    for (int i=0;i<20;++i) j = WELLRNG1024a();
}

/*public member for printing seeds on file*/
void Random::print_seeds_on_file(){
    
    string mystringseedinit; //Name of the file where data will be written
    ostringstream seedsite_str;
    ostringstream nsite_str;
    seedsite_str<<_seed;
    nsite_str<<_NSITES;
    
    const string titleseed="Seed_init_N=";
    const string seedstring="seed=";
    const string dat=".txt";
    mystringseedinit=directory;
    mystringseedinit+=titleseed;
    mystringseedinit+=nsite_str.str();
    mystringseedinit+=seedstring;
    mystringseedinit+=seedsite_str.str();
    mystringseedinit+=dat;
    
    ofstream outfilewriteseed(mystringseedinit.c_str(), ios_base::app );
    for (int ci=0; ci<INITLEN; ++ci) {
        outfilewriteseed<<_vecinittoprint[ci]<<endl;
    }
}

#endif /* Random_h */
