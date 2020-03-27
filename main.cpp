//
//  main.cpp
// 
//
//  Created by Raffaele Marino on 2018-08-25.
//  Copyright © 2018 Raffaele Marino. All rights reserved.
//

/*
  This code describes the greedy algorithm presented in the paper "Large independent set for random d-regular graphs with small d"
 For any bug present, please send an email to:
 raffaele.marino@epfl.ch
 marino.raffaele@mail.huji.ac.il
 marinoraffaele.nunziatella@gmail.com
 */

#include "Header.h"
#include "Random.hpp"
unsigned int D=0;

/*Class graph: it contains all information about the graph*/
class Graph{
public:
    
    unsigned int d=D;
    /*constructor*/
    Graph(){};
    /*get values for the graph*/
    void Get_Values_Graph(long unsigned int &N, int &_d){
        order_graph=N;
    };
    /*order of the graph*/
    long unsigned int _N(){return order_graph;};
    /*connectivity of each single node*/
    int _d(){return d;};
    /*print the order of the graph*/
    void Order_Graph(){
        cout<<"Order Graph is: "<<order_graph<<endl;
    };
    /*print the size of the graph*/
    void Size_Graph(){
        cout<<"Size Graph is: "<<d*order_graph<<endl;
    };
    /*print the connectivity of the graph*/
    void Connectivity_Graph(){
        cout<<"Connectivity Graph is: "<<d<<endl;
    };
    

private:
    long unsigned int order_graph;
    
};


/*Class Vertex: it is a derived class from Graph class. This class contains all information about the node */
class Vertex:public Graph{
public:
    /*constructor*/
    Vertex(){
        _refdeg=0;
        _ref_char_val=' ';
        _who_am_I=' ';
        _I_am_a_two=false;
        _I_am_contracted=false;
        _I_am_the_contractor=false;
        _I_am_a_possible_2=false;
        My_contractor_is=NULL;
        _I_am_a_three=false;
        _I_am_in_vec_list=false;
    };
    
    /*print the index of the vertex*/
    void vertex_i(){
        cout<<"Vertex labeled: "<<_vertex_i<<endl;
    };
    
    /*operator vertex class*/
    Vertex* operator = (Vertex* V){
        _vertex_i=V->_vertex_i;
        _refdeg=V->_refdeg;
        _ref_char_val=V->_ref_char_val;
        _list_ngb=V->_list_ngb;
        _who_am_I=V->_who_am_I;
        it_list2=V->it_list2;
        it_poss2=V->it_poss2;
        _I_am_a_possible_2=V->_I_am_a_possible_2;
        _I_am_a_three=V->_I_am_a_three;
        _I_am_contracted=V->_I_am_contracted;
        _I_am_the_contractor=V->_I_am_the_contractor;
        My_contractor_is=V->My_contractor_is;
        _I_am_in_vec_list=V-> _I_am_in_vec_list;
        it_vec_list=V->it_vec_list;
        return this;
    };
    
    /*operator for printing the object*/
    friend ostream &operator<<(ostream &out, Vertex *V){
        out<<V->_vertex_i<<" "<<V->_refdeg<<" "<<V->_ref_char_val<<" "<<V->_who_am_I<<" "<<V->_I_am_the_contractor<<" "<<V->_I_am_contracted<<endl;
        return out;
    }
    
    /*operator for printing the object*/
    friend ostream &operator<<(ostream &out, Vertex V){
        out<<V._vertex_i<<" "<<V._refdeg<<" "<<V._ref_char_val<<" "<<V._who_am_I<<" "<<V._I_am_the_contractor<<" "<<V._I_am_contracted<<endl;
        return out;
    }
    
    /*operator [] vertex object*/
    long long int operator[](long unsigned int i){
        return _list_ngb[i];
    }
    
    vector<long long int> _list_ngb; //vector of ngbs of i
    list<Vertex>::iterator it_list2; //iterator for list of P sites
    list<Vertex *>::iterator it_list3;// iterator for list of three
    list<Vertex *>::iterator it_poss2;
    list<Vertex *>::iterator it_vec_list; // iterator vector list
    Vertex * My_contractor_is; // pointer to a P site
    long unsigned int _vertex_i; //label vertex i
    int _refdeg; // degree of a node
    char _ref_char_val; // label of the node
    char _who_am_I; // final label
    bool _I_am_a_two;
    bool _I_am_contracted;
    bool _I_am_the_contractor;
    bool _I_am_a_possible_2;
    bool _I_am_a_three;
    bool _I_am_in_vec_list;
   
};

/*This class allows to merge together different nodes and create a virtual node*/
/* With contractor we identify on P sites that becomes the centre of the object to which the other nodes are referring to.*/
class Contraction:public Graph{
    public:
    // constructor
    Contraction(){
        _connectivity=0;
        _i=0;
        _P=0;
        _C=0;
        _I_am_in=false;
        _I_am_in_the_list=false;
        _I_have_a_loop_inside=false;
    };
    
    Vertex Give_Me_Contractor(){
        return *_Vref;
    };
    
    Vertex* Give_Me_Contractor_ptr(){
        return _Vref;
    };
    
    Vertex* Give_Me_Contracts_ptr(int &i){
        return _fill_vec[i];
    };

    void Get_Contractor(Vertex  &V1){
        _Vref= &V1;
    };
    
    Vertex Give_Me_Vertex(){
        return *_fill_vec[--_i];
    };
    
    Vertex* Give_Me_Vertex(long unsigned int &i){
        return _fill_vec[i];
    };
    
    void check(){

        for (unsigned int i=0; i<_fill_vec.size(); ++i) {
            if(_fill_vec[i]->My_contractor_is->_vertex_i!=_Vref->My_contractor_is->_vertex_i){
                cout<<"ERROR"<<endl;
                cout<<_fill_vec[i]->My_contractor_is->_vertex_i<<" "<<_Vref->My_contractor_is->_vertex_i<<endl;
                exit(-1);
            }
        }
    }
    
    void Clean(){
        _fill_vec.clear();
        _i=0;
        _connectivity=0;
        _P=0;
        _C=0;
        _I_am_in=false;
        _I_am_in_the_list=false;
        _I_have_a_loop_inside=false;
    }
    
    void Contract(Vertex &V2){
        _fill_vec.push_back(&V2);
        _i=_fill_vec.size();
    };
    
    void update_connectivity(vector<list<Contraction *> > &listsofcontractsites){
        if(_I_am_in_the_list){
            listsofcontractsites[_connectivity].erase(ptr_vec_list);
        }
        _connectivity=0;
        for (int i=0; i<_fill_vec.size(); ++i) {
            _connectivity+=d-((*_fill_vec[i])._refdeg);
        }
        listsofcontractsites[_connectivity].push_front(this);
        ptr_vec_list=listsofcontractsites[_connectivity].begin();
        _I_am_in_the_list=true;
    }
    
    long unsigned int _fill_vec_size(){
        return _fill_vec.size();
    };
    
    void _insert(long unsigned int &_j, Vertex  &V2){
        _fill_vec[_j]=&V2;
    }
    
    void swap_P_C(){
        _P=0;
        _C=0;
        for (int i=0; i<_fill_vec_size(); ++i) {
            if(_fill_vec[i]->_ref_char_val=='P'){
                _fill_vec[i]->_ref_char_val='C';
                ++_C;
            }else{
                _fill_vec[i]->_ref_char_val='P';
                ++_P;
            }
        }
        if(_Vref->_ref_char_val=='P'){
            _Vref->_ref_char_val='C';
            ++_C;
        }else{
            _Vref->_ref_char_val='P';
            ++_P;
        }
    }
    
    long unsigned int Find_V_to_fill(long unsigned int &N){
        Vertex * V=NULL;
        
        for (int i=0; i<_fill_vec_size(); ++i) {
            if(_fill_vec[i]->_ref_char_val=='P' && _fill_vec[i]->_refdeg!=d){
                V=_fill_vec[i];
                break;
            }
        }
        if(V==NULL){
            cout<<"WARNING 2"<<endl;
            return N+1;
//            exit(-1);
        }
        return V->_vertex_i;
    }
    
    
    list<long unsigned int> Find_V_to_fill_vec(){
        list<long unsigned int> v;
        for (int i=0; i<_fill_vec_size(); ++i) {
            if(_fill_vec[i]->_ref_char_val=='P' && _fill_vec[i]->_refdeg!=d){
                v.push_back(_fill_vec[i]->_vertex_i);
            }
        }
        if(v.empty()){
            cout<<"WARNING 3"<<endl;
//            exit(-1);
        }
        return v;
    }
    
    void erase(long unsigned int J){
        for (long unsigned int i=0; i<_fill_vec.size(); ++i) {
            if(_fill_vec[i]->_vertex_i==J){
                swap(_fill_vec[i], _fill_vec[_fill_vec.size()-1]);
                break;
            }
        }
        
        if(!_fill_vec.empty()){
            _fill_vec.pop_back();
            _i=_fill_vec.size();
        }
        else _i=0;
    }
    
    void fix_IS_VC(){
        for (int i=0; i<_fill_vec.size(); ++i) {
            if(_fill_vec[i]->_ref_char_val=='P'){
                _fill_vec[i]->_who_am_I='I';
            }else{
                _fill_vec[i]->_who_am_I='V';
                if(_fill_vec[i]->_refdeg!=d)_VC_to_fill.push_front(_fill_vec[i]);
            }
        }
        if(_Vref->_ref_char_val=='P'){
            _Vref->_who_am_I='I';
        }else{
             _Vref->_who_am_I='V';
            if(_Vref->_refdeg!=d)_VC_to_fill.push_front(_Vref);
        }
    }
    
    
    void find_sites_to_fill(){
        for (int i=0; i<_fill_vec_size(); ++i) {
            if(_fill_vec[i]->_refdeg!=d)_VC_to_fill.push_front(_fill_vec[i]);
        }
            if(_Vref->_refdeg!=d)_VC_to_fill.push_front(_Vref);
        
    }
    
    
    Contraction* _fusion_c(Contraction &C){
        for (int i=0;i<C._fill_vec_size(); ++i) {
            _fill_vec.push_back(C._fill_vec[i]);
            C._fill_vec[i]->My_contractor_is=Give_Me_Contractor_ptr();
        }
        _fill_vec.push_back(C.Give_Me_Contractor_ptr());
        C.Give_Me_Contractor_ptr()->My_contractor_is=Give_Me_Contractor_ptr();
        _i=_fill_vec.size();
        return this;
    };
    
    friend ostream &operator<<(ostream &out, Contraction *V){
        out<<V->_Vref;
        for (int i=0; i<V->_i; ++i) {
            out<<V->_fill_vec[i];
            
        }
        return out;
    }
    
    friend ostream &operator<<(ostream &out, Contraction V){
        out<<"/*******************************/"<<endl;
        out<<V._Vref->_vertex_i<<" "<<V._Vref->_who_am_I<<" "<<V._Vref->_ref_char_val<<endl;
        out<<"/*******************************/"<<endl;
        for (int i=0; i<V._i; ++i) {
            out<<V._fill_vec[i]->_vertex_i<<" "<<V._fill_vec[i]->_who_am_I<<" "<<V._fill_vec[i]->_ref_char_val<<endl;
        }
        return out;
    }
    
    bool operator !=(Contraction &C){
        if((_Vref->_vertex_i)==(C._Vref->_vertex_i))return false;
        else return true;
    };
    
    Contraction & operator = (Contraction &C){
        _Vref=C._Vref;
        _fill_vec=C._fill_vec;
        _P=C._P;
        _C=C._C;
        _i=C._i;
        _I_am_in_the_list=C._I_am_in_the_list;
        ptr_vec_list=C.ptr_vec_list;
        _connectivity=C._connectivity;
        _I_have_a_loop_inside=C._I_have_a_loop_inside;
        return *this;
    };
    
    Vertex operator [](int &i){
        return *_fill_vec[i];
    };
    
    bool _I_am_in;
    bool _I_am_in_the_list;
    bool _I_have_a_loop_inside;
    long unsigned int _connectivity;
    list<Contraction *> ::iterator ptr_vec_list;
    long unsigned int _i;
    long unsigned int _P, _C;
    list<Vertex *> _VC_to_fill;
private:
    Vertex *_Vref;
    vector<Vertex *> _fill_vec;
};

// end declaration and definition of classes

// declaration and definition of functions used for splitting and storing information of the graph.
void _myswap(long unsigned int &i, int &d, vector<vector<long unsigned int> > &BACKLIST, vector<long unsigned int> &IJLINKS, long unsigned int &I, long unsigned int &J, long unsigned int &iflag){
    long unsigned int _j=IJLINKS.size();
    int k1=0, k2=0;
    for(k1=0;k1<d;++k1) {if (BACKLIST[J][k1]>=i) break;}
    _j= BACKLIST[J][k1];
    swap(IJLINKS[i], IJLINKS[_j]);
    BACKLIST[J][k1]=i;
    for (k2=0;k2<d;++k2) {if (BACKLIST[I][k2]==i) break;}
    BACKLIST[I][k2] = _j;
    I=IJLINKS[i];
    }

void _myswap(long unsigned int &i, int &d, vector<vector<long unsigned int> > &BACKLIST, vector<long unsigned int> &IJLINKS, long unsigned int &I, long unsigned int &J, long unsigned int &iflag, long unsigned int &_z){
    long unsigned int _j=IJLINKS.size();\
    int k1=0, k2=0;
    for(k1=0;k1<d;++k1) {if (BACKLIST[J][k1]>=i) break;}
    _j= BACKLIST[J][k1];
    swap(IJLINKS[_z], IJLINKS[_j]);
    BACKLIST[J][k1]=_z;
    for (k2=0;k2<d;++k2) {if (BACKLIST[I][k2]==_z) break;}
    BACKLIST[I][k2] = _j;
    I=IJLINKS[_z];
}

void _backlist(long unsigned int &i, int &d, vector<vector<long unsigned int> > &BACKLIST, vector<long unsigned int> &IJLINKS, long unsigned int &I, long unsigned int &J, long unsigned int &K, long unsigned int &iflag){
    long unsigned int k=0;
    long unsigned int _j=IJLINKS.size();
    int k1=0, k2=0;
    for ( k1 = 0;k1<d;++k1) if (BACKLIST[I][k1]>i) {
        iflag += 2;++k;/*k counts number of connections to fill*/
        if(iflag>=IJLINKS.size())break;
        /*swap this connection into position i+2k, */
        _j = BACKLIST[I][k1];
        J = IJLINKS[_j];/*now swap J with whatever's at location i + k*2 */

        K = IJLINKS[i+2*k];
        for (k2=0;k2<d;++k2){if (BACKLIST[K][k2] == i+2*k) break;};
        BACKLIST[K][k2] = _j;
      /*  if(k2<d-1){
            if(BACKLIST[K][k2]>BACKLIST[K][k2+1])swap(BACKLIST[K][k2],BACKLIST[K][k2+1]);
        }*/
        BACKLIST[J][k1] = i+2*k;
        if(_j>BACKLIST.size()*d)iflag=IJLINKS.size();
        IJLINKS[i+2*k] = I;IJLINKS[_j]=K;
        
    }

   
}

void _backlist(long unsigned int &i, int &d, vector<vector<long unsigned int> > &BACKLIST, vector<long unsigned int> &IJLINKS, long unsigned int &I, long unsigned int &J, long unsigned int &K, long unsigned int &iflag, long unsigned int &_z){
    long unsigned int k=0;
    long unsigned int _j=IJLINKS.size();
    long unsigned int cpiflag=iflag;
    int k1=0, k2=0;
    for ( k1 = 0;k1<d;++k1) if (BACKLIST[I][k1]>_z) {
      
        /*swap this connection into position i+2k, */
        _j = BACKLIST[I][k1];
        if(_j!=cpiflag-2){
        iflag += 2;++k;/*k counts number of connections to fill*/
        if(iflag>=IJLINKS.size())break;
        J = IJLINKS[_j];/*now swap J with whatever's at location i + k*2 */
        
        K = IJLINKS[i+2*k];
        for (k2=0;k2<d;++k2){if (BACKLIST[K][k2] == i+2*k) break;};
        BACKLIST[K][k2] = _j;BACKLIST[J][k1] = i+2*k;
        if(_j>BACKLIST.size()*d)iflag=IJLINKS.size();
        IJLINKS[i+2*k] = I;IJLINKS[_j]=K;
    }
    }
    
    
}


void update_backlist(long unsigned int &i, int &d, vector<vector<long unsigned int> > &BACKLIST, long unsigned int &temp, long unsigned int &J, long unsigned int &_j){
    int k1=0, k2=0;
    for(k1=0;k1<d;++k1) {if (BACKLIST[temp][k1]==(i+1)) break;}
    BACKLIST[temp][k1] = _j;
    for(k2=0;k2<d;++k2) {if (BACKLIST[J][k2]==_j) break;}
    BACKLIST[J][k2] = i+1;

}

void update_connectivity(vector<Contraction> &ptc, vector<list<Contraction *> > &listsofcontractsites, long unsigned int &I, long int &_max_in_list){
    ptc[I].update_connectivity(listsofcontractsites);
    if(_max_in_list<static_cast<int>(ptc[I]._connectivity)){
        _max_in_list=ptc[I]._connectivity;
    }
}

void _swap_single_element(long unsigned int &i, int &d, vector<vector<long unsigned int> > &BACKLIST, vector<long unsigned int> &IJLINKS, long unsigned int &I, long unsigned int &J, long unsigned int &K, long unsigned int &iflag){
    long unsigned int _z=i+2;
    long unsigned int _z2=i+2;
    
    I=IJLINKS[_z];
    iflag +=2;
    
    _myswap(_z2, d, BACKLIST, IJLINKS, I, J, iflag,_z);
    _backlist(_z, d, BACKLIST, IJLINKS, I, J, K, iflag, _z2);
    
}


// declaration and definition of functions for statistics.
template<typename T>
double _mean(vector<T> &a){ /*mean of the elements of the vector a */
    T sum=0;
    for (int i=0; i<a.size(); ++i) {
        sum+=a[i];
    }
    const double mean=(double)sum/(double)a.size();
    return mean;
}

template<typename T>
double _stdv(vector<T> &a, const double &mean){ /* standard deviation and standard deviation of the mean of the vector a */
    double sum=0.;
    double stdv=0.;
    for (int i=0; i<a.size(); ++i) {
        sum+=((double)a[i]-mean)*((double)a[i]-mean);
    }
    stdv=sqrt(sum/(double)(a.size()-1));
    //stdv=stdv/sqrt((double)(a.size()));
    return stdv;
}



// Definition of the algorithm

int main(int argc, const char * argv[]) {

    
    long unsigned int N=atoi(argv[argc-3]);
    int d=atoi(argv[argc-2]);
    ITER=atoi(argv[argc-1]);
    D=d;
    long unsigned int k=0, nlinks=0, I, Itemp=0, Jtemp=0, J = 0, K, _Jr=0, iflag = 0, connection=0, temp=-1;
    long long int nleft=0;
    long unsigned int loops = 0, hlinks = 0, _max_in_list;
    long unsigned int _j=-1;
    long double IS;
    long double VC;
    vector<long unsigned int> ADJLIST; /*vector which is the adjacency matrix*/
    vector<long unsigned int> IS_vec;
    Random MyRand;
    MyRand.print_seed(); /* seed random number generator */
    Graph myGraph;
    myGraph.Get_Values_Graph(N, d);
    vector<long unsigned int> FIRST, LAST, IJLINKS;
    vector<vector<long unsigned int> > BACKLIST;
    list<Vertex> listof2;
    list<Vertex *> listof3;
    vector<list<Contraction *> > listsofcontractsites;
    vector<Vertex> ptv;
    vector<Contraction> ptc;
    list<Vertex *> listtocomplete;
    vector<long double> v;
    list<Vertex *> _list_ngb_tuched;
    list<Vertex *> possible_2;
    v.resize(ITER);
    Vertex *ptrV;
    list<Contraction *> _list_of_C;
    vector<list<Vertex *> > vec_list;
    bool flag_list=false;
    bool flag_contractor=false;
    bool flag_max_in_list=false;
    bool flag_max_in_list1=false;
    bool flag2=false;
    bool flagb4=false;
    bool flagI=false;
    bool flag_I_met_a_contractor=false;
    bool flag_door=false;
    unsigned int H=0;
    long unsigned int MAX=0;
    long unsigned int counter__=0;
    ostringstream n,_d;
    n<<N;
    _d<<d;
    string title="Mean";
    string Nt="N=";
    string txt=".txt";
    string _str;
    _str=directory;
    _str+=title;
    _str+=Nt;
    _str+=n.str();
    _str+=txt;
    ofstream outfile(_str.c_str(),ios_base::app);
    string title_anal="Anal_d=";
    string _str_anal;
    _str_anal=directory;
    _str_anal+=title_anal;
    _str_anal+=_d.str();
    _str_anal+=Nt;
    _str_anal+=n.str();
    _str_anal+=txt;
    
    ofstream outfile_anal(_str_anal.c_str(),ios_base::app);
    // for loop over the number of graphs analyzed.
    for (int s=0; s<ITER; ++s) {
START:
      
    cout<<"I START THE ANALYSIS"<<endl;
    // here the algorithm initialises all variables and vectors.
    ptrV=NULL;
    flag_max_in_list=false;
    counter__=0;
    IS=0;
    VC=0;
    myGraph.Get_Values_Graph(N, d);
    IS_vec.clear();
    listof3.clear();
    _list_ngb_tuched.clear();
    possible_2.clear();
    listtocomplete.clear();
    IS_vec.reserve(N);
    listof2.clear();
    if(!listsofcontractsites.empty())
        for (unsigned int i=0; i<MAX_LIST; ++i) {
            listsofcontractsites[i].clear();
        }
        if(d>=4){
        if(!vec_list.empty())
        for (unsigned int i=0; i<d; ++i) {
            vec_list[i].clear();
        }
    vec_list.clear();
    vec_list.resize(d);
         }
    listsofcontractsites.clear();
    
    FIRST.clear();
    LAST.clear();
    IJLINKS.clear();
        if(!BACKLIST.empty()){
            for (long unsigned int i=0; i<BACKLIST.size(); ++i) {
                BACKLIST[i].clear();
            }
        }
    
        if(!ADJLIST.empty())ADJLIST.clear();
       
        if(!ptc.empty()){
            ptc.clear();
        }
        if(!ptv.empty())ptv.clear();
    listsofcontractsites.resize(MAX_LIST);
    FIRST.resize(N+1);
    LAST.resize(N+1);
    IJLINKS.resize(d*N);
    BACKLIST.resize(N);
    ADJLIST.resize(N*d);
    ptv.resize(N);
    ptc.resize(N);
    _list_of_C.clear();
    _j=-1;
    _Jr=0;
    loops = 0;
    hlinks = 0;
    _max_in_list=0;
    k=0;
    H=0;
    nlinks=0; nleft=0; iflag = 0; connection=0; temp=-1;

    flag_contractor=false;
    flag_max_in_list=false;
    flag_max_in_list1=false;
    flag2=false;
    flag_list=false;
    flagb4=false;
    flag_I_met_a_contractor=false;
    flag_door=false;
    for (long unsigned int i=0; i<N; ++i) {
        
        ptv[i]._vertex_i=i; // I am labeling the elements
        ptv[i]._list_ngb.resize(d); // I am preparing the list of ngb;
        FIRST[i] = LAST[i] = i*d;
        BACKLIST[i].resize(d);
        for (long unsigned int j=0;j<d;++j) {
            IJLINKS[k] = i;BACKLIST[i][j]=k;
            ptv[i]._list_ngb[j]=-1;
            ++k;
        }
    }

    FIRST[N] = d*N;
    nlinks =d*N;
    nleft = nlinks;
        
   // here the algorithm starts
   if(d==3)ptv[0]._who_am_I='V';
   if(d>=4)ptv[0]._who_am_I='I';
   // build the graph and find the maximal independent set
    for (long unsigned int i=0;i<nlinks-2;i+=2){
        flagI=true;
        I = IJLINKS[i];
        if (i<iflag) goto CONNECT;
        _list_ngb_tuched.clear();
        iflag +=2;
        if(I==0){
            ptv[I]._vertex_i=I;
            _backlist(i, d, BACKLIST, IJLINKS, I, J, K, iflag);
        }else{
            if(!listof2.empty()){
                long unsigned int F=N+1;
                for (list<Vertex>::iterator it=listof2.begin(); it!=listof2.end(); ++it) {
                    if(it->_who_am_I==' ' and it->_ref_char_val==' '){
                        F=it->_vertex_i;
                        listof2.pop_front();
                        ptv[F]._I_am_a_two=false;
                        break;
                    }else{
                        if(listof2.empty())break;
                       listof2.pop_front();
                    }
                }
                J=F;
                if(!ptv[J]._I_am_contracted){
                    ptv[J]._ref_char_val='P';
                    ptv[J]._I_am_contracted=true;
                    ptv[J]._I_am_the_contractor=true;
                    ptv[J].My_contractor_is=&ptv[J];
                    ptc[J].Get_Contractor(ptv[J]);
                }
                _myswap(i, d, BACKLIST, IJLINKS, I, J, iflag);
                _backlist(i, d, BACKLIST, IJLINKS, I, J, K, iflag);
                if(iflag>=IJLINKS.size()) {
                    iflag-=2;
                }
            }
        }
       
        
    CONNECT:
        // connect the site with d ngbs
        /*START CHECK VALUE J*/
        nleft = nlinks - iflag;
        if(nleft==0){
            break;
        }
        if (loops + hlinks>2*N*d) {
            cout<<"error "<<loops<<" loops, "<<hlinks<<" hlinks"<<endl;
            goto START;
        }

        _j =  (MyRand.random_number()*nleft)+iflag;
        
        if(_j==IJLINKS.size()){
            goto START;
            
        }
        J = IJLINKS[_j];
        if(I == J) {++loops;
            if(loops>3*N) goto START;
            goto CONNECT;};
        
        for (k=FIRST[I];k<LAST[I];++k) if (ADJLIST[k] == J){
            ++hlinks;
            if(hlinks>3*N) goto START;
            goto CONNECT;};
        
        if(ptv[J]._refdeg==d){
            ++counter__;
            if(counter__>3*N) goto START;
            goto CONNECT;
        }
        if(ptv[I]._refdeg==d){
            cout<<"1 "<<i<<endl;
            cout<<ptv[I];
            exit(-1);
        }
        /*END CHECK VALUE J*/
        
        if(ptv[I]._I_am_a_three){
            listof3.erase(ptv[I].it_list3);
            ptv[I]._I_am_a_three=false;
        }
        if(ptv[J]._I_am_a_three){
            listof3.erase(ptv[J].it_list3);
            ptv[J]._I_am_a_three=false;
        }
        ADJLIST[LAST[I]]=J;
        ++LAST[I];
        ADJLIST[LAST[J]]=I;
        ++LAST[J];
        
        _list_ngb_tuched.push_front(&ptv[J]);
        connection+=2;
        if(d>=4){
        if(ptv[I]._I_am_in_vec_list){
            ptv[I]._I_am_in_vec_list=false;
            vec_list[d-ptv[I]._refdeg].erase(ptv[I].it_vec_list);
        }
        
        if(ptv[J]._I_am_in_vec_list){
            ptv[J]._I_am_in_vec_list=false;
            vec_list[d-ptv[J]._refdeg].erase(ptv[J].it_vec_list);
        }
       }
        ptv[I]._list_ngb[ptv[I]._refdeg++]=J;
        ptv[J]._list_ngb[ptv[J]._refdeg++]=I;
        ptv[J]._vertex_i=J;

        

        
       
        /*Update connectivity contractor*/
        if(ptv[I]._who_am_I=='V' and ptv[J]._I_am_contracted and ptv[J]._who_am_I==' '){
            ptc[ptv[J].My_contractor_is->_vertex_i].update_connectivity(listsofcontractsites);
        }
        
        /*select two's*/
        if(ptv[I]._who_am_I=='V' and (d-ptv[J]._refdeg)<=2 and !ptv[J]._I_am_contracted){
            if(ptv[J]._I_am_a_two){
                listof2.erase(ptv[J].it_list2);
                    if((d-ptv[J]._refdeg)==0){
                        ptv[J]._who_am_I='I';
                        flagI=false;
                        ptv[J]._I_am_a_two=false;
                    }else{
                        listof2.push_front(ptv[J]);
                        ptv[J]._I_am_a_two=true;
                        ptv[J].it_list2=listof2.begin();
                    }
            }else{
                ptv[J]._I_am_a_two=true;
                listof2.push_front(ptv[J]);
                ptv[J].it_list2=listof2.begin();
            }
        }
        /*Contract sites and make a new sites*/
        if(ptv[I]._ref_char_val=='P' && ptv[I]._who_am_I==' '){
            if(ptv[J]._I_am_a_two) {
                ptv[J]._I_am_a_two=false;
                listof2.erase(ptv[J].it_list2);
            }
            if(!ptv[J]._I_am_contracted){
                ptv[J]._I_am_contracted=true;
                ptv[J]._ref_char_val='C';
                ptv[J].My_contractor_is=ptv[I].My_contractor_is;
                ptc[ptv[J].My_contractor_is->_vertex_i].Contract(ptv[J]);
            }else{
                if(ptv[J].My_contractor_is->_vertex_i!=ptv[I].My_contractor_is->_vertex_i){
                K=ptv[I].My_contractor_is->_vertex_i;
                ptc[ptv[J].My_contractor_is->_vertex_i]._fusion_c(ptc[K]);
                    if(ptc[K]._I_am_in_the_list){
                        ptc[K]._I_am_in_the_list=false;
                        listsofcontractsites[ptc[K]._connectivity].erase(ptc[K].ptr_vec_list);
                    }
                ptc[K].Clean();
                ptc[ptv[J].My_contractor_is->_vertex_i].check();
                ptv[I]._I_am_the_contractor=false;
                ptv[K]._I_am_the_contractor=false;
                ptv[I]._I_am_contracted=true;
                ptv[I].My_contractor_is=ptv[J].My_contractor_is;
                ptv[K].My_contractor_is=ptv[J].My_contractor_is;
                ptc[ptv[J].My_contractor_is->_vertex_i].update_connectivity(listsofcontractsites);
                }else{
                    if(ptv[J]._ref_char_val=='P'){
                        if(ptv[J]._I_am_a_two){
                            ptv[J]._I_am_a_two=false;
                            listof2.erase(ptv[J].it_list2);
                        }
                    ptc[ptv[J].My_contractor_is->_vertex_i].swap_P_C();
                    ptc[ptv[J].My_contractor_is->_vertex_i].fix_IS_VC();
                    ptv[ptv[J].My_contractor_is->_vertex_i]._I_am_the_contractor=false;
                    ptc[ptv[J].My_contractor_is->_vertex_i].Clean();
                        ptc[ptv[J].My_contractor_is->_vertex_i].check();
                    ptv[I]._I_am_the_contractor=false;
                    ptv[I]._I_am_contracted=false;
                    flagI=false;
                    }else{
                        ptc[ptv[I].My_contractor_is->_vertex_i].update_connectivity(listsofcontractsites);
                        if(ptc[ptv[J].My_contractor_is->_vertex_i]._connectivity==0){
                            ptc[ptv[I].My_contractor_is->_vertex_i]._I_am_in_the_list=false;
                            listsofcontractsites[0].erase( ptc[ptv[I].My_contractor_is->_vertex_i].ptr_vec_list);
                        ptc[ptv[J].My_contractor_is->_vertex_i].fix_IS_VC();
                        ptv[ptv[J].My_contractor_is->_vertex_i]._I_am_the_contractor=false;
                        ptc[ptv[J].My_contractor_is->_vertex_i].Clean();
                            ptc[ptv[J].My_contractor_is->_vertex_i].check();
                        ptv[I]._I_am_the_contractor=false;
                        ptv[I]._I_am_contracted=false;
                        flagI=false;
                        }else{
                        ptv[J]._I_am_contracted=true;
                        ptv[J]._ref_char_val='C';
                        }
                    }
        
                }
            }
            if(ptv[I]._refdeg==d && ptv[I]._I_am_contracted){
                ptc[ptv[I].My_contractor_is->_vertex_i].update_connectivity(listsofcontractsites);
                if(_max_in_list<ptc[ptv[I].My_contractor_is->_vertex_i]._connectivity){
                    _max_in_list=ptc[ptv[I].My_contractor_is->_vertex_i]._connectivity;
                }
            }
        }
        

        

        

        Itemp=I;
        Jtemp=J;
        
     temp = IJLINKS[i+1];IJLINKS[i+1] = J;IJLINKS[_j] = temp;
     update_backlist(i, d, BACKLIST, temp, J, _j);
        
        if(d>=4){
            if(!ptv[J]._I_am_in_vec_list and ptv[J]._ref_char_val==' '){
                ptv[J]._I_am_in_vec_list=true;
                vec_list[d-ptv[J]._refdeg].push_front(&ptv[J]);
                ptv[J].it_vec_list=vec_list[d-ptv[J]._refdeg].begin();
            }
            
            if(!ptv[I]._I_am_in_vec_list and ptv[I]._ref_char_val==' '){
                ptv[I]._I_am_in_vec_list=true;
                vec_list[d-ptv[I]._refdeg].push_front(&ptv[I]);
                ptv[I].it_vec_list=vec_list[d-ptv[I]._refdeg].begin();
            }
        }


        if(ptv[J]._I_am_contracted){
                ptc[ptv[J].My_contractor_is->_vertex_i].update_connectivity(listsofcontractsites);
        }
        
        
        if(ptv[I]._who_am_I=='I' and flagI){
            ptv[J]._who_am_I='V';
            if(ptv[J]._I_am_a_two){
                ptv[J]._I_am_a_two=false;
                listof2.erase(ptv[J].it_list2);
            }
            if(ptv[J]._ref_char_val!='C')ptv[J]._ref_char_val='C';
            if(ptv[J]._I_am_contracted){
                if(ptv[J]._I_am_the_contractor)exit(-1);
                ptc[ptv[J].My_contractor_is->_vertex_i].erase(J);
                ptc[ptv[J].My_contractor_is->_vertex_i].update_connectivity(listsofcontractsites);
                ptv[J]._I_am_contracted=false;
            }
            if(ptv[J]._refdeg!=d){
            listtocomplete.push_front(&ptv[J]);
            }
        }
        
        
       
        if(flagb4 and (i==iflag-2 and iflag<nlinks-2)){
            
          
            bool _I_have_met_a_two=false;
                ptv[I]._who_am_I='I';
                ptv[I]._refdeg='F';
                ptv[I]._I_am_in_vec_list=false;
                vec_list[0].erase(ptv[I].it_vec_list);
                 for (list<Vertex *>::iterator _l_ngb=_list_ngb_tuched.begin(); _l_ngb!=_list_ngb_tuched.end(); ++_l_ngb) {
                     if((*_l_ngb)->_who_am_I!='I'){
                         (*_l_ngb)->_who_am_I='V';
                         if((*_l_ngb)->_I_am_a_two){
                             (*_l_ngb)->_I_am_a_two=false;
                             listof2.erase((*_l_ngb)->it_list2);
                         }
                         if((*_l_ngb)->_ref_char_val!='C')(*_l_ngb)->_ref_char_val='C';
                         if((*_l_ngb)->_I_am_contracted){
                             if((*_l_ngb)->_I_am_the_contractor)exit(-1);
                             ptc[(*_l_ngb)->My_contractor_is->_vertex_i].erase((*_l_ngb)->_vertex_i);
                             ptc[(*_l_ngb)->My_contractor_is->_vertex_i].update_connectivity(listsofcontractsites);
                             (*_l_ngb)->_I_am_contracted=false;
                         }
                         if((*_l_ngb)->_refdeg!=d){
                             listtocomplete.push_front((*_l_ngb));
                         }
                         
                     }else{
                         cout<<"I met an I"<<endl;
                         exit(-1);
                     }
                 }

            flagb4=false;
             _list_ngb_tuched.clear();
        }
        
        if(i==iflag-2 and iflag<nlinks-2 and listof2.empty()){
            
            I=Itemp;
            J=Jtemp;
        FINDANDFIX:
            for (list<Vertex *>::iterator _l_ngb=_list_ngb_tuched.begin(); _l_ngb!=_list_ngb_tuched.end(); ++_l_ngb) {
                if((*_l_ngb)->_I_am_contracted){
                    if(ptc[(*_l_ngb)->My_contractor_is->_vertex_i]._I_am_in_the_list){
                        ptc[(*_l_ngb)->My_contractor_is->_vertex_i].update_connectivity(listsofcontractsites);
                    }
                }
            }
            _list_ngb_tuched.clear();
            if(listtocomplete.empty()){
            bool flag_max_d4=true;
            list<Contraction *>::iterator _it;
            if(!listsofcontractsites[0].empty()){
              
            for (_it=listsofcontractsites[0].begin(); _it!=listsofcontractsites[0].end(); ++_it) {
                if((*_it)->_C>(*_it)->_P){
                    if((*_it)->_I_have_a_loop_inside)exit(-1);
                    (*_it)->swap_P_C();
                }
                (*_it)->fix_IS_VC();
                (*_it)->Clean();
                }
          
            listsofcontractsites[0].clear();
            }
            
            if(!listsofcontractsites[1].empty()){
              
                 flag_max_d4=false;
                _it=listsofcontractsites[1].begin();
                  long unsigned int _z=i+2;
                 long unsigned int _z2=i+2;
                if((*_it)->_I_have_a_loop_inside)exit(-1);
                (*_it)->swap_P_C();
                (*_it)->fix_IS_VC();
                J=(*_it)->Find_V_to_fill(N);
                if(J==N+1){
                    listsofcontractsites[1].erase((*_it)->ptr_vec_list);
                    (*_it)->Clean();
                    goto FINDANDFIX;
                }
                I=IJLINKS[_z];
    
                iflag +=2;
                _myswap(_z2, d, BACKLIST, IJLINKS, I, J, iflag,_z);
                _backlist(_z, d, BACKLIST, IJLINKS, I, J, K, iflag, _z2);
                (*_it)->_I_am_in_the_list=false;
                listsofcontractsites[1].erase((*_it)->ptr_vec_list);
                (*_it)->Clean();
                if(iflag>=IJLINKS.size()) {
                    iflag -=2;
                }
            }else if(!listsofcontractsites[2].empty()){
            
                 flag_max_d4=false;
                _it=listsofcontractsites[2].begin();
                if((*_it)->_I_have_a_loop_inside)exit(-1);
                (*_it)->swap_P_C();
                list<long unsigned int> l_J;
                l_J=(*_it)->Find_V_to_fill_vec();
                if(l_J.empty()){
                    (*_it)->fix_IS_VC();
                    listsofcontractsites[2].erase((*_it)->ptr_vec_list);
                    (*_it)->Clean();
                       goto FINDANDFIX;
                }
            
                
                long unsigned int _z=i+2;
                     long unsigned int _z2=i+2;
                if(_z==IJLINKS.size())goto END;
                for (list<long unsigned int>::iterator lit=l_J.begin(); lit!=l_J.end(); ++lit) {
                    J=*lit;
                    I=IJLINKS[_z];
                    iflag +=2;// questo e' il punto pericoloso perche' troppe operazioni faccio
                    //devo migliorare questo punto
                    _myswap(_z2, d, BACKLIST, IJLINKS, I, J, iflag, _z);
                    _backlist(_z, d, BACKLIST, IJLINKS, I, J, K, iflag, _z2);//qui ho il bug perche' sto
                    if(iflag>=IJLINKS.size()) {
                        iflag-=2;
                    }
                    _z=iflag;
                }
            }else{
                
            while (listsofcontractsites[_max_in_list].empty()){
                if(_max_in_list==0)break;
                --_max_in_list;
            }
                if(_max_in_list!=0){
                  
                    if(d==3 or d>=4){
                        
                     flag_max_d4=false;
                        flag_door=false;
                    _it=listsofcontractsites[_max_in_list].begin();
                    (*_it)->fix_IS_VC();
                        Vertex *cp;
                        cp=(*_it)->Give_Me_Contractor_ptr();
                        const long unsigned int K_1=cp->_vertex_i;
                       // cout<<*cp<<endl;
                    long unsigned int _z=i+2;
                        long unsigned int _z2=i+2;
                        if(_z==IJLINKS.size())goto END;
                    for (list<Vertex *>::iterator lit=(*_it)->_VC_to_fill.begin(); lit!=(*_it)->_VC_to_fill.end(); ++lit) {
                         if(_z==IJLINKS.size())goto END;
                            J=(*lit)->_vertex_i;
                            I=IJLINKS[_z];
                            iflag +=2;
                         _myswap(_z2, d, BACKLIST, IJLINKS, I, J, iflag,_z);
                        _backlist(_z, d, BACKLIST, IJLINKS, I, J, K, iflag, _z2);
                        if(iflag>=IJLINKS.size()){
                            iflag -=2;
                        }
                            _z=iflag;
                       
                        
                 }
                (*_it)->_I_am_in_the_list=false;
                listsofcontractsites[_max_in_list].erase((*_it)->ptr_vec_list);
                ptc[K_1].Clean();
                    }
                }
            }
                
                if(d>=4 and flag_max_d4){
                    list<Vertex *>::iterator it;
                    if(!vec_list[0].empty()){
                        while (!vec_list[0].empty()) {
                            it=vec_list[0].begin();
                            bool flag_for=true;
                            for (long unsigned int h=0; h<d; ++h) {
                                if(ptv[(*it)->_list_ngb[h]]._who_am_I!='V'){
                                    flag_for=false;
                                    break;
                                }
                            }
                            if(flag_for){
                                (*it)->_who_am_I='I';
                            }
                            (*it)->_I_am_in_vec_list=false;
                            vec_list[0].erase(it);
                        }
                    }
                    int counteer_empty=0;
                    for (unsigned int i_=0; i_<d; ++i_) {
                        if(!vec_list[i_].empty()){
                            it=vec_list[i_].begin();
                            K=(*it)->_vertex_i;
                            ptv[K]._I_am_in_vec_list=false;
                            vec_list[i_].erase(it);
                            if(ptv[K]._ref_char_val==' ')break;
                        }else{
                            ++counteer_empty;
                        }
                    }
                    if (counteer_empty==d) {
                        flag_door=true;
                        goto FINDANDFIX;
                    }
                    _swap_single_element(i, d, BACKLIST, IJLINKS, I, K, J, iflag);
                   flagb4=true;
                }

            }else{
                long unsigned int _z=i+2;
                long unsigned int _z2=i+2;
                for (list<Vertex *>::iterator lit=listtocomplete.begin(); lit!=listtocomplete.end(); ++lit) {
                    if(_z==IJLINKS.size())goto END;
                    J=(*lit)->_vertex_i;
                    I=IJLINKS[_z];
                    iflag +=2;
                    _myswap(_z2, d, BACKLIST, IJLINKS, I, J, iflag, _z);
                    _backlist(_z, d, BACKLIST, IJLINKS, I, J, K, iflag,_z2);
                    if(iflag>=IJLINKS.size()) {
                        iflag -=2;
                    }
                    _z=iflag;
                }
                listtocomplete.clear();
            }
        }
    
        
    }
    END:
    
    I = IJLINKS[nlinks-2];J = IJLINKS[nlinks-1];
    ADJLIST[LAST[I]]=J;++LAST[I];
    ADJLIST[LAST[J]]=I;++LAST[J];
    connection+=2;
        
        IS=0;
        for (long unsigned int i=0; i<N; ++i) {
            if(ptv[i]._who_am_I=='I')IS++;
        
        }
        
        // print the number of IS found
        cout<<IS<<endl;
        // check that the list of virtual nodes is empty, if not empty then fix all nodes
        list<Contraction *>::iterator _it;
        for (long unsigned int i=0; i<MAX_LIST; ++i) {
            while (!listsofcontractsites[i].empty()) {
                list<Contraction *>::iterator _it;
                _it=listsofcontractsites[i].begin();
                (*_it)->fix_IS_VC();
                listsofcontractsites[i].pop_front();
                
            }
        }
        
    if (I==J){ printf("\nloop on last pair");
        cout<<endl;
        cout<<"WARNING"<<endl;
        
    }else{
        ptv[I]._list_ngb[ptv[I]._refdeg++]=J;
        ptv[J]._list_ngb[ptv[J]._refdeg++]=I;
        ptv[J]._vertex_i=J;
    }
    
    cout<<"nloops: "<<loops<<" h-links"<<hlinks<<endl;
    
    
    for (int i=1;i<d;++i) {
        if (LAST[i-1]!=FIRST[i])
            cout<<"last != first "<<i<<" "<<LAST[i-1]<<" "<<FIRST[i]<<" "<<LAST[i]<<endl;
        
    };
    
    for (k = 0;k<N;++k) {
        for (int k1 = 0;k1<d;++k1) {if (IJLINKS[ BACKLIST[k][k1]] != k){
            cout<<" k, k1, BACKLIST[k][k1]"<<k<<" "<<k1<<" "<<BACKLIST[k][k1]<<endl;
             goto START;
        }
            };};
    
        IS=0;
        
        VC=0;
        IS_vec.clear();
        IS_vec.reserve(N);
        for (long unsigned int i=0; i<N; ++i) {
            if(ptv[i]._who_am_I=='I'){
                ++IS;
                IS_vec.push_back(i);
            }else if(ptv[i]._who_am_I=='V'){
                ++VC;
            }else if(ptv[i]._who_am_I==' '){
                for (int a=0; a<d;++a) {
                    if(ptv[ptv[i]._list_ngb[a]]._vertex_i!=-1){
                    if(ptv[ptv[i]._list_ngb[a]]._who_am_I=='I'){
                        ptv[i]._who_am_I='V';
                        ++VC;
                        break;
                    }
                    }else{
                        ptv[i]._who_am_I='V';
                        ++VC;
                        break;
                    }
                    
                }
                if(ptv[i]._who_am_I==' '){
                    ptv[i]._who_am_I='I';
                    ++IS;
                    IS_vec.push_back(i);
                }
            }

            
        }
    // end part checking that all nodes are labelled V or I
        
    vector<long unsigned int> (IS_vec).swap(IS_vec);
    // check that the solution found is a solution, in other words we check that all sites labbeled with letter I are covered by only V sites.
    for (long unsigned int i=0; i<IS_vec.size(); ++i) {
        for (long unsigned int j=0; j<d; ++j) {
            if(ptv[IS_vec[i]]._list_ngb[j]!=-1){
            if(ptv[ptv[IS_vec[i]]._list_ngb[j]]._who_am_I=='I'){
                cout<<IS_vec[i]<<" "<<ptv[IS_vec[i]]._list_ngb[j]<<" "<<ptv[IS_vec[i]]._who_am_I<<" "<<ptv[ptv[IS_vec[i]]._list_ngb[j]]._who_am_I<<endl;
                cout<<"NO SOLUTION"<<endl;
                cout<<"I restart"<<endl;
                goto START;
            }
            }else{
                if(IS_vec[i]!=I||IS_vec[i]!=J){
                    cout<<"WARNING LOOPS"<<endl;
                    goto START;
                }
            }
        }
    }
    cout<<"YES SOLUTION"<<endl;
        cout<<s<<" "<<IS<<" "<<VC<<" "<<IS+VC<<" "<<IS/(IS+VC)<<endl;
        v[s]=IS/(IS+VC);
    }
    const double mean=_mean<long double>(v);
    const double stdv=_stdv<long double>(v, mean);
    cout<<"MEAN: "<<mean<<" STD: "<<stdv<<endl;
    outfile<<N<<" "<<d<<" "<<mean<<" "<<stdv<<endl;
    cout<<"END"<<endl;
    return 0;
}

