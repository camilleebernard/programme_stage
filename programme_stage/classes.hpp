//
//  classes.hpp
//  programme_stage
//
//  Created by Camille Bernard on 22/07/2021.
//

#ifndef classes_hpp
#define classes_hpp


#include <iostream>
#include <complex>
using namespace std;

#include <math.h>

#include <armadillo>
//#include<boost/math/special_functions/bessel.hpp>
//using namespace arma;

#include<TApplication.h>
#include<TCanvas.h>
#include<TLine.h>
#include<TEllipse.h>
#include<TMarker.h>
#include<TBox.h>
#include<TText.h>
#include<TPad.h>
#include<TAxis.h>
#include<TGraph.h>
#include<TGraph2D.h>
#include<TMultiGraph.h>
#include<TF2.h>
#include<TPolyLine3D.h> //Graphe/lignes de niveaux d'une fonction

#include<TLatex.h>
#include<TMathText.h>

#include<TVector.h>
#include<TStyle.h>
#include<TH2F.h>
#include <TRandom3.h>

#include<TSlider.h>
#include<TSystem.h>
#include<TROOT.h>
#include<TBenchmark.h>


class ensemble {
private:
    string nom;
public:
    //constructeur
    ensemble(const string);
    //constructeur par recopie
    ensemble(const ensemble &);
    //affectation (obligé de définir ça pour faire des vecteurs d'ensembles)
    ensemble &  operator = (const ensemble &);
    //indicatrice (par defaut ensemble vide)
    virtual bool indicatrice(void);
};

class Reel: ensemble{
private:
    string nom;
public:
    //constructeur
    Reel(const string);
    //indicatrice
    virtual bool indicatrice(double);
};

class ouvert {
private:
    string nom;
public:
    //constructeur
    ouvert(const string);
    //constructeur par recopie
    ouvert(const ouvert &);
    //affectation
    ouvert &  operator = (const ouvert &);
    //destructeur
    //~ouvert(void);
};


class interval_ouvert : ouvert {
private:
    string nom;
    pair<double,double> bornes;
public:
    //constructeur
    interval_ouvert(const string, pair<double,double>);
    //constructeur par recopie
    interval_ouvert(const interval_ouvert &);
    //affectation
    interval_ouvert &  operator = (const interval_ouvert &);
    //destructeur
    //~ouvert(void);
};

class application {
protected:
    string nom;
    ouvert domaine_def;
    vector<double> parametres_reels;
    vector<complex<double> > parametres_complexes;
public:
    //constructeur
    application(const string,const ouvert);
    //constructeur par recopie
    application(const application &);
   //affectation
    application & operator= (const application &);
    //Arithmetique
    //friend application &  operator+ (const application &,const application &);
    //friend application &  operator- (const application &,const application &);
    //friend application &  operator* (const application &,const application &);
    //Evaluation
    virtual complex<double> operator() (complex<double>);
};

//Arithmetique
//application &  operator+ (const application &,const application &);
//application &  operator- (const application &,const application &);
//application &  operator* (const application &,const application &);

class un: application {
public:
    //constructeur
    un(const string,const ouvert);
    //Evaluation
    complex<double>  operator() (complex<double>);
};

class zero: application {
public:
    //constructeur
    zero(const string,const ouvert);
    //Evaluation
    complex<double>  operator() (complex<double>);
};

/*
 class carte: public ouvert {
protected:
    string name;
    vector<pair<ouvert,application> >;
    
public:
    //constructeur
    carte(const string);
    //constructeur par recopie
    carte(const carte &);
   //affectation
    carte &  operator = (const carte &);
    //destructeur
    ~carte(void);
};
*/




#endif /* classes_hpp */
