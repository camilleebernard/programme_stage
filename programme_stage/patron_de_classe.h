//  patron_de_classe.h
//  programme_stage
//
//  Created by Camille Bernard on 27/07/2021.
//

//Voir comment utiliser export pour séparer la déclaration de la définition des patrons. (Delannoy p369)

#ifndef patron_de_classe_h
#define patron_de_classe_h

#include <iostream>
#include <complex>
using namespace std;

#include <math.h>
#include<float.h>

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

#define H_P_COOD 0
#define N_N_PRIM_COOD 1

//def du parametre expression C
#define BASE 0
#define PLUS 1
#define MOINS 2
#define MULT 3
#define DIV 4



template <class A,class B>
class application;
template <class A,class B>
class un;
template <class A,class B>
class constante;


/*
template <class A,class B>
application<A,B>  operator+ (const application<A,B> &,const application<A,B> &);
template <class A,class B>
application<A,B>  operator- (const application<A,B> &,const application<A,B> &);
template <class A,class B>
application<A,B>  operator* (const application<A,B> &,const application<A,B> &);
template <class A,class B>
application<A,B>  operator/ (const application<A,B> &,const application<A,B> &);
*/

template <class A,class B> class application {
protected:
    string nom;
    pair<const application<A,B>*,const application<A,B>*> entrees;
    int cat;
    //vector<double> parametres_reels;
    //vector<complex<double> > parametres_complexes;
public:
    //constructeur pour une conversion B -> application<A,B> (fct constante)
    //application(B);
    //constructeur
    application(const string);
    //application(const string,pair<un<A,B>*,un<A,B>*>,int);
    application(const string,pair<const application<A,B>*,const application<A,B>*>,int);
    
    //constructeur par recopie
    application(const application<A,B> &);
   //affectation
    application<A,B> & operator= (const application<A,B> &);
    
    //domaine_def
    virtual bool domaine_def(A);
    
    //friend application<A,B>  operator+<> (un<A,B> &,un<A,B> &);
    
    //friend application<A,B>  operator+<> (const application<A,B> &,const application<A,B> &);
    //friend application<A,B>  operator-<> (const application<A,B> &,const application<A,B> &);
    //friend application<A,B>  operator*<> (const application<A,B> &,const application<A,B> &);
    //friend application<A,B>  operator/<> (const application<A,B> &,const application<A,B> &);
    
    //obligé de déclarer les opérateur inline pour que les conversions implicites marchent
    
    inline friend application<A,B>  operator+ (const application<A,B> & f,const application<A,B> & g){
        pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
        application<A,B> h("h",entrees,PLUS);
        return h;
    }

    inline friend application<A,B>  operator- (const application<A,B> & f,const application<A,B> & g){
        pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
        application<A,B> h("f",entrees,MOINS);
        return h;
    }

    inline friend application<A,B>  operator* (const application<A,B> & f,const application<A,B> & g){
        pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
        application<A,B> h("f",entrees,MULT);
        return h;
    }
    
    inline friend application<A,B>  operator* (const constante<A,B> & f,const application<A,B> & g){
        cout<< "mult SCA_f" <<endl;
        const application<A,B>* adr_f = &f;//pour qu'il n'y ait pas de conversion implicite de f en application
        cout << "(*adr_f)(5)"<<(*adr_f)(5) << endl;
        pair<const application<A,B>*,const application<A,B>*> entrees(adr_f,&g);
        application<A,B> h("f",entrees,MULT);
        return h;
    }
    
    inline friend application<A,B>  operator/ (const application<A,B> & f,const application<A,B> & g){
        pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
        application<A,B> h("f",entrees,DIV);
        return h;
    }
    

    //Evaluation
    virtual B operator() (A) const; //fonction virtuelle et peut agir sur des objets constants
};

//Problème quand je mets les définitions des méthodes dans de le fichier .cpp ça fait des erreurs de liens.
// Début de solution ici: https://www.cs.technion.ac.il/users/yechiel/c++-faq/separate-template-class-defn-from-decl.html

template <class A, class B> application<A,B>::application(const string input_nom){
    nom = input_nom;
    cat = BASE;
}

/*
template <class A, class B> application<A,B>::application(B y){
    cout << "conv implicite" << endl;
    constante<A,B> f(y);
    *this = f;
}
*/

template <class A, class B> application<A,B>::application(const string input_nom,pair<const application<A,B>*, const application<A,B>*> input_entrees,int input_cat){
    nom = input_nom;
    cat = input_cat;
    entrees = input_entrees;
}

template <class A, class B> application<A,B>::application(const application<A,B> & input_application){
    nom = input_application.nom;
    cat = input_application.cat;
    entrees = input_application.entrees;
    //parametres_reels = input_application.parametres_reels;
    //parametres_complexes = input_application.parametres_complexes;
}
//rerediger
template <class A, class B> application<A,B> & application<A,B>::operator= (const application<A,B> & input_application){
    if(this!= &input_application){
        entrees = input_application.entrees;
        cat = input_application.cat;
    }
    return *this;
}

template <class A, class B> bool application<A,B>::domaine_def(A x){
    //Par defaut tout R ou tout C
    return true;
}

//possible de ne définir que des specifications de cette fonction -> erreur de compiation si on ne l'utilise pas avec le bon type -> plus difficile à reperer
// Je choisi de la définir mais il y aura une erreur si le constructeur du type/class B requiert des paramètres
template <class A,class B> B application<A,B>::operator() (A x) const {
    B result;
    cout << "nom"<< nom << endl;
    cout << "catégorie"<< cat << endl;
    switch(cat){
        case PLUS : result = (*entrees.first)(x)+ (*entrees.second)(x);
            return result;
        case MOINS : result = (*entrees.first)(x)- (*entrees.second)(x);
            return result;
        case MULT : result = ((*entrees.first)(x))* ((*entrees.second)(x));
            return result;
        case DIV :result = ((*entrees.first)(x))/ ((*entrees.second)(x));
            return result;
    }
    /*
    cout << "Erreur: le patron application<A,B> a été instancié avec des types A,B inconnus." << endl;
    B y = 1;
    return y;
     */
    cout << "Problème d'évaluation" << endl;
    B y = 1;
    return y;
}
/*
template <class A,class B> application<A,B>  operator+ (const application<A,B> & f,const application<A,B> & g){
    pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
    application<A,B> h("h",entrees,PLUS);
    return h;
}

template <class A,class B> application<A,B>  operator- (const application<A,B> & f,const application<A,B> & g){
    pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
    application<A,B> h("f",entrees,MOINS);
    return h;
}

template <class A,class B> application<A,B>  operator* (const application<A,B> & f,const application<A,B> & g){
    pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
    application<A,B> h("f",entrees,MULT);
    return h;
}

template <class A,class B> application<A,B>  operator/ (const application<A,B> & f,const application<A,B> & g){
    pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
    application<A,B> h("f",entrees,DIV);
    return h;
}
 */
/*----------------------------------------*/
//Un peu lourd: pas possible de specialiser une methode juste pour un paramètre de type ?

template <class A, class B> class un: public application<A,B>{
public:
    //constructeur
    un(const string);
    //Evaluation
    B  operator() (A) const;
    
};

template <class A> class un<A,double>: public application<A,double>{
public:
    //constructeur
    un(const string);
    
    //Evaluation
    double  operator() (A) const;
};

template <class A> class un<A,complex<double> >: public application<A,complex<double> >{
public:
    //constructeur
    un(const string);
    //Evaluation
    complex<double>  operator() (A) const;
};

template <class A,class B> un<A,B>::un(const string input_nom): application<A,B>(input_nom){}



template <class A> un<A,double>::un(const string input_nom): application<A,double>(input_nom){}
template<class A> double  un<A,double>::operator() (A z) const{
    double a = 1.;
    return a;
}


template <class A> un<A,complex<double> >::un(const string input_nom): application<A,complex<double> >(input_nom){}
template<class A> complex<double>  un<A, complex<double> >::operator() (A z) const{
    complex<double> a(1,0);
    return a;
}
/*----------------------------------------*/
//Cette classe sert uniquement à faire la conversion  B -> application<A,B>
//Il faut que le constructeur d'un objet de B ne requiert pas de paramètres.
template <class A, class B> class constante: public application<A,B>{
private:
    B y;
public:
    //constructeur
    //constante(const string,B);
    //Pour une conversion B -> constante<A,B>
    constante(B);
    //Evaluation
    B  operator() (A) const;
    
};
/*
template <class A,class B> constante<A,B>::constante(const string input_nom,B input_x): application<A,B>(input_nom){
    x = input_x;
}
 */
template <class A,class B> constante<A,B>::constante(B input_y): application<A,B>("fonction_constante"){
    cout << "conversion implicite" << endl;
    y = input_y;
}
template<class A,class B> B  constante<A,B>::operator() (A x) const{
    return y;
}

/*----------------------------------------*/

//En principe pas besoin de 0
/*
template <class A,class B> class zero: public application<A,B> {
public:
    //constructeur
    zero(const string);
    //Evaluation
    B  operator() (A) const;
};

template <class A> class zero<A,double>: public application<A,double> {
public:
    //constructeur
    zero(const string);
    //Evaluation
    double  operator() (A) const;
};
template <class A> class zero<A,complex<double> >: public application<A,complex<double> > {
public:
    //constructeur
    zero(const string);
    //Evaluation
    complex<double>  operator() (A) const;
};


template <class A,class B> zero<A,B>::zero(const string input_nom): application<A,B>(input_nom){}

template<class A> double  zero<A,double>::operator() (A z) const{
    double a = 0.;
    return a;
}

template<class A> complex<double>  zero<A,complex<double> >::operator() (A z) const{
    complex<double> a(0,0);
    return a;
}
*/
/*----------------------------------------*/
//mettre en parametre le type d'enveloppe


class paquet_d_onde: public application<double,complex<double> > {
    double mu_x; //moyenne en position
    double sigma_x;
    double mu_xi; //moyenne en freq
    public:
        //constructeur
        paquet_d_onde(const string,double,double,double);
        //Evaluation
        complex<double>  operator() (double) const;
};


paquet_d_onde::paquet_d_onde(const string input_nom,double input_mu_x,double input_sigma_x,double input_mu_xi): application<double,complex<double> >(input_nom){
    mu_x = input_mu_x;
    sigma_x = input_sigma_x;
    mu_xi = input_mu_xi;
}

complex<double>  paquet_d_onde::operator() (double x) const{
    complex<double> I(0,1);
    return (1./sqrt(sigma_x*sqrt(M_PI)))*exp(-(1./2)*pow((x-mu_x)/sigma_x,2)+I*mu_xi*(x-mu_x));
}

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

#endif /* patron_de_classe_h */

