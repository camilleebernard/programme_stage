//
//  patron_de_classe_1.h
//  programme_stage
//
//  Created by Camille Bernard on 01/08/2021.
//

#ifndef patron_de_classe_1_h
#define patron_de_classe_1_h

#include <iostream>
#include <complex>
using namespace std;

#include<boost/type_index.hpp> //pour convertir paramètre de type->string

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


//def du parametre cat
/*
#define BASE 0
#define PLUS 1
#define MOINS 2
#define MULT 3
#define DIV 4
*/
//La norme préfère les constantes aux macro en C++ (par rapport au C)
/*
const int BASE = 0;
const int PLUS = 1;
const int MOINS = 2;
const int MULT = 3;
const int DIV = 4;
const int COMP = 5;
 */
//Obligé de les mettre dans une enum sinon le constructeur de application<A,B> est utilisé pour faire des conversions implicites int->application<A,B>
enum cat_appli{BUG,BASE,PLUS,MOINS,MULT,DIV,COMP};


template <class A,class B>
class application;
template <class A,class C,class B>
class compo_externe;
template <class A,class B>
class compo_interne;
template <class A,class B>
class constante;

//chercher une fonction qu'on puisse déclarer comme virtuelle pure pour rendre cette classe abstraite.

template <class A,class B> class application{
protected:
    cat_appli cat;
    application<A,B> * adr_def = NULL;
    //Ce constructeur n'est pas censé être utilisé par l'utilisateur.
    application(application<A,B> *, cat_appli);
public:
    //constructeur pour une conversion B -> application<A,B> (fct constante)
    //application(B);
    inline application(void){
#ifdef COMMENTAIRES
        cout << "-----------------------------------" << endl;
        cout << "Création d'une application de " << boost::typeindex::type_id<A>().pretty_name() << " vers "<< boost::typeindex::type_id<B>().pretty_name()<< endl;
        cout << "  à l'adresse " << this << endl;
        cout << "Cette application n'a pas de définition." << endl;
        cout << "-----------------------------------" << endl;
#endif
        cat = BUG;
        adr_def = NULL;
    }
    //constructeur par recopie
    application(const application<A,B> &);
    //affectation
    application<A,B> & operator= (const application<A,B> &);
    //domaine_def
    virtual bool domaine_def(A);
    
    //obligé de déclarer les opérateur inline pour que les conversions implicites marchent
    
    inline friend compo_interne<A,B> & operator+ (const application<A,B> & f,const application<A,B> & g){
        pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
        // création dynamique de l'objet qui est retourné par référence
        // il faut prévoir de le suprimer dans le main ou dans une autre fonction
        compo_interne<A,B> * adr_h;
        adr_h = new compo_interne<A,B>(entrees,PLUS);
        return *adr_h;
    }

    inline friend compo_interne<A,B> & operator- (const application<A,B> & f,const application<A,B> & g){
        pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
        // création dynamique de l'objet qui est retourné par référence
        // il faut prévoir de le suprimer dans le main ou dans une autre fonction
        compo_interne<A,B> * adr_h;
        adr_h = new compo_interne<A,B>(entrees,MOINS);
        return *adr_h;
    }

    inline friend compo_interne<A,B> & operator* (const application<A,B> & f,const application<A,B> & g){
        pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
        // création dynamique de l'objet qui est retourné par référence
        // il faut prévoir de le suprimer dans le main ou dans une autre fonction
        compo_interne<A,B> * adr_h;
        adr_h = new compo_interne<A,B>(entrees,MULT);
        return *adr_h;
    }
    
    inline friend compo_interne<A,B> & operator* (const constante<A,B> & f,const application<A,B> & g){
        cout<< "mult SCA_f" <<endl;
        const application<A,B>* adr_f = &f;//pour qu'il n'y ait pas de conversion implicite de f en application
        pair<const application<A,B>*,const application<A,B>*> entrees(adr_f,&g);
        // création dynamique de l'objet qui est retourné par référence
        // il faut prévoir de le suprimer dans le main ou dans une autre fonction
        compo_interne<A,B> * adr_h;
        adr_h = new compo_interne<A,B>(entrees,MULT);
        return *adr_h;
    }
    
    inline friend compo_interne<A,B> & operator/ (const application<A,B> & f,const application<A,B> & g){
        pair<const application<A,B>*,const application<A,B>*> entrees(&f,&g);
        // création dynamique de l'objet qui est retourné par référence
        // il faut prévoir de le suprimer dans le main ou dans une autre fonction
        compo_interne<A,B> * adr_h;
        adr_h = new compo_interne<A,B>(entrees,DIV);
        return *adr_h;
    }
    

    //Evaluation
    inline virtual B operator() (A a) const {
        //Voir si ce test ralenti beaucoup l'évaluation
#ifdef COMMENTAIRES
        cout << "Evaluation de l'application<" << boost::typeindex::type_id<A>().pretty_name()<< ","<<boost::typeindex::type_id<B>().pretty_name() << "> à d'adresse " << this << endl;
#endif
        if(adr_def!=NULL){
            return (*adr_def)(a);
        }else{
            cout << "Problème d'évaluation: adr_def==NULL" << endl;
            B b;
            return b;
        }
    };
    //affichage
    //Obligé de creer une fonction info parcequ'on ne peut pas utiliser le typage dynamique avec une fonction amie
    inline virtual ostream & info(ostream & os) const {
        //Boost a une fonction qui permet  de mieux afficher les paramètres de type
        //A tester si ça marche pour des types défini par l'utilisateur
        
        //os << ": "<< typeid(A).name() << "->"<< typeid(B).name() << endl << "cat : " << f.cat << endl << "entrees : (" << f.entrees.first << "," << f.entrees.second << ")";
        os << ": "<< boost::typeindex::type_id<A>().pretty_name()  << "->"<< boost::typeindex::type_id<B>().pretty_name()  << endl;
        return os;
    }
    /*
    inline friend ostream & operator<<(ostream & os,const application<A,B> & f){
        ostream is = os;
        os << f.info(is);
        return os;
    }
     */
};

//Problème quand je mets les définitions des méthodes dans de le fichier .cpp ça fait des erreurs de liens.
// Début de solution ici: https://www.cs.technion.ac.il/users/yechiel/c++-faq/separate-template-class-defn-from-decl.html

template <class A, class B> application<A,B>::application(application<A,B> * input_adr_def, cat_appli input_cat){
#ifdef COMMENTAIRES
    cout << "Création d'une application<" << boost::typeindex::type_id<A>().pretty_name()<< ","<<boost::typeindex::type_id<B>().pretty_name() << "> à l'adresse :  " << this << endl;
    cout << "Avec les paramètres : " << endl;
    cout << "cat " << cat << endl;
    cout << "adr_def " << adr_def << endl;
#endif
    cat = input_cat;
    adr_def = input_adr_def;
}

/*
template <class A, class B> application<A,B>::application(B y){
    cout << "conv implicite" << endl;
    constante<A,B> f(y);
    *this = f;
}
*/

template <class A, class B> application<A,B>::application(const application<A,B> & input_application){
    cat = input_application.cat;
    adr_def = input_application.adr_def;
}
//rerediger
template <class A, class B> application<A,B> & application<A,B>::operator= (const application<A,B> & input_application){
#ifdef COMMENTAIRES
    cout << "Appel de application<" << boost::typeindex::type_id<A>().pretty_name()<< ","<<boost::typeindex::type_id<B>().pretty_name() << "::operator= par" <<this << endl;
#endif
    if(this!= &input_application){
        cat = input_application.cat;
        adr_def = input_application.adr_def;
    }
    return *this;
}

template <class A, class B> bool application<A,B>::domaine_def(A x){
    //Par defaut tout R ou tout C
    return true;
}

/*---------------------------------------*/
template <class A,class C, class B> class compo_externe : public application<A,B>{
protected:
    pair<const application<A,C> *,const application<C,B>*> entrees;
public:
    
    //constructeur
    //Ce type d'application n'est destiné à être construit que lors d'une composition (le mettre en privé ?)
    compo_externe(pair<const application<A,C> *,const application<C,B>*>,cat_appli);
    //constructeur par recopie
    compo_externe(const compo_externe<A,C,B> &);
    //affectation
    compo_externe<A,C,B> & operator= (const compo_externe<A,C,B> &);
    
    //essayer ne pas mettre cette def inline
    //friend ostream & operator<< (ostream &,const application<A,B> &);
    inline virtual ostream & info(ostream & os) const {
        return os;
    }
    //Evaluation
    virtual B operator() (A) const; //fonction virtuelle et peut agir sur des objets constants
};


template <class A,class C, class B> compo_externe<A,C,B>::compo_externe(pair<const application<A,C> *,const application<C,B>*> input_entrees,cat_appli input_cat): application<A,B>(this,input_cat){
    entrees = input_entrees;
}

template <class A,class C, class B>  compo_externe<A,C,B>::compo_externe(const compo_externe<A,C,B> & input_application): application<A,B>(input_application){
    entrees = input_application.entrees;
}

template <class A,class C, class B> compo_externe<A,C,B> & compo_externe<A,C,B>::operator= (const compo_externe<A,C,B> & input_application){
    if(this!= &input_application){
        application<A,B> * ad1;
        const application<A,B> * ad2;
        ad1 = this;
        ad2 = &input_application;
        *ad1 = *ad2;
        entrees = input_application.entrees;
    }
    return *this;
}

//possible de ne définir que des specifications de cette fonction -> erreur de compiation si on ne l'utilise pas avec le bon type -> plus difficile à reperer
// Je choisi de la définir mais il y aura une erreur si le constructeur du type/class B requiert des paramètres
template <class A,class C,class B> B compo_externe<A,C,B>::operator() (A x) const {
    B result;
    cout << "Appel evaluation d'application" << endl;
    switch((*this).cat){
        case COMP : result = (*entrees.first)((*entrees.second)(x));
            return result;
    }
    cout << "Problème d'évaluation" << endl;
    B y = 1;
    return y;
}


/*----------------------------------------*/
template <class A,class B> class compo_interne : public application<A,B>{
protected:
    pair<const application<A,B> *,const application<A,B> *> entrees;
public:
    
    //constructeur
    //Ce type d'application n'est destiné à être construit que lors d'une composition (le mettre en privé ?)
    compo_interne(pair<const application<A,B> *,const application<A,B> *>,cat_appli);
    
    //constructeur par recopie
    compo_interne(const compo_interne<A,B> &);
    //affectation
    compo_interne<A,B> & operator= (const compo_interne<A,B> &);
    
    //essayer ne pas mettre cette def inline
    //friend ostream & operator<< (ostream &,const application<A,B> &);
    inline virtual ostream & info(ostream & os) const {
        return os;
    }
    //Evaluation
    virtual B operator() (A) const; //fonction virtuelle et peut agir sur des objets constants
};

template <class A, class B>  compo_interne<A,B>::compo_interne(pair<const application<A,B> *,const application<A,B> *> input_entrees,cat_appli input_cat): application<A,B>(this,input_cat) {
#ifdef COMMENTAIRES
    cout << "Création d'une compo_interne<" << boost::typeindex::type_id<A>().pretty_name()<< ","<<boost::typeindex::type_id<B>().pretty_name() << "> à l'adresse :  " << this << endl;
#endif
    entrees = input_entrees;
}

template <class A, class B>  compo_interne<A,B>::compo_interne(const compo_interne<A,B> & input_application): application<A,B>(input_application){
    entrees = input_application.entrees;
}

template <class A, class B>  compo_interne<A,B> & compo_interne<A,B>::operator= (const compo_interne<A,B> & input_application){
    if(this!= &input_application){
        application<A,B> * ad1;
        const application<A,B> * ad2;
        ad1 = this;
        ad2 = &input_application;
        *ad1 = *ad2;
        entrees = input_application.entrees;
    }
    return *this;
}

//possible de ne définir que des specifications de cette fonction -> erreur de compiation si on ne l'utilise pas avec le bon type -> plus difficile à reperer
// Je choisi de la définir mais il y aura une erreur si le constructeur du type/class B requiert des paramètres
template <class A, class B>  B compo_interne<A,B>::operator() (A x) const {
    B result;
    //cout << "catégorie"<< cat << endl;
    switch((*this).cat){
        case PLUS : result = (*entrees.first)(x)+ (*entrees.second)(x);
            return result;
        case MOINS : result = (*entrees.first)(x)- (*entrees.second)(x);
            return result;
        case MULT : result = (*entrees.first)(x)* (*entrees.second)(x);
            return result;
        case DIV :result = (*entrees.first)(x)/ (*entrees.second)(x);
            return result;
    }
    cout << "Problème d'évaluation" << endl;
    B y = 1;
    return y;
}



/*---------------------------------------*/
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
    
    inline virtual ostream & info(ostream & os) const  {
        os << "Application constante égale à "<< y <<endl;
        return os;
    }
    /*
    inline friend ostream & operator<<(ostream & os,const constante & f){
        //Boost a une fonction qui permet  de mieux afficher les paramètres de type
        //A tester si ça marche pour des types défini par l'utilisateur
        
        //os << ": "<< typeid(A).name() << "->"<< typeid(B).name() << endl << "cat : " << f.cat << endl << "entrees : (" << f.entrees.first << "," << f.entrees.second << ")";
        os << ": "<< boost::typeindex::type_id<A>().pretty_name()  << "->"<< boost::typeindex::type_id<B>().pretty_name()  << endl;
        os << "cat : BASE(" << f.cat << ")" << endl;
        os << "Application constante égale à "<< f.y <<endl;
        return os;
    }
    */
};
/*
template <class A,class B> constante<A,B>::constante(const string input_nom,B input_x): application<A,B>(input_nom){
    x = input_x;
}
 */
template <class A,class B> constante<A,B>::constante(B input_y): application<A,B>(this,BASE){
#ifdef COMMENTAIRES
    cout << "Création d'une constante<" << boost::typeindex::type_id<A>().pretty_name()<< ","<<boost::typeindex::type_id<B>().pretty_name() << "> à l'adresse :  " << this << endl;
#endif
    y = input_y;
}
template<class A,class B> B  constante<A,B>::operator() (A x) const{
    return y;
}


/*----------------------------------------*/
//mettre en parametre le type d'enveloppe

class paquet_d_onde: public application<double,complex<double> > {
    double mu_x; //moyenne en position
    double sigma_x;
    double mu_xi; //moyenne en freq
    public:
        //constructeur
        paquet_d_onde(double,double,double);
        //Evaluation
        complex<double>  operator() (double) const;
    
    
    inline virtual ostream & info(ostream & os) const  {
        os << "paquet_d_onde de parametres : "<<endl;
        os << "mu_x : "<< mu_x <<endl;
        os << "sigma_x : "<< sigma_x <<endl;
        os << "mu_xi : "<< mu_xi <<endl;
        return os;
    }
    /*
    inline friend ostream & operator<<(ostream & os,const paquet_d_onde & f){
        //Boost a une fonction qui permet  de mieux afficher les paramètres de type
        //A tester si ça marche pour des types défini par l'utilisateur
        
        //os << ": "<< typeid(A).name() << "->"<< typeid(B).name() << endl << "cat : " << f.cat << endl << "entrees : (" << f.entrees.first << "," << f.entrees.second << ")";
        os << ": "<< boost::typeindex::type_id<double>().pretty_name()  << "->"<< boost::typeindex::type_id<complex<double> >().pretty_name()  << endl;
        os << "cat : " << f.cat << endl;
        os << "paquet_d_onde de parametres : "<<endl;
        os << "mu_x : "<< f.mu_x <<endl;
        os << "sigma_x : "<< f.sigma_x <<endl;
        os << "mu_xi : "<< f.mu_xi <<endl;
        return os;
    }
     */
};

//Le constructeur de application<double,complex<double> > est appelé implicitement ?
paquet_d_onde::paquet_d_onde(double input_mu_x,double input_sigma_x,double input_mu_xi):application(this,BASE) {
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

#endif /* patron_de_classe_1_h */
