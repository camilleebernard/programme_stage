//
//  parties_du_main_2.hpp
//  programme_stage
//
//  Created by Camille Bernard on 01/08/2021.
//

//Pour afficher une description du programme.
#define COMMENTAIRES

#ifndef parties_du_main_2_hpp
#define parties_du_main_2_hpp

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

#define H_P_COOD 0  //Il faut mieux utiliser "const int H_P_COOD=5;" d'apres Delannoy chap 31
#define N_N_PRIM_COOD 1

//D'après SaeidMo7 sur:
// https://stackoverflow.com/questions/3386861/converting-a-variable-name-to-a-string-in-c
//Donne la chaîne de caractères associé au nom de la variable.
#define NOM(Variable) (#Variable)

// voir plus précisément comment (#f) fonctionne
#define INFO(f) cout <<"------------------------------------"<< endl << NOM(f) << f << endl \
<< "------------------------------------" << endl;
//backslash pour mettre le macro sur 2 lignes


void partie_0(void);
void partie_3(void);
void partie_5(void);
void partie_6(void);

complex<double> psi_0_X_psi_1_wp(double x_0, double p_0, double sigma_0,double x_1, double p_1, double sigma_1);
complex<double> psi_0_etX_psi_1_wp(double x_0, double p_0, double sigma_0,double x_1, double p_1, double sigma_1,double t);


//Pseudospectre
class resolvant {
private:
    arma::cx_mat A; //the resolvant of A will be (z -> (zI-A)^{-1}
public:
   //constructor
    resolvant(const arma::cx_mat &);
    
   //copy constructor
    resolvant(const resolvant &);
    
    void Update(const arma::cx_mat &);
    
    double norm(double *, double *);
};

bool pseudospectre_1(const arma::cx_mat *, const complex<double>,const float);
void affiche_pseudospectre_1(const arma::cx_mat *,double,double,double,double);


void affiche_vp(arma::cx_vec);
void affiche_vp(arma::vec eigval);


#endif /* parties_du_main_2_hpp */
