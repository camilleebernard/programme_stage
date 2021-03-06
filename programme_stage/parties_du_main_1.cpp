//
//  parties_du_main_1.cpp
//  programme_stage
//
//  Created by Camille Bernard on 22/07/2021.
//

#include "parties_du_main_1.hpp"
#include "patron_de_classe_1.h"

void partie_0(void){
    //Test pseudospectre de B = 1*E(1,1)+j*E(2,2)+ j^2*E(3,3)
    //Définition de A
    arma::Mat<arma::cx_double> B;
    B.eye(3,3);
    complex<double> I(0.,1.),z;
    double r,theta;
    r = 2.;
    theta = M_PI/6. ;
    z = r*exp(I*theta);
    B(0,0) = z;
    B(1,1) = z*exp(2.*I*M_PI/3.);
    B(2,2) = z*exp(4.*I*M_PI/3.);
    affiche_pseudospectre_1(&B,-4.,4.,-4.,4.);
}


void partie_3(void){
    //Essai paquet d'ondes
    
    //Définition de B
    int N = 200;  //40

    arma::Mat<arma::cx_double> B;
    B.zeros(2*N+1,2*N+1);
    
    
    
    double p_0 = 1;
    double lambda_0 = 2*M_PI/p_0;
    
    //espacement des paquets d'ondes a->espace / b-> fréquence
    double a = 0.1;
    double b = 0.1;
    
    complex<double> bnn_prim;
    
    for(int n= -N;n<=N;n++){
        for(int n_prim= -N;n_prim<=N;n_prim++){
            bnn_prim = psi_0_etX_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1,0);
            B(n+N,n_prim+N) = arma::cx_double(real(bnn_prim),imag(bnn_prim));
        }
    }
    //B.print();
    
    //diagonalise B
    arma::vec eigval_B;
    arma::cx_mat eigvec_B;
    
    cout << "B.is_hermitian()" << B.is_hermitian() << endl;
    arma::eig_sym(eigval_B, eigvec_B, B,"std");
    //affiche_vp(eigval_B);
    
    
    //Selection des directions avec vp assez grandes
    // Calcul des vp gardées
    double epsilon = 5;// pour a = b = 1 :1.0000795; // 1/epsilon  = 10^(-2)
    
    arma::cx_mat sqrt_B;
    cout << "min(vp): " << min(eigval_B) << endl;
    cout << "max(vp): " << max(eigval_B) << endl;
    cout << "epsilon: " << epsilon << endl;
    cout << "nb de directions propres gardées: " << size(find(eigval_B > epsilon))(0) << endl;
    sqrt_B = eigvec_B.cols(find(eigval_B > epsilon))*arma::diagmat(sqrt(eigval_B(find(eigval_B > epsilon))));
    //*(eigvec_B.cols(find(eigval_B > epsilon))).t();
    //affiche_vp(eigval_B);
    
    //Définition de A, t = 0
    arma::Mat<arma::cx_double> A;
    A.zeros(2*N+1,2*N+1);
    
    /*
    for(int n= -N;n<=N;n++){
        for(int n_prim= -N;n_prim<=N;n_prim++){
            bnn_prim = psi_0_etX_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1,0);
            A(n+N,n_prim+N) = arma::cx_double(real(bnn_prim),imag(bnn_prim));
        }
    }
    */
    //Définition de A_W
    arma::Mat<arma::cx_double> A_W;
    //A_W = sqrt_B.t()*A*sqrt_B; // .t() pour cx_mat = trans conj
    
    
    //Crée le GIF
    
    //Le temps va varier de sur [0,T] avec un pas de 1/NT
    const double minT = -5;
    const double maxT = 5;
    //const double zoom = 10;
    const double NT = 100;
    const double xmin = -50;
    const double xmax = 50;
    const double ymin = -50;
    const double ymax = 50;
    
    //nombre de paramètres dans TF2
    const int npar = 0;
    
    
    //Supprime l'ancien gif
    gSystem->Unlink("/Users/camillebernard/Desktop/Stage/Test_numériques/programme_stage/Gif/gif_pseudospectre1_flot.gif");
    
    auto c2 = new TCanvas("c2","Pseudospectre");
    
    //On construit un ojet resolvant qui contient la matrice A_W qu'on va mettre à jour après chaque incrément de temps.
    //Une méthode de cette classe permet de caculer la norme de la résolvante.
//https://root.cern.ch/doc/master/classTF2.html
    //https://root.cern/root/html528/TF1.html
    resolvant * adr_res = new resolvant(A_W);

    
    //Cercle r = 1
    double x[100], y[100];
    for(int i=0;i<100; i++){
        x[i] = cos(2*M_PI*((double)i/100.));
        y[i] = sin(2*M_PI*((double)i/100.));
    }
    
    TGraph *cercle = new TGraph(100,x,y);
    cercle->SetName("cercle");
    cercle->SetLineColor(2);
    
    
    for(int nt = 0;nt<=NT;nt++){
        //Update A
        A.zeros(2*N+1,2*N+1);
        for(int n= -N;n<=N;n++){
            for(int n_prim= -N;n_prim<=N;n_prim++){
                bnn_prim = psi_0_etX_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1,minT+ nt*(maxT-minT)/NT);
                //bnn_prim = psi_0_etX_psi_1(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1,fct_zoom(nt/NT, zoom, minT, maxT));
                A(n+N,n_prim+N) = arma::cx_double(real(bnn_prim),imag(bnn_prim));
            }
        }
        
        //Update A_W
        A_W = sqrt_B.t()*A*sqrt_B;
        
        
        //Calcule le spectre de A_W
        arma::cx_vec eigval_A_W;
        arma::eig_gen(eigval_A_W, A_W);
        
        //Trace le spectre de A_W
        TGraph *spectre = new TGraph();
        spectre -> SetMarkerStyle(7);
        spectre -> SetMarkerColor(2);
        spectre -> SetMarkerSize(15);
        for(int i= 0;i<size(eigval_A_W)(0);i++){
            spectre->SetPoint(i,eigval_A_W(i).real(),eigval_A_W(i).imag());
        }
        
        
        //Pseudospectre
        (*adr_res).Update(A_W);
        //Pour construire un obj TF2 avec une classe définie préalablement https://root.cern.ch/doc/master/classTF2.html
        //https://root.cern/root/html528/TF1.html
        TF2 *f = new TF2("Pseudospectre",adr_res,&resolvant::norm,xmin,xmax,ymin,ymax,npar,2);
        f -> Draw("CONT1"); // doc: https://root.cern.ch/doc/master/draw2dopt_8C.html
        //CONT1 SURF1
        
        //Titre
        f-> SetTitle(Form("Pseudospectre t = %f",minT+ nt*(maxT-minT)/NT));
        
        //Ajoute le cercle unité
        cercle->Draw("L");
        //Ajoute le spectre
        spectre->Draw("P");
        
        //Ajouter le dessin au Gif
        c2->Print("/Users/camillebernard/Desktop/Stage/Test_numériques/programme_stage/Gif/gif_pseudospectre1_flot.gif+");
        //Suppression de l'objet TF2
        delete f;
        //Suppression de l'objet TGraph
        delete spectre;
    }
}

void partie_4(void){
    /*
    application<int,double>  f("f");
    cout << f(10) <<endl;
    application<complex<double>,complex<double> >  g("g");
    complex<double> I(0,1);
    cout << g(I) <<endl;
    
    pair<double,double> a(45.,67.);
    cout << a.first <<endl;
    cout << a.second <<endl;
    */
    complex<double> I(0,1);
    un<complex<double>,complex<double> > h("h");
    cout << "h(I)"<< h(I) << endl;
    un<complex<double>,complex<double> > i("i");
    cout << "i(I)"<< i(I) << endl;
    un<complex<double>,complex<double> > k("k");
    cout << "k(I)"<< k(3.*I) << endl;
    application<complex<double>,complex<double> > j("j");
    j = h+i+k;
    cout << "j(I)"<< j(I) << endl;
    application<complex<double>,complex<double> > l("l");
    l = k-j;
    cout << "l(I)"<< l(I) << endl;
    //test paquet d'onde
    paquet_d_onde psi_1("psi_1",1.,2.,3.);
    cout << "psi_1(7)"<< psi_1(7.) <<endl;
    paquet_d_onde psi_2("psi_2",8.,4.,1.);
    cout << "psi_2(7)"<< psi_2(7.) <<endl;
    application<double,complex<double> > n("n");
    n = psi_1+psi_2;
    cout << "n(7)"<< n(7.) <<endl;
    cout << "psi_1(7)+psi_2(7)"<< psi_1(7.)+psi_2(7.) <<endl;
    
    //test conversion constante->fct constante
    constante<double, double> m(5.);
    cout << "m(10.)" <<m(10.) << endl;
    application<double, double> o(m);
    cout << "m(10.)" <<o(10.) << endl;
    //test fonction virtuelle
    application<double, double> * adr_appli_1;
    adr_appli_1 = & m;
    cout << "(*adr_appli_1)(10.)" << (*adr_appli_1)(10.) << endl;
    const application<double, double> * adr_appli_2 = &m;
    cout << "(*adr_appli_2)(10.)" << (*adr_appli_2)(10.) << endl;
    
    
    //test mult sca*f
    application<complex<double>, complex<double> > p("p");
    p = I*l;
    cout << "p(I)" << p(I)<< endl;
    
    /*
    //test mult sca*f
    cout << "j(I) avant" << j(I)<< endl;
    j = I*l;
    cout << "j(I) après" << j(I)<< endl;
    */
}
void partie_5(void){
    application<double,double> a;
    paquet_d_onde b(1.,2.,3.,6.);    
    
}


complex<double> psi_0_X_psi_1_wp(double x_0, double p_0, double sigma_0,double x_1, double p_1, double sigma_1){
    complex<double> result, I(0.,1.),alpha,A,C;
    double a,b,c;
    a= (pow(sigma_0,2)+pow(sigma_1,2))/(pow(sigma_0,2)*pow(sigma_1,2));
    b = (pow(sigma_0,2)*x_1 + pow(sigma_1,2)*x_0)/(pow(sigma_0,2)+pow(sigma_1,2));
    c = pow(x_0-x_1,2)/(pow(sigma_0,2)+pow(sigma_1,2));
    alpha = x_1 + I * pow(sigma_1,2)*p_1;
    A = -c/2. - I*(p_1*x_1-p_0*x_0)- pow(p_1-p_0,2)/(2.*a)+ 2. * I* (p_1-p_0)*b;
    C = 1./(pow(sigma_1,2)*sqrt(sigma_0*sigma_1))*sqrt(2)/pow(a,2)*(1.+a*b*(b-alpha)+(p_1-p_0)*(I*(2.*b-alpha)-(p_1-p_0)/a));
    result = C*exp(A);
    return result;
}

complex<double> psi_0_etX_psi_1_wp(double x_0, double p_0, double sigma_0,double x_1, double p_1, double sigma_1,double t){
    complex<double> result, I(0.,1.);
    double a,lambda;
    lambda = exp(-t);
    a = (pow(sigma_1,2)+pow(sigma_0,2)*pow(lambda,2))/(pow(sigma_0,2)*pow(sigma_1,2));
    result = (1/(pow(sigma_0,2)*pow(sigma_1,2)))*sqrt(2/a)*exp(-pow(p_1*lambda-p_0,2)/(2*a)-pow(x_0*lambda-x_1,2)/(2*pow(sigma_0,2)*pow(sigma_1,2)*a)+I*((x_0*pow(sigma_1,2)+x_1*pow(sigma_0,2)*lambda)*(p_1*exp(-t)-p_0)/(pow(sigma_0,2)*pow(sigma_1,2)*a)-(p_1*x_1-p_0*x_0)));
    return result;
}



//Pseudospectre

//arma::cx_mat = arma::Mat<arma::cx_double>
bool pseudospectre_1(const arma::cx_mat * adr_matrice, const complex<double> z, const float epsilon){
    // A = zI-(*adr_matrice)
    arma::cx_mat A;
    A.eye(size(*adr_matrice));
    A = z*A-(*adr_matrice);
    
    arma::cx_mat U;
    arma::vec s;
    arma::cx_mat V;
    
    arma::svd(U,s,V,A);
    cout << "arma::min(s)" << arma::min(s) << endl;
    
    if(arma::min(s)<epsilon){
        return true;
    }else{
        return false;
    }
}

//La fonction vient de la page: https://root.cern.ch/doc/master/annotation3d_8C.html
Double_t func(Double_t *val, Double_t *par)
{
   Float_t x = val[0];
   Float_t y = val[1];
   Double_t f = x*x-y*y;
   return f;
}

void affiche_pseudospectre_1(const arma::cx_mat * adr_matrice, double x_min, double x_max, double y_min, double y_max){
    
    
    TApplication theApp("App",nullptr,nullptr);
    
    //test avec la fonction func
    //TF2 *f = new TF2("f",func,-1,1,-1,1);
    
    //Pour construire un obj TF2 avec une classe définie préalablement https://root.cern.ch/doc/master/classTF2.html
    //https://root.cern/root/html528/TF1.html
    resolvant * adr_res = new resolvant(*adr_matrice);
    int npar = 0;
    TF2 *f = new TF2("Pseudospectre",adr_res,&resolvant::norm,x_min,x_max,y_min,y_max,npar,2);
    f -> Draw("CONT1"); // doc: https://root.cern.ch/doc/master/draw2dopt_8C.html
    
    //Cercle r = 1
    double x[1000], y[1000];
    for(int i=0;i<1000; i++){
        x[i] = cos(2*M_PI*((double)i/1000.));
        y[i] = sin(2*M_PI*((double)i/1000.));
    }
    
    TGraph *cercle = new TGraph(1000,x,y);
    cercle->SetName("cercle");
    cercle->SetLineColor(2);
    //cercle->Fit(f);
    cercle->Draw("same");
    
    theApp.Run();
}

//Class resolvant

resolvant::resolvant(const arma::cx_mat & input_A){
    A = input_A;
}
resolvant::resolvant(const resolvant & input_R){
    A = input_R.A;
}

double resolvant::norm(double * val, double * par){
    arma::cx_mat B;
    B.eye(size(A));
    B = arma::cx_double(val[0],val[1])*B-A;
    
    arma::vec s = arma::svd(B);
    //cout << "arma::min(s)" << arma::min(s) << endl;
    return ((double)arma::min(s));
}

void resolvant::Update(const arma::cx_mat & input_A){
    A = input_A;
}


void affiche_spectre(arma::cx_vec eigval){
    int N = arma::size(eigval)(0);
    //cout <<"N = " << N <<endl;
    
    auto c2 = new TCanvas("c2","Resonances de Ruelles");
    c2 -> SetCanvasSize(700, 700);
    c2 -> SetWindowSize(700, 700);
    c2 -> SetGrid(1,1); //set the grid on both axis https://root.cern/doc/v612/classTGaxis.html
    
    
    double x[1000], y[1000];
    for(int i=0;i<1000; i++){
        x[i] = cos(2*M_PI*((double)i/1000.));
        y[i] = sin(2*M_PI*((double)i/1000.));
    }
    
    TGraph *cercle = new TGraph(1000,x,y);
    cercle->SetName("cercle");
    cercle->SetTitle("Norme de l'opérateur");
    cercle->SetLineColor(3);
    cercle->SetLineWidth(1);
    
    TGraph * spectre = new TGraph();
    spectre->SetName("spectre");
    spectre-> SetTitle("Spetre de l'opérateur");
    spectre -> SetMarkerStyle(7);
    spectre -> SetMarkerColor(2);
    spectre -> SetMarkerSize(15);
    spectre ->GetXaxis()->SetRange(-1.1,1.1);
    spectre ->GetYaxis()->SetRange(-1.1,1.1);
    
    //Crée un multigraphe avec les deux graph précédents
    TMultiGraph * mg = new TMultiGraph();
    mg -> SetTitle("Spectre de A_W;\\mathcal{Re};\\mathcal{Im}");
    mg -> Add(cercle,"L");
    mg -> Add(spectre,"P");
    
    gSystem->Unlink("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif"); //Suprime l'ancien gif
    
    //Calcul des valeurs et vecteurs propres

    const int kUPDATE = 2;
    for (int nt = 2;nt<=N;nt++){
        int i = 0;
        while(i<nt){
            //cout << "nt = " << nt <<endl;
            //cout << "eigenval i" << eigval(i) <<endl;
            spectre->SetPoint(i,eigval(i).real(),eigval(i).imag());
            i++;
        }
       if (nt && (nt%kUPDATE) == 0) {
           if (nt == kUPDATE) {
            mg->Draw("A"); //par dessus l'autre
            c2->Update();
          }
           c2->Modified();
           c2->Update();
           if (gROOT->IsBatch()) {
               c2->Print("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif+");
               //printf("nt = %d\n", nt);
           } else {
               if (gSystem->ProcessEvents()) break;
          }
       }
    }
    // make infinite animation by adding "++" to the file name
    if (gROOT->IsBatch()) c2->Print("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif++");
}


void affiche_vp(arma::cx_vec eigval){
    int N = arma::size(eigval)(0);
    //cout <<"N = " << N <<endl;
    
    auto c2 = new TCanvas("c2","Resonances de Ruelles");
    c2 -> SetCanvasSize(700, 700);
    c2 -> SetWindowSize(700, 700);
    c2 -> SetGrid(1,1); //set the grid on both axis https://root.cern/doc/v612/classTGaxis.html
    
    
    double x[1000], y[1000];
    for(int i=0;i<1000; i++){
        x[i] = cos(2*M_PI*((double)i/1000.));
        y[i] = sin(2*M_PI*((double)i/1000.));
    }
    
    TGraph *cercle = new TGraph(1000,x,y);
    cercle->SetName("cercle");
    cercle->SetTitle("Norme de l'opérateur");
    cercle->SetLineColor(3);
    cercle->SetLineWidth(1);
    
    TGraph * spectre = new TGraph();
    spectre->SetName("spectre");
    spectre-> SetTitle("Spetre de l'opérateur");
    spectre -> SetMarkerStyle(7);
    spectre -> SetMarkerColor(2);
    spectre -> SetMarkerSize(15);
    spectre ->GetXaxis()->SetRange(-1.1,1.1);
    spectre ->GetYaxis()->SetRange(-1.1,1.1);
    
    //Crée un multigraphe avec les deux graph précédents
    TMultiGraph * mg = new TMultiGraph();
    mg -> SetTitle("Spectre de A_W;\\mathcal{Re};\\mathcal{Im}");
    mg -> Add(cercle,"L");
    mg -> Add(spectre,"P");
    
    gSystem->Unlink("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif"); //Suprime l'ancien gif
    
    //Calcul des valeurs et vecteurs propres

    const int kUPDATE = 2;
    for (int nt = 2;nt<=N;nt++){
        int i = 0;
        while(i<nt){
            //cout << "nt = " << nt <<endl;
            //cout << "eigenval i" << eigval(i) <<endl;
            spectre->SetPoint(i,eigval(i).real(),eigval(i).imag());
            i++;
        }
       if (nt && (nt%kUPDATE) == 0) {
           if (nt == kUPDATE) {
            mg->Draw("A"); //par dessus l'autre
            c2->Update();
          }
           c2->Modified();
           c2->Update();
           if (gROOT->IsBatch()) {
               c2->Print("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif+");
               //printf("nt = %d\n", nt);
           } else {
               if (gSystem->ProcessEvents()) break;
          }
       }
    }
    // make infinite animation by adding "++" to the file name
    if (gROOT->IsBatch()) c2->Print("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif++");
}


void affiche_vp(arma::vec eigval){
    int N = arma::size(eigval)(0);
    cout <<"N = " << N <<endl;
    
    auto c2 = new TCanvas("c2","Resonances de Ruelles");
    c2 -> SetCanvasSize(700, 700);
    c2 -> SetWindowSize(700, 700);
    c2 -> SetGrid(1,1); //set the grid on both axis https://root.cern/doc/v612/classTGaxis.html
    
    
    double x[1000], y[1000];
    for(int i=0;i<1000; i++){
        x[i] = cos(2*M_PI*((double)i/1000.));
        y[i] = sin(2*M_PI*((double)i/1000.));
    }
    
    TGraph *cercle = new TGraph(1000,x,y);
    cercle->SetName("cercle");
    cercle->SetTitle("Norme de l'opérateur");
    cercle->SetLineColor(3);
    cercle->SetLineWidth(1);
    
    TGraph * spectre = new TGraph();
    spectre->SetName("spectre");
    spectre-> SetTitle("Spetre de l'opérateur");
    spectre -> SetMarkerStyle(7);
    spectre -> SetMarkerColor(2);
    spectre -> SetMarkerSize(15);
    spectre ->GetXaxis()->SetRange(-1.1,1.1);
    spectre ->GetYaxis()->SetRange(-1.1,1.1);
    
    //Crée un multigraphe avec les deux graph précédents
    TMultiGraph * mg = new TMultiGraph();
    mg -> SetTitle("Spectre de A_W;\\mathcal{Re};\\mathcal{Im}");
    mg -> Add(cercle,"L");
    mg -> Add(spectre,"P");
    
    gSystem->Unlink("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif"); //Suprime l'ancien gif
    
    //Calcul des valeurs et vecteurs propres

    const int kUPDATE = 1; // on crée une image tout les kUPDATE changement
    for (int nt = 1;nt<=N;nt++){
        int i = 0;
        while(i<nt){
            //cout << "nt = " << nt <<endl;
            //cout << "eigenval i" << eigval(i) <<endl;
            spectre->SetPoint(i,eigval(i),0);
            i++;
        }
       if (nt && (nt%kUPDATE) == 0) {
           cout << "nt%kUPDATE" << nt%kUPDATE << endl;
           if (nt == kUPDATE) { //première fois qu'on affiche l'image
            mg->Draw("A"); //par dessus l'autre
            c2->Update();
          }
           c2->Modified();
           c2->Update();
           if (gROOT->IsBatch()) {
               c2->Print("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif+");
               //printf("nt = %d\n", nt);
           } else {
               if (gSystem->ProcessEvents()) break;
          }
       }
    }
    // make infinite animation by adding "++" to the file name
    if (gROOT->IsBatch()){ c2->Print("/Users/camillebernard/Desktop/Stage/Test_numériques/new_prog1/Gif/gif_vp.gif++");
    }
}
