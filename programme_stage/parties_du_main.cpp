//
//  parties_du_main.cpp
//  programme_stage
//
//  Created by Camille Bernard on 08/07/2021.
//

#include "parties_du_main.hpp"


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
void partie_1(void){
    //Essai paquet d'ondes
    
    //Définition de B
    int N = 40;
    arma::Mat<arma::cx_double> B;
    B.zeros(2*N+1,2*N+1);
    //B.zeros(N,N);
    double p_0 = 1;
    double lambda_0 = 2*M_PI/p_0;
    
    
    //espacement des paquets d'ondes a->espace / b-> fréquence
    double a = 1;
    double b = 1;
    
    complex<double> bnn_prim;
    
    for(int n= -N;n<=N;n++){
        for(int n_prim= -N;n_prim<=N;n_prim++){
            bnn_prim = psi_0_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1);
            B(n+N,n_prim+N) = arma::cx_double(real(bnn_prim),imag(bnn_prim));
        }
    }
    
    //B.print();
    
    //Définition de A
    
    arma::Mat<arma::cx_double> A;
    A.zeros(2*N+1,2*N+1);
    //A.zeros(N,N);
    
    
    for(int n= -N;n<=N;n++){
        for(int n_prim= -N;n_prim<=N;n_prim++){
            bnn_prim = psi_0_X_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1);
            A(n+N,n_prim+N) = arma::cx_double(real(bnn_prim),imag(bnn_prim));
        }
    }
    
    //A.print();
    
    //Diagonalise B
    arma::vec eigval_B;
    arma::cx_mat eigvec_B;
    
    
    cout << "B.is_hermitian()" << B.is_hermitian() << endl;
    arma::eig_sym(eigval_B, eigvec_B, B,"std");
    //affiche_vp(eigval_B);
    
    //Selection des directions avec vp assez grandes
    
    // Calcul des vp gardées
    double epsilon = 0.99; // 1/epsilon  = 10^(-2) 0.99
    int M = 0; //nb de vp gardées
    
    for (int i = -N;i<=N;i++){
        cout << "eigval_B(i+N)" << i+N <<"   " <<eigval_B(i+N) << endl;
        if(eigval_B(i+N)>epsilon){
            cout << "retenue" << endl;
            M++;
        }
    }
    
    cout << " nb_eigval_gardee " << M << endl;
    
    //Retient les indices ok
    int ind[M];
    int ind1 = 0;
    for (int i = -N;i<=N;i++){
        if(eigval_B(i+N)>epsilon){
            ind[ind1] =i+N;
            cout << "ind " << ind1 <<" : " << ind[ind1] << endl;
            ind1++;
        }
    }
    
    //P_s
    arma::cx_mat P_s(2*N+1,M);
    for(int i = -N; i<=N;i++){
        for(int j = 0; j<=M-1;j++){
            P_s(i+N,j) = eigvec_B(i+N,ind[j]);
        }
    }
    
    
    //Pour tester
    //arma::cx_mat D(M,M);
    //D = P_s.t()*B*P_s;
    //cout << "D" << D << endl;
    
    
    arma::cx_mat D;
    D.zeros(M,M);
    for(int i = 0; i<=M-1;i++){
        D(i,i) = eigval_B(ind[i]);
    }
    //cout << "D" << D << endl;
    
    arma::cx_mat sqrt_B;
    sqrt_B = P_s*D*P_s.t();  // .t() pour cx_mat = trans conj
    cout << "sqrt_B" << size(sqrt_B) << endl;
    
    
    //A_W de taille M*M
    arma::Mat<arma::cx_double> A_W;
    
    /*
    A_W.zeros(M,M);
    for (int j = 0;j <= M-1;j++){
        for (int k = 0;k <= M-1;k++){
            
            for (int i = -N;i<=N;i++){
                for (int l = -N;l<=N;l++){
                    A_W(j,k) += 1./sqrt(eigval_B(ind[j])*eigval_B(ind[k]))*conj(eigvec_B(i+N,ind[j]))*eigvec_B(l+N,ind[k])*A(i+N,l+N);
                }
            }
        }
    }
    */
    
    A_W = (sqrt_B.i()).t()*A*sqrt_B.i(); //.i -> inverse
    arma::cx_vec eigval_A_W = arma::eig_gen(A_W);
    affiche_pseudospectre_1(&A_W,-1.,1.,-1.,1.);
    //affiche_vp(eigval_A_W);
    //cout << C << endl;
    
    //tuto_gif_vp_1();
    
}


//ancienne version inutile
void partie_2(void){
    //Essai paquet d'ondes
    
    //Définition de B
    int N = 40;
    arma::Mat<arma::cx_double> B;
    B.zeros(2*N+1,2*N+1);
    //B.zeros(N,N);
    double p_0 = 1;
    double lambda_0 = 2*M_PI/p_0;
    
    //espacement des paquets d'ondes a->espace / b-> fréquence
    double a = 1;
    double b = 1;
    
    complex<double> bnn_prim;
    
    for(int n= -N;n<=N;n++){
        for(int n_prim= -N;n_prim<=N;n_prim++){
            bnn_prim = psi_0_etX_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1,0);
            //bnn_prim = psi_0_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1);
            B(n+N,n_prim+N) = arma::cx_double(real(bnn_prim),imag(bnn_prim));
        }
    }
    //B.print();
    
    //Diagonalise B
    arma::vec eigval_B;
    arma::cx_mat eigvec_B;
    
    cout << "B.is_hermitian()" << B.is_hermitian() << endl;
    arma::eig_sym(eigval_B, eigvec_B, B,"std");
    //affiche_vp(eigval_B);
    
    //Selection des directions avec vp assez grandes
    
    // Calcul des vp gardées
    double epsilon = 0.99; // 1/epsilon  = 10^(-2)
    int M = 0; //nb de vp gardées
    
    for (int i = -N;i<=N;i++){
        //cout << "eigval_B(i+N)" << i+N <<"   " <<eigval_B(i+N) << endl;
        if(eigval_B(i+N)>epsilon){
            //cout << "retenue" << endl;
            M++;
        }
    }
    cout << " nb_eigval_gardee " << M << endl;
    
    //Retient les indices ok
    int ind[M];
    int ind1 = 0;
    for (int i = -N;i<=N;i++){
        if(eigval_B(i+N)>epsilon){
            ind[ind1] =i+N;
            //cout << "ind " << ind1 <<" : " << ind[ind1] << endl;
            ind1++;
        }
    }
    
    //P_s
    arma::cx_mat P_s(2*N+1,M);
    for(int i = -N; i<=N;i++){
        for(int j = 0; j<=M-1;j++){
            P_s(i+N,j) = eigvec_B(i+N,ind[j]);
        }
    }
    
    
    //Pour tester
    //arma::cx_mat D(M,M);
    //D = P_s.t()*B*P_s;
    //cout << "D" << D << endl;
    
    
    arma::cx_mat D;
    D.zeros(M,M);
    for(int i = 0; i<=M-1;i++){
        D(i,i) = eigval_B(ind[i]);
    }
    //cout << "D" << D << endl;
    
    arma::cx_mat sqrt_B;
    sqrt_B = P_s*sqrt(D)*P_s.t();  // .t() pour cx_mat = trans conj
    //cout << "sqrt_B" << size(sqrt_B) << endl;
    
    //Définition de A, t = 0
    
    arma::Mat<arma::cx_double> A;
    A.zeros(2*N+1,2*N+1);
    //A.zeros(N,N);
    
    
    for(int n= -N;n<=N;n++){
        for(int n_prim= -N;n_prim<=N;n_prim++){
            bnn_prim = psi_0_etX_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1,0);
            A(n+N,n_prim+N) = arma::cx_double(real(bnn_prim),imag(bnn_prim));
        }
    }
    
    //A.print();
    
    //A_W de taille M*M
    arma::Mat<arma::cx_double> A_W;
    
    /*
    A_W.zeros(M,M);
    for (int j = 0;j <= M-1;j++){
        for (int k = 0;k <= M-1;k++){
            
            for (int i = -N;i<=N;i++){
                for (int l = -N;l<=N;l++){
                    A_W(j,k) += 1./sqrt(eigval_B(ind[j])*eigval_B(ind[k]))*conj(eigvec_B(i+N,ind[j]))*eigvec_B(l+N,ind[k])*A(i+N,l+N);
                }
            }
        }
    }
    */
    
    A_W = (sqrt_B.i()).t()*A*sqrt_B.i(); //.i -> inverse
    //arma::cx_vec eigval_A_W = arma::eig_gen(A_W);
    //affiche_pseudospectre_1(&A_W);
    
    //Crée le GIF
    
    //Le temps va varier de sur [0,T] avec un pas de 1/NT
    const double minT = -5;
    const double maxT = 5;
    //const double zoom = 10;
    const double NT = 50;
    const double xmin = -10;
    const double xmax = 10;
    const double ymin = -10;
    const double ymax = 10;
    
    //nombre de paramètres dans TF2
    const int npar = 0;
    
    //Supprime l'ancien gif
    gSystem->Unlink("/Users/camillebernard/Desktop/Stage/Test_numériques/programme_stage/Gif/gif_pseudospectre1_flot.gif");
    
    //Update A
    A.zeros(2*N+1,2*N+1);
    for(int n= -N;n<=N;n++){
        for(int n_prim= -N;n_prim<=N;n_prim++){
            bnn_prim = psi_0_etX_psi_1_wp(n*a*lambda_0,n*b,1,n_prim*a*lambda_0,n_prim*b,1,0);
            A(n+N,n_prim+N) = arma::cx_double(real(bnn_prim),imag(bnn_prim));
        }
    }
    
    //Update A_W
    A_W = (sqrt_B.i()).t()*A*sqrt_B.i();
    
    auto c2 = new TCanvas("c2","Pseudospectre");
    
    //Pour construire un obj TF2 avec une classe définie préalablement https://root.cern.ch/doc/master/classTF2.html
    //https://root.cern/root/html528/TF1.html
    
    resolvant * adr_res = new resolvant(A_W);
    //TF2 *f = new TF2("Pseudospectre",adr_res,&resolvant::norm,xmin,xmax,ymin,ymax,npar,2);
    //f -> Draw("CONT1"); // doc: https://root.cern.ch/doc/master/draw2dopt_8C.html
    
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
    //cercle->Draw("same");
    
    //c2->Print("/Users/camillebernard/Desktop/Stage/Test_numériques/programme_stage/Gif/gif_pseudospectre1_flot.gif+");
    
    //delete f;
    //delete adr_res;
    
    for (int nt = 0;nt<=NT;nt++){
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
        A_W = (sqrt_B.i()).t()*A*sqrt_B.i();
        
        //Pour construire un obj TF2 avec une classe définie préalablement https://root.cern.ch/doc/master/classTF2.html
        //https://root.cern/root/html528/TF1.html
        (*adr_res).Update(A_W);
        TF2 *f = new TF2("Pseudospectre",adr_res,&resolvant::norm,xmin,xmax,ymin,ymax,npar,2);
        f -> Draw("CONT1"); // doc: https://root.cern.ch/doc/master/draw2dopt_8C.html
        //CONT1 SURF1
        
        f-> SetTitle(Form("Pseudospectre t = %f",minT+ nt*(maxT-minT)/NT));
        cercle->Draw("same");
        
        
        c2 -> Print("/Users/camillebernard/Desktop/Stage/Test_numériques/programme_stage/Gif/gif_pseudospectre1_flot.gif+");
        
        delete f;
    }
    
    //TF2 *f = new TF2("Pseudospectre",adr_res,&resolvant::norm,-1,1,-1,1,npar,2);
    //f -> Draw("CONT1"); // doc: https://root.cern.ch/doc/master/draw2dopt_8C.html
    //cercle->Draw("same");
    
    //c2 -> Print("/Users/camillebernard/Desktop/Stage/Test_numériques/programme_stage/Gif/gif_pseudospectre1_flot.gif+");
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


//Pour le calcul d'elements de matrice
complex<double> psi_0_psi_1_wp(double x_0, double p_0, double sigma_0,double x_1, double p_1, double sigma_1){
    complex<double> result, I(0.,1.);
    //complex<double> alpha;
    double a,b,c;
    //double b;
    a= (pow(sigma_0,2)+pow(sigma_1,2))/(pow(sigma_0,2)*pow(sigma_1,2));
    b = (pow(sigma_0,2)*x_1 + pow(sigma_1,2)*x_0)/(pow(sigma_0,2)+pow(sigma_1,2));
    c = pow(x_0-x_1,2)/(pow(sigma_0,2)+pow(sigma_1,2));
    //b = (pow(sigma_0,2)*x_1 + pow(sigma_1,2)*x_0)/(pow(sigma_0,2)+pow(sigma_1,2));
    //alpha = x_1 + I * pow(sigma_1,2)*p_1;
    
    //result = sqrt(2./(sigma_0*sigma_1*a))*exp(-(pow(x_0-x_1,2)+pow(sigma_0,2)*pow(sigma_1,2)*pow(p_1-p_0,2))/(2*(pow(sigma_0,2)+pow(sigma_1,2))))*exp(I*(pow(sigma_0,2)*p_0+pow(sigma_1,2)*p_1)/(pow(sigma_0,2)+pow(sigma_1,2))*(x_0-x_1));
    result = sqrt(2./(sigma_0*sigma_1*a))*exp(-pow(p_1-p_0,2)/(2*a)-c/2.+I*(p_1-p_0)*b-I*(p_1*x_1-p_0*x_0));
    return result;
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


void reconditionne(arma::mat * adr_M,double epsilon){
    
    //svd
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U,s,V,*adr_M);
    
    cout << "Avant reconditionnement : " << endl;
    
    cout << "* adr_M" << * adr_M <<endl;
    
    cout << "s_min : " << arma::min(s) << endl;
    cout << "s_max : " << arma::max(s) << endl;
    cout << "kappa : " << arma::max(s)/arma::min(s) << endl;
    
    cout << "*adr_M" << (*adr_M) <<endl;
    cout << "U*arma::diagmat(s)*V.t()" << U*arma::diagmat(s)*V.t() <<endl;
    
    arma::vec eigval_M;
    arma::mat eigvec_M;
    
    cout << "M.is_hermitian()" << (*adr_M).is_hermitian() << endl;
    arma::eig_sym(eigval_M, eigvec_M, *adr_M,"std");
    cout << eigval_M << endl;
    cout << find(eigval_M < epsilon) << endl;
    
    
    //Selection des directions avec vp assez grandes
    *adr_M = eigvec_M.cols(find(eigval_M < epsilon))*arma::diagmat(eigval_M(find(eigval_M < epsilon)))*(eigvec_M.cols(find(eigval_M < epsilon))).t();  // .t() pour cx_mat = trans conj
    
    
    s = arma::svd(*adr_M);
    cout << "Après reconditionnement : " << endl;
    
    cout << "*adr_M" << (*adr_M) << endl;
    
    cout << "s_min : " << arma::min(s) <<endl;
    cout << "s_max : " << arma::max(s) <<endl;
    cout << "kappa : " << arma::max(s)/arma::min(s) << endl;
    
}


void reconditionne(arma::cx_mat * adr_M,double epsilon){
    //svd
    arma::cx_mat U;
    arma::vec s;
    arma::cx_mat V;
    arma::svd(U,s,V,*adr_M);
    
    cout << "Avant reconditionnement : " << endl;
    cout << "s_min : " << arma::min(s) << endl;
    cout << "s_max : " << arma::max(s) << endl;
    cout << "kappa : " << arma::max(s)/arma::min(s) << endl;

    
    cout << "*adr_M" << *adr_M <<endl;
    cout << "U*arma::diagmat(s)*V.t()" << U*arma::diagmat(s)*V.t() <<endl;
    
    arma::mat D = arma::diagmat(s);
    
    //Supprime les directions avec des vp trop petites
    for (int i=0;i<size(D)(0);i++){
        cout << "s("<<i<<")"<< s(i) << endl;
        if(s(i)<epsilon){
            cout << "trop petite" << endl;
            D(i,i) = 0;
        }
    }
    
    *adr_M = U*D*V.t();
    
    s = arma::svd(*adr_M);
    
    cout << "Après reconditionnement : " << endl;
    cout << "s_min : " << arma::min(s) <<endl;
    cout << "s_max : " << arma::max(s) <<endl;
    cout << "kappa : " << arma::max(s)/arma::min(s) << endl;
    
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





void test_graph_surface(void){
    TApplication theApp("App",nullptr,nullptr);
    
    //Le code vient de la page: https://root.cern.ch/doc/master/annotation3d_8C.html
    
    TCanvas *c = new TCanvas("c", "c", 600, 600);
    c->SetTheta(30);
    c->SetPhi(50);
    gStyle->SetOptStat(0);
    gStyle->SetHistTopMargin(0);
    gStyle->SetOptTitle(kFALSE);
    
    // Define and draw a surface
    TF2 *f = new TF2("f", "[0]*cos(x)*cos(y)", -1, 1, -1, 1);
    f->SetParameter(0, 1);
    double s = 1./f->Integral(-1, 1, -1, 1);
    f->SetParameter(0, s);
    f->SetNpx(50);
    f->SetNpy(50);
    
    f->GetXaxis()->SetTitle("x");
    f->GetXaxis()->SetTitleOffset(1.4);
    f->GetXaxis()->SetTitleSize(0.04);
    f->GetXaxis()->CenterTitle();
    f->GetXaxis()->SetNdivisions(505);
    f->GetXaxis()->SetTitleOffset(1.3);
    f->GetXaxis()->SetLabelSize(0.03);
    f->GetXaxis()->ChangeLabel(2,-1,-1,-1,kRed,-1,"X_{0}");
    
    f->GetYaxis()->SetTitle("y");
    f->GetYaxis()->CenterTitle();
    f->GetYaxis()->SetTitleOffset(1.4);
    f->GetYaxis()->SetTitleSize(0.04);
    f->GetYaxis()->SetTitleOffset(1.3);
    f->GetYaxis()->SetNdivisions(505);
    f->GetYaxis()->SetLabelSize(0.03);
    
    f->GetZaxis()->SetTitle("dP/dx");
    f->GetZaxis()->CenterTitle();
    f->GetZaxis()->SetTitleOffset(1.3);
    f->GetZaxis()->SetNdivisions(505);
    f->GetZaxis()->SetTitleSize(0.04);
    f->GetZaxis()->SetLabelSize(0.03);
    
    f->SetLineWidth(1);
    f->SetLineColorAlpha(kAzure-2, 0.3);
    
    f->Draw("surf1 fb");
    
    // Lines for 3D annotation
    double x[11] = {-0.500, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.500};
    double y[11] = {-0.985, -0.8, -0.6, -0.4, -0.2,  0.0,  0.2,  0.4,  0.6,  0.8,  0.985};
    double z[11];
    for (int i = 0; i < 11; ++i) z[i] = s*cos(x[i])*cos(y[i]);
    TPolyLine3D *g2 = new TPolyLine3D(11, x, y, z);
    
    double xx[2] = {-0.5, -0.5};
    double yy[2] = {-0.985, -0.985};
    double zz[2] = {0.11, s*cos(-0.5)*cos(-0.985)};
    TPolyLine3D *l2 = new TPolyLine3D(2, xx, yy, zz);
    
    g2->SetLineColor(kRed);
    g2->SetLineWidth(3);
    g2->Draw();
    
    l2->SetLineColor(kRed);
    l2->SetLineStyle(2);
    l2->SetLineWidth(1);
    l2->Draw();
    
    // Draw text Annotations
    TLatex *txt = new TLatex(0.05, 0, "f(y,x_{0})");
    txt->SetTextFont(42);
    txt->SetTextColor(kRed);
    txt->Draw();
    
    TLatex *txt1 = new TLatex(0.12, 0.52, "f(x,y)");
    txt1->SetTextColor(kBlue);
    txt1->SetTextFont(42);
    txt1->Draw();
    
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

//Test
void test_decompositions(void){
    //SVD decomposition
    arma::mat X(5, 5, arma::fill::randu);
    
    arma::mat U;
    arma::vec s;
    arma::mat V;
    
    arma::svd(U,s,V,X);
    
    cout << "X = " << X << endl;
    cout << "U = " << U << endl;
    cout << "s = " << s << endl;
    cout << "V = " << V << endl;
    
    //vecteurs propres valeurs propres

    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    arma::eig_gen(eigval, eigvec, X);
    
    cout << "X = " << X << endl;
    cout << "eigvec = " << eigvec << endl;
    cout << "eigval = " << eigval << endl;
}

void hsumanim(void){
    auto c1 = new TCanvas("c1","The HSUM example",200,10,600,400);
    c1->SetGrid();
    gBenchmark->Start("hsum");
    
    // Create some histograms.
    auto total  = new TH1F("total","This is the total distribution",100,-4,4);
    auto main   = new TH1F("main","Main contributor",100,-4,4);
    auto s1     = new TH1F("s1","This is the first signal",100,-4,4);
    auto s2     = new TH1F("s2","This is the second signal",100,-4,4);
    total->Sumw2();   // this makes sure that the sum of squares of weights will be stored
    total->SetMarkerStyle(21);
    total->SetMarkerSize(0.7);
    main->SetFillColor(16);
    s1->SetFillColor(42);
    s2->SetFillColor(46);
    TSlider *slider = 0;
    gSystem->Unlink("hsumanim.gif"); // delete old file
    
    // Fill histograms randomly
    gRandom->SetSeed();
    const Int_t kUPDATE = 500;
    Float_t xs1, xs2, xmain;
    Int_t gifcnt = 0;
    
    for ( Int_t i=0; i<10000; i++) {
        xmain = gRandom->Gaus(-1,1.5);
        xs1   = gRandom->Gaus(-0.5,0.5);
        xs2   = gRandom->Landau(1,0.15);
        main->Fill(xmain);
        s1->Fill(xs1,0.3);
        s2->Fill(xs2,0.2);
        total->Fill(xmain);
        total->Fill(xs1,0.3);
        total->Fill(xs2,0.2);
        if (i && (i%kUPDATE) == 0) {
            if (i == kUPDATE) {
                total->Draw("e1p");
                main->Draw("same");
                s1->Draw("same");
                s2->Draw("same");
                c1->Update();
                slider = new TSlider("slider","test",4.2,0,4.6,total->GetMaximum(),38);
                slider->SetFillColor(46);
            }
            if (slider) slider->SetRange(0,Float_t(i)/10000.);
            c1->Modified();
            c1->Update();
            if (gROOT->IsBatch()) {
                c1->Print("hsumanim.gif+");
                printf("i = %d\n", i);
                
            }else{
                if (gSystem->ProcessEvents())
                    break;
            }
        }
    }
    
    slider->SetRange(0,1);
    total->Draw("sameaxis"); // to redraw axis hidden by the fill area
    c1->Modified();
    // make infinite animation by adding "++" to the file name
    if (gROOT->IsBatch()) c1->Print("hsumanim.gif++");
    gBenchmark->Show("hsum");
}

//double fct_zoom(double t,double zoom, double t_min , double t_max){
//    return (t_min-t_max)/(M_PI*zoom)*tan(M_PI*(t-0.5))+t_min;
//}

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

