//
//  classes.cpp
//  programme_stage
//
//  Created by Camille Bernard on 22/07/2021.
//

#include "classes.hpp"


ensemble::ensemble(const string input_nom){
    nom = input_nom;
}
ensemble::ensemble(const ensemble & input_ensemble){
    nom = input_ensemble.nom;
}
ensemble & ensemble::operator= (const ensemble & input_ensemble){
    nom = input_ensemble.nom;
    return *this;
}
bool ensemble::indicatrice(void){
    return false;
}

Reel::Reel(const string input_nom): ensemble(input_nom){}

bool Reel::indicatrice(double){
        return true;
}


ouvert::ouvert(const string input_nom){
    nom = input_nom;
}
ouvert::ouvert(const ouvert & input_ouvert){
    nom = input_ouvert.nom;
}
ouvert & ouvert::operator= (const ouvert & input_ouvert){
    nom = input_ouvert.nom;
    return *this;
}


interval_ouvert::interval_ouvert(const string input_nom,pair<double,double> input_bornes): ouvert(input_nom){
    bornes = input_bornes;
}
interval_ouvert::interval_ouvert(const interval_ouvert & input_ouvert): ouvert(input_ouvert){
    bornes = input_ouvert.bornes;
}
interval_ouvert & interval_ouvert::operator= (const interval_ouvert & input_ouvert){
    nom = input_ouvert.nom;
    bornes = input_ouvert.bornes;
    return *this;
}

application::application(const string input_nom, const ouvert input_domaine_def): domaine_def(input_domaine_def){
    nom = input_nom;
}

application::application(const application & input_application): domaine_def(input_application.domaine_def){
    nom = input_application.nom;
}

application & application::operator= (const application & input_ouvert){
    nom = input_ouvert.nom;
    domaine_def = input_ouvert.domaine_def;
    return *this;
}

complex<double>  application::operator() (complex<double> z){
    complex<double> zero(0,1);
    return zero;
}


un::un(const string input_nom, const ouvert input_domaine_def): application(input_nom,input_domaine_def){
    nom = input_nom;
}
complex<double>  un::operator() (complex<double> z){
    complex<double> un(1,0);
    return un;
}

zero::zero(const string input_nom, const ouvert input_domaine_def): application(input_nom,input_domaine_def){
    nom = input_nom;
}
complex<double>  zero::operator() (complex<double> z){
    complex<double> zero(0,0);
    return zero;
}

