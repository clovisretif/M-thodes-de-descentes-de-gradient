#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

class Stockage{

    public:
        Stockage(int N, double mu_, double dx_);
        const double& unlock_dt();
        const double& unlock_dx();
        std::vector<double> operator*(const std::vector<double>& x) const;

    private:
        int taille, nb_non_zero;
        double dt, dx;
        std::vector<int> lignes, colonnes;
        std::vector<double> valeurs;
};

//Autre outil

std::vector<std::vector<double>> matrice(int N, double mu_, double dx_); //forme la matrice sans le stockage

double fonctionnelle(const std::vector<double>& x,const std::vector<double>& b,const std::vector<std::vector<double>>& A); //calcul la fonctionnelle

double fonctionnelle(const std::vector<double>& x, const std::vector<double>& b, const Stockage& A); //fonctionnelle pour la classe Stockage

//outil pour calcul:

std::vector<double> multiplication(const std::vector<double>& M, double rho);

std::vector<double> soustraction(const std::vector<double>& M1,const std::vector<double>& M2);

std::vector<double> multiplicationMat(const std::vector<std::vector<double>>& M,const std::vector<double>& x);

std::vector<double> multiplicationMat(const Stockage& S, const std::vector<double>& x); //multiplicationMat pour la classe Stockage

double produitScal(const std::vector<double>& A, const std::vector<double>& B);

double norm(const std::vector<double>& a);

std::vector<double> operator+(const std::vector<double>& u, const std::vector<double>& v);

std::vector<double> operator-(const std::vector<double>& u, const std::vector<double>& v);

//erreur:

double solution_exacte(double t, double x, double mu);

double erreur_max(const std::vector<double>& U, double tf, double mu, int N);

double condit_init(double x);