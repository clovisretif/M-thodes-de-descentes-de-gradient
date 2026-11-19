#include "ressource.hpp"

//constructeur de la matrice avec le stockage:

Stockage::Stockage(int N, double dt_, double dx_){
    //Initialise les composantes connues du stockage à l'aide de fonction source pour dt et dx
    taille = N+2;
    double mu = 0.3;
    dx = dx_;
    dt = dt_;

    double alpha = 2*mu*dt / (dx*dx);

    //initialisation à 0 de chaque composante inconnues du stockage pour ensuite les remplirs
    nb_non_zero = 0;
    valeurs = {};
    lignes = {};
    colonnes = {};

    double diag = 1+alpha;
    double non_diag = -0.5*alpha;

    //gestion des indices 0,0 et N+1,N+1

    valeurs.push_back(1.0);
    lignes.push_back(0);
    colonnes.push_back(0);
    nb_non_zero++;

    valeurs.push_back(1.0);
    lignes.push_back(taille-1);
    colonnes.push_back(taille-1);
    nb_non_zero++;

    for(int i=1; i<taille-1; ++i){
        nb_non_zero += 1;
        valeurs.push_back(diag);
        lignes.push_back(i);
        colonnes.push_back(i);
    }
    for(int i = 2; i < taille - 1; ++i){
        //sauvegarde des données au dessus de la diagonale
        nb_non_zero +=1;
        valeurs.push_back(non_diag);
        lignes.push_back(i-1);
        colonnes.push_back(i);

        //sauvegarde des données en dessous de la diagonale
        nb_non_zero += 1;
        valeurs.push_back(non_diag);
        lignes.push_back(i);
        colonnes.push_back(i-1);
    }


}

//On récupère dt et dx car ces deux variables privées sont comme "bloquées":

const double& Stockage::unlock_dt(){
    return dt;
}

const double& Stockage::unlock_dx(){
    return dx;
}

//Fonction qui surcharge l'opérateur * pour la classe stockage 

std::vector<double> Stockage::operator*(const std::vector<double>& x) const{
    if(x.size()!=taille){
        std::cout<<"x admet une taille différente que celle de la matrice: erreur"<<std::endl;
        exit(0);
    }

    std::vector<double> y(taille,0.);
    for (int i=0; i<valeurs.size(); ++i) {
        y[lignes[i]]+=valeurs[i]*x[colonnes[i]];
    }

    return y;
}

//construit la matrice avec le stockage "classique":

std::vector<std::vector<double>> matrice(int N, double dt_, double dx_){
    int taille = N+2;
    double mu = 0.3;
    double alpha = 2*mu*dt_/(dx_*dx_);
    std::vector<std::vector<double>> A(taille, std::vector<double>(taille, 0.0));
    A[0][0] = 1;
    A[taille-1][taille-1] = 1;
    for(int i=1; i<taille-1; ++i){
        A[i][i]=1+alpha;
    }
    for(int i=2; i<taille-1; ++i){
        A[i-1][i]= -0.5*alpha;
        A[i][i-1]= -0.5*alpha;
    }
    return A;
}

//surcharge des opérateur: 

std::vector<double> operator-(const std::vector<double>& u, const std::vector<double>& v) {
    std::vector<double> resultat(u.size());
    for(size_t i=0; i<u.size(); ++i){
        resultat[i]=u[i]-v[i];
    } 
    return resultat;
}

std::vector<double> operator+(const std::vector<double>& u, const std::vector<double>& v) {
    std::vector<double> resultat(u.size());
    for(size_t i=0; i<u.size(); ++i){
        resultat[i]=u[i]+v[i];
    } 
    return resultat;
}

//Autre fonctions:

std::vector<double> multiplication(const std::vector<double>& M, double rho){
    std::vector<double> R = M;
    for (int i = 0; i < R.size(); i = i + 1){
        R[i] = R[i] * rho;
    }
    return R;
}

std::vector<double> multiplicationMat(const std::vector<std::vector<double>>& M,const std::vector<double>& x){
    std::vector<double> y(x.size(), 0.0);
    for (int i = 0; i < x.size(); i = i + 1){
        for (int j = 0; j < x.size(); j = j + 1){
            y[i] = y[i] + M[i][j] * x[j];
        }
    }
    return y;
}

std::vector<double> multiplicationMat(const Stockage& M, const std::vector<double>& x){ //fonction de pour la classe stockage
    return M*x;
}

double produitScal(const std::vector<double>& A, const std::vector<double>& B){
    double prod = 0.0;
    for (int i = 0; i < A.size(); i = i + 1){
        prod = prod + A[i] * B[i];
    }
    return prod;
}

double fonctionnelle(const std::vector<double>& x,const std::vector<double>& b,const std::vector<std::vector<double>>& A){
    std::vector<double> M = multiplicationMat(A, x);
    double scal1 = produitScal(M, x);
    double scal2 = produitScal(b, x);
    double J = 0.5 * scal1 - scal2;
    return J;
}

double fonctionnelle(const std::vector<double>& x, const std::vector<double>& b, const Stockage& A){ //fonctionnelle de la classe Stockage
    std::vector<double> Ax = A*x;
    double scal1 = produitScal(Ax,x);
    double scal2 = produitScal(b,x);
    return 0.5*scal1-scal2;
}

double norm(const std::vector<double>& a){
    return std::sqrt(produitScal(a,a));
}