#pragma once

#include "ressource.hpp"

template <typename T>
void gradientPasFixe(const T& A, const std::vector<double>& b, std::vector<double>& x, double rho, int N){
    
    std::vector<double> grad(x.size());

    for (int k = 0; k < N; k++){

        grad = multiplicationMat(A, x);
        grad = grad-b;

        for (int i=0; i<x.size(); i=i+1){
            x[i] = x[i]-rho*grad[i];
        }
    }
}

template <typename T>
void gradient_optimal(const T& A, const std::vector<double>& b, std::vector<double>& x, int iter, double epsilon){

    std::vector<double>Ax= multiplicationMat(A,x);
    std::vector<double>r = b-Ax;

    for (int k = 0; k < iter; ++k) {
        
        double norm_r = norm(r);
        if (fabs(norm_r-epsilon)<1e-3){
            return;
        }

        std::vector<double> Ar = multiplicationMat(A,r);
        double numerateur = produitScal(r,r);
        double denominateur = produitScal(Ar,r);

        if (denominateur==0){

            break;
        }
        double tau = numerateur / denominateur;
        x = x + multiplication(r,tau);
        r = r - multiplication(Ar,tau);
    }
}

template <typename T>
void gradientConjugue(const T& A, const std::vector<double>& b, std::vector<double>& x, int N) {
    std::vector<double> r=b-multiplicationMat(A,x);
    std::vector<double> p=r;
    double scal=produitScal(r,r);

    for (int k=0; k<N; k++) {
        std::vector<double> stock = multiplicationMat(A, p);
        double alpha = scal / produitScal(p, stock);
        for (int i=0; i<x.size(); i++) {
            x[i] = x[i]+alpha*p[i];
            r[i] = r[i]-alpha*stock[i];
        }
        double scal_2=produitScal(r, r);
        if(sqrt(scal_2)<1e-10){
            break; 
        }
        for(int i=0; i<x.size(); i++){
            p[i] = r[i]+(scal_2/scal)*p[i];
        }
        scal=scal_2;
    }
}