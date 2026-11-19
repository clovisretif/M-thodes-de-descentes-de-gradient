#include "ressource.hpp"

double condit_init(double x){ //calcul de la condition initiale pour la première itération
    return std::exp(-4096.0*((x-0.5)*(x-0.5)));
}

double solution_exacte(double t, double x, double mu) {  //calcul solution exacte à partir de celle donnée dasn la partie théorique
    double denom = 1.0 + 16384.0 * mu * t;
    return (1.0 / std::sqrt(denom)) * std::exp(-4096.0 * (x - 0.5) * (x - 0.5) / denom);
}

double erreur_max(const std::vector<double>& U, double tf, double mu, int N) { //calcul de l'erreur max
    double h = 1.0 / (N + 1);
    double erreur = 0.0;

    for(int i = 0; i < N + 2; ++i) {
        double xi = i * h;
        double u_exact = solution_exacte(tf, xi, mu);
        double diff = std::abs(U[i] - u_exact);
        if(diff > erreur){
            erreur = diff;
        }
    }

    return erreur;
}