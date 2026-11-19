#include "ressource.hpp"
#include "gradient.hpp"

int main() {

    //Variables à changer pour les tests: 
    int N = 32*2*2*2*2;
    double alpha = 0.5; // alpha = 1 le plus optimal
    double Tf = 0.02;
    int iter = 1000;
    
    //Autres initialisation:
    int dim = N;
    double mu = 0.3;
    double dx = 1.0/(N+1);
    double dt = alpha *(dx*dx)/(2.0*mu);
    double t =0.0;
    int choix_mat;
    int choix_grad;
    std::vector<std::vector<double>> A_classique;
    Stockage* A_stock = nullptr;

    std::cout << "Quelle méthode choisir ? (0)stockage creux, (1)stockage classique";
    std::cin >> choix_mat;

    std::cout << std::endl;

    switch (choix_mat){
        case 0:{

            //génération avec le stockage creux

            A_stock = new Stockage(dim, dt, dx);
            break;
        }

        case 1:{

            //génération avec stockage classique

            A_classique = matrice(dim, dt, dx);
            break;
        }
        
        default:{
            std::cout << "problème de choix stockage" << std::endl;
            abort();
        }
    }

    // Initialisation de la solution

    std::vector<double> x(dim + 2);//vecteur de la température 

    for(int i=0; i<dim+2; ++i) {
        double pos = i*dx;
        x[i] = condit_init(pos);
    }
    x[0] = 0.;
    x[dim+1] = 0.;

    std::cout << "Quel gradient voulez-vous choisir ? (0) Pas fixe, (1) Optimal, (2) Conjugue";
    std::cin >> choix_grad;
    while(t<Tf){
        std::vector<double>b=x;

        switch (choix_grad){
            case 0:{
                //test grad pas fixe
                double rho = -0.0000001;

                if(choix_mat == 0){
                    gradientPasFixe(*A_stock, b, x, rho, iter);
                    break;
                }

                gradientPasFixe(A_classique, b, x, rho, iter);
                break;
                
            }
            case 1:{
                //test grad optimal
                double epsilon = 0.01;

                if(choix_mat==0){
                    gradient_optimal(*A_stock, b, x, iter, epsilon);
                    break;
                }

                gradient_optimal(A_classique, b, x, iter, epsilon);
                break;
                
            }
            case 2:{
                //test grad conj

                if(choix_mat==0){
                    gradientConjugue(*A_stock, b, x, iter);
                    break;
                }

                gradientConjugue(A_classique, b, x, iter);
                break;
                
            }
            default:{
                std::cout << "problème de choix gradient" << std::endl;
                abort();
            }
        }
        t+=dt;
    }

    std::ofstream fichier("resultats.txt");
    if(fichier){
        for (int i=0; i<dim+2; i++){
            fichier << (i*dx) << " " << x[i] << "\n";
        }
        fichier.close();
    }else{
        std::cout << "Erreur" << std::endl;
    }

    double erreur_finale = erreur_max(x, t, mu, dim);

    std::cout << std::endl;

    std::cout << "terminée au temps t = " << t << std::endl;
    std::cout << "Erreur (Max): " << erreur_finale << std::endl;

    delete A_stock;

    return 0;
}