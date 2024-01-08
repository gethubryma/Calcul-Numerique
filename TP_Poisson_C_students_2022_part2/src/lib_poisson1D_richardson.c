/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/

#include "lib_poisson1D.h"
#define M_PI 3.14159265358979323846

void eig_poisson1D(double* eigval, int *la){
int n = *la;
for (int k = 0; k < n; ++k) {
  double lambda = 2.0 - (2.0 * cos(M_PI * k * (n + 1)));
  eigval[k] = lambda;
  
}

}

double eigmax_poisson1D(int *la){
  //lambda N
  int n = *la;
  double eigval[n];
  eig_poisson1D(eigval, la);

  double max = eigval[0];
  for (int i = 1; i < n; ++i) {
      if (eigval[i] > max) {
          max = eigval[i];
      }
  }

  return max;
}
  

double eigmin_poisson1D(int *la){
  //lambda 0
  int n = *la;
  double eigval[n];
  eig_poisson1D(eigval, la);

  double min = eigval[0];
  for (int i = 1; i < n; ++i) {
      if (eigval[i] < min) {
          min = eigval[i];
      }
  }

  return min;
}

double richardson_alpha_opt(int *la){
  //alpha_opti=2/(lambda 0 + lambda N)
  double lambda_min = eigmin_poisson1D(la);
  double lambda_max = eigmax_poisson1D(la);

  double alpha_opti = 2.0 / (lambda_min + lambda_max);

  return alpha_opti;
}


//une fonction qui calcule lerreur  par rapport a la solution analytique 

void calculate_error_vector(double* error_vector, double* numerical_solution, double* analytical_solution, int n) {
    for (int i = 0; i < n; i++) {
        error_vector[i] = fabs(numerical_solution[i] - analytical_solution[i]);
    }
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
 //x^k+1=x^k + alpha(b-A*x^k)
    int n = *la;
    int ld = *ku + *kl + 1;  

    //la méthode de Richardson
    for (int k = 0; k < *maxit; k++) {
        //matrice-vecteur A * X
        for (int i = 0; i < n; i++) {
            resvec[i] = 0;
            for (int j = fmax(0, i - 1); j <= fmin(n - 1, i + 1); j++) {
                resvec[i] += AB[1 + i - j + j * ld] * X[j];
            }
        }

        //résidu R = b - A * X
        for (int i = 0; i < n; i++) {
            resvec[i] = RHS[i] - resvec[i];
        }
        for (int i = 0; i < n; i++) {
            X[i] += *alpha_rich * resvec[i];
        }

        //la norme
        double norm = 0;
        for (int i = 0; i < n; i++) {
            norm += resvec[i] * resvec[i];
        }
        norm = sqrt(norm);

        resvec[k] = norm;

        //Sauvegarde de la convergence dans un fichier pour tracer la courbe de convergence 
        FILE *file = fopen("convergence_data.dat", "a"); 
        if (file != NULL) {
        fprintf(file, "%d %lf\n", k + 1, resvec[k]);
        fclose(file);
        }

        if (norm < *tol) {
            *nbite = k + 1;  
            break;
        }

    }
}


void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
      int n = *la;

    // Parcourir la matriceAB
    for (int i = 0; i < n; i++) {
        
    
      for (int j = i - *kl; j <= i + *ku; ++j) {
        if (j >= 0 && j < n) {
          // Extraire 
          if (i == j) {
              MB[i + j * n] = AB[*ku + i - j + i * (*lab)];
          }else if (i - j == *ku) {
              MB[i + j * n] = AB[*ku + i - j + i * (*lab)];
          }else if (j - i == *kl) {
              MB[i + j * n] = AB[*ku + i - j + i * (*lab)];
          } else {
              MB[i + j * n] = 0;
          }
        }
      }
    }
}


void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  int n = *la;

  // Parcourir AB
  for (int i = 0; i < n; ++i) {
      for (int j = i - *kl; j <= i + *ku; ++j) {
          if (j >= 0 && j < n) {
              
              if (i == j || i - j == *ku || j - i == *kl) {
                  MB[i + j * n] = AB[*ku + i - j + i * (*lab)];
              } else {
                  MB[i + j * n] = 0;
              }
          }
      }
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  int n = *la;

  for (int k = 0; k < *maxit; k++) {
    //matrice-vecteur MB * X
    for (int i = 0; i < n; i++) {
        resvec[i] = 0;
        for (int j = fmax(0, i - 1); j <= fmin(n - 1, i + 1); j++) {
            resvec[i] += MB[i + j * n] * X[j];
        }
    }

    // R = b - MB * X
    for (int i = 0; i < n; i++) {
        resvec[i] = RHS[i] - resvec[i];
    }
    for (int i = 0; i < n; i++) {
        X[i] += resvec[i] / MB[i + i * n];
    }

    double norm = 0;
    for (int i = 0; i < n; i++) {
        norm += resvec[i] * resvec[i];
    }
    norm = sqrt(norm);

    resvec[k] = norm;

    //tolerence 
    if (norm < *tol) {
      *nbite = k + 1; 
    }

  }
}




