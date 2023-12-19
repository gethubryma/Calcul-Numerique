/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
   int n = *la;

    // Parcourir la matrice en utilisant une seule boucle i (lignes)
    for (int i = 0; i < n; ++i) {

        if (i > 0) {
            AB[*kv + 1 - 1 + i * (*lab) ] = -1.0;
        }

        if (i < n - 1) {
            
            AB[*kv + 1 + 1 + i * (*lab) ] = -1.0;
        }

        AB[*kv + 1 + i * (*lab) ] = 2.0;
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int n = *la;

    for (int i = 0; i < n; ++i) {

        if (i > 0) {
            AB[*kv + 1 - 1 + i * (*lab)] = -1.0;
        }

        if (i < n - 1) {
            AB[*kv + 1 + 1 + i * (*lab)] = -1.0;
        }

        AB[*kv + 1 + i * (*lab)] = 2.0;
    }

    for (int i = 0; i < n; ++i) {
        AB[*kv + 1 + i * (*lab)] += 1.0;
    }

}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
    //resolution d'un systeme lineaire Ax = b cest la fonction qui permet de creer b
    int n=*la;
    RHS[0] = *BC0;
    
    for (int i = 1; i <n- 1; i++) {
        RHS[i] = 0.0; 
    }
    RHS[n-1] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
    // Calcul de la solution analytique T(x) = T0 + x(T1 - T0)
    int n=*la;
     for (int i = 0; i < n; i++) {
        EX_SOL[i] = *BC0 + X[i] * (*BC1 - *BC0);
    }
}  

void set_grid_points_1D(double* x, int* la){
    int n = *la;
    double dx = 1.0 / (n - 1);  // h

    for (int i = 0; i < n; ++i) {
        x[i] = i * dx;
    }
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

int indexABCol(int i, int j, int *lab){
  return 0;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
