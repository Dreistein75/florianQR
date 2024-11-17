#include <iostream>
#include <string>
#include <cstdlib>
#include <limits>
#include <cmath>

#include "Matrix.cpp"

void headline(string);
void print_matrix(int);
void pause();
bool yes_no_question(string);
Matrix* giveMatrix();
int sgn(double);
double scalarProduct(Matrix*, Matrix*);
double norm(Matrix*);
//Berechne Householder-Vektoren  v = a - beta * e_1;    calcBeta = -sgn(a_1) norm(a)
double calcBeta(Matrix *A);
Matrix* betaVector(int, double);
Matrix* calculate_householder_vector(Matrix*, double);
int doStep(Matrix*, double*);
int householder(double**, double*, int , int);
int householder(Matrix*, double*);
int rw_subst(double**, double*, int, double*);