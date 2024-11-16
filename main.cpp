#include <iostream>
#include <string>
#include <cstdlib>
#include <limits>
#include <cmath>

//#include "Matrix.h"
#include "Matrix.cpp"

using namespace std;

void headline(string);
void print_matrix(int);
void pause();
bool yes_no_question(string);
Matrix* giveMatrix();
Matrix* buildVector(int, double*);
Matrix* extractColN(Matrix*, int);
Matrix* extractFirstCol(Matrix*);
int sgn(double);
double scalarProduct(Matrix*, Matrix*);
double norm(Matrix*);
//Berechne Householder-Vektoren  v = a - beta * e_1;    calcBeta = -sgn(a_1) norm(a)
double calcBeta(Matrix *A);
Matrix* betaVector(int, double);
Matrix* v_h(Matrix*);
int doStep(Matrix**, double*);
int householder(double**, double*, int , int);

int main() {
    Matrix* A = giveMatrix();
    Matrix* a = extractFirstCol(A);

    cout << A->print();
    cout << a->print();
    cout << norm(a) << endl;
    cout << calcBeta(a) << endl << endl;

    double** M = new double*[5];
    for (int i = 0; i < 5; i++) {
        M[i] = new double[3];
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 3; j++) {
            M[i][j] = *(A->getEntry(i,j));
        }
    }

    double* alpha = new double[5];
    for (int i = 0; i < 5; i++) {
        alpha[i] = 1;
    }


    cout << householder(M, alpha, 5, 3);

    delete alpha;
    delete *M;
    delete M;
    delete A;
    delete a;
    return 0;
}

int householder(double** M, double* alpha, int n, int m){
    Matrix* A = new Matrix(n, m, M);
    Matrix* res = new Matrix(*A);
    double beta = 0;

    for(int i = 0; i < m; i++){
        if(doStep(&A, &beta) == -1){
            return 1;
        }

        alpha[i] = beta;
        res->replace(i, A);
        Matrix* tmp = A->cancelRowAndCol(0,0);
        delete A;
        A = tmp;
    }
    delete A;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            M[i][j] = *(res->getEntry(i,j));
        }
    }
    delete res;
    return 0;

}

int doStep(Matrix** A, double* beta) {
    *beta = calcBeta(*A);
    if (*beta == 0) {
        return -1;
    }

    Matrix* v = v_h(*A);

    if ((*A)->replaceColumn(0,*v) == -1) {
        return -1;
    }

    for (int i = 1; i < (*A)->getCols(); i++) {
        Matrix* x = extractColN(*A, i);
        double factor = 2 * scalarProduct(v, x) / scalarProduct(v,v);
        (*A)->subtractCol(i, 0, factor);
        delete x;
    }

    return 0;
}

Matrix* v_h(Matrix* A) {
    Matrix* columnOfA = extractFirstCol(A);
    Matrix* tmp = betaVector(A->getRows(), calcBeta(A));
    Matrix* res = (*columnOfA) - (*tmp);

    delete columnOfA;
    delete tmp;
    return res;
}

Matrix* betaVector(int n, double value){
    double tmp[n];
    tmp[0] = value;
    for (int i = 1; i < n; i++) {
        tmp[i] = 0;
    }
    return buildVector(n, tmp);
}

double calcBeta(Matrix* A){
    return -sgn(*(A->getEntry(0,0))) * norm(A);;
}

double norm(Matrix* A){
    return sqrt(scalarProduct(A, A));
}

double scalarProduct(Matrix* A, Matrix* B){
    if(A->getCols() > 1 || A->getCols() != B->getCols()) {
        cout << "Skalarprodukt nicht berechenbar \n";
        return -1;
    }
    double sum = 0;
    for (int i = 0; i < A->getRows(); i++) {
        double* a = A->getEntry(i, 0);
        double* b = B->getEntry(i, 0);
        if (a == nullptr || b == nullptr) {
            cout << "Skalarprodukt nicht berechenbar (kann hier eh nicht passieren) \n";
            return -1;
        }
        sum += (*a) * (*b);
    }

    return sum;
}

int sgn(double x) {
    if(x >= 0) {
        return 1;
    } else {
        return -1;
    }
}

/*Matrix* extractFirstCol(Matrix* A){
    double* entries = A->getCol(0);
    Matrix* res = buildVector(A->getRows(), entries);
    delete[] entries;
    return res;
}*/
Matrix* extractFirstCol(Matrix* A) {
    return extractColN(A, 0);
}

Matrix* extractColN(Matrix* A, int N){
    double* entries = A->getCol(N);
    Matrix* res = buildVector(A->getRows(), entries);
    delete[] entries;
    return res;
}

Matrix* buildVector(int n, double* entries) {
    double** a_i = new double*[n];
    for(int i = 0; i < n; i++) {
        a_i[i] = new double[1];
        a_i[i][0] = entries[i];
    }

    return new Matrix(n, 1, a_i);
}

Matrix* giveMatrix() {
    double** coeff = new double*[5];
    for(int i = 0; i < 5; i++) {
        coeff[i] = new double[3];
    }
    coeff[0][0] = 2;
    coeff[0][1] = sqrt(5);
    coeff[0][2] = 0;
    coeff[1][0] = 2;
    coeff[1][1] = 3;
    coeff[1][2] = -1;
    coeff[2][0] = 2;
    coeff[2][1] = 3;
    coeff[2][2] = -1;
    coeff[3][0] = 2;
    coeff[3][1] = 0;
    coeff[3][2] = -1;
    coeff[4][0] = 0;
    coeff[4][1] = sqrt(7);
    coeff[4][2] = 15/sqrt(7);

    return new Matrix(5, 3, coeff);
}

void headline(string text) {
    system("clear");                    //funktioniert nur auf linux (Windows cls) loescht screen einmal
    cout << "\n\t" << text << "\n\t";

    for (int i = 0; i<text.length(); i++) {     //unterstrichen
        cout << "=";
    }

    cout << "\n\n";
}

void print_matrix(Matrix* A) {
    headline("Ergebnis");

    cout << A->print() << "\n\n\n";          //Matrix ausgeben

    pause();
}

void pause() {
    cout << endl << "Druecke <Enter>, um weiterzumachen...";            //Ziel ist, dass einfach in Ruhe der Output gelesen werden kann.
    cin.ignore(numeric_limits<streamsize>::max(), '\n');        //Man muss also zuerst Enter druecken, um weiterzumachen
    cin.get();
}

bool yes_no_question(string question) {
    char answer;
    do {
        cout << "\n\n\n" << question << " (y/n)?\t";        //man soll einfach nur mit y/n ein ja oder ein nein an das system weitergeben
        cin >> answer;
    } while(answer != 'y' && answer != 'n');

    return answer == 'y';
}

int rw_subst(double** A, double* alpha, int n, double* b)
{
    // bau dir hier eine Matrix R   <--- wirst du wahrscheinlich nicht brauchen
    // Löse mit Rückwärtssubstitution  Rx = b
    // Überschreibe b mit der Lösung


}