#include <iostream>
#include <string>
#include <cstdlib>
#include <limits>
#include <cmath>

#include "Matrix.cpp"

using namespace std;

/*
 * void headline(string);
void print_matrix(int);
void pause();
bool yes_no_question(string);
 */

Matrix* giveMatrix();
int sgn(double);
double scalarProduct(Matrix*, Matrix*);
double norm(Matrix*);
double calcBeta(Matrix *A);
Matrix* betaVector(int, double);
Matrix* calculate_householder_vector(Matrix*, double);
int doStep(Matrix*, double*);
int householder(double**, double*, int , int);
int rw_subst(double**, double*, int, double*);
int qtb(double**, double*, int, double*)


int main() {


    Matrix* A = giveMatrix();

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
        alpha[i] = 0;
    }

    cout << householder(M, alpha, 3, 5);

    auto foo = new Matrix(5,3,M);
    cout << foo->print();
    return 0;
}

int householder(double** M, double* alpha, int n, int m) {
    /*Strategie:
     * Kopieren EingabeMatrix in A und in res
     * fuehre den ersten Schritt durch
     * Hier wird erste Spalte von A durch Householdervektor ersetzt
     * Dann wird von den restlichen Spalten in A das nach der Formel jeweils richtige Vielfache abgezogen
     * Dann wird res als die veraenderte von A gespeichert
     * nun wird die erste Zeile und erste Spalte von A gestrichen
     * Dann wird analog bei der Veraenderung der neuen Matrix A vorgegangen
     * Nun soll die erste Zeile und erste Spalte von der neuen wiederum veraenderten Matrix in res gespeichert werden
     * Dazu ueberschreibt man genau die Zeile i und Spalte i von res ab dem Diagonalelement (i,i) (so wird dann auch richtige Groesse sichergestellt
     * Der letzte Schritt ist dann genau wie in der VL diesen Prozess einmal noch durchzufuehren bei der Matrix-Groesse 2xk
     * (fuer k der Groesse entsprechend der Reduktion durch Algorithmus)
     * */

    Matrix* A = new Matrix(n, m, M);
    Matrix* res = new Matrix(*A);
    double beta = 0;

    for(int i = 0; i < m - 1; i++) {    // endet bei vorletzter Spalte da in letzter Spalte nichts zu tun
        if(doStep(A, &beta) == -1) {
            return 1;
        }

        alpha[i] = beta;
        res->replace(i, A);     // ersetzt in Matrix "res" ab dem Diagonaleintrag (i,i) die Restzeile und Restspalte durch die erste Zeile/Spalte von Matrix A
        Matrix* tmp = A->cancelRowAndCol(0,0);
        delete A;
        A = tmp;
    }

    // if(doStep(A, &beta) == -1) {   // da wir aus der for-Schleife raus sind: A hat jetzt nur noch 1 Spalte und mit der rufen wir doStep auf
    //     return 1;
    // }
    // alpha[m - 1] = beta;
    // res->replace(m - 1, A);
    // Berechne beta für letzte Spalte
    // Streiche A[0,0] und überschreibe

    delete A;

    // res ist jetzt das Ergebnis. Müssen aber M überschreiben, also mach ich das jetzt:
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            M[i][j] = *(res->getEntry(i,j));
        }
    }
    delete res;

    return 0;
}

int doStep(Matrix* A, double* beta) {
    *beta = calcBeta(A);
    if (*beta == 0) {       // beta == 0 heißt wir haben ne Nullspalte --> Problem
        return -1;
    }

    Matrix* householder_vector = calculate_householder_vector(A, *beta);

    if (A->replaceColumn(0,*householder_vector) == -1) {
        return -1;
    }

    for (int i = 1; i < A->getCols(); i++) {
        Matrix* x = A->extractColN(i);
        double factor = 2 * scalarProduct(householder_vector, x) / scalarProduct(householder_vector,householder_vector);   // v != 0 weil Householder is nie Null
        A->subtractCol(i, 0, factor);
        delete x;
    }

    delete householder_vector;

    return 0;
}

Matrix* calculate_householder_vector(Matrix* A, double beta) {
    Matrix* firstColumnOfA = A->extractColN(0);
    int n = A->getRows();
    Matrix* tmp = Matrix::unit_vector(n).scale(beta);
    Matrix* res = (*firstColumnOfA) - (*tmp);

    delete firstColumnOfA;
    delete tmp;

    return res;
}

double calcBeta(Matrix* A){
    // Ich berechne den Faktor der im Algoritmus eigentlich \alpha heißt (geht hier nicht weil alpha besetzt)
    Matrix* firstColumnOfA = A->extractColN(0);

    double res = -sgn(*(A->getEntry(0,0))) * norm(firstColumnOfA);
    delete firstColumnOfA;

    return res;
}

double norm(Matrix* A){
    // Matrix A ist ein Vektor!
    return sqrt(scalarProduct(A, A));
}

double scalarProduct(Matrix* A, Matrix* B){
    // Matrizen A und B müssen Vektoren sein. Teste das:
    if(A->getCols() > 1 || B->getCols() > 1) {
        cout << "Skalarprodukt nicht berechenbar \n";
        return -1;
    }
    // Außerdem müssen A und B dieselbe Dimension haben:
    if (A->getRows() != B->getRows())
    {
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
    }

    return -1;
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


int rw_subst(double** A, double* alpha, int n, double* b)
{
    for(int i = n - 1 ; i >= 0; i--)
    {
        double sum = 0;
        if(alpha[i] == 0)   // Wenn Diagonale == 0 --> nicht lösbar
        {
            return 1;
        }
        for(int j = n - 1; j > i; j--)
        {
            sum += A[i][j] * b[j];
        }
        b[i] = (b[i] - sum) / alpha[i];
    }

    return 0;
}

int qtb(double** M, double* alpha, int n, double* b) {
    double** b_etrs = new double*[n];
    for (int i = 0; i < n; i++) {
        b_etrs[i] = new double[1];
        b_etrs[i][0] = b[i];
    }

    Matrix* A = new Matrix(n, n, M);
    Matrix* b_tmp = new Matrix(n, 1, b_etrs);
    Matrix* res = new Matrix(*b_tmp);
    for(int i = 0; i < n - 1; i++) {
        Matrix* v_house = A->extractColN(0);
        double factor = 2 * scalarProduct(v_house, b_tmp) / norm(v_house);
        b_tmp = b_tmp - v_house->scale(factor);

        res->replace(i, b_tmp);     //ersetzt alles ab dem Eintrag i nach unten
        Matrix* tmp = b_tmp->cancelRow(0);          //streicht die erste Zeile also den ersten Eintrag, damit wir den kleineren Vektor b_tmp betrachten koennen zum iterieren
        delete b_tmp;
        b_tmp = tmp;

        Matrix* A_smaller = A->cancelRowAndCol(0,0);        //Matrix wird kleiner, damit wieder die Householdervektoren der richtigen Groesse extracted werden koennen
        delete A;
        A = A_smaller;
    }
    delete b_tmp;
    delete A;

    //res ist jetzt daqs Ergebnis. Muessen aber b ueberschreiben, als mahc ihc das jetzt:
    for(int i = 0; i < n; i++) {
        b[i] = *(res->getEntry(i,0));
    }

    delete res;
    delete A;

    //Habe definitiv noch memory leaks mit b_etrs

    return 0;
}






/*
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
 */