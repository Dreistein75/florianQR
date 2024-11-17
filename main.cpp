#include "main.h"

using namespace std;

int aufgabe1()
{
    double** M = new double*[5];
    for (int i = 0; i < 5; i++) {
        M[i] = new double[3];
    }

    M[0][0] = 2;
    M[0][1] = sqrt(5);
    M[0][2] = 0;
    M[1][0] = 2;
    M[1][1] = 3;
    M[1][2] = -1;
    M[2][0] = 2;
    M[2][1] = 3;
    M[2][2] = -1;
    M[3][0] = 2;
    M[3][1] = 0;
    M[3][2] = -1;
    M[4][0] = 0;
    M[4][1] = sqrt(7);
    M[4][2] = 15/sqrt(7);

    double* alpha = new double[5];
    for (int i = 0; i < 5; i++) {
        alpha[i] = 0;
    }

    if (householder(M, alpha, 3, 5) == 1)
    {
        return 1;
    }

    headline("Aufgabe 1 * Matrix R");
    cout << fixed << setprecision(3);
    for (int i = 0; i < 5; i++, cout << endl)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i > j) { cout  << 0. << "\t\t"; }
            if (i == j) { cout << alpha[i] << " \t\t"; }
            if (i < j) { cout << M[i][j] << " \t\t"; }
        }
    }

    pause("Fuer die H-Vektoren zweimal Enter bitte");
    headline("Aufgabe 1 * Householders");

    for (int i = 0; i < 5; i++, cout << endl)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i >= j) { cout  << M[i][j] << "\t\t"; }
            if (i < j) { cout << "  \t\t"; }
        }
    }

    for (int i = 0; i < 5; i++)
    {
        delete[] M[i];
    }
    delete[] M;
    delete[] alpha;

    return 0;
}

int aufgabe2()
{
    double** M = new double*[4];
    for (int i = 0; i < 4; i++) {
        M[i] = new double[4];
    }

    M[0][0] = 2;
    M[0][1] = -1;
    M[0][2] = 3;
    M[0][3] = 0;
    M[1][0] = 7;
    M[1][1] = 0;
    M[1][2] = sqrt(2);
    M[1][3] = 1;
    M[2][0] = 1;
    M[2][1] = 1;
    M[2][2] = 1;
    M[2][3] = 1;
    M[3][0] = 0;
    M[3][1] = 10;
    M[3][2] = 3;
    M[3][3] = 2;

    double* alpha = new double[4];
    for (int i = 0; i < 4; i++) {
        alpha[i] = 0;
    }

    double* b = new double[4];
    b[0] = 1;
    for (int i = 1; i < 4; i++) {
        b[i] = 0;
    }

    if (solve_QR(M, alpha, 4, b) == 1)
    {
        return 1;
    }

    headline("Aufgabe 2 * Matrix R");
    cout << fixed << setprecision(3);
    for (int i = 0; i < 4; i++, cout << endl)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i > j) { cout  << 0. << "\t\t"; }
            if (i == j) { cout << alpha[i] << " \t\t"; }
            if (i < j) { cout << M[i][j] << " \t\t"; }
        }
    }

    pause("Fuer die H-Vektoren zweimal Enter bitte");
    headline("Aufgabe 2 * Householders");

    for (int i = 0; i < 4; i++, cout << endl)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i >= j) { cout  << M[i][j] << "\t\t"; }
            if (i < j) { cout << "  \t\t"; }
        }
    }

    pause("Fuer den Loesungsvektor x zweimal Enter bitte");
    headline("Aufgabe 2 * Loesungsvektor x");

    for (int i = 0; i < 4; i++, cout << endl)
    {
        cout  << b[i] << "\t\t";

    }

    for (int i = 0; i < 4; i++)
    {
        delete[] M[i];
    }
    delete[] M;
    delete[] alpha;
    delete[] b;

    return 0;

}

int main() {
    return aufgabe1();
}

int householder(Matrix* M, double* alpha)
{
    /*
     * Strategie: Ein Arbeitsschritt arbeitet mit Zeile 1 und Spalte 1. Wenn das fertig ist, schneiden wir
     * die erste Zeile und Spalte weg, erhalten eine kleinere Matrix und fangen wieder von vorne an.
     * Damit wir durch das wegschneiden nix verlieren: Ich erstelle eine Matrix "working" die eine Kopie von M ist mithilfe
     * des Kopierkonstruktors. Ich führe dann alle Arbeitsschritte mit "work" durch anstatt mit M und direkt
     * vor dem Wegschneiden von Spalte 1 udn Zeile 1 von work überschreibe ich meine Hauptmatrix (ist auch die Ergebnismatrix) M
     * indem ich die entsprechende Spalte/Zeile von "work" an die richtige Stelle von M schreibe.
     * Genauerer Kommentar zum überschreiben kommt an der Stelle wo ich überschreibe.
     * */

    Matrix* work = new Matrix(*M);
    const int numberOfCols = M->getCols();
    double beta;    // Variable für den Diagonaleintrag im i-ten Schritt (s.u.); wird nach Berechnung in Vektor alpha gespeichert

    for (int i = 0; i < numberOfCols; i++)
    {
        // Folgende Zeile macht den Arbeitsschritt des Algorithmus angewandt auf "work" und überschreibt dabei "work" und überschreibt beta
        // deswegen wird auch die Adresse von beta übergeben, nicht beta als Wert. Das "stepResult" ist -1 bei Fehlern
        // weil ich das so gewohnt bin MINUS 1 zu benutzen und sonst 0

        int stepResult = doStep(work, &beta);
        if (stepResult == -1)
        {
            return 1;   // Laut Aufgabenblatt müssen wir bei Fehlern eine Plus 1 returnen
        }

        alpha[i] = beta;

        // work wurde in doStep() überschrieben. Speichern wollen wir die erste Spalte/Zeile von work bevor wir sie dann
        // wegschneiden und die for-Schleife von vorne beginnen mit einem work das kleiner ist weil wir weggeschnitten
        // haben.

        if (M->replace(i, work) == -1)
        {
            // die replace()-Methode überschreibt in Matrix "M" folgende Einträge:
            // Vom Diagonaleintrag (i,i) senkrecht runter wird ersetzt durch die erste Spalte von work
            // und vom Diagonaleintrag (i,i) wird waagrecht nach rechts ersetzt durch erste Zeile von work
            // Der Rückgabewert ist -1 wenn was schief läuft (was es bei sinnvollen Eingaben nicht tut!)

            return 1;
        }

        if (i < numberOfCols - 1)   // Streichen der ersten Zeile/Spalte nur Sinnvoll wenn noch mindestens 2 Spalten da sind.
        {
            if (work->cancelRowAndCol(0,0) == -1)  // Streiche in work Spalte 1 und Zeile 1 (haben Index 0); Rückgabe = -1 nur wenn Fehler (hier: Nie)
            {
                return 1;
            }
        }
    }

    delete work;

    // M und alpha sind jetzt in der gewünschten Form

    return 0;
}

int householder(double** M, double* alpha, int n, int m) {  // suboptimale Funktionsparameter Reihenfolge da jetzt erst Spalten-, dann Zeilen-Zahl übergeben wird
    // Da ich mit Matrizen arbeiten möchte (war meine Vertiefung im C-Kurs), wandle ich das double** array
    // in eine Matrix um, mache dann die ganze Arbeit und wandle dann zurück um die Signatur dieser
    // Funktion hier so zu lassen wie sie laut Angabe sein soll

    Matrix* A = new Matrix(m, n, M);    // erstelle (mxn)-Matrix mit Einträgen gespeichert aus M

    // Rufe nun die householder(Matrix*, double*) Funktion auf. Da alles mit pointern gelöst wird überschreibt sie bereits
    // das double**-Array M wie gewünscht, also ich danach nichts mehr zu tun und wir leiten den Rückgabewert einfach weiter

    return householder(A, alpha);
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

void headline(string text) {
    system("clear");                    //funktioniert nur auf linux (Windows cls) loescht screen einmal
    cout << "\n\t" << text << "\n\t";

    for (int i = 0; i<text.length(); i++) {     //unterstrichen
        cout << "=";
    }

    cout << "\n\n";
}

void pause(string msg) {
    cout << endl << msg;            //Ziel ist, dass einfach in Ruhe der Output gelesen werden kann.
    cin.ignore(numeric_limits<streamsize>::max(), '\n');        //Man muss also zuerst Enter druecken, um weiterzumachen
    cin.get();
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
        b[i] = (b[i] - sum) / alpha[i];     // Nenner ungleich Null
    }

    return 0;
}

int qtb(double** M, double* alpha, int n, double* b) {  // wir brauchen kein alpha aber laut Blatt muss ich es in der Signatur haben: okay.
    // mache aus b einen Vektor (aka eine nx1-Matrix):
    double** b_etrs = new double*[n];
    for (int i = 0; i < n; i++) {
        b_etrs[i] = new double[1];
        b_etrs[i][0] = b[i];
    }
    Matrix* b_vect = new Matrix(n, 1, b_etrs);

    Matrix* A = new Matrix(n, n, M);
    Matrix* B = new Matrix(*A);     // Eine Kopie von A damit wir beim Streichen von Zeilen/Spalten unser M nicht kaputt machen --> wir arbeiten nur mit B

    // Matrix* res = new Matrix(*b_vect);
    for(int i = 0; i < n; i++) {
        Matrix* v_house = B->extractColN(0);
        double factor = 2 * scalarProduct(v_house, b_vect) / scalarProduct(v_house, v_house);       // Nenner ungleich Null

        //mach das schöner:
        Matrix* scaled_householder = v_house->scale(factor);
        Matrix* tmp_res = (*b_vect) - (*scaled_householder);

        delete v_house;
        delete scaled_householder;
        delete b_vect;
        b_vect = tmp_res;

        b[i] = *(b_vect->getEntry(0,0));

        b_vect->cancelFirstRow();          //streicht die erste Zeile also den ersten Eintrag, damit wir den kleineren Vektor b_tmp betrachten koennen zum iterieren

        B->cancelRowAndCol(0,0);        //Matrix wird kleiner, damit wieder die Householdervektoren der richtigen Groesse extracted werden können
    }
    delete b_vect;
    delete B;

    return 0;
}

int solve_QR(double** M, double* alpha, int n, double* b)
{
    if (householder(M, alpha, n, n) == 1)
    {
        return 1;
    }

    if (qtb(M, alpha, n, b) == 1)
    {
        return 1;
    }

    if (rw_subst(M, alpha, n, b) == 1)
    {
        return 1;
    }

    return 0;


    // Fun fact: Diese gesamte Funktion is äquivalent zu:
    return (householder(M, alpha, n, n) == 1) || (qtb(M, alpha, n, b) == 1) || (rw_subst(M, alpha, n, b) == 1);
}