#include "Matrix.h"

Matrix::Matrix(int r, int c, double** etrs) {               //konstruktor mit anzahl reihen und spalten und den Eintraegen
    rows = r;
    cols = c;
    entries = etrs;
}

Matrix::~Matrix() {                         //Destruktor, der dann im Code verantwortlich ist, die nicht mehr genutzten Matrizen zu loeschen
    cleanUpEntries();
}

void Matrix::cleanUpEntries()
{
    for (int i = 0; i < rows; i++) {
        delete[] entries[i];                //zuerst Inhalte (der Spalten) loeschen
    }
    delete entries;                         //dann auch den Rest (die Zeilen) loeschen
}


Matrix::Matrix(const Matrix& src) {                 //Kopierkonstruktor
    rows = src.rows;                                //gleiche Anzahl an Reihen und Spalten
    cols = src.cols;

    entries = new double*[rows];                    //Matrix der entsprechenden Groesse erstellen
    for (int i = 0; i < rows; i++) {
        entries[i] = new double[cols];
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {            //Matrix dann mit den Source-Eintraegen befuellen
            entries[i][j] = src.entries[i][j];
        }
    }
}


Matrix::Matrix(int size, double diagEntry) {        //Konstruktor fuer eine Diagonalmatrix von bestimmter Groesse und mit einem bestimmten double Eintrag
    rows = cols = size;                             //muss quadratisch sein
    entries = new double*[size];
    for (int i = 0; i < size; i++) {                //Matrix der bestimmten Groesse erstellen
        entries[i] = new double[size];
    }

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (i == j)                             //Matrix auf der Diagonale mit dem bestimmten Diagonaleintrag und
                entries[i][j] = diagEntry;          //auf der Nicht-Hauptdiagonalen der Matrix(den Rest) mit 0 befuellen
            else
                entries[i][j] = 0;
}

string Matrix::print() {
    string res;                                             //Matrix wird als String angesehen

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res += to_string(entries[i][j]) + "\t";     //haengt an den String den i-j-ten Eintrag als String gesehen an und laesst ein bisschen Platz zum naechsten
        }
        res += "\n";                                        //nach jeder fertigen Zeile soll er auch in die naechste wechseln
    }

    return res;
}

int Matrix::getRows() {
    return rows;
}


int Matrix::getCols() {
    return cols;
}
                                                //nun folgt eine Reihe von Operator ueberladen

Matrix & Matrix::operator=(Matrix src) {                //Zuweisungsoperator
    if (this != &src) {                                 //Algorithmus soll nur ausgefuehrt werden, wenn die Matrizen nicht eh schon uebereinnstimmen
        rows = src.rows;
        cols = src.cols;                                //gleiche Anzahl an Reihen und Spalten
        for (int i = 0; i < rows; i++) {
            delete[] entries[i];                        //vorherige Eintraege loeschen
        }
        delete[] entries;

        entries = new double*[rows];                    //Matrix der passenden Groesse wieder zuweisen
        for (int i = 0; i < rows; i++) {
            entries[i] = new double[cols];
        }

        for (int i = 0; i < rows; i++) {                //Matrix mit den Eintraegen aus Source befuellen
            for (int j = 0; j < cols; j++) {
                entries[i][j] = src.entries[i][j];
            }
        }
    }

    return *this;
}


Matrix* Matrix::operator+(Matrix rhs) {
    if(rows != rhs.rows || cols != rhs.cols) {                  //nur Matrizen mit gleicher Zahl an Spalten und gleicher Zahl an Reihen koennen addiert werden
        return nullptr;                                         //falls also dies nicht der fall ist wird Nullpointer zurueckgegeben (keine neue Info)
    }

    double** etrs = new double*[rows];                          //neue Matrix erstellen mit gleicher Anzahl an Reihen und Spalten
    for (int i = 0; i < rows; i++) {
        etrs[i] = new double[cols];
    }

    for (int i = 0; i < rows; i++) {                            //Jeden Eintrag entprechend der math. Formel berechnen
        for (int j = 0; j < cols; j++) {
            etrs[i][j] = entries[i][j] + rhs.entries[i][j];
        }
    }

    return new Matrix(rows, cols, etrs);
}



Matrix* Matrix::scale(double factor) {
    return *this * Matrix(cols, factor);            //Matrix mit einer Diagonalmatrix der entprechenden Laenge (dass Regel der Multiplikation nicht
}                                                                 //verletzt wird) und den Diagonaleintraegen entsprechen dem Skalenfaktor multiplizieren



Matrix* Matrix::operator-(Matrix rhs) {
    return *this + *(rhs.scale(-1));                        //Subtraktion ist Addition mit dem additiven inversen
}



Matrix* Matrix::operator*(Matrix rhs) {
    if (cols != rhs.rows) {
        return nullptr;                                             //math. Bedingung, dass Multiplikation funktionieren kann
    }                                                               //Anzahl Spalten der ersten und Anzahl Zeilen der zweiten muessen gleich sein

    int new_rows = rows;                                            //neue Matrix mit logischerweise Anzahl der Zeilen der ersten und Anzahl der Spalten der zweiten Matrix
    int new_cols = rhs.cols;
    auto** new_entries = new double*[new_rows];

    for (int i = 0; i < new_rows; i++) {                            //Matrix der entsprechenden Groesse wie immer
        new_entries[i] = new double[new_cols];
    }

    for (int i = 0; i < new_rows; i++) {
        for (int j = 0; j < new_cols; j++) {
            double tmp = 0;  // ChangePoint1                        //temporaere double als Zwischenergebnis, dass fuer jeden Eintrag i-j auf 0 zurueckgesetzt wird
            for (int k = 0; k < cols; k++) {
                tmp += entries[i][k] * rhs.entries[k][j];           //Formel fuer den i-j-ten Eintrag fuer Matrizenmultiplikation (vgl. z.B. Wikipedia)
            }

            new_entries[i][j] = tmp;
        }
    }

    return new Matrix(new_rows, new_cols, new_entries);
}



Matrix* Matrix::operator/(Matrix rhs) {
    return *this * *rhs.inverse();                      //Matrizendivision ist Matrizenmultiplikation mit dem Inversen
}

// Matrix* Matrix::cancelRowAndCol(int r, int c) {
//     int new_rows = rows - 1;                                //neue Matrix rows und cols 1 niedriger wegen wegstreichen einer Zeile bzw. Spalte
//     int new_cols = cols - 1;
//
//     double** new_entries = new double*[new_rows];
//     for (int i = 0; i < new_rows; i++) {                    //Matrix der passenden Groesse erstellen
//         new_entries[i] = new double[cols];
//     }
//
//     for (int i = 0; i < rows; i++) {
//         for (int j = 0; j < cols; j++) {                    //alle faelle werden abgedeckt (bei i = r oder j = c wird nichts gemacht also diese eben ausgelassen wie gewuenscht
//             if (i < r && j < c) {
//                 new_entries[i][j] = entries[i][j];          //bis 1 vor gestrichener Zeile und Spalte copy paste alte Matrix
//             }
//
//             if (i < r && j > c) {
//                 new_entries[i][j - 1] = entries[i][j];      //falls aktuelle Spalte von alter Matrix groesser als gestrichene soll das der naechste Eintrag sein in der neuen
//             }
//
//             if (i > r && j < c) {
//                 new_entries[i - 1][j] = entries[i][j];      //gleiches mit Zeilen
//             }
//
//             if (i > r && j > c) {
//                 new_entries[i - 1][j - 1] = entries[i][j];  //gleiches kombiniert Reihen und Zeilen
//             }
//                                                             //insgesamt wird also die i-j-ten Eintraege ausgelassen und damit auch die Matrix verkleinert (Eintraege kopiert)
//         }
//     }
//
//     return new Matrix(new_rows, new_cols, new_entries);
// }

int Matrix::cancelRowAndCol(int r, int c) {
    // wir erstellen temporär die kleinere Matrix und überschreiben uns dann am Ende mit der kleineren

    if (r < 0 || c < 0 || r >= rows || c >= cols)
    {
        return -1;
    }

    int new_rows = rows - 1;
    int new_cols = cols - 1;

    double** new_entries = new double*[new_rows];
    for (int i = 0; i < new_rows; i++) {
        new_entries[i] = new double[cols];
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {                    //alle faelle werden abgedeckt (bei i = r oder j = c wird nichts gemacht also diese eben ausgelassen wie gewuenscht
            if (i < r && j < c) {
                new_entries[i][j] = entries[i][j];          //bis 1 vor gestrichener Zeile und Spalte copy paste alte Matrix
            }

            if (i < r && j > c) {
                new_entries[i][j - 1] = entries[i][j];      //falls aktuelle Spalte von alter Matrix groesser als gestrichene soll das der naechste Eintrag sein in der neuen
            }

            if (i > r && j < c) {
                new_entries[i - 1][j] = entries[i][j];      //gleiches mit Zeilen
            }

            if (i > r && j > c) {
                new_entries[i - 1][j - 1] = entries[i][j];  //gleiches kombiniert Reihen und Zeilen
            }
            //insgesamt wird also die i-j-ten Eintraege ausgelassen und damit auch die Matrix verkleinert (Eintraege kopiert)
        }
    }

    // bevor wir den pointer "entries" überschreiben müssen wir den Speicherplatz dort freigeben:
    cleanUpEntries();

    rows = new_rows;
    cols = new_cols;
    entries = new_entries;

    return 0;
}

int plus_minus(int x) {
    if (x % 2 == 0) {
        return 1;               //wird benoetigt fuer Vorzeichen bei Determinante/Adjungierte etc. vgl Formel (-1)^(i+j)*...
    }
    return -1;
}


double* Matrix::det() {
    if (rows != cols) {
        return nullptr;                 //Matrix muss quadratisch sein
    }

    if (rows == 1) {
        return &entries[0][0];          //1x1 Matrix hat Determinante genau dessen Eintrag
    }

    auto* res = new double(0);  //ChangePoint1

    for (int i = 0; i < rows; i++) {                //Laplacescher Entwicklungssatz
        //ich entwickle hier nach der ersten Spalte (jeweils)
        //entries[i][0] gibt dann jeweiligen Eintrag an (dabei wird bei der Matrix diese Zeile und diese Spalte gestrichen und dessen Determinante berechnet
        // -> rekursiv Bedingung von oben fuer 1x1 greift)
        //insgesamt also Determinante von der gestrichenen Matrix multipliziert mit dem jeweiligen Eintrag und dem jeweiligen Vorzeichen
        // *res += *(cancelRowAndCol(i, 0)->det()) * plus_minus(i) * entries[i][0];
    }

    return res;
}


double* Matrix::trace() {
    if (rows != cols) {
        return nullptr;
    }

    auto* res = new double(0);

    for (int i = 0; i < rows; i++) {        //analog zu det()
        *res += entries[i][i];              //res berechnet sicher aber nur aus Summe der Diagonaleintraege als Unterschied zu det()
    }

    return res;
}


Matrix* Matrix::transpose() {
    int new_rows = cols;
    int new_cols = rows;                //neue Matrix mit Zeilen bzw. Spalten entsprechend dem jeweils anderen der alten Matrix

    double** new_entries =  new double*[new_rows];
    for (int i = 0; i < new_rows; i++) {
        new_entries[i] = new double[new_cols];
    }

    for (int i = 0; i < new_rows; i++) {
        for (int j = 0; j < new_cols; j++) {            //Eintraege einfach ueber Permutation der Indizes befuellen
            new_entries[i][j] = entries[j][i];
        }
    }

    return new Matrix(new_rows, new_cols, new_entries);
}


Matrix* Matrix::adjoint() {     //Adjungierte A^(ad) ist diejenige Matrix zu A, sodass A*A^(ad) = det(A) * E_n (Cramersche Regel)
    if (rows != cols) {         //nur fuer quadratische
        return nullptr;
    }

    double** new_entries = new double*[rows];
    for (int i = 0; i < rows; i++) {                //Matrix erstellen
        new_entries[i] = new double[cols];
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            //jeden Eintrag ueber Regel einfach bestimmen
            //i-te Zeile und j-te Spalte streichen; Determinante davon berechnen; Vorzeichen bestimmen
            // new_entries[i][j] = (*cancelRowAndCol(i, j)->det() * plus_minus(i + j));
        }
    }
    auto tmp = new Matrix(rows, cols, new_entries);
    auto result = tmp->transpose();

    delete tmp;             //schliesslich muss nach der Regel noch transponiert werden und dann wird sie auch gleich so zurueckgegeben

    return result;
}


Matrix* Matrix::inverse() {
    double* d = this->det();

    if (d == nullptr || *d == 0) {          //prueft ob Matrix invertierbar ist
        return nullptr;
    }

    double preFactor = 1/(*d);

    auto adjointMatrix = this->adjoint();                      //bildet Adjungierte
    auto result = adjointMatrix->scale(preFactor);       //skaliert um 1/Determinante nach Cramerscher Regel

    delete adjointMatrix;           //gibt speicher von Adjungierter Matrix wieder frei

    return result;
    /*
    if (rows != cols || this->det() == 0) { // spaeter:  rows != cols || !this->det()->isInvertible() (also wenn man noch die anderen Klassen zugelassen haette)
        return nullptr;
    }

    double d = *(this->det());

    double preFactor = 1/d;                    //spaeter: *(d.inverse()) (also wenn man noch die anderen Klassen zugelassen haette)


    return this->adjoint()->scale(preFactor);
     */
}


void Matrix::swapRows(int i, int j) {
    if (i == j) {
        return;             //falls gleiche Zeile angegeben wurde, soll nichts passieren
    }

    double tmp;             //Zwischenspeicher double

    for (int k = 0; k < cols; k++) {
        tmp = entries[i][k];
        entries[i][k] = entries[j][k];      //klassischer Tausch der Eintraege ueber Zwischenspeicher
        entries[j][k] = tmp;
    }
}


void Matrix::swapCols(int i, int j) {
    if (i == j) {
        return;
    }

    double tmp;                                 //analog

    for (int k = 0; k < cols; k++) {
        tmp = entries[k][i];
        entries[k][i] = entries[k][j];
        entries[k][j] = tmp;
    }
}


void Matrix::subtractRow(int targetRow, int sourceRow, double factor) {
    if (targetRow == sourceRow) {
        return;                                 //falls gleiche Zeile angegeben wurde, soll nichts passieren
    }

    for (int i = 0; i < cols; i++) {
        entries[targetRow][i] -= factor * entries[sourceRow][i];        //von der Ziel-Reihe wird 'factor' mal die Source-Reihe abgezogen
    }
}


void Matrix::subtractCol(int targetCol, int sourceCol, double factor) {
    if (targetCol == sourceCol) {
        return;                                                     //analog
    }

    for (int i = 0; i < rows; i++) {
        entries[i][targetCol] -= factor * entries[i][sourceCol];
    }
}


Matrix* Matrix::attach(int r, int c, double entry) {
    int new_rows = rows + 1;
    int new_cols = cols + 1;

    double** new_entries = new double*[new_rows];
    for (int i = 0; i < new_rows; i++) {
        new_entries[i] = new double[new_cols];
    }

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {               //Matrix copy paste bis zur gewuenschten Zeile/Spalte die ergaenzt werden
            new_entries[i][j] = entries[i][j];
        }
    }

    for (int i = 0; i < new_rows; i++) {            //die ergaenzte Spalte bzw. Zeile wird mit 0 ausgefuellt
        new_entries[i][c] = 0;
    }

    for (int i = 0; i < new_cols; i++) {
        new_entries[r][i] = 0;
    }

    new_entries[r][c] = entry;                      //die Ueberschneidung von Zeile und Spalte wird nun nochmal ueberschrieben mit dem gewuenschten Eintrag

    for (int i = r + 1; i < new_rows; i++) {        //der Rest wird wieder mit copy paste unter Beachtung der nun um im vergleich zur alten Matrix um 1 niedrigeren
        for (int j = c + 1; j < new_cols; j++) {    //Zeile/Spalte ausgefuellt
            new_entries[i][j] = entries[i - 1][j - 1];
        }
    }

    return new Matrix(new_rows, new_cols, new_entries);
}



bool Matrix::isZsf() {
    // gib als bool zurück ob du (also "this") in ZSF ist.
    // in jeder Zeile müssen am Anfang der Zeile mehr Nullen stehen als in allen vorherigen Zeilen.
    int anzahlNull = 0;

    for (int i = 0; i < rows; i++) {
        int aktuellNull = 0;                //wird bei jedem neuem i nochmal neu auf 0 zurueckgesetzt
        for (int j = 0; j < cols; j++) {
            if (entries[i][j] == 0) {
                aktuellNull++;
            }
        }
        if (aktuellNull <= anzahlNull) {
            return false;
        }

        anzahlNull = aktuellNull;

    }

    return true;
}

bool Matrix::isZeroCol(int c, int s) {      //nimmt, um welche Spalte es sich handelt und ab welcher Zeile das ganze gewertet werden soll
    for (int j = s; j < rows; j++) {
        if (entries[j][c] != 0) {           //wenn einer nicht 0 ist dann keine Null-Spalte (Beachte den oben beschriebenen Kontext)
            return false;
        }
    }
    return true;
}

Matrix* Matrix::zsf() {
    auto* result = new Matrix(*this);
    if (rows < 2 || this->isZsf())          //wenn Matrix 1x1 ist, dann fertig; Wenn Matrix schon in ZSF ist dann auch fertig und einfach nur noch ausgeben
        return result;

    int stufe = 0;                      //soll sichern auf welcher Stufe man sich befindet, wenn nun eine Stufe maximal vereinfacht ist und auch wirklich einer
                                        //stufe entspricht soll in die naechste gewechselt werden, um dort auch diese zu vereinfachen und zu einer wirklichen Stufe zu machen

    for (int i = 0; i < cols; i++) {
        if (result->isZeroCol(i, stufe)) {
            continue;                                   //Nullspalten (im Anwendungskontext natuerlich gesprochen) werden ausgelassen
        }
        int rowIndex;
        for (rowIndex = stufe; rowIndex < rows; rowIndex++) {
            if (result->entries[rowIndex][i] != 0) {                //soll bei aktueller Spalte einen Eintrag != 0 auf die Hoehe 'Stufe' tauschen
                break;
            }
        }
        result->swapRows(stufe, rowIndex);

        for (int k = stufe + 1; k < rows; k++) {
            result->subtractRow(k, stufe, result->entries[k][i] / result->entries[stufe][i]);
            //nun sollen die unter Stufe liegenden Zeilen von der mit einem Faktor skalierten Stufenzeile abgezogen werden, sodass in der aktuellen Spalte nur noch
            //Nullen unter Stufe stehen
        }

        stufe++;
        if (stufe == rows - 1) {        //falls nun die vorletzte Stufe erreicht wurde, dann gilt nichts mehr zu tun, da eben in der aktuellen Spalte unterhalb Stufe schon 0
            break;                      //steht und damit die letzte Stufe automatisch schon fertig ist
        }
    }

    return result;
}


int Matrix::rank() {
    auto* myZsf = this->zsf();
    int result = rows - myZsf->countZeroLines();

    delete myZsf;               //gibt speicher von der entstandenen Matrix in ZSF wieder frei
    return result;
}


int Matrix::countZeroLines() {
    int res = 0;
    for (int i = 0; i < rows; i++) {
        double tmp = 0;  // ChangePoint1
        for (int j = 0; j < cols; j++)
            tmp += entries[i][j] * entries[i][j];

        res += tmp == 0;            //wenn tmp == 0 dann ist Ergebnis wahr also 1; wenn nicht dann falsch also 0
                                    //Dementsprechend wird entweder 1 hinzuaddiert, wenn eben eine 0 Zeile vorliegt oder auch nicht
    }

    return res;
}

double* Matrix::getCol(int n){
    double* col = new double[rows];         //muss noch geloescht werden nach Aufruf
    for(int i = 0; i < rows; i++) {
        col[i] = entries[i][n];
    }
    return col;
}

double* Matrix::getEntry(int r, int c) {
    if(r >= rows || c >= cols || r < 0 || c < 0){
        return nullptr;
    }
    return &entries[r][c];
}


int Matrix::replaceColumn(int targetCol, Matrix A) {
    if(targetCol < 0 || targetCol >= cols) {
        return -1;
    }
    for (int i = 0; i < rows; i++) {
        entries[i][targetCol] = *(A.getEntry(i,targetCol));
    }

    return 0;
}

int Matrix::replace(int index, Matrix* A) {
    if(index < 0 || index >= cols) {    // kann nicht vorkommen bei meiner Implementierung (und wir wissen dass rows >= cols!)
        return -1;
    }
    if(A->getRows() != rows - index || A->getCols() != cols - index) {    // kann nicht vorkommen bei meiner Implementierung
        return -1;
    }

    for (int i = index; i < rows; i++) {
        entries[i][index] = *(A->getEntry(i - index,0));
    }

    for (int i = index; i < cols; i++) {
        entries[index][i] = *(A->getEntry(0,i - index));
    }

    return 0;
}

Matrix* Matrix::extractColN(int n)
{
    double** new_entries = new double*[rows];
    for (int i = 0; i < rows; i++)
    {
        new_entries[i] = new double[1];
        new_entries[i][0] = entries[i][n];
    }

    return new Matrix(rows, 1, new_entries);
}

