#include <string>

using namespace std;

class Matrix {
private:                //alles was nur innerhalb der Klasse Matrix benutzt wird (also eben nicht in main)
    int rows;
    int cols;
    double** entries;

    bool isZsf();
    bool isZeroCol(int, int);

    int countZeroLines();

public:
    Matrix(int, int, double**);
    Matrix(const Matrix&);
    Matrix(int, double);

    ~Matrix();
    string print();

    Matrix* operator+(Matrix);
    Matrix* operator-(Matrix);
    Matrix* operator*(Matrix);
    Matrix* operator/(Matrix);
    Matrix & operator=(Matrix);

    double* det();
    double* trace();
    int rank();
    int getRows();
    int getCols();

    void swapRows(int, int );
    void swapCols(int, int );               //wuerde fuer elementarteilerform benutzt werden
    void subtractRow(int, int, double);
    void subtractCol(int, int, double);

    Matrix* attach(int, int, double);       //wuerde fuer elementarteilerform benutzt werden
    Matrix* cancelRowAndCol(int, int);
    Matrix* zsf();
    Matrix* elementary_form();              //Code nicht fertig ausgereift, da leider keine Zeit mehr und es hier bei double als Koerper eh keinen Sinn ergibt sie zu bilden
    Matrix* inverse();
    Matrix* scale(double);
    Matrix* transpose();
    Matrix* adjoint();

    double* getCol(int);
    double* getEntry(int, int);

    int replaceColumn(int, Matrix);
    void replace(int, Matrix*);

};

