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
    Matrix* cancleRow(int);
    Matrix* zsf();
    Matrix* elementary_form();              //Code nicht fertig ausgereift, da leider keine Zeit mehr und es hier bei double als Koerper eh keinen Sinn ergibt sie zu bilden
    Matrix* inverse();
    Matrix* scale(double);
    Matrix* transpose();
    Matrix* adjoint();
    Matrix* extractColN(int);

    double* getCol(int);
    double* getEntry(int, int);

    int replaceColumn(int, Matrix);
    void replace(int, Matrix*);

    static Matrix unit_vector(int n) {
        double** etr = new double*[n];
        for (int i = 0; i < n; i++)
        {
            etr[i] = new double[1];
            etr[i][0] = 0;
        }

        etr[0][0] = 1;

        return Matrix(n, 1, etr);
    }
};

