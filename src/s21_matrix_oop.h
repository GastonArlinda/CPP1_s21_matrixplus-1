#ifndef _S21_MATRIX_OOP
#define _S21_MATRIX_OOP

#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  // Atributes
  int rows_, cols_;
  double** matrix_;

 public:
  // Methods
  S21Matrix();   // Constructor
  ~S21Matrix();  // Destructor

  S21Matrix(int rows, int cols);      // Constructor with arguments
  S21Matrix(const S21Matrix& other);  // Сopy constructor
  S21Matrix(S21Matrix&& other);       // Transfer constructor

  // Functions
  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  // Accessor & mutator
  int GetRows() const;
  int GetCols() const;
  void SetRows(int rows);
  void SetCols(int cols);

  // Additional Methods
  void AlocMatrix(int rows, int cols, double*** matrix);
  void DelMatrix(double** matrix);
  void create_determinant_matrix(S21Matrix* determ, int minor_row,
                                 int minor_col);
  void algebraic_complements(S21Matrix* result);

  // Operations
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double& num);
  bool operator==(const S21Matrix& other) const;
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(const double& num);
  double& operator()(int i, int j);
};

#endif