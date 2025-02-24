#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(1), cols_(1) {
  AlocMatrix(rows_, cols_, &matrix_);
}

S21Matrix::~S21Matrix() {
  if (matrix_ != nullptr) {
    DelMatrix(matrix_);
  }
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  AlocMatrix(rows_, cols_, &matrix_);
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(std::move(other.rows_)),
      cols_(std::move(other.cols_)),
      matrix_(std::move(other.matrix_)) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

void S21Matrix::AlocMatrix(int rows, int cols, double*** matrix) {
  if (rows < 1) {
    throw std::invalid_argument("Invalid number of rows " +
                                std::to_string(rows));
  }
  if (cols < 1) {
    throw std::invalid_argument("Invalid number of cols " +
                                std::to_string(cols));
  }

  *matrix = new double*[rows]();
  for (int i = 0; i < rows; i++) {
    (*matrix)[i] = new double[cols]();
  }
}

void S21Matrix::DelMatrix(double** matrix) {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < rows_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) {
        return false;
      }
    }
  }

  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument(
        "Invalid argument, different matrix dimensions");
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < rows_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument(
        "Invalid argument, different matrix dimensions");
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < rows_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < rows_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument(
        "Invalid argument, the number of columns of the first matrix is not "
        "equal to the number of rows of the second matrix");
  }

  S21Matrix res(rows_, other.cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        res.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = std::move(res);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix transp_(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      transp_.matrix_[j][i] = matrix_[i][j];
    }
  }
  return transp_;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::invalid_argument("Invalid argument, matrix is not square");
  }

  S21Matrix res(rows_, cols_);

  if (rows_ < 2) {
    res.matrix_[0][0] = matrix_[0][0];
  } else {
    S21Matrix determ(rows_ - 1, cols_ - 1);
    for (int minor_row = 0; minor_row < rows_ && rows_ > 1; minor_row++) {
      for (int minor_col = 0; minor_col < cols_; minor_col++) {
        create_determinant_matrix(&determ, minor_row, minor_col);
        res.matrix_[minor_row][minor_col] = determ.Determinant();
      }
    }
    algebraic_complements(&res);
  }

  return res;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::invalid_argument("Invalid argument, matrix is not square");
  }

  S21Matrix copy(*this);

  int not_zero = 0, all_zero = 1, sign = 1;
  double* tmp = NULL;
  double result = 1.;

  for (int rows = 0; (rows < copy.rows_ - 1) && result; rows++) {
    for (not_zero = rows, all_zero = 1; not_zero < copy.rows_; not_zero++) {
      if (fabs(copy.matrix_[not_zero][rows])) {
        tmp = copy.matrix_[not_zero];
        all_zero = 0;
        break;
      }
    }

    if (copy.matrix_[rows][rows] == 0. && tmp && !all_zero) {
      copy.matrix_[not_zero] = copy.matrix_[rows];
      copy.matrix_[rows] = tmp;
      sign = -sign;
    }

    if (!all_zero) {
      for (int m = rows + 1; m < copy.rows_; m++) {
        if (copy.matrix_[m][rows] && copy.matrix_[rows][rows]) {
          double tmp = copy.matrix_[m][rows] / copy.matrix_[rows][rows];
          for (int n = rows; n < copy.cols_; n++) {
            copy.matrix_[m][n] -= copy.matrix_[rows][n] * tmp;
          }
        }
      }
    } else {
      result = 0.;
    }
  }

  for (int k = 0; k < rows_; k++) {
    result *= copy.matrix_[k][k];
  }
  result *= sign;

  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determ = 0.;
  determ = Determinant();
  S21Matrix result(*this);

  if (fabs(determ) > 1e-7) {
    if (rows_ == 1) {
      result.matrix_[0][0] = 1 / matrix_[0][0];
      return result;
    }

    result = result.CalcComplements();
    result = result.Transpose();

    determ = 1 / determ;
    result.MulNumber(determ);
  } else {
    throw std::runtime_error(
        "Error: determinant is zero, inverse matrix doesn't exist");
  }
  return result;
}

void S21Matrix::create_determinant_matrix(S21Matrix* determ, int minor_row,
                                          int minor_col) {
  for (int src_row = 0, det_rows = 0; src_row < rows_; src_row++) {
    for (int src_col = 0, det_col = 0; src_col < cols_; src_col++) {
      if (src_row != minor_row && src_col != minor_col) {
        determ->matrix_[det_rows][det_col] = matrix_[src_row][src_col];
        det_col++;

        if (det_col == determ->cols_) {
          det_rows++;
          det_col = 0;
        }
      }
    }
  }
}

void S21Matrix::algebraic_complements(S21Matrix* result) {
  for (int r = 0; r < rows_; r++) {
    for (int c = 0; c < cols_; c++) {
      if ((r + c) % 2) {
        result->matrix_[r][c] *= -1;
      }
    }
  }
}

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows == rows_) {
    return;
  }

  double** newMatrix;

  AlocMatrix(rows, cols_, &newMatrix);

  for (int i = 0; i < rows_ && i < rows; i++) {
    for (int j = 0; j < cols_; j++) {
      newMatrix[i][j] = matrix_[i][j];
    }
  }

  DelMatrix(matrix_);
  rows_ = rows;
  matrix_ = newMatrix;
}

void S21Matrix::SetCols(int cols) {
  if (cols == cols_) {
    return;
  }

  double** newMatrix;

  AlocMatrix(rows_, cols, &newMatrix);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_ && j < cols; j++) {
      newMatrix[i][j] = matrix_[i][j];
    }
  }

  DelMatrix(matrix_);
  cols_ = cols;
  matrix_ = newMatrix;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  S21Matrix(other.rows_, other.cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix res = *this;

  res.SumMatrix(other);

  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix res = *this;

  res.SubMatrix(other);

  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res = *this;

  res.MulMatrix(other);

  return res;
}

S21Matrix S21Matrix::operator*(const double& num) {
  S21Matrix res = *this;

  res.MulNumber(num);

  return res;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}

void S21Matrix::operator+=(const S21Matrix& other) { SumMatrix(other); }

void S21Matrix::operator-=(const S21Matrix& other) { SubMatrix(other); }

void S21Matrix::operator*=(const S21Matrix& other) { MulMatrix(other); }

void S21Matrix::operator*=(const double& num) { MulNumber(num); }

double& S21Matrix::operator()(int i, int j) {
  if (i < 0 || i > rows_ - 1) {
    throw std::out_of_range("Invalid argument: i - is out of range");
  }

  if (j < 0 || j > cols_ - 1) {
    throw std::out_of_range("Invalid argument: j - is out of range");
  }

  return matrix_[i][j];
}
