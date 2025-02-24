#include "gtest/gtest.h"
#include "s21_matrix_oop.h"

void fillMatrix(double* num, S21Matrix& other) {
  int n = 0;

  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      other(i, j) = num[n++];
    }
  }
}

TEST(Matrix, Constructors) {
  S21Matrix first;
  first(0, 0) = 73.895;

  S21Matrix second;
  second(0, 0) = 73.895;

  EXPECT_EQ(first, second);
}

TEST(Matrix, Copy) {
  S21Matrix source(2, 2);

  double num[] = {11.11, 12.12, 21.21, 31.31};
  fillMatrix(num, source);

  S21Matrix copy(source);
  EXPECT_EQ(source, copy);
}

TEST(Matrix, Move) {
  S21Matrix source(2, 2);

  double num[] = {11.11, 12.12, 21.21, 31.31};
  fillMatrix(num, source);

  S21Matrix copy(source);
  EXPECT_EQ(source, copy);

  S21Matrix move(std::move(source));
  EXPECT_FALSE(move == source);
}

TEST(Matrix, SetterGetter) {
  S21Matrix right(123, 911);

  EXPECT_EQ(right.GetCols(), 911);
  EXPECT_EQ(right.GetRows(), 123);

  right.SetRows(2);
  EXPECT_EQ(right.GetRows(), 2);

  right.SetCols(20);
  EXPECT_EQ(right.GetCols(), 20);

  right.SetRows(2);
  right.SetCols(20);
}

TEST(Matrix, Operations) {
  S21Matrix basic(3, 3);
  double num[] = {1, 2, 4, 3, 3, 5, 2, 4, 4};
  fillMatrix(num, basic);

  S21Matrix one(3, 3);
  double num_one[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  fillMatrix(num_one, one);

  S21Matrix calc(3, 3);
  calc = basic.InverseMatrix();
  EXPECT_TRUE(calc == basic.InverseMatrix());

  calc.MulMatrix(basic);
  EXPECT_EQ(one, calc);

  S21Matrix one_index = S21Matrix();
  one_index(0, 0) = 8121980;
  one_index.CalcComplements();
  EXPECT_EQ(one_index(0, 0), 8121980);

  one_index.InverseMatrix();
  EXPECT_EQ(one_index(0, 0), 8121980);
}

TEST(Matrix, Determinant) {
  S21Matrix other = S21Matrix();
  other(0, 0) = 214214321.4325452;
  EXPECT_EQ(other.Determinant(), other(0, 0));

  other.SetCols(4);
  other.SetRows(4);
  double seq[] = {1, 2, 3, 5, 6, 1, 12, 3, 0, 11, 2, 33, 7, 68, 9, 71};

  fillMatrix(seq, other);
  double det = -6984.;

  EXPECT_DOUBLE_EQ(det, other.Determinant());

  other.SetCols(80);
  other.SetRows(80);
  EXPECT_EQ(other.Determinant(), 0.0);

  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      other(i, j) = (random() % 1000) * 0.01;
    }
  }

  det = other.Determinant();
  EXPECT_EQ(det, other.Determinant());

  S21Matrix zero_row(3, 3);
  double zer[] = {0, 2, 3, 4, 0, 6, 7, 8, 0};
  fillMatrix(zer, zero_row);
  EXPECT_EQ(zero_row.Determinant(), 180);
}

TEST(Matrix, Operators) {
  S21Matrix basic(3, 3);
  double seq[] = {11.01, 22.02, 33.03, 44.04, 55.05,
                  66.06, 77.07, 88.08, 99.09};
  fillMatrix(seq, basic);

  S21Matrix sum(3, 3);
  sum = basic + basic;
  basic = basic * 2.0;
  EXPECT_EQ(sum, basic);

  for (int i = 0; i < 3; i++) {
    sum += basic;
    sum *= 2.0;
    EXPECT_FALSE(sum == basic);
    sum *= 0.5;
    sum -= basic;
  }
  EXPECT_EQ(sum, basic);

  sum = sum * sum;
  basic = basic * basic;
  EXPECT_EQ(basic, sum);

  basic = basic - sum;
  sum = sum - basic;
  EXPECT_FALSE(basic == sum);
}

TEST(Matrix, Errors) {
  EXPECT_THROW(S21Matrix errors(-123, 0), std::invalid_argument);
  EXPECT_THROW(S21Matrix errors(1, 0), std::invalid_argument);

  S21Matrix errorLeft(1, 55);
  S21Matrix errorRight(2, 55);

  EXPECT_THROW(errorLeft -= errorRight, std::invalid_argument);
  EXPECT_THROW(errorLeft += errorRight, std::invalid_argument);
  EXPECT_THROW(errorLeft *= errorRight, std::invalid_argument);

  S21Matrix basic(4, 4);

  EXPECT_THROW(errorLeft.Determinant(), std::invalid_argument);
  EXPECT_THROW(errorLeft.CalcComplements(), std::invalid_argument);
  EXPECT_THROW(basic.InverseMatrix(), std::runtime_error);
  EXPECT_THROW(basic(-12, 1) = 1, std::out_of_range);
  EXPECT_THROW(basic(1, 123) = 14, std::out_of_range);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}