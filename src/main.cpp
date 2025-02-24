#include "s21_matrix_oop.h"

void fillMatrix(double* num, S21Matrix& other) {
  int n = 0;

  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      other(i, j) = num[n++];
    }
  }
}

int main() {
  S21Matrix basic(3, 3);
  double seq[] = {11.01, 22.02, 33.03, 44.04, 55.05,
                  66.06, 77.07, 88.08, 99.09};
  fillMatrix(seq, basic);

  S21Matrix sum(3, 3);
  sum = basic + basic;
  basic = basic * 2.0;

  for (int i = 0; i < 3; i++) {
    sum += basic;
    sum *= 2.0;
    sum *= 0.5;
    sum -= basic;
  }

  sum = sum * sum;
  basic = basic * basic;

  //   basic = basic - sum;
  sum = sum - basic;

  for (int i = 0; i < 3; i++) {
    for (int a = 0; a < 3; a++) {
      printf("%lf ", (basic)(i, a));
    }
    printf("\n");
  }

  for (int i = 0; i < 3; i++) {
    for (int a = 0; a < 3; a++) {
      printf("%lf ", (sum)(i, a));
    }
    printf("\n");
  }

  return 0;
}