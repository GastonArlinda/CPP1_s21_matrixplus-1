# CPP1_s21_matrixplus-1

## Описание проекта
**CPP1_s21_matrixplus-1** – учебный проект, направленный на реализацию библиотеки для работы с матрицами на языке C++ с применением объектно-ориентированного программирования.  
Библиотека реализует основные операции над матрицами и перегружает соответствующие операторы.  

Проект выполнен в соответствии с:
- Стандартом C++17
- Google Style
- Принципами объектно-ориентированного программирования (ООП)

## Реализованные функции
Все методы разработаны в классе `S21Matrix` и соответствуют стандартным операциям с матрицами:

- `bool EqMatrix(const S21Matrix& other)` – проверка на равенство матриц
- `void SumMatrix(const S21Matrix& other)` – сложение матриц
- `void SubMatrix(const S21Matrix& other)` – вычитание матриц
- `void MulNumber(const double num)` – умножение матрицы на число
- `void MulMatrix(const S21Matrix& other)` – умножение матриц
- `S21Matrix Transpose()` – транспонирование матрицы
- `S21Matrix CalcComplements()` – вычисление алгебраических дополнений
- `double Determinant()` – вычисление определителя
- `S21Matrix InverseMatrix()` – вычисление обратной матрицы

## Перегруженные операторы
- `+` – сложение матриц  
- `-` – вычитание матриц  
- `*` – умножение матриц и умножение на число  
- `==` – проверка на равенство  
- `=` – присваивание  
- `+=` – сложение с присваиванием  
- `-=` – вычитание с присваиванием  
- `*=` – умножение с присваиванием  
- `()` – доступ к элементу матрицы по индексу  

## Сборка и использование
### Сборка библиотеки:
```sh
make
```
### Тесты
```sh
make tests
```
### Проверка покрытия
```sh
make gcov_report
```
