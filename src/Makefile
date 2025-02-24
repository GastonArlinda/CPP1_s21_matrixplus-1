CFLAGS=-Wall -Werror -Wextra --std=c++17 -lm
GCOVFLAGS=--coverage -fprofile-arcs -ftest-coverage 

all: clean s21_matrix_oop.a

s21_matrix_oop.a:
	g++ $(CFLAGS) $(GCOVFLAGS) s21_matrix_oop.cpp -c
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

tests: clean s21_matrix_oop.a
	g++ tests.cpp s21_matrix_oop.a -o test $(CFLAGS) -lgtest -lgcov -o tests
	./tests

gcov_report: tests
	gcov *.gcda
	gcovr
	gcovr --html-details -o report.html

clean:
	rm -rf *.o *.a *.gcda *.gcno *.gcov *.gch rep.info *.html *.css tests

format:
	clang-format -i *.cpp *.h