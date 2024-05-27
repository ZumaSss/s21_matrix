CC = g++
FLAGS = -Wall -Wextra -Werror -std=c++17

LDFLAGS = -lcheck -lsubunit -lrt -lpthread --std=c++17 -lstdc++ -lm -lgtest -lgcov
TEST_FLAGS = --coverage $(LDFLAGS)
OPEN_CMD = xdg-open


GCOV_DIR = report

LIB_NAME = s21_matrix_oop

SRCS=$(wildcard s21_*.cpp)
OBJS=$(SRCS:.c=.o)

all: $(LIB_NAME).a test clean

$(LIB_NAME).a:
	$(CC) $(FLAGS) -c s21_matrix_oop.cpp
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

test: clean $(LIB_NAME).a
	$(CC) tests.cpp s21_matrix_oop.a $(LDFLAGS) -o test
	./test

gcov_report:
	make clean
	mkdir -p $(GCOV_DIR)
	$(CC) tests.cpp s21_matrix_oop.cpp -o greport $(TEST_FLAGS)
	./greport
	lcov  -t "s21_matrix_oop" -o s21_matrix_oop.info -c -d . --no-external
	genhtml -o ./$(GCOV_DIR) s21_matrix_oop.info
	open ./$(GCOV_DIR)/index.html


clean:
	@-rm -rf $(OBJ_DIR)/*.o
	@-rm -rf $(EXECUTABLE)
	@-rm -rf *.o && rm -rf *.gcno
	@-rm -rf *.a && rm -rf *.gcda
	@-rm -rf *.info && rm -rf *.gcov
	@-rm -rf ./test && rm -rf ./gcov_report
	@-rm -rf ./report/
	@-rm -rf s21_test greport
	@-rm -rf RESULT_VALGRIND.txt



style: clean
	@cp ../materials/linters/.clang-format ../src/.clang-format
	@clang-format -n *.h *.cpp
	@rm .clang-format

dev:
	$(CC) $(FLAGS) -g -fsanitize=address $(SRCS)
	./a.out