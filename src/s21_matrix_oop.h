#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_;
  bool isValid() const;
  void CalcMinor(S21Matrix& m, int i, int j);

 public:
  S21Matrix();
  S21Matrix(const int rows, const int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  S21Matrix& operator=(const S21Matrix& other);
  bool operator==(const S21Matrix& other);
  double& operator()(int i, int j);

  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  int GetRows();
  int GetCols();
  void SetRows(const int rows);
  void SetCols(const int cols);
};

S21Matrix operator+(const S21Matrix& left, const S21Matrix& right);
S21Matrix operator-(const S21Matrix& left, const S21Matrix& right);
S21Matrix operator*(const S21Matrix& left, const S21Matrix& right);
S21Matrix operator*(const S21Matrix& m, const double num);
bool operator==(const S21Matrix& left, const S21Matrix& right);

#endif  // S21_MATRIX_OOP_H