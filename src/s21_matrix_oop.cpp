#include "s21_matrix_oop.h"

#include <cmath>
#include <iostream>

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(const int rows, const int cols) {
  if (rows < 1 || cols < 1)
    throw std::invalid_argument("Invalid number of cols or rows!");
  rows_ = rows;
  cols_ = cols;

  matrix_ = new double*[rows_]();
  if (matrix_ == nullptr) throw "Could not allocate memory!";
  for (int i = 0; i < rows_; i++) {
    this->matrix_[i] = new double[cols_]();
    if (matrix_[i] == nullptr) throw "Could not allocate memory!";
  }
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] = 0;
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double*[rows_]();
  if (matrix_ == nullptr) throw "Could not allocate memory!";
  for (int i = 0; i < rows_; i++) {
    this->matrix_[i] = new double[cols_]();
    if (matrix_[i] == nullptr) throw "Could not allocate memory!";
  }

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] = other.matrix_[i][j];
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept {
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = other.matrix_;

  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  if (isValid()) {
    for (int i = 0; i < rows_; i++) delete[] matrix_[i];
    delete[] matrix_;
  }
  rows_ = cols_ = 0;
  matrix_ = nullptr;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  if (!other.isValid()) throw std::logic_error("Incorrect Matrix!");

  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  if (!other.isValid()) throw std::logic_error("Incorrect Matrix!");

  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  if (!other.isValid()) throw std::logic_error("Incorrect Matrix!");

  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& other) {
  return this->EqMatrix(other);
}

bool operator==(const S21Matrix& left, const S21Matrix& right) {
  S21Matrix copy{left};
  return copy.EqMatrix(right);
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  this->~S21Matrix();
  rows_ = other.rows_;
  cols_ = other.cols_;

  matrix_ = new double*[rows_]();
  if (matrix_ == nullptr) throw "Could not allocate memory!";
  for (int i = 0; i < rows_; i++) {
    this->matrix_[i] = new double[cols_]();
    if (matrix_[i] == nullptr) throw "Could not allocate memory!";
  }

  for (int i = 0; i < other.rows_; i++) {
    for (int j = 0; j < other.cols_; j++) matrix_[i][j] = other.matrix_[i][j];
  }
  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0)
    throw std::logic_error("Out of range!");
  return matrix_[i][j];
}

S21Matrix operator+(const S21Matrix& left, const S21Matrix& right) {
  S21Matrix res;
  res = left;
  res += right;
  return res;
}

S21Matrix operator-(const S21Matrix& left, const S21Matrix& right) {
  S21Matrix res;
  res = left;
  res -= right;
  return res;
}

S21Matrix operator*(const S21Matrix& left, const S21Matrix& right) {
  S21Matrix res;
  res = left;
  res *= right;
  return res;
}

S21Matrix operator*(const S21Matrix& m, const double num) {
  S21Matrix res;
  res = m;
  res *= num;
  return res;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool res = true;
  int flag = 0;

  if (this->cols_ != other.cols_ || this->rows_ != other.rows_ ||
      !this->isValid() || !other.isValid() || this->cols_ <= 0 ||
      this->rows_ <= 0 || other.cols_ <= 0 || other.rows_ <= 0)
    res = false;
  else {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        if (fabs(this->matrix_[i][j] - other.matrix_[i][j]) > 1e-07) {
          res = false;
          flag = 1;
          break;
        }
        if (flag) break;
      }
    }
  }

  return res;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (!this->isValid() || !other.isValid() || this->rows_ <= 0 ||
      this->cols_ <= 0 || other.rows_ <= 0 || other.cols_ <= 0)
    throw std::logic_error("Matrices are not valid!");
  else if (this->cols_ != other.cols_ || this->rows_ != other.rows_)
    throw std::logic_error("Matrices are not size equal!");
  else {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] += other.matrix_[i][j];
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (!this->isValid() || !other.isValid() || this->rows_ <= 0 ||
      this->cols_ <= 0 || other.rows_ <= 0 || other.cols_ <= 0)
    throw std::logic_error("Matrices are not valid!");
  else if (this->cols_ != other.cols_ || this->rows_ != other.rows_)
    throw std::logic_error("Matrices are not size equal!");
  else {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] -= other.matrix_[i][j];
  }
}

void S21Matrix::MulNumber(const double num) {
  if (!this->isValid() || this->rows_ <= 0 || this->cols_ <= 0)
    throw std::logic_error("Matrix is not valid!");
  else {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++) this->matrix_[i][j] *= num;
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (!this->isValid() || !other.isValid() || this->rows_ <= 0 ||
      this->cols_ <= 0 || other.rows_ <= 0 || other.cols_ <= 0)
    throw std::logic_error("Matrices are not valid!");
  else if (this->cols_ != other.rows_)
    throw std::logic_error(
        "Number of 1st matrix cols is not equal number of 2nd matrix rows!");
  else {
    S21Matrix res = S21Matrix(rows_, other.cols_);

    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < other.cols_; j++) {
        res.matrix_[i][j] = 0;
        for (int k = 0; k < this->cols_; k++)
          res.matrix_[i][j] += this->matrix_[i][k] * other.matrix_[k][j];
      }

    *this = res;
  }
}

S21Matrix S21Matrix::Transpose() {
  if (!isValid() || rows_ <= 0 || cols_ <= 0)
    throw std::logic_error("Invalid matrix!");

  S21Matrix res(cols_, rows_);

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) res.matrix_[j][i] = matrix_[i][j];

  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (!isValid() || rows_ <= 0 || cols_ <= 0)
    throw std::logic_error("Invalid matrix!");
  else if (rows_ != cols_)
    throw std::logic_error("Matrix is not square!");

  S21Matrix res(rows_, cols_);
  if (rows_ == 1)
    res.matrix_[0][0] = 1;
  else {
    double det;
    S21Matrix temp(rows_ - 1, cols_ - 1);

    for (int i = 0; i < rows_; i++)
      for (int j = 0; j < cols_; j++) {
        CalcMinor(temp, i, j);
        det = temp.Determinant();
        if ((i + j) % 2 == 1)
          res.matrix_[i][j] = -det;
        else
          res.matrix_[i][j] = det;
      }
  }
  return res;
}

double S21Matrix::Determinant() {
  double res;
  if (!isValid() || rows_ <= 0 || cols_ <= 0)
    throw std::logic_error("Invalid matrix!");
  else if (rows_ != cols_)
    throw std::logic_error("Matrix is not square!");
  else {
    if (rows_ == 1)
      res = matrix_[0][0];
    else if (rows_ == 2)
      res = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
    else {
      res = 0;
      double det;
      S21Matrix temp(rows_ - 1, cols_ - 1);

      for (int i = 0; i < rows_; i++) {
        CalcMinor(temp, 0, i);
        if (i % 2 == 0) {
          det = temp.Determinant();
          res += matrix_[0][i] * det;
        } else {
          det = temp.Determinant();
          res -= matrix_[0][i] * det;
        }
      }
    }
  }
  return res;
}

void S21Matrix::CalcMinor(S21Matrix& m, int i, int j) {
  if (!this->isValid() || this->rows_ <= 0 || this->cols_ <= 0 ||
      !m.isValid() || m.rows_ <= 0 || m.cols_ <= 0)
    throw std::logic_error("Invalid matrix!");
  else {
    m.matrix_[0][0] = this->matrix_[0][0];

    int row_t = 0, col_t = 0;
    for (int k = 0; k < this->rows_; k++)
      for (int l = 0; l < this->cols_; l++) {
        if (k != i && l != j) {
          m.matrix_[row_t][col_t] = this->matrix_[k][l];
          if (col_t != m.cols_ - 1)
            col_t++;
          else {
            row_t++;
            col_t = 0;
          }
        }
      }
  }
}

S21Matrix S21Matrix::InverseMatrix() {
  if (!isValid() || rows_ <= 0 || cols_ <= 0)
    throw std::logic_error("Invalid matrix!");
  else if (rows_ != cols_)
    throw std::logic_error("Matrix is not square!");

  double det;
  det = Determinant();
  S21Matrix res;

  if (det == 0)
    throw std::logic_error("Determinant equal 0!");
  else {
    if (rows_ == 1) {
      res = S21Matrix(1, 1);
      res.matrix_[0][0] = 1.0 / matrix_[0][0];
    } else {
      res = CalcComplements();
      res = res.Transpose();
      for (int i = 0; i < res.rows_; i++)
        for (int j = 0; j < res.cols_; j++) res.matrix_[i][j] /= det;
    }
  }
  return res;
}

bool S21Matrix::isValid() const { return matrix_ != nullptr; }

int S21Matrix::GetRows() { return rows_; }

int S21Matrix::GetCols() { return cols_; }

void S21Matrix::SetRows(const int rows) {
  if (rows != rows_) {
    S21Matrix temp(rows, cols_);
    if (rows > rows_)
      for (int i = 0; i < rows_; i++)
        for (int j = 0; j < cols_; j++) temp(i, j) = matrix_[i][j];
    else
      for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols_; j++) temp(i, j) = matrix_[i][j];
    *this = temp;
  }
}

void S21Matrix::SetCols(const int cols) {
  if (cols != cols_) {
    S21Matrix temp(rows_, cols);
    if (cols > cols_)
      for (int i = 0; i < rows_; i++)
        for (int j = 0; j < cols_; j++) temp(i, j) = matrix_[i][j];
    else
      for (int i = 0; i < rows_; i++)
        for (int j = 0; j < cols; j++) temp(i, j) = matrix_[i][j];
    *this = temp;
  }
}