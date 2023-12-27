
// Copyright 2022 Miheev Ivan

#include <iostream>
#include <random>
#include <vector>

class SprMatCCS {
 private:
  int dim;                  // number of matrix dimension
  int cap;                  // number of elements other than zero
  std::vector<double> val;  // vector that contains values of elements
  std::vector<int> rows;    // vector that contains position in a column
  std::vector<int> ptr;     // vector of pointers of positions

 public:
  SprMatCCS(int _dim = 0, int _cap = 0,
            const std::vector<double>& _val = std::vector<double>(),
            const std::vector<int>& _rows = std::vector<int>(),
            const std::vector<int>& _ptr = std::vector<int>())
      : dim(_dim), cap(_cap), val(_val), rows(_rows), ptr(_ptr) {}
  SprMatCCS(const SprMatCCS& mat)
      : dim(mat.dim),
        cap(mat.cap),
        val(mat.val),
        rows(mat.rows),
        ptr(mat.ptr) {}
  ~SprMatCCS();

  int getCapacity() { return this->cap; }
  int getDimension() { return this->dim; }

  std::vector<double> getValues() { return this->val; }
  std::vector<int> getRows() { return this->rows; }
  std::vector<int> getPtr() { return this->ptr; }

  SprMatCCS& randMat(int size, int cnt = -1);
  SprMatCCS& clrMat();  // function to clear Matrix
  SprMatCCS transMat();

  SprMatCCS& operator=(const SprMatCCS& mat);
  SprMatCCS operator*(SprMatCCS mat);
  bool operator==(SprMatCCS mat);
  bool operator!=(const SprMatCCS& mat) { return !(*this == mat); }

  void shwVal();
};

