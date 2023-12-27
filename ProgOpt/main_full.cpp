#include <omp.h>

#include <iostream>
#include <random>
#include <vector>

// Utils.h
inline bool isZero(const double num) { return (std::abs(num) < 0.00000001); }
inline bool isEqual(double x, double y) {
  return std::fabs(x - y) < 0.00000001;
}

// SprMatCCSOMP
class SprMatCCSOMP {
 private:
  int dim;                  // number of matrix dimension
  int cap;                  // number of elements other than zero
  std::vector<double> val;  // vector that contains values of elements
  std::vector<int> rows;    // vector that contains position in a column
  std::vector<int> ptr;     // vector of pointers of positions

 public:
  SprMatCCSOMP(int _dim = 0, int _cap = 0,
               const std::vector<double>& _val = std::vector<double>(),
               const std::vector<int>& _rows = std::vector<int>(),
               const std::vector<int>& _ptr = std::vector<int>())
      : dim(_dim), cap(_cap), val(_val), rows(_rows), ptr(_ptr) {}
  SprMatCCSOMP(const SprMatCCSOMP& mat)
      : dim(mat.dim),
        cap(mat.cap),
        val(mat.val),
        rows(mat.rows),
        ptr(mat.ptr) {}
  ~SprMatCCSOMP() {
    this->val.clear();
    this->rows.clear();
    this->ptr.clear();
  };

  int getCapacity() { return this->cap; }
  int getDimension() { return this->dim; }

  std::vector<double> getValues() { return this->val; }
  std::vector<int> getRows() { return this->rows; }
  std::vector<int> getPtr() { return this->ptr; }

  SprMatCCSOMP& randMat(int size, int spr = -1) {
    std::random_device rd;
    std::mt19937 gen(rd());
    if (spr <= 0) {
      spr = 1;
    }

    this->clrMat();
    this->dim = size;
    this->cap = spr * size;
    this->val.resize(spr * size);
    this->rows.resize(spr * size);
    this->ptr.resize(size + 1);

    for (int i = 0; i < size; i++) {
      for (int j = 0; j < spr; j++) {
        bool b;
        do {
          this->rows[i * spr + j] = gen() % size + 1;
          b = false;
          for (int k = 0; k < j; k++) {
            if (this->rows[i * spr + j] == this->rows[i * spr + k]) {
              b = true;
            }
          }
        } while (b);
      }

      for (int k = 0; k < spr - 1; k++) {
        if (this->rows[i * spr + k] > this->rows[i * spr + k + 1]) {
          int tmp = this->rows[i * spr + k];
          this->rows[i * spr + k] = this->rows[i * spr + k + 1];
          this->rows[i * spr + k + 1] = tmp;
        }
      }
    }

    int rng = 100;  // upper bound of values
    for (int i = 0; i < spr * size; i++) {
      this->val[i] = (int)gen() % rng;
    }

    int sum = 1;
    for (int i = 0; i < size + 1; i++) {
      this->ptr[i] = sum;
      sum += spr;
    }

    return *this;
  };
  SprMatCCSOMP& clrMat() {
    this->cap = 0;
    this->dim = 0;
    this->val.clear();
    this->rows.clear();
    this->ptr.clear();

    return *this;
  };  // function to clear Matrix
  SprMatCCSOMP transMat() {
    SprMatCCSOMP res;
    res.dim = this->dim;
    res.cap = this->cap;
    res.val.resize(this->cap);
    res.rows.resize(this->cap);
    res.ptr.resize(this->dim + 1);

    for (int i = 0; i < this->cap; i++) {
      res.ptr[this->rows[i] - 1]++;
    }

    int tmp, sum = 1;
    for (int i = 0; i < this->dim + 1; i++) {
      tmp = res.ptr[i];
      res.ptr[i] = sum;
      sum += tmp;
    }
    std::vector<int> _ptr = res.ptr;
    for (int i = 0; i < this->dim; i++) {
      for (int j = this->ptr[i]; j < this->ptr[i + 1]; j++) {
        int r_ind = this->rows[j - 1];
        int i_ind = _ptr[r_ind - 1];
        res.rows[i_ind - 1] = i + 1;
        res.val[i_ind - 1] = this->val[j - 1];
        _ptr[r_ind - 1]++;
      }
    }
    return res;
  };

  SprMatCCSOMP& operator=(const SprMatCCSOMP& mat) {
    this->cap = mat.cap;
    this->dim = mat.dim;
    this->val = mat.val;
    this->rows = mat.rows;
    this->ptr = mat.ptr;

    return *this;
  };
  SprMatCCSOMP operator*(SprMatCCSOMP mat) {
    if (this->dim != mat.dim) throw "wrong sizes";
    SprMatCCSOMP trans = this->transMat();

    std::vector<double> val_res;
    std::vector<int> rows_res;
    std::vector<int> ptr_res{1};

    int capacity = 1;
    double sum;
    for (int i = 0; i < mat.dim; i++) {
      std::vector<int> ip(trans.dim, 0);
      for (int j = mat.ptr[i]; j < mat.ptr[i + 1]; j++) {
        ip[mat.rows[j - 1] - 1] = j;
      }
      for (int j = 0; j < trans.dim; j++) {
        sum = 0;
        for (int k = trans.ptr[j]; k < trans.ptr[j + 1]; k++) {
          int p = ip[trans.rows[k - 1] - 1];
          if (p) {
            sum += mat.val[p - 1] * trans.val[k - 1];
          }
        }
        if (!isZero(sum)) {
          rows_res.push_back(j + 1);
          val_res.push_back(sum);
          capacity++;
        }
      }
      ptr_res.push_back(capacity);
    }
    SprMatCCSOMP res(trans.dim, capacity - 1, val_res, rows_res, ptr_res);
    return res;
  };
  bool operator==(SprMatCCSOMP mat) {
    if (this->dim != mat.dim || this->cap != mat.cap ||
        this->rows != mat.rows || this->ptr != mat.ptr)
      return false;
    for (int i = 0; i < this->cap; i++)
      if (!isEqual(val[i], mat.val[i])) return false;
    return true;
  };
  bool operator!=(const SprMatCCSOMP& mat) { return !(*this == mat); }

  SprMatCCSOMP ParallelMult(SprMatCCSOMP mat) {
    if (this->dim != mat.dim) throw "wrong sizes";
    SprMatCCSOMP trans = this->transMat();

    const int num_threads = 7;

    std::vector<std::vector<double>> val_shared(num_threads);
    std::vector<std::vector<int>> rows_shared(num_threads);
    std::vector<int> ptr_counter(trans.dim);

#pragma omp parallel num_threads(num_threads)
    {
      int ind = omp_get_thread_num();

#pragma omp for
      for (int i = 0; i < trans.dim; i++) {
        std::vector<int> ip(trans.dim, 0);
        int cnt = 0;
        for (int j = mat.ptr[i]; j < mat.ptr[i + 1]; j++) {
          ip[mat.rows[j - 1] - 1] = j;
        }
        for (int j = 0; j < trans.dim; j++) {
          double sum = 0;
          for (int k = trans.ptr[j]; k < trans.ptr[j + 1]; k++) {
            int irow = trans.rows[k - 1];
            int p = ip[irow - 1];
            if (p) {
              sum += mat.val[p - 1] * trans.val[k - 1];
            }
          }
          if (!isZero(sum)) {
            val_shared[ind].push_back(sum);
            rows_shared[ind].push_back(j + 1);
            cnt++;
          }
        }
        ptr_counter[i] += cnt;
      }
    }

    std::vector<double> val_res;
    std::vector<int> row_res;
    std::vector<int> ptr_res{1};

    for (int i = 0; i < num_threads; i++) {
      val_res.insert(val_res.end(), val_shared[i].begin(), val_shared[i].end());
      row_res.insert(row_res.end(), rows_shared[i].begin(),
                     rows_shared[i].end());
    }
    int sum = 1;
    for (int i = 0; i < trans.dim; i++) {
      sum += ptr_counter[i];
      ptr_res.push_back(sum);
    }

    SprMatCCSOMP res(trans.dim, val_res.size(), val_res, row_res, ptr_res);
    return res;
  };

  void shwVal() {
    for (int i = 0; i < this->cap; i++) {
      std::cout << val[i] << " ";
    }
    std::cout << "\n";
    for (int i = 0; i < this->cap; i++) {
      std::cout << rows[i] << " ";
    }
    std::cout << "\n";
    for (int i = 0; i < static_cast<int>(ptr.size()); i++) {
      std::cout << ptr[i] << " ";
    }
    std::cout << "\n";
  };
};

// MAIN
int main() {
  clock_t start, end;
  int dim = 15000;
  int spr = 50;
  SprMatCCSOMP A, B;

  start = clock();
  A.randMat(dim, spr);
  B.randMat(dim, spr);
  end = clock();
  std::cout << "Time required -> " << (end - start + .0) / CLOCKS_PER_SEC
            << " <-" << std::endl;

  start = clock();
  A* B;
  end = clock();
  std::cout << "Time required Sequ-> " << (end - start + .0) / CLOCKS_PER_SEC
            << " <-" << std::endl;

  start = clock();
  A.ParallelMult(A);
  end = clock();

  std::cout << "Time required OMP -> " << (end - start + .0) / CLOCKS_PER_SEC
            << " <-" << std::endl;
  return 0;
};