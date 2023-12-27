#include "matrix_ccs.h"
#include "matrix_ccs_omp.h"

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
