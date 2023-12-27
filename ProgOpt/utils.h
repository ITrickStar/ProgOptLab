#ifndef UTILS_H
#define UTILS_H

#include <iostream>

inline bool isZero(const double num) { return (std::abs(num) < 0.00000001); }
inline bool isEqual(double x, double y) { return std::fabs(x - y) < 0.00000001; }

#endif  // UTILS_H