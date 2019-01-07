#ifndef _MATHEMATICADEFINITIONS_H_
#define _MATHEMATICADEFINITIONS_H_

#include <math.h>

namespace hylc {
namespace mathematica {

typedef double Real;
inline Real Power(const Real &A, const Real &B) { return std::pow(A, B); }
// inline Real Power(const Real &A, const int &B) { return std::pow(A, B); }
inline Real Power(const Real &A, const int &B) {
  if (abs(B) > 8) // arbitary threshold
    return std::pow(A, B);
  Real out = 1.0;
  if (B >= 0) {
    for (int i = 0; i < B; ++i)
      out *= A;
  } else {
    Real invA = 1.0 / A;
    for (int i = 0; i < -B; ++i)
      out *= invA;
  }
  return out;
}
inline Real Abs(const Real &A) { return std::abs(A); }
inline Real Sqrt(const Real &A) { return std::sqrt(A); }
inline Real Exp(const Real &A) { return std::exp(A); }
inline Real Sin(const Real &A) { return std::sin(A); }
inline Real Cos(const Real &A) { return std::cos(A); }
inline Real ArcCos(const Real &A) { return std::acos(A); }
inline Real ArcTan(const Real &X, const Real &Y) { return std::atan2(Y, X); }
inline Real Tan(const Real &X) { return std::tan(X); }
inline Real Sec(const Real &X) { return 1.0/std::cos(X); }

} // namespace Mathematica
} // namespace HYLC

#endif /* end of include guard: _MATHEMATICADEFINITIONS_H_ */
