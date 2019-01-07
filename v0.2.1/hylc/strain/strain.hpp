#ifndef _STRAIN_HPP_
#define _STRAIN_HPP_

#include "../../vectors.hpp"
#include "../MathematicaDefinitions.h"

#include <utility>
#include <tuple>
#include <vector>


// (xlocal)(\d+) -> $1($2)
// (invDm)(\d)(\d) -> $1($2,$3)

namespace hylc {
  typedef Mat<18, 18> Mat18x18;
  typedef Vec<18> Vec18;
  typedef Vec<6> Vec6;
  typedef Mat<6, 18> Mat6x18;
namespace mathematica {
Vec6 ek(const Vec18 &xlocal, const Mat2x2 &invDm, const Real &A, const Real &l0,
        const Real &l1, const Real &l2, const Vec2 &t0, const Vec2 &t1,
        const Vec2 &t2);

std::tuple<std::vector<Mat18x18>, Mat6x18, Vec6>
ek_valdrv(const Vec18 &xlocal, const Mat2x2 &invDm, const Real &A, const Real &l0,
      const Real &l1, const Real &l2, const Vec2 &t0, const Vec2 &t1,
      const Vec2 &t2);
} // namespace mathematica
} // namespace hylc

#endif /* end of include guard: _STRAIN_HPP_ */
