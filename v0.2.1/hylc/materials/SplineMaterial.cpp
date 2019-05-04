#include "SplineMaterial.hpp"
#include <memory>
#include <vector>

using namespace hylc;
using namespace fitpackpp;

SplineMaterial::SplineMaterial() {}

double SplineMaterial::psi(const Vec6 &strain) {
  double val = 0;
  if (!initialized)
    return val;

  // clamp / barrier
  // Vec6 X;
  // for (int i = 0; i < 6;
  //      i++) {  // TODO clamp only bending for 1D but all for 2D?
  //              // if (i > 2)
  //   X(i) =
  //       std::min(std::max(strain(i), this->strain_min(i)), this->strain_max(i));
  //   // X(i) = strain(i);
  //   // else
  //   // X(i) = strain(i);
  // }

  // if (use_barrier) {
  //   val += psi_barrier(strain, X);
  // }
  Vec6 X = strain;

  // normalize
  for (int i = 0; i < 6; i++)
    X(i)     = (X(i) - this->strainshift(i)) / this->strainscale(i);

  // const
  val += C0;

  // 1D
  for (auto &spline1d : splines_1d) {
    double x=X(spline1d.k);
    spline1d.clamp(x);
    val += spline1d.spline->eval(x);
  }

  // 2D splines old
  for (auto &spline2d : splines_2d) {
    val += spline2d.spline->eval(X(spline2d.k0), X(spline2d.k1));
  }

  // 2D polys
  for (auto &poly : polys_2d) {
    double x = X(poly.k0);
    double y = X(poly.k1);
    poly.clamp(x,y);
    val += poly.eval(x, y);
  }

  return val;
}

Vec6 SplineMaterial::psi_grad(const Vec6 &strain) {
  Vec6 grad(0);
  if (!initialized)
    return grad;

  // clamp / barrier
  // Vec6 X;
  // for (int i = 0; i < 6;
  //      i++) {  // TODO clamp only bending for 1D but all for 2D?
  //              // if (i > 2)
  //   X(i) =
  //       std::min(std::max(strain(i), this->strain_min(i)), this->strain_max(i));
  //   // X(i) = strain(i);
  //   // else
  //   // X(i) = strain(i);
  // }

  // if (use_barrier) {
  //   grad += grad_barrier(strain, X);
  // }
  Vec6 X = strain;

  // normalize
  for (int i = 0; i < 6; i++)
    X(i)     = (X(i) - this->strainshift(i)) / this->strainscale(i);

  // 1D
  for (auto &spline1d : splines_1d) {
    double x=X(spline1d.k);
    spline1d.clamp(x);
    grad(spline1d.k) +=
        spline1d.spline->der(x, 1) / this->strainscale(spline1d.k);
  }

  // 2D splines old
  for (auto &spline2d : splines_2d) {
    double x = X(spline2d.k0);
    double y = X(spline2d.k1);
    grad(spline2d.k0) +=
        spline2d.spline->der(x, y, 1, 0) / this->strainscale(spline2d.k0);
    grad(spline2d.k1) +=
        spline2d.spline->der(x, y, 0, 1) / this->strainscale(spline2d.k1);
  }

  // 2D polys
  for (auto &poly : polys_2d) {
    double x = X(poly.k0);
    double y = X(poly.k1);
    poly.clamp(x,y);
    grad(poly.k0) += poly.dx(x, y) / this->strainscale(poly.k0);
    grad(poly.k1) += poly.dy(x, y) / this->strainscale(poly.k1);
  }

  return grad;
}

std::pair<Mat6x6, Vec6> SplineMaterial::psi_drv(const Vec6 &strain) {
  Vec6 grad(0);
  Mat6x6 hess(0);
  if (!initialized)
    return std::make_pair(hess, grad);

  // // clamp / barrier
  // Vec6 X;
  // for (int i = 0; i < 6;
  //      i++) {  // TODO clamp only bending for 1D but all for 2D?
  //              // if (i > 2)
  //   X(i) =
  //       std::min(std::max(strain(i), this->strain_min(i)), this->strain_max(i));
  //   // X(i) = strain(i);
  //   // else
  //   // X(i) = strain(i);
  // }

  // if (use_barrier) {
  //   auto ghB = gradhess_barrier(strain, X);
  //   hess += ghB.first;
  //   grad += ghB.second;
  // }
  Vec6 X = strain;

  // normalize
  for (int i = 0; i < 6; i++)
    X(i)     = (X(i) - this->strainshift(i)) / this->strainscale(i);

  // 1D
  for (auto &spline1d : splines_1d) {
    // if (spline1d.k > 2)
    // continue;
    // if (spline1d.k == 4)
    // continue; // DEBUG

    double x=X(spline1d.k);
    spline1d.clamp(x);

    grad(spline1d.k) +=
        spline1d.spline->der(x, 1) / this->strainscale(spline1d.k);
    hess(spline1d.k, spline1d.k) +=
        spline1d.spline->der(x, 2) /
        (this->strainscale(spline1d.k) * this->strainscale(spline1d.k));
  }

  // 2D splines old
  for (auto &spline2d : splines_2d) {
    continue;
    // TODO CLAMP SHIT
    // if (!(spline2d.k0==0 && spline2d.k1 == 2))
    //   continue;
    // if (!(spline2d.k0<3 && spline2d.k1 <3))
    //   continue;
    double x = X(spline2d.k0);
    double y = X(spline2d.k1);
    grad(spline2d.k0) +=
        spline2d.spline->der(x, y, 1, 0) / this->strainscale(spline2d.k0);
    grad(spline2d.k1) +=
        spline2d.spline->der(x, y, 0, 1) / this->strainscale(spline2d.k1);

    hess(spline2d.k0, spline2d.k0) +=
        spline2d.spline->der(x, y, 2, 0) /
        (this->strainscale(spline2d.k0) * this->strainscale(spline2d.k0));
    double dxdy =
        spline2d.spline->der(x, y, 1, 1) /
        (this->strainscale(spline2d.k0) * this->strainscale(spline2d.k1));
    hess(spline2d.k0, spline2d.k1) += dxdy;
    hess(spline2d.k1, spline2d.k0) += dxdy;
    hess(spline2d.k1, spline2d.k1) +=
        spline2d.spline->der(x, y, 0, 2) /
        (this->strainscale(spline2d.k1) * this->strainscale(spline2d.k1));
  }

  // 2D polys
  for (auto &poly : polys_2d) {
    // continue; // DEBUG
    // if (poly.k0 == 4 || poly.k1 == 4)
    //   continue;
    // if (poly.k0 > 2 || poly.k1 > 2)
    //   continue;
    double x = X(poly.k0);
    double y = X(poly.k1);
    poly.clamp(x,y);
    grad(poly.k0) += poly.dx(x, y) / this->strainscale(poly.k0);
    grad(poly.k1) += poly.dy(x, y) / this->strainscale(poly.k1);

    hess(poly.k0, poly.k0) += poly.dxdx(x, y) / (this->strainscale(poly.k0) *
                                                 this->strainscale(poly.k0));
    double dxdy = poly.dxdy(x, y) /
                  (this->strainscale(poly.k0) * this->strainscale(poly.k1));
    hess(poly.k0, poly.k1) += dxdy;
    hess(poly.k1, poly.k0) += dxdy;
    hess(poly.k1, poly.k1) += poly.dydy(x, y) / (this->strainscale(poly.k1) *
                                                 this->strainscale(poly.k1));
  }

  // grad(0) = 1.0;
  return std::make_pair(hess, grad);
}

using namespace hylc::mathematica;
typedef Vec6 strain_type;
typedef double value_type;
typedef Vec6 grad_type;
typedef std::pair<Mat6x6, Vec6> gradhess_type;

#define BARRIERPOWER 6  // should be even integer!
value_type SplineMaterial::psi_barrier(const strain_type &straincpp,
                                       const strain_type &strainclamped) {
  Real copt22  = strainscale(0);
  Real copt30  = 1 / copt22;
  Real copt47  = strainscale(1);
  Real copt48  = 1 / copt47;
  Real copt65  = strainscale(2);
  Real copt66  = 1 / copt65;
  Real copt82  = strainscale(3);
  Real copt90  = 1 / copt82;
  Real copt101 = strainscale(4);
  Real copt102 = 1 / copt101;
  Real copt109 = strainscale(5);
  Real copt110 = 1 / copt109;
  return bscale * (Power(bspeed * copt30 * (-strainclamped(0) + straincpp(0)),
                         BARRIERPOWER) +
                   Power(bspeed * copt48 * (-strainclamped(1) + straincpp(1)),
                         BARRIERPOWER) +
                   Power(bspeed * copt66 * (-strainclamped(2) + straincpp(2)),
                         BARRIERPOWER) +
                   Power(bspeed * copt90 * (-strainclamped(3) + straincpp(3)),
                         BARRIERPOWER) +
                   Power(bspeed * copt102 * (-strainclamped(4) + straincpp(4)),
                         BARRIERPOWER) +
                   Power(bspeed * copt110 * (-strainclamped(5) + straincpp(5)),
                         BARRIERPOWER));
}

grad_type SplineMaterial::grad_barrier(const strain_type &straincpp,
                                       const strain_type &strainclamped) {
  Vec6 out(0);
  Real copt22  = strainscale(0);
  Real copt30  = 1 / copt22;
  Real copt50  = strainscale(1);
  Real copt51  = 1 / copt50;
  Real copt33  = -1 + BARRIERPOWER;
  Real copt75  = strainscale(2);
  Real copt77  = 1 / copt75;
  Real copt95  = strainscale(3);
  Real copt96  = 1 / copt95;
  Real copt106 = strainscale(4);
  Real copt107 = 1 / copt106;
  Real copt115 = strainscale(5);
  Real copt116 = 1 / copt115;
  out(0)       = BARRIERPOWER * bscale * bspeed * copt30 *
           Power(bspeed * copt30 * (-strainclamped(0) + straincpp(0)), copt33);
  out(1) = BARRIERPOWER * bscale * bspeed * copt51 *
           Power(bspeed * copt51 * (-strainclamped(1) + straincpp(1)), copt33);
  out(2) = BARRIERPOWER * bscale * bspeed * copt77 *
           Power(bspeed * copt77 * (-strainclamped(2) + straincpp(2)), copt33);
  out(3) = BARRIERPOWER * bscale * bspeed * copt96 *
           Power(bspeed * copt96 * (-strainclamped(3) + straincpp(3)), copt33);
  out(4) = BARRIERPOWER * bscale * bspeed * copt107 *
           Power(bspeed * copt107 * (-strainclamped(4) + straincpp(4)), copt33);
  out(5) = BARRIERPOWER * bscale * bspeed * copt116 *
           Power(bspeed * copt116 * (-strainclamped(5) + straincpp(5)), copt33);

  return out;
}

gradhess_type SplineMaterial::gradhess_barrier(
    const strain_type &straincpp, const strain_type &strainclamped) {
  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real copt22  = strainscale(0);
  Real copt30  = 1 / copt22;
  Real copt50  = strainscale(1);
  Real copt51  = 1 / copt50;
  Real copt33  = -1 + BARRIERPOWER;
  Real copt75  = strainscale(2);
  Real copt77  = 1 / copt75;
  Real copt95  = strainscale(3);
  Real copt96  = 1 / copt95;
  Real copt106 = strainscale(4);
  Real copt107 = 1 / copt106;
  Real copt115 = strainscale(5);
  Real copt116 = 1 / copt115;
  Real copt120 = Power(bspeed, 2);
  Real copt7   = strainclamped(0);
  Real copt15  = -copt7;
  Real copt18  = straincpp(0);
  Real copt20  = copt15 + copt18;
  Real copt32  = bspeed * copt20 * copt30;
  Real copt123 = Power(copt22, 2);
  Real copt124 = 1 / copt123;
  Real copt37  = strainclamped(1);
  Real copt45  = -copt37;
  Real copt47  = straincpp(1);
  Real copt48  = copt45 + copt47;
  Real copt52  = bspeed * copt48 * copt51;
  Real copt121 = -2 + BARRIERPOWER;
  Real copt127 = Power(copt50, 2);
  Real copt128 = 1 / copt127;
  Real copt63  = strainclamped(2);
  Real copt65  = -copt63;
  Real copt66  = straincpp(2);
  Real copt67  = copt65 + copt66;
  Real copt78  = bspeed * copt67 * copt77;
  Real copt131 = Power(copt75, 2);
  Real copt132 = 1 / copt131;
  Real copt82  = strainclamped(3);
  Real copt90  = -copt82;
  Real copt92  = straincpp(3);
  Real copt93  = copt90 + copt92;
  Real copt97  = bspeed * copt93 * copt96;
  Real copt135 = Power(copt95, 2);
  Real copt136 = 1 / copt135;
  Real copt102 = strainclamped(4);
  Real copt103 = -copt102;
  Real copt104 = straincpp(4);
  Real copt105 = copt103 + copt104;
  Real copt108 = bspeed * copt105 * copt107;
  Real copt141 = Power(copt106, 2);
  Real copt142 = 1 / copt141;
  Real copt111 = strainclamped(5);
  Real copt112 = -copt111;
  Real copt113 = straincpp(5);
  Real copt114 = copt112 + copt113;
  Real copt117 = bspeed * copt114 * copt116;
  Real copt145 = Power(copt115, 2);
  Real copt146 = 1 / copt145;
  out1(0) = BARRIERPOWER * bscale * bspeed * copt30 * Power(copt32, copt33);
  out1(1) = BARRIERPOWER * bscale * bspeed * copt51 * Power(copt52, copt33);
  out1(2) = BARRIERPOWER * bscale * bspeed * copt77 * Power(copt78, copt33);
  out1(3) = BARRIERPOWER * bscale * bspeed * copt96 * Power(copt97, copt33);
  out1(4) = BARRIERPOWER * bscale * bspeed * copt107 * Power(copt108, copt33);
  out1(5) = BARRIERPOWER * bscale * bspeed * copt116 * Power(copt117, copt33);
  out2(0, 0) = BARRIERPOWER * bscale * copt120 * copt124 *
               Power(copt32, copt121) * copt33;
  out2(0, 1) = 0;
  out2(0, 2) = 0;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) = BARRIERPOWER * bscale * copt120 * copt128 * copt33 *
               Power(copt52, copt121);
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) = BARRIERPOWER * bscale * copt120 * copt132 * copt33 *
               Power(copt78, copt121);
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = BARRIERPOWER * bscale * copt120 * copt136 * copt33 *
               Power(copt97, copt121);
  out2(3, 4) = 0;
  out2(3, 5) = 0;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = 0;
  out2(4, 4) = BARRIERPOWER * bscale * Power(copt108, copt121) * copt120 *
               copt142 * copt33;
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) = BARRIERPOWER * bscale * Power(copt117, copt121) * copt120 *
               copt146 * copt33;

  return std::make_pair(hess, grad);
}
