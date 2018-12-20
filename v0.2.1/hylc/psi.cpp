#include "hylc_mm.hpp"

using namespace hylc;
using namespace hylc::mathematica;

// (ek)(\d+) -> $1($2)

Real hylc::mathematica::psi(const Vec6 &_ek, double a0, double a1, double b0,
                            double b1) {
  // define input
  auto ek = [&](int i) -> const Real & { return _ek[i - 1]; };
  return a1 *
             (Power(-1 + ek(1), 2) + 2 * ek(2) * ek(2) + Power(-1 + ek(3), 2)) +
         a0 * Power(-2 + ek(1) + ek(3), 2) + b0 * Power(ek(4) + ek(6), 2) +
         b1 * Power(ek(5) * ek(5) - ek(4) * ek(6), 2);
}

std::pair<Mat6x6, Vec6> hylc::mathematica::psi_drv(const Vec6 &_ek, double a0,
                                                   double a1, double b0,
                                                   double b1) {
  // define input
  auto ek = [&](int i) -> const Real & { return _ek[i - 1]; };

  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real c3 = -2 + ek(1) + ek(3);
  Real c4 = 2 * a0 * c3;
  Real c14 = ek(5) * ek(5);
  Real c11 = ek(4) + ek(6);
  Real c12 = 2 * b0 * c11;
  Real c15 = -c14;
  Real c16 = ek(4) * ek(6);
  Real c18 = c15 + c16;
  Real c32 = 2 * a0;
  Real c29 = a0 + a1;
  Real c30 = 2 * c29;
  Real c42 = -4 * b1 * ek(5) * ek(6);
  Real c22 = -(ek(4) * ek(6));
  Real c46 = -(b1 * c14);
  Real c47 = 2 * b1 * ek(4) * ek(6);
  Real c48 = b0 + c46 + c47;
  Real c49 = 2 * c48;
  Real c54 = -4 * b1 * ek(4) * ek(5);
  out1(0) = 2 * a1 * (-1 + ek(1)) + c4;
  out1(1) = 4 * a1 * ek(2);
  out1(2) = 2 * a1 * (-1 + ek(3)) + c4;
  out1(3) = c12 + 2 * b1 * ek(6) * c18;
  out1(4) = 4 * b1 * ek(5) * (c14 + c22);
  out1(5) = c12 + 2 * b1 * ek(4) * c18;
  out2(0, 0) = c30;
  out2(0, 1) = 0;
  out2(0, 2) = c32;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) = 4 * a1;
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = c32;
  out2(2, 1) = 0;
  out2(2, 2) = c30;
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = 2 * (b0 + b1 * ek(6) * ek(6));
  out2(3, 4) = c42;
  out2(3, 5) = c49;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = c42;
  out2(4, 4) = 4 * b1 * (3 * c14 + c22);
  out2(4, 5) = c54;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = c49;
  out2(5, 4) = c54;
  out2(5, 5) = 2 * (b0 + b1 * ek(4) * ek(4));

  return std::make_pair(hess, grad);
}
