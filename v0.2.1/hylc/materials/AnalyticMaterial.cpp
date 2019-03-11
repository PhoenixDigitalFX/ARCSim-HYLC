#include "AnalyticMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

double AnalyticMaterial::psi(const Vec6 &ekdof) {
  Real c1 = ekdof(0);
  Real c7 = ekdof(2);
  Real c16 = ekdof(3);
  Real c17 = ekdof(5);
  return b0 * Power(c16 + c17, 2) + a0 * Power(-2 + c1 + c7, 2) +
         a1 * (Power(-1 + c1, 2) + Power(-1 + c7, 2) + 2 * Power(ekdof(1), 2)) +
         b1 * Power(-(c16 * c17) + Power(ekdof(4), 2), 2);
}

Vec6 AnalyticMaterial::psi_grad(const Vec6 &ekdof) {
  Vec6 out(0);

  Real c1 = ekdof(0);
  Real c4 = ekdof(2);
  Real c5 = -2 + c1 + c4;
  Real c6 = 2 * a0 * c5;
  Real c15 = ekdof(5);
  Real c14 = ekdof(3);
  Real c18 = ekdof(4);
  Real c19 = Power(c18, 2);
  Real c16 = c14 + c15;
  Real c17 = 2 * b0 * c16;
  Real c21 = -c19;
  Real c22 = c14 * c15;
  Real c23 = c21 + c22;
  out(0) = 2 * a1 * (-1 + c1) + c6;
  out(1) = 4 * a1 * ekdof(1);
  out(2) = 2 * a1 * (-1 + c4) + c6;
  out(3) = c17 + 2 * b1 * c15 * c23;
  out(4) = 4 * b1 * c18 * (-(c14 * c15) + c19);
  out(5) = c17 + 2 * b1 * c14 * c23;

  return out;
}

std::pair<Mat6x6, Vec6> AnalyticMaterial::psi_drv(const Vec6 &ekdof) {

  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real c1 = ekdof(0);
  Real c4 = ekdof(2);
  Real c5 = -2 + c1 + c4;
  Real c6 = 2 * a0 * c5;
  Real c15 = ekdof(5);
  Real c14 = ekdof(3);
  Real c18 = ekdof(4);
  Real c19 = Power(c18, 2);
  Real c16 = c14 + c15;
  Real c17 = 2 * b0 * c16;
  Real c21 = -c19;
  Real c22 = c14 * c15;
  Real c23 = c21 + c22;
  Real c33 = 2 * a0;
  Real c31 = a0 + a1;
  Real c32 = 2 * c31;
  Real c39 = -4 * b1 * c15 * c18;
  Real c26 = -(c14 * c15);
  Real c40 = -(b1 * c19);
  Real c42 = 2 * b1 * c14 * c15;
  Real c43 = b0 + c40 + c42;
  Real c44 = 2 * c43;
  Real c48 = -4 * b1 * c14 * c18;
  out1(0) = 2 * a1 * (-1 + c1) + c6;
  out1(1) = 4 * a1 * ekdof(1);
  out1(2) = 2 * a1 * (-1 + c4) + c6;
  out1(3) = c17 + 2 * b1 * c15 * c23;
  out1(4) = 4 * b1 * c18 * (c19 + c26);
  out1(5) = c17 + 2 * b1 * c14 * c23;
  out2(0, 0) = c32;
  out2(0, 1) = 0;
  out2(0, 2) = c33;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) = 4 * a1;
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = c33;
  out2(2, 1) = 0;
  out2(2, 2) = c32;
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = 2 * (b0 + b1 * Power(c15, 2));
  out2(3, 4) = c39;
  out2(3, 5) = c44;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = c39;
  out2(4, 4) = 4 * b1 * (3 * c19 + c26);
  out2(4, 5) = c48;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = c44;
  out2(5, 4) = c48;
  out2(5, 5) = 2 * (b0 + b1 * Power(c14, 2));

  return std::make_pair(hess, grad);
}
