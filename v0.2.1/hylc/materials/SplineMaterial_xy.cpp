#define PC_CURVATURE
#ifndef PC_CURVATURE

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

  // only normalize inplane components, bc bending normalization depends on
  // PC
  Vec6 S;
  for (int i = 0; i < 6; i++)
    S(i) = (strain(i) - this->strainshift[i]) / this->strainscale[i];

  // const
  val += C0;

  // 1D
  for (auto &s : hsplines_1d) {
    double x = S(s.k);
    val += s.eval(x);
  }

  // 2D
  for (auto &s : hsplines_2d) {
    // // double x = X(s.k0);
    // // double y = X(s.k1);
    // // if (!(s.k0 == 0 && s.k0 == 2)) {
    // //   if (s.k0 == 0 || s.k0 == 2)
    // //     if (x < 0)
    // //       continue;

    // //   if (s.k1 == 0 || s.k1 == 2)
    // //     if (y < 0)
    // //       continue;
    // // }

    double x = S(s.k0);
    double y = S(s.k1);
    val += s.eval(x, y);
  }

  return val;
}

Vec6 SplineMaterial::psi_grad(const Vec6 &strain) {
  Vec6 grad(0);
  if (!initialized)
    return grad;

  Vec6 S;
  for (int i = 0; i < 6; i++)
    S(i) = (strain(i) - this->strainshift[i]) / this->strainscale[i];

  // 1D
  for (auto &s : hsplines_1d) {
    double invsc = 1.0 / this->strainscale[s.k];
      double x = S(s.k);
      grad(s.k) += s.dx(x) * invsc;
  }

  // 2D
  for (auto &s : hsplines_2d) {
    continue;  // DEBUG

    // if (!(s.k0 == 0 && s.k0 == 2)) {
    //   if (s.k0 == 0) 
    //     if (S(s.k0) < 0)
    //       continue;

    //   if (s.k1 == 2)
    //     if (S(s.k1) < 0)
    //       continue;
    // }

    double invsc0 = 1.0 / this->strainscale[s.k0];
    double invsc1 = 1.0 / this->strainscale[s.k1];

    double x = S(s.k0);
    double y = S(s.k1);
    grad(s.k0) += s.dx(x, y) * invsc0;
    grad(s.k1) += s.dy(x, y) * invsc1;
  }
  return grad;
}

std::pair<Mat6x6, Vec6> SplineMaterial::psi_drv(const Vec6 &strain) {
  Vec6 grad(0);
  Mat6x6 hess(0);
  if (!initialized)
    return std::make_pair(hess, grad);

  Vec6 S;
  for (int i = 0; i < 6; i++)
    S(i) = (strain(i) - this->strainshift[i]) / this->strainscale[i];

  // 1D
  for (auto &s : hsplines_1d) {
    double invsc = 1.0 / this->strainscale[s.k];
      double x = S(s.k);
      grad(s.k) += s.dx(x) * invsc;
      hess(s.k, s.k) += s.dxdx(x) * invsc * invsc;
  }

  // 2D
  for (auto &s : hsplines_2d) {
    // continue;  // DEBUG

    // if (!(s.k0 == 0 && s.k0 == 2)) {
    //   if (s.k0 == 0) 
    //     if (S(s.k0) < 0)
    //       continue;

    //   if (s.k1 == 2)
    //     if (S(s.k1) < 0)
    //       continue;
    // }

    double invsc0 = 1.0 / this->strainscale[s.k0];
    double invsc1 = 1.0 / this->strainscale[s.k1];

    double x = S(s.k0);
    double y = S(s.k1);
    grad(s.k0) += s.dx(x, y) * invsc0;
    grad(s.k1) += s.dy(x, y) * invsc1;
    hess(s.k0, s.k0) += s.dxdx(x, y) * invsc0 * invsc0;
    double dxdy = s.dxdy(x, y) * invsc0 * invsc1;
    hess(s.k0, s.k1) += dxdy;
    hess(s.k1, s.k0) += dxdy;
    hess(s.k1, s.k1) += s.dydy(x, y) * invsc1 * invsc1;

  }
  return std::make_pair(hess, grad);
}

#endif //!PC_CURVATURE
