#include "SplineMaterial.hpp"
#include <memory>
#include <vector>

using namespace hylc;
using namespace fitpackpp;

SplineMaterial::SplineMaterial() {}

double SplineMaterial::psi(const Vec6 &strain) {
  if (!initialized)
    return 0;

  // normalize strain
  Vec6 X;
  for(int i = 0; i < 6; i++)
    X(i) = (strain(i) - this->strainshift(i))/this->strainscale(i);
  


  // const
  double val = C0;

  // 1D
  for(auto & spline1d : splines_1d) {
    spline1d.spline->eval(X(spline1d.k));
  }
  
  // ...

  // 2D
  // ...


  return val;  // + 1d + 2d
}

Vec6 SplineMaterial::psi_grad(const Vec6 &strain) {
  Vec6 grad(0);
  if (!initialized)
    return grad;

  // normalize strain
  Vec6 X;
  for(int i = 0; i < 6; i++)
    X(i) = (strain(i) - this->strainshift(i))/this->strainscale(i);

  // 1D
  for(auto & spline1d : splines_1d) {
    grad(spline1d.k) += spline1d.spline->der(X(spline1d.k), 1);
  }

  // 2D
  // ...

  return grad;
}

std::pair<Mat6x6, Vec6> SplineMaterial::psi_drv(const Vec6 &strain) {
  Vec6 grad(0);
  Mat6x6 hess(0);
  if (!initialized)
    return std::make_pair(hess, grad);

  // normalize strain
  Vec6 X;
  for(int i = 0; i < 6; i++)
    if (i > 2)
      X(i) = std::min(std::max(strain(i),this->strain_min(i)),this->strain_max(i));
      
  for(int i = 0; i < 6; i++)
    X(i) = (X(i) - this->strainshift(i))/this->strainscale(i);

  // for(int i = 0; i < 6; i++)
  //   printf("   %.2f", strain(i));


  // 1D
  for(auto & spline1d : splines_1d) {
    grad(spline1d.k) += spline1d.spline->der(X(spline1d.k), 1);
    hess(spline1d.k, spline1d.k) += spline1d.spline->der(X(spline1d.k), 2);
    // grad(spline1d.k) += 2 * X(spline1d.k);
    // hess(spline1d.k, spline1d.k) += 2;
  }

  // 2D
  // ...

  grad(0) = 1.0;
  return std::make_pair(hess, grad);
}
