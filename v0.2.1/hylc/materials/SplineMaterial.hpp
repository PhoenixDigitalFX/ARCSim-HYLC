#ifndef _SPLINEMATERIAL_HPP_
#define _SPLINEMATERIAL_HPP_

#include "BaseMaterial.hpp"
#include <vector>
#include <memory>

#include <BSplineCurve.h>
#include <BSplineSurface.h>

namespace hylc {
class SplineMaterial : public BaseMaterial {
public:
  SplineMaterial();

  virtual double psi(const Vec6 &strain);
  virtual Vec6 psi_grad(const Vec6 &strain);
  virtual std::pair<Mat6x6, Vec6> psi_drv(const Vec6 &strain);

  struct Spline1D {
    int k;
    std::shared_ptr<fitpackpp::BSplineCurve> spline;
  };

  struct Spline2D {
    int k0, k1;
    std::shared_ptr<fitpackpp::BSplineSurface> spline;
  };

  double C0=0;
  std::vector<Spline1D> splines_1d;
  std::vector<Spline2D> splines_2d;
  
  Vec6 strainshift;
  Vec6 strainscale;
  Vec6 strain_min, strain_max;

  bool initialized;

private:

};

} // namespace hylc
#endif /* end of include guard: _SPLINEMATERIAL_HPP_ */
