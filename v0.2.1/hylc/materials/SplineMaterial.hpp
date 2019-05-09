#ifndef _SPLINEMATERIAL_HPP_
#define _SPLINEMATERIAL_HPP_

#include <memory>
#include <vector>
#include "BaseMaterial.hpp"

#include <BSplineCurve.h>
#include <BSplineSurface.h>
#include "HermiteSpline1D.hpp"
#include "HermiteSpline2D.hpp"
#include "Poly1D.hpp"
#include "Poly2D.hpp"

namespace hylc {
class SplineMaterial : public BaseMaterial {
 public:
  SplineMaterial();

  virtual double psi(const Vec6 &strain);
  virtual Vec6 psi_grad(const Vec6 &strain);
  virtual std::pair<Mat6x6, Vec6> psi_drv(const Vec6 &strain);

  struct Spline1D {
    int k;
    double xmin, xmax;
    bool clampx;
    inline void clamp(double &x) {
      if (clampx)
        x = std::min(std::max(x, xmin), xmax);
    }
    std::shared_ptr<fitpackpp::BSplineCurve> spline;
  };

  struct Spline2D {
    int k0, k1;
    std::shared_ptr<fitpackpp::BSplineSurface> spline;
  };

  struct HSpline2Das1D {
    int k0, k1;
    // HermiteSpline1D fun;
    Poly1D fun;
  };

  double C0 = 0;
  std::vector<HermiteSpline1D> hsplines_1d;
  std::vector<HermiteSpline2D> hsplines_2d;
  std::vector<HSpline2Das1D> hsplines_2d1d;
  std::vector<Spline1D> splines_1d;
  std::vector<Spline2D> splines_2d;
  std::vector<Poly1D> polys_1d;
  std::vector<Poly2D> polys_2d;

  Vec6 strainshift;
  Vec6 strainscale;
  // Vec6 strain_min, strain_max;

  bool initialized;

  // double TST(Vec6 &x) {
  //   double a = 0;

  //   for (auto &spline1d : splines_1d)
  //     a += spline1d.spline->der(x(spline1d.k), 1);
  //   return a;
  // }

 private:
  double bspeed = 1e0, bscale = 1e2;
  bool use_barrier = true;
  double psi_barrier(const Vec6 &strain, const Vec6 &strainclamped);
  Vec6 grad_barrier(const Vec6 &strain, const Vec6 &strainclamped);
  std::pair<Mat6x6, Vec6> gradhess_barrier(const Vec6 &strain,
                                           const Vec6 &strainclamped);
};

}  // namespace hylc
#endif /* end of include guard: _SPLINEMATERIAL_HPP_ */
