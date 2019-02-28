#ifndef _ANALYTICMATERIAL_HPP_
#define _ANALYTICMATERIAL_HPP_

#include "BaseMaterial.hpp"
namespace hylc {
class AnalyticMaterial : public BaseMaterial {
public:
  AnalyticMaterial(double a0, double a1, double b0, double b1) {
    this->a0 = a0;
    this->a1 = a1;
    this->b0 = b0;
    this->b1 = b1;
    this->density = 0.187; // copied from gray-interlock
  }
  virtual double psi(const Vec6 &ek);
  virtual std::pair<Mat6x6, Vec6> psi_drv(const Vec6 &ek);

private:
  double a0, a1, b0, b1;
};

} // namespace hylc
#endif /* end of include guard: _ANALYTICMATERIAL_HPP_ */
