#ifndef _FITTEDMATERIAL_HPP_
#define _FITTEDMATERIAL_HPP_

#include "BaseMaterial.hpp"
#include <vector>
#include <tuple>

namespace hylc {
class FittedMaterial : public BaseMaterial {
public:
  FittedMaterial(int type);

  virtual double psi(const Vec6 &ek);
  virtual std::pair<Mat6x6, Vec6> psi_drv(const Vec6 &ek);
  virtual Vec6 psi_grad(const Vec6 &ek);

private:
  // base::density
  // coefficients
  Vec6 ekscale;
  double C0, C01, C02, C03, C04, C11, C12, C13, C14, C21, C22, C23, C24, C31,
      C32, C33, C34, C41, C42, C43, C44, C51, C52, C53, C54, C0111, C0112,
      C0121, C0113, C0122, C0131, C0211, C0212, C0221, C0213, C0222, C0231,
      C0311, C0312, C0321, C0313, C0322, C0331, C0411, C0412, C0421, C0413,
      C0422, C0431, C0511, C0512, C0521, C0513, C0522, C0531, C2311, C2312,
      C2321, C2313, C2322, C2331, C2411, C2412, C2421, C2413, C2422, C2431,
      C2511, C2512, C2521, C2513, C2522, C2531, C3411, C3412, C3421, C3413,
      C3422, C3431, C3511, C3512, C3521, C3513, C3522, C3531;
  std::vector<double> ek_min, ek_max; // TODO from pydata

  void clamp_strains(Vec6 &ek, std::vector<int> &clamped_coords,
                     std::vector<double> &dek);
  double min_taylor_grad = 0.0;
  double min_taylor_hess = 1e-5;


  // barrier
  double s;
  double psi_barrier(const Vec6 &ek,const Vec6 &ekclamped);
  Vec6 grad_barrier(const Vec6 &ek,const Vec6 &ekclamped);
  std::pair<Mat6x6, Vec6> gradhess_barrier(const Vec6 &ek,const Vec6 &ekclamped);

  // 0th derivative, i.e. actual values
  double psi_taylor_0(const Vec6 &ek);
  Vec6 grad_taylor_0(const Vec6 &ek);
  std::pair<Mat6x6, Vec6> gradhess_taylor_0(const Vec6 &ek);
  // 1st and 2nd derivative per strain coord for taylor
  std::pair<double, double> psi_taylor_12_i(const Vec6 &ek, int i);
  std::pair<Vec6, Vec6> grad_taylor_12_i(const Vec6 &ek, int i);
  std::pair<std::pair<Mat6x6, Vec6>, std::pair<Mat6x6, Vec6>>
  gradhess_taylor_12_i(const Vec6 &ek, int i);


  std::pair<double, double> psi_taylor_12_0(const Vec6 &ek);
  std::pair<double, double> psi_taylor_12_1(const Vec6 &ek);
  std::pair<double, double> psi_taylor_12_2(const Vec6 &ek);
  std::pair<double, double> psi_taylor_12_3(const Vec6 &ek);
  std::pair<double, double> psi_taylor_12_4(const Vec6 &ek);
  std::pair<double, double> psi_taylor_12_5(const Vec6 &ek);
  std::pair<Vec6, Vec6> grad_taylor_12_0(const Vec6 &ek);
  std::pair<Vec6, Vec6> grad_taylor_12_1(const Vec6 &ek);
  std::pair<Vec6, Vec6> grad_taylor_12_2(const Vec6 &ek);
  std::pair<Vec6, Vec6> grad_taylor_12_3(const Vec6 &ek);
  std::pair<Vec6, Vec6> grad_taylor_12_4(const Vec6 &ek);
  std::pair<Vec6, Vec6> grad_taylor_12_5(const Vec6 &ek);
  std::pair<std::pair<Mat6x6, Vec6>, std::pair<Mat6x6, Vec6>>
  gradhess_taylor_12_0(const Vec6 &ek);
  std::pair<std::pair<Mat6x6, Vec6>, std::pair<Mat6x6, Vec6>>
  gradhess_taylor_12_1(const Vec6 &ek);
  std::pair<std::pair<Mat6x6, Vec6>, std::pair<Mat6x6, Vec6>>
  gradhess_taylor_12_2(const Vec6 &ek);
  std::pair<std::pair<Mat6x6, Vec6>, std::pair<Mat6x6, Vec6>>
  gradhess_taylor_12_3(const Vec6 &ek);
  std::pair<std::pair<Mat6x6, Vec6>, std::pair<Mat6x6, Vec6>>
  gradhess_taylor_12_4(const Vec6 &ek);
  std::pair<std::pair<Mat6x6, Vec6>, std::pair<Mat6x6, Vec6>>
  gradhess_taylor_12_5(const Vec6 &ek);
};

} // namespace hylc
#endif /* end of include guard: _FITTEDMATERIAL_HPP_ */
