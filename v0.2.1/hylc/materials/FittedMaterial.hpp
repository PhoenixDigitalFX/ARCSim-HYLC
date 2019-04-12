#ifndef _FITTEDMATERIAL_HPP_
#define _FITTEDMATERIAL_HPP_

#include "BaseMaterial.hpp"
#include <tuple>
#include <vector>

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
  double C0 = 0, C01 = 0, C02 = 0, C03 = 0, C04 = 0, C11 = 0, C12 = 0, C13 = 0,
         C14 = 0, C21 = 0, C22 = 0, C23 = 0, C24 = 0, C31 = 0, C32 = 0, C33 = 0,
         C34 = 0, C41 = 0, C42 = 0, C43 = 0, C44 = 0, C51 = 0, C52 = 0, C53 = 0,
         C54 = 0, C0111 = 0, C0112 = 0, C0121 = 0, C0113 = 0, C0122 = 0,
         C0131 = 0, C0211 = 0, C0212 = 0, C0221 = 0, C0213 = 0, C0222 = 0,
         C0231 = 0, C0311 = 0, C0312 = 0, C0321 = 0, C0313 = 0, C0322 = 0,
         C0331 = 0, C0411 = 0, C0412 = 0, C0421 = 0, C0413 = 0, C0422 = 0,
         C0431 = 0, C0511 = 0, C0512 = 0, C0521 = 0, C0513 = 0, C0522 = 0,
         C0531 = 0, C1211 = 0, C1212 = 0, C1221 = 0, C1213 = 0, C1222 = 0,
         C1231 = 0, C1311 = 0, C1312 = 0, C1321 = 0, C1313 = 0, C1322 = 0,
         C1331 = 0, C1411 = 0, C1412 = 0, C1421 = 0, C1413 = 0, C1422 = 0,
         C1431 = 0, C1511 = 0, C1512 = 0, C1521 = 0, C1513 = 0, C1522 = 0,
         C1531 = 0, C2311 = 0, C2312 = 0, C2321 = 0, C2313 = 0, C2322 = 0,
         C2331 = 0, C2411 = 0, C2412 = 0, C2421 = 0, C2413 = 0, C2422 = 0,
         C2431 = 0, C2511 = 0, C2512 = 0, C2521 = 0, C2513 = 0, C2522 = 0,
         C2531 = 0, C3411 = 0, C3412 = 0, C3421 = 0, C3413 = 0, C3422 = 0,
         C3431 = 0, C3511 = 0, C3512 = 0, C3521 = 0, C3513 = 0, C3522 = 0,
         C3531 = 0, C4511 = 0, C4512 = 0, C4521 = 0, C4513 = 0, C4522 = 0,
         C4531 = 0;
  std::vector<double> ek_min, ek_max; // TODO from pydata

  bool clamp_strains(Vec6 &ek, std::vector<int> &clamped_coords,
                     std::vector<double> &dek);
  double min_taylor_grad = 0.0;
  double min_taylor_hess = 1e-5;

  // barrier
  double bspeed, bscale;
  double psi_barrier(const Vec6 &ek, const Vec6 &ekclamped);
  Vec6 grad_barrier(const Vec6 &ek, const Vec6 &ekclamped);
  std::pair<Mat6x6, Vec6> gradhess_barrier(const Vec6 &ek,
                                           const Vec6 &ekclamped);

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
