#include "FittedMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

#include <iostream>
FittedMaterial::FittedMaterial(int type) {


  bspeed = 1e1; // barrier strength // TODO, also maybe barrier not 6th order TODO
  bscale = 1e2;
  // if barrier is 4th order and bspeed is 10 but strain is normalized same way then it should grow 10^4 times faster ?

  std::cout << "Loading Material " << type << "\n";
  // TODO switch or subclass with other constructor
  density = 4.8530851068e-01;
  ek_min = {0.3, -0.8, 0.3, -70, -40, -70}; // TODO measure bounds vs barrier?
  ek_max = {4.0, 0.8, 4.0, 70, 40, 70};
  // ek_min = {0.1, -2, 0.1, -300, -100, -300}; // TODO <--- confidence
  // ek_max = {7.0, 2, 7.0, 300, 100, 300};
  // ek_max = {0.1, 2, 0.1, 300, 100, 300};// DEBUG
  // ek_min = {0.0, -10, 0.0, -1000, -1000, -1000}; // TODO <--- inf
  // ek_max = {10.0, 10, 10.0, 1000, 1000, 1000};
  ekscale[0] = 3.000000000000000e+00;
  ekscale[1] = 1.600000000000000e+00;
  ekscale[2] = 3.000000000000000e+00;
  ekscale[3] = 1.399864101000000e+02;
  ekscale[4] = 3.546223147500000e+01;
  ekscale[5] = 1.399860782000000e+02;
  double bendmult = 1;//1e2;
  C0    = -2.377419535288772e+00;
  C01   = -6.353491632007823e+00;
  C02   = 2.592702012390881e+02;
  C03   = -4.250956455651524e+02;
  C04   = 2.332829716578757e+02;
  C11   = 8.629258760186682e-01;
  C12   = 3.072934731374260e+01;
  C13   = -2.139756725664153e-01;
  C14   = 1.314536054295262e+01;
  C21   = -9.183306964792202e+00;
  C22   = 1.884982131987582e+02;
  C23   = -2.650059360980511e+02;
  C24   = 1.438318992599692e+02;
  C31   = bendmult*4.826001482726263e-01;
  C32   = bendmult*1.812581445312619e+01;
  C33   = bendmult*-4.203978104893830e-01;
  C34   = bendmult*1.093866293185757e-07;
  C41   = bendmult*-4.743176262124552e-02;
  C42   = bendmult*5.976623978295525e+00;
  C43   = bendmult*1.227988766966033e-01;
  C44   = bendmult*1.000000000000001e-10;
  C51   = bendmult*1.253486154329573e+00;
  C52   = bendmult*1.879349915716852e+01;
  C53   = bendmult*5.514722220845923e-01;
  C54   = bendmult*2.810403860094979e-10;
  C0111 = -2.448591773295877e+00;
  C0112 = -1.035849737305243e+02;
  C0121 = -2.379095122747371e+00;
  C0113 = 1.341023677501257e+00;
  C0122 = 6.895544677906587e+01;
  C0131 = 3.527446260668672e+00;
  C1211 = -8.655932889920050e+00;
  C1221 = -7.594521842691807e+01;
  C1212 = 5.911568344154627e+00;
  C1231 = 1.000000023160279e-10;
  C1222 = 5.290854575246620e+01;
  C1213 = 3.270235067773814e-10;
  C0211 = 2.609697882234321e+02;
  C0212 = -2.183919890524269e+02;
  C0221 = -1.438762325694731e+02;
  C0213 = 5.974622381041176e+01;
  C0222 = 2.108509821644536e+02;
  C0231 = 1.329545128607632e+01;
  C0311 = 1.049208549298532e+01;
  C0312 = -6.889265572956786e+01;
  C0321 = -1.060916969807765e+01;
  C0313 = 2.492782701540272e-01;
  C0322 = 5.097943832059515e+01;
  C0331 = 2.553620077984803e+00;
  C2511 = -9.528691969971483e+00;
  C2512 = -6.663247849733990e+01;
  C2521 = 8.441030251921262e+00;
  C2513 = 1.000000014669581e-10;
  C2522 = 5.097116884220844e+01;
  C2531 = 1.000000000763121e-10;
  C3511 = 4.186805370886005e-01;
  C3512 = -1.995054403583357e-01;
  C3521 = -7.077077308421980e+00;
  C3513 = 1.000000000077098e-10;
  C3522 = 3.107621786093128e+02;
  C3531 = 1.000469198978339e-10;
}

bool FittedMaterial::clamp_strains(Vec6 &ek, std::vector<int> &clamped_coords,
                                   std::vector<double> &dek) {
  bool clamped=false;
  for (int i = 0; i < 6; ++i) {
    if (ek[i] < ek_min[i]) {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_min[i]); // (ek - ekclamped)
      ek[i] = ek_min[i];
      clamped=true;
      // std::cout<<"clamping min "<<i<<" w distance "<<ek[i] - ek_min[i]<<"\n";
    } else if (ek[i] > ek_max[i]) {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_max[i]); // (ek - ekclamped)
      ek[i] = ek_max[i];
      clamped=true;
      // std::cout<<"clamping max "<<i<<" w distance "<<ek[i] - ek_min[i]<<"\n";
    }
  }
  return clamped;
}

#include <iostream>
double FittedMaterial::psi(const Vec6 &ek) {
  Vec6 ekclamped = ek;
  std::vector<int> clamped_coords;
  clamped_coords.reserve(6);
  std::vector<double> dek; // (ei - clamped(ei)) Taylor distance
  dek.reserve(6);
  clamp_strains(ekclamped, clamped_coords, dek);

  double out = 0;
  // std::cout<<ekclamped[0]<<" from "<<ek[0]<<"\n\n";

  // if unclamped dek is 0 and only actual energy remains
  // psi(ek + dek) = psi(ek) + dpsid_i(ek) * dek_i + 0.5 * dpsidd_i(ek) *
  // dek_i^2
  out = psi_taylor_0(
      ekclamped); //  const taylor part, actual energy if not clamped
  // for (size_t i = 0; i < clamped_coords.size(); i++) {
  //   int coord = clamped_coords[i];
  //   double d = dek[i];
  //
  //   //
//   // std::cout<<"clamping "<<i<<" w distance "<<d<<"\n"
  //   std::pair<double, double> psi_D = psi_taylor_12_i(ekclamped, coord);
  //
  //   double tgrad = psi_D.first;
  //   if (d < 0)
  //     tgrad = std::min(-min_taylor_grad, tgrad);
  //   else
  //     tgrad = std::max(min_taylor_grad, tgrad);
  //   double crv = psi_D.second;
  //   crv = std::max(min_taylor_hess, crv);
  //   out += tgrad * d + crv * (0.5 * d * d);
  // }
  out += psi_barrier(ek, ekclamped);

  return out;
}

Vec6 FittedMaterial::psi_grad(const Vec6 &ek) {
  Vec6 ekclamped = ek;
  std::vector<int> clamped_coords;
  clamped_coords.reserve(6);
  std::vector<double> dek; // (ei - clamped(ei)) Taylor distance
  dek.reserve(6);
  clamp_strains(ekclamped, clamped_coords, dek);

  Vec6 out(0);
  out = grad_taylor_0(
      ekclamped); //  const taylor part, actual energy if not clamped
  // for (size_t i = 0; i < clamped_coords.size(); i++) {
  //   int coord = clamped_coords[i];
  //   double d = dek[i];
  //
  //   std::pair<Vec6, Vec6> grad_D = grad_taylor_12_i(ekclamped, coord);
  //
  //   Vec6 tgrad = grad_D.first;
  //   if (d < 0)
  //     for (int i = 0; i < 6; ++i)
  //       tgrad(i) = std::min(-min_taylor_grad, tgrad(i));
  //   else
  //     for (int i = 0; i < 6; ++i)
  //       tgrad(i) = std::max(min_taylor_grad, tgrad(i));
  //   Vec6 crv = grad_D.second;
  //   for (int i = 0; i < 6; ++i)
  //     crv(i) = std::max(min_taylor_hess, crv(i));
  //   out += tgrad * d + crv * (0.5 * d * d);
  // }
  out += grad_barrier(ek, ekclamped);

  return out;
}

std::pair<Mat6x6, Vec6> FittedMaterial::psi_drv(const Vec6 &ek) {

  Vec6 ekclamped = ek;
  std::vector<int> clamped_coords;
  clamped_coords.reserve(6);
  std::vector<double> dek; // (ei - clamped(ei)) Taylor distance
  dek.reserve(6);
  clamp_strains(ekclamped, clamped_coords, dek);

  Mat6x6 hess(0);
  Vec6 grad(0);

  auto gh0 = gradhess_taylor_0(
      ekclamped); //  const taylor part, actual energy if not clamped
  hess += gh0.first;
  grad += gh0.second;
  // for (size_t i = 0; i < clamped_coords.size(); i++) {
  //   int coord = clamped_coords[i];
  //   double d = dek[i];
  //
  //   // {{hD,gD},{hDD,gDD}}
  //   auto gh_D = gradhess_taylor_12_i(ekclamped, coord);
  //
  //   Vec6 tgradgrad = gh_D.first.second;
  //   Mat6x6 tgradhess = gh_D.first.first;
  //   if (d < 0) {
  //     for (int i = 0; i < 6; ++i)
  //       tgradgrad(i) = std::min(-min_taylor_grad, tgradgrad(i));
  //     for (int i = 0; i < 6; ++i)
  //       for (int j = 0; j < 6; ++i)
  //         tgradhess(i, j) = std::min(-min_taylor_grad, tgradhess(i, j));
  //   } else {
  //     for (int i = 0; i < 6; ++i)
  //       tgradgrad(i) = std::max(min_taylor_grad, tgradgrad(i));
  //     for (int i = 0; i < 6; ++i)
  //       for (int j = 0; j < 6; ++i)
  //         tgradhess(i, j) = std::max(min_taylor_grad, tgradhess(i, j));
  //   }
  //
  //   Vec6 crvgrad = gh_D.second.second;
  //   for (int i = 0; i < 6; ++i)
  //     crvgrad(i) = std::max(min_taylor_hess, crvgrad(i));
  //   Mat6x6 crvhess = gh_D.second.first;
  //   for (int i = 0; i < 6; ++i)
  //     for (int j = 0; j < 6; ++i)
  //       crvhess(i, j) = std::max(min_taylor_hess, crvhess(i, j));
  //
  //   hess += tgradhess * d + crvhess * (0.5 * d * d);
  //   grad += tgradgrad * d + crvgrad * (0.5 * d * d);
  // }
  auto ghB = gradhess_barrier(
      ek,
      ekclamped); //  const taylor part, actual energy if not clamped
  hess += ghB.first;
  grad += ghB.second;

  return std::make_pair(hess, grad);
}
