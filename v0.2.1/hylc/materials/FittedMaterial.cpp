#include "FittedMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

#include <iostream>
FittedMaterial::FittedMaterial(int type) {

  bspeed =
      1e1; // barrier strength // TODO, also maybe barrier not 6th order TODO
  bscale = 1e2;
  // if barrier is 4th order and bspeed is 10 but strain is normalized same way
  // then it should grow 10^4 times faster ?

  std::cout << "Loading Material " << type << "\n";
  // TODO switch or subclass with other constructor
  density = 4.8530851068e-01;
  ek_min = {0.0, -0.8, 0.0, -70, -40, -70}; // TODO measure bounds vs barrier?
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
  double bendmult = 1; // 1e2;
  C0 = -3.173365863561212e+00;
  C01 = 1.274173003488781e+01;
  C02 = 4.747231875296465e+01;
  C03 = -4.480497102275712e-01;
  C04 = 1.000000062449521e-10;
  C11 = -8.254626670431638e-02;
  C12 = 1.694450101740435e+01;
  C13 = 8.856231617824880e-01;
  C14 = 1.249142869700673e+01;
  C21 = 1.293384050109903e+01;
  C22 = 2.304301037900553e+01;
  C23 = 2.279546080079567e+01;
  C24 = 1.000000002700019e-10;
  C31 = 4.271820882958070e-01;
  C32 = 7.105730259356392e+00;
  C33 = -4.512893307033707e-01;
  C34 = 3.346153921544564e-03;
  C41 = -4.596599884973425e-02;
  C42 = 9.500414130012255e-01;
  C43 = 1.210041988288516e-01;
  C44 = 1.000000000000354e-10;
  C51 = 2.692952968668278e-01;
  C52 = 3.674778093960921e+00;
  C53 = 7.009053310225889e-01;
  C54 = 9.843922790388412e-01;
  C0111 = 1.228169641433861e+00;
  C0112 = -4.178450046639179e+01;
  C0121 = -8.272436924352622e+00;
  C0113 = 1.722200262355749e-01;
  C0122 = 1.968737309207855e+01;
  C0131 = 6.832757413035157e+00;
  C1211 = -6.091033484877289e+00;
  C1221 = -2.184350879510441e+01;
  C1212 = 3.540093247569986e+00;
  C1231 = 1.000000424511114e-10;
  C1222 = 1.031121344118416e+01;
  C1213 = 1.077358510766908e-10;
  C0211 = 1.295350875922079e+01;
  C0212 = 8.264307218082685e+01;
  C0221 = 1.171109934061652e+02;
  C0213 = 1.000031758144436e-10;
  C0222 = 1.000000216137947e-10;
  C0231 = 1.000003346836002e-10;
  C0311 = 1.198363764520588e+01;
  C0312 = -2.174293040442834e+01;
  C0321 = -1.487177753497701e+01;
  C0313 = 1.814448778738220e-01;
  C0322 = 1.300367900051241e+01;
  C0331 = 5.491209974476847e+00;
  C2511 = -5.298133748971060e+00;
  C2512 = -1.424386818243591e+01;
  C2521 = 5.424247859338142e+00;
  C2513 = 1.656680848353157e-10;
  C2522 = 1.104873452274377e+01;
  C2531 = 1.055998868359542e-10;
  C3511 = 3.179424187234492e-01;
  C3512 = 1.829232842834747e-01;
  C3521 = 6.183366848611576e+00;
  C3513 = 1.000000021229392e-10;
  C3522 = 7.392664856532826e+01;
  C3531 = 6.904746092521139e+00;
  Ccompr0 = 4.078976900696096e+00;
  Ccompr1 = 1.518191854741584e-02;
}

bool FittedMaterial::clamp_strains(Vec6 &ek, std::vector<int> &clamped_coords,
                                   std::vector<double> &dek) {
  bool clamped = false;
  for (int i = 0; i < 6; ++i) {
    if (ek[i] < ek_min[i]) {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_min[i]); // (ek - ekclamped)
      ek[i] = ek_min[i];
      clamped = true;
      // std::cout<<"clamping min "<<i<<" w distance "<<ek[i] - ek_min[i]<<"\n";
    } else if (ek[i] > ek_max[i]) {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_max[i]); // (ek - ekclamped)
      ek[i] = ek_max[i];
      clamped = true;
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
  out += psi_compr(ek);

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
  out += grad_compr(ek);

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
      ekclamped);
  hess += ghB.first;
  grad += ghB.second;


  auto ghC = gradhess_compr(
      ek);
  hess += ghC.first;
  grad += ghC.second;

  return std::make_pair(hess, grad);
}
