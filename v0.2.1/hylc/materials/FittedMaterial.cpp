#include "FittedMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

#include <iostream>
FittedMaterial::FittedMaterial(int type) {


  s = 1e5; // barrier strength // TODO, also maybe barrier not 6th order TODO

  std::cout << "Loading Material " << type << "\n";
  // TODO switch or subclass with other constructor
  density = 4.8530851068e-01;
  ek_min = {0.3, -0.8, 0.3, -70, -40, -70}; // TODO measure bounds vs barrier?
  ek_max = {4.0, 0.8, 4.0, 70, 40, 70};
  ek_min = {0.1, -2, 0.1, -300, -100, -300}; // TODO <--- confidence
  ek_max = {0.1, 2, 0.1, 300, 100, 300};
  // ek_min = {0.0, -10, 0.0, -1000, -1000, -1000}; // TODO <--- inf
  // ek_max = {10.0, 10, 10.0, 1000, 1000, 1000};
  ekscale[0] = 3.000000000000000e+00;
  ekscale[1] = 1.600000000000000e+00;
  ekscale[2] = 3.000000000000000e+00;
  ekscale[3] = 1.399864101000000e+02;
  ekscale[4] = 3.546223147500000e+01;
  ekscale[5] = 1.399860782000000e+02;
  C0 = -1.426222500212410e-01;
  C01 = -6.451431039609062e+00;
  C02 = 1.325357790095562e+02;
  C03 = -1.172661360597956e+02;
  C04 = 4.841848036614542e+01;
  C11 = 8.629328630779201e-01;
  C12 = 3.025567997317892e+01;
  C13 = -2.139581021504493e-01;
  C14 = 1.543626911175780e+01;
  C21 = -7.971583629075118e+00;
  C22 = 9.405330104166907e+01;
  C23 = -4.937026138467384e+01;
  C24 = 1.988613061812202e+01;
  C31 = 4.825996100640121e-01;
  C32 = 1.319675253165803e+01;
  C33 = -4.203836423948742e-01;
  C34 = 2.414150917335131e+00;
  C41 = -4.743033764878506e-02;
  C42 = 2.469405901769147e+00;
  C43 = 1.227969035067547e-01;
  C44 = 1.004060441757639e-10;
  C51 = 1.253477268891114e+00;
  C52 = 1.020543882286709e+01;
  C53 = 5.514749282422567e-01;
  C54 = 5.159079000954860e+00;
  C0111 = -2.448662000077315e+00;
  C0112 = -1.033129390075163e+02;
  C0121 = -2.379014523949212e+00;
  C0113 = 1.341013409687814e+00;
  C0122 = 6.692090709052944e+01;
  C0131 = 3.527426975535120e+00;
  C1211 = -8.655992590829319e+00;
  C1221 = -7.607626048827098e+01;
  C1212 = 5.911613288903617e+00;
  C1231 = 1.000000262813812e-10;
  C1222 = 4.928446997946500e+01;
  C1213 = 1.000019635401658e-10;
  C0211 = 1.140914558437495e+02;
  C0212 = 1.685100813465264e+01;
  C0221 = 4.063907930831111e+01;
  C0213 = 1.000805688736867e-10;
  C0222 = 1.524843112638780e-10;
  C0231 = 1.000000010929201e-10;
  C0311 = 1.049202433328033e+01;
  C0312 = -5.519173763402622e+01;
  C0321 = -1.060899811260474e+01;
  C0313 = 2.492529307844227e-01;
  C0322 = 3.958724193594173e+01;
  C0331 = 2.553515202020149e+00;
  C2511 = -9.528630021715076e+00;
  C2512 = -4.554089437944626e+01;
  C2521 = 8.440968537305412e+00;
  C2513 = 1.300904972154314e-10;
  C2522 = 3.144551170927630e+01;
  C2531 = 1.182396571688130e-10;
  C3511 = 1.049825812874267e+00;
  C3512 = -1.995315235047130e-01;
  C3521 = -7.076971262814091e+00;
  C3513 = 1.000001601509200e-10;
  C3522 = 6.738986289577487e+01;
  C3531 = 2.499568968977100e-01;
}

void FittedMaterial::clamp_strains(Vec6 &ek, std::vector<int> &clamped_coords,
                                   std::vector<double> &dek) {
  // return;
  for (int i = 0; i < 6; ++i) {
    if (ek[i] < ek_min[i]) {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_min[i]); // (ek - ekclamped)
      ek[i] = ek_min[i];
      // std::cout<<"clamping min "<<i<<" w distance "<<ek[i] - ek_min[i]<<"\n";
    } else if (ek[i] > ek_max[i]) {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_max[i]); // (ek - ekclamped)
      ek[i] = ek_max[i];
      // std::cout<<"clamping max "<<i<<" w distance "<<ek[i] - ek_min[i]<<"\n";
    }
  }
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
