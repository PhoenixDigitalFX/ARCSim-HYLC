#include "FittedMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

#include <iostream>
FittedMaterial::FittedMaterial(int type)
{
  std::cout << "Loading Material " << type << "\n";

  bspeed =
      1e1; // barrier strength // TODO, also maybe barrier not 6th order TODO
  bscale = 1e2;
  // if barrier is 4th order and bspeed is 10 but strain is normalized same way
  // then it should grow 10^4 times faster ?
  if (type == 0)
  {
    ek_min = {0.0, -0.8, 0.0, -70, -40, -70}; // TODO measure bounds vs barrier?
    ek_max = {4.0, 0.8, 4.0, 70, 40, 70};
    density = 4.1562049890e-01;
    ekscale[0] = 3.000000000000000e+00;
    ekscale[1] = 1.600000000000000e+00;
    ekscale[2] = 3.000000000000000e+00;
    ekscale[3] = 1.399920440000000e+02;
    ekscale[4] = 3.695105080000000e+01;
    ekscale[5] = 1.399918732000000e+02;
    C0 = -1.383312731430795e+01;
    C01 = 4.762220787822391e+01;
    C02 = 6.122951639018105e+01;
    C03 = -9.461277161058091e+01;
    C04 = 6.115528143343278e+01;
    C11 = 1.448438300117031e+00;
    C12 = 4.717967449161662e+01;
    C13 = -3.702742987568058e+00;
    C14 = 2.811012476573370e+01;
    C21 = 6.101152400769500e+01;
    C22 = -4.968584329432040e+01;
    C23 = 6.259984771585656e+01;
    C24 = 1.000033181737568e-10;
    C31 = 9.773343437313715e-01;
    C32 = 1.465787392525264e+01;
    C33 = 3.993696939285856e-01;
    C34 = 7.470584424894682e+00;
    C41 = -6.520808987623057e-02;
    C42 = 4.651187786120465e+00;
    C43 = 5.307754645070940e-02;
    C44 = 1.000000000000225e-10;
    C51 = 1.987251408773745e+00;
    C52 = 1.808780140616588e+01;
    C53 = 7.852803401673418e-01;
    C54 = 8.148243513289176e+00;
    C0111 = 9.114906991559929e-01;
    C0112 = -2.041338724630685e+02;
    C0121 = -9.638701359558235e+00;
    C0113 = 5.940823391337779e+00;
    C0122 = 1.423373880483361e+02;
    C0131 = 6.158802442149463e+00;
    C1211 = -9.441862544700703e+00;
    C1221 = -1.644296277128158e+02;
    C1212 = 7.819464785325664e+00;
    C1231 = 1.000004763154176e-10;
    C1222 = 1.130766793318817e+02;
    C1213 = 6.575870791622929e-08;
    C0211 = -5.797676770607567e+01;
    C0212 = 1.004644230268284e+02;
    C0221 = 1.597103946106755e+02;
    C0213 = 1.000011114437434e-10;
    C0222 = 1.000061228019741e-10;
    C0231 = 1.000000077898358e-10;
    C0311 = 3.975386385315942e+00;
    C0312 = -8.913289975691512e+01;
    C0321 = -2.717351741841032e+00;
    C0313 = 1.001742632894281e-10;
    C0322 = 6.881917013359228e+01;
    C0331 = 1.000004506030324e-10;
    C2511 = -1.314465484142953e+01;
    C2512 = -9.030255356049773e+01;
    C2521 = 1.242018284461180e+01;
    C2513 = 3.922546834972385e-10;
    C2522 = 6.633114137452904e+01;
    C2531 = 1.000000101429552e-10;
    C3511 = -3.885759169184497e-01;
    C3512 = -5.916779466773975e+00;
    C3521 = -2.020823664793698e+01;
    C3513 = 1.000002274119113e-10;
    C3522 = 1.072198415859202e+02;
    C3531 = 1.000000437504055e-10;
    Ccompr0 = 1.221488122826855e+01;
    Ccompr1 = 1.000000000000004e-10;
  }

  else if (type == 1)
  {
    ek_min = {0.0, -0.8, 0.0, -100, -60, -100}; // TODO measure bounds vs barrier?
    ek_max = {4.0, 0.8, 4.0, 100, 60, 100};
    density = 6.0205506498e-01;
    ekscale[0] = 3.000000000000000e+00;
    ekscale[1] = 1.600000000000000e+00;
    ekscale[2] = 3.000000000000000e+00;
    ekscale[3] = 1.399937629000000e+02;
    ekscale[4] = 5.298510940500000e+01;
    ekscale[5] = 1.399941661000000e+02;
    C0 = -2.237123353487991e+02;
    C01 = 1.856552356350204e+03;
    C02 = -1.149045410101043e+03;
    C03 = 4.572064757220791e+02;
    C04 = 1.000000000000000e-10;
    C11 = -4.313486213787195e+02;
    C12 = 1.194367599760837e+03;
    C13 = -1.856168649235895e+02;
    C14 = 3.031395362344056e+02;
    C21 = 1.494755753128106e+03;
    C22 = -1.330209247230039e+02;
    C23 = 1.791063248489721e+01;
    C24 = 2.181783107250892e-08;
    C31 = 1.086084852910490e+02;
    C32 = 7.832489473787348e+01;
    C33 = -1.580750041936218e+01;
    C34 = 1.000000000001496e-10;
    C41 = 2.877892776033187e+01;
    C42 = 6.628269041901426e+01;
    C43 = 1.015537865995712e+01;
    C44 = 1.000000000000000e-10;
    C51 = 7.705066661334546e+01;
    C52 = 1.036728755915088e+02;
    C53 = -1.706292860175709e+01;
    C54 = 1.000000000000000e-10;
    C0111 = 2.460924408691637e+02;
    C0112 = -2.796668416107188e+03;
    C0121 = -2.948667062958633e+02;
    C0113 = 1.560531694311601e+02;
    C0122 = 1.652379097660342e+03;
    C0131 = 1.636240201580273e+02;
    C1211 = 6.094482576135239e+01;
    C1221 = -3.090647098248120e+03;
    C1212 = 5.596589465665287e+00;
    C1231 = 1.511329920555722e+02;
    C1222 = 1.798014220060494e+03;
    C1213 = 1.000000000036182e-10;
    C0211 = -1.937384027778640e+02;
    C0212 = 2.830967002095807e+02;
    C0221 = -1.134129877834272e+03;
    C0213 = 1.002436401457279e-10;
    C0222 = 2.974468360059241e-10;
    C0231 = 1.180680643446769e+03;
    C0311 = -3.134614968436791e+02;
    C0312 = -2.711619096689856e+02;
    C0321 = 2.250697041662561e+02;
    C0313 = 1.830245801971105e+01;
    C0322 = 2.158439748145468e+02;
    C0331 = 1.000000005192369e-10;
    C2511 = -1.157351630003389e+02;
    C2512 = -2.774563625429581e+02;
    C2521 = 5.290984449634400e+01;
    C2513 = 2.374381124560871e+01;
    C2522 = 1.640019919356215e+02;
    C2531 = 1.000000271822735e-10;
    C3511 = 4.588771291977308e+00;
    C3512 = -1.814107521758856e+02;
    C3521 = -8.537296250185464e+01;
    C3513 = 1.000000441651963e-10;
    C3522 = 1.718474286135270e+03;
    C3531 = 1.084361830906294e+02;
    Ccompr0 = 1.270618777149060e+03;
    Ccompr1 = 4.510877350994876e-01;
  }
  else if (type == 2)
  {
    ek_min = {0.0, -0.8, 0.0, -150, -100, -150}; // TODO measure bounds vs barrier?
    ek_max = {4.0, 0.8, 4.0, 150, 100, 150};
    density = 3.1553164892e-01;
    ekscale[0] = 3.000000000000000e+00;
    ekscale[1] = 1.600000000000000e+00;
    ekscale[2] = 3.000000000000000e+00;
    ekscale[3] = 1.499743416000000e+02;
    ekscale[4] = 7.744316544000000e+01;
    ekscale[5] = 1.499745997000000e+02;
    C0 = -2.161800502701954e+02;
    C01 = 7.349978508294802e+02;
    C02 = -4.764606392301475e+02;
    C03 = 4.222804946897930e+02;
    C04 = 4.935943641668013e-10;
    C11 = 2.749487249681770e+00;
    C12 = 5.029959981949225e+02;
    C13 = 3.584301930890814e-01;
    C14 = 2.056886507688828e+02;
    C21 = 6.061936092370763e+02;
    C22 = -1.114129036372289e+02;
    C23 = 1.506730115729834e+02;
    C24 = 6.173224465767273e-07;
    C31 = 1.389364655800102e-01;
    C32 = 2.888318260750401e+01;
    C33 = -1.498095094187991e-01;
    C34 = 2.786230803523323e-10;
    C41 = -5.663110230820932e-02;
    C42 = 2.610281617868567e+01;
    C43 = 9.299717065949033e-02;
    C44 = 1.000699526693422e-10;
    C51 = 2.026322104668477e-03;
    C52 = 2.642095834119355e+01;
    C53 = -3.045407453585553e-03;
    C54 = 1.063550677583818e-10;
    C0111 = 3.143556106934323e+01;
    C0112 = -1.487080005126687e+03;
    C0121 = -1.102475588253000e+02;
    C0113 = 1.852393971803182e-10;
    C0122 = 9.152746138862564e+02;
    C0131 = 7.782877378723002e+01;
    C1211 = -2.451659125160949e+01;
    C1221 = -1.535931406766860e+03;
    C1212 = 2.470261842882782e+01;
    C1231 = 1.011893635288147e-10;
    C1222 = 9.673853234531789e+02;
    C1213 = 1.001171398616361e-10;
    C0211 = -8.801348524068304e+02;
    C0212 = 4.524435043734538e+02;
    C0221 = 3.939863694117493e+02;
    C0213 = 1.635681305524175e-10;
    C0222 = 1.325322397792186e-09;
    C0231 = 1.919430648834414e-09;
    C0311 = -6.424536673856824e-01;
    C0312 = -1.250521527596394e+02;
    C0321 = 5.485466573484021e-01;
    C0313 = 1.765820249103039e-01;
    C0322 = 8.257825518038538e+01;
    C0331 = 3.555862482756407e-10;
    C2511 = -7.313580505978607e-03;
    C2512 = -1.071740001877660e+02;
    C2521 = 9.939121989008142e-04;
    C2513 = 9.654481160453256e-03;
    C2522 = 6.951637017003624e+01;
    C2531 = 1.030976672908113e-03;
    C3511 = -2.691691636444475e+00;
    C3512 = -1.603864361964612e+00;
    C3521 = -5.349224439521996e-02;
    C3513 = 2.761851512903444e-10;
    C3522 = 2.041718060257125e+03;
    C3531 = 1.332451913911396e-09;
    Ccompr0 = 2.523471117547179e+02;
    Ccompr1 = 1.777509093910859e-01;
  }
}

bool FittedMaterial::clamp_strains(Vec6 &ek, std::vector<int> &clamped_coords,
                                   std::vector<double> &dek)
{
  bool clamped = false;
  for (int i = 0; i < 6; ++i)
  {
    if (ek[i] < ek_min[i])
    {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_min[i]); // (ek - ekclamped)
      ek[i] = ek_min[i];
      clamped = true;
      // std::cout<<"clamping min "<<i<<" w distance "<<ek[i] - ek_min[i]<<"\n";
    }
    else if (ek[i] > ek_max[i])
    {
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
double FittedMaterial::psi(const Vec6 &ek)
{
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

Vec6 FittedMaterial::psi_grad(const Vec6 &ek)
{
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

std::pair<Mat6x6, Vec6> FittedMaterial::psi_drv(const Vec6 &ek)
{

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
