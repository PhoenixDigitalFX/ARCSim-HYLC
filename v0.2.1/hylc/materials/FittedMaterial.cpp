#include "FittedMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

#include <iostream>
FittedMaterial::FittedMaterial(int type) {
  std::cout << "Loading Material " << type << "\n";

  use_barrier = true;
  bspeed =
      1e0;  // barrier strength // TODO, also maybe barrier not 6th order TODO
  bscale = 1e2;
  // if barrier is 4th order and bspeed is 10 but strain is normalized same way
  // then it should grow 10^4 times faster ?
  if (type == 0) {  // slipstitchhoney
    ek_min = {0.0, -0.8, 0.0,
              -70, -40,  -70};  // TODO measure bounds vs barrier?
    ek_max = {4.0, 0.8, 4.0, 70, 40, 70};
    density = 4.1562049890e-01;
    ekscale[0] = 3.000000000000000e+00;
    ekscale[1] = 1.600000000000000e+00;
    ekscale[2] = 3.000000000000000e+00;
    ekscale[3] = 1.399920440000000e+02;
    ekscale[4] = 3.695105080000000e+01;
    ekscale[5] = 1.399918732000000e+02;
    C0 = -4.169937212160532e+00;
    C01 = 1.612808100445473e+01;
    C02 = 3.687606885814537e+01;
    C03 = 1.630694747600710e+01;
    C04 = 2.382401937689732e-10;
    C11 = -2.696929653981711e-01;
    C12 = 1.298474597931686e+01;
    C13 = 1.221015944845739e-01;
    C14 = 3.152803934637734e+01;
    C21 = 1.599013775207774e+01;
    C22 = 1.225140311416962e+01;
    C23 = 3.861195234079322e+01;
    C24 = 1.000000413344519e-10;
    C31 = 5.807048954209569e-01;
    C32 = 4.369084235726388e+00;
    C33 = -4.464239909360339e-01;
    C34 = 1.011055968390358e-10;
    C41 = -6.447241577011194e-02;
    C42 = 4.188815046414406e-01;
    C43 = 5.195768032876835e-02;
    C44 = 6.311461870581086e-10;
    C51 = 1.267692381948303e-01;
    C52 = 2.187005715058245e+00;
    C53 = 1.258602278927919e+00;
    C54 = 1.048079313749996e-08;
    C0111 = 1.819202064110867e+00;
    C0112 = -5.683547733791996e+01;
    C0121 = -8.703459205692377e+00;
    C0113 = 3.267696622189011e+00;
    C0122 = 2.316225729156167e+01;
    C0131 = 5.866995788428657e+00;
    C1211 = -8.663158087722103e+00;
    C1221 = -3.223459879222595e+01;
    C1212 = 6.170149957063725e+00;
    C1231 = 1.000000428847936e-10;
    C1222 = 1.110134248858557e+01;
    C1213 = 9.022681097841101e-04;
    C0211 = 3.603806796864024e+00;
    C0212 = 8.823914585655265e+01;
    C0221 = 1.442580909415224e+02;
    C0213 = 2.389283651131812e-10;
    C0222 = 7.241728050110792e+01;
    C0231 = 1.191287180634241e-10;
    C0311 = 9.701199240542470e+00;
    C0312 = -1.475752983966270e+01;
    C0321 = -1.289791719295748e+01;
    C0313 = 7.357423004545709e-01;
    C0322 = 9.135863266119332e+00;
    C0331 = 4.813927766264184e+00;
    C2511 = -7.477484862441563e+00;
    C2512 = -1.375353050806992e+01;
    C2521 = 1.082450464144567e+01;
    C2513 = 1.007176468148605e-10;
    C2522 = 1.650443367623757e+01;
    C2531 = 2.214301532827918e-10;
    C3511 = 2.484734681109843e-02;
    C3512 = -3.131788317313266e-01;
    C3521 = 3.026474438411697e+00;
    C3513 = 1.006060475324269e-10;
    C3522 = 2.777458751933085e+01;
    C3531 = 2.127733857341799e+00;
    Ccompr0 = 5.176879294184138e+00;
    Ccompr1 = 1.000086011346457e-10;
  } else if (type == 1) {  // slipstitchrib
    ek_min = {0.0,  -0.8, 0.0,
              -100, -60,  -100};  // TODO measure bounds vs barrier?
    ek_max = {4.0, 0.8, 4.0, 100, 60, 100};
    density = 6.0205506498e-01;
    ekscale[0] = 3.000000000000000e+00;
    ekscale[1] = 1.600000000000000e+00;
    ekscale[2] = 3.000000000000000e+00;
    ekscale[3] = 1.399937629000000e+02;
    ekscale[4] = 5.298510940500000e+01;
    ekscale[5] = 1.399941661000000e+02;
    C0 = 1.029232426866122e+02;
    C01 = 1.191745935979640e+03;
    C02 = -3.173938433105700e+02;
    C03 = 4.399282264251593e+01;
    C04 = 3.364269846107152e+01;
    C11 = -4.470402892408852e+02;
    C12 = 7.540800969630250e+02;
    C13 = -9.150615610164708e+01;
    C14 = 5.227934359829638e+02;
    C21 = 9.125387758749604e+02;
    C22 = 3.683619809727810e+02;
    C23 = -1.017761358871539e+02;
    C24 = 1.000000000218192e-10;
    C31 = 9.604282556903969e+01;
    C32 = 2.547566846228422e+01;
    C33 = 5.284782560282889e+00;
    C34 = 1.000000118406923e-10;
    C41 = 2.821465794057206e+01;
    C42 = 3.021184935572633e+01;
    C43 = 1.140401662571048e+01;
    C44 = 1.204346351549734e+01;
    C51 = 7.574720772853283e+01;
    C52 = 5.548120773849816e+01;
    C53 = -1.495175834219085e+01;
    C54 = 8.007829115015808e+00;
    C0111 = 1.658182483171825e+02;
    C0112 = -1.631396473587136e+03;
    C0121 = -4.805468481292214e+01;
    C0113 = 1.000000000016789e-10;
    C0122 = 7.240187142005940e+02;
    C0131 = 4.905703906455661e+01;
    C1211 = -2.634668910332831e+01;
    C1221 = -1.856686497864073e+03;
    C1212 = 2.059028413983595e+02;
    C1231 = 1.000000000001183e-10;
    C1222 = 7.665110215815808e+02;
    C1213 = 1.000000000000019e-10;
    C0211 = 3.917006731462708e+01;
    C0212 = 1.300930198181952e+02;
    C0221 = -5.615214350694951e+02;
    C0213 = 1.009432376385894e-10;
    C0222 = 1.000000010878550e-10;
    C0231 = 6.372973017250500e+02;
    C0311 = -2.927192195731965e+02;
    C0312 = -4.066032507997035e+01;
    C0321 = 2.225617919586869e+02;
    C0313 = 1.000000000066990e-10;
    C0322 = 4.375726429873126e+00;
    C0331 = 1.000000000000000e-10;
    C2511 = 3.466756349471953e+01;
    C2512 = -1.493363714184856e+02;
    C2521 = -4.921919401160355e+02;
    C2513 = 1.582546416049004e+01;
    C2522 = 7.459736872601562e+01;
    C2531 = 4.330321074414447e+02;
    C3511 = 7.243261136329862e+00;
    C3512 = -4.526147854470868e+01;
    C3521 = -7.455695002581437e+01;
    C3513 = 1.000000000001624e-10;
    C3522 = 3.068881056832164e+02;
    C3531 = 1.038977627184902e+02;
    Ccompr0 = 8.076289473979385e+02;
    Ccompr1 = 4.269709335102413e-01;
  } else if (type == 2) {  // basket2_2
    ek_min = {0.0,  -0.8, 0.0,
              -150, -100, -150};  // TODO measure bounds vs barrier?
    ek_max = {4.0, 0.8, 4.0, 150, 100, 150};
    density = 3.1553164892e-01;
    ekscale[0] = 3.000000000000000e+00;
    ekscale[1] = 1.600000000000000e+00;
    ekscale[2] = 3.000000000000000e+00;
    ekscale[3] = 1.499743416000000e+02;
    ekscale[4] = 7.744316544000000e+01;
    ekscale[5] = 1.499745997000000e+02;
    C0 = -4.510815718245586e+01;
    C01 = 1.953928325010298e+02;
    C02 = 4.120239656646683e+02;
    C03 = -5.282259110786783e+01;
    C04 = 1.000000000000000e-10;
    C11 = 5.979915676677089e-01;
    C12 = 1.991971299802532e+02;
    C13 = 1.409442526284929e-01;
    C14 = 4.280894091876473e+02;
    C21 = 1.729497528776026e+02;
    C22 = 4.367151166161981e+02;
    C23 = -8.401030727402784e+01;
    C24 = 1.000000000000000e-10;
    C31 = 4.020491323756263e-02;
    C32 = 7.745567172076953e+00;
    C33 = -2.477904583056343e-02;
    C34 = 1.000000000000000e-10;
    C41 = -4.564811508612415e-02;
    C42 = -2.235361032959430e+00;
    C43 = 7.033703148767138e-02;
    C44 = 1.201845468916092e+01;
    C51 = 5.756769720987511e-04;
    C52 = 7.847812170045665e+00;
    C53 = -8.444760759873492e-04;
    C54 = 1.000000000000000e-10;
    C0111 = -7.384709978660982e+00;
    C0112 = -8.253801497268571e+02;
    C0121 = 8.656278156852242e+00;
    C0113 = 1.997682447370252e+00;
    C0122 = 4.299856967505843e+02;
    C0131 = 1.000000000001371e-10;
    C1211 = 6.527951700635992e-01;
    C1221 = -8.407059558960646e+02;
    C1212 = -1.423731586176511e+01;
    C1231 = 6.099095377064363e+00;
    C1222 = 4.432667442845294e+02;
    C1213 = 1.258341621202148e+01;
    C0211 = -3.131727274519760e+02;
    C0212 = 2.040024233126967e+02;
    C0221 = 1.546442418516037e+02;
    C0213 = 1.000000000000000e-10;
    C0222 = 2.008849199032312e+01;
    C0231 = 1.000000000000000e-10;
    C0311 = -7.312993982949361e-01;
    C0312 = -3.333315231355029e+01;
    C0321 = 1.124414982931318e+00;
    C0313 = 1.598859872498802e-10;
    C0322 = 1.022980186682059e+01;
    C0331 = 1.000000143957413e-10;
    C2511 = -6.060814165801954e-03;
    C2512 = -2.963254257879662e+01;
    C2521 = 8.478667052160816e-03;
    C2513 = 1.290568221449840e-03;
    C2522 = 8.517601936631651e+00;
    C2531 = 2.621216946376302e-06;
    C3511 = -3.671722850048618e-01;
    C3512 = -3.563367171359974e-01;
    C3521 = -4.575887163026612e-02;
    C3513 = 1.000000000229910e-10;
    C3522 = 2.186750707073139e+02;
    C3531 = 1.000000000020131e-10;
    Ccompr0 = 5.928807828394078e+01;
    Ccompr1 = 7.530058687210900e-02;
  } else if (type == 3) {  // left diagonal
    ek_min = {0.0,  -0.8, 0.0,
              -100, -70, -100};  // TODO measure bounds vs barrier?
    ek_max = {4.0, 0.8, 4.0, 100, 70, 100}; // sampling bounds here..
    density = 4.5347352190e-01;
    ekscale[0] = 3.000000000000000e+00;
    ekscale[1] = 1.600000000000000e+00;
    ekscale[2] = 3.000000000000000e+00;
    ekscale[3] = 1.799623856000000e+02;
    ekscale[4] = 4.992170441500000e+01;
    ekscale[5] = 1.799645233000000e+02;
    C0 = -3.988314402795971e+00;
    C01 = 2.588641064361665e+01;
    C02 = 6.006611215906834e+01;
    C03 = -4.503855122966828e+00;
    C04 = 1.000000395115077e-10;
    C11 = 7.661379806894944e+00;
    C12 = 8.055135622827743e+01;
    C13 = -1.472037145052538e+00;
    C14 = 4.626967892827536e+01;
    C21 = 2.912391962941827e+01;
    C22 = 1.414853244971180e+02;
    C23 = -4.439170017645874e+01;
    C24 = 1.000001295856314e-10;
    C31 = -4.355681951400870e-01;
    C32 = 1.275804573166941e+01;
    C33 = -1.829104182730263e+00;
    C34 = 2.141613460438668e-01;
    C41 = -1.093504188006281e+00;
    C42 = 2.365366870816034e+00;
    C43 = 9.261533397746690e-01;
    C44 = 3.776556775768014e+00;
    C51 = 5.202602498051429e+00;
    C52 = 1.215716688913655e+01;
    C53 = 3.701555707504530e-01;
    C54 = 1.031991302312095e+00;
    C0111 = 1.461280107967107e+01;
    C0112 = -1.585261848759336e+02;
    C0121 = -3.406364210033444e+01;
    C0113 = 2.441975378111478e+00;
    C0122 = 7.023262677436136e+01;
    C0131 = 2.163594770363614e+01;
    C1211 = -3.065157313067541e+01;
    C1221 = -1.728051591108353e+02;
    C1212 = 1.481058361540472e+01;
    C1231 = 1.174778024351352e+01;
    C1222 = 9.451482914568031e+01;
    C1213 = 1.178953197703341e-10;
    C0211 = 3.266109709072965e+01;
    C0212 = 5.676760483324583e+01;
    C0221 = 1.162741632715654e+02;
    C0213 = 1.008799989763360e-10;
    C0222 = 4.655962802622021e+01;
    C0231 = 1.000576632684764e-10;
    C0311 = 8.329765257751728e+00;
    C0312 = -3.929480837273046e+01;
    C0321 = -1.836904509206550e+01;
    C0313 = 1.717292848693538e+00;
    C0322 = 2.837528172710601e+01;
    C0331 = 8.742948531856596e+00;
    C2511 = -9.663883369995627e+00;
    C2512 = -4.913933607927721e+01;
    C2521 = 1.528351249504757e+01;
    C2513 = 1.618531654567264e-08;
    C2522 = 4.127971141018681e+01;
    C2531 = 1.028434412497740e-10;
    C3511 = -1.081663676630940e+01;
    C3512 = -6.960296231228806e+00;
    C3521 = -1.131219737297940e+01;
    C3513 = 1.018060365818484e-10;
    C3522 = 4.389777240956852e+02;
    C3531 = 1.000283574855876e-10;
    Ccompr0 = 6.875353632997077e+00;
    Ccompr1 = 1.081150724220022e-10;
  }
}

bool FittedMaterial::clamp_strains(Vec6 &ek, std::vector<int> &clamped_coords,
                                   std::vector<double> &dek) {
  bool clamped = false;
  for (int i = 0; i < 6; ++i) {
    if (ek[i] < ek_min[i]) {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_min[i]);  // (ek - ekclamped)
      ek[i] = ek_min[i];
      clamped = true;
      // std::cout<<"clamping min "<<i<<" w distance "<<ek[i] - ek_min[i]<<"\n";
    } else if (ek[i] > ek_max[i]) {
      clamped_coords.push_back(i);
      dek.push_back(ek[i] - ek_max[i]);  // (ek - ekclamped)
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
  std::vector<double> dek;  // (ei - clamped(ei)) Taylor distance
  dek.reserve(6);
  clamp_strains(ekclamped, clamped_coords, dek);

  double out = 0;
  // std::cout<<ekclamped[0]<<" from "<<ek[0]<<"\n\n";

  // if unclamped dek is 0 and only actual energy remains
  // psi(ek + dek) = psi(ek) + dpsid_i(ek) * dek_i + 0.5 * dpsidd_i(ek) *
  // dek_i^2
  out = psi_taylor_0(
      use_barrier ? ekclamped
                  : ek);  //  const taylor part, actual energy if not clamped
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
  if (use_barrier) out += psi_barrier(ek, ekclamped);
  out += psi_compr(ek);

  return out;
}

Vec6 FittedMaterial::psi_grad(const Vec6 &ek) {
  Vec6 ekclamped = ek;
  std::vector<int> clamped_coords;
  clamped_coords.reserve(6);
  std::vector<double> dek;  // (ei - clamped(ei)) Taylor distance
  dek.reserve(6);
  clamp_strains(ekclamped, clamped_coords, dek);

  Vec6 out(0);
  out = grad_taylor_0(
      use_barrier ? ekclamped
                  : ek);  //  const taylor part, actual energy if not clamped
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
  if (use_barrier) out += grad_barrier(ek, ekclamped);
  out += grad_compr(ek);

  return out;
}

std::pair<Mat6x6, Vec6> FittedMaterial::psi_drv(const Vec6 &ek) {
  Vec6 ekclamped = ek;
  std::vector<int> clamped_coords;
  clamped_coords.reserve(6);
  std::vector<double> dek;  // (ei - clamped(ei)) Taylor distance
  dek.reserve(6);
  clamp_strains(ekclamped, clamped_coords, dek);

  Mat6x6 hess(0);
  Vec6 grad(0);

  auto gh0 = gradhess_taylor_0(
      use_barrier ? ekclamped
                  : ek);  //  const taylor part, actual energy if not clamped
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
  if (use_barrier) {
    auto ghB = gradhess_barrier(ek, ekclamped);
    hess += ghB.first;
    grad += ghB.second;
  }
  auto ghC = gradhess_compr(ek);
  hess += ghC.first;
  grad += ghC.second;

  return std::make_pair(hess, grad);
}
