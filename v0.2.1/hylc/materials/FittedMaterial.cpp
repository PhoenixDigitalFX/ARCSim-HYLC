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
  if (type == 0) {  // stockinette
    strain_min     = {0.4,  -0.8, 0.4,
                  -200, -100, -250};  // TODO measure bounds vs barrier?
    strain_max     = {2.0, 0.8, 2.0, 200, 100, 250};
    density        = 3.2949044794e-01;
    strainscale[0] = 1.000000000000000e+00;
    strainscale[1] = 7.505106984103984e-01;
    strainscale[2] = 1.000000000000000e+00;
    strainscale[3] = 3.998970322000000e+02;
    strainscale[4] = 1.050047476000000e+02;
    strainscale[5] = 4.997989018000000e+02;
    C0             = -1.003818644749069e+01;
    C01            = 7.816830717858565e+00;
    C02            = 8.136266447997481e-01;
    C03            = 3.816425231857323e+00;
    C04            = 1.000000000000000e-10;
    C11            = -3.438291685847008e-02;
    C12            = 1.833470497060952e-01;
    C13            = 7.110645208475560e-03;
    C14            = 1.103235967349711e+00;
    C21            = 7.841967484398829e+00;
    C22            = 4.090128178556819e+00;
    C23            = -1.018778362581240e+00;
    C24            = 2.822261419714615e+00;
    C31            = -7.619133154361377e-01;
    C32            = 1.857687857070609e+00;
    C33            = -1.348992387575448e+00;
    C34            = -1.458056659230248e-01;
    C41            = -8.892665534985685e-03;
    C42            = 5.896720647659408e-02;
    C43            = -1.666882621413648e-02;
    C44            = 3.076289004215776e-01;
    C51            = 9.899500846260979e-01;
    C52            = 3.149893126915625e+00;
    C53            = 2.879352421143433e+00;
    C54            = 1.147999667676190e+00;
    C0111          = -1.021214417322628e-01;
    C0112          = 9.740590829140793e-01;
    C0121          = 2.250637988511391e-01;
    C0113          = 1.000000000000000e-10;
    C0122          = 1.000000000000000e-10;
    C0131          = 1.000000000000000e-10;
    C0211          = 6.722905769602443e-01;
    C0212          = 7.788688595759385e+00;
    C0221          = 1.042145593617905e+01;
    C0213          = 1.000000000000000e-10;
    C0222          = 5.558806219882548e+00;
    C0231          = 1.000000000000000e-10;
    C0311          = -3.858563468771280e-02;
    C0312          = -4.773333016197464e+00;
    C0321          = -5.248956158762187e+00;
    C0313          = 1.975866449386912e+00;
    C0322          = 3.645711722125020e+00;
    C0331          = 2.883357864029369e+00;
    C0411          = 1.580061372212539e-01;
    C0412          = -2.404475071934213e+00;
    C0421          = -5.126701087752591e-01;
    C0413          = 1.321829786299092e-02;
    C0422          = 1.000000000000000e-10;
    C0431          = 4.081663311156233e-01;
    C0511          = 1.416825391742409e+00;
    C0512          = -4.292946776283403e+00;
    C0521          = 8.705478931085874e+00;
    C0513          = 1.000000000000000e-10;
    C0522          = 9.024782777172598e+00;
    C0531          = 1.000000000000000e-10;
    C1211          = -2.299279451167805e-01;
    C1212          = 2.949702425863953e-01;
    C1221          = -1.840049316482592e+00;
    C1213          = 1.000000000000000e-10;
    C1222          = 3.101559964156303e+00;
    C1231          = 1.570238944369045e-01;
    C1311          = -1.907108387993680e-01;
    C1312          = 2.370755735413338e-02;
    C1321          = -1.149241667465178e+00;
    C1313          = 1.001715701559122e+00;
    C1322          = 3.760848541522382e+00;
    C1331          = 4.728774166570899e-01;
    C1411          = -4.620479781200793e-01;
    C1412          = 3.013334809222201e-02;
    C1421          = 9.203145956748006e-02;
    C1413          = 1.000000000000000e-10;
    C1422          = 1.814628014033463e+00;
    C1431          = 1.000000000000000e-10;
    C1511          = -1.593541508869527e-01;
    C1512          = -1.461569636128655e-01;
    C1521          = 2.885515807745717e+00;
    C1513          = 1.000000000000000e-10;
    C1522          = 5.310539498030637e+00;
    C1531          = 1.000000000000000e-10;
    C2311          = -3.901330307657239e+00;
    C2312          = 1.683676963209165e+00;
    C2321          = -1.049063776208416e+01;
    C2313          = 1.854063787334209e+00;
    C2322          = 1.000000000000000e-10;
    C2331          = 5.497803868745760e+00;
    C2411          = 1.574566484247732e-01;
    C2412          = -2.359146170292670e+00;
    C2421          = -3.455493157278433e-01;
    C2413          = 1.325978764881029e-03;
    C2422          = 2.433177259818043e-01;
    C2431          = 2.324936408195044e-01;
    C2511          = -1.958934109071078e+00;
    C2512          = -1.072121616422114e+01;
    C2521          = 3.971175207449036e+00;
    C2513          = 1.000000000000000e-10;
    C2522          = 8.997133583318595e+00;
    C2531          = 1.000000000000000e-10;
    C3411          = -2.811469646512593e-02;
    C3412          = -2.606702204469629e-01;
    C3421          = 2.491820856706644e-02;
    C3413          = 1.000000000000000e-10;
    C3422          = 1.000000000000000e-10;
    C3431          = 2.644997296536375e-02;
    C3511          = -6.967006024893937e+00;
    C3512          = -8.961944644965266e+00;
    C3521          = 4.424557077813521e+00;
    C3513          = 1.000000000000000e-10;
    C3522          = 1.318346780299234e+02;
    C3531          = 1.711103275246939e+01;
    C4511          = -5.805110347941198e-02;
    C4512          = 2.102420770185180e-01;
    C4521          = 9.416336966911285e-01;
    C4513          = 4.084829486416055e-01;
    C4522          = 1.180848801623229e+00;
    C4531          = 1.580686201576425e-01;
    Ccompr0        = 1.425726086604097e+01;
    Ccompr1        = 3.954034266555811e-01;
  }
}

bool FittedMaterial::clamp_strains(Vec6 &strain,
                                   std::vector<int> &clamped_coords,
                                   std::vector<double> &dstrain) {
  bool clamped = false;
  for (int i = 0; i < 6; ++i) {
    if (strain[i] < strain_min[i]) {
      clamped_coords.push_back(i);
      dstrain.push_back(strain[i] - strain_min[i]);  // (strain - strainclamped)
      strain[i] = strain_min[i];
      clamped   = true;
      // std::cout<<"clamping min "<<i<<" w distance "<<strain[i] -
      // strain_min[i]<<"\n";
    } else if (strain[i] > strain_max[i]) {
      clamped_coords.push_back(i);
      dstrain.push_back(strain[i] - strain_max[i]);  // (strain - strainclamped)
      strain[i] = strain_max[i];
      clamped   = true;
      // std::cout<<"clamping max "<<i<<" w distance "<<strain[i] -
      // strain_min[i]<<"\n";
    }
  }
  return clamped;
}

#include <iostream>
double FittedMaterial::psi(const Vec6 &strain) {
  Vec6 strainclamped = strain;
  std::vector<int> clamped_coords;
  clamped_coords.reserve(6);
  std::vector<double> dstrain;  // (ei - clamped(ei)) Taylor distance
  dstrain.reserve(6);

  // TODO input to
  // TODO input here
  clamp_strains(strainclamped, clamped_coords, dstrain);

  double out = 0;

  // if unclamped dstrain is 0 and only actual energy remains
  // psi(strain + dstrain) = psi(strain) + dpsid_i(strain) * dstrain_i + 0.5 *
  // dpsidd_i(strain) * dstrain_i^2
  out = psi_taylor_0(
      use_barrier
          ? strainclamped
          : strain);  //  const taylor part, actual energy if not clamped

  for (int i = 3; i < 6; i++)
    strainclamped(i) =
        strain(i);  // pretend that distance is 0 for bending strains

  if (use_barrier)
    out += psi_barrier(strain, strainclamped);

  return out;
}

Vec6 FittedMaterial::psi_grad(const Vec6 &strain) {
  Vec6 strainclamped = strain;
  std::vector<int> clamped_coords;
  clamped_coords.reserve(6);
  std::vector<double> dstrain;  // (ei - clamped(ei)) Taylor distance
  dstrain.reserve(6);
  clamp_strains(strainclamped, clamped_coords, dstrain);

  Vec6 out(0);
  out = grad_taylor_0(
      use_barrier
          ? strainclamped
          : strain);  //  const taylor part, actual energy if not clamped

  for (int i = 3; i < 6; i++)
    strainclamped(i) =
        strain(i);  // pretend that distance is 0 for bending strains
  if (use_barrier)
    out += grad_barrier(strain, strainclamped);

  return out;
}

std::pair<Mat6x6, Vec6> FittedMaterial::psi_drv(const Vec6 &strain) {
  Vec6 strainclamped = strain;
  std::vector<int> clamped_coords;
  clamped_coords.reserve(6);
  std::vector<double> dstrain;  // (ei - clamped(ei)) Taylor distance
  dstrain.reserve(6);
  clamp_strains(strainclamped, clamped_coords, dstrain);

  Mat6x6 hess(0);
  Vec6 grad(0);

  auto gh0 = gradhess_taylor_0(
      use_barrier
          ? strainclamped
          : strain);  //  const taylor part, actual energy if not clamped
  hess += gh0.first;
  grad += gh0.second;

  for (int i = 3; i < 6; i++)
    strainclamped(i) =
        strain(i);  // pretend that distance is 0 for bending strains
  if (use_barrier) {
    auto ghB = gradhess_barrier(strain, strainclamped);
    hess += ghB.first;
    grad += ghB.second;
  }

  return std::make_pair(hess, grad);
}
