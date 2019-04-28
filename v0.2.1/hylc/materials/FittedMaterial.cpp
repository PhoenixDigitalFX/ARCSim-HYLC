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
    strain_min     = {0.0,  -0.75, 0.0,
                  -200, -90, -250};  // TODO measure bounds vs barrier?
    strain_max     = {2.0, 0.75, 2.0, 200, 90, 250};
    density        = 3.2949044794e-01;
    strainscale[0] = 1.000000000000000e+00;
    strainscale[1] = 7.505306426861681e-01;
    strainscale[2] = 1.000000000000000e+00;
    strainscale[3] = 3.998970322000000e+02;
    strainscale[4] = 1.024603321500000e+02;
    strainscale[5] = 4.997989018000000e+02;
    C0             = -2.583689260735773e+00;
    C01            = 3.761998302468429e+00;
    C02            = 1.740168952451396e+00;
    C03            = 4.489256381779957e+00;
    C04            = 1.054010877541357e-10;
    C11            = -3.740007076789804e-02;
    C12            = 8.673820771314596e-02;
    C13            = 1.308112329237473e-02;
    C14            = 1.306703020353469e+00;
    C21            = 3.693537888883060e+00;
    C22            = 4.414594428294858e+00;
    C23            = 2.674128212960473e+00;
    C24            = 1.000035774816541e-10;
    C31            = -7.665672381935169e-01;
    C32            = 1.707426361480081e+00;
    C33            = -1.463662514797240e+00;
    C34            = 3.559798313215725e-01;
    C41            = -5.010268532813047e-04;
    C42            = 3.319576164293321e-02;
    C43            = -2.430839910646501e-02;
    C44            = 2.757919236804862e-01;
    C51            = 8.873374771626198e-01;
    C52            = 3.230308144865501e+00;
    C53            = 3.427972875943248e+00;
    C54            = 6.152449009696280e-01;
    C0111          = 4.171584267102617e-01;
    C0112          = 1.707173457448742e+00;
    C0121          = -3.386373075515804e-01;
    C0113          = -6.861322837804851e-01;
    C0122          = -1.706680139522157e+00;
    C0131          = 5.403988524739924e-01;
    C0211          = -1.676924921842356e+00;
    C0212          = 1.553907961666330e+01;
    C0221          = 2.501019875607848e+01;
    C0213          = -5.779103669932293e+00;
    C0222          = 4.866524976267702e-01;
    C0231          = -1.252292735785256e+01;
    C0311          = 8.710451148684276e-03;
    C0312          = -5.846396639234984e+00;
    C0321          = -3.519291027498746e+00;
    C0313          = 2.408471660966941e+00;
    C0322          = 4.679118555622282e+00;
    C0331          = 2.300914829448796e-01;
    C0411          = 7.571057228137787e-02;
    C0412          = -1.135001499828273e+00;
    C0421          = -8.304328388402165e-01;
    C0413          = 6.967091269622792e-02;
    C0422          = -2.306424774794670e+00;
    C0431          = 9.722686039759629e-01;
    C0511          = 1.210295385025024e+00;
    C0512          = -3.253424471700641e+00;
    C0521          = 8.900554409149258e+00;
    C0513          = -1.936622304480658e-01;
    C0522          = 7.862913386105775e+00;
    C0531          = 6.417624882033456e-01;
    C1211          = -1.044743758907113e-01;
    C1212          = 3.349252888781991e-01;
    C1221          = -2.020358761638191e+00;
    C1213          = -5.693881409578808e-01;
    C1222          = 3.736251110367001e+00;
    C1231          = 2.967851450520432e-01;
    C1311          = -3.625090064607752e-01;
    C1312          = 2.152549469332077e-02;
    C1321          = -1.177239526713471e+00;
    C1313          = 2.041932592961278e+00;
    C1322          = 4.395507855368609e+00;
    C1331          = 4.341999128686170e-01;
    C1411          = -3.160782989337617e-01;
    C1412          = -9.950815010491254e-03;
    C1421          = -1.294973189910542e-02;
    C1413          = -3.503993179176550e-02;
    C1422          = 1.930533646871038e+00;
    C1431          = -3.240298894532069e-01;
    C1511          = -8.490941940971408e-02;
    C1512          = -5.270133972996015e-02;
    C1521          = 3.027133189885599e+00;
    C1513          = 1.684908653096500e-01;
    C1522          = 4.694754874836040e+00;
    C1531          = -1.234829743702751e-01;
    C2311          = -3.171620480348300e+00;
    C2312          = 1.043480421999083e+00;
    C2321          = -1.104474670395795e+01;
    C2313          = 1.564664498103869e-01;
    C2322          = 2.975236001868980e+00;
    C2331          = 4.630370343453738e+00;
    C2411          = 1.247635578149544e-01;
    C2412          = -3.421494265040022e+00;
    C2421          = 2.642946175133767e-01;
    C2413          = -2.469996598738923e-01;
    C2422          = 2.751615183082597e+00;
    C2431          = -3.346316923561914e-01;
    C2511          = -2.236360814922877e+00;
    C2512          = -5.277448714562856e+00;
    C2521          = 1.785007439983631e+01;
    C2513          = -5.232772757609895e+00;
    C2522          = 5.376176398426669e-01;
    C2531          = -1.438460026624184e+01;
    C3411          = 1.878743851371322e-02;
    C3412          = -2.785027024493537e-01;
    C3421          = -2.348556415167639e-03;
    C3413          = -1.076615566473933e-01;
    C3422          = 3.040568516771601e-01;
    C3431          = 1.309279110309722e-01;
    C3511          = -2.335739637330160e+00;
    C3512          = -1.205410846394313e+01;
    C3521          = 7.522665607148684e+00;
    C3513          = -3.898998809125761e+01;
    C3522          = 1.475879939631389e+02;
    C3531          = -1.189414787344278e+01;
    C4511          = -4.845292930072063e-02;
    C4512          = 2.076441744044206e-01;
    C4521          = 8.614382976789762e-01;
    C4513          = 3.512397517658248e-01;
    C4522          = 1.367499350232036e+00;
    C4531          = 1.667597665425031e-01;
    Ccompr0        = 2.369562846086346e+00;
    Ccompr1        = -1.436011891220309e-01;
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
