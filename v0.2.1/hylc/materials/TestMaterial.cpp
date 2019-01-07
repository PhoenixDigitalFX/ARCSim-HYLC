#include "TestMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

double TestMaterial::psi(const Vec6 &_ek) {
  // define input
  auto ek = [&](int i) -> const Real & { return _ek[i]; };
  auto C = [&](int i) -> const Real & { return m_C[i]; };
  auto ekstd = [&](int i) -> const Real & { return m_ekstd[i]; };
  auto ekmean = [&](int i) -> const Real & { return m_ekmean[i]; };
  Real c9 = ek(0);
  Real c10 = ekmean(0);
  Real c11 = -c10;
  Real c12 = c11 + c9;
  Real c15 = ekstd(0);
  Real c31 = ek(1);
  Real c32 = ekmean(1);
  Real c33 = -c32;
  Real c34 = c31 + c33;
  Real c36 = ekstd(1);
  Real c52 = ek(2);
  Real c53 = ekmean(2);
  Real c54 = -c53;
  Real c55 = c52 + c54;
  Real c57 = ekstd(2);
  Real c72 = ek(3);
  Real c73 = ekmean(3);
  Real c74 = -c73;
  Real c75 = c72 + c74;
  Real c77 = ekstd(3);
  Real c92 = ek(4);
  Real c93 = ekmean(4);
  Real c94 = -c93;
  Real c95 = c92 + c94;
  Real c97 = ekstd(4);
  Real c112 = ek(5);
  Real c113 = ekmean(5);
  Real c114 = -c113;
  Real c115 = c112 + c114;
  Real c117 = ekstd(5);
  return s *
         (C(0) + C(1) + (c12 * C(2)) / c15 + (c12 * c12 * C(3)) / c15 * c15 +
          (Power(c12, 3) * C(4)) / Power(c15, 3) +
          (Power(c12, 4) * C(5)) / Power(c15, 4) + C(6) + (c34 * C(7)) / c36 +
          (c34 * c34 * C(8)) / c36 * c36 +
          (Power(c34, 3) * C(9)) / Power(c36, 3) +
          (Power(c34, 4) * C(10)) / Power(c36, 4) + C(11) +
          (c55 * C(12)) / c57 + (c55 * c55 * C(13)) / c57 * c57 +
          (Power(c55, 3) * C(14)) / Power(c57, 3) +
          (Power(c55, 4) * C(15)) / Power(c57, 4) + C(16) +
          (c75 * C(17)) / c77 + (c75 * c75 * C(18)) / c77 * c77 +
          (Power(c75, 3) * C(19)) / Power(c77, 3) +
          (Power(c75, 4) * C(20)) / Power(c77, 4) + C(21) +
          (c95 * C(22)) / c97 + (c95 * c95 * C(23)) / c97 * c97 +
          (Power(c95, 3) * C(24)) / Power(c97, 3) +
          (Power(c95, 4) * C(25)) / Power(c97, 4) + C(26) +
          (c115 * C(27)) / c117 + (c115 * c115 * C(28)) / c117 * c117 +
          (Power(c115, 3) * C(29)) / Power(c117, 3) +
          (Power(c115, 4) * C(30)) / Power(c117, 4));
}

std::pair<Mat6x6, Vec6> TestMaterial::psi_drv(const Vec6 &_ek) {
  // define input
  auto ek = [&](int i) -> const Real & { return _ek[i]; };
  auto C = [&](int i) -> const Real & { return m_C[i]; };
  auto ekstd = [&](int i) -> const Real & { return m_ekstd[i]; };
  auto ekmean = [&](int i) -> const Real & { return m_ekmean[i]; };

  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real c1 = ekstd(0);
  Real c4 = ek(0);
  Real c5 = ekmean(0);
  Real c6 = -c5;
  Real c7 = c4 + c6;
  Real c25 = ekstd(1);
  Real c28 = ek(1);
  Real c29 = ekmean(1);
  Real c30 = -c29;
  Real c31 = c28 + c30;
  Real c48 = ekstd(2);
  Real c51 = ek(2);
  Real c52 = ekmean(2);
  Real c53 = -c52;
  Real c54 = c51 + c53;
  Real c70 = ekstd(3);
  Real c73 = ek(3);
  Real c74 = ekmean(3);
  Real c75 = -c74;
  Real c76 = c73 + c75;
  Real c92 = ekstd(4);
  Real c95 = ek(4);
  Real c96 = ekmean(4);
  Real c97 = -c96;
  Real c98 = c95 + c97;
  Real c114 = ekstd(5);
  Real c117 = ek(5);
  Real c118 = ekmean(5);
  Real c119 = -c118;
  Real c120 = c117 + c119;
  Real c2 = Power(c1, -4);
  Real c3 = C(5);
  Real c11 = c7 * c7;
  Real c10 = C(4);
  Real c14 = C(3);
  Real c26 = Power(c25, -4);
  Real c27 = C(10);
  Real c35 = c31 * c31;
  Real c34 = C(9);
  Real c37 = C(8);
  Real c49 = Power(c48, -4);
  Real c50 = C(15);
  Real c58 = c54 * c54;
  Real c57 = C(14);
  Real c60 = C(13);
  Real c71 = Power(c70, -4);
  Real c72 = C(20);
  Real c80 = c76 * c76;
  Real c79 = C(19);
  Real c82 = C(18);
  Real c93 = Power(c92, -4);
  Real c94 = C(25);
  Real c102 = c98 * c98;
  Real c101 = C(24);
  Real c104 = C(23);
  Real c115 = Power(c114, -4);
  Real c116 = C(30);
  Real c124 = c120 * c120;
  Real c123 = C(29);
  Real c126 = C(28);
  out1(0) = s * c2 *
            (4 * c3 * Power(c7, 3) +
             c1 * (3 * c10 * c11 + c1 * (2 * c14 * c7 + c1 * C(2))));
  out1(1) = s * c26 *
            (4 * c27 * Power(c31, 3) +
             c25 * (3 * c34 * c35 + c25 * (2 * c31 * c37 + c25 * C(7))));
  out1(2) = s * c49 *
            (4 * c50 * Power(c54, 3) +
             c48 * (3 * c57 * c58 + c48 * (2 * c54 * c60 + c48 * C(12))));
  out1(3) = s * c71 *
            (4 * c72 * Power(c76, 3) +
             c70 * (3 * c79 * c80 + c70 * (2 * c76 * c82 + c70 * C(17))));
  out1(4) = s * c93 *
            (4 * c94 * Power(c98, 3) +
             c92 * (3 * c101 * c102 + c92 * (2 * c104 * c98 + c92 * C(22))));
  out1(5) =
      s * c115 *
      (4 * c116 * Power(c120, 3) +
       c114 * (3 * c123 * c124 + c114 * (2 * c120 * c126 + c114 * C(27))));
  out2(0, 0) = 2 * s * c2 * (6 * c11 * c3 + c1 * (c1 * c14 + 3 * c10 * c7));
  out2(0, 1) = 0;
  out2(0, 2) = 0;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) =
      2 * s * c26 * (6 * c27 * c35 + c25 * (3 * c31 * c34 + c25 * c37));
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) =
      2 * s * c49 * (6 * c50 * c58 + c48 * (3 * c54 * c57 + c48 * c60));
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) =
      2 * s * c71 * (6 * c72 * c80 + c70 * (3 * c76 * c79 + c70 * c82));
  out2(3, 4) = 0;
  out2(3, 5) = 0;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = 0;
  out2(4, 4) =
      2 * s * c93 * (6 * c102 * c94 + c92 * (c104 * c92 + 3 * c101 * c98));
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) =
      2 * s * c115 * (6 * c116 * c124 + c114 * (3 * c120 * c123 + c114 * c126));

  return std::make_pair(hess, grad);
}
