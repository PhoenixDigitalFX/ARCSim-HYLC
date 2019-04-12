#include "FittedMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

typedef Vec6 strain_type;
typedef double value_type;
typedef Vec6 grad_type;
typedef std::pair<Mat6x6, Vec6> gradhess_type;

// 0th derivative of psi, grad, tuple<grad,hess>
value_type FittedMaterial::psi_taylor_0(const strain_type &ek) {
  Real c7 = ek(0);
  Real c9 = -1 + c7;
  Real c14 = ekscale(0);
  Real c46 = ek(1);
  Real c50 = ekscale(1);
  Real c63 = Power(c46, 3);
  Real c40 = 1 / c14;
  Real c66 = Power(c50, -3);
  Real c28 = Power(c9, 2);
  Real c76 = Power(c46, 2);
  Real c30 = Power(c14, -2);
  Real c80 = Power(c50, -2);
  Real c22 = Power(c9, 3);
  Real c25 = Power(c14, -3);
  Real c89 = 1 / c50;
  Real c106 = ek(2);
  Real c108 = -1 + c106;
  Real c118 = ekscale(2);
  Real c125 = Power(c108, 3);
  Real c127 = Power(c118, -3);
  Real c134 = Power(c108, 2);
  Real c136 = Power(c118, -2);
  Real c149 = 1 / c118;
  Real c161 = ek(3);
  Real c165 = ekscale(3);
  Real c171 = Power(c161, 3);
  Real c173 = Power(c165, -3);
  Real c183 = Power(c161, 2);
  Real c185 = Power(c165, -2);
  Real c201 = 1 / c165;
  Real c222 = ek(4);
  Real c226 = ekscale(4);
  Real c231 = Power(c222, 3);
  Real c233 = Power(c226, -3);
  Real c240 = Power(c222, 2);
  Real c241 = Power(c226, -2);
  Real c256 = 1 / c226;
  Real c283 = ek(5);
  Real c287 = ekscale(5);
  Real c293 = Power(c283, 3);
  Real c295 = Power(c287, -3);
  Real c307 = Power(c283, 2);
  Real c308 = Power(c287, -2);
  Real c321 = 1 / c287;
  return C0 + (C24 * Power(c108, 4)) / Power(c118, 4) + C23 * c125 * c127 +
         C22 * c134 * c136 + C21 * c108 * c149 +
         (C34 * Power(c161, 4)) / Power(c165, 4) + C33 * c171 * c173 +
         C2313 * c108 * c149 * c171 * c173 + C32 * c183 * c185 +
         C2322 * c134 * c136 * c183 * c185 + C2312 * c108 * c149 * c183 * c185 +
         C31 * c161 * c201 + C2331 * c125 * c127 * c161 * c201 +
         C2321 * c134 * c136 * c161 * c201 + C2311 * c108 * c149 * c161 * c201 +
         (C44 * Power(c222, 4)) / Power(c226, 4) + C43 * c231 * c233 +
         C2413 * c108 * c149 * c231 * c233 + C3413 * c161 * c201 * c231 * c233 +
         C42 * c240 * c241 + C2422 * c134 * c136 * c240 * c241 +
         C2412 * c108 * c149 * c240 * c241 + C3422 * c183 * c185 * c240 * c241 +
         C3412 * c161 * c201 * c240 * c241 + C03 * c22 * c25 +
         C0231 * c108 * c149 * c22 * c25 + C0331 * c161 * c201 * c22 * c25 +
         C41 * c222 * c256 + C2431 * c125 * c127 * c222 * c256 +
         C2421 * c134 * c136 * c222 * c256 + C2411 * c108 * c149 * c222 * c256 +
         C3431 * c171 * c173 * c222 * c256 + C3421 * c183 * c185 * c222 * c256 +
         C3411 * c161 * c201 * c222 * c256 + C0431 * c22 * c222 * c25 * c256 +
         (C54 * Power(c283, 4)) / Power(c287, 4) + C53 * c293 * c295 +
         C2513 * c108 * c149 * c293 * c295 + C3513 * c161 * c201 * c293 * c295 +
         C4513 * c222 * c256 * c293 * c295 + C02 * c28 * c30 +
         C0222 * c134 * c136 * c28 * c30 + C0221 * c108 * c149 * c28 * c30 +
         C0322 * c183 * c185 * c28 * c30 + C0321 * c161 * c201 * c28 * c30 +
         C0422 * c240 * c241 * c28 * c30 + C0421 * c222 * c256 * c28 * c30 +
         C52 * c307 * c308 + C2522 * c134 * c136 * c307 * c308 +
         C2512 * c108 * c149 * c307 * c308 + C3522 * c183 * c185 * c307 * c308 +
         C3512 * c161 * c201 * c307 * c308 + C4522 * c240 * c241 * c307 * c308 +
         C4512 * c222 * c256 * c307 * c308 + C0522 * c28 * c30 * c307 * c308 +
         C51 * c283 * c321 + C2531 * c125 * c127 * c283 * c321 +
         C2521 * c134 * c136 * c283 * c321 + C2511 * c108 * c149 * c283 * c321 +
         C3531 * c171 * c173 * c283 * c321 + C3521 * c183 * c185 * c283 * c321 +
         C3511 * c161 * c201 * c283 * c321 + C4531 * c231 * c233 * c283 * c321 +
         C4521 * c240 * c241 * c283 * c321 + C0531 * c22 * c25 * c283 * c321 +
         C4511 * c222 * c256 * c283 * c321 + C0521 * c28 * c283 * c30 * c321 +
         (C14 * Power(c46, 4)) / Power(c50, 4) + C13 * c63 * c66 +
         C1231 * c108 * c149 * c63 * c66 + C1331 * c161 * c201 * c63 * c66 +
         C1431 * c222 * c256 * c63 * c66 + C1531 * c283 * c321 * c63 * c66 +
         C12 * c76 * c80 + C1222 * c134 * c136 * c76 * c80 +
         C1221 * c108 * c149 * c76 * c80 + C1322 * c183 * c185 * c76 * c80 +
         C1321 * c161 * c201 * c76 * c80 + C1422 * c240 * c241 * c76 * c80 +
         C1421 * c222 * c256 * c76 * c80 + C0122 * c28 * c30 * c76 * c80 +
         C1522 * c307 * c308 * c76 * c80 + C1521 * c283 * c321 * c76 * c80 +
         C11 * c46 * c89 + C1213 * c125 * c127 * c46 * c89 +
         C1212 * c134 * c136 * c46 * c89 + C1211 * c108 * c149 * c46 * c89 +
         C1313 * c171 * c173 * c46 * c89 + C1312 * c183 * c185 * c46 * c89 +
         C1311 * c161 * c201 * c46 * c89 + C1413 * c231 * c233 * c46 * c89 +
         C1412 * c240 * c241 * c46 * c89 + C0131 * c22 * c25 * c46 * c89 +
         C1411 * c222 * c256 * c46 * c89 + C1513 * c293 * c295 * c46 * c89 +
         C0121 * c28 * c30 * c46 * c89 + C1512 * c307 * c308 * c46 * c89 +
         C1511 * c283 * c321 * c46 * c89 + C01 * c40 * c9 +
         C0213 * c125 * c127 * c40 * c9 + C0212 * c134 * c136 * c40 * c9 +
         C0211 * c108 * c149 * c40 * c9 + C0313 * c171 * c173 * c40 * c9 +
         C0312 * c183 * c185 * c40 * c9 + C0311 * c161 * c201 * c40 * c9 +
         C0413 * c231 * c233 * c40 * c9 + C0412 * c240 * c241 * c40 * c9 +
         C0411 * c222 * c256 * c40 * c9 + C0513 * c293 * c295 * c40 * c9 +
         C0512 * c307 * c308 * c40 * c9 + C0511 * c283 * c321 * c40 * c9 +
         C0113 * c40 * c63 * c66 * c9 + C0112 * c40 * c76 * c80 * c9 +
         C0111 * c40 * c46 * c89 * c9 + (C04 * Power(c9, 4)) / Power(c14, 4);
}

grad_type FittedMaterial::grad_taylor_0(const strain_type &ek) {
  Vec6 out(0);

  Real c11 = ek(0);
  Real c14 = -1 + c11;
  Real c7 = ekscale(0);
  Real c30 = Power(c7, 3);
  Real c40 = ek(1);
  Real c27 = Power(c7, 2);
  Real c46 = ekscale(1);
  Real c57 = Power(c40, 2);
  Real c61 = Power(c46, -2);
  Real c22 = Power(c14, 2);
  Real c68 = 1 / c46;
  Real c82 = ek(2);
  Real c85 = -1 + c82;
  Real c89 = ekscale(2);
  Real c101 = Power(c85, 2);
  Real c104 = Power(c89, -2);
  Real c114 = 1 / c89;
  Real c125 = ek(3);
  Real c128 = ekscale(3);
  Real c134 = Power(c125, 2);
  Real c136 = Power(c128, -2);
  Real c142 = 1 / c128;
  Real c151 = ek(4);
  Real c154 = ekscale(4);
  Real c159 = Power(c151, 2);
  Real c160 = Power(c154, -2);
  Real c165 = 1 / c154;
  Real c173 = ek(5);
  Real c177 = ekscale(5);
  Real c183 = Power(c173, 2);
  Real c185 = Power(c177, -2);
  Real c191 = 1 / c177;
  Real c43 = Power(c40, 3);
  Real c211 = Power(c46, 2);
  Real c209 = 1 / c7;
  Real c17 = Power(c14, 3);
  Real c222 = Power(c46, 3);
  Real c216 = Power(c7, -2);
  Real c87 = Power(c85, 3);
  Real c95 = Power(c89, -3);
  Real c127 = Power(c125, 3);
  Real c129 = Power(c128, -3);
  Real c153 = Power(c151, 3);
  Real c155 = Power(c154, -3);
  Real c175 = Power(c173, 3);
  Real c179 = Power(c177, -3);
  Real c283 = Power(c89, 2);
  Real c226 = Power(c7, -3);
  Real c295 = Power(c89, 3);
  Real c48 = Power(c46, -3);
  Real c345 = Power(c128, 2);
  Real c360 = Power(c128, 3);
  Real c400 = Power(c154, 2);
  Real c419 = Power(c154, 3);
  Real c457 = Power(c177, 2);
  Real c469 = Power(c177, 3);
  out(0) =
      (4 * C04 * c17 + 2 * C02 * c14 * c27 +
       2 * C0222 * c101 * c104 * c14 * c27 +
       2 * C0322 * c134 * c136 * c14 * c27 +
       2 * C0321 * c125 * c14 * c142 * c27 +
       2 * C0422 * c14 * c159 * c160 * c27 +
       2 * C0421 * c14 * c151 * c165 * c27 +
       2 * C0522 * c14 * c183 * c185 * c27 +
       2 * C0521 * c14 * c173 * c191 * c27 + C01 * c30 +
       C0212 * c101 * c104 * c30 + C0313 * c127 * c129 * c30 +
       C0312 * c134 * c136 * c30 + C0311 * c125 * c142 * c30 +
       C0413 * c153 * c155 * c30 + C0412 * c159 * c160 * c30 +
       C0411 * c151 * c165 * c30 + C0513 * c175 * c179 * c30 +
       C0512 * c183 * c185 * c30 + C0511 * c173 * c191 * c30 +
       C0113 * c30 * c43 * c48 + 2 * C0122 * c14 * c27 * c57 * c61 +
       C0112 * c30 * c57 * c61 + 2 * C0121 * c14 * c27 * c40 * c68 +
       C0111 * c30 * c40 * c68 + 3 * C03 * c22 * c7 +
       3 * C0331 * c125 * c142 * c22 * c7 + 3 * C0431 * c151 * c165 * c22 * c7 +
       3 * C0531 * c173 * c191 * c22 * c7 + 3 * C0131 * c22 * c40 * c68 * c7 +
       2 * C0221 * c114 * c14 * c27 * c85 + C0211 * c114 * c30 * c85 +
       3 * C0231 * c114 * c22 * c7 * c85 + C0213 * c30 * c87 * c95) /
      Power(c7, 4);
  out(1) = (C11 * c222 + C1212 * c101 * c104 * c222 +
            C1313 * c127 * c129 * c222 + C1312 * c134 * c136 * c222 +
            C1311 * c125 * c142 * c222 + C1413 * c153 * c155 * c222 +
            C1412 * c159 * c160 * c222 + C1411 * c151 * c165 * c222 +
            C1513 * c175 * c179 * c222 + C1512 * c183 * c185 * c222 +
            C1511 * c173 * c191 * c222 + C0111 * c14 * c209 * c222 +
            C0121 * c216 * c22 * c222 + C0131 * c17 * c222 * c226 +
            2 * C12 * c211 * c40 + 2 * C1222 * c101 * c104 * c211 * c40 +
            2 * C1322 * c134 * c136 * c211 * c40 +
            2 * C1321 * c125 * c142 * c211 * c40 +
            2 * C1422 * c159 * c160 * c211 * c40 +
            2 * C1421 * c151 * c165 * c211 * c40 +
            2 * C1522 * c183 * c185 * c211 * c40 +
            2 * C1521 * c173 * c191 * c211 * c40 +
            2 * C0112 * c14 * c209 * c211 * c40 +
            2 * C0122 * c211 * c216 * c22 * c40 + 4 * C14 * c43 +
            3 * C13 * c46 * c57 + 3 * C1331 * c125 * c142 * c46 * c57 +
            3 * C1431 * c151 * c165 * c46 * c57 +
            3 * C1531 * c173 * c191 * c46 * c57 +
            3 * C0113 * c14 * c209 * c46 * c57 + C1211 * c114 * c222 * c85 +
            2 * C1221 * c114 * c211 * c40 * c85 +
            3 * C1231 * c114 * c46 * c57 * c85 + C1213 * c222 * c87 * c95) /
           Power(c46, 4);
  out(2) = (C21 * c295 + C2313 * c127 * c129 * c295 +
            C2312 * c134 * c136 * c295 + C2311 * c125 * c142 * c295 +
            C2413 * c153 * c155 * c295 + C2412 * c159 * c160 * c295 +
            C2411 * c151 * c165 * c295 + C2513 * c175 * c179 * c295 +
            C2512 * c183 * c185 * c295 + C2511 * c173 * c191 * c295 +
            C0211 * c14 * c209 * c295 + C0221 * c216 * c22 * c295 +
            C0231 * c17 * c226 * c295 + C1231 * c295 * c43 * c48 +
            C1221 * c295 * c57 * c61 + C1211 * c295 * c40 * c68 +
            2 * C22 * c283 * c85 + 2 * C2322 * c134 * c136 * c283 * c85 +
            2 * C2321 * c125 * c142 * c283 * c85 +
            2 * C2422 * c159 * c160 * c283 * c85 +
            2 * C2421 * c151 * c165 * c283 * c85 +
            2 * C2522 * c183 * c185 * c283 * c85 +
            2 * C2521 * c173 * c191 * c283 * c85 +
            2 * C0212 * c14 * c209 * c283 * c85 +
            2 * C0222 * c216 * c22 * c283 * c85 +
            2 * C1222 * c283 * c57 * c61 * c85 +
            2 * C1212 * c283 * c40 * c68 * c85 + 4 * C24 * c87 +
            3 * C23 * c101 * c89 + 3 * C2331 * c101 * c125 * c142 * c89 +
            3 * C2431 * c101 * c151 * c165 * c89 +
            3 * C2531 * c101 * c173 * c191 * c89 +
            3 * C0213 * c101 * c14 * c209 * c89 +
            3 * C1213 * c101 * c40 * c68 * c89) /
           Power(c89, 4);
  out(3) = (4 * C34 * c127 + 3 * C33 * c128 * c134 +
            3 * C3431 * c128 * c134 * c151 * c165 +
            3 * C3531 * c128 * c134 * c173 * c191 +
            3 * C0313 * c128 * c134 * c14 * c209 + 2 * C32 * c125 * c345 +
            2 * C2322 * c101 * c104 * c125 * c345 +
            2 * C3422 * c125 * c159 * c160 * c345 +
            2 * C3421 * c125 * c151 * c165 * c345 +
            2 * C3522 * c125 * c183 * c185 * c345 +
            2 * C3521 * c125 * c173 * c191 * c345 +
            2 * C0312 * c125 * c14 * c209 * c345 +
            2 * C0322 * c125 * c216 * c22 * c345 + C31 * c360 +
            C2321 * c101 * c104 * c360 + C3413 * c153 * c155 * c360 +
            C3412 * c159 * c160 * c360 + C3411 * c151 * c165 * c360 +
            C3513 * c175 * c179 * c360 + C3512 * c183 * c185 * c360 +
            C3511 * c173 * c191 * c360 + C0311 * c14 * c209 * c360 +
            C0321 * c216 * c22 * c360 + C0331 * c17 * c226 * c360 +
            C1331 * c360 * c43 * c48 + 2 * C1322 * c125 * c345 * c57 * c61 +
            C1321 * c360 * c57 * c61 + 3 * C1313 * c128 * c134 * c40 * c68 +
            2 * C1312 * c125 * c345 * c40 * c68 + C1311 * c360 * c40 * c68 +
            3 * C2313 * c114 * c128 * c134 * c85 +
            2 * C2312 * c114 * c125 * c345 * c85 + C2311 * c114 * c360 * c85 +
            C2331 * c360 * c87 * c95) /
           Power(c128, 4);
  out(4) = (4 * C44 * c153 + 3 * C43 * c154 * c159 +
            3 * C3413 * c125 * c142 * c154 * c159 +
            3 * C4531 * c154 * c159 * c173 * c191 +
            3 * C0413 * c14 * c154 * c159 * c209 + 2 * C42 * c151 * c400 +
            2 * C2422 * c101 * c104 * c151 * c400 +
            2 * C3422 * c134 * c136 * c151 * c400 +
            2 * C3412 * c125 * c142 * c151 * c400 +
            2 * C4522 * c151 * c183 * c185 * c400 +
            2 * C4521 * c151 * c173 * c191 * c400 +
            2 * C0412 * c14 * c151 * c209 * c400 +
            2 * C0422 * c151 * c216 * c22 * c400 + C41 * c419 +
            C2421 * c101 * c104 * c419 + C3431 * c127 * c129 * c419 +
            C3421 * c134 * c136 * c419 + C3411 * c125 * c142 * c419 +
            C4513 * c175 * c179 * c419 + C4512 * c183 * c185 * c419 +
            C4511 * c173 * c191 * c419 + C0411 * c14 * c209 * c419 +
            C0421 * c216 * c22 * c419 + C0431 * c17 * c226 * c419 +
            C1431 * c419 * c43 * c48 + 2 * C1422 * c151 * c400 * c57 * c61 +
            C1421 * c419 * c57 * c61 + 3 * C1413 * c154 * c159 * c40 * c68 +
            2 * C1412 * c151 * c40 * c400 * c68 + C1411 * c40 * c419 * c68 +
            3 * C2413 * c114 * c154 * c159 * c85 +
            2 * C2412 * c114 * c151 * c400 * c85 + C2411 * c114 * c419 * c85 +
            C2431 * c419 * c87 * c95) /
           Power(c154, 4);
  out(5) = (4 * C54 * c175 + 3 * C53 * c177 * c183 +
            3 * C3513 * c125 * c142 * c177 * c183 +
            3 * C4513 * c151 * c165 * c177 * c183 +
            3 * C0513 * c14 * c177 * c183 * c209 + 2 * C52 * c173 * c457 +
            2 * C2522 * c101 * c104 * c173 * c457 +
            2 * C3522 * c134 * c136 * c173 * c457 +
            2 * C3512 * c125 * c142 * c173 * c457 +
            2 * C4522 * c159 * c160 * c173 * c457 +
            2 * C4512 * c151 * c165 * c173 * c457 +
            2 * C0512 * c14 * c173 * c209 * c457 +
            2 * C0522 * c173 * c216 * c22 * c457 + C51 * c469 +
            C2521 * c101 * c104 * c469 + C3531 * c127 * c129 * c469 +
            C3521 * c134 * c136 * c469 + C3511 * c125 * c142 * c469 +
            C4531 * c153 * c155 * c469 + C4521 * c159 * c160 * c469 +
            C4511 * c151 * c165 * c469 + C0511 * c14 * c209 * c469 +
            C0521 * c216 * c22 * c469 + C0531 * c17 * c226 * c469 +
            C1531 * c43 * c469 * c48 + 2 * C1522 * c173 * c457 * c57 * c61 +
            C1521 * c469 * c57 * c61 + 3 * C1513 * c177 * c183 * c40 * c68 +
            2 * C1512 * c173 * c40 * c457 * c68 + C1511 * c40 * c469 * c68 +
            3 * C2513 * c114 * c177 * c183 * c85 +
            2 * C2512 * c114 * c173 * c457 * c85 + C2511 * c114 * c469 * c85 +
            C2531 * c469 * c87 * c95) /
           Power(c177, 4);

  return out;
}
gradhess_type FittedMaterial::gradhess_taylor_0(const strain_type &ek) {
  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real c11 = ek(0);
  Real c14 = -1 + c11;
  Real c7 = ekscale(0);
  Real c30 = Power(c7, 3);
  Real c40 = ek(1);
  Real c27 = Power(c7, 2);
  Real c46 = ekscale(1);
  Real c57 = Power(c40, 2);
  Real c61 = Power(c46, -2);
  Real c22 = Power(c14, 2);
  Real c68 = 1 / c46;
  Real c82 = ek(2);
  Real c85 = -1 + c82;
  Real c89 = ekscale(2);
  Real c101 = Power(c85, 2);
  Real c104 = Power(c89, -2);
  Real c114 = 1 / c89;
  Real c125 = ek(3);
  Real c128 = ekscale(3);
  Real c134 = Power(c125, 2);
  Real c136 = Power(c128, -2);
  Real c142 = 1 / c128;
  Real c151 = ek(4);
  Real c154 = ekscale(4);
  Real c159 = Power(c151, 2);
  Real c160 = Power(c154, -2);
  Real c165 = 1 / c154;
  Real c173 = ek(5);
  Real c177 = ekscale(5);
  Real c183 = Power(c173, 2);
  Real c185 = Power(c177, -2);
  Real c191 = 1 / c177;
  Real c43 = Power(c40, 3);
  Real c211 = Power(c46, 2);
  Real c209 = 1 / c7;
  Real c17 = Power(c14, 3);
  Real c222 = Power(c46, 3);
  Real c216 = Power(c7, -2);
  Real c87 = Power(c85, 3);
  Real c95 = Power(c89, -3);
  Real c127 = Power(c125, 3);
  Real c129 = Power(c128, -3);
  Real c153 = Power(c151, 3);
  Real c155 = Power(c154, -3);
  Real c175 = Power(c173, 3);
  Real c179 = Power(c177, -3);
  Real c283 = Power(c89, 2);
  Real c226 = Power(c7, -3);
  Real c295 = Power(c89, 3);
  Real c48 = Power(c46, -3);
  Real c345 = Power(c128, 2);
  Real c360 = Power(c128, 3);
  Real c400 = Power(c154, 2);
  Real c419 = Power(c154, 3);
  Real c457 = Power(c177, 2);
  Real c469 = Power(c177, 3);
  Real c9 = Power(c7, -4);
  Real c578 = 3 * C0113 * c27 * c57;
  Real c579 = 4 * C0122 * c14 * c40 * c7;
  Real c580 = 2 * C0112 * c27 * c40;
  Real c581 = 3 * C0131 * c22;
  Real c582 = 2 * C0121 * c14;
  Real c583 = C0111 * c7;
  Real c584 = c582 + c583;
  Real c585 = c584 * c7;
  Real c586 = c581 + c585;
  Real c587 = c46 * c586;
  Real c588 = c579 + c580 + c587;
  Real c589 = c46 * c588;
  Real c638 = c578 + c589;
  Real c639 = c226 * c48 * c638;
  Real c203 = Power(c46, -4);
  Real c640 = 3 * C0213 * c101 * c27;
  Real c641 = 4 * C0222 * c14 * c7 * c85;
  Real c642 = 2 * C0212 * c27 * c85;
  Real c643 = 3 * C0231 * c22;
  Real c644 = 2 * C0221 * c14;
  Real c645 = C0211 * c7;
  Real c646 = c644 + c645;
  Real c647 = c646 * c7;
  Real c648 = c643 + c647;
  Real c649 = c648 * c89;
  Real c650 = c641 + c642 + c649;
  Real c651 = c650 * c89;
  Real c652 = c640 + c651;
  Real c653 = c226 * c652 * c95;
  Real c794 = 3 * C1213 * c101 * c211;
  Real c795 = 4 * C1222 * c40 * c46 * c85;
  Real c796 = 2 * C1212 * c211 * c85;
  Real c797 = 3 * C1231 * c57;
  Real c798 = 2 * C1221 * c40 * c46;
  Real c799 = C1211 * c211;
  Real c800 = c797 + c798 + c799;
  Real c801 = c800 * c89;
  Real c802 = c795 + c796 + c801;
  Real c803 = c802 * c89;
  Real c804 = c794 + c803;
  Real c805 = c48 * c804 * c95;
  Real c274 = Power(c89, -4);
  Real c654 = 3 * C0313 * c134 * c27;
  Real c655 = 4 * C0322 * c125 * c14 * c7;
  Real c656 = 2 * C0312 * c125 * c27;
  Real c657 = 3 * C0331 * c22;
  Real c658 = 2 * C0321 * c14;
  Real c659 = C0311 * c7;
  Real c660 = c658 + c659;
  Real c661 = c660 * c7;
  Real c662 = c657 + c661;
  Real c663 = c128 * c662;
  Real c706 = c655 + c656 + c663;
  Real c707 = c128 * c706;
  Real c708 = c654 + c707;
  Real c709 = c129 * c226 * c708;
  Real c806 = 3 * C1313 * c134 * c211;
  Real c807 = 4 * C1322 * c125 * c40 * c46;
  Real c808 = 2 * C1312 * c125 * c211;
  Real c809 = 3 * C1331 * c57;
  Real c810 = 2 * C1321 * c40 * c46;
  Real c811 = C1311 * c211;
  Real c812 = c809 + c810 + c811;
  Real c813 = c128 * c812;
  Real c814 = c807 + c808 + c813;
  Real c815 = c128 * c814;
  Real c816 = c806 + c815;
  Real c817 = c129 * c48 * c816;
  Real c862 = 3 * C2313 * c134 * c283;
  Real c863 = 4 * C2322 * c125 * c85 * c89;
  Real c864 = 2 * C2312 * c125 * c283;
  Real c865 = 3 * C2331 * c101;
  Real c866 = 2 * C2321 * c85;
  Real c867 = C2311 * c89;
  Real c868 = c866 + c867;
  Real c869 = c868 * c89;
  Real c870 = c865 + c869;
  Real c871 = c128 * c870;
  Real c872 = c863 + c864 + c871;
  Real c873 = c128 * c872;
  Real c874 = c862 + c873;
  Real c875 = c129 * c874 * c95;
  Real c333 = Power(c128, -4);
  Real c710 = 3 * C0413 * c159 * c27;
  Real c711 = 4 * C0422 * c14 * c151 * c7;
  Real c712 = 2 * C0412 * c151 * c27;
  Real c713 = 3 * C0431 * c22;
  Real c714 = 2 * C0421 * c14;
  Real c715 = C0411 * c7;
  Real c716 = c714 + c715;
  Real c717 = c7 * c716;
  Real c718 = c713 + c717;
  Real c719 = c154 * c718;
  Real c720 = c711 + c712 + c719;
  Real c721 = c154 * c720;
  Real c722 = c710 + c721;
  Real c723 = c155 * c226 * c722;
  Real c818 = 3 * C1413 * c159 * c211;
  Real c819 = 4 * C1422 * c151 * c40 * c46;
  Real c820 = 2 * C1412 * c151 * c211;
  Real c821 = 3 * C1431 * c57;
  Real c822 = 2 * C1421 * c40 * c46;
  Real c823 = C1411 * c211;
  Real c824 = c821 + c822 + c823;
  Real c825 = c154 * c824;
  Real c826 = c819 + c820 + c825;
  Real c827 = c154 * c826;
  Real c828 = c818 + c827;
  Real c829 = c155 * c48 * c828;
  Real c876 = 3 * C2413 * c159 * c283;
  Real c877 = 4 * C2422 * c151 * c85 * c89;
  Real c878 = 2 * C2412 * c151 * c283;
  Real c879 = 3 * C2431 * c101;
  Real c880 = 2 * C2421 * c85;
  Real c881 = C2411 * c89;
  Real c882 = c880 + c881;
  Real c883 = c882 * c89;
  Real c884 = c879 + c883;
  Real c885 = c154 * c884;
  Real c886 = c877 + c878 + c885;
  Real c887 = c154 * c886;
  Real c888 = c876 + c887;
  Real c889 = c155 * c888 * c95;
  Real c924 = 3 * C3413 * c159 * c345;
  Real c925 = 4 * C3422 * c125 * c128 * c151;
  Real c926 = 2 * C3412 * c151 * c345;
  Real c927 = 3 * C3431 * c134;
  Real c928 = 2 * C3421 * c125 * c128;
  Real c929 = C3411 * c345;
  Real c930 = c927 + c928 + c929;
  Real c931 = c154 * c930;
  Real c932 = c925 + c926 + c931;
  Real c933 = c154 * c932;
  Real c934 = c924 + c933;
  Real c935 = c129 * c155 * c934;
  Real c392 = Power(c154, -4);
  Real c724 = 3 * C0513 * c183 * c27;
  Real c725 = 4 * C0522 * c14 * c173 * c7;
  Real c726 = 2 * C0512 * c173 * c27;
  Real c727 = 3 * C0531 * c22;
  Real c728 = 2 * C0521 * c14;
  Real c729 = C0511 * c7;
  Real c730 = c728 + c729;
  Real c731 = c7 * c730;
  Real c756 = c727 + c731;
  Real c757 = c177 * c756;
  Real c758 = c725 + c726 + c757;
  Real c759 = c177 * c758;
  Real c760 = c724 + c759;
  Real c761 = c179 * c226 * c760;
  Real c830 = 3 * C1513 * c183 * c211;
  Real c831 = 4 * C1522 * c173 * c40 * c46;
  Real c832 = 2 * C1512 * c173 * c211;
  Real c833 = 3 * C1531 * c57;
  Real c834 = 2 * C1521 * c40 * c46;
  Real c835 = C1511 * c211;
  Real c836 = c833 + c834 + c835;
  Real c837 = c177 * c836;
  Real c838 = c831 + c832 + c837;
  Real c839 = c177 * c838;
  Real c840 = c830 + c839;
  Real c841 = c179 * c48 * c840;
  Real c890 = 3 * C2513 * c183 * c283;
  Real c891 = 4 * C2522 * c173 * c85 * c89;
  Real c892 = 2 * C2512 * c173 * c283;
  Real c893 = 3 * C2531 * c101;
  Real c894 = 2 * C2521 * c85;
  Real c895 = C2511 * c89;
  Real c896 = c894 + c895;
  Real c897 = c89 * c896;
  Real c898 = c893 + c897;
  Real c899 = c177 * c898;
  Real c900 = c891 + c892 + c899;
  Real c901 = c177 * c900;
  Real c902 = c890 + c901;
  Real c903 = c179 * c902 * c95;
  Real c936 = 3 * C3513 * c183 * c345;
  Real c937 = 4 * C3522 * c125 * c128 * c173;
  Real c938 = 2 * C3512 * c173 * c345;
  Real c939 = 3 * C3531 * c134;
  Real c940 = 2 * C3521 * c125 * c128;
  Real c941 = C3511 * c345;
  Real c942 = c939 + c940 + c941;
  Real c943 = c177 * c942;
  Real c944 = c937 + c938 + c943;
  Real c945 = c177 * c944;
  Real c946 = c936 + c945;
  Real c947 = c129 * c179 * c946;
  Real c968 = 3 * C4513 * c183 * c400;
  Real c969 = 4 * C4522 * c151 * c154 * c173;
  Real c970 = 2 * C4512 * c173 * c400;
  Real c971 = 3 * C4531 * c159;
  Real c972 = 2 * C4521 * c151 * c154;
  Real c973 = C4511 * c400;
  Real c974 = c971 + c972 + c973;
  Real c975 = c177 * c974;
  Real c976 = c969 + c970 + c975;
  Real c977 = c177 * c976;
  Real c978 = c968 + c977;
  Real c979 = c155 * c179 * c978;
  Real c448 = Power(c177, -4);
  out1(0) =
      c9 *
      (4 * C04 * c17 + 2 * C02 * c14 * c27 +
       2 * C0222 * c101 * c104 * c14 * c27 +
       2 * C0322 * c134 * c136 * c14 * c27 +
       2 * C0321 * c125 * c14 * c142 * c27 +
       2 * C0422 * c14 * c159 * c160 * c27 +
       2 * C0421 * c14 * c151 * c165 * c27 +
       2 * C0522 * c14 * c183 * c185 * c27 +
       2 * C0521 * c14 * c173 * c191 * c27 + C01 * c30 +
       C0212 * c101 * c104 * c30 + C0313 * c127 * c129 * c30 +
       C0312 * c134 * c136 * c30 + C0311 * c125 * c142 * c30 +
       C0413 * c153 * c155 * c30 + C0412 * c159 * c160 * c30 +
       C0411 * c151 * c165 * c30 + C0513 * c175 * c179 * c30 +
       C0512 * c183 * c185 * c30 + C0511 * c173 * c191 * c30 +
       C0113 * c30 * c43 * c48 + 2 * C0122 * c14 * c27 * c57 * c61 +
       C0112 * c30 * c57 * c61 + 2 * C0121 * c14 * c27 * c40 * c68 +
       C0111 * c30 * c40 * c68 + 3 * C03 * c22 * c7 +
       3 * C0331 * c125 * c142 * c22 * c7 + 3 * C0431 * c151 * c165 * c22 * c7 +
       3 * C0531 * c173 * c191 * c22 * c7 + 3 * C0131 * c22 * c40 * c68 * c7 +
       2 * C0221 * c114 * c14 * c27 * c85 + C0211 * c114 * c30 * c85 +
       3 * C0231 * c114 * c22 * c7 * c85 + C0213 * c30 * c87 * c95);
  out1(1) =
      c203 * (C11 * c222 + C1212 * c101 * c104 * c222 +
              C1313 * c127 * c129 * c222 + C1312 * c134 * c136 * c222 +
              C1311 * c125 * c142 * c222 + C1413 * c153 * c155 * c222 +
              C1412 * c159 * c160 * c222 + C1411 * c151 * c165 * c222 +
              C1513 * c175 * c179 * c222 + C1512 * c183 * c185 * c222 +
              C1511 * c173 * c191 * c222 + C0111 * c14 * c209 * c222 +
              C0121 * c216 * c22 * c222 + C0131 * c17 * c222 * c226 +
              2 * C12 * c211 * c40 + 2 * C1222 * c101 * c104 * c211 * c40 +
              2 * C1322 * c134 * c136 * c211 * c40 +
              2 * C1321 * c125 * c142 * c211 * c40 +
              2 * C1422 * c159 * c160 * c211 * c40 +
              2 * C1421 * c151 * c165 * c211 * c40 +
              2 * C1522 * c183 * c185 * c211 * c40 +
              2 * C1521 * c173 * c191 * c211 * c40 +
              2 * C0112 * c14 * c209 * c211 * c40 +
              2 * C0122 * c211 * c216 * c22 * c40 + 4 * C14 * c43 +
              3 * C13 * c46 * c57 + 3 * C1331 * c125 * c142 * c46 * c57 +
              3 * C1431 * c151 * c165 * c46 * c57 +
              3 * C1531 * c173 * c191 * c46 * c57 +
              3 * C0113 * c14 * c209 * c46 * c57 + C1211 * c114 * c222 * c85 +
              2 * C1221 * c114 * c211 * c40 * c85 +
              3 * C1231 * c114 * c46 * c57 * c85 + C1213 * c222 * c87 * c95);
  out1(2) =
      c274 * (C21 * c295 + C2313 * c127 * c129 * c295 +
              C2312 * c134 * c136 * c295 + C2311 * c125 * c142 * c295 +
              C2413 * c153 * c155 * c295 + C2412 * c159 * c160 * c295 +
              C2411 * c151 * c165 * c295 + C2513 * c175 * c179 * c295 +
              C2512 * c183 * c185 * c295 + C2511 * c173 * c191 * c295 +
              C0211 * c14 * c209 * c295 + C0221 * c216 * c22 * c295 +
              C0231 * c17 * c226 * c295 + C1231 * c295 * c43 * c48 +
              C1221 * c295 * c57 * c61 + C1211 * c295 * c40 * c68 +
              2 * C22 * c283 * c85 + 2 * C2322 * c134 * c136 * c283 * c85 +
              2 * C2321 * c125 * c142 * c283 * c85 +
              2 * C2422 * c159 * c160 * c283 * c85 +
              2 * C2421 * c151 * c165 * c283 * c85 +
              2 * C2522 * c183 * c185 * c283 * c85 +
              2 * C2521 * c173 * c191 * c283 * c85 +
              2 * C0212 * c14 * c209 * c283 * c85 +
              2 * C0222 * c216 * c22 * c283 * c85 +
              2 * C1222 * c283 * c57 * c61 * c85 +
              2 * C1212 * c283 * c40 * c68 * c85 + 4 * C24 * c87 +
              3 * C23 * c101 * c89 + 3 * C2331 * c101 * c125 * c142 * c89 +
              3 * C2431 * c101 * c151 * c165 * c89 +
              3 * C2531 * c101 * c173 * c191 * c89 +
              3 * C0213 * c101 * c14 * c209 * c89 +
              3 * C1213 * c101 * c40 * c68 * c89);
  out1(3) =
      c333 * (4 * C34 * c127 + 3 * C33 * c128 * c134 +
              3 * C3431 * c128 * c134 * c151 * c165 +
              3 * C3531 * c128 * c134 * c173 * c191 +
              3 * C0313 * c128 * c134 * c14 * c209 + 2 * C32 * c125 * c345 +
              2 * C2322 * c101 * c104 * c125 * c345 +
              2 * C3422 * c125 * c159 * c160 * c345 +
              2 * C3421 * c125 * c151 * c165 * c345 +
              2 * C3522 * c125 * c183 * c185 * c345 +
              2 * C3521 * c125 * c173 * c191 * c345 +
              2 * C0312 * c125 * c14 * c209 * c345 +
              2 * C0322 * c125 * c216 * c22 * c345 + C31 * c360 +
              C2321 * c101 * c104 * c360 + C3413 * c153 * c155 * c360 +
              C3412 * c159 * c160 * c360 + C3411 * c151 * c165 * c360 +
              C3513 * c175 * c179 * c360 + C3512 * c183 * c185 * c360 +
              C3511 * c173 * c191 * c360 + C0311 * c14 * c209 * c360 +
              C0321 * c216 * c22 * c360 + C0331 * c17 * c226 * c360 +
              C1331 * c360 * c43 * c48 + 2 * C1322 * c125 * c345 * c57 * c61 +
              C1321 * c360 * c57 * c61 + 3 * C1313 * c128 * c134 * c40 * c68 +
              2 * C1312 * c125 * c345 * c40 * c68 + C1311 * c360 * c40 * c68 +
              3 * C2313 * c114 * c128 * c134 * c85 +
              2 * C2312 * c114 * c125 * c345 * c85 + C2311 * c114 * c360 * c85 +
              C2331 * c360 * c87 * c95);
  out1(4) =
      c392 * (4 * C44 * c153 + 3 * C43 * c154 * c159 +
              3 * C3413 * c125 * c142 * c154 * c159 +
              3 * C4531 * c154 * c159 * c173 * c191 +
              3 * C0413 * c14 * c154 * c159 * c209 + 2 * C42 * c151 * c400 +
              2 * C2422 * c101 * c104 * c151 * c400 +
              2 * C3422 * c134 * c136 * c151 * c400 +
              2 * C3412 * c125 * c142 * c151 * c400 +
              2 * C4522 * c151 * c183 * c185 * c400 +
              2 * C4521 * c151 * c173 * c191 * c400 +
              2 * C0412 * c14 * c151 * c209 * c400 +
              2 * C0422 * c151 * c216 * c22 * c400 + C41 * c419 +
              C2421 * c101 * c104 * c419 + C3431 * c127 * c129 * c419 +
              C3421 * c134 * c136 * c419 + C3411 * c125 * c142 * c419 +
              C4513 * c175 * c179 * c419 + C4512 * c183 * c185 * c419 +
              C4511 * c173 * c191 * c419 + C0411 * c14 * c209 * c419 +
              C0421 * c216 * c22 * c419 + C0431 * c17 * c226 * c419 +
              C1431 * c419 * c43 * c48 + 2 * C1422 * c151 * c400 * c57 * c61 +
              C1421 * c419 * c57 * c61 + 3 * C1413 * c154 * c159 * c40 * c68 +
              2 * C1412 * c151 * c40 * c400 * c68 + C1411 * c40 * c419 * c68 +
              3 * C2413 * c114 * c154 * c159 * c85 +
              2 * C2412 * c114 * c151 * c400 * c85 + C2411 * c114 * c419 * c85 +
              C2431 * c419 * c87 * c95);
  out1(5) =
      c448 * (4 * C54 * c175 + 3 * C53 * c177 * c183 +
              3 * C3513 * c125 * c142 * c177 * c183 +
              3 * C4513 * c151 * c165 * c177 * c183 +
              3 * C0513 * c14 * c177 * c183 * c209 + 2 * C52 * c173 * c457 +
              2 * C2522 * c101 * c104 * c173 * c457 +
              2 * C3522 * c134 * c136 * c173 * c457 +
              2 * C3512 * c125 * c142 * c173 * c457 +
              2 * C4522 * c159 * c160 * c173 * c457 +
              2 * C4512 * c151 * c165 * c173 * c457 +
              2 * C0512 * c14 * c173 * c209 * c457 +
              2 * C0522 * c173 * c216 * c22 * c457 + C51 * c469 +
              C2521 * c101 * c104 * c469 + C3531 * c127 * c129 * c469 +
              C3521 * c134 * c136 * c469 + C3511 * c125 * c142 * c469 +
              C4531 * c153 * c155 * c469 + C4521 * c159 * c160 * c469 +
              C4511 * c151 * c165 * c469 + C0511 * c14 * c209 * c469 +
              C0521 * c216 * c22 * c469 + C0531 * c17 * c226 * c469 +
              C1531 * c43 * c469 * c48 + 2 * C1522 * c173 * c457 * c57 * c61 +
              C1521 * c469 * c57 * c61 + 3 * C1513 * c177 * c183 * c40 * c68 +
              2 * C1512 * c173 * c40 * c457 * c68 + C1511 * c40 * c469 * c68 +
              3 * C2513 * c114 * c177 * c183 * c85 +
              2 * C2512 * c114 * c173 * c457 * c85 + C2511 * c114 * c469 * c85 +
              C2531 * c469 * c87 * c95);
  out2(0, 0) =
      2 *
      (6 * C04 * c22 + C02 * c27 + C0222 * c101 * c104 * c27 +
       C0322 * c134 * c136 * c27 + C0321 * c125 * c142 * c27 +
       C0422 * c159 * c160 * c27 + C0421 * c151 * c165 * c27 +
       C0522 * c183 * c185 * c27 + C0521 * c173 * c191 * c27 +
       C0122 * c27 * c57 * c61 + C0121 * c27 * c40 * c68 + 3 * C03 * c14 * c7 +
       3 * C0331 * c125 * c14 * c142 * c7 + 3 * C0431 * c14 * c151 * c165 * c7 +
       3 * C0531 * c14 * c173 * c191 * c7 + 3 * C0131 * c14 * c40 * c68 * c7 +
       C0221 * c114 * c27 * c85 + 3 * C0231 * c114 * c14 * c7 * c85) *
      c9;
  out2(0, 1) = c639;
  out2(0, 2) = c653;
  out2(0, 3) = c709;
  out2(0, 4) = c723;
  out2(0, 5) = c761;
  out2(1, 0) = c639;
  out2(1, 1) = 2 * c203 *
               (C12 * c211 + C1222 * c101 * c104 * c211 +
                C1322 * c134 * c136 * c211 + C1321 * c125 * c142 * c211 +
                C1422 * c159 * c160 * c211 + C1421 * c151 * c165 * c211 +
                C1522 * c183 * c185 * c211 + C1521 * c173 * c191 * c211 +
                C0112 * c14 * c209 * c211 + C0122 * c211 * c216 * c22 +
                3 * C13 * c40 * c46 + 3 * C1331 * c125 * c142 * c40 * c46 +
                3 * C1431 * c151 * c165 * c40 * c46 +
                3 * C1531 * c173 * c191 * c40 * c46 +
                3 * C0113 * c14 * c209 * c40 * c46 + 6 * C14 * c57 +
                C1221 * c114 * c211 * c85 + 3 * C1231 * c114 * c40 * c46 * c85);
  out2(1, 2) = c805;
  out2(1, 3) = c817;
  out2(1, 4) = c829;
  out2(1, 5) = c841;
  out2(2, 0) = c653;
  out2(2, 1) = c805;
  out2(2, 2) =
      2 * c274 *
      (6 * C24 * c101 + C22 * c283 + C2322 * c134 * c136 * c283 +
       C2321 * c125 * c142 * c283 + C2422 * c159 * c160 * c283 +
       C2421 * c151 * c165 * c283 + C2522 * c183 * c185 * c283 +
       C2521 * c173 * c191 * c283 + C0212 * c14 * c209 * c283 +
       C0222 * c216 * c22 * c283 + C1222 * c283 * c57 * c61 +
       C1212 * c283 * c40 * c68 + 3 * C23 * c85 * c89 +
       3 * C2331 * c125 * c142 * c85 * c89 +
       3 * C2431 * c151 * c165 * c85 * c89 +
       3 * C2531 * c173 * c191 * c85 * c89 +
       3 * C0213 * c14 * c209 * c85 * c89 + 3 * C1213 * c40 * c68 * c85 * c89);
  out2(2, 3) = c875;
  out2(2, 4) = c889;
  out2(2, 5) = c903;
  out2(3, 0) = c709;
  out2(3, 1) = c817;
  out2(3, 2) = c875;
  out2(3, 3) =
      2 * c333 *
      (3 * C33 * c125 * c128 + 6 * C34 * c134 +
       3 * C3431 * c125 * c128 * c151 * c165 +
       3 * C3531 * c125 * c128 * c173 * c191 +
       3 * C0313 * c125 * c128 * c14 * c209 + C32 * c345 +
       C2322 * c101 * c104 * c345 + C3422 * c159 * c160 * c345 +
       C3421 * c151 * c165 * c345 + C3522 * c183 * c185 * c345 +
       C3521 * c173 * c191 * c345 + C0312 * c14 * c209 * c345 +
       C0322 * c216 * c22 * c345 + C1322 * c345 * c57 * c61 +
       3 * C1313 * c125 * c128 * c40 * c68 + C1312 * c345 * c40 * c68 +
       3 * C2313 * c114 * c125 * c128 * c85 + C2312 * c114 * c345 * c85);
  out2(3, 4) = c935;
  out2(3, 5) = c947;
  out2(4, 0) = c723;
  out2(4, 1) = c829;
  out2(4, 2) = c889;
  out2(4, 3) = c935;
  out2(4, 4) =
      2 * c392 *
      (3 * C43 * c151 * c154 + 3 * C3413 * c125 * c142 * c151 * c154 +
       6 * C44 * c159 + 3 * C4531 * c151 * c154 * c173 * c191 +
       3 * C0413 * c14 * c151 * c154 * c209 + C42 * c400 +
       C2422 * c101 * c104 * c400 + C3422 * c134 * c136 * c400 +
       C3412 * c125 * c142 * c400 + C4522 * c183 * c185 * c400 +
       C4521 * c173 * c191 * c400 + C0412 * c14 * c209 * c400 +
       C0422 * c216 * c22 * c400 + C1422 * c400 * c57 * c61 +
       3 * C1413 * c151 * c154 * c40 * c68 + C1412 * c40 * c400 * c68 +
       3 * C2413 * c114 * c151 * c154 * c85 + C2412 * c114 * c400 * c85);
  out2(4, 5) = c979;
  out2(5, 0) = c761;
  out2(5, 1) = c841;
  out2(5, 2) = c903;
  out2(5, 3) = c947;
  out2(5, 4) = c979;
  out2(5, 5) =
      2 * c448 *
      (3 * C53 * c173 * c177 + 3 * C3513 * c125 * c142 * c173 * c177 +
       3 * C4513 * c151 * c165 * c173 * c177 + 6 * C54 * c183 +
       3 * C0513 * c14 * c173 * c177 * c209 + C52 * c457 +
       C2522 * c101 * c104 * c457 + C3522 * c134 * c136 * c457 +
       C3512 * c125 * c142 * c457 + C4522 * c159 * c160 * c457 +
       C4512 * c151 * c165 * c457 + C0512 * c14 * c209 * c457 +
       C0522 * c216 * c22 * c457 + C1522 * c457 * c57 * c61 +
       3 * C1513 * c173 * c177 * c40 * c68 + C1512 * c40 * c457 * c68 +
       3 * C2513 * c114 * c173 * c177 * c85 + C2512 * c114 * c457 * c85);

  return std::make_pair(hess, grad);
}

#define BARRIERPOWER 4 // should be even integer!
value_type FittedMaterial::psi_barrier(const strain_type &ek,
                                       const strain_type &ekclamped) {
  Real copt106 = ekscale(0);
  Real copt107 = 1 / copt106;
  Real copt120 = ekscale(1);
  Real copt121 = 1 / copt120;
  Real copt132 = ekscale(2);
  Real copt133 = 1 / copt132;
  Real copt148 = ekscale(3);
  Real copt149 = 1 / copt148;
  Real copt158 = ekscale(4);
  Real copt159 = 1 / copt158;
  Real copt166 = ekscale(5);
  Real copt167 = 1 / copt166;
  return bscale *
         (Power(bspeed * copt107 * (ek(0) - ekclamped(0)), BARRIERPOWER) +
          Power(bspeed * copt121 * (ek(1) - ekclamped(1)), BARRIERPOWER) +
          Power(bspeed * copt133 * (ek(2) - ekclamped(2)), BARRIERPOWER) +
          Power(bspeed * copt149 * (ek(3) - ekclamped(3)), BARRIERPOWER) +
          Power(bspeed * copt159 * (ek(4) - ekclamped(4)), BARRIERPOWER) +
          Power(bspeed * copt167 * (ek(5) - ekclamped(5)), BARRIERPOWER));
}
grad_type FittedMaterial::grad_barrier(const strain_type &ek,
                                       const strain_type &ekclamped) {
  Vec6 out(0);
  Real copt106 = ekscale(0);
  Real copt107 = 1 / copt106;
  Real copt122 = ekscale(1);
  Real copt123 = 1 / copt122;
  int copt109 = -1 + BARRIERPOWER;
  Real copt139 = ekscale(2);
  Real copt140 = 1 / copt139;
  Real copt152 = ekscale(3);
  Real copt153 = 1 / copt152;
  Real copt163 = ekscale(4);
  Real copt164 = 1 / copt163;
  Real copt172 = ekscale(5);
  Real copt173 = 1 / copt172;
  out(0) = BARRIERPOWER * bscale * bspeed * copt107 *
           Power(bspeed * copt107 * (ek(0) - ekclamped(0)), copt109);
  out(1) = BARRIERPOWER * bscale * bspeed * copt123 *
           Power(bspeed * copt123 * (ek(1) - ekclamped(1)), copt109);
  out(2) = BARRIERPOWER * bscale * bspeed * copt140 *
           Power(bspeed * copt140 * (ek(2) - ekclamped(2)), copt109);
  out(3) = BARRIERPOWER * bscale * bspeed * copt153 *
           Power(bspeed * copt153 * (ek(3) - ekclamped(3)), copt109);
  out(4) = BARRIERPOWER * bscale * bspeed * copt164 *
           Power(bspeed * copt164 * (ek(4) - ekclamped(4)), copt109);
  out(5) = BARRIERPOWER * bscale * bspeed * copt173 *
           Power(bspeed * copt173 * (ek(5) - ekclamped(5)), copt109);

  return out;
}
gradhess_type FittedMaterial::gradhess_barrier(const strain_type &ek,
                                               const strain_type &ekclamped) {
  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };
  Real copt26 = ekscale(0);
  Real copt29 = 1 / copt26;
  Real copt38 = ekscale(1);
  Real copt39 = 1 / copt38;
  int copt31 = -1 + BARRIERPOWER;
  Real copt47 = ekscale(2);
  Real copt48 = 1 / copt47;
  Real copt56 = ekscale(3);
  Real copt57 = 1 / copt56;
  Real copt65 = ekscale(4);
  Real copt66 = 1 / copt65;
  Real copt74 = ekscale(5);
  Real copt75 = 1 / copt74;
  Real copt79 = Power(bspeed, 2);
  Real copt1 = ek(0);
  Real copt5 = ekclamped(0);
  Real copt18 = -copt5;
  Real copt23 = copt1 + copt18;
  Real copt30 = bspeed * copt23 * copt29;
  Real copt82 = Power(copt26, 2);
  Real copt83 = 1 / copt82;
  Real copt34 = ek(1);
  Real copt35 = ekclamped(1);
  Real copt36 = -copt35;
  Real copt37 = copt34 + copt36;
  Real copt40 = bspeed * copt37 * copt39;
  int copt80 = -2 + BARRIERPOWER;
  Real copt86 = Power(copt38, 2);
  Real copt87 = 1 / copt86;
  Real copt43 = ek(2);
  Real copt44 = ekclamped(2);
  Real copt45 = -copt44;
  Real copt46 = copt43 + copt45;
  Real copt49 = bspeed * copt46 * copt48;
  Real copt90 = Power(copt47, 2);
  Real copt91 = 1 / copt90;
  Real copt52 = ek(3);
  Real copt53 = ekclamped(3);
  Real copt54 = -copt53;
  Real copt55 = copt52 + copt54;
  Real copt58 = bspeed * copt55 * copt57;
  Real copt94 = Power(copt56, 2);
  Real copt95 = 1 / copt94;
  Real copt61 = ek(4);
  Real copt62 = ekclamped(4);
  Real copt63 = -copt62;
  Real copt64 = copt61 + copt63;
  Real copt67 = bspeed * copt64 * copt66;
  Real copt98 = Power(copt65, 2);
  Real copt99 = 1 / copt98;
  Real copt70 = ek(5);
  Real copt71 = ekclamped(5);
  Real copt72 = -copt71;
  Real copt73 = copt70 + copt72;
  Real copt76 = bspeed * copt73 * copt75;
  Real copt102 = Power(copt74, 2);
  Real copt103 = 1 / copt102;
  out1(0) = BARRIERPOWER * bscale * bspeed * copt29 * Power(copt30, copt31);
  out1(1) = BARRIERPOWER * bscale * bspeed * copt39 * Power(copt40, copt31);
  out1(2) = BARRIERPOWER * bscale * bspeed * copt48 * Power(copt49, copt31);
  out1(3) = BARRIERPOWER * bscale * bspeed * copt57 * Power(copt58, copt31);
  out1(4) = BARRIERPOWER * bscale * bspeed * copt66 * Power(copt67, copt31);
  out1(5) = BARRIERPOWER * bscale * bspeed * copt75 * Power(copt76, copt31);
  out2(0, 0) =
      BARRIERPOWER * bscale * Power(copt30, copt80) * copt31 * copt79 * copt83;
  out2(0, 1) = 0;
  out2(0, 2) = 0;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) =
      BARRIERPOWER * bscale * copt31 * Power(copt40, copt80) * copt79 * copt87;
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) =
      BARRIERPOWER * bscale * copt31 * Power(copt49, copt80) * copt79 * copt91;
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) =
      BARRIERPOWER * bscale * copt31 * Power(copt58, copt80) * copt79 * copt95;
  out2(3, 4) = 0;
  out2(3, 5) = 0;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = 0;
  out2(4, 4) =
      BARRIERPOWER * bscale * copt31 * Power(copt67, copt80) * copt79 * copt99;
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) =
      BARRIERPOWER * bscale * copt103 * copt31 * Power(copt76, copt80) * copt79;
  return std::make_pair(hess, grad);
}
