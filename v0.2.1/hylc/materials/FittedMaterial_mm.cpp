#include "FittedMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

typedef Vec6 strain_type;
typedef double value_type;
typedef Vec6 grad_type;
typedef std::pair<Mat6x6, Vec6> gradhess_type;

// 0th derivative of psi, grad, tuple<grad,hess>
value_type FittedMaterial::psi_taylor_0(const strain_type &ek) {
  Real c1 = ek(0);
  Real c2 = -1 + c1;
  Real c10 = ekscale(0);
  Real c23 = ek(1);
  Real c25 = ekscale(1);
  Real c28 = Power(c23, 3);
  Real c21 = 1 / c10;
  Real c29 = Power(c25, -3);
  Real c17 = Power(c2, 2);
  Real c32 = Power(c23, 2);
  Real c18 = Power(c10, -2);
  Real c33 = Power(c25, -2);
  Real c14 = Power(c2, 3);
  Real c15 = Power(c10, -3);
  Real c37 = 1 / c25;
  Real c43 = ek(2);
  Real c44 = -1 + c43;
  Real c46 = ekscale(2);
  Real c49 = Power(c44, 3);
  Real c50 = Power(c46, -3);
  Real c53 = Power(c44, 2);
  Real c54 = Power(c46, -2);
  Real c58 = 1 / c46;
  Real c63 = ek(3);
  Real c65 = ekscale(3);
  Real c68 = Power(c63, 3);
  Real c69 = Power(c65, -3);
  Real c73 = Power(c63, 2);
  Real c74 = Power(c65, -2);
  Real c80 = 1 / c65;
  Real c88 = ek(4);
  Real c90 = ekscale(4);
  Real c93 = Power(c88, 3);
  Real c94 = Power(c90, -3);
  Real c99 = Power(c88, 2);
  Real c100 = Power(c90, -2);
  Real c108 = 1 / c90;
  Real c119 = ek(5);
  Real c121 = ekscale(5);
  Real c124 = Power(c119, 3);
  Real c125 = Power(c121, -3);
  Real c130 = Power(c119, 2);
  Real c131 = Power(c121, -2);
  Real c139 = 1 / c121;
  return C0 + (C54 * Power(c119, 4)) / Power(c121, 4) + C53 * c124 * c125 +
         C52 * c130 * c131 + C51 * c119 * c139 + C03 * c14 * c15 +
         C0531 * c119 * c139 * c14 * c15 + C02 * c17 * c18 +
         C0522 * c130 * c131 * c17 * c18 + C0521 * c119 * c139 * c17 * c18 +
         (C04 * Power(c2, 4)) / Power(c10, 4) + C01 * c2 * c21 +
         C0513 * c124 * c125 * c2 * c21 + C0512 * c130 * c131 * c2 * c21 +
         C0511 * c119 * c139 * c2 * c21 +
         (C14 * Power(c23, 4)) / Power(c25, 4) + C13 * c28 * c29 +
         C0113 * c2 * c21 * c28 * c29 + C12 * c32 * c33 +
         C0122 * c17 * c18 * c32 * c33 + C0112 * c2 * c21 * c32 * c33 +
         C11 * c23 * c37 + C0131 * c14 * c15 * c23 * c37 +
         C0121 * c17 * c18 * c23 * c37 + C0111 * c2 * c21 * c23 * c37 +
         (C24 * Power(c44, 4)) / Power(c46, 4) + C23 * c49 * c50 +
         C2531 * c119 * c139 * c49 * c50 + C0213 * c2 * c21 * c49 * c50 +
         C22 * c53 * c54 + C2522 * c130 * c131 * c53 * c54 +
         C2521 * c119 * c139 * c53 * c54 + C0222 * c17 * c18 * c53 * c54 +
         C0212 * c2 * c21 * c53 * c54 + C21 * c44 * c58 +
         C2513 * c124 * c125 * c44 * c58 + C2512 * c130 * c131 * c44 * c58 +
         C2511 * c119 * c139 * c44 * c58 + C0231 * c14 * c15 * c44 * c58 +
         C0221 * c17 * c18 * c44 * c58 + C0211 * c2 * c21 * c44 * c58 +
         (C34 * Power(c63, 4)) / Power(c65, 4) + C33 * c68 * c69 +
         C3531 * c119 * c139 * c68 * c69 + C0313 * c2 * c21 * c68 * c69 +
         C2313 * c44 * c58 * c68 * c69 + C32 * c73 * c74 +
         C3522 * c130 * c131 * c73 * c74 + C3521 * c119 * c139 * c73 * c74 +
         C0322 * c17 * c18 * c73 * c74 + C0312 * c2 * c21 * c73 * c74 +
         C2322 * c53 * c54 * c73 * c74 + C2312 * c44 * c58 * c73 * c74 +
         C31 * c63 * c80 + C3513 * c124 * c125 * c63 * c80 +
         C3512 * c130 * c131 * c63 * c80 + C3511 * c119 * c139 * c63 * c80 +
         C0331 * c14 * c15 * c63 * c80 + C0321 * c17 * c18 * c63 * c80 +
         C0311 * c2 * c21 * c63 * c80 + C2331 * c49 * c50 * c63 * c80 +
         C2321 * c53 * c54 * c63 * c80 + C2311 * c44 * c58 * c63 * c80 +
         C41 * c108 * c88 + C0431 * c108 * c14 * c15 * c88 +
         C0421 * c108 * c17 * c18 * c88 + C0411 * c108 * c2 * c21 * c88 +
         C2431 * c108 * c49 * c50 * c88 + C2421 * c108 * c53 * c54 * c88 +
         C2411 * c108 * c44 * c58 * c88 + C3431 * c108 * c68 * c69 * c88 +
         C3421 * c108 * c73 * c74 * c88 + C3411 * c108 * c63 * c80 * c88 +
         (C44 * Power(c88, 4)) / Power(c90, 4) + C43 * c93 * c94 +
         C0413 * c2 * c21 * c93 * c94 + C2413 * c44 * c58 * c93 * c94 +
         C3413 * c63 * c80 * c93 * c94 + C42 * c100 * c99 +
         C0422 * c100 * c17 * c18 * c99 + C0412 * c100 * c2 * c21 * c99 +
         C2422 * c100 * c53 * c54 * c99 + C2412 * c100 * c44 * c58 * c99 +
         C3422 * c100 * c73 * c74 * c99 + C3412 * c100 * c63 * c80 * c99;
}
grad_type FittedMaterial::grad_taylor_0(const strain_type &ek) {
  Vec6 out(0);

  Real c9 = ek(0);
  Real c10 = -1 + c9;
  Real c1 = ekscale(0);
  Real c18 = Power(c1, 3);
  Real c21 = ek(1);
  Real c16 = Power(c1, 2);
  Real c23 = ekscale(1);
  Real c26 = Power(c21, 2);
  Real c27 = Power(c23, -2);
  Real c14 = Power(c10, 2);
  Real c30 = 1 / c23;
  Real c34 = ek(2);
  Real c35 = -1 + c34;
  Real c37 = ekscale(2);
  Real c40 = Power(c35, 2);
  Real c42 = Power(c37, -2);
  Real c45 = 1 / c37;
  Real c49 = ek(3);
  Real c51 = ekscale(3);
  Real c54 = Power(c49, 2);
  Real c55 = Power(c51, -2);
  Real c58 = 1 / c51;
  Real c62 = ek(4);
  Real c64 = ekscale(4);
  Real c67 = Power(c62, 2);
  Real c68 = Power(c64, -2);
  Real c71 = 1 / c64;
  Real c75 = ek(5);
  Real c77 = ekscale(5);
  Real c80 = Power(c75, 2);
  Real c81 = Power(c77, -2);
  Real c84 = 1 / c77;
  Real c22 = Power(c21, 3);
  Real c95 = Power(c23, 2);
  Real c93 = 1 / c1;
  Real c11 = Power(c10, 3);
  Real c100 = Power(c23, 3);
  Real c97 = Power(c1, -2);
  Real c36 = Power(c35, 3);
  Real c112 = Power(c37, 2);
  Real c102 = Power(c1, -3);
  Real c116 = Power(c37, 3);
  Real c50 = Power(c49, 3);
  Real c52 = Power(c51, -3);
  Real c63 = Power(c62, 3);
  Real c65 = Power(c64, -3);
  Real c76 = Power(c75, 3);
  Real c78 = Power(c77, -3);
  Real c146 = Power(c51, 2);
  Real c152 = Power(c51, 3);
  Real c38 = Power(c37, -3);
  Real c180 = Power(c64, 2);
  Real c188 = Power(c64, 3);
  Real c207 = Power(c77, 2);
  Real c215 = Power(c77, 3);
  out(0) =
      (4 * C04 * c11 + 3 * C03 * c1 * c14 + 2 * C02 * c10 * c16 + C01 * c18 +
       (C0113 * c18 * c22) / Power(c23, 3) + 2 * C0122 * c10 * c16 * c26 * c27 +
       C0112 * c18 * c26 * c27 + 3 * C0131 * c1 * c14 * c21 * c30 +
       2 * C0121 * c10 * c16 * c21 * c30 + C0111 * c18 * c21 * c30 +
       C0213 * c18 * c36 * c38 + 2 * C0222 * c10 * c16 * c40 * c42 +
       C0212 * c18 * c40 * c42 + 3 * C0231 * c1 * c14 * c35 * c45 +
       2 * C0221 * c10 * c16 * c35 * c45 + C0211 * c18 * c35 * c45 +
       C0313 * c18 * c50 * c52 + 2 * C0322 * c10 * c16 * c54 * c55 +
       C0312 * c18 * c54 * c55 + 3 * C0331 * c1 * c14 * c49 * c58 +
       2 * C0321 * c10 * c16 * c49 * c58 + C0311 * c18 * c49 * c58 +
       C0413 * c18 * c63 * c65 + 2 * C0422 * c10 * c16 * c67 * c68 +
       C0412 * c18 * c67 * c68 + 3 * C0431 * c1 * c14 * c62 * c71 +
       2 * C0421 * c10 * c16 * c62 * c71 + C0411 * c18 * c62 * c71 +
       C0513 * c18 * c76 * c78 + 2 * C0522 * c10 * c16 * c80 * c81 +
       C0512 * c18 * c80 * c81 + 3 * C0531 * c1 * c14 * c75 * c84 +
       2 * C0521 * c10 * c16 * c75 * c84 + C0511 * c18 * c75 * c84) /
      Power(c1, 4);
  out(1) = (C11 * c100 + C0131 * c100 * c102 * c11 + 4 * C14 * c22 +
            3 * C13 * c23 * c26 + C0111 * c10 * c100 * c93 +
            3 * C0113 * c10 * c23 * c26 * c93 + 2 * C12 * c21 * c95 +
            2 * C0112 * c10 * c21 * c93 * c95 + C0121 * c100 * c14 * c97 +
            2 * C0122 * c14 * c21 * c95 * c97) /
           Power(c23, 4);
  out(2) =
      (C21 * c116 + C0231 * c102 * c11 * c116 + 2 * C22 * c112 * c35 +
       4 * C24 * c36 + 3 * C23 * c37 * c40 + C2313 * c116 * c50 * c52 +
       C2312 * c116 * c54 * c55 + 2 * C2322 * c112 * c35 * c54 * c55 +
       C2311 * c116 * c49 * c58 + 2 * C2321 * c112 * c35 * c49 * c58 +
       3 * C2331 * c37 * c40 * c49 * c58 + C2413 * c116 * c63 * c65 +
       C2412 * c116 * c67 * c68 + 2 * C2422 * c112 * c35 * c67 * c68 +
       C2411 * c116 * c62 * c71 + 2 * C2421 * c112 * c35 * c62 * c71 +
       3 * C2431 * c37 * c40 * c62 * c71 + C2513 * c116 * c76 * c78 +
       C2512 * c116 * c80 * c81 + 2 * C2522 * c112 * c35 * c80 * c81 +
       C2511 * c116 * c75 * c84 + 2 * C2521 * c112 * c35 * c75 * c84 +
       3 * C2531 * c37 * c40 * c75 * c84 + C0211 * c10 * c116 * c93 +
       2 * C0212 * c10 * c112 * c35 * c93 + 3 * C0213 * c10 * c37 * c40 * c93 +
       C0221 * c116 * c14 * c97 + 2 * C0222 * c112 * c14 * c35 * c97) /
      Power(c37, 4);
  out(3) =
      (C31 * c152 + C0331 * c102 * c11 * c152 + C2331 * c152 * c36 * c38 +
       C2321 * c152 * c40 * c42 + C2311 * c152 * c35 * c45 +
       2 * C32 * c146 * c49 + 2 * C2322 * c146 * c40 * c42 * c49 +
       2 * C2312 * c146 * c35 * c45 * c49 + 4 * C34 * c50 +
       3 * C33 * c51 * c54 + 3 * C2313 * c35 * c45 * c51 * c54 +
       C3413 * c152 * c63 * c65 + C3412 * c152 * c67 * c68 +
       2 * C3422 * c146 * c49 * c67 * c68 + C3411 * c152 * c62 * c71 +
       2 * C3421 * c146 * c49 * c62 * c71 + 3 * C3431 * c51 * c54 * c62 * c71 +
       C3513 * c152 * c76 * c78 + C3512 * c152 * c80 * c81 +
       2 * C3522 * c146 * c49 * c80 * c81 + C3511 * c152 * c75 * c84 +
       2 * C3521 * c146 * c49 * c75 * c84 + 3 * C3531 * c51 * c54 * c75 * c84 +
       C0311 * c10 * c152 * c93 + 2 * C0312 * c10 * c146 * c49 * c93 +
       3 * C0313 * c10 * c51 * c54 * c93 + C0321 * c14 * c152 * c97 +
       2 * C0322 * c14 * c146 * c49 * c97) /
      Power(c51, 4);
  out(4) =
      (C41 * c188 + C0431 * c102 * c11 * c188 + C2431 * c188 * c36 * c38 +
       C2421 * c188 * c40 * c42 + C2411 * c188 * c35 * c45 +
       C3431 * c188 * c50 * c52 + C3421 * c188 * c54 * c55 +
       C3411 * c188 * c49 * c58 + 2 * C42 * c180 * c62 +
       2 * C2422 * c180 * c40 * c42 * c62 + 2 * C2412 * c180 * c35 * c45 * c62 +
       2 * C3422 * c180 * c54 * c55 * c62 + 2 * C3412 * c180 * c49 * c58 * c62 +
       4 * C44 * c63 + 3 * C43 * c64 * c67 + 3 * C2413 * c35 * c45 * c64 * c67 +
       3 * C3413 * c49 * c58 * c64 * c67 + C0411 * c10 * c188 * c93 +
       2 * C0412 * c10 * c180 * c62 * c93 + 3 * C0413 * c10 * c64 * c67 * c93 +
       C0421 * c14 * c188 * c97 + 2 * C0422 * c14 * c180 * c62 * c97) /
      Power(c64, 4);
  out(5) =
      (C51 * c215 + C0531 * c102 * c11 * c215 + C2531 * c215 * c36 * c38 +
       C2521 * c215 * c40 * c42 + C2511 * c215 * c35 * c45 +
       C3531 * c215 * c50 * c52 + C3521 * c215 * c54 * c55 +
       C3511 * c215 * c49 * c58 + 2 * C52 * c207 * c75 +
       2 * C2522 * c207 * c40 * c42 * c75 + 2 * C2512 * c207 * c35 * c45 * c75 +
       2 * C3522 * c207 * c54 * c55 * c75 + 2 * C3512 * c207 * c49 * c58 * c75 +
       4 * C54 * c76 + 3 * C53 * c77 * c80 + 3 * C2513 * c35 * c45 * c77 * c80 +
       3 * C3513 * c49 * c58 * c77 * c80 + C0511 * c10 * c215 * c93 +
       2 * C0512 * c10 * c207 * c75 * c93 + 3 * C0513 * c10 * c77 * c80 * c93 +
       C0521 * c14 * c215 * c97 + 2 * C0522 * c14 * c207 * c75 * c97) /
      Power(c77, 4);

  return out;
}
gradhess_type FittedMaterial::gradhess_taylor_0(const strain_type &ek) {
  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real c9 = ek(0);
  Real c10 = -1 + c9;
  Real c1 = ekscale(0);
  Real c18 = Power(c1, 3);
  Real c21 = ek(1);
  Real c16 = Power(c1, 2);
  Real c23 = ekscale(1);
  Real c26 = Power(c21, 2);
  Real c27 = Power(c23, -2);
  Real c14 = Power(c10, 2);
  Real c30 = 1 / c23;
  Real c34 = ek(2);
  Real c35 = -1 + c34;
  Real c37 = ekscale(2);
  Real c40 = Power(c35, 2);
  Real c42 = Power(c37, -2);
  Real c45 = 1 / c37;
  Real c49 = ek(3);
  Real c51 = ekscale(3);
  Real c54 = Power(c49, 2);
  Real c55 = Power(c51, -2);
  Real c58 = 1 / c51;
  Real c62 = ek(4);
  Real c64 = ekscale(4);
  Real c67 = Power(c62, 2);
  Real c68 = Power(c64, -2);
  Real c71 = 1 / c64;
  Real c75 = ek(5);
  Real c77 = ekscale(5);
  Real c80 = Power(c75, 2);
  Real c81 = Power(c77, -2);
  Real c84 = 1 / c77;
  Real c22 = Power(c21, 3);
  Real c95 = Power(c23, 2);
  Real c93 = 1 / c1;
  Real c11 = Power(c10, 3);
  Real c100 = Power(c23, 3);
  Real c97 = Power(c1, -2);
  Real c36 = Power(c35, 3);
  Real c112 = Power(c37, 2);
  Real c102 = Power(c1, -3);
  Real c116 = Power(c37, 3);
  Real c50 = Power(c49, 3);
  Real c52 = Power(c51, -3);
  Real c63 = Power(c62, 3);
  Real c65 = Power(c64, -3);
  Real c76 = Power(c75, 3);
  Real c78 = Power(c77, -3);
  Real c146 = Power(c51, 2);
  Real c152 = Power(c51, 3);
  Real c38 = Power(c37, -3);
  Real c180 = Power(c64, 2);
  Real c188 = Power(c64, 3);
  Real c207 = Power(c77, 2);
  Real c215 = Power(c77, 3);
  Real c2 = Power(c1, -4);
  Real c24 = Power(c23, -3);
  Real c248 = 3 * C0113 * c16 * c26;
  Real c249 = 4 * C0122 * c1 * c10 * c21;
  Real c250 = 2 * C0112 * c16 * c21;
  Real c251 = 3 * C0131 * c14;
  Real c252 = 2 * C0121 * c10;
  Real c253 = C0111 * c1;
  Real c254 = c252 + c253;
  Real c255 = c1 * c254;
  Real c256 = c251 + c255;
  Real c257 = c23 * c256;
  Real c258 = c249 + c250 + c257;
  Real c259 = c23 * c258;
  Real c260 = c248 + c259;
  Real c261 = c102 * c24 * c260;
  Real c90 = Power(c23, -4);
  Real c262 = 3 * C0213 * c16 * c40;
  Real c263 = 4 * C0222 * c1 * c10 * c35;
  Real c264 = 2 * C0212 * c16 * c35;
  Real c265 = 3 * C0231 * c14;
  Real c266 = 2 * C0221 * c10;
  Real c267 = C0211 * c1;
  Real c268 = c266 + c267;
  Real c269 = c1 * c268;
  Real c270 = c265 + c269;
  Real c271 = c270 * c37;
  Real c272 = c263 + c264 + c271;
  Real c273 = c272 * c37;
  Real c274 = c262 + c273;
  Real c275 = c102 * c274 * c38;
  Real c108 = Power(c37, -4);
  Real c276 = 3 * C0313 * c16 * c54;
  Real c277 = 4 * C0322 * c1 * c10 * c49;
  Real c278 = 2 * C0312 * c16 * c49;
  Real c279 = 3 * C0331 * c14;
  Real c280 = 2 * C0321 * c10;
  Real c281 = C0311 * c1;
  Real c282 = c280 + c281;
  Real c283 = c1 * c282;
  Real c284 = c279 + c283;
  Real c285 = c284 * c51;
  Real c286 = c277 + c278 + c285;
  Real c287 = c286 * c51;
  Real c288 = c276 + c287;
  Real c289 = c102 * c288 * c52;
  Real c349 = 3 * C2313 * c112 * c54;
  Real c350 = 4 * C2322 * c35 * c37 * c49;
  Real c351 = 2 * C2312 * c112 * c49;
  Real c352 = 3 * C2331 * c40;
  Real c353 = 2 * C2321 * c35;
  Real c354 = C2311 * c37;
  Real c355 = c353 + c354;
  Real c356 = c355 * c37;
  Real c357 = c352 + c356;
  Real c358 = c357 * c51;
  Real c359 = c350 + c351 + c358;
  Real c360 = c359 * c51;
  Real c361 = c349 + c360;
  Real c362 = c361 * c38 * c52;
  Real c141 = Power(c51, -4);
  Real c290 = 3 * C0413 * c16 * c67;
  Real c291 = 4 * C0422 * c1 * c10 * c62;
  Real c292 = 2 * C0412 * c16 * c62;
  Real c293 = 3 * C0431 * c14;
  Real c294 = 2 * C0421 * c10;
  Real c295 = C0411 * c1;
  Real c296 = c294 + c295;
  Real c297 = c1 * c296;
  Real c298 = c293 + c297;
  Real c299 = c298 * c64;
  Real c300 = c291 + c292 + c299;
  Real c301 = c300 * c64;
  Real c302 = c290 + c301;
  Real c303 = c102 * c302 * c65;
  Real c363 = 3 * C2413 * c112 * c67;
  Real c364 = 4 * C2422 * c35 * c37 * c62;
  Real c365 = 2 * C2412 * c112 * c62;
  Real c366 = 3 * C2431 * c40;
  Real c367 = 2 * C2421 * c35;
  Real c368 = C2411 * c37;
  Real c369 = c367 + c368;
  Real c370 = c369 * c37;
  Real c371 = c366 + c370;
  Real c372 = c371 * c64;
  Real c373 = c364 + c365 + c372;
  Real c374 = c373 * c64;
  Real c375 = c363 + c374;
  Real c376 = c375 * c38 * c65;
  Real c408 = 3 * C3413 * c146 * c67;
  Real c409 = 4 * C3422 * c49 * c51 * c62;
  Real c410 = 2 * C3412 * c146 * c62;
  Real c411 = 3 * C3431 * c54;
  Real c412 = 2 * C3421 * c49 * c51;
  Real c413 = C3411 * c146;
  Real c414 = c411 + c412 + c413;
  Real c415 = c414 * c64;
  Real c416 = c409 + c410 + c415;
  Real c417 = c416 * c64;
  Real c418 = c408 + c417;
  Real c419 = c418 * c52 * c65;
  Real c174 = Power(c64, -4);
  Real c304 = 3 * C0513 * c16 * c80;
  Real c305 = 4 * C0522 * c1 * c10 * c75;
  Real c306 = 2 * C0512 * c16 * c75;
  Real c307 = 3 * C0531 * c14;
  Real c308 = 2 * C0521 * c10;
  Real c309 = C0511 * c1;
  Real c310 = c308 + c309;
  Real c311 = c1 * c310;
  Real c312 = c307 + c311;
  Real c313 = c312 * c77;
  Real c314 = c305 + c306 + c313;
  Real c315 = c314 * c77;
  Real c316 = c304 + c315;
  Real c317 = c102 * c316 * c78;
  Real c377 = 3 * C2513 * c112 * c80;
  Real c378 = 4 * C2522 * c35 * c37 * c75;
  Real c379 = 2 * C2512 * c112 * c75;
  Real c380 = 3 * C2531 * c40;
  Real c381 = 2 * C2521 * c35;
  Real c382 = C2511 * c37;
  Real c383 = c381 + c382;
  Real c384 = c37 * c383;
  Real c385 = c380 + c384;
  Real c386 = c385 * c77;
  Real c387 = c378 + c379 + c386;
  Real c388 = c387 * c77;
  Real c389 = c377 + c388;
  Real c390 = c38 * c389 * c78;
  Real c420 = 3 * C3513 * c146 * c80;
  Real c421 = 4 * C3522 * c49 * c51 * c75;
  Real c422 = 2 * C3512 * c146 * c75;
  Real c423 = 3 * C3531 * c54;
  Real c424 = 2 * C3521 * c49 * c51;
  Real c425 = C3511 * c146;
  Real c426 = c423 + c424 + c425;
  Real c427 = c426 * c77;
  Real c428 = c421 + c422 + c427;
  Real c429 = c428 * c77;
  Real c430 = c420 + c429;
  Real c431 = c430 * c52 * c78;
  Real c201 = Power(c77, -4);
  out1(0) = c2 * (4 * C04 * c11 + 3 * C03 * c1 * c14 + 2 * C02 * c10 * c16 +
                  C01 * c18 + C0113 * c18 * c22 * c24 +
                  2 * C0122 * c10 * c16 * c26 * c27 + C0112 * c18 * c26 * c27 +
                  3 * C0131 * c1 * c14 * c21 * c30 +
                  2 * C0121 * c10 * c16 * c21 * c30 + C0111 * c18 * c21 * c30 +
                  C0213 * c18 * c36 * c38 + 2 * C0222 * c10 * c16 * c40 * c42 +
                  C0212 * c18 * c40 * c42 + 3 * C0231 * c1 * c14 * c35 * c45 +
                  2 * C0221 * c10 * c16 * c35 * c45 + C0211 * c18 * c35 * c45 +
                  C0313 * c18 * c50 * c52 + 2 * C0322 * c10 * c16 * c54 * c55 +
                  C0312 * c18 * c54 * c55 + 3 * C0331 * c1 * c14 * c49 * c58 +
                  2 * C0321 * c10 * c16 * c49 * c58 + C0311 * c18 * c49 * c58 +
                  C0413 * c18 * c63 * c65 + 2 * C0422 * c10 * c16 * c67 * c68 +
                  C0412 * c18 * c67 * c68 + 3 * C0431 * c1 * c14 * c62 * c71 +
                  2 * C0421 * c10 * c16 * c62 * c71 + C0411 * c18 * c62 * c71 +
                  C0513 * c18 * c76 * c78 + 2 * C0522 * c10 * c16 * c80 * c81 +
                  C0512 * c18 * c80 * c81 + 3 * C0531 * c1 * c14 * c75 * c84 +
                  2 * C0521 * c10 * c16 * c75 * c84 + C0511 * c18 * c75 * c84);
  out1(1) =
      c90 * (C11 * c100 + C0131 * c100 * c102 * c11 + 4 * C14 * c22 +
             3 * C13 * c23 * c26 + C0111 * c10 * c100 * c93 +
             3 * C0113 * c10 * c23 * c26 * c93 + 2 * C12 * c21 * c95 +
             2 * C0112 * c10 * c21 * c93 * c95 + C0121 * c100 * c14 * c97 +
             2 * C0122 * c14 * c21 * c95 * c97);
  out1(2) =
      c108 *
      (C21 * c116 + C0231 * c102 * c11 * c116 + 2 * C22 * c112 * c35 +
       4 * C24 * c36 + 3 * C23 * c37 * c40 + C2313 * c116 * c50 * c52 +
       C2312 * c116 * c54 * c55 + 2 * C2322 * c112 * c35 * c54 * c55 +
       C2311 * c116 * c49 * c58 + 2 * C2321 * c112 * c35 * c49 * c58 +
       3 * C2331 * c37 * c40 * c49 * c58 + C2413 * c116 * c63 * c65 +
       C2412 * c116 * c67 * c68 + 2 * C2422 * c112 * c35 * c67 * c68 +
       C2411 * c116 * c62 * c71 + 2 * C2421 * c112 * c35 * c62 * c71 +
       3 * C2431 * c37 * c40 * c62 * c71 + C2513 * c116 * c76 * c78 +
       C2512 * c116 * c80 * c81 + 2 * C2522 * c112 * c35 * c80 * c81 +
       C2511 * c116 * c75 * c84 + 2 * C2521 * c112 * c35 * c75 * c84 +
       3 * C2531 * c37 * c40 * c75 * c84 + C0211 * c10 * c116 * c93 +
       2 * C0212 * c10 * c112 * c35 * c93 + 3 * C0213 * c10 * c37 * c40 * c93 +
       C0221 * c116 * c14 * c97 + 2 * C0222 * c112 * c14 * c35 * c97);
  out1(3) =
      c141 *
      (C31 * c152 + C0331 * c102 * c11 * c152 + C2331 * c152 * c36 * c38 +
       C2321 * c152 * c40 * c42 + C2311 * c152 * c35 * c45 +
       2 * C32 * c146 * c49 + 2 * C2322 * c146 * c40 * c42 * c49 +
       2 * C2312 * c146 * c35 * c45 * c49 + 4 * C34 * c50 +
       3 * C33 * c51 * c54 + 3 * C2313 * c35 * c45 * c51 * c54 +
       C3413 * c152 * c63 * c65 + C3412 * c152 * c67 * c68 +
       2 * C3422 * c146 * c49 * c67 * c68 + C3411 * c152 * c62 * c71 +
       2 * C3421 * c146 * c49 * c62 * c71 + 3 * C3431 * c51 * c54 * c62 * c71 +
       C3513 * c152 * c76 * c78 + C3512 * c152 * c80 * c81 +
       2 * C3522 * c146 * c49 * c80 * c81 + C3511 * c152 * c75 * c84 +
       2 * C3521 * c146 * c49 * c75 * c84 + 3 * C3531 * c51 * c54 * c75 * c84 +
       C0311 * c10 * c152 * c93 + 2 * C0312 * c10 * c146 * c49 * c93 +
       3 * C0313 * c10 * c51 * c54 * c93 + C0321 * c14 * c152 * c97 +
       2 * C0322 * c14 * c146 * c49 * c97);
  out1(4) =
      c174 *
      (C41 * c188 + C0431 * c102 * c11 * c188 + C2431 * c188 * c36 * c38 +
       C2421 * c188 * c40 * c42 + C2411 * c188 * c35 * c45 +
       C3431 * c188 * c50 * c52 + C3421 * c188 * c54 * c55 +
       C3411 * c188 * c49 * c58 + 2 * C42 * c180 * c62 +
       2 * C2422 * c180 * c40 * c42 * c62 + 2 * C2412 * c180 * c35 * c45 * c62 +
       2 * C3422 * c180 * c54 * c55 * c62 + 2 * C3412 * c180 * c49 * c58 * c62 +
       4 * C44 * c63 + 3 * C43 * c64 * c67 + 3 * C2413 * c35 * c45 * c64 * c67 +
       3 * C3413 * c49 * c58 * c64 * c67 + C0411 * c10 * c188 * c93 +
       2 * C0412 * c10 * c180 * c62 * c93 + 3 * C0413 * c10 * c64 * c67 * c93 +
       C0421 * c14 * c188 * c97 + 2 * C0422 * c14 * c180 * c62 * c97);
  out1(5) =
      c201 *
      (C51 * c215 + C0531 * c102 * c11 * c215 + C2531 * c215 * c36 * c38 +
       C2521 * c215 * c40 * c42 + C2511 * c215 * c35 * c45 +
       C3531 * c215 * c50 * c52 + C3521 * c215 * c54 * c55 +
       C3511 * c215 * c49 * c58 + 2 * C52 * c207 * c75 +
       2 * C2522 * c207 * c40 * c42 * c75 + 2 * C2512 * c207 * c35 * c45 * c75 +
       2 * C3522 * c207 * c54 * c55 * c75 + 2 * C3512 * c207 * c49 * c58 * c75 +
       4 * C54 * c76 + 3 * C53 * c77 * c80 + 3 * C2513 * c35 * c45 * c77 * c80 +
       3 * C3513 * c49 * c58 * c77 * c80 + C0511 * c10 * c215 * c93 +
       2 * C0512 * c10 * c207 * c75 * c93 + 3 * C0513 * c10 * c77 * c80 * c93 +
       C0521 * c14 * c215 * c97 + 2 * C0522 * c14 * c207 * c75 * c97);
  out2(0, 0) = 2 * c2 *
               (3 * C03 * c1 * c10 + 6 * C04 * c14 + C02 * c16 +
                C0122 * c16 * c26 * c27 + 3 * C0131 * c1 * c10 * c21 * c30 +
                C0121 * c16 * c21 * c30 + C0222 * c16 * c40 * c42 +
                3 * C0231 * c1 * c10 * c35 * c45 + C0221 * c16 * c35 * c45 +
                C0322 * c16 * c54 * c55 + 3 * C0331 * c1 * c10 * c49 * c58 +
                C0321 * c16 * c49 * c58 + C0422 * c16 * c67 * c68 +
                3 * C0431 * c1 * c10 * c62 * c71 + C0421 * c16 * c62 * c71 +
                C0522 * c16 * c80 * c81 + 3 * C0531 * c1 * c10 * c75 * c84 +
                C0521 * c16 * c75 * c84);
  out2(0, 1) = c261;
  out2(0, 2) = c275;
  out2(0, 3) = c289;
  out2(0, 4) = c303;
  out2(0, 5) = c317;
  out2(1, 0) = c261;
  out2(1, 1) = 2 * c90 *
               (6 * C14 * c26 +
                c23 *
                    (3 * C0113 * c1 * c10 * c21 + 3 * C13 * c16 * c21 +
                     (c1 * (C12 * c1 + C0112 * c10) + C0122 * c14) * c23) *
                    c97);
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = c275;
  out2(2, 1) = 0;
  out2(2, 2) = 2 * c108 *
               (C22 * c112 + 3 * C23 * c35 * c37 + 6 * C24 * c40 +
                C2322 * c112 * c54 * c55 + C2321 * c112 * c49 * c58 +
                3 * C2331 * c35 * c37 * c49 * c58 + C2422 * c112 * c67 * c68 +
                C2421 * c112 * c62 * c71 + 3 * C2431 * c35 * c37 * c62 * c71 +
                C2522 * c112 * c80 * c81 + C2521 * c112 * c75 * c84 +
                3 * C2531 * c35 * c37 * c75 * c84 + C0212 * c10 * c112 * c93 +
                3 * C0213 * c10 * c35 * c37 * c93 + C0222 * c112 * c14 * c97);
  out2(2, 3) = c362;
  out2(2, 4) = c376;
  out2(2, 5) = c390;
  out2(3, 0) = c289;
  out2(3, 1) = 0;
  out2(3, 2) = c362;
  out2(3, 3) =
      2 * c141 *
      (C32 * c146 + C2322 * c146 * c40 * c42 + C2312 * c146 * c35 * c45 +
       3 * C33 * c49 * c51 + 3 * C2313 * c35 * c45 * c49 * c51 + 6 * C34 * c54 +
       C3422 * c146 * c67 * c68 + C3421 * c146 * c62 * c71 +
       3 * C3431 * c49 * c51 * c62 * c71 + C3522 * c146 * c80 * c81 +
       C3521 * c146 * c75 * c84 + 3 * C3531 * c49 * c51 * c75 * c84 +
       C0312 * c10 * c146 * c93 + 3 * C0313 * c10 * c49 * c51 * c93 +
       C0322 * c14 * c146 * c97);
  out2(3, 4) = c419;
  out2(3, 5) = c431;
  out2(4, 0) = c303;
  out2(4, 1) = 0;
  out2(4, 2) = c376;
  out2(4, 3) = c419;
  out2(4, 4) =
      2 * c174 *
      (C42 * c180 + C2422 * c180 * c40 * c42 + C2412 * c180 * c35 * c45 +
       C3422 * c180 * c54 * c55 + C3412 * c180 * c49 * c58 +
       3 * C43 * c62 * c64 + 3 * C2413 * c35 * c45 * c62 * c64 +
       3 * C3413 * c49 * c58 * c62 * c64 + 6 * C44 * c67 +
       C0412 * c10 * c180 * c93 + 3 * C0413 * c10 * c62 * c64 * c93 +
       C0422 * c14 * c180 * c97);
  out2(4, 5) = 0;
  out2(5, 0) = c317;
  out2(5, 1) = 0;
  out2(5, 2) = c390;
  out2(5, 3) = c431;
  out2(5, 4) = 0;
  out2(5, 5) =
      2 * c201 *
      (C52 * c207 + C2522 * c207 * c40 * c42 + C2512 * c207 * c35 * c45 +
       C3522 * c207 * c54 * c55 + C3512 * c207 * c49 * c58 +
       3 * C53 * c75 * c77 + 3 * C2513 * c35 * c45 * c75 * c77 +
       3 * C3513 * c49 * c58 * c75 * c77 + 6 * C54 * c80 +
       C0512 * c10 * c207 * c93 + 3 * C0513 * c10 * c75 * c77 * c93 +
       C0522 * c14 * c207 * c97);

  return std::make_pair(hess, grad);
}

value_type FittedMaterial::psi_barrier(const strain_type &ek,
                                       const Vec6 &ekclamped) {
  return s * (Power(ek(0) - ekclamped(0), 6) / Power(ekscale(0), 6) +
              Power(ek(1) - ekclamped(1), 6) / Power(ekscale(1), 6) +
              Power(ek(2) - ekclamped(2), 6) / Power(ekscale(2), 6) +
              Power(ek(3) - ekclamped(3), 6) / Power(ekscale(3), 6) +
              Power(ek(4) - ekclamped(4), 6) / Power(ekscale(4), 6) +
              Power(ek(5) - ekclamped(5), 6) / Power(ekscale(5), 6));
}
grad_type FittedMaterial::grad_barrier(const strain_type &ek,
                                       const Vec6 &ekclamped) {
  Vec6 out(0);
  out(0) = (6 * s * Power(ek(0) - ekclamped(0), 5)) / Power(ekscale(0), 6);
  out(1) = (6 * s * Power(ek(1) - ekclamped(1), 5)) / Power(ekscale(1), 6);
  out(2) = (6 * s * Power(ek(2) - ekclamped(2), 5)) / Power(ekscale(2), 6);
  out(3) = (6 * s * Power(ek(3) - ekclamped(3), 5)) / Power(ekscale(3), 6);
  out(4) = (6 * s * Power(ek(4) - ekclamped(4), 5)) / Power(ekscale(4), 6);
  out(5) = (6 * s * Power(ek(5) - ekclamped(5), 5)) / Power(ekscale(5), 6);

  return out;
}
gradhess_type FittedMaterial::gradhess_barrier(const strain_type &ek,
                                               const Vec6 &ekclamped) {
  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real c1 = ek(0);
  Real c2 = ekclamped(0);
  Real c4 = -c2;
  Real c5 = c1 + c4;
  Real c7 = ekscale(0);
  Real c8 = Power(c7, -6);
  Real c10 = ek(1);
  Real c11 = ekclamped(1);
  Real c12 = -c11;
  Real c14 = c10 + c12;
  Real c16 = ekscale(1);
  Real c17 = Power(c16, -6);
  Real c19 = ek(2);
  Real c21 = ekclamped(2);
  Real c22 = -c21;
  Real c23 = c19 + c22;
  Real c25 = ekscale(2);
  Real c26 = Power(c25, -6);
  Real c28 = ek(3);
  Real c29 = ekclamped(3);
  Real c30 = -c29;
  Real c31 = c28 + c30;
  Real c33 = ekscale(3);
  Real c34 = Power(c33, -6);
  Real c36 = ek(4);
  Real c37 = ekclamped(4);
  Real c38 = -c37;
  Real c39 = c36 + c38;
  Real c42 = ekscale(4);
  Real c43 = Power(c42, -6);
  Real c45 = ek(5);
  Real c46 = ekclamped(5);
  Real c47 = -c46;
  Real c48 = c45 + c47;
  Real c50 = ekscale(5);
  Real c51 = Power(c50, -6);
  out1(0) = 6 * s * Power(c5, 5) * c8;
  out1(1) = 6 * s * Power(c14, 5) * c17;
  out1(2) = 6 * s * Power(c23, 5) * c26;
  out1(3) = 6 * s * Power(c31, 5) * c34;
  out1(4) = 6 * s * Power(c39, 5) * c43;
  out1(5) = 6 * s * Power(c48, 5) * c51;
  out2(0, 0) = 30 * s * Power(c5, 4) * c8;
  out2(0, 1) = 0;
  out2(0, 2) = 0;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) = 30 * s * Power(c14, 4) * c17;
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) = 30 * s * Power(c23, 4) * c26;
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = 30 * s * Power(c31, 4) * c34;
  out2(3, 4) = 0;
  out2(3, 5) = 0;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = 0;
  out2(4, 4) = 30 * s * Power(c39, 4) * c43;
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) = 30 * s * Power(c48, 4) * c51;
  return std::make_pair(hess, grad);
}

// 1st and 2nd derivative of psi, grad, tuple<grad,hess> for 2nd order taylor
// wrt ek[0 ... 5]
// -------------- values -------------- //
std::pair<value_type, value_type>
FittedMaterial::psi_taylor_12_0(const strain_type &ek) {
  auto _out = std::make_pair<value_type, value_type>(0, 0);
  auto out = [&](int i) -> Real & {
    return (i == 0) ? _out.first : _out.second;
  };
  Real c1 = ek(0);
  Real c2 = -1 + c1;
  Real c11 = ekscale(0);
  Real c21 = 1 / c11;
  Real c23 = ek(1);
  Real c18 = Power(c11, -2);
  Real c25 = ekscale(1);
  Real c28 = Power(c23, 2);
  Real c29 = Power(c25, -2);
  Real c15 = Power(c2, 2);
  Real c16 = Power(c11, -3);
  Real c32 = 1 / c25;
  Real c36 = ek(2);
  Real c37 = -1 + c36;
  Real c39 = ekscale(2);
  Real c43 = Power(c37, 2);
  Real c44 = Power(c39, -2);
  Real c47 = 1 / c39;
  Real c51 = ek(3);
  Real c53 = ekscale(3);
  Real c56 = Power(c51, 2);
  Real c57 = Power(c53, -2);
  Real c60 = 1 / c53;
  Real c64 = ek(4);
  Real c66 = ekscale(4);
  Real c69 = Power(c64, 2);
  Real c70 = Power(c66, -2);
  Real c73 = 1 / c66;
  Real c77 = ek(5);
  Real c79 = ekscale(5);
  Real c82 = Power(c77, 2);
  Real c83 = Power(c79, -2);
  Real c86 = 1 / c79;
  Real c12 = Power(c11, -4);
  out(0) =
      3 * C03 * c15 * c16 + 2 * C02 * c18 * c2 + 4 * C04 * c12 * Power(c2, 3) +
      C01 * c21 + (C0113 * c21 * Power(c23, 3)) / Power(c25, 3) +
      2 * C0122 * c18 * c2 * c28 * c29 + C0112 * c21 * c28 * c29 +
      3 * C0131 * c15 * c16 * c23 * c32 + 2 * C0121 * c18 * c2 * c23 * c32 +
      C0111 * c21 * c23 * c32 + (C0213 * c21 * Power(c37, 3)) / Power(c39, 3) +
      2 * C0222 * c18 * c2 * c43 * c44 + C0212 * c21 * c43 * c44 +
      3 * C0231 * c15 * c16 * c37 * c47 + 2 * C0221 * c18 * c2 * c37 * c47 +
      C0211 * c21 * c37 * c47 + (C0313 * c21 * Power(c51, 3)) / Power(c53, 3) +
      2 * C0322 * c18 * c2 * c56 * c57 + C0312 * c21 * c56 * c57 +
      3 * C0331 * c15 * c16 * c51 * c60 + 2 * C0321 * c18 * c2 * c51 * c60 +
      C0311 * c21 * c51 * c60 + (C0413 * c21 * Power(c64, 3)) / Power(c66, 3) +
      2 * C0422 * c18 * c2 * c69 * c70 + C0412 * c21 * c69 * c70 +
      3 * C0431 * c15 * c16 * c64 * c73 + 2 * C0421 * c18 * c2 * c64 * c73 +
      C0411 * c21 * c64 * c73 + (C0513 * c21 * Power(c77, 3)) / Power(c79, 3) +
      2 * C0522 * c18 * c2 * c82 * c83 + C0512 * c21 * c82 * c83 +
      3 * C0531 * c15 * c16 * c77 * c86 + 2 * C0521 * c18 * c2 * c77 * c86 +
      C0511 * c21 * c77 * c86;
  out(1) = 12 * C04 * c12 * c15 + 2 * C02 * c18 + 6 * C03 * c16 * c2 +
           2 * C0122 * c18 * c28 * c29 + 2 * C0121 * c18 * c23 * c32 +
           6 * C0131 * c16 * c2 * c23 * c32 + 2 * C0222 * c18 * c43 * c44 +
           2 * C0221 * c18 * c37 * c47 + 6 * C0231 * c16 * c2 * c37 * c47 +
           2 * C0322 * c18 * c56 * c57 + 2 * C0321 * c18 * c51 * c60 +
           6 * C0331 * c16 * c2 * c51 * c60 + 2 * C0422 * c18 * c69 * c70 +
           2 * C0421 * c18 * c64 * c73 + 6 * C0431 * c16 * c2 * c64 * c73 +
           2 * C0522 * c18 * c82 * c83 + 2 * C0521 * c18 * c77 * c86 +
           6 * C0531 * c16 * c2 * c77 * c86;
  return _out;
}
std::pair<value_type, value_type>
FittedMaterial::psi_taylor_12_1(const strain_type &ek) {
  auto _out = std::make_pair<value_type, value_type>(0, 0);
  auto out = [&](int i) -> Real & {
    return (i == 0) ? _out.first : _out.second;
  };
  Real c1 = ek(1);
  Real c10 = ekscale(1);
  Real c14 = Power(c1, 2);
  Real c15 = Power(c10, -3);
  Real c17 = ek(0);
  Real c18 = -1 + c17;
  Real c19 = ekscale(0);
  Real c23 = Power(c10, -2);
  Real c21 = 1 / c19;
  Real c29 = 1 / c10;
  Real c25 = Power(c18, 2);
  Real c26 = Power(c19, -2);
  Real c11 = Power(c10, -4);
  out(0) = 4 * C14 * Power(c1, 3) * c11 + 3 * C13 * c14 * c15 +
           3 * C0113 * c14 * c15 * c18 * c21 + 2 * C12 * c1 * c23 +
           2 * C0112 * c1 * c18 * c21 * c23 + 2 * C0122 * c1 * c23 * c25 * c26 +
           C11 * c29 + (C0131 * Power(c18, 3) * c29) / Power(c19, 3) +
           C0111 * c18 * c21 * c29 + C0121 * c25 * c26 * c29;
  out(1) = 12 * C14 * c11 * c14 + 6 * C13 * c1 * c15 +
           6 * C0113 * c1 * c15 * c18 * c21 + 2 * C12 * c23 +
           2 * C0112 * c18 * c21 * c23 + 2 * C0122 * c23 * c25 * c26;
  return _out;
}
std::pair<value_type, value_type>
FittedMaterial::psi_taylor_12_2(const strain_type &ek) {
  auto _out = std::make_pair<value_type, value_type>(0, 0);
  auto out = [&](int i) -> Real & {
    return (i == 0) ? _out.first : _out.second;
  };
  Real c1 = ek(2);
  Real c2 = -1 + c1;
  Real c11 = ekscale(2);
  Real c15 = Power(c2, 2);
  Real c16 = Power(c11, -3);
  Real c18 = ek(0);
  Real c19 = -1 + c18;
  Real c21 = ekscale(0);
  Real c24 = Power(c11, -2);
  Real c22 = 1 / c21;
  Real c30 = 1 / c11;
  Real c26 = Power(c19, 2);
  Real c27 = Power(c21, -2);
  Real c37 = ek(3);
  Real c39 = ekscale(3);
  Real c43 = Power(c37, 2);
  Real c44 = Power(c39, -2);
  Real c47 = 1 / c39;
  Real c51 = ek(4);
  Real c53 = ekscale(4);
  Real c56 = Power(c51, 2);
  Real c57 = Power(c53, -2);
  Real c60 = 1 / c53;
  Real c64 = ek(5);
  Real c66 = ekscale(5);
  Real c69 = Power(c64, 2);
  Real c70 = Power(c66, -2);
  Real c73 = 1 / c66;
  Real c12 = Power(c11, -4);
  out(0) =
      3 * C23 * c15 * c16 + 4 * C24 * c12 * Power(c2, 3) +
      3 * C0213 * c15 * c16 * c19 * c22 + 2 * C22 * c2 * c24 +
      2 * C0212 * c19 * c2 * c22 * c24 + 2 * C0222 * c2 * c24 * c26 * c27 +
      C21 * c30 + (C0231 * Power(c19, 3) * c30) / Power(c21, 3) +
      C0211 * c19 * c22 * c30 + C0221 * c26 * c27 * c30 +
      (C2313 * c30 * Power(c37, 3)) / Power(c39, 3) +
      2 * C2322 * c2 * c24 * c43 * c44 + C2312 * c30 * c43 * c44 +
      3 * C2331 * c15 * c16 * c37 * c47 + 2 * C2321 * c2 * c24 * c37 * c47 +
      C2311 * c30 * c37 * c47 + (C2413 * c30 * Power(c51, 3)) / Power(c53, 3) +
      2 * C2422 * c2 * c24 * c56 * c57 + C2412 * c30 * c56 * c57 +
      3 * C2431 * c15 * c16 * c51 * c60 + 2 * C2421 * c2 * c24 * c51 * c60 +
      C2411 * c30 * c51 * c60 + (C2513 * c30 * Power(c64, 3)) / Power(c66, 3) +
      2 * C2522 * c2 * c24 * c69 * c70 + C2512 * c30 * c69 * c70 +
      3 * C2531 * c15 * c16 * c64 * c73 + 2 * C2521 * c2 * c24 * c64 * c73 +
      C2511 * c30 * c64 * c73;
  out(1) = 12 * C24 * c12 * c15 + 6 * C23 * c16 * c2 +
           6 * C0213 * c16 * c19 * c2 * c22 + 2 * C22 * c24 +
           2 * C0212 * c19 * c22 * c24 + 2 * C0222 * c24 * c26 * c27 +
           2 * C2322 * c24 * c43 * c44 + 6 * C2331 * c16 * c2 * c37 * c47 +
           2 * C2321 * c24 * c37 * c47 + 2 * C2422 * c24 * c56 * c57 +
           6 * C2431 * c16 * c2 * c51 * c60 + 2 * C2421 * c24 * c51 * c60 +
           2 * C2522 * c24 * c69 * c70 + 6 * C2531 * c16 * c2 * c64 * c73 +
           2 * C2521 * c24 * c64 * c73;
  return _out;
}
std::pair<value_type, value_type>
FittedMaterial::psi_taylor_12_3(const strain_type &ek) {
  auto _out = std::make_pair<value_type, value_type>(0, 0);
  auto out = [&](int i) -> Real & {
    return (i == 0) ? _out.first : _out.second;
  };
  Real c1 = ek(3);
  Real c10 = ekscale(3);
  Real c14 = Power(c1, 2);
  Real c15 = Power(c10, -3);
  Real c17 = ek(0);
  Real c18 = -1 + c17;
  Real c19 = ekscale(0);
  Real c28 = Power(c10, -2);
  Real c21 = 1 / c19;
  Real c23 = ek(2);
  Real c24 = -1 + c23;
  Real c25 = ekscale(2);
  Real c26 = 1 / c25;
  Real c38 = 1 / c10;
  Real c30 = Power(c18, 2);
  Real c31 = Power(c19, -2);
  Real c34 = Power(c24, 2);
  Real c35 = Power(c25, -2);
  Real c51 = ek(4);
  Real c53 = ekscale(4);
  Real c56 = Power(c51, 2);
  Real c57 = Power(c53, -2);
  Real c60 = 1 / c53;
  Real c64 = ek(5);
  Real c66 = ekscale(5);
  Real c69 = Power(c64, 2);
  Real c70 = Power(c66, -2);
  Real c73 = 1 / c66;
  Real c11 = Power(c10, -4);
  out(0) =
      4 * C34 * Power(c1, 3) * c11 + 3 * C33 * c14 * c15 +
      3 * C0313 * c14 * c15 * c18 * c21 + 3 * C2313 * c14 * c15 * c24 * c26 +
      2 * C32 * c1 * c28 + 2 * C0312 * c1 * c18 * c21 * c28 +
      2 * C2312 * c1 * c24 * c26 * c28 + 2 * C0322 * c1 * c28 * c30 * c31 +
      2 * C2322 * c1 * c28 * c34 * c35 + C31 * c38 +
      (C0331 * Power(c18, 3) * c38) / Power(c19, 3) + C0311 * c18 * c21 * c38 +
      (C2331 * Power(c24, 3) * c38) / Power(c25, 3) + C2311 * c24 * c26 * c38 +
      C0321 * c30 * c31 * c38 + C2321 * c34 * c35 * c38 +
      (C3413 * c38 * Power(c51, 3)) / Power(c53, 3) +
      2 * C3422 * c1 * c28 * c56 * c57 + C3412 * c38 * c56 * c57 +
      3 * C3431 * c14 * c15 * c51 * c60 + 2 * C3421 * c1 * c28 * c51 * c60 +
      C3411 * c38 * c51 * c60 + (C3513 * c38 * Power(c64, 3)) / Power(c66, 3) +
      2 * C3522 * c1 * c28 * c69 * c70 + C3512 * c38 * c69 * c70 +
      3 * C3531 * c14 * c15 * c64 * c73 + 2 * C3521 * c1 * c28 * c64 * c73 +
      C3511 * c38 * c64 * c73;
  out(1) = 12 * C34 * c11 * c14 + 6 * C33 * c1 * c15 +
           6 * C0313 * c1 * c15 * c18 * c21 + 6 * C2313 * c1 * c15 * c24 * c26 +
           2 * C32 * c28 + 2 * C0312 * c18 * c21 * c28 +
           2 * C2312 * c24 * c26 * c28 + 2 * C0322 * c28 * c30 * c31 +
           2 * C2322 * c28 * c34 * c35 + 2 * C3422 * c28 * c56 * c57 +
           6 * C3431 * c1 * c15 * c51 * c60 + 2 * C3421 * c28 * c51 * c60 +
           2 * C3522 * c28 * c69 * c70 + 6 * C3531 * c1 * c15 * c64 * c73 +
           2 * C3521 * c28 * c64 * c73;
  return _out;
}
std::pair<value_type, value_type>
FittedMaterial::psi_taylor_12_4(const strain_type &ek) {
  auto _out = std::make_pair<value_type, value_type>(0, 0);
  auto out = [&](int i) -> Real & {
    return (i == 0) ? _out.first : _out.second;
  };
  Real c1 = ek(4);
  Real c10 = ekscale(4);
  Real c14 = Power(c1, 2);
  Real c15 = Power(c10, -3);
  Real c17 = ek(0);
  Real c18 = -1 + c17;
  Real c19 = ekscale(0);
  Real c32 = Power(c10, -2);
  Real c21 = 1 / c19;
  Real c23 = ek(2);
  Real c24 = -1 + c23;
  Real c25 = ekscale(2);
  Real c26 = 1 / c25;
  Real c28 = ek(3);
  Real c29 = ekscale(3);
  Real c30 = 1 / c29;
  Real c47 = 1 / c10;
  Real c34 = Power(c18, 2);
  Real c35 = Power(c19, -2);
  Real c38 = Power(c24, 2);
  Real c39 = Power(c25, -2);
  Real c43 = Power(c28, 2);
  Real c44 = Power(c29, -2);
  Real c11 = Power(c10, -4);
  out(0) =
      4 * C44 * Power(c1, 3) * c11 + 3 * C43 * c14 * c15 +
      3 * C0413 * c14 * c15 * c18 * c21 + 3 * C2413 * c14 * c15 * c24 * c26 +
      3 * C3413 * c14 * c15 * c28 * c30 + 2 * C42 * c1 * c32 +
      2 * C0412 * c1 * c18 * c21 * c32 + 2 * C2412 * c1 * c24 * c26 * c32 +
      2 * C3412 * c1 * c28 * c30 * c32 + 2 * C0422 * c1 * c32 * c34 * c35 +
      2 * C2422 * c1 * c32 * c38 * c39 + 2 * C3422 * c1 * c32 * c43 * c44 +
      C41 * c47 + (C0431 * Power(c18, 3) * c47) / Power(c19, 3) +
      C0411 * c18 * c21 * c47 + (C2431 * Power(c24, 3) * c47) / Power(c25, 3) +
      C2411 * c24 * c26 * c47 + (C3431 * Power(c28, 3) * c47) / Power(c29, 3) +
      C3411 * c28 * c30 * c47 + C0421 * c34 * c35 * c47 +
      C2421 * c38 * c39 * c47 + C3421 * c43 * c44 * c47;
  out(1) = 12 * C44 * c11 * c14 + 6 * C43 * c1 * c15 +
           6 * C0413 * c1 * c15 * c18 * c21 + 6 * C2413 * c1 * c15 * c24 * c26 +
           6 * C3413 * c1 * c15 * c28 * c30 + 2 * C42 * c32 +
           2 * C0412 * c18 * c21 * c32 + 2 * C2412 * c24 * c26 * c32 +
           2 * C3412 * c28 * c30 * c32 + 2 * C0422 * c32 * c34 * c35 +
           2 * C2422 * c32 * c38 * c39 + 2 * C3422 * c32 * c43 * c44;
  return _out;
}
std::pair<value_type, value_type>
FittedMaterial::psi_taylor_12_5(const strain_type &ek) {
  auto _out = std::make_pair<value_type, value_type>(0, 0);
  auto out = [&](int i) -> Real & {
    return (i == 0) ? _out.first : _out.second;
  };
  Real c1 = ek(5);
  Real c10 = ekscale(5);
  Real c14 = Power(c1, 2);
  Real c15 = Power(c10, -3);
  Real c17 = ek(0);
  Real c18 = -1 + c17;
  Real c19 = ekscale(0);
  Real c32 = Power(c10, -2);
  Real c21 = 1 / c19;
  Real c23 = ek(2);
  Real c24 = -1 + c23;
  Real c25 = ekscale(2);
  Real c26 = 1 / c25;
  Real c28 = ek(3);
  Real c29 = ekscale(3);
  Real c30 = 1 / c29;
  Real c47 = 1 / c10;
  Real c34 = Power(c18, 2);
  Real c35 = Power(c19, -2);
  Real c38 = Power(c24, 2);
  Real c39 = Power(c25, -2);
  Real c43 = Power(c28, 2);
  Real c44 = Power(c29, -2);
  Real c11 = Power(c10, -4);
  out(0) =
      4 * C54 * Power(c1, 3) * c11 + 3 * C53 * c14 * c15 +
      3 * C0513 * c14 * c15 * c18 * c21 + 3 * C2513 * c14 * c15 * c24 * c26 +
      3 * C3513 * c14 * c15 * c28 * c30 + 2 * C52 * c1 * c32 +
      2 * C0512 * c1 * c18 * c21 * c32 + 2 * C2512 * c1 * c24 * c26 * c32 +
      2 * C3512 * c1 * c28 * c30 * c32 + 2 * C0522 * c1 * c32 * c34 * c35 +
      2 * C2522 * c1 * c32 * c38 * c39 + 2 * C3522 * c1 * c32 * c43 * c44 +
      C51 * c47 + (C0531 * Power(c18, 3) * c47) / Power(c19, 3) +
      C0511 * c18 * c21 * c47 + (C2531 * Power(c24, 3) * c47) / Power(c25, 3) +
      C2511 * c24 * c26 * c47 + (C3531 * Power(c28, 3) * c47) / Power(c29, 3) +
      C3511 * c28 * c30 * c47 + C0521 * c34 * c35 * c47 +
      C2521 * c38 * c39 * c47 + C3521 * c43 * c44 * c47;
  out(1) = 12 * C54 * c11 * c14 + 6 * C53 * c1 * c15 +
           6 * C0513 * c1 * c15 * c18 * c21 + 6 * C2513 * c1 * c15 * c24 * c26 +
           6 * C3513 * c1 * c15 * c28 * c30 + 2 * C52 * c32 +
           2 * C0512 * c18 * c21 * c32 + 2 * C2512 * c24 * c26 * c32 +
           2 * C3512 * c28 * c30 * c32 + 2 * C0522 * c32 * c34 * c35 +
           2 * C2522 * c32 * c38 * c39 + 2 * C3522 * c32 * c43 * c44;
  return _out;
}
// -------------- grads -------------- //
std::pair<grad_type, grad_type>
FittedMaterial::grad_taylor_12_0(const strain_type &ek) {
  auto out1 = grad_type(0);
  auto out2 = grad_type(0);
  Real c3 = ek(0);
  Real c4 = -1 + c3;
  Real c1 = ekscale(0);
  Real c8 = Power(c1, 2);
  Real c10 = ek(1);
  Real c12 = ekscale(1);
  Real c16 = 1 / c12;
  Real c19 = ek(2);
  Real c21 = -1 + c19;
  Real c23 = ekscale(2);
  Real c26 = 1 / c23;
  Real c29 = ek(3);
  Real c31 = ekscale(3);
  Real c34 = 1 / c31;
  Real c37 = ek(4);
  Real c39 = ekscale(4);
  Real c43 = 1 / c39;
  Real c46 = ek(5);
  Real c48 = ekscale(5);
  Real c51 = 1 / c48;
  Real c11 = Power(c10, 2);
  Real c57 = 1 / c1;
  Real c60 = Power(c12, 2);
  Real c5 = Power(c4, 2);
  Real c59 = Power(c1, -2);
  Real c64 = Power(c12, 3);
  Real c22 = Power(c21, 2);
  Real c72 = Power(c23, 2);
  Real c63 = Power(c1, -3);
  Real c75 = Power(c23, 3);
  Real c30 = Power(c29, 2);
  Real c83 = Power(c31, 2);
  Real c86 = Power(c31, 3);
  Real c38 = Power(c37, 2);
  Real c94 = Power(c39, 2);
  Real c97 = Power(c39, 3);
  Real c47 = Power(c46, 2);
  Real c105 = Power(c48, 2);
  Real c108 = Power(c48, 3);
  Real c2 = Power(c1, -4);
  Real c56 = Power(c12, -4);
  Real c70 = Power(c23, -4);
  Real c81 = Power(c31, -4);
  Real c92 = Power(c39, -4);
  Real c103 = Power(c48, -4);
  out1(0) =
      c2 *
      (6 * C03 * c1 * c4 + 6 * C0131 * c1 * c10 * c16 * c4 +
       6 * C0231 * c1 * c21 * c26 * c4 + 6 * C0331 * c1 * c29 * c34 * c4 +
       6 * C0431 * c1 * c37 * c4 * c43 + 12 * C04 * c5 +
       6 * C0531 * c1 * c4 * c46 * c51 + 2 * C02 * c8 +
       (2 * C0122 * c11 * c8) / Power(c12, 2) + 2 * C0121 * c10 * c16 * c8 +
       (2 * C0222 * c22 * c8) / Power(c23, 2) + 2 * C0221 * c21 * c26 * c8 +
       (2 * C0322 * c30 * c8) / Power(c31, 2) + 2 * C0321 * c29 * c34 * c8 +
       (2 * C0422 * c38 * c8) / Power(c39, 2) + 2 * C0421 * c37 * c43 * c8 +
       (2 * C0522 * c47 * c8) / Power(c48, 2) + 2 * C0521 * c46 * c51 * c8);
  out1(1) = c56 * (3 * C0113 * c11 * c12 * c57 + 2 * C0112 * c10 * c57 * c60 +
                   4 * C0122 * c10 * c4 * c59 * c60 + C0111 * c57 * c64 +
                   2 * C0121 * c4 * c59 * c64 + 3 * C0131 * c5 * c63 * c64);
  out1(2) = c70 * (3 * C0213 * c22 * c23 * c57 + 2 * C0212 * c21 * c57 * c72 +
                   4 * C0222 * c21 * c4 * c59 * c72 + C0211 * c57 * c75 +
                   2 * C0221 * c4 * c59 * c75 + 3 * C0231 * c5 * c63 * c75);
  out1(3) = c81 * (3 * C0313 * c30 * c31 * c57 + 2 * C0312 * c29 * c57 * c83 +
                   4 * C0322 * c29 * c4 * c59 * c83 + C0311 * c57 * c86 +
                   2 * C0321 * c4 * c59 * c86 + 3 * C0331 * c5 * c63 * c86);
  out1(4) = c92 * (3 * C0413 * c38 * c39 * c57 + 2 * C0412 * c37 * c57 * c94 +
                   4 * C0422 * c37 * c4 * c59 * c94 + C0411 * c57 * c97 +
                   2 * C0421 * c4 * c59 * c97 + 3 * C0431 * c5 * c63 * c97);
  out1(5) =
      c103 * (C0511 * c108 * c57 + 2 * C0512 * c105 * c46 * c57 +
              3 * C0513 * c47 * c48 * c57 + 2 * C0521 * c108 * c4 * c59 +
              4 * C0522 * c105 * c4 * c46 * c59 + 3 * C0531 * c108 * c5 * c63);
  out2(0) = c2 * (6 * C03 * c1 + 6 * C0131 * c1 * c10 * c16 +
                  6 * C0231 * c1 * c21 * c26 + 6 * C0331 * c1 * c29 * c34 +
                  24 * C04 * c4 + 6 * C0431 * c1 * c37 * c43 +
                  6 * C0531 * c1 * c46 * c51);
  out2(1) = c56 * (4 * C0122 * c10 * c59 * c60 + 2 * C0121 * c59 * c64 +
                   6 * C0131 * c4 * c63 * c64);
  out2(2) = c70 * (4 * C0222 * c21 * c59 * c72 + 2 * C0221 * c59 * c75 +
                   6 * C0231 * c4 * c63 * c75);
  out2(3) = c81 * (4 * C0322 * c29 * c59 * c83 + 2 * C0321 * c59 * c86 +
                   6 * C0331 * c4 * c63 * c86);
  out2(4) = c92 * (4 * C0422 * c37 * c59 * c94 + 2 * C0421 * c59 * c97 +
                   6 * C0431 * c4 * c63 * c97);
  out2(5) = c103 * (2 * C0521 * c108 * c59 + 4 * C0522 * c105 * c46 * c59 +
                    6 * C0531 * c108 * c4 * c63);
  return std::make_pair(out1, out2);
}
std::pair<grad_type, grad_type>
FittedMaterial::grad_taylor_12_1(const strain_type &ek) {
  auto out1 = grad_type(0);
  auto out2 = grad_type(0);
  Real c1 = ekscale(0);
  Real c3 = ek(1);
  Real c6 = ekscale(1);
  Real c5 = Power(c1, 3);
  Real c12 = Power(c6, -2);
  Real c9 = ek(0);
  Real c10 = -1 + c9;
  Real c11 = Power(c1, 2);
  Real c17 = 1 / c6;
  Real c4 = Power(c3, 2);
  Real c16 = Power(c10, 2);
  Real c29 = Power(c6, 2);
  Real c27 = 1 / c1;
  Real c2 = Power(c1, -4);
  Real c7 = Power(c6, -3);
  Real c24 = Power(c6, -4);
  out1(0) = c2 * (2 * C0121 * c10 * c11 * c17 + 3 * C0131 * c1 * c16 * c17 +
                  4 * C0122 * c10 * c11 * c12 * c3 + C0111 * c17 * c5 +
                  2 * C0112 * c12 * c3 * c5 + 3 * C0113 * c4 * c5 * c7);
  out1(1) = c24 * (2 * C12 * c29 + (2 * C0122 * c16 * c29) / Power(c1, 2) +
                   2 * C0112 * c10 * c27 * c29 + 12 * C14 * c4 +
                   6 * C13 * c3 * c6 + 6 * C0113 * c10 * c27 * c3 * c6);
  out1(2) = 0;
  out1(3) = 0;
  out1(4) = 0;
  out1(5) = 0;
  out2(0) = c2 * (4 * C0122 * c10 * c11 * c12 + 2 * C0112 * c12 * c5 +
                  6 * C0113 * c3 * c5 * c7);
  out2(1) = c24 * (24 * C14 * c3 + 6 * C13 * c6 + 6 * C0113 * c10 * c27 * c6);
  out2(2) = 0;
  out2(3) = 0;
  out2(4) = 0;
  out2(5) = 0;
  return std::make_pair(out1, out2);
}
std::pair<grad_type, grad_type>
FittedMaterial::grad_taylor_12_2(const strain_type &ek) {
  auto out1 = grad_type(0);
  auto out2 = grad_type(0);
  Real c1 = ekscale(0);
  Real c3 = ek(2);
  Real c4 = -1 + c3;
  Real c7 = ekscale(2);
  Real c6 = Power(c1, 3);
  Real c14 = Power(c7, -2);
  Real c10 = ek(0);
  Real c11 = -1 + c10;
  Real c12 = Power(c1, 2);
  Real c18 = 1 / c7;
  Real c5 = Power(c4, 2);
  Real c17 = Power(c11, 2);
  Real c30 = Power(c7, 2);
  Real c28 = 1 / c1;
  Real c35 = ek(3);
  Real c37 = ekscale(3);
  Real c40 = 1 / c37;
  Real c44 = ek(4);
  Real c46 = ekscale(4);
  Real c49 = 1 / c46;
  Real c52 = ek(5);
  Real c54 = ekscale(5);
  Real c57 = 1 / c54;
  Real c36 = Power(c35, 2);
  Real c64 = Power(c37, 2);
  Real c8 = Power(c7, -3);
  Real c67 = Power(c37, 3);
  Real c45 = Power(c44, 2);
  Real c75 = Power(c46, 2);
  Real c78 = Power(c46, 3);
  Real c53 = Power(c52, 2);
  Real c86 = Power(c54, 2);
  Real c89 = Power(c54, 3);
  Real c2 = Power(c1, -4);
  Real c25 = Power(c7, -4);
  Real c62 = Power(c37, -4);
  Real c73 = Power(c46, -4);
  Real c84 = Power(c54, -4);
  out1(0) = c2 * (2 * C0221 * c11 * c12 * c18 + 3 * C0231 * c1 * c17 * c18 +
                  4 * C0222 * c11 * c12 * c14 * c4 + C0211 * c18 * c6 +
                  2 * C0212 * c14 * c4 * c6 + 3 * C0213 * c5 * c6 * c8);
  out1(1) = 0;
  out1(2) =
      c25 *
      (2 * C22 * c30 + (2 * C0222 * c17 * c30) / Power(c1, 2) +
       2 * C0212 * c11 * c28 * c30 + (2 * C2322 * c30 * c36) / Power(c37, 2) +
       2 * C2321 * c30 * c35 * c40 + (2 * C2422 * c30 * c45) / Power(c46, 2) +
       2 * C2421 * c30 * c44 * c49 + 12 * C24 * c5 +
       (2 * C2522 * c30 * c53) / Power(c54, 2) + 2 * C2521 * c30 * c52 * c57 +
       6 * C23 * c4 * c7 + 6 * C0213 * c11 * c28 * c4 * c7 +
       6 * C2331 * c35 * c4 * c40 * c7 + 6 * C2431 * c4 * c44 * c49 * c7 +
       6 * C2531 * c4 * c52 * c57 * c7);
  out1(3) = c62 * (3 * C2313 * c18 * c36 * c37 + 2 * C2312 * c18 * c35 * c64 +
                   4 * C2322 * c14 * c35 * c4 * c64 + C2311 * c18 * c67 +
                   2 * C2321 * c14 * c4 * c67 + 3 * C2331 * c5 * c67 * c8);
  out1(4) = c73 * (3 * C2413 * c18 * c45 * c46 + 2 * C2412 * c18 * c44 * c75 +
                   4 * C2422 * c14 * c4 * c44 * c75 + C2411 * c18 * c78 +
                   2 * C2421 * c14 * c4 * c78 + 3 * C2431 * c5 * c78 * c8);
  out1(5) = c84 * (3 * C2513 * c18 * c53 * c54 + 2 * C2512 * c18 * c52 * c86 +
                   4 * C2522 * c14 * c4 * c52 * c86 + C2511 * c18 * c89 +
                   2 * C2521 * c14 * c4 * c89 + 3 * C2531 * c5 * c8 * c89);
  out2(0) = c2 * (4 * C0222 * c11 * c12 * c14 + 2 * C0212 * c14 * c6 +
                  6 * C0213 * c4 * c6 * c8);
  out2(1) = 0;
  out2(2) = c25 * (24 * C24 * c4 + 6 * C23 * c7 + 6 * C0213 * c11 * c28 * c7 +
                   6 * C2331 * c35 * c40 * c7 + 6 * C2431 * c44 * c49 * c7 +
                   6 * C2531 * c52 * c57 * c7);
  out2(3) = c62 * (4 * C2322 * c14 * c35 * c64 + 2 * C2321 * c14 * c67 +
                   6 * C2331 * c4 * c67 * c8);
  out2(4) = c73 * (4 * C2422 * c14 * c44 * c75 + 2 * C2421 * c14 * c78 +
                   6 * C2431 * c4 * c78 * c8);
  out2(5) = c84 * (4 * C2522 * c14 * c52 * c86 + 2 * C2521 * c14 * c89 +
                   6 * C2531 * c4 * c8 * c89);
  return std::make_pair(out1, out2);
}
std::pair<grad_type, grad_type>
FittedMaterial::grad_taylor_12_3(const strain_type &ek) {
  auto out1 = grad_type(0);
  auto out2 = grad_type(0);
  Real c1 = ekscale(0);
  Real c3 = ek(3);
  Real c6 = ekscale(3);
  Real c5 = Power(c1, 3);
  Real c12 = Power(c6, -2);
  Real c9 = ek(0);
  Real c10 = -1 + c9;
  Real c11 = Power(c1, 2);
  Real c17 = 1 / c6;
  Real c4 = Power(c3, 2);
  Real c24 = ekscale(2);
  Real c7 = Power(c6, -3);
  Real c26 = Power(c24, 3);
  Real c28 = ek(2);
  Real c29 = -1 + c28;
  Real c30 = Power(c24, 2);
  Real c16 = Power(c10, 2);
  Real c47 = Power(c6, 2);
  Real c43 = 1 / c1;
  Real c33 = Power(c29, 2);
  Real c45 = 1 / c24;
  Real c55 = ek(4);
  Real c57 = ekscale(4);
  Real c60 = 1 / c57;
  Real c63 = ek(5);
  Real c65 = ekscale(5);
  Real c68 = 1 / c65;
  Real c56 = Power(c55, 2);
  Real c75 = Power(c57, 2);
  Real c78 = Power(c57, 3);
  Real c64 = Power(c63, 2);
  Real c86 = Power(c65, 2);
  Real c89 = Power(c65, 3);
  Real c2 = Power(c1, -4);
  Real c25 = Power(c24, -4);
  Real c39 = Power(c6, -4);
  Real c73 = Power(c57, -4);
  Real c84 = Power(c65, -4);
  out1(0) = c2 * (2 * C0321 * c10 * c11 * c17 + 3 * C0331 * c1 * c16 * c17 +
                  4 * C0322 * c10 * c11 * c12 * c3 + C0311 * c17 * c5 +
                  2 * C0312 * c12 * c3 * c5 + 3 * C0313 * c4 * c5 * c7);
  out1(1) = 0;
  out1(2) =
      c25 * (C2311 * c17 * c26 + 2 * C2312 * c12 * c26 * c3 +
             2 * C2321 * c17 * c29 * c30 + 4 * C2322 * c12 * c29 * c3 * c30 +
             3 * C2331 * c17 * c24 * c33 + 3 * C2313 * c26 * c4 * c7);
  out1(3) =
      c39 *
      (12 * C34 * c4 + 2 * C32 * c47 + (2 * C0322 * c16 * c47) / Power(c1, 2) +
       (2 * C2322 * c33 * c47) / Power(c24, 2) + 2 * C0312 * c10 * c43 * c47 +
       2 * C2312 * c29 * c45 * c47 + (2 * C3422 * c47 * c56) / Power(c57, 2) +
       6 * C33 * c3 * c6 + 6 * C0313 * c10 * c3 * c43 * c6 +
       6 * C2313 * c29 * c3 * c45 * c6 + 2 * C3421 * c47 * c55 * c60 +
       6 * C3431 * c3 * c55 * c6 * c60 +
       (2 * C3522 * c47 * c64) / Power(c65, 2) + 2 * C3521 * c47 * c63 * c68 +
       6 * C3531 * c3 * c6 * c63 * c68);
  out1(4) = c73 * (3 * C3413 * c17 * c56 * c57 + 2 * C3412 * c17 * c55 * c75 +
                   4 * C3422 * c12 * c3 * c55 * c75 + C3411 * c17 * c78 +
                   2 * C3421 * c12 * c3 * c78 + 3 * C3431 * c4 * c7 * c78);
  out1(5) = c84 * (3 * C3513 * c17 * c64 * c65 + 2 * C3512 * c17 * c63 * c86 +
                   4 * C3522 * c12 * c3 * c63 * c86 + C3511 * c17 * c89 +
                   2 * C3521 * c12 * c3 * c89 + 3 * C3531 * c4 * c7 * c89);
  out2(0) = c2 * (4 * C0322 * c10 * c11 * c12 + 2 * C0312 * c12 * c5 +
                  6 * C0313 * c3 * c5 * c7);
  out2(1) = 0;
  out2(2) = c25 * (2 * C2312 * c12 * c26 + 4 * C2322 * c12 * c29 * c30 +
                   6 * C2313 * c26 * c3 * c7);
  out2(3) = c39 * (24 * C34 * c3 + 6 * C33 * c6 + 6 * C0313 * c10 * c43 * c6 +
                   6 * C2313 * c29 * c45 * c6 + 6 * C3431 * c55 * c6 * c60 +
                   6 * C3531 * c6 * c63 * c68);
  out2(4) = c73 * (4 * C3422 * c12 * c55 * c75 + 2 * C3421 * c12 * c78 +
                   6 * C3431 * c3 * c7 * c78);
  out2(5) = c84 * (4 * C3522 * c12 * c63 * c86 + 2 * C3521 * c12 * c89 +
                   6 * C3531 * c3 * c7 * c89);
  return std::make_pair(out1, out2);
}
std::pair<grad_type, grad_type>
FittedMaterial::grad_taylor_12_4(const strain_type &ek) {
  auto out1 = grad_type(0);
  auto out2 = grad_type(0);
  Real c1 = ekscale(0);
  Real c3 = ek(4);
  Real c6 = ekscale(4);
  Real c5 = Power(c1, 3);
  Real c12 = Power(c6, -2);
  Real c9 = ek(0);
  Real c10 = -1 + c9;
  Real c11 = Power(c1, 2);
  Real c17 = 1 / c6;
  Real c4 = Power(c3, 2);
  Real c24 = ekscale(2);
  Real c7 = Power(c6, -3);
  Real c26 = Power(c24, 3);
  Real c28 = ek(2);
  Real c29 = -1 + c28;
  Real c30 = Power(c24, 2);
  Real c39 = ekscale(3);
  Real c42 = Power(c39, 3);
  Real c44 = ek(3);
  Real c45 = Power(c39, 2);
  Real c16 = Power(c10, 2);
  Real c63 = Power(c6, 2);
  Real c57 = 1 / c1;
  Real c33 = Power(c29, 2);
  Real c59 = 1 / c24;
  Real c48 = Power(c44, 2);
  Real c61 = 1 / c39;
  Real c2 = Power(c1, -4);
  Real c25 = Power(c24, -4);
  Real c40 = Power(c39, -4);
  Real c54 = Power(c6, -4);
  out1(0) = c2 * (2 * C0421 * c10 * c11 * c17 + 3 * C0431 * c1 * c16 * c17 +
                  4 * C0422 * c10 * c11 * c12 * c3 + C0411 * c17 * c5 +
                  2 * C0412 * c12 * c3 * c5 + 3 * C0413 * c4 * c5 * c7);
  out1(1) = 0;
  out1(2) =
      c25 * (C2411 * c17 * c26 + 2 * C2412 * c12 * c26 * c3 +
             2 * C2421 * c17 * c29 * c30 + 4 * C2422 * c12 * c29 * c3 * c30 +
             3 * C2431 * c17 * c24 * c33 + 3 * C2413 * c26 * c4 * c7);
  out1(3) =
      c40 * (C3411 * c17 * c42 + 2 * C3412 * c12 * c3 * c42 +
             2 * C3421 * c17 * c44 * c45 + 4 * C3422 * c12 * c3 * c44 * c45 +
             3 * C3431 * c17 * c39 * c48 + 3 * C3413 * c4 * c42 * c7);
  out1(4) =
      c54 *
      (12 * C44 * c4 + 6 * C43 * c3 * c6 + 6 * C0413 * c10 * c3 * c57 * c6 +
       6 * C2413 * c29 * c3 * c59 * c6 + 6 * C3413 * c3 * c44 * c6 * c61 +
       2 * C42 * c63 + (2 * C0422 * c16 * c63) / Power(c1, 2) +
       (2 * C2422 * c33 * c63) / Power(c24, 2) +
       (2 * C3422 * c48 * c63) / Power(c39, 2) + 2 * C0412 * c10 * c57 * c63 +
       2 * C2412 * c29 * c59 * c63 + 2 * C3412 * c44 * c61 * c63);
  out1(5) = 0;
  out2(0) = c2 * (4 * C0422 * c10 * c11 * c12 + 2 * C0412 * c12 * c5 +
                  6 * C0413 * c3 * c5 * c7);
  out2(1) = 0;
  out2(2) = c25 * (2 * C2412 * c12 * c26 + 4 * C2422 * c12 * c29 * c30 +
                   6 * C2413 * c26 * c3 * c7);
  out2(3) = c40 * (2 * C3412 * c12 * c42 + 4 * C3422 * c12 * c44 * c45 +
                   6 * C3413 * c3 * c42 * c7);
  out2(4) = c54 * (24 * C44 * c3 + 6 * C43 * c6 + 6 * C0413 * c10 * c57 * c6 +
                   6 * C2413 * c29 * c59 * c6 + 6 * C3413 * c44 * c6 * c61);
  out2(5) = 0;
  return std::make_pair(out1, out2);
}
std::pair<grad_type, grad_type>
FittedMaterial::grad_taylor_12_5(const strain_type &ek) {
  auto out1 = grad_type(0);
  auto out2 = grad_type(0);
  Real c1 = ekscale(0);
  Real c3 = ek(5);
  Real c6 = ekscale(5);
  Real c5 = Power(c1, 3);
  Real c12 = Power(c6, -2);
  Real c9 = ek(0);
  Real c10 = -1 + c9;
  Real c11 = Power(c1, 2);
  Real c17 = 1 / c6;
  Real c4 = Power(c3, 2);
  Real c24 = ekscale(2);
  Real c7 = Power(c6, -3);
  Real c26 = Power(c24, 3);
  Real c28 = ek(2);
  Real c29 = -1 + c28;
  Real c30 = Power(c24, 2);
  Real c39 = ekscale(3);
  Real c42 = Power(c39, 3);
  Real c44 = ek(3);
  Real c45 = Power(c39, 2);
  Real c16 = Power(c10, 2);
  Real c63 = Power(c6, 2);
  Real c57 = 1 / c1;
  Real c33 = Power(c29, 2);
  Real c59 = 1 / c24;
  Real c48 = Power(c44, 2);
  Real c61 = 1 / c39;
  Real c2 = Power(c1, -4);
  Real c25 = Power(c24, -4);
  Real c40 = Power(c39, -4);
  Real c54 = Power(c6, -4);
  out1(0) = c2 * (2 * C0521 * c10 * c11 * c17 + 3 * C0531 * c1 * c16 * c17 +
                  4 * C0522 * c10 * c11 * c12 * c3 + C0511 * c17 * c5 +
                  2 * C0512 * c12 * c3 * c5 + 3 * C0513 * c4 * c5 * c7);
  out1(1) = 0;
  out1(2) =
      c25 * (C2511 * c17 * c26 + 2 * C2512 * c12 * c26 * c3 +
             2 * C2521 * c17 * c29 * c30 + 4 * C2522 * c12 * c29 * c3 * c30 +
             3 * C2531 * c17 * c24 * c33 + 3 * C2513 * c26 * c4 * c7);
  out1(3) =
      c40 * (C3511 * c17 * c42 + 2 * C3512 * c12 * c3 * c42 +
             2 * C3521 * c17 * c44 * c45 + 4 * C3522 * c12 * c3 * c44 * c45 +
             3 * C3531 * c17 * c39 * c48 + 3 * C3513 * c4 * c42 * c7);
  out1(4) = 0;
  out1(5) =
      c54 *
      (12 * C54 * c4 + 6 * C53 * c3 * c6 + 6 * C0513 * c10 * c3 * c57 * c6 +
       6 * C2513 * c29 * c3 * c59 * c6 + 6 * C3513 * c3 * c44 * c6 * c61 +
       2 * C52 * c63 + (2 * C0522 * c16 * c63) / Power(c1, 2) +
       (2 * C2522 * c33 * c63) / Power(c24, 2) +
       (2 * C3522 * c48 * c63) / Power(c39, 2) + 2 * C0512 * c10 * c57 * c63 +
       2 * C2512 * c29 * c59 * c63 + 2 * C3512 * c44 * c61 * c63);
  out2(0) = c2 * (4 * C0522 * c10 * c11 * c12 + 2 * C0512 * c12 * c5 +
                  6 * C0513 * c3 * c5 * c7);
  out2(1) = 0;
  out2(2) = c25 * (2 * C2512 * c12 * c26 + 4 * C2522 * c12 * c29 * c30 +
                   6 * C2513 * c26 * c3 * c7);
  out2(3) = c40 * (2 * C3512 * c12 * c42 + 4 * C3522 * c12 * c44 * c45 +
                   6 * C3513 * c3 * c42 * c7);
  out2(4) = 0;
  out2(5) = c54 * (24 * C54 * c3 + 6 * C53 * c6 + 6 * C0513 * c10 * c57 * c6 +
                   6 * C2513 * c29 * c59 * c6 + 6 * C3513 * c44 * c6 * c61);
  return std::make_pair(out1, out2);
}
// -------------- gradhesss -------------- //
std::pair<gradhess_type, gradhess_type>
FittedMaterial::gradhess_taylor_12_0(const strain_type &ek) {
  Mat6x6 hessD(0), hessDD(0);
  Vec6 gradD(0), gradDD(0);
  // MM order: {gradd, hessd, graddd, hessdd}
  auto out1 = [&](int i) -> Real & { return gradD[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hessD(i, j); };
  auto out3 = [&](int i) -> Real & { return gradDD[i]; };
  auto out4 = [&](int i, int j) -> Real & { return hessDD(i, j); };
  Real c10 = ek(0);
  Real c11 = -1 + c10;
  Real c1 = ekscale(0);
  Real c16 = Power(c1, 2);
  Real c18 = ek(1);
  Real c21 = ekscale(1);
  Real c24 = 1 / c21;
  Real c27 = ek(2);
  Real c28 = -1 + c27;
  Real c30 = ekscale(2);
  Real c33 = 1 / c30;
  Real c36 = ek(3);
  Real c38 = ekscale(3);
  Real c42 = 1 / c38;
  Real c45 = ek(4);
  Real c47 = ekscale(4);
  Real c50 = 1 / c47;
  Real c53 = ek(5);
  Real c55 = ekscale(5);
  Real c58 = 1 / c55;
  Real c19 = Power(c18, 2);
  Real c64 = 1 / c1;
  Real c67 = Power(c21, 2);
  Real c12 = Power(c11, 2);
  Real c66 = Power(c1, -2);
  Real c71 = Power(c21, 3);
  Real c29 = Power(c28, 2);
  Real c79 = Power(c30, 2);
  Real c70 = Power(c1, -3);
  Real c82 = Power(c30, 3);
  Real c37 = Power(c36, 2);
  Real c90 = Power(c38, 2);
  Real c93 = Power(c38, 3);
  Real c46 = Power(c45, 2);
  Real c101 = Power(c47, 2);
  Real c104 = Power(c47, 3);
  Real c54 = Power(c53, 2);
  Real c112 = Power(c55, 2);
  Real c115 = Power(c55, 3);
  Real c2 = Power(c1, -4);
  Real c22 = Power(c21, -2);
  Real c31 = Power(c30, -2);
  Real c39 = Power(c38, -2);
  Real c48 = Power(c47, -2);
  Real c56 = Power(c55, -2);
  Real c130 = 4 * C0122 * c1 * c18;
  Real c131 = 6 * C0131 * c11;
  Real c132 = 2 * C0121 * c1;
  Real c133 = c131 + c132;
  Real c134 = c133 * c21;
  Real c135 = c130 + c134;
  Real c136 = c135 * c22 * c70;
  Real c137 = 4 * C0222 * c1 * c28;
  Real c138 = 6 * C0231 * c11;
  Real c139 = 2 * C0221 * c1;
  Real c140 = c138 + c139;
  Real c141 = c140 * c30;
  Real c142 = c137 + c141;
  Real c143 = c142 * c31 * c70;
  Real c77 = Power(c30, -4);
  Real c144 = 4 * C0322 * c1 * c36;
  Real c145 = 6 * C0331 * c11;
  Real c146 = 2 * C0321 * c1;
  Real c147 = c145 + c146;
  Real c148 = c147 * c38;
  Real c149 = c144 + c148;
  Real c150 = c149 * c39 * c70;
  Real c88 = Power(c38, -4);
  Real c151 = 4 * C0422 * c1 * c45;
  Real c152 = 6 * C0431 * c11;
  Real c153 = 2 * C0421 * c1;
  Real c154 = c152 + c153;
  Real c155 = c154 * c47;
  Real c156 = c151 + c155;
  Real c157 = c156 * c48 * c70;
  Real c99 = Power(c47, -4);
  Real c158 = 4 * C0522 * c1 * c53;
  Real c159 = 6 * C0531 * c11;
  Real c160 = 2 * C0521 * c1;
  Real c161 = c159 + c160;
  Real c162 = c161 * c55;
  Real c163 = c158 + c162;
  Real c164 = c163 * c56 * c70;
  Real c110 = Power(c55, -4);
  Real c63 = Power(c21, -4);
  Real c228 = 6 * C0131 * c24 * c70;
  Real c229 = 6 * C0231 * c33 * c70;
  Real c230 = 6 * C0331 * c42 * c70;
  Real c231 = 6 * C0431 * c50 * c70;
  Real c232 = 6 * C0531 * c58 * c70;
  out1(0) =
      c2 * (6 * C03 * c1 * c11 + 12 * C04 * c12 + 2 * C02 * c16 +
            2 * C0122 * c16 * c19 * c22 + 6 * C0131 * c1 * c11 * c18 * c24 +
            2 * C0121 * c16 * c18 * c24 + 2 * C0222 * c16 * c29 * c31 +
            6 * C0231 * c1 * c11 * c28 * c33 + 2 * C0221 * c16 * c28 * c33 +
            2 * C0322 * c16 * c37 * c39 + 6 * C0331 * c1 * c11 * c36 * c42 +
            2 * C0321 * c16 * c36 * c42 + 2 * C0422 * c16 * c46 * c48 +
            6 * C0431 * c1 * c11 * c45 * c50 + 2 * C0421 * c16 * c45 * c50 +
            2 * C0522 * c16 * c54 * c56 + 6 * C0531 * c1 * c11 * c53 * c58 +
            2 * C0521 * c16 * c53 * c58);
  out1(1) = c63 * (3 * C0113 * c19 * c21 * c64 + 2 * C0112 * c18 * c64 * c67 +
                   4 * C0122 * c11 * c18 * c66 * c67 + C0111 * c64 * c71 +
                   2 * C0121 * c11 * c66 * c71 + 3 * C0131 * c12 * c70 * c71);
  out1(2) = c77 * (3 * C0213 * c29 * c30 * c64 + 2 * C0212 * c28 * c64 * c79 +
                   4 * C0222 * c11 * c28 * c66 * c79 + C0211 * c64 * c82 +
                   2 * C0221 * c11 * c66 * c82 + 3 * C0231 * c12 * c70 * c82);
  out1(3) = c88 * (3 * C0313 * c37 * c38 * c64 + 2 * C0312 * c36 * c64 * c90 +
                   4 * C0322 * c11 * c36 * c66 * c90 + C0311 * c64 * c93 +
                   2 * C0321 * c11 * c66 * c93 + 3 * C0331 * c12 * c70 * c93);
  out1(4) =
      (C0411 * c104 * c64 + 2 * C0412 * c101 * c45 * c64 +
       3 * C0413 * c46 * c47 * c64 + 2 * C0421 * c104 * c11 * c66 +
       4 * C0422 * c101 * c11 * c45 * c66 + 3 * C0431 * c104 * c12 * c70) *
      c99;
  out1(5) = c110 *
            (C0511 * c115 * c64 + 2 * C0512 * c112 * c53 * c64 +
             3 * C0513 * c54 * c55 * c64 + 2 * C0521 * c11 * c115 * c66 +
             4 * C0522 * c11 * c112 * c53 * c66 + 3 * C0531 * c115 * c12 * c70);
  out2(0, 0) = 2 * c2 *
               (3 * C03 * c1 + 12 * C04 * c11 + 3 * C0131 * c1 * c18 * c24 +
                3 * C0231 * c1 * c28 * c33 + 3 * C0331 * c1 * c36 * c42 +
                3 * C0431 * c1 * c45 * c50 + 3 * C0531 * c1 * c53 * c58);
  out2(0, 1) = c136;
  out2(0, 2) = c143;
  out2(0, 3) = c150;
  out2(0, 4) = c157;
  out2(0, 5) = c164;
  out2(1, 0) = c136;
  out2(1, 1) =
      (2 * (3 * C0113 * c1 * c18 + (C0112 * c1 + 2 * C0122 * c11) * c21) *
       c66) /
      Power(c21, 3);
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = c143;
  out2(2, 1) = 0;
  out2(2, 2) = 2 * c77 *
               (3 * C0213 * c28 * c30 * c64 + C0212 * c64 * c79 +
                2 * C0222 * c11 * c66 * c79);
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = c150;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = 2 * c88 *
               (3 * C0313 * c36 * c38 * c64 + C0312 * c64 * c90 +
                2 * C0322 * c11 * c66 * c90);
  out2(3, 4) = 0;
  out2(3, 5) = 0;
  out2(4, 0) = c157;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = 0;
  out2(4, 4) = 2 *
               (C0412 * c101 * c64 + 3 * C0413 * c45 * c47 * c64 +
                2 * C0422 * c101 * c11 * c66) *
               c99;
  out2(4, 5) = 0;
  out2(5, 0) = c164;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) = 2 * c110 *
               (C0512 * c112 * c64 + 3 * C0513 * c53 * c55 * c64 +
                2 * C0522 * c11 * c112 * c66);
  out3(0) = c2 * (6 * C03 * c1 + 24 * C04 * c11 + 6 * C0131 * c1 * c18 * c24 +
                  6 * C0231 * c1 * c28 * c33 + 6 * C0331 * c1 * c36 * c42 +
                  6 * C0431 * c1 * c45 * c50 + 6 * C0531 * c1 * c53 * c58);
  out3(1) = c63 * (4 * C0122 * c18 * c66 * c67 + 2 * C0121 * c66 * c71 +
                   6 * C0131 * c11 * c70 * c71);
  out3(2) = c77 * (4 * C0222 * c28 * c66 * c79 + 2 * C0221 * c66 * c82 +
                   6 * C0231 * c11 * c70 * c82);
  out3(3) = c88 * (4 * C0322 * c36 * c66 * c90 + 2 * C0321 * c66 * c93 +
                   6 * C0331 * c11 * c70 * c93);
  out3(4) = (2 * C0421 * c104 * c66 + 4 * C0422 * c101 * c45 * c66 +
             6 * C0431 * c104 * c11 * c70) *
            c99;
  out3(5) = c110 * (2 * C0521 * c115 * c66 + 4 * C0522 * c112 * c53 * c66 +
                    6 * C0531 * c11 * c115 * c70);
  out4(0, 0) = 24 * C04 * c2;
  out4(0, 1) = c228;
  out4(0, 2) = c229;
  out4(0, 3) = c230;
  out4(0, 4) = c231;
  out4(0, 5) = c232;
  out4(1, 0) = c228;
  out4(1, 1) = 4 * C0122 * c22 * c66;
  out4(1, 2) = 0;
  out4(1, 3) = 0;
  out4(1, 4) = 0;
  out4(1, 5) = 0;
  out4(2, 0) = c229;
  out4(2, 1) = 0;
  out4(2, 2) = 4 * C0222 * c31 * c66;
  out4(2, 3) = 0;
  out4(2, 4) = 0;
  out4(2, 5) = 0;
  out4(3, 0) = c230;
  out4(3, 1) = 0;
  out4(3, 2) = 0;
  out4(3, 3) = 4 * C0322 * c39 * c66;
  out4(3, 4) = 0;
  out4(3, 5) = 0;
  out4(4, 0) = c231;
  out4(4, 1) = 0;
  out4(4, 2) = 0;
  out4(4, 3) = 0;
  out4(4, 4) = 4 * C0422 * c48 * c66;
  out4(4, 5) = 0;
  out4(5, 0) = c232;
  out4(5, 1) = 0;
  out4(5, 2) = 0;
  out4(5, 3) = 0;
  out4(5, 4) = 0;
  out4(5, 5) = 4 * C0522 * c56 * c66;
  return std::make_pair(std::make_pair(hessD, gradD),
                        std::make_pair(hessDD, gradDD));
}
std::pair<gradhess_type, gradhess_type>
FittedMaterial::gradhess_taylor_12_1(const strain_type &ek) {
  Mat6x6 hessD(0), hessDD(0);
  Vec6 gradD(0), gradDD(0);
  // MM order: {gradd, hessd, graddd, hessdd}
  auto out1 = [&](int i) -> Real & { return gradD[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hessD(i, j); };
  auto out3 = [&](int i) -> Real & { return gradDD[i]; };
  auto out4 = [&](int i, int j) -> Real & { return hessDD(i, j); };
  Real c1 = ekscale(0);
  Real c10 = ek(1);
  Real c14 = ekscale(1);
  Real c12 = Power(c1, 3);
  Real c21 = Power(c14, -2);
  Real c17 = ek(0);
  Real c18 = -1 + c17;
  Real c19 = Power(c1, 2);
  Real c25 = 1 / c14;
  Real c11 = Power(c10, 2);
  Real c24 = Power(c18, 2);
  Real c36 = Power(c14, 2);
  Real c34 = 1 / c1;
  Real c2 = Power(c1, -4);
  Real c15 = Power(c14, -3);
  Real c49 = Power(c1, -3);
  Real c50 = 6 * C0113 * c10 * c19;
  Real c51 = 4 * C0122 * c1 * c18;
  Real c52 = 2 * C0112 * c19;
  Real c53 = c51 + c52;
  Real c54 = c14 * c53;
  Real c55 = c50 + c54;
  Real c56 = c15 * c49 * c55;
  Real c31 = Power(c14, -4);
  Real c38 = Power(c1, -2);
  Real c75 = 6 * C0113 * c15 * c34;
  out1(0) = c2 * (3 * C0113 * c11 * c12 * c15 + 2 * C0112 * c10 * c12 * c21 +
                  4 * C0122 * c10 * c18 * c19 * c21 + C0111 * c12 * c25 +
                  2 * C0121 * c18 * c19 * c25 + 3 * C0131 * c1 * c24 * c25);
  out1(1) = c31 * (12 * C14 * c11 + 6 * C13 * c10 * c14 +
                   6 * C0113 * c10 * c14 * c18 * c34 + 2 * C12 * c36 +
                   2 * C0112 * c18 * c34 * c36 + 2 * C0122 * c24 * c36 * c38);
  out1(2) = 0;
  out1(3) = 0;
  out1(4) = 0;
  out1(5) = 0;
  out2(0, 0) = 2 * c2 *
               (2 * C0122 * c10 * c19 * c21 + 3 * C0131 * c1 * c18 * c25 +
                C0121 * c19 * c25);
  out2(0, 1) = c56;
  out2(0, 2) = 0;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = c56;
  out2(1, 1) =
      2 * c31 *
      (12 * C14 * c10 + c14 * (3 * C0113 * c1 * c18 + 3 * C13 * c19) * c38);
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) = 0;
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = 0;
  out2(3, 4) = 0;
  out2(3, 5) = 0;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = 0;
  out2(4, 4) = 0;
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) = 0;
  out3(0) = c2 * (6 * C0113 * c10 * c12 * c15 + 2 * C0112 * c12 * c21 +
                  4 * C0122 * c18 * c19 * c21);
  out3(1) =
      c31 * (24 * C14 * c10 + 6 * C13 * c14 + 6 * C0113 * c14 * c18 * c34);
  out3(2) = 0;
  out3(3) = 0;
  out3(4) = 0;
  out3(5) = 0;
  out4(0, 0) = 4 * C0122 * c21 * c38;
  out4(0, 1) = c75;
  out4(0, 2) = 0;
  out4(0, 3) = 0;
  out4(0, 4) = 0;
  out4(0, 5) = 0;
  out4(1, 0) = c75;
  out4(1, 1) = 24 * C14 * c31;
  out4(1, 2) = 0;
  out4(1, 3) = 0;
  out4(1, 4) = 0;
  out4(1, 5) = 0;
  out4(2, 0) = 0;
  out4(2, 1) = 0;
  out4(2, 2) = 0;
  out4(2, 3) = 0;
  out4(2, 4) = 0;
  out4(2, 5) = 0;
  out4(3, 0) = 0;
  out4(3, 1) = 0;
  out4(3, 2) = 0;
  out4(3, 3) = 0;
  out4(3, 4) = 0;
  out4(3, 5) = 0;
  out4(4, 0) = 0;
  out4(4, 1) = 0;
  out4(4, 2) = 0;
  out4(4, 3) = 0;
  out4(4, 4) = 0;
  out4(4, 5) = 0;
  out4(5, 0) = 0;
  out4(5, 1) = 0;
  out4(5, 2) = 0;
  out4(5, 3) = 0;
  out4(5, 4) = 0;
  out4(5, 5) = 0;
  return std::make_pair(std::make_pair(hessD, gradD),
                        std::make_pair(hessDD, gradDD));
}
std::pair<gradhess_type, gradhess_type>
FittedMaterial::gradhess_taylor_12_2(const strain_type &ek) {
  Mat6x6 hessD(0), hessDD(0);
  Vec6 gradD(0), gradDD(0);
  // MM order: {gradd, hessd, graddd, hessdd}
  auto out1 = [&](int i) -> Real & { return gradD[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hessD(i, j); };
  auto out3 = [&](int i) -> Real & { return gradDD[i]; };
  auto out4 = [&](int i, int j) -> Real & { return hessDD(i, j); };
  Real c1 = ekscale(0);
  Real c10 = ek(2);
  Real c11 = -1 + c10;
  Real c15 = ekscale(2);
  Real c14 = Power(c1, 3);
  Real c22 = Power(c15, -2);
  Real c18 = ek(0);
  Real c19 = -1 + c18;
  Real c21 = Power(c1, 2);
  Real c26 = 1 / c15;
  Real c12 = Power(c11, 2);
  Real c25 = Power(c19, 2);
  Real c37 = Power(c15, 2);
  Real c35 = 1 / c1;
  Real c43 = ek(3);
  Real c45 = ekscale(3);
  Real c48 = 1 / c45;
  Real c51 = ek(4);
  Real c53 = ekscale(4);
  Real c56 = 1 / c53;
  Real c59 = ek(5);
  Real c61 = ekscale(5);
  Real c64 = 1 / c61;
  Real c44 = Power(c43, 2);
  Real c71 = Power(c45, 2);
  Real c16 = Power(c15, -3);
  Real c74 = Power(c45, 3);
  Real c52 = Power(c51, 2);
  Real c82 = Power(c53, 2);
  Real c85 = Power(c53, 3);
  Real c60 = Power(c59, 2);
  Real c93 = Power(c61, 2);
  Real c96 = Power(c61, 3);
  Real c2 = Power(c1, -4);
  Real c107 = Power(c1, -3);
  Real c108 = 6 * C0213 * c11 * c21;
  Real c109 = 4 * C0222 * c1 * c19;
  Real c110 = 2 * C0212 * c21;
  Real c111 = c109 + c110;
  Real c112 = c111 * c15;
  Real c113 = c108 + c112;
  Real c114 = c107 * c113 * c16;
  Real c32 = Power(c15, -4);
  Real c46 = Power(c45, -2);
  Real c54 = Power(c53, -2);
  Real c62 = Power(c61, -2);
  Real c123 = 4 * C2322 * c15 * c43;
  Real c124 = 6 * C2331 * c11;
  Real c125 = 2 * C2321 * c15;
  Real c126 = c124 + c125;
  Real c127 = c126 * c45;
  Real c128 = c123 + c127;
  Real c129 = c128 * c16 * c46;
  Real c69 = Power(c45, -4);
  Real c130 = 4 * C2422 * c15 * c51;
  Real c131 = 6 * C2431 * c11;
  Real c132 = 2 * C2421 * c15;
  Real c133 = c131 + c132;
  Real c134 = c133 * c53;
  Real c135 = c130 + c134;
  Real c136 = c135 * c16 * c54;
  Real c80 = Power(c53, -4);
  Real c137 = 4 * C2522 * c15 * c59;
  Real c138 = 6 * C2531 * c11;
  Real c139 = 2 * C2521 * c15;
  Real c140 = c138 + c139;
  Real c141 = c140 * c61;
  Real c142 = c137 + c141;
  Real c143 = c142 * c16 * c62;
  Real c91 = Power(c61, -4);
  Real c39 = Power(c1, -2);
  Real c188 = 6 * C0213 * c16 * c35;
  Real c190 = 6 * C2331 * c16 * c48;
  Real c191 = 6 * C2431 * c16 * c56;
  Real c192 = 6 * C2531 * c16 * c64;
  out1(0) = c2 * (3 * C0213 * c12 * c14 * c16 + 2 * C0212 * c11 * c14 * c22 +
                  4 * C0222 * c11 * c19 * c21 * c22 + C0211 * c14 * c26 +
                  2 * C0221 * c19 * c21 * c26 + 3 * C0231 * c1 * c25 * c26);
  out1(1) = 0;
  out1(2) =
      c32 * (12 * C24 * c12 + 6 * C23 * c11 * c15 +
             6 * C0213 * c11 * c15 * c19 * c35 + 2 * C22 * c37 +
             2 * C0212 * c19 * c35 * c37 + 2 * C0222 * c25 * c37 * c39 +
             2 * C2322 * c37 * c44 * c46 + 6 * C2331 * c11 * c15 * c43 * c48 +
             2 * C2321 * c37 * c43 * c48 + 2 * C2422 * c37 * c52 * c54 +
             6 * C2431 * c11 * c15 * c51 * c56 + 2 * C2421 * c37 * c51 * c56 +
             2 * C2522 * c37 * c60 * c62 + 6 * C2531 * c11 * c15 * c59 * c64 +
             2 * C2521 * c37 * c59 * c64);
  out1(3) =
      c69 * (3 * C2313 * c26 * c44 * c45 + 4 * C2322 * c11 * c22 * c43 * c71 +
             2 * C2312 * c26 * c43 * c71 + 3 * C2331 * c12 * c16 * c74 +
             2 * C2321 * c11 * c22 * c74 + C2311 * c26 * c74);
  out1(4) =
      c80 * (3 * C2413 * c26 * c52 * c53 + 4 * C2422 * c11 * c22 * c51 * c82 +
             2 * C2412 * c26 * c51 * c82 + 3 * C2431 * c12 * c16 * c85 +
             2 * C2421 * c11 * c22 * c85 + C2411 * c26 * c85);
  out1(5) =
      c91 * (3 * C2513 * c26 * c60 * c61 + 4 * C2522 * c11 * c22 * c59 * c93 +
             2 * C2512 * c26 * c59 * c93 + 3 * C2531 * c12 * c16 * c96 +
             2 * C2521 * c11 * c22 * c96 + C2511 * c26 * c96);
  out2(0, 0) = 2 * c2 *
               (2 * C0222 * c11 * c21 * c22 + 3 * C0231 * c1 * c19 * c26 +
                C0221 * c21 * c26);
  out2(0, 1) = 0;
  out2(0, 2) = c114;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) = 0;
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = c114;
  out2(2, 1) = 0;
  out2(2, 2) = 2 * c32 *
               (12 * C24 * c11 + 3 * C23 * c15 + 3 * C0213 * c15 * c19 * c35 +
                3 * C2331 * c15 * c43 * c48 + 3 * C2431 * c15 * c51 * c56 +
                3 * C2531 * c15 * c59 * c64);
  out2(2, 3) = c129;
  out2(2, 4) = c136;
  out2(2, 5) = c143;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = c129;
  out2(3, 3) = 2 * c69 *
               (3 * C2313 * c26 * c43 * c45 + 2 * C2322 * c11 * c22 * c71 +
                C2312 * c26 * c71);
  out2(3, 4) = 0;
  out2(3, 5) = 0;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = c136;
  out2(4, 3) = 0;
  out2(4, 4) = 2 * c80 *
               (3 * C2413 * c26 * c51 * c53 + 2 * C2422 * c11 * c22 * c82 +
                C2412 * c26 * c82);
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = c143;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) = 2 * c91 *
               (3 * C2513 * c26 * c59 * c61 + 2 * C2522 * c11 * c22 * c93 +
                C2512 * c26 * c93);
  out3(0) = c2 * (6 * C0213 * c11 * c14 * c16 + 2 * C0212 * c14 * c22 +
                  4 * C0222 * c19 * c21 * c22);
  out3(1) = 0;
  out3(2) = c32 * (24 * C24 * c11 + 6 * C23 * c15 +
                   6 * C0213 * c15 * c19 * c35 + 6 * C2331 * c15 * c43 * c48 +
                   6 * C2431 * c15 * c51 * c56 + 6 * C2531 * c15 * c59 * c64);
  out3(3) = c69 * (4 * C2322 * c22 * c43 * c71 + 6 * C2331 * c11 * c16 * c74 +
                   2 * C2321 * c22 * c74);
  out3(4) = c80 * (4 * C2422 * c22 * c51 * c82 + 6 * C2431 * c11 * c16 * c85 +
                   2 * C2421 * c22 * c85);
  out3(5) = c91 * (4 * C2522 * c22 * c59 * c93 + 6 * C2531 * c11 * c16 * c96 +
                   2 * C2521 * c22 * c96);
  out4(0, 0) = 4 * C0222 * c22 * c39;
  out4(0, 1) = 0;
  out4(0, 2) = c188;
  out4(0, 3) = 0;
  out4(0, 4) = 0;
  out4(0, 5) = 0;
  out4(1, 0) = 0;
  out4(1, 1) = 0;
  out4(1, 2) = 0;
  out4(1, 3) = 0;
  out4(1, 4) = 0;
  out4(1, 5) = 0;
  out4(2, 0) = c188;
  out4(2, 1) = 0;
  out4(2, 2) = 24 * C24 * c32;
  out4(2, 3) = c190;
  out4(2, 4) = c191;
  out4(2, 5) = c192;
  out4(3, 0) = 0;
  out4(3, 1) = 0;
  out4(3, 2) = c190;
  out4(3, 3) = 4 * C2322 * c22 * c46;
  out4(3, 4) = 0;
  out4(3, 5) = 0;
  out4(4, 0) = 0;
  out4(4, 1) = 0;
  out4(4, 2) = c191;
  out4(4, 3) = 0;
  out4(4, 4) = 4 * C2422 * c22 * c54;
  out4(4, 5) = 0;
  out4(5, 0) = 0;
  out4(5, 1) = 0;
  out4(5, 2) = c192;
  out4(5, 3) = 0;
  out4(5, 4) = 0;
  out4(5, 5) = 4 * C2522 * c22 * c62;
  return std::make_pair(std::make_pair(hessD, gradD),
                        std::make_pair(hessDD, gradDD));
}
std::pair<gradhess_type, gradhess_type>
FittedMaterial::gradhess_taylor_12_3(const strain_type &ek) {
  Mat6x6 hessD(0), hessDD(0);
  Vec6 gradD(0), gradDD(0);
  // MM order: {gradd, hessd, graddd, hessdd}
  auto out1 = [&](int i) -> Real & { return gradD[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hessD(i, j); };
  auto out3 = [&](int i) -> Real & { return gradDD[i]; };
  auto out4 = [&](int i, int j) -> Real & { return hessDD(i, j); };
  Real c1 = ekscale(0);
  Real c10 = ek(3);
  Real c14 = ekscale(3);
  Real c12 = Power(c1, 3);
  Real c21 = Power(c14, -2);
  Real c17 = ek(0);
  Real c18 = -1 + c17;
  Real c19 = Power(c1, 2);
  Real c25 = 1 / c14;
  Real c11 = Power(c10, 2);
  Real c31 = ekscale(2);
  Real c15 = Power(c14, -3);
  Real c33 = Power(c31, 3);
  Real c35 = ek(2);
  Real c36 = -1 + c35;
  Real c37 = Power(c31, 2);
  Real c24 = Power(c18, 2);
  Real c54 = Power(c14, 2);
  Real c50 = 1 / c1;
  Real c40 = Power(c36, 2);
  Real c52 = 1 / c31;
  Real c62 = ek(4);
  Real c64 = ekscale(4);
  Real c67 = 1 / c64;
  Real c70 = ek(5);
  Real c72 = ekscale(5);
  Real c75 = 1 / c72;
  Real c63 = Power(c62, 2);
  Real c82 = Power(c64, 2);
  Real c85 = Power(c64, 3);
  Real c71 = Power(c70, 2);
  Real c93 = Power(c72, 2);
  Real c96 = Power(c72, 3);
  Real c2 = Power(c1, -4);
  Real c32 = Power(c31, -4);
  Real c107 = Power(c1, -3);
  Real c108 = 6 * C0313 * c10 * c19;
  Real c109 = 4 * C0322 * c1 * c18;
  Real c110 = 2 * C0312 * c19;
  Real c111 = c109 + c110;
  Real c112 = c111 * c14;
  Real c113 = c108 + c112;
  Real c114 = c107 * c113 * c15;
  Real c120 = Power(c31, -3);
  Real c121 = 6 * C2313 * c10 * c37;
  Real c122 = 4 * C2322 * c31 * c36;
  Real c123 = 2 * C2312 * c37;
  Real c124 = c122 + c123;
  Real c125 = c124 * c14;
  Real c126 = c121 + c125;
  Real c127 = c120 * c126 * c15;
  Real c47 = Power(c14, -4);
  Real c65 = Power(c64, -2);
  Real c73 = Power(c72, -2);
  Real c136 = 4 * C3422 * c14 * c62;
  Real c137 = 6 * C3431 * c10;
  Real c138 = 2 * C3421 * c14;
  Real c139 = c137 + c138;
  Real c140 = c139 * c64;
  Real c141 = c136 + c140;
  Real c142 = c141 * c15 * c65;
  Real c80 = Power(c64, -4);
  Real c143 = 4 * C3522 * c14 * c70;
  Real c144 = 6 * C3531 * c10;
  Real c145 = 2 * C3521 * c14;
  Real c146 = c144 + c145;
  Real c147 = c146 * c72;
  Real c148 = c143 + c147;
  Real c149 = c148 * c15 * c73;
  Real c91 = Power(c72, -4);
  Real c56 = Power(c1, -2);
  Real c59 = Power(c31, -2);
  Real c189 = 6 * C0313 * c15 * c50;
  Real c191 = 6 * C2313 * c15 * c52;
  Real c193 = 6 * C3431 * c15 * c67;
  Real c194 = 6 * C3531 * c15 * c75;
  out1(0) = c2 * (3 * C0313 * c11 * c12 * c15 + 2 * C0312 * c10 * c12 * c21 +
                  4 * C0322 * c10 * c18 * c19 * c21 + C0311 * c12 * c25 +
                  2 * C0321 * c18 * c19 * c25 + 3 * C0331 * c1 * c24 * c25);
  out1(1) = 0;
  out1(2) = c32 * (3 * C2313 * c11 * c15 * c33 + 2 * C2312 * c10 * c21 * c33 +
                   C2311 * c25 * c33 + 4 * C2322 * c10 * c21 * c36 * c37 +
                   2 * C2321 * c25 * c36 * c37 + 3 * C2331 * c25 * c31 * c40);
  out1(3) =
      c47 * (12 * C34 * c11 + 6 * C33 * c10 * c14 +
             6 * C0313 * c10 * c14 * c18 * c50 +
             6 * C2313 * c10 * c14 * c36 * c52 + 2 * C32 * c54 +
             2 * C0312 * c18 * c50 * c54 + 2 * C2312 * c36 * c52 * c54 +
             2 * C0322 * c24 * c54 * c56 + 2 * C2322 * c40 * c54 * c59 +
             2 * C3422 * c54 * c63 * c65 + 6 * C3431 * c10 * c14 * c62 * c67 +
             2 * C3421 * c54 * c62 * c67 + 2 * C3522 * c54 * c71 * c73 +
             6 * C3531 * c10 * c14 * c70 * c75 + 2 * C3521 * c54 * c70 * c75);
  out1(4) =
      c80 * (3 * C3413 * c25 * c63 * c64 + 4 * C3422 * c10 * c21 * c62 * c82 +
             2 * C3412 * c25 * c62 * c82 + 3 * C3431 * c11 * c15 * c85 +
             2 * C3421 * c10 * c21 * c85 + C3411 * c25 * c85);
  out1(5) =
      c91 * (3 * C3513 * c25 * c71 * c72 + 4 * C3522 * c10 * c21 * c70 * c93 +
             2 * C3512 * c25 * c70 * c93 + 3 * C3531 * c11 * c15 * c96 +
             2 * C3521 * c10 * c21 * c96 + C3511 * c25 * c96);
  out2(0, 0) = 2 * c2 *
               (2 * C0322 * c10 * c19 * c21 + 3 * C0331 * c1 * c18 * c25 +
                C0321 * c19 * c25);
  out2(0, 1) = 0;
  out2(0, 2) = 0;
  out2(0, 3) = c114;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) = 0;
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) = 2 * c32 *
               (3 * C2331 * c25 * c31 * c36 + 2 * C2322 * c10 * c21 * c37 +
                C2321 * c25 * c37);
  out2(2, 3) = c127;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = c114;
  out2(3, 1) = 0;
  out2(3, 2) = c127;
  out2(3, 3) = 2 * c47 *
               (12 * C34 * c10 + 3 * C33 * c14 + 3 * C0313 * c14 * c18 * c50 +
                3 * C2313 * c14 * c36 * c52 + 3 * C3431 * c14 * c62 * c67 +
                3 * C3531 * c14 * c70 * c75);
  out2(3, 4) = c142;
  out2(3, 5) = c149;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = c142;
  out2(4, 4) = 2 * c80 *
               (3 * C3413 * c25 * c62 * c64 + 2 * C3422 * c10 * c21 * c82 +
                C3412 * c25 * c82);
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = c149;
  out2(5, 4) = 0;
  out2(5, 5) = 2 * c91 *
               (3 * C3513 * c25 * c70 * c72 + 2 * C3522 * c10 * c21 * c93 +
                C3512 * c25 * c93);
  out3(0) = c2 * (6 * C0313 * c10 * c12 * c15 + 2 * C0312 * c12 * c21 +
                  4 * C0322 * c18 * c19 * c21);
  out3(1) = 0;
  out3(2) = c32 * (6 * C2313 * c10 * c15 * c33 + 2 * C2312 * c21 * c33 +
                   4 * C2322 * c21 * c36 * c37);
  out3(3) = c47 * (24 * C34 * c10 + 6 * C33 * c14 +
                   6 * C0313 * c14 * c18 * c50 + 6 * C2313 * c14 * c36 * c52 +
                   6 * C3431 * c14 * c62 * c67 + 6 * C3531 * c14 * c70 * c75);
  out3(4) = c80 * (4 * C3422 * c21 * c62 * c82 + 6 * C3431 * c10 * c15 * c85 +
                   2 * C3421 * c21 * c85);
  out3(5) = c91 * (4 * C3522 * c21 * c70 * c93 + 6 * C3531 * c10 * c15 * c96 +
                   2 * C3521 * c21 * c96);
  out4(0, 0) = 4 * C0322 * c21 * c56;
  out4(0, 1) = 0;
  out4(0, 2) = 0;
  out4(0, 3) = c189;
  out4(0, 4) = 0;
  out4(0, 5) = 0;
  out4(1, 0) = 0;
  out4(1, 1) = 0;
  out4(1, 2) = 0;
  out4(1, 3) = 0;
  out4(1, 4) = 0;
  out4(1, 5) = 0;
  out4(2, 0) = 0;
  out4(2, 1) = 0;
  out4(2, 2) = 4 * C2322 * c21 * c59;
  out4(2, 3) = c191;
  out4(2, 4) = 0;
  out4(2, 5) = 0;
  out4(3, 0) = c189;
  out4(3, 1) = 0;
  out4(3, 2) = c191;
  out4(3, 3) = 24 * C34 * c47;
  out4(3, 4) = c193;
  out4(3, 5) = c194;
  out4(4, 0) = 0;
  out4(4, 1) = 0;
  out4(4, 2) = 0;
  out4(4, 3) = c193;
  out4(4, 4) = 4 * C3422 * c21 * c65;
  out4(4, 5) = 0;
  out4(5, 0) = 0;
  out4(5, 1) = 0;
  out4(5, 2) = 0;
  out4(5, 3) = c194;
  out4(5, 4) = 0;
  out4(5, 5) = 4 * C3522 * c21 * c73;
  return std::make_pair(std::make_pair(hessD, gradD),
                        std::make_pair(hessDD, gradDD));
}
std::pair<gradhess_type, gradhess_type>
FittedMaterial::gradhess_taylor_12_4(const strain_type &ek) {
  Mat6x6 hessD(0), hessDD(0);
  Vec6 gradD(0), gradDD(0);
  // MM order: {gradd, hessd, graddd, hessdd}
  auto out1 = [&](int i) -> Real & { return gradD[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hessD(i, j); };
  auto out3 = [&](int i) -> Real & { return gradDD[i]; };
  auto out4 = [&](int i, int j) -> Real & { return hessDD(i, j); };
  Real c1 = ekscale(0);
  Real c10 = ek(4);
  Real c14 = ekscale(4);
  Real c12 = Power(c1, 3);
  Real c21 = Power(c14, -2);
  Real c17 = ek(0);
  Real c18 = -1 + c17;
  Real c19 = Power(c1, 2);
  Real c25 = 1 / c14;
  Real c11 = Power(c10, 2);
  Real c31 = ekscale(2);
  Real c15 = Power(c14, -3);
  Real c33 = Power(c31, 3);
  Real c35 = ek(2);
  Real c36 = -1 + c35;
  Real c37 = Power(c31, 2);
  Real c47 = ekscale(3);
  Real c49 = Power(c47, 3);
  Real c51 = ek(3);
  Real c52 = Power(c47, 2);
  Real c24 = Power(c18, 2);
  Real c70 = Power(c14, 2);
  Real c64 = 1 / c1;
  Real c40 = Power(c36, 2);
  Real c66 = 1 / c31;
  Real c55 = Power(c51, 2);
  Real c68 = 1 / c47;
  Real c2 = Power(c1, -4);
  Real c32 = Power(c31, -4);
  Real c48 = Power(c47, -4);
  Real c88 = Power(c1, -3);
  Real c89 = 6 * C0413 * c10 * c19;
  Real c90 = 4 * C0422 * c1 * c18;
  Real c91 = 2 * C0412 * c19;
  Real c92 = c90 + c91;
  Real c93 = c14 * c92;
  Real c94 = c89 + c93;
  Real c95 = c15 * c88 * c94;
  Real c101 = Power(c31, -3);
  Real c102 = 6 * C2413 * c10 * c37;
  Real c103 = 4 * C2422 * c31 * c36;
  Real c104 = 2 * C2412 * c37;
  Real c105 = c103 + c104;
  Real c106 = c105 * c14;
  Real c107 = c102 + c106;
  Real c108 = c101 * c107 * c15;
  Real c114 = Power(c47, -3);
  Real c115 = 6 * C3413 * c10 * c52;
  Real c116 = 4 * C3422 * c47 * c51;
  Real c117 = 2 * C3412 * c52;
  Real c118 = c116 + c117;
  Real c119 = c118 * c14;
  Real c120 = c115 + c119;
  Real c121 = c114 * c120 * c15;
  Real c61 = Power(c14, -4);
  Real c72 = Power(c1, -2);
  Real c75 = Power(c31, -2);
  Real c78 = Power(c47, -2);
  Real c152 = 6 * C0413 * c15 * c64;
  Real c154 = 6 * C2413 * c15 * c66;
  Real c156 = 6 * C3413 * c15 * c68;
  out1(0) = c2 * (3 * C0413 * c11 * c12 * c15 + 2 * C0412 * c10 * c12 * c21 +
                  4 * C0422 * c10 * c18 * c19 * c21 + C0411 * c12 * c25 +
                  2 * C0421 * c18 * c19 * c25 + 3 * C0431 * c1 * c24 * c25);
  out1(1) = 0;
  out1(2) = c32 * (3 * C2413 * c11 * c15 * c33 + 2 * C2412 * c10 * c21 * c33 +
                   C2411 * c25 * c33 + 4 * C2422 * c10 * c21 * c36 * c37 +
                   2 * C2421 * c25 * c36 * c37 + 3 * C2431 * c25 * c31 * c40);
  out1(3) = c48 * (3 * C3413 * c11 * c15 * c49 + 2 * C3412 * c10 * c21 * c49 +
                   C3411 * c25 * c49 + 4 * C3422 * c10 * c21 * c51 * c52 +
                   2 * C3421 * c25 * c51 * c52 + 3 * C3431 * c25 * c47 * c55);
  out1(4) = c61 * (12 * C44 * c11 + 6 * C43 * c10 * c14 +
                   6 * C0413 * c10 * c14 * c18 * c64 +
                   6 * C2413 * c10 * c14 * c36 * c66 +
                   6 * C3413 * c10 * c14 * c51 * c68 + 2 * C42 * c70 +
                   2 * C0412 * c18 * c64 * c70 + 2 * C2412 * c36 * c66 * c70 +
                   2 * C3412 * c51 * c68 * c70 + 2 * C0422 * c24 * c70 * c72 +
                   2 * C2422 * c40 * c70 * c75 + 2 * C3422 * c55 * c70 * c78);
  out1(5) = 0;
  out2(0, 0) = 2 * c2 *
               (2 * C0422 * c10 * c19 * c21 + 3 * C0431 * c1 * c18 * c25 +
                C0421 * c19 * c25);
  out2(0, 1) = 0;
  out2(0, 2) = 0;
  out2(0, 3) = 0;
  out2(0, 4) = c95;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) = 0;
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) = 2 * c32 *
               (3 * C2431 * c25 * c31 * c36 + 2 * C2422 * c10 * c21 * c37 +
                C2421 * c25 * c37);
  out2(2, 3) = 0;
  out2(2, 4) = c108;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = 2 * c48 *
               (3 * C3431 * c25 * c47 * c51 + 2 * C3422 * c10 * c21 * c52 +
                C3421 * c25 * c52);
  out2(3, 4) = c121;
  out2(3, 5) = 0;
  out2(4, 0) = c95;
  out2(4, 1) = 0;
  out2(4, 2) = c108;
  out2(4, 3) = c121;
  out2(4, 4) = 2 * c61 *
               (12 * C44 * c10 + 3 * C43 * c14 + 3 * C0413 * c14 * c18 * c64 +
                3 * C2413 * c14 * c36 * c66 + 3 * C3413 * c14 * c51 * c68);
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) = 0;
  out3(0) = c2 * (6 * C0413 * c10 * c12 * c15 + 2 * C0412 * c12 * c21 +
                  4 * C0422 * c18 * c19 * c21);
  out3(1) = 0;
  out3(2) = c32 * (6 * C2413 * c10 * c15 * c33 + 2 * C2412 * c21 * c33 +
                   4 * C2422 * c21 * c36 * c37);
  out3(3) = c48 * (6 * C3413 * c10 * c15 * c49 + 2 * C3412 * c21 * c49 +
                   4 * C3422 * c21 * c51 * c52);
  out3(4) =
      c61 * (24 * C44 * c10 + 6 * C43 * c14 + 6 * C0413 * c14 * c18 * c64 +
             6 * C2413 * c14 * c36 * c66 + 6 * C3413 * c14 * c51 * c68);
  out3(5) = 0;
  out4(0, 0) = 4 * C0422 * c21 * c72;
  out4(0, 1) = 0;
  out4(0, 2) = 0;
  out4(0, 3) = 0;
  out4(0, 4) = c152;
  out4(0, 5) = 0;
  out4(1, 0) = 0;
  out4(1, 1) = 0;
  out4(1, 2) = 0;
  out4(1, 3) = 0;
  out4(1, 4) = 0;
  out4(1, 5) = 0;
  out4(2, 0) = 0;
  out4(2, 1) = 0;
  out4(2, 2) = 4 * C2422 * c21 * c75;
  out4(2, 3) = 0;
  out4(2, 4) = c154;
  out4(2, 5) = 0;
  out4(3, 0) = 0;
  out4(3, 1) = 0;
  out4(3, 2) = 0;
  out4(3, 3) = 4 * C3422 * c21 * c78;
  out4(3, 4) = c156;
  out4(3, 5) = 0;
  out4(4, 0) = c152;
  out4(4, 1) = 0;
  out4(4, 2) = c154;
  out4(4, 3) = c156;
  out4(4, 4) = 24 * C44 * c61;
  out4(4, 5) = 0;
  out4(5, 0) = 0;
  out4(5, 1) = 0;
  out4(5, 2) = 0;
  out4(5, 3) = 0;
  out4(5, 4) = 0;
  out4(5, 5) = 0;
  return std::make_pair(std::make_pair(hessD, gradD),
                        std::make_pair(hessDD, gradDD));
}
std::pair<gradhess_type, gradhess_type>
FittedMaterial::gradhess_taylor_12_5(const strain_type &ek) {
  Mat6x6 hessD(0), hessDD(0);
  Vec6 gradD(0), gradDD(0);
  // MM order: {gradd, hessd, graddd, hessdd}
  auto out1 = [&](int i) -> Real & { return gradD[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hessD(i, j); };
  auto out3 = [&](int i) -> Real & { return gradDD[i]; };
  auto out4 = [&](int i, int j) -> Real & { return hessDD(i, j); };
  Real c1 = ekscale(0);
  Real c10 = ek(5);
  Real c14 = ekscale(5);
  Real c12 = Power(c1, 3);
  Real c21 = Power(c14, -2);
  Real c17 = ek(0);
  Real c18 = -1 + c17;
  Real c19 = Power(c1, 2);
  Real c25 = 1 / c14;
  Real c11 = Power(c10, 2);
  Real c31 = ekscale(2);
  Real c15 = Power(c14, -3);
  Real c33 = Power(c31, 3);
  Real c35 = ek(2);
  Real c36 = -1 + c35;
  Real c37 = Power(c31, 2);
  Real c47 = ekscale(3);
  Real c49 = Power(c47, 3);
  Real c51 = ek(3);
  Real c52 = Power(c47, 2);
  Real c24 = Power(c18, 2);
  Real c70 = Power(c14, 2);
  Real c64 = 1 / c1;
  Real c40 = Power(c36, 2);
  Real c66 = 1 / c31;
  Real c55 = Power(c51, 2);
  Real c68 = 1 / c47;
  Real c2 = Power(c1, -4);
  Real c32 = Power(c31, -4);
  Real c48 = Power(c47, -4);
  Real c88 = Power(c1, -3);
  Real c89 = 6 * C0513 * c10 * c19;
  Real c90 = 4 * C0522 * c1 * c18;
  Real c91 = 2 * C0512 * c19;
  Real c92 = c90 + c91;
  Real c93 = c14 * c92;
  Real c94 = c89 + c93;
  Real c95 = c15 * c88 * c94;
  Real c101 = Power(c31, -3);
  Real c102 = 6 * C2513 * c10 * c37;
  Real c103 = 4 * C2522 * c31 * c36;
  Real c104 = 2 * C2512 * c37;
  Real c105 = c103 + c104;
  Real c106 = c105 * c14;
  Real c107 = c102 + c106;
  Real c108 = c101 * c107 * c15;
  Real c114 = Power(c47, -3);
  Real c115 = 6 * C3513 * c10 * c52;
  Real c116 = 4 * C3522 * c47 * c51;
  Real c117 = 2 * C3512 * c52;
  Real c118 = c116 + c117;
  Real c119 = c118 * c14;
  Real c120 = c115 + c119;
  Real c121 = c114 * c120 * c15;
  Real c61 = Power(c14, -4);
  Real c72 = Power(c1, -2);
  Real c75 = Power(c31, -2);
  Real c78 = Power(c47, -2);
  Real c152 = 6 * C0513 * c15 * c64;
  Real c154 = 6 * C2513 * c15 * c66;
  Real c156 = 6 * C3513 * c15 * c68;
  out1(0) = c2 * (3 * C0513 * c11 * c12 * c15 + 2 * C0512 * c10 * c12 * c21 +
                  4 * C0522 * c10 * c18 * c19 * c21 + C0511 * c12 * c25 +
                  2 * C0521 * c18 * c19 * c25 + 3 * C0531 * c1 * c24 * c25);
  out1(1) = 0;
  out1(2) = c32 * (3 * C2513 * c11 * c15 * c33 + 2 * C2512 * c10 * c21 * c33 +
                   C2511 * c25 * c33 + 4 * C2522 * c10 * c21 * c36 * c37 +
                   2 * C2521 * c25 * c36 * c37 + 3 * C2531 * c25 * c31 * c40);
  out1(3) = c48 * (3 * C3513 * c11 * c15 * c49 + 2 * C3512 * c10 * c21 * c49 +
                   C3511 * c25 * c49 + 4 * C3522 * c10 * c21 * c51 * c52 +
                   2 * C3521 * c25 * c51 * c52 + 3 * C3531 * c25 * c47 * c55);
  out1(4) = 0;
  out1(5) = c61 * (12 * C54 * c11 + 6 * C53 * c10 * c14 +
                   6 * C0513 * c10 * c14 * c18 * c64 +
                   6 * C2513 * c10 * c14 * c36 * c66 +
                   6 * C3513 * c10 * c14 * c51 * c68 + 2 * C52 * c70 +
                   2 * C0512 * c18 * c64 * c70 + 2 * C2512 * c36 * c66 * c70 +
                   2 * C3512 * c51 * c68 * c70 + 2 * C0522 * c24 * c70 * c72 +
                   2 * C2522 * c40 * c70 * c75 + 2 * C3522 * c55 * c70 * c78);
  out2(0, 0) = 2 * c2 *
               (2 * C0522 * c10 * c19 * c21 + 3 * C0531 * c1 * c18 * c25 +
                C0521 * c19 * c25);
  out2(0, 1) = 0;
  out2(0, 2) = 0;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = c95;
  out2(1, 0) = 0;
  out2(1, 1) = 0;
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) = 2 * c32 *
               (3 * C2531 * c25 * c31 * c36 + 2 * C2522 * c10 * c21 * c37 +
                C2521 * c25 * c37);
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = c108;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = 2 * c48 *
               (3 * C3531 * c25 * c47 * c51 + 2 * C3522 * c10 * c21 * c52 +
                C3521 * c25 * c52);
  out2(3, 4) = 0;
  out2(3, 5) = c121;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = 0;
  out2(4, 4) = 0;
  out2(4, 5) = 0;
  out2(5, 0) = c95;
  out2(5, 1) = 0;
  out2(5, 2) = c108;
  out2(5, 3) = c121;
  out2(5, 4) = 0;
  out2(5, 5) = 2 * c61 *
               (12 * C54 * c10 + 3 * C53 * c14 + 3 * C0513 * c14 * c18 * c64 +
                3 * C2513 * c14 * c36 * c66 + 3 * C3513 * c14 * c51 * c68);
  out3(0) = c2 * (6 * C0513 * c10 * c12 * c15 + 2 * C0512 * c12 * c21 +
                  4 * C0522 * c18 * c19 * c21);
  out3(1) = 0;
  out3(2) = c32 * (6 * C2513 * c10 * c15 * c33 + 2 * C2512 * c21 * c33 +
                   4 * C2522 * c21 * c36 * c37);
  out3(3) = c48 * (6 * C3513 * c10 * c15 * c49 + 2 * C3512 * c21 * c49 +
                   4 * C3522 * c21 * c51 * c52);
  out3(4) = 0;
  out3(5) =
      c61 * (24 * C54 * c10 + 6 * C53 * c14 + 6 * C0513 * c14 * c18 * c64 +
             6 * C2513 * c14 * c36 * c66 + 6 * C3513 * c14 * c51 * c68);
  out4(0, 0) = 4 * C0522 * c21 * c72;
  out4(0, 1) = 0;
  out4(0, 2) = 0;
  out4(0, 3) = 0;
  out4(0, 4) = 0;
  out4(0, 5) = c152;
  out4(1, 0) = 0;
  out4(1, 1) = 0;
  out4(1, 2) = 0;
  out4(1, 3) = 0;
  out4(1, 4) = 0;
  out4(1, 5) = 0;
  out4(2, 0) = 0;
  out4(2, 1) = 0;
  out4(2, 2) = 4 * C2522 * c21 * c75;
  out4(2, 3) = 0;
  out4(2, 4) = 0;
  out4(2, 5) = c154;
  out4(3, 0) = 0;
  out4(3, 1) = 0;
  out4(3, 2) = 0;
  out4(3, 3) = 4 * C3522 * c21 * c78;
  out4(3, 4) = 0;
  out4(3, 5) = c156;
  out4(4, 0) = 0;
  out4(4, 1) = 0;
  out4(4, 2) = 0;
  out4(4, 3) = 0;
  out4(4, 4) = 0;
  out4(4, 5) = 0;
  out4(5, 0) = c152;
  out4(5, 1) = 0;
  out4(5, 2) = c154;
  out4(5, 3) = c156;
  out4(5, 4) = 0;
  out4(5, 5) = 24 * C54 * c61;
  return std::make_pair(std::make_pair(hessD, gradD),
                        std::make_pair(hessDD, gradDD));
}

// helper for selecting correct 1d derivative
std::pair<value_type, value_type>
FittedMaterial::psi_taylor_12_i(const strain_type &ek, int i) {
  if (i == 0)
    return psi_taylor_12_0(ek);
  else if (i == 1)
    return psi_taylor_12_1(ek);
  else if (i == 2)
    return psi_taylor_12_2(ek);
  else if (i == 3)
    return psi_taylor_12_3(ek);
  else if (i == 4)
    return psi_taylor_12_4(ek);
  else
    return psi_taylor_12_5(ek);
}

std::pair<grad_type, grad_type>
FittedMaterial::grad_taylor_12_i(const strain_type &ek, int i) {
  if (i == 0)
    return grad_taylor_12_0(ek);
  else if (i == 1)
    return grad_taylor_12_1(ek);
  else if (i == 2)
    return grad_taylor_12_2(ek);
  else if (i == 3)
    return grad_taylor_12_3(ek);
  else if (i == 4)
    return grad_taylor_12_4(ek);
  else
    return grad_taylor_12_5(ek);
}

std::pair<gradhess_type, gradhess_type>
FittedMaterial::gradhess_taylor_12_i(const strain_type &ek, int i) {
  if (i == 0)
    return gradhess_taylor_12_0(ek);
  else if (i == 1)
    return gradhess_taylor_12_1(ek);
  else if (i == 2)
    return gradhess_taylor_12_2(ek);
  else if (i == 3)
    return gradhess_taylor_12_3(ek);
  else if (i == 4)
    return gradhess_taylor_12_4(ek);
  else
    return gradhess_taylor_12_5(ek);
}
