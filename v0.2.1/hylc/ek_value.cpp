#include "hylc_mm.hpp"

using namespace hylc;
using namespace hylc::mathematica;

Vec6 hylc::mathematica::ek(const Vec18 &_xlocal, const Mat2x2 &_invDm,
                           const Real &A, const Real &l0, const Real &l1,
                           const Real &l2, const Vec2 &t0, const Vec2 &t1,
                           const Vec2 &t2) {
  // define input
  auto xlocal = [&](int i) -> const Real & {
    assert(i >= 1 && i <= 18);
    return _xlocal[i - 1];
  };
  auto invDm = [&](int i, int j) -> const Real & {
    assert(i >= 1 && i <= 2 && j >= 1 && j <= 2);
    return _invDm(i - 1, j - 1);
  };
  const Real &t01 = t0[0];
  const Real &t02 = t0[1];
  const Real &t11 = t1[0];
  const Real &t12 = t1[1];
  const Real &t21 = t2[0];
  const Real &t22 = t2[1];

  // define output
  Vec6 _out;
  auto out = [&](int i) -> Real & {
    assert(i >= 0 && i < 6);
    return _out[i];
  };

  Real c1 = -xlocal(4);
  Real c2 = xlocal(1) + c1;
  Real c3 = invDm(1, 1) * c2;
  Real c4 = -xlocal(7);
  Real c5 = xlocal(1) + c4;
  Real c6 = invDm(2, 1) * c5;
  Real c7 = c3 + c6;
  Real c9 = -xlocal(5);
  Real c10 = xlocal(2) + c9;
  Real c11 = invDm(1, 1) * c10;
  Real c12 = -xlocal(8);
  Real c14 = xlocal(2) + c12;
  Real c15 = invDm(2, 1) * c14;
  Real c16 = c11 + c15;
  Real c18 = -xlocal(6);
  Real c19 = xlocal(3) + c18;
  Real c21 = invDm(1, 1) * c19;
  Real c22 = -xlocal(9);
  Real c23 = xlocal(3) + c22;
  Real c24 = invDm(2, 1) * c23;
  Real c25 = c21 + c24;
  Real c28 = invDm(1, 2) * c2;
  Real c29 = invDm(2, 2) * c5;
  Real c30 = c28 + c29;
  Real c32 = invDm(1, 2) * c10;
  Real c33 = invDm(2, 2) * c14;
  Real c34 = c32 + c33;
  Real c36 = invDm(1, 2) * c19;
  Real c37 = invDm(2, 2) * c23;
  Real c38 = c36 + c37;
  Real c54 = xlocal(2) * xlocal(2);
  Real c52 = xlocal(3) * xlocal(3);
  Real c61 = xlocal(7) * xlocal(7);
  Real c79 = xlocal(8) * xlocal(8);
  Real c100 = xlocal(9) * xlocal(9);
  Real c114 = xlocal(4) + c4;
  Real c120 = -2 * xlocal(4);
  Real c121 = xlocal(7) + c120;
  Real c119 = -(xlocal(5) * xlocal(7) * xlocal(8));
  Real c123 = xlocal(4) * c79;
  Real c124 = -(xlocal(6) * xlocal(7) * xlocal(9));
  Real c126 = xlocal(4) * c100;
  Real c105 = xlocal(1) * xlocal(1);
  Real c141 = xlocal(2) * c114;
  Real c192 = xlocal(4) * xlocal(4);
  Real c199 = xlocal(5) * xlocal(5);
  Real c207 = xlocal(6) * xlocal(6);
  Real c58 = xlocal(4) * xlocal(7) * c54;
  Real c60 = xlocal(4) * xlocal(7) * c52;
  Real c73 = xlocal(5) * xlocal(8) * c52;
  Real c76 = -(xlocal(2) * xlocal(3) * xlocal(6) * xlocal(8));
  Real c108 = -(xlocal(5) * xlocal(8));
  Real c109 = xlocal(6) + c22;
  Real c88 = -(xlocal(2) * xlocal(3) * xlocal(5) * xlocal(9));
  Real c90 = xlocal(6) * xlocal(9) * c54;
  Real c115 = c114 * c54;
  Real c116 = c114 * c52;
  Real c117 = xlocal(2) * xlocal(5) * xlocal(7);
  Real c245 = -(xlocal(7) * c199);
  Real c118 = xlocal(3) * xlocal(6) * xlocal(7);
  Real c247 = -(xlocal(7) * c207);
  Real c248 = xlocal(4) * xlocal(5) * xlocal(8);
  Real c170 = xlocal(5) * xlocal(7);
  Real c251 = xlocal(4) * xlocal(6) * xlocal(9);
  Real c272 = xlocal(3) * xlocal(4);
  Real c274 = -(xlocal(3) * xlocal(7));
  Real c275 = xlocal(6) * xlocal(7);
  Real c213 = -(xlocal(2) * xlocal(4) * xlocal(5) * xlocal(7));
  Real c215 = -(xlocal(3) * xlocal(4) * xlocal(6) * xlocal(7));
  Real c67 = xlocal(2) * xlocal(5) * c61;
  Real c69 = xlocal(3) * xlocal(6) * c61;
  Real c219 = xlocal(2) * xlocal(8) * c192;
  Real c224 = -(xlocal(3) * xlocal(5) * xlocal(6) * xlocal(8));
  Real c226 = xlocal(2) * xlocal(8) * c207;
  Real c78 = -(xlocal(2) * xlocal(4) * xlocal(7) * xlocal(8));
  Real c83 = xlocal(3) * xlocal(6) * c79;
  Real c236 = xlocal(3) * xlocal(9) * c192;
  Real c240 = xlocal(3) * xlocal(9) * c199;
  Real c243 = -(xlocal(2) * xlocal(5) * xlocal(6) * xlocal(9));
  Real c92 = -(xlocal(3) * xlocal(4) * xlocal(7) * xlocal(9));
  Real c97 = -(xlocal(3) * xlocal(5) * xlocal(8) * xlocal(9));
  Real c99 = -(xlocal(2) * xlocal(6) * xlocal(8) * xlocal(9));
  Real c104 = xlocal(2) * xlocal(5) * c100;
  Real c106 = xlocal(5) + c12;
  Real c361 = -2 * xlocal(4) * xlocal(7);
  Real c362 = -2 * xlocal(5) * xlocal(8);
  Real c387 = -2 * xlocal(6) * xlocal(9);
  Real c163 = -(xlocal(3) * xlocal(5) * xlocal(7));
  Real c164 = xlocal(2) * xlocal(6) * xlocal(7);
  Real c166 = xlocal(3) * xlocal(4) * xlocal(8);
  Real c167 = -(xlocal(1) * xlocal(6) * xlocal(8));
  Real c169 = -(xlocal(1) * xlocal(5));
  Real c171 = xlocal(1) * xlocal(8);
  Real c172 = -(xlocal(4) * xlocal(8));
  Real c173 = c141 + c169 + c170 + c171 + c172;
  Real c176 = -(xlocal(2) * xlocal(4) * xlocal(9));
  Real c177 = xlocal(1) * xlocal(5) * xlocal(9);
  Real c179 = -(xlocal(3) * xlocal(4));
  Real c180 = xlocal(1) * xlocal(6);
  Real c181 = xlocal(3) * xlocal(7);
  Real c182 = -(xlocal(6) * xlocal(7));
  Real c183 = -(xlocal(1) * xlocal(9));
  Real c184 = xlocal(4) * xlocal(9);
  Real c185 = c179 + c180 + c181 + c182 + c183 + c184;
  Real c46 = 1 / A;
  Real c47 = 1 / l0;
  Real c48 = 1 / l1;
  Real c49 = 1 / l2;
  Real c51 = xlocal(18) * xlocal(2) * xlocal(3) * xlocal(5);
  Real c53 = -(xlocal(17) * xlocal(5) * c52);
  Real c55 = -(xlocal(18) * xlocal(6) * c54);
  Real c56 = xlocal(17) * xlocal(2) * xlocal(3) * xlocal(6);
  Real c57 = -(xlocal(17) * xlocal(2) * xlocal(4) * xlocal(7));
  Real c59 = -(xlocal(18) * xlocal(3) * xlocal(4) * xlocal(7));
  Real c62 = xlocal(17) * xlocal(2) * c61;
  Real c63 = -(c54 * c61);
  Real c64 = xlocal(18) * xlocal(3) * c61;
  Real c65 = -(c52 * c61);
  Real c66 = -(xlocal(17) * xlocal(5) * c61);
  Real c68 = -(xlocal(18) * xlocal(6) * c61);
  Real c70 = -(xlocal(18) * xlocal(2) * xlocal(3) * xlocal(8));
  Real c71 = xlocal(17) * xlocal(8) * c52;
  Real c72 = -(xlocal(18) * xlocal(3) * xlocal(5) * xlocal(8));
  Real c74 = 2 * xlocal(18) * xlocal(2) * xlocal(6) * xlocal(8);
  Real c75 = -(xlocal(17) * xlocal(3) * xlocal(6) * xlocal(8));
  Real c77 = xlocal(17) * xlocal(4) * xlocal(7) * xlocal(8);
  Real c80 = xlocal(18) * xlocal(3) * c79;
  Real c81 = -(c52 * c79);
  Real c82 = -(xlocal(18) * xlocal(6) * c79);
  Real c84 = xlocal(18) * xlocal(9) * c54;
  Real c85 = -(xlocal(17) * xlocal(2) * xlocal(3) * xlocal(9));
  Real c86 = -(xlocal(18) * xlocal(2) * xlocal(5) * xlocal(9));
  Real c87 = 2 * xlocal(17) * xlocal(3) * xlocal(5) * xlocal(9);
  Real c89 = -(xlocal(17) * xlocal(2) * xlocal(6) * xlocal(9));
  Real c91 = xlocal(18) * xlocal(4) * xlocal(7) * xlocal(9);
  Real c93 = -(xlocal(18) * xlocal(2) * xlocal(8) * xlocal(9));
  Real c94 = -(xlocal(17) * xlocal(3) * xlocal(8) * xlocal(9));
  Real c95 = 2 * xlocal(2) * xlocal(3) * xlocal(8) * xlocal(9);
  Real c96 = xlocal(18) * xlocal(5) * xlocal(8) * xlocal(9);
  Real c98 = xlocal(17) * xlocal(6) * xlocal(8) * xlocal(9);
  Real c101 = xlocal(17) * xlocal(2) * c100;
  Real c102 = -(c100 * c54);
  Real c103 = -(xlocal(17) * xlocal(5) * c100);
  Real c107 = xlocal(17) * c106;
  Real c110 = xlocal(18) * c109;
  Real c111 = -(xlocal(6) * xlocal(9));
  Real c112 = c100 + c107 + c108 + c110 + c111 + c79;
  Real c113 = -(c105 * c112);
  Real c122 = xlocal(2) * xlocal(8) * c121;
  Real c125 = xlocal(3) * xlocal(9) * c121;
  Real c127 =
      c115 + c116 + c117 + c118 + c119 + c122 + c123 + c124 + c125 + c126;
  Real c128 = -(xlocal(16) * c127);
  Real c129 = xlocal(18) * xlocal(3) * xlocal(4);
  Real c130 = xlocal(16) * xlocal(2) * xlocal(5);
  Real c131 = xlocal(16) * xlocal(3) * xlocal(6);
  Real c132 = -(xlocal(18) * xlocal(3) * xlocal(7));
  Real c133 = -(xlocal(2) * xlocal(5) * xlocal(7));
  Real c134 = 2 * xlocal(18) * xlocal(6) * xlocal(7);
  Real c135 = -(xlocal(3) * xlocal(6) * xlocal(7));
  Real c136 = -(xlocal(16) * xlocal(2) * xlocal(8));
  Real c137 = -(xlocal(2) * xlocal(4) * xlocal(8));
  Real c138 = -(xlocal(16) * xlocal(5) * xlocal(8));
  Real c139 = 2 * xlocal(2) * xlocal(7) * xlocal(8);
  Real c140 = xlocal(16) * c79;
  Real c142 = 2 * xlocal(5) * xlocal(7);
  Real c143 = xlocal(4) + xlocal(7);
  Real c144 = -(xlocal(8) * c143);
  Real c145 = c141 + c142 + c144;
  Real c146 = xlocal(17) * c145;
  Real c147 = -(xlocal(16) * xlocal(3) * xlocal(9));
  Real c148 = -(xlocal(18) * xlocal(4) * xlocal(9));
  Real c149 = -(xlocal(3) * xlocal(4) * xlocal(9));
  Real c150 = -(xlocal(16) * xlocal(6) * xlocal(9));
  Real c151 = -(xlocal(18) * xlocal(7) * xlocal(9));
  Real c152 = 2 * xlocal(3) * xlocal(7) * xlocal(9);
  Real c153 = xlocal(16) * c100;
  Real c154 = c119 + c123 + c124 + c126 + c129 + c130 + c131 + c132 + c133 +
              c134 + c135 + c136 + c137 + c138 + c139 + c140 + c146 + c147 +
              c148 + c149 + c150 + c151 + c152 + c153;
  Real c155 = xlocal(1) * c154;
  Real c156 = c101 + c102 + c103 + c104 + c113 + c128 + c155 + c51 + c53 + c55 +
              c56 + c57 + c58 + c59 + c60 + c62 + c63 + c64 + c65 + c66 + c67 +
              c68 + c69 + c70 + c71 + c72 + c73 + c74 + c75 + c76 + c77 + c78 +
              c80 + c81 + c82 + c83 + c84 + c85 + c86 + c87 + c88 + c89 + c90 +
              c91 + c92 + c93 + c94 + c95 + c96 + c97 + c98 + c99;
  Real c157 = -2 * xlocal(1) * xlocal(7);
  Real c158 = -2 * xlocal(2) * xlocal(8);
  Real c159 = -2 * xlocal(3) * xlocal(9);
  Real c160 = c100 + c105 + c157 + c158 + c159 + c52 + c54 + c61 + c79;
  Real c161 = xlocal(16) * xlocal(3) * xlocal(5);
  Real c162 = -(xlocal(16) * xlocal(2) * xlocal(6));
  Real c165 = -(xlocal(16) * xlocal(3) * xlocal(8));
  Real c168 = xlocal(16) * xlocal(6) * xlocal(8);
  Real c174 = xlocal(18) * c173;
  Real c175 = xlocal(16) * xlocal(2) * xlocal(9);
  Real c178 = -(xlocal(16) * xlocal(5) * xlocal(9));
  Real c186 = xlocal(17) * c185;
  Real c187 = c161 + c162 + c163 + c164 + c165 + c166 + c167 + c168 + c174 +
              c175 + c176 + c177 + c178 + c186;
  Real c188 = c160 * c187;
  Real c189 = ArcTan(c156, c188);
  Real c193 = xlocal(11) * xlocal(2) * c192;
  Real c194 = -(c192 * c54);
  Real c195 = xlocal(12) * xlocal(3) * c192;
  Real c196 = -(c192 * c52);
  Real c197 = -(xlocal(12) * xlocal(2) * xlocal(3) * xlocal(5));
  Real c198 = xlocal(11) * xlocal(5) * c52;
  Real c200 = xlocal(12) * xlocal(3) * c199;
  Real c201 = -(c199 * c52);
  Real c202 = xlocal(12) * xlocal(6) * c54;
  Real c203 = -(xlocal(11) * xlocal(2) * xlocal(3) * xlocal(6));
  Real c204 = -(xlocal(12) * xlocal(2) * xlocal(5) * xlocal(6));
  Real c205 = -(xlocal(11) * xlocal(3) * xlocal(5) * xlocal(6));
  Real c206 = 2 * xlocal(2) * xlocal(3) * xlocal(5) * xlocal(6);
  Real c208 = xlocal(11) * xlocal(2) * c207;
  Real c209 = -(c207 * c54);
  Real c210 = -(xlocal(11) * xlocal(2) * xlocal(4) * xlocal(7));
  Real c211 = -(xlocal(12) * xlocal(3) * xlocal(4) * xlocal(7));
  Real c212 = xlocal(11) * xlocal(4) * xlocal(5) * xlocal(7);
  Real c214 = xlocal(12) * xlocal(4) * xlocal(6) * xlocal(7);
  Real c216 = xlocal(12) * xlocal(2) * xlocal(3) * xlocal(8);
  Real c217 = -(xlocal(11) * xlocal(8) * c52);
  Real c218 = -(xlocal(11) * xlocal(8) * c192);
  Real c220 = -(xlocal(12) * xlocal(3) * xlocal(5) * xlocal(8));
  Real c221 = -(xlocal(12) * xlocal(2) * xlocal(6) * xlocal(8));
  Real c222 = 2 * xlocal(11) * xlocal(3) * xlocal(6) * xlocal(8);
  Real c223 = xlocal(12) * xlocal(5) * xlocal(6) * xlocal(8);
  Real c225 = -(xlocal(11) * xlocal(8) * c207);
  Real c227 = xlocal(8) + c9;
  Real c228 = xlocal(11) * c227;
  Real c229 = xlocal(12) + c18;
  Real c230 = -(c109 * c229);
  Real c231 = c108 + c199 + c228 + c230;
  Real c232 = -(c105 * c231);
  Real c233 = -(xlocal(12) * xlocal(9) * c54);
  Real c234 = xlocal(11) * xlocal(2) * xlocal(3) * xlocal(9);
  Real c235 = -(xlocal(12) * xlocal(9) * c192);
  Real c237 = 2 * xlocal(12) * xlocal(2) * xlocal(5) * xlocal(9);
  Real c238 = -(xlocal(11) * xlocal(3) * xlocal(5) * xlocal(9));
  Real c239 = -(xlocal(12) * xlocal(9) * c199);
  Real c241 = -(xlocal(11) * xlocal(2) * xlocal(6) * xlocal(9));
  Real c242 = xlocal(11) * xlocal(5) * xlocal(6) * xlocal(9);
  Real c244 = 2 * xlocal(2) * xlocal(5) * xlocal(7);
  Real c246 = 2 * xlocal(3) * xlocal(6) * xlocal(7);
  Real c249 = xlocal(5) + xlocal(8);
  Real c250 = -(xlocal(2) * xlocal(4) * c249);
  Real c252 = xlocal(6) + xlocal(9);
  Real c253 = -(xlocal(3) * xlocal(4) * c252);
  Real c254 =
      c115 + c116 + c244 + c245 + c246 + c247 + c248 + c250 + c251 + c253;
  Real c255 = xlocal(10) * c254;
  Real c256 = xlocal(10) * xlocal(2) * xlocal(5);
  Real c257 = -2 * xlocal(2) * xlocal(4) * xlocal(5);
  Real c258 = -(xlocal(10) * c199);
  Real c259 = xlocal(10) * xlocal(3) * xlocal(6);
  Real c260 = -2 * xlocal(3) * xlocal(4) * xlocal(6);
  Real c261 = -(xlocal(10) * c207);
  Real c262 = -(xlocal(10) * xlocal(2) * xlocal(8));
  Real c263 = xlocal(2) * xlocal(4) * xlocal(8);
  Real c264 = xlocal(10) * xlocal(5) * xlocal(8);
  Real c265 = xlocal(4) * xlocal(5);
  Real c266 = -2 * xlocal(4) * xlocal(8);
  Real c267 = c141 + c170 + c265 + c266;
  Real c268 = xlocal(11) * c267;
  Real c269 = -(xlocal(10) * xlocal(3) * xlocal(9));
  Real c270 = xlocal(3) * xlocal(4) * xlocal(9);
  Real c271 = xlocal(10) * xlocal(6) * xlocal(9);
  Real c273 = xlocal(4) * xlocal(6);
  Real c276 = -2 * xlocal(4) * xlocal(9);
  Real c277 = c272 + c273 + c274 + c275 + c276;
  Real c278 = xlocal(12) * c277;
  Real c279 = c117 + c118 + c245 + c247 + c248 + c251 + c256 + c257 + c258 +
              c259 + c260 + c261 + c262 + c263 + c264 + c268 + c269 + c270 +
              c271 + c278;
  Real c280 = -(xlocal(1) * c279);
  Real c281 = c193 + c194 + c195 + c196 + c197 + c198 + c200 + c201 + c202 +
              c203 + c204 + c205 + c206 + c208 + c209 + c210 + c211 + c212 +
              c213 + c214 + c215 + c216 + c217 + c218 + c219 + c220 + c221 +
              c222 + c223 + c224 + c225 + c226 + c232 + c233 + c234 + c235 +
              c236 + c237 + c238 + c239 + c240 + c241 + c242 + c243 + c255 +
              c280 + c58 + c60 + c73 + c76 + c88 + c90;
  Real c282 = -2 * xlocal(1) * xlocal(4);
  Real c283 = -2 * xlocal(2) * xlocal(5);
  Real c284 = -2 * xlocal(3) * xlocal(6);
  Real c285 = c105 + c192 + c199 + c207 + c282 + c283 + c284 + c52 + c54;
  Real c286 = -(xlocal(10) * xlocal(3) * xlocal(5));
  Real c287 = xlocal(10) * xlocal(2) * xlocal(6);
  Real c288 = xlocal(3) * xlocal(5) * xlocal(7);
  Real c289 = -(xlocal(2) * xlocal(6) * xlocal(7));
  Real c290 = xlocal(10) * xlocal(3) * xlocal(8);
  Real c291 = -(xlocal(3) * xlocal(4) * xlocal(8));
  Real c292 = xlocal(1) * xlocal(6) * xlocal(8);
  Real c293 = -(xlocal(10) * xlocal(6) * xlocal(8));
  Real c294 = xlocal(1) * xlocal(5);
  Real c295 = -(xlocal(5) * xlocal(7));
  Real c296 = xlocal(7) + c1;
  Real c297 = xlocal(2) * c296;
  Real c298 = -(xlocal(1) * xlocal(8));
  Real c299 = xlocal(4) * xlocal(8);
  Real c300 = c294 + c295 + c297 + c298 + c299;
  Real c301 = xlocal(12) * c300;
  Real c302 = -(xlocal(10) * xlocal(2) * xlocal(9));
  Real c303 = xlocal(2) * xlocal(4) * xlocal(9);
  Real c304 = -(xlocal(1) * xlocal(5) * xlocal(9));
  Real c305 = xlocal(10) * xlocal(5) * xlocal(9);
  Real c306 = -(xlocal(1) * xlocal(6));
  Real c307 = xlocal(1) * xlocal(9);
  Real c308 = -(xlocal(4) * xlocal(9));
  Real c309 = c272 + c274 + c275 + c306 + c307 + c308;
  Real c310 = xlocal(11) * c309;
  Real c311 = c286 + c287 + c288 + c289 + c290 + c291 + c292 + c293 + c301 +
              c302 + c303 + c304 + c305 + c310;
  Real c312 = -(c285 * c311);
  Real c313 = ArcTan(c281, c312);
  Real c316 = xlocal(13) * xlocal(2) * xlocal(4) * xlocal(5);
  Real c317 = -(xlocal(1) * xlocal(13) * c199);
  Real c318 = xlocal(13) * xlocal(3) * xlocal(4) * xlocal(6);
  Real c319 = -(xlocal(1) * xlocal(13) * c207);
  Real c320 = -(xlocal(13) * xlocal(2) * xlocal(5) * xlocal(7));
  Real c321 = xlocal(1) * xlocal(7) * c199;
  Real c322 = xlocal(13) * xlocal(7) * c199;
  Real c323 = -(xlocal(13) * xlocal(3) * xlocal(6) * xlocal(7));
  Real c324 = xlocal(1) * xlocal(7) * c207;
  Real c325 = xlocal(13) * xlocal(7) * c207;
  Real c326 = -(c199 * c61);
  Real c327 = -(c207 * c61);
  Real c328 = -(xlocal(13) * xlocal(2) * xlocal(4) * xlocal(8));
  Real c329 = 2 * xlocal(1) * xlocal(13) * xlocal(5) * xlocal(8);
  Real c330 = -(xlocal(1) * xlocal(4) * xlocal(5) * xlocal(8));
  Real c331 = -(xlocal(13) * xlocal(4) * xlocal(5) * xlocal(8));
  Real c332 = xlocal(13) * xlocal(2) * xlocal(7) * xlocal(8);
  Real c333 = -(xlocal(1) * xlocal(5) * xlocal(7) * xlocal(8));
  Real c334 = -(xlocal(13) * xlocal(5) * xlocal(7) * xlocal(8));
  Real c335 = 2 * xlocal(4) * xlocal(5) * xlocal(7) * xlocal(8);
  Real c336 = -(xlocal(1) * xlocal(13) * c79);
  Real c337 = xlocal(1) * xlocal(4) * c79;
  Real c338 = xlocal(13) * xlocal(4) * c79;
  Real c339 = -(c192 * c79);
  Real c340 = -(c207 * c79);
  Real c341 = -(xlocal(13) * xlocal(3) * xlocal(4) * xlocal(9));
  Real c342 = 2 * xlocal(1) * xlocal(13) * xlocal(6) * xlocal(9);
  Real c343 = -(xlocal(1) * xlocal(4) * xlocal(6) * xlocal(9));
  Real c344 = -(xlocal(13) * xlocal(4) * xlocal(6) * xlocal(9));
  Real c345 = xlocal(13) * xlocal(3) * xlocal(7) * xlocal(9);
  Real c346 = -(xlocal(1) * xlocal(6) * xlocal(7) * xlocal(9));
  Real c347 = -(xlocal(13) * xlocal(6) * xlocal(7) * xlocal(9));
  Real c348 = 2 * xlocal(4) * xlocal(6) * xlocal(7) * xlocal(9);
  Real c349 = 2 * xlocal(5) * xlocal(6) * xlocal(8) * xlocal(9);
  Real c350 = -(xlocal(1) * xlocal(13) * c100);
  Real c351 = xlocal(1) * xlocal(4) * c100;
  Real c352 = xlocal(13) * xlocal(4) * c100;
  Real c353 = -(c100 * c192);
  Real c354 = -(c100 * c199);
  Real c355 = -(xlocal(2) * xlocal(5) * xlocal(6));
  Real c356 = xlocal(4) * xlocal(6) * xlocal(7);
  Real c357 = -(xlocal(6) * c61);
  Real c358 = xlocal(2) * xlocal(6) * xlocal(8);
  Real c359 = xlocal(5) * xlocal(6) * xlocal(8);
  Real c360 = -(xlocal(6) * c79);
  Real c363 = c192 + c199 + c361 + c362 + c61 + c79;
  Real c364 = xlocal(3) * c363;
  Real c365 = -(xlocal(1) * c109 * c114);
  Real c366 = -(xlocal(9) * c192);
  Real c367 = xlocal(2) * xlocal(5) * xlocal(9);
  Real c368 = -(xlocal(9) * c199);
  Real c369 = xlocal(4) * xlocal(7) * xlocal(9);
  Real c370 = -(xlocal(2) * xlocal(8) * xlocal(9));
  Real c371 = xlocal(5) * xlocal(8) * xlocal(9);
  Real c372 = c355 + c356 + c357 + c358 + c359 + c360 + c364 + c365 + c366 +
              c367 + c368 + c369 + c370 + c371;
  Real c373 = -(xlocal(15) * c372);
  Real c374 = -(xlocal(3) * xlocal(5) * xlocal(6));
  Real c375 = xlocal(4) * xlocal(5) * xlocal(7);
  Real c376 = -(xlocal(5) * c61);
  Real c377 = -(xlocal(1) * c106 * c114);
  Real c378 = -(xlocal(8) * c192);
  Real c379 = xlocal(3) * xlocal(6) * xlocal(8);
  Real c380 = -(xlocal(8) * c207);
  Real c381 = xlocal(4) * xlocal(7) * xlocal(8);
  Real c382 = xlocal(3) * xlocal(5) * xlocal(9);
  Real c383 = xlocal(5) * xlocal(6) * xlocal(9);
  Real c384 = -(xlocal(3) * xlocal(8) * xlocal(9));
  Real c385 = xlocal(6) * xlocal(8) * xlocal(9);
  Real c386 = -(xlocal(5) * c100);
  Real c388 = c100 + c192 + c207 + c361 + c387 + c61;
  Real c389 = xlocal(2) * c388;
  Real c390 = c374 + c375 + c376 + c377 + c378 + c379 + c380 + c381 + c382 +
              c383 + c384 + c385 + c386 + c389;
  Real c391 = -(xlocal(14) * c390);
  Real c392 = c104 + c213 + c215 + c219 + c224 + c226 + c236 + c240 + c243 +
              c316 + c317 + c318 + c319 + c320 + c321 + c322 + c323 + c324 +
              c325 + c326 + c327 + c328 + c329 + c330 + c331 + c332 + c333 +
              c334 + c335 + c336 + c337 + c338 + c339 + c340 + c341 + c342 +
              c343 + c344 + c345 + c346 + c347 + c348 + c349 + c350 + c351 +
              c352 + c353 + c354 + c373 + c391 + c67 + c69 + c78 + c83 + c92 +
              c97 + c99;
  Real c393 = c100 + c192 + c199 + c207 + c361 + c362 + c387 + c61 + c79;
  Real c394 = xlocal(13) * xlocal(3) * xlocal(5);
  Real c395 = -(xlocal(13) * xlocal(2) * xlocal(6));
  Real c396 = -(xlocal(13) * xlocal(3) * xlocal(8));
  Real c397 = xlocal(13) * xlocal(6) * xlocal(8);
  Real c398 = xlocal(15) * c173;
  Real c399 = xlocal(13) * xlocal(2) * xlocal(9);
  Real c400 = -(xlocal(13) * xlocal(5) * xlocal(9));
  Real c401 = xlocal(14) * c185;
  Real c402 = c163 + c164 + c166 + c167 + c176 + c177 + c394 + c395 + c396 +
              c397 + c398 + c399 + c400 + c401;
  Real c403 = c393 * c402;
  Real c404 = ArcTan(c392, c403);
  out(0) = c16 * c16 + c25 * c25 + c7 * c7;
  out(1) = c16 * c34 + c25 * c38 + c30 * c7;
  out(2) = c30 * c30 + c34 * c34 + c38 * c38;
  out(3) = ((l0 * l1 * t21 * t21 * c189 + l1 * l2 * t01 * t01 * c313 +
             l0 * l2 * t11 * t11 * c404) *
            c46 * c47 * c48 * c49) /
           2.;
  out(4) = ((l0 * l1 * t21 * t22 * c189 + l1 * l2 * t01 * t02 * c313 +
             l0 * l2 * t11 * t12 * c404) *
            c46 * c47 * c48 * c49) /
           2.;
  out(5) = ((l0 * l1 * t22 * t22 * c189 + l1 * l2 * t02 * t02 * c313 +
             l0 * l2 * t12 * t12 * c404) *
            c46 * c47 * c48 * c49) /
           2.;

  return _out;
}
