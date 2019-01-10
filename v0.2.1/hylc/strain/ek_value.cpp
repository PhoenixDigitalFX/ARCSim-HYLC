#include "strain.hpp"

using namespace hylc;
using namespace hylc::mathematica;

// (xlocal)(\d+) -> $1($2)
// (invDm)(\d)(\d) -> $1($2,$3)

Vec6 hylc::mathematica::ek(const Vec18 &_xlocal, const Mat2x2 &_invDm,
                           const Real &A, const Real &l0, const Real &l1,
                           const Real &l2, const Vec2 &t0, const Vec2 &t1,
                           const Vec2 &t2) {
  // define input
  const Real &xlocal1 = _xlocal[0];
  const Real &xlocal2 = _xlocal[1];
  const Real &xlocal3 = _xlocal[2];
  const Real &xlocal4 = _xlocal[3];
  const Real &xlocal5 = _xlocal[4];
  const Real &xlocal6 = _xlocal[5];
  const Real &xlocal7 = _xlocal[6];
  const Real &xlocal8 = _xlocal[7];
  const Real &xlocal9 = _xlocal[8];
  const Real &xlocal10 = _xlocal[9];
  const Real &xlocal11 = _xlocal[10];
  const Real &xlocal12 = _xlocal[11];
  const Real &xlocal13 = _xlocal[12];
  const Real &xlocal14 = _xlocal[13];
  const Real &xlocal15 = _xlocal[14];
  const Real &xlocal16 = _xlocal[15];
  const Real &xlocal17 = _xlocal[16];
  const Real &xlocal18 = _xlocal[17];
  const Real &invDm11 = _invDm(0, 0);
  const Real &invDm12 = _invDm(0, 1);
  const Real &invDm21 = _invDm(1, 0);
  const Real &invDm22 = _invDm(1, 1);
  // auto xlocal = [&](int i) -> const Real & {
  //   assert(i >= 1 && i <= 18);
  //   return _xlocal[i - 1];
  // };
  // auto invDm = [&](int i, int j) -> const Real & {
  //   assert(i >= 1 && i <= 2 && j >= 1 && j <= 2);
  //   return _invDm(i - 1, j - 1);
  // };
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
  Real c1 = -xlocal4;
  Real c2 = xlocal1 + c1;
  Real c3 = invDm11 * c2;
  Real c4 = -xlocal7;
  Real c5 = xlocal1 + c4;
  Real c6 = invDm21 * c5;
  Real c7 = c3 + c6;
  Real c9 = -xlocal5;
  Real c10 = xlocal2 + c9;
  Real c11 = invDm11 * c10;
  Real c12 = -xlocal8;
  Real c14 = xlocal2 + c12;
  Real c15 = invDm21 * c14;
  Real c16 = c11 + c15;
  Real c18 = -xlocal6;
  Real c19 = xlocal3 + c18;
  Real c21 = invDm11 * c19;
  Real c22 = -xlocal9;
  Real c23 = xlocal3 + c22;
  Real c24 = invDm21 * c23;
  Real c25 = c21 + c24;
  Real c28 = invDm12 * c2;
  Real c29 = invDm22 * c5;
  Real c30 = c28 + c29;
  Real c32 = invDm12 * c10;
  Real c33 = invDm22 * c14;
  Real c34 = c32 + c33;
  Real c36 = invDm12 * c19;
  Real c37 = invDm22 * c23;
  Real c38 = c36 + c37;
  Real c54 = xlocal2 * xlocal2;
  Real c52 = xlocal3 * xlocal3;
  Real c61 = xlocal7 * xlocal7;
  Real c79 = xlocal8 * xlocal8;
  Real c100 = xlocal9 * xlocal9;
  Real c114 = xlocal4 + c4;
  Real c120 = -2 * xlocal4;
  Real c121 = xlocal7 + c120;
  Real c119 = -(xlocal5 * xlocal7 * xlocal8);
  Real c123 = xlocal4 * c79;
  Real c124 = -(xlocal6 * xlocal7 * xlocal9);
  Real c126 = xlocal4 * c100;
  Real c105 = xlocal1 * xlocal1;
  Real c141 = xlocal2 * c114;
  Real c194 = xlocal4 * xlocal4;
  Real c201 = xlocal5 * xlocal5;
  Real c209 = xlocal6 * xlocal6;
  Real c58 = xlocal4 * xlocal7 * c54;
  Real c60 = xlocal4 * xlocal7 * c52;
  Real c73 = xlocal5 * xlocal8 * c52;
  Real c76 = -(xlocal2 * xlocal3 * xlocal6 * xlocal8);
  Real c108 = -(xlocal5 * xlocal8);
  Real c109 = xlocal6 + c22;
  Real c88 = -(xlocal2 * xlocal3 * xlocal5 * xlocal9);
  Real c90 = xlocal6 * xlocal9 * c54;
  Real c115 = c114 * c54;
  Real c116 = c114 * c52;
  Real c117 = xlocal2 * xlocal5 * xlocal7;
  Real c247 = -(xlocal7 * c201);
  Real c118 = xlocal3 * xlocal6 * xlocal7;
  Real c249 = -(xlocal7 * c209);
  Real c250 = xlocal4 * xlocal5 * xlocal8;
  Real c170 = xlocal5 * xlocal7;
  Real c253 = xlocal4 * xlocal6 * xlocal9;
  Real c274 = xlocal3 * xlocal4;
  Real c276 = -(xlocal3 * xlocal7);
  Real c277 = xlocal6 * xlocal7;
  Real c215 = -(xlocal2 * xlocal4 * xlocal5 * xlocal7);
  Real c217 = -(xlocal3 * xlocal4 * xlocal6 * xlocal7);
  Real c67 = xlocal2 * xlocal5 * c61;
  Real c69 = xlocal3 * xlocal6 * c61;
  Real c221 = xlocal2 * xlocal8 * c194;
  Real c226 = -(xlocal3 * xlocal5 * xlocal6 * xlocal8);
  Real c228 = xlocal2 * xlocal8 * c209;
  Real c78 = -(xlocal2 * xlocal4 * xlocal7 * xlocal8);
  Real c83 = xlocal3 * xlocal6 * c79;
  Real c238 = xlocal3 * xlocal9 * c194;
  Real c242 = xlocal3 * xlocal9 * c201;
  Real c245 = -(xlocal2 * xlocal5 * xlocal6 * xlocal9);
  Real c92 = -(xlocal3 * xlocal4 * xlocal7 * xlocal9);
  Real c97 = -(xlocal3 * xlocal5 * xlocal8 * xlocal9);
  Real c99 = -(xlocal2 * xlocal6 * xlocal8 * xlocal9);
  Real c104 = xlocal2 * xlocal5 * c100;
  Real c106 = xlocal5 + c12;
  Real c365 = -2 * xlocal4 * xlocal7;
  Real c366 = -2 * xlocal5 * xlocal8;
  Real c391 = -2 * xlocal6 * xlocal9;
  Real c163 = -(xlocal3 * xlocal5 * xlocal7);
  Real c164 = xlocal2 * xlocal6 * xlocal7;
  Real c166 = xlocal3 * xlocal4 * xlocal8;
  Real c167 = -(xlocal1 * xlocal6 * xlocal8);
  Real c169 = -(xlocal1 * xlocal5);
  Real c171 = xlocal1 * xlocal8;
  Real c172 = -(xlocal4 * xlocal8);
  Real c173 = c141 + c169 + c170 + c171 + c172;
  Real c176 = -(xlocal2 * xlocal4 * xlocal9);
  Real c177 = xlocal1 * xlocal5 * xlocal9;
  Real c179 = -(xlocal3 * xlocal4);
  Real c180 = xlocal1 * xlocal6;
  Real c181 = xlocal3 * xlocal7;
  Real c182 = -(xlocal6 * xlocal7);
  Real c183 = -(xlocal1 * xlocal9);
  Real c184 = xlocal4 * xlocal9;
  Real c185 = c179 + c180 + c181 + c182 + c183 + c184;
  Real c46 = 1 / A;
  Real c47 = 1 / l0;
  Real c48 = 1 / l1;
  Real c49 = 1 / l2;
  Real c51 = xlocal18 * xlocal2 * xlocal3 * xlocal5;
  Real c53 = -(xlocal17 * xlocal5 * c52);
  Real c55 = -(xlocal18 * xlocal6 * c54);
  Real c56 = xlocal17 * xlocal2 * xlocal3 * xlocal6;
  Real c57 = -(xlocal17 * xlocal2 * xlocal4 * xlocal7);
  Real c59 = -(xlocal18 * xlocal3 * xlocal4 * xlocal7);
  Real c62 = xlocal17 * xlocal2 * c61;
  Real c63 = -(c54 * c61);
  Real c64 = xlocal18 * xlocal3 * c61;
  Real c65 = -(c52 * c61);
  Real c66 = -(xlocal17 * xlocal5 * c61);
  Real c68 = -(xlocal18 * xlocal6 * c61);
  Real c70 = -(xlocal18 * xlocal2 * xlocal3 * xlocal8);
  Real c71 = xlocal17 * xlocal8 * c52;
  Real c72 = -(xlocal18 * xlocal3 * xlocal5 * xlocal8);
  Real c74 = 2 * xlocal18 * xlocal2 * xlocal6 * xlocal8;
  Real c75 = -(xlocal17 * xlocal3 * xlocal6 * xlocal8);
  Real c77 = xlocal17 * xlocal4 * xlocal7 * xlocal8;
  Real c80 = xlocal18 * xlocal3 * c79;
  Real c81 = -(c52 * c79);
  Real c82 = -(xlocal18 * xlocal6 * c79);
  Real c84 = xlocal18 * xlocal9 * c54;
  Real c85 = -(xlocal17 * xlocal2 * xlocal3 * xlocal9);
  Real c86 = -(xlocal18 * xlocal2 * xlocal5 * xlocal9);
  Real c87 = 2 * xlocal17 * xlocal3 * xlocal5 * xlocal9;
  Real c89 = -(xlocal17 * xlocal2 * xlocal6 * xlocal9);
  Real c91 = xlocal18 * xlocal4 * xlocal7 * xlocal9;
  Real c93 = -(xlocal18 * xlocal2 * xlocal8 * xlocal9);
  Real c94 = -(xlocal17 * xlocal3 * xlocal8 * xlocal9);
  Real c95 = 2 * xlocal2 * xlocal3 * xlocal8 * xlocal9;
  Real c96 = xlocal18 * xlocal5 * xlocal8 * xlocal9;
  Real c98 = xlocal17 * xlocal6 * xlocal8 * xlocal9;
  Real c101 = xlocal17 * xlocal2 * c100;
  Real c102 = -(c100 * c54);
  Real c103 = -(xlocal17 * xlocal5 * c100);
  Real c107 = xlocal17 * c106;
  Real c110 = xlocal18 * c109;
  Real c111 = -(xlocal6 * xlocal9);
  Real c112 = c100 + c107 + c108 + c110 + c111 + c79;
  Real c113 = -(c105 * c112);
  Real c122 = xlocal2 * xlocal8 * c121;
  Real c125 = xlocal3 * xlocal9 * c121;
  Real c127 =
      c115 + c116 + c117 + c118 + c119 + c122 + c123 + c124 + c125 + c126;
  Real c128 = -(xlocal16 * c127);
  Real c129 = xlocal18 * xlocal3 * xlocal4;
  Real c130 = xlocal16 * xlocal2 * xlocal5;
  Real c131 = xlocal16 * xlocal3 * xlocal6;
  Real c132 = -(xlocal18 * xlocal3 * xlocal7);
  Real c133 = -(xlocal2 * xlocal5 * xlocal7);
  Real c134 = 2 * xlocal18 * xlocal6 * xlocal7;
  Real c135 = -(xlocal3 * xlocal6 * xlocal7);
  Real c136 = -(xlocal16 * xlocal2 * xlocal8);
  Real c137 = -(xlocal2 * xlocal4 * xlocal8);
  Real c138 = -(xlocal16 * xlocal5 * xlocal8);
  Real c139 = 2 * xlocal2 * xlocal7 * xlocal8;
  Real c140 = xlocal16 * c79;
  Real c142 = 2 * xlocal5 * xlocal7;
  Real c143 = xlocal4 + xlocal7;
  Real c144 = -(xlocal8 * c143);
  Real c145 = c141 + c142 + c144;
  Real c146 = xlocal17 * c145;
  Real c147 = -(xlocal16 * xlocal3 * xlocal9);
  Real c148 = -(xlocal18 * xlocal4 * xlocal9);
  Real c149 = -(xlocal3 * xlocal4 * xlocal9);
  Real c150 = -(xlocal16 * xlocal6 * xlocal9);
  Real c151 = -(xlocal18 * xlocal7 * xlocal9);
  Real c152 = 2 * xlocal3 * xlocal7 * xlocal9;
  Real c153 = xlocal16 * c100;
  Real c154 = c119 + c123 + c124 + c126 + c129 + c130 + c131 + c132 + c133 +
              c134 + c135 + c136 + c137 + c138 + c139 + c140 + c146 + c147 +
              c148 + c149 + c150 + c151 + c152 + c153;
  Real c155 = xlocal1 * c154;
  Real c156 = c101 + c102 + c103 + c104 + c113 + c128 + c155 + c51 + c53 + c55 +
              c56 + c57 + c58 + c59 + c60 + c62 + c63 + c64 + c65 + c66 + c67 +
              c68 + c69 + c70 + c71 + c72 + c73 + c74 + c75 + c76 + c77 + c78 +
              c80 + c81 + c82 + c83 + c84 + c85 + c86 + c87 + c88 + c89 + c90 +
              c91 + c92 + c93 + c94 + c95 + c96 + c97 + c98 + c99;
  Real c157 = -2 * xlocal1 * xlocal7;
  Real c158 = -2 * xlocal2 * xlocal8;
  Real c159 = -2 * xlocal3 * xlocal9;
  Real c160 = c100 + c105 + c157 + c158 + c159 + c52 + c54 + c61 + c79;
  Real c161 = xlocal16 * xlocal3 * xlocal5;
  Real c162 = -(xlocal16 * xlocal2 * xlocal6);
  Real c165 = -(xlocal16 * xlocal3 * xlocal8);
  Real c168 = xlocal16 * xlocal6 * xlocal8;
  Real c174 = xlocal18 * c173;
  Real c175 = xlocal16 * xlocal2 * xlocal9;
  Real c178 = -(xlocal16 * xlocal5 * xlocal9);
  Real c186 = xlocal17 * c185;
  Real c187 = c161 + c162 + c163 + c164 + c165 + c166 + c167 + c168 + c174 +
              c175 + c176 + c177 + c178 + c186;
  Real c188 = c160 * c187;
  Real c189 = ArcTan(c156, c188);
  Real c190 = c189 / 2.;
  Real c191 = Tan(c190);
  Real c195 = xlocal11 * xlocal2 * c194;
  Real c196 = -(c194 * c54);
  Real c197 = xlocal12 * xlocal3 * c194;
  Real c198 = -(c194 * c52);
  Real c199 = -(xlocal12 * xlocal2 * xlocal3 * xlocal5);
  Real c200 = xlocal11 * xlocal5 * c52;
  Real c202 = xlocal12 * xlocal3 * c201;
  Real c203 = -(c201 * c52);
  Real c204 = xlocal12 * xlocal6 * c54;
  Real c205 = -(xlocal11 * xlocal2 * xlocal3 * xlocal6);
  Real c206 = -(xlocal12 * xlocal2 * xlocal5 * xlocal6);
  Real c207 = -(xlocal11 * xlocal3 * xlocal5 * xlocal6);
  Real c208 = 2 * xlocal2 * xlocal3 * xlocal5 * xlocal6;
  Real c210 = xlocal11 * xlocal2 * c209;
  Real c211 = -(c209 * c54);
  Real c212 = -(xlocal11 * xlocal2 * xlocal4 * xlocal7);
  Real c213 = -(xlocal12 * xlocal3 * xlocal4 * xlocal7);
  Real c214 = xlocal11 * xlocal4 * xlocal5 * xlocal7;
  Real c216 = xlocal12 * xlocal4 * xlocal6 * xlocal7;
  Real c218 = xlocal12 * xlocal2 * xlocal3 * xlocal8;
  Real c219 = -(xlocal11 * xlocal8 * c52);
  Real c220 = -(xlocal11 * xlocal8 * c194);
  Real c222 = -(xlocal12 * xlocal3 * xlocal5 * xlocal8);
  Real c223 = -(xlocal12 * xlocal2 * xlocal6 * xlocal8);
  Real c224 = 2 * xlocal11 * xlocal3 * xlocal6 * xlocal8;
  Real c225 = xlocal12 * xlocal5 * xlocal6 * xlocal8;
  Real c227 = -(xlocal11 * xlocal8 * c209);
  Real c229 = xlocal8 + c9;
  Real c230 = xlocal11 * c229;
  Real c231 = xlocal12 + c18;
  Real c232 = -(c109 * c231);
  Real c233 = c108 + c201 + c230 + c232;
  Real c234 = -(c105 * c233);
  Real c235 = -(xlocal12 * xlocal9 * c54);
  Real c236 = xlocal11 * xlocal2 * xlocal3 * xlocal9;
  Real c237 = -(xlocal12 * xlocal9 * c194);
  Real c239 = 2 * xlocal12 * xlocal2 * xlocal5 * xlocal9;
  Real c240 = -(xlocal11 * xlocal3 * xlocal5 * xlocal9);
  Real c241 = -(xlocal12 * xlocal9 * c201);
  Real c243 = -(xlocal11 * xlocal2 * xlocal6 * xlocal9);
  Real c244 = xlocal11 * xlocal5 * xlocal6 * xlocal9;
  Real c246 = 2 * xlocal2 * xlocal5 * xlocal7;
  Real c248 = 2 * xlocal3 * xlocal6 * xlocal7;
  Real c251 = xlocal5 + xlocal8;
  Real c252 = -(xlocal2 * xlocal4 * c251);
  Real c254 = xlocal6 + xlocal9;
  Real c255 = -(xlocal3 * xlocal4 * c254);
  Real c256 =
      c115 + c116 + c246 + c247 + c248 + c249 + c250 + c252 + c253 + c255;
  Real c257 = xlocal10 * c256;
  Real c258 = xlocal10 * xlocal2 * xlocal5;
  Real c259 = -2 * xlocal2 * xlocal4 * xlocal5;
  Real c260 = -(xlocal10 * c201);
  Real c261 = xlocal10 * xlocal3 * xlocal6;
  Real c262 = -2 * xlocal3 * xlocal4 * xlocal6;
  Real c263 = -(xlocal10 * c209);
  Real c264 = -(xlocal10 * xlocal2 * xlocal8);
  Real c265 = xlocal2 * xlocal4 * xlocal8;
  Real c266 = xlocal10 * xlocal5 * xlocal8;
  Real c267 = xlocal4 * xlocal5;
  Real c268 = -2 * xlocal4 * xlocal8;
  Real c269 = c141 + c170 + c267 + c268;
  Real c270 = xlocal11 * c269;
  Real c271 = -(xlocal10 * xlocal3 * xlocal9);
  Real c272 = xlocal3 * xlocal4 * xlocal9;
  Real c273 = xlocal10 * xlocal6 * xlocal9;
  Real c275 = xlocal4 * xlocal6;
  Real c278 = -2 * xlocal4 * xlocal9;
  Real c279 = c274 + c275 + c276 + c277 + c278;
  Real c280 = xlocal12 * c279;
  Real c281 = c117 + c118 + c247 + c249 + c250 + c253 + c258 + c259 + c260 +
              c261 + c262 + c263 + c264 + c265 + c266 + c270 + c271 + c272 +
              c273 + c280;
  Real c282 = -(xlocal1 * c281);
  Real c283 = c195 + c196 + c197 + c198 + c199 + c200 + c202 + c203 + c204 +
              c205 + c206 + c207 + c208 + c210 + c211 + c212 + c213 + c214 +
              c215 + c216 + c217 + c218 + c219 + c220 + c221 + c222 + c223 +
              c224 + c225 + c226 + c227 + c228 + c234 + c235 + c236 + c237 +
              c238 + c239 + c240 + c241 + c242 + c243 + c244 + c245 + c257 +
              c282 + c58 + c60 + c73 + c76 + c88 + c90;
  Real c284 = -2 * xlocal1 * xlocal4;
  Real c285 = -2 * xlocal2 * xlocal5;
  Real c286 = -2 * xlocal3 * xlocal6;
  Real c287 = c105 + c194 + c201 + c209 + c284 + c285 + c286 + c52 + c54;
  Real c288 = -(xlocal10 * xlocal3 * xlocal5);
  Real c289 = xlocal10 * xlocal2 * xlocal6;
  Real c290 = xlocal3 * xlocal5 * xlocal7;
  Real c291 = -(xlocal2 * xlocal6 * xlocal7);
  Real c292 = xlocal10 * xlocal3 * xlocal8;
  Real c293 = -(xlocal3 * xlocal4 * xlocal8);
  Real c294 = xlocal1 * xlocal6 * xlocal8;
  Real c295 = -(xlocal10 * xlocal6 * xlocal8);
  Real c296 = xlocal1 * xlocal5;
  Real c297 = -(xlocal5 * xlocal7);
  Real c298 = xlocal7 + c1;
  Real c299 = xlocal2 * c298;
  Real c300 = -(xlocal1 * xlocal8);
  Real c301 = xlocal4 * xlocal8;
  Real c302 = c296 + c297 + c299 + c300 + c301;
  Real c303 = xlocal12 * c302;
  Real c304 = -(xlocal10 * xlocal2 * xlocal9);
  Real c305 = xlocal2 * xlocal4 * xlocal9;
  Real c306 = -(xlocal1 * xlocal5 * xlocal9);
  Real c307 = xlocal10 * xlocal5 * xlocal9;
  Real c308 = -(xlocal1 * xlocal6);
  Real c309 = xlocal1 * xlocal9;
  Real c310 = -(xlocal4 * xlocal9);
  Real c311 = c274 + c276 + c277 + c308 + c309 + c310;
  Real c312 = xlocal11 * c311;
  Real c313 = c288 + c289 + c290 + c291 + c292 + c293 + c294 + c295 + c303 +
              c304 + c305 + c306 + c307 + c312;
  Real c314 = -(c287 * c313);
  Real c315 = ArcTan(c283, c314);
  Real c316 = c315 / 2.;
  Real c317 = Tan(c316);
  Real c320 = xlocal13 * xlocal2 * xlocal4 * xlocal5;
  Real c321 = -(xlocal1 * xlocal13 * c201);
  Real c322 = xlocal13 * xlocal3 * xlocal4 * xlocal6;
  Real c323 = -(xlocal1 * xlocal13 * c209);
  Real c324 = -(xlocal13 * xlocal2 * xlocal5 * xlocal7);
  Real c325 = xlocal1 * xlocal7 * c201;
  Real c326 = xlocal13 * xlocal7 * c201;
  Real c327 = -(xlocal13 * xlocal3 * xlocal6 * xlocal7);
  Real c328 = xlocal1 * xlocal7 * c209;
  Real c329 = xlocal13 * xlocal7 * c209;
  Real c330 = -(c201 * c61);
  Real c331 = -(c209 * c61);
  Real c332 = -(xlocal13 * xlocal2 * xlocal4 * xlocal8);
  Real c333 = 2 * xlocal1 * xlocal13 * xlocal5 * xlocal8;
  Real c334 = -(xlocal1 * xlocal4 * xlocal5 * xlocal8);
  Real c335 = -(xlocal13 * xlocal4 * xlocal5 * xlocal8);
  Real c336 = xlocal13 * xlocal2 * xlocal7 * xlocal8;
  Real c337 = -(xlocal1 * xlocal5 * xlocal7 * xlocal8);
  Real c338 = -(xlocal13 * xlocal5 * xlocal7 * xlocal8);
  Real c339 = 2 * xlocal4 * xlocal5 * xlocal7 * xlocal8;
  Real c340 = -(xlocal1 * xlocal13 * c79);
  Real c341 = xlocal1 * xlocal4 * c79;
  Real c342 = xlocal13 * xlocal4 * c79;
  Real c343 = -(c194 * c79);
  Real c344 = -(c209 * c79);
  Real c345 = -(xlocal13 * xlocal3 * xlocal4 * xlocal9);
  Real c346 = 2 * xlocal1 * xlocal13 * xlocal6 * xlocal9;
  Real c347 = -(xlocal1 * xlocal4 * xlocal6 * xlocal9);
  Real c348 = -(xlocal13 * xlocal4 * xlocal6 * xlocal9);
  Real c349 = xlocal13 * xlocal3 * xlocal7 * xlocal9;
  Real c350 = -(xlocal1 * xlocal6 * xlocal7 * xlocal9);
  Real c351 = -(xlocal13 * xlocal6 * xlocal7 * xlocal9);
  Real c352 = 2 * xlocal4 * xlocal6 * xlocal7 * xlocal9;
  Real c353 = 2 * xlocal5 * xlocal6 * xlocal8 * xlocal9;
  Real c354 = -(xlocal1 * xlocal13 * c100);
  Real c355 = xlocal1 * xlocal4 * c100;
  Real c356 = xlocal13 * xlocal4 * c100;
  Real c357 = -(c100 * c194);
  Real c358 = -(c100 * c201);
  Real c359 = -(xlocal2 * xlocal5 * xlocal6);
  Real c360 = xlocal4 * xlocal6 * xlocal7;
  Real c361 = -(xlocal6 * c61);
  Real c362 = xlocal2 * xlocal6 * xlocal8;
  Real c363 = xlocal5 * xlocal6 * xlocal8;
  Real c364 = -(xlocal6 * c79);
  Real c367 = c194 + c201 + c365 + c366 + c61 + c79;
  Real c368 = xlocal3 * c367;
  Real c369 = -(xlocal1 * c109 * c114);
  Real c370 = -(xlocal9 * c194);
  Real c371 = xlocal2 * xlocal5 * xlocal9;
  Real c372 = -(xlocal9 * c201);
  Real c373 = xlocal4 * xlocal7 * xlocal9;
  Real c374 = -(xlocal2 * xlocal8 * xlocal9);
  Real c375 = xlocal5 * xlocal8 * xlocal9;
  Real c376 = c359 + c360 + c361 + c362 + c363 + c364 + c368 + c369 + c370 +
              c371 + c372 + c373 + c374 + c375;
  Real c377 = -(xlocal15 * c376);
  Real c378 = -(xlocal3 * xlocal5 * xlocal6);
  Real c379 = xlocal4 * xlocal5 * xlocal7;
  Real c380 = -(xlocal5 * c61);
  Real c381 = -(xlocal1 * c106 * c114);
  Real c382 = -(xlocal8 * c194);
  Real c383 = xlocal3 * xlocal6 * xlocal8;
  Real c384 = -(xlocal8 * c209);
  Real c385 = xlocal4 * xlocal7 * xlocal8;
  Real c386 = xlocal3 * xlocal5 * xlocal9;
  Real c387 = xlocal5 * xlocal6 * xlocal9;
  Real c388 = -(xlocal3 * xlocal8 * xlocal9);
  Real c389 = xlocal6 * xlocal8 * xlocal9;
  Real c390 = -(xlocal5 * c100);
  Real c392 = c100 + c194 + c209 + c365 + c391 + c61;
  Real c393 = xlocal2 * c392;
  Real c394 = c378 + c379 + c380 + c381 + c382 + c383 + c384 + c385 + c386 +
              c387 + c388 + c389 + c390 + c393;
  Real c395 = -(xlocal14 * c394);
  Real c396 = c104 + c215 + c217 + c221 + c226 + c228 + c238 + c242 + c245 +
              c320 + c321 + c322 + c323 + c324 + c325 + c326 + c327 + c328 +
              c329 + c330 + c331 + c332 + c333 + c334 + c335 + c336 + c337 +
              c338 + c339 + c340 + c341 + c342 + c343 + c344 + c345 + c346 +
              c347 + c348 + c349 + c350 + c351 + c352 + c353 + c354 + c355 +
              c356 + c357 + c358 + c377 + c395 + c67 + c69 + c78 + c83 + c92 +
              c97 + c99;
  Real c397 = c100 + c194 + c201 + c209 + c365 + c366 + c391 + c61 + c79;
  Real c398 = xlocal13 * xlocal3 * xlocal5;
  Real c399 = -(xlocal13 * xlocal2 * xlocal6);
  Real c400 = -(xlocal13 * xlocal3 * xlocal8);
  Real c401 = xlocal13 * xlocal6 * xlocal8;
  Real c402 = xlocal15 * c173;
  Real c403 = xlocal13 * xlocal2 * xlocal9;
  Real c404 = -(xlocal13 * xlocal5 * xlocal9);
  Real c405 = xlocal14 * c185;
  Real c406 = c163 + c164 + c166 + c167 + c176 + c177 + c398 + c399 + c400 +
              c401 + c402 + c403 + c404 + c405;
  Real c407 = c397 * c406;
  Real c408 = ArcTan(c396, c407);
  Real c409 = c408 / 2.;
  Real c410 = Tan(c409);
  out(0) = c16 * c16 + c25 * c25 + c7 * c7;
  out(1) = c16 * c34 + c25 * c38 + c30 * c7;
  out(2) = c30 * c30 + c34 * c34 + c38 * c38;
  out(3) = (l0 * l1 * t21 * t21 * c191 + l1 * l2 * t01 * t01 * c317 +
            l0 * l2 * t11 * t11 * c410) *
           c46 * c47 * c48 * c49;
  out(4) = (l0 * l1 * t21 * t22 * c191 + l1 * l2 * t01 * t02 * c317 +
            l0 * l2 * t11 * t12 * c410) *
           c46 * c47 * c48 * c49;
  out(5) = (l0 * l1 * t22 * t22 * c191 + l1 * l2 * t02 * t02 * c317 +
            l0 * l2 * t12 * t12 * c410) *
           c46 * c47 * c48 * c49;

  return _out;
}