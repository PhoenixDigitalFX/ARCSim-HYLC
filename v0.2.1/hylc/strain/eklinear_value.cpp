#include "strain.hpp"

using namespace hylc;
using namespace hylc::mathematica;

Vec6 hylc::mathematica::eklinear(const Vec18 &xloc, const Mat2x2 &invDm,
                                 const Real &A, const Real &thetarest0,
                                 const Real &thetarest1, const Real &thetarest2,
                                 const Real &l0, const Real &l1, const Real &l2,
                                 const Vec2 &t0, const Vec2 &t1,
                                 const Vec2 &t2) {
  // define output
  Vec6 out;

  Real copt1 = invDm(0, 0);
  Real copt2 = xloc(0);
  Real copt3 = -copt2;
  Real copt4 = xloc(3);
  Real copt5 = copt3 + copt4;
  Real copt6 = copt1 * copt5;
  Real copt7 = invDm(1, 0);
  Real copt8 = xloc(6);
  Real copt9 = copt3 + copt8;
  Real copt10 = copt7 * copt9;
  Real copt11 = copt10 + copt6;
  Real copt12 = Power(copt11, 2);
  Real copt13 = xloc(1);
  Real copt14 = -copt13;
  Real copt15 = xloc(4);
  Real copt16 = copt14 + copt15;
  Real copt17 = copt1 * copt16;
  Real copt18 = xloc(7);
  Real copt19 = copt14 + copt18;
  Real copt20 = copt19 * copt7;
  Real copt21 = copt17 + copt20;
  Real copt22 = Power(copt21, 2);
  Real copt23 = xloc(2);
  Real copt24 = -copt23;
  Real copt25 = xloc(5);
  Real copt26 = copt24 + copt25;
  Real copt27 = copt1 * copt26;
  Real copt28 = xloc(8);
  Real copt29 = copt24 + copt28;
  Real copt30 = copt29 * copt7;
  Real copt31 = copt27 + copt30;
  Real copt32 = Power(copt31, 2);
  Real copt34 = invDm(0, 1);
  Real copt36 = invDm(1, 1);
  Real copt35 = copt34 * copt5;
  Real copt37 = copt36 * copt9;
  Real copt38 = copt35 + copt37;
  Real copt49 = Power(copt38, 2);
  Real copt40 = copt16 * copt34;
  Real copt41 = copt19 * copt36;
  Real copt42 = copt40 + copt41;
  Real copt50 = Power(copt42, 2);
  Real copt44 = copt26 * copt34;
  Real copt45 = copt29 * copt36;
  Real copt46 = copt44 + copt45;
  Real copt51 = Power(copt46, 2);
  Real copt53 = 1 / A;
  Real copt54 = 1 / l0;
  Real copt55 = 1 / l1;
  Real copt56 = 1 / l2;
  Real copt57 = t0(0);
  Real copt58 = Power(copt57, 2);
  Real copt60 = -copt4;
  Real copt61 = copt2 + copt60;
  Real copt62 = Power(copt61, 2);
  Real copt63 = -copt15;
  Real copt64 = copt13 + copt63;
  Real copt65 = Power(copt64, 2);
  Real copt66 = -copt25;
  Real copt67 = copt23 + copt66;
  Real copt68 = Power(copt67, 2);
  Real copt69 = copt62 + copt65 + copt68;
  Real copt70 = Sqrt(copt69);
  Real copt71 = Power(copt2, 2);
  Real copt72 = Power(copt15, 2);
  Real copt74 = Power(copt25, 2);
  Real copt82 = xloc(9);
  Real copt91 = xloc(10);
  Real copt98 = Power(copt4, 2);
  Real copt102 = Power(copt23, 2);
  Real copt111 = xloc(11);
  Real copt121 = Power(copt13, 2);
  Real copt103 = copt8 * copt82;
  Real copt104 = copt8 + copt82;
  Real copt105 = -(copt104 * copt4);
  Real copt123 = copt111 + copt28;
  Real copt137 = copt104 * copt25;
  Real copt107 = copt18 + copt91;
  Real copt143 = -2 * copt25;
  Real copt144 = copt111 + copt143 + copt28;
  Real copt161 = -(copt18 * copt82);
  Real copt136 = copt28 * copt82;
  Real copt142 = -(copt111 * copt8);
  Real copt219 = t1(0);
  Real copt220 = Power(copt219, 2);
  Real copt194 = -copt8;
  Real copt222 = copt194 + copt4;
  Real copt223 = Power(copt222, 2);
  Real copt224 = -copt18;
  Real copt225 = copt15 + copt224;
  Real copt226 = Power(copt225, 2);
  Real copt210 = -copt28;
  Real copt227 = copt210 + copt25;
  Real copt228 = Power(copt227, 2);
  Real copt229 = copt223 + copt226 + copt228;
  Real copt230 = Sqrt(copt229);
  Real copt232 = Power(copt8, 2);
  Real copt238 = Power(copt18, 2);
  Real copt248 = Power(copt28, 2);
  Real copt251 = xloc(12);
  Real copt264 = xloc(13);
  Real copt277 = xloc(14);
  Real copt295 = -copt251;
  Real copt296 = copt295 + copt8;
  Real copt304 = -copt264;
  Real copt305 = copt18 + copt304;
  Real copt312 = -copt277;
  Real copt313 = copt28 + copt312;
  Real copt186 = -(copt18 * copt2 * copt25);
  Real copt187 = copt15 * copt2 * copt28;
  Real copt330 = copt305 * copt4;
  Real copt378 = t2(0);
  Real copt379 = Power(copt378, 2);
  Real copt381 = copt194 + copt2;
  Real copt382 = Power(copt381, 2);
  Real copt383 = copt13 + copt224;
  Real copt384 = Power(copt383, 2);
  Real copt385 = copt210 + copt23;
  Real copt386 = Power(copt385, 2);
  Real copt387 = copt382 + copt384 + copt386;
  Real copt388 = Sqrt(copt387);
  Real copt397 = xloc(15);
  Real copt411 = xloc(16);
  Real copt406 = -copt232;
  Real copt407 = -copt397;
  Real copt408 = copt407 + copt8;
  Real copt409 = copt4 * copt408;
  Real copt410 = copt397 * copt8;
  Real copt426 = xloc(17);
  Real copt412 = -copt411;
  Real copt413 = copt18 + copt412;
  Real copt427 = -copt426;
  Real copt428 = copt28 + copt427;
  Real copt458 = copt4 * copt413;
  Real copt459 = copt411 * copt8;
  Real copt483 = copt25 * copt408;
  Real copt484 = copt28 * copt397;
  Real copt73 = copt71 * copt72;
  Real copt75 = copt71 * copt74;
  Real copt76 = -(copt2 * copt72 * copt8);
  Real copt77 = -(copt2 * copt74 * copt8);
  Real copt78 = -(copt15 * copt18 * copt71);
  Real copt79 = copt15 * copt18 * copt2 * copt4;
  Real copt80 = -(copt25 * copt28 * copt71);
  Real copt81 = copt2 * copt25 * copt28 * copt4;
  Real copt83 = -(copt2 * copt72 * copt82);
  Real copt84 = -(copt2 * copt74 * copt82);
  Real copt85 = copt72 * copt8 * copt82;
  Real copt86 = copt74 * copt8 * copt82;
  Real copt87 = copt15 * copt18 * copt2 * copt82;
  Real copt88 = -(copt15 * copt18 * copt4 * copt82);
  Real copt89 = copt2 * copt25 * copt28 * copt82;
  Real copt90 = -(copt25 * copt28 * copt4 * copt82);
  Real copt92 = -(copt15 * copt71 * copt91);
  Real copt93 = copt15 * copt2 * copt4 * copt91;
  Real copt94 = copt15 * copt2 * copt8 * copt91;
  Real copt95 = -(copt15 * copt4 * copt8 * copt91);
  Real copt96 = copt18 * copt71 * copt91;
  Real copt97 = -2 * copt18 * copt2 * copt4 * copt91;
  Real copt99 = copt18 * copt91 * copt98;
  Real copt100 = copt18 * copt74 * copt91;
  Real copt101 = -(copt15 * copt25 * copt28 * copt91);
  Real copt106 = copt18 * copt91;
  Real copt108 = -(copt107 * copt15);
  Real copt109 = copt103 + copt105 + copt106 + copt108 + copt72 + copt98;
  Real copt110 = copt102 * copt109;
  Real copt112 = -(copt111 * copt25 * copt71);
  Real copt113 = copt111 * copt2 * copt25 * copt4;
  Real copt114 = copt111 * copt2 * copt25 * copt8;
  Real copt115 = -(copt111 * copt25 * copt4 * copt8);
  Real copt116 = -(copt111 * copt15 * copt18 * copt25);
  Real copt117 = copt111 * copt28 * copt71;
  Real copt118 = -2 * copt111 * copt2 * copt28 * copt4;
  Real copt119 = copt111 * copt28 * copt98;
  Real copt120 = copt111 * copt28 * copt72;
  Real copt122 = copt111 * copt28;
  Real copt124 = -(copt123 * copt25);
  Real copt125 = copt103 + copt105 + copt122 + copt124 + copt74 + copt98;
  Real copt126 = copt121 * copt125;
  Real copt127 = copt15 * copt18 * copt25;
  Real copt128 = -(copt28 * copt72);
  Real copt129 = -2 * copt25 * copt8 * copt82;
  Real copt130 = copt15 * copt25 * copt91;
  Real copt131 = -2 * copt18 * copt25 * copt91;
  Real copt132 = copt15 * copt28 * copt91;
  Real copt133 = -(copt111 * copt72);
  Real copt134 = copt111 * copt15 * copt18;
  Real copt135 = -(copt123 * copt98);
  Real copt138 = copt111 * copt8;
  Real copt139 = copt136 + copt137 + copt138;
  Real copt140 = copt139 * copt4;
  Real copt141 = -(copt28 * copt82);
  Real copt145 = copt144 * copt4;
  Real copt146 = copt137 + copt141 + copt142 + copt145;
  Real copt147 = copt146 * copt2;
  Real copt148 = copt127 + copt128 + copt129 + copt130 + copt131 + copt132 +
                 copt133 + copt134 + copt135 + copt140 + copt147;
  Real copt149 = copt148 * copt23;
  Real copt150 = copt15 * copt4 * copt8;
  Real copt151 = -(copt18 * copt98);
  Real copt152 = -(copt18 * copt74);
  Real copt153 = copt15 * copt25 * copt28;
  Real copt154 = copt15 * copt4 * copt82;
  Real copt155 = -2 * copt15 * copt8 * copt82;
  Real copt156 = copt18 * copt4 * copt82;
  Real copt157 = -(copt91 * copt98);
  Real copt158 = -(copt74 * copt91);
  Real copt159 = copt4 * copt8 * copt91;
  Real copt160 = copt25 * copt28 * copt91;
  Real copt162 = copt104 * copt15;
  Real copt163 = -(copt8 * copt91);
  Real copt164 = -2 * copt15;
  Real copt165 = copt164 + copt18 + copt91;
  Real copt166 = copt165 * copt4;
  Real copt167 = copt161 + copt162 + copt163 + copt166;
  Real copt168 = copt167 * copt2;
  Real copt169 = copt111 * copt15 * copt25;
  Real copt170 = copt111 * copt18 * copt25;
  Real copt171 = -2 * copt111 * copt15 * copt28;
  Real copt172 = -(copt28 * copt91);
  Real copt173 = copt107 * copt25;
  Real copt174 = -(copt111 * copt18);
  Real copt175 = copt144 * copt15;
  Real copt176 = copt172 + copt173 + copt174 + copt175;
  Real copt177 = copt176 * copt23;
  Real copt178 = copt150 + copt151 + copt152 + copt153 + copt154 + copt155 +
                 copt156 + copt157 + copt158 + copt159 + copt160 + copt168 +
                 copt169 + copt170 + copt171 + copt177;
  Real copt179 = copt13 * copt178;
  Real copt180 = copt100 + copt101 + copt110 + copt112 + copt113 + copt114 +
                 copt115 + copt116 + copt117 + copt118 + copt119 + copt120 +
                 copt126 + copt149 + copt179 + copt73 + copt75 + copt76 +
                 copt77 + copt78 + copt79 + copt80 + copt81 + copt83 + copt84 +
                 copt85 + copt86 + copt87 + copt88 + copt89 + copt90 + copt92 +
                 copt93 + copt94 + copt95 + copt96 + copt97 + copt99;
  Real copt181 = -(copt180 * copt70);
  Real copt182 = -2 * copt2 * copt4;
  Real copt183 = -2 * copt13 * copt15;
  Real copt184 = -2 * copt23 * copt25;
  Real copt185 = copt102 + copt121 + copt182 + copt183 + copt184 + copt71 +
                 copt72 + copt74 + copt98;
  Real copt188 = copt18 * copt25 * copt82;
  Real copt189 = -(copt15 * copt28 * copt82);
  Real copt190 = copt2 * copt25 * copt91;
  Real copt191 = -(copt25 * copt8 * copt91);
  Real copt192 = -(copt2 * copt28 * copt91);
  Real copt193 = copt28 * copt4 * copt91;
  Real copt195 = copt194 + copt82;
  Real copt196 = copt15 * copt195;
  Real copt197 = -copt91;
  Real copt198 = copt18 + copt197;
  Real copt199 = copt198 * copt4;
  Real copt200 = copt8 * copt91;
  Real copt201 = copt161 + copt196 + copt199 + copt200;
  Real copt202 = copt201 * copt23;
  Real copt203 = -(copt111 * copt15 * copt2);
  Real copt204 = copt111 * copt15 * copt8;
  Real copt205 = copt111 * copt18 * copt2;
  Real copt206 = -(copt111 * copt18 * copt4);
  Real copt207 = -copt82;
  Real copt208 = copt207 + copt8;
  Real copt209 = copt208 * copt25;
  Real copt211 = copt111 + copt210;
  Real copt212 = copt211 * copt4;
  Real copt213 = copt136 + copt142 + copt209 + copt212;
  Real copt214 = copt13 * copt213;
  Real copt215 = copt186 + copt187 + copt188 + copt189 + copt190 + copt191 +
                 copt192 + copt193 + copt202 + copt203 + copt204 + copt205 +
                 copt206 + copt214;
  Real copt216 = copt185 * copt215;
  Real copt217 = ArcTan(copt181, copt216);
  Real copt523 = t0(1);
  Real copt231 = -(copt23 * copt25 * copt4 * copt8);
  Real copt233 = -(copt232 * copt72);
  Real copt234 = copt23 * copt232 * copt25;
  Real copt235 = -(copt232 * copt74);
  Real copt236 = -(copt15 * copt18 * copt23 * copt25);
  Real copt237 = 2 * copt15 * copt18 * copt4 * copt8;
  Real copt239 = -(copt238 * copt98);
  Real copt240 = copt23 * copt238 * copt25;
  Real copt241 = -(copt238 * copt74);
  Real copt242 = copt23 * copt28 * copt98;
  Real copt243 = copt23 * copt28 * copt72;
  Real copt244 = -(copt23 * copt28 * copt4 * copt8);
  Real copt245 = 2 * copt25 * copt28 * copt4 * copt8;
  Real copt246 = -(copt15 * copt18 * copt23 * copt28);
  Real copt247 = 2 * copt15 * copt18 * copt25 * copt28;
  Real copt249 = -(copt248 * copt98);
  Real copt250 = -(copt248 * copt72);
  Real copt252 = copt23 * copt25 * copt251 * copt4;
  Real copt253 = copt251 * copt72 * copt8;
  Real copt254 = -(copt23 * copt25 * copt251 * copt8);
  Real copt255 = copt251 * copt74 * copt8;
  Real copt256 = -(copt15 * copt18 * copt251 * copt4);
  Real copt257 = -(copt15 * copt18 * copt251 * copt8);
  Real copt258 = copt238 * copt251 * copt4;
  Real copt259 = -(copt23 * copt251 * copt28 * copt4);
  Real copt260 = -(copt25 * copt251 * copt28 * copt4);
  Real copt261 = copt23 * copt251 * copt28 * copt8;
  Real copt262 = -(copt25 * copt251 * copt28 * copt8);
  Real copt263 = copt248 * copt251 * copt4;
  Real copt265 = copt15 * copt23 * copt25 * copt264;
  Real copt266 = -(copt15 * copt264 * copt4 * copt8);
  Real copt267 = copt15 * copt232 * copt264;
  Real copt268 = copt18 * copt264 * copt98;
  Real copt269 = -(copt18 * copt23 * copt25 * copt264);
  Real copt270 = copt18 * copt264 * copt74;
  Real copt271 = -(copt18 * copt264 * copt4 * copt8);
  Real copt272 = -(copt15 * copt23 * copt264 * copt28);
  Real copt273 = -(copt15 * copt25 * copt264 * copt28);
  Real copt274 = copt18 * copt23 * copt264 * copt28;
  Real copt275 = -(copt18 * copt25 * copt264 * copt28);
  Real copt276 = copt15 * copt248 * copt264;
  Real copt278 = -(copt23 * copt277 * copt98);
  Real copt279 = -(copt23 * copt277 * copt72);
  Real copt280 = 2 * copt23 * copt277 * copt4 * copt8;
  Real copt281 = -(copt25 * copt277 * copt4 * copt8);
  Real copt282 = -(copt23 * copt232 * copt277);
  Real copt283 = copt232 * copt25 * copt277;
  Real copt284 = 2 * copt15 * copt18 * copt23 * copt277;
  Real copt285 = -(copt15 * copt18 * copt25 * copt277);
  Real copt286 = -(copt23 * copt238 * copt277);
  Real copt287 = copt238 * copt25 * copt277;
  Real copt288 = copt277 * copt28 * copt98;
  Real copt289 = copt277 * copt28 * copt72;
  Real copt290 = -(copt277 * copt28 * copt4 * copt8);
  Real copt291 = -(copt15 * copt18 * copt277 * copt28);
  Real copt292 = copt18 * copt74;
  Real copt293 = -(copt18 * copt25 * copt28);
  Real copt294 = copt18 * copt251 * copt8;
  Real copt297 = copt15 * copt296;
  Real copt298 = copt18 * copt251;
  Real copt299 = -2 * copt264;
  Real copt300 = copt18 + copt299;
  Real copt301 = copt300 * copt8;
  Real copt302 = copt297 + copt298 + copt301;
  Real copt303 = -(copt302 * copt4);
  Real copt306 = copt305 * copt98;
  Real copt307 = -(copt264 * copt74);
  Real copt308 = -(copt232 * copt264);
  Real copt309 = 2 * copt25 * copt264 * copt28;
  Real copt310 = -(copt248 * copt264);
  Real copt311 = -(copt251 * copt8);
  Real copt314 = -(copt227 * copt313);
  Real copt315 = copt232 + copt311 + copt314;
  Real copt316 = copt15 * copt315;
  Real copt317 = -(copt18 * copt25 * copt277);
  Real copt318 = copt18 * copt277 * copt28;
  Real copt319 = copt292 + copt293 + copt294 + copt303 + copt306 + copt307 +
                 copt308 + copt309 + copt310 + copt316 + copt317 + copt318;
  Real copt320 = copt13 * copt319;
  Real copt321 = copt238 * copt4;
  Real copt322 = copt248 * copt4;
  Real copt323 = copt296 * copt72;
  Real copt324 = copt296 * copt74;
  Real copt325 = -(copt238 * copt251);
  Real copt326 = -(copt248 * copt251);
  Real copt327 = -(copt18 * copt264 * copt4);
  Real copt328 = copt18 * copt264 * copt8;
  Real copt329 = -2 * copt18 * copt251;
  Real copt331 = copt18 + copt264;
  Real copt332 = copt331 * copt8;
  Real copt333 = copt329 + copt330 + copt332;
  Real copt334 = -(copt15 * copt333);
  Real copt335 = -(copt277 * copt28 * copt4);
  Real copt336 = copt277 * copt28 * copt8;
  Real copt337 = -2 * copt251 * copt28;
  Real copt338 = copt313 * copt4;
  Real copt339 = copt277 + copt28;
  Real copt340 = copt339 * copt8;
  Real copt341 = copt337 + copt338 + copt340;
  Real copt342 = -(copt25 * copt341);
  Real copt343 = copt321 + copt322 + copt323 + copt324 + copt325 + copt326 +
                 copt327 + copt328 + copt334 + copt335 + copt336 + copt342;
  Real copt344 = copt2 * copt343;
  Real copt345 = copt231 + copt233 + copt234 + copt235 + copt236 + copt237 +
                 copt239 + copt240 + copt241 + copt242 + copt243 + copt244 +
                 copt245 + copt246 + copt247 + copt249 + copt250 + copt252 +
                 copt253 + copt254 + copt255 + copt256 + copt257 + copt258 +
                 copt259 + copt260 + copt261 + copt262 + copt263 + copt265 +
                 copt266 + copt267 + copt268 + copt269 + copt270 + copt271 +
                 copt272 + copt273 + copt274 + copt275 + copt276 + copt278 +
                 copt279 + copt280 + copt281 + copt282 + copt283 + copt284 +
                 copt285 + copt286 + copt287 + copt288 + copt289 + copt290 +
                 copt291 + copt320 + copt344;
  Real copt346 = copt230 * copt345;
  Real copt347 = -2 * copt4 * copt8;
  Real copt348 = -2 * copt15 * copt18;
  Real copt349 = -2 * copt25 * copt28;
  Real copt350 = copt232 + copt238 + copt248 + copt347 + copt348 + copt349 +
                 copt72 + copt74 + copt98;
  Real copt351 = copt18 * copt25 * copt251;
  Real copt352 = -(copt15 * copt251 * copt28);
  Real copt353 = copt2 * copt25 * copt264;
  Real copt354 = -(copt25 * copt264 * copt8);
  Real copt355 = -(copt2 * copt264 * copt28);
  Real copt356 = copt264 * copt28 * copt4;
  Real copt357 = -(copt18 * copt251);
  Real copt358 = copt194 + copt251;
  Real copt359 = copt15 * copt358;
  Real copt360 = copt264 * copt8;
  Real copt361 = copt330 + copt357 + copt359 + copt360;
  Real copt362 = copt23 * copt361;
  Real copt363 = -(copt15 * copt2 * copt277);
  Real copt364 = copt15 * copt277 * copt8;
  Real copt365 = copt18 * copt2 * copt277;
  Real copt366 = -(copt18 * copt277 * copt4);
  Real copt367 = copt25 * copt296;
  Real copt368 = copt251 * copt28;
  Real copt369 = -(copt277 * copt8);
  Real copt370 = copt210 + copt277;
  Real copt371 = copt370 * copt4;
  Real copt372 = copt367 + copt368 + copt369 + copt371;
  Real copt373 = copt13 * copt372;
  Real copt374 = copt186 + copt187 + copt351 + copt352 + copt353 + copt354 +
                 copt355 + copt356 + copt362 + copt363 + copt364 + copt365 +
                 copt366 + copt373;
  Real copt375 = copt350 * copt374;
  Real copt376 = ArcTan(copt346, copt375);
  Real copt526 = t1(1);
  Real copt389 = copt15 * copt18 * copt71;
  Real copt390 = -(copt15 * copt18 * copt2 * copt8);
  Real copt391 = -(copt238 * copt71);
  Real copt392 = copt2 * copt238 * copt4;
  Real copt393 = copt25 * copt28 * copt71;
  Real copt394 = -(copt2 * copt25 * copt28 * copt8);
  Real copt395 = -(copt248 * copt71);
  Real copt396 = copt2 * copt248 * copt4;
  Real copt398 = -(copt15 * copt18 * copt2 * copt397);
  Real copt399 = copt15 * copt18 * copt397 * copt8;
  Real copt400 = copt2 * copt238 * copt397;
  Real copt401 = -(copt238 * copt397 * copt4);
  Real copt402 = -(copt2 * copt25 * copt28 * copt397);
  Real copt403 = copt25 * copt28 * copt397 * copt8;
  Real copt404 = copt2 * copt248 * copt397;
  Real copt405 = -(copt248 * copt397 * copt4);
  Real copt414 = copt225 * copt413;
  Real copt415 = copt406 + copt409 + copt410 + copt414;
  Real copt416 = copt102 * copt415;
  Real copt417 = -(copt15 * copt411 * copt71);
  Real copt418 = 2 * copt15 * copt2 * copt411 * copt8;
  Real copt419 = -(copt15 * copt232 * copt411);
  Real copt420 = copt18 * copt411 * copt71;
  Real copt421 = -(copt18 * copt2 * copt4 * copt411);
  Real copt422 = -(copt18 * copt2 * copt411 * copt8);
  Real copt423 = copt18 * copt4 * copt411 * copt8;
  Real copt424 = copt18 * copt25 * copt28 * copt411;
  Real copt425 = -(copt15 * copt248 * copt411);
  Real copt429 = copt227 * copt428;
  Real copt430 = copt406 + copt409 + copt410 + copt429;
  Real copt431 = copt121 * copt430;
  Real copt432 = -(copt25 * copt426 * copt71);
  Real copt433 = 2 * copt2 * copt25 * copt426 * copt8;
  Real copt434 = -(copt232 * copt25 * copt426);
  Real copt435 = -(copt238 * copt25 * copt426);
  Real copt436 = copt28 * copt426 * copt71;
  Real copt437 = -(copt2 * copt28 * copt4 * copt426);
  Real copt438 = -(copt2 * copt28 * copt426 * copt8);
  Real copt439 = copt28 * copt4 * copt426 * copt8;
  Real copt440 = copt15 * copt18 * copt28 * copt426;
  Real copt441 = -(copt15 * copt232);
  Real copt442 = copt18 * copt23 * copt25;
  Real copt443 = copt18 * copt4 * copt8;
  Real copt444 = -2 * copt18 * copt23 * copt28;
  Real copt445 = copt18 * copt25 * copt28;
  Real copt446 = copt15 * copt397 * copt8;
  Real copt447 = -2 * copt18 * copt397 * copt4;
  Real copt448 = copt18 * copt397 * copt8;
  Real copt449 = -(copt23 * copt25 * copt411);
  Real copt450 = copt4 * copt411 * copt8;
  Real copt451 = -(copt232 * copt411);
  Real copt452 = copt23 * copt28 * copt411;
  Real copt453 = copt25 * copt28 * copt411;
  Real copt454 = -(copt248 * copt411);
  Real copt455 = -2 * copt18 * copt8;
  Real copt456 = copt15 * copt408;
  Real copt457 = copt18 * copt397;
  Real copt460 = copt455 + copt456 + copt457 + copt458 + copt459;
  Real copt461 = copt2 * copt460;
  Real copt462 = copt15 * copt385 * copt428;
  Real copt463 = copt18 * copt23 * copt426;
  Real copt464 = -2 * copt18 * copt25 * copt426;
  Real copt465 = copt18 * copt28 * copt426;
  Real copt466 = copt441 + copt442 + copt443 + copt444 + copt445 + copt446 +
                 copt447 + copt448 + copt449 + copt450 + copt451 + copt452 +
                 copt453 + copt454 + copt461 + copt462 + copt463 + copt464 +
                 copt465;
  Real copt467 = -(copt13 * copt466);
  Real copt468 = copt28 * copt4 * copt8;
  Real copt469 = copt15 * copt18 * copt28;
  Real copt470 = -2 * copt28 * copt397 * copt4;
  Real copt471 = copt28 * copt397 * copt8;
  Real copt472 = -2 * copt15 * copt28 * copt411;
  Real copt473 = copt18 * copt28 * copt411;
  Real copt474 = -(copt397 * copt8);
  Real copt475 = -(copt18 * copt411);
  Real copt476 = copt232 + copt238 + copt474 + copt475;
  Real copt477 = -(copt25 * copt476);
  Real copt478 = copt4 * copt426 * copt8;
  Real copt479 = -(copt232 * copt426);
  Real copt480 = copt15 * copt18 * copt426;
  Real copt481 = -(copt238 * copt426);
  Real copt482 = -2 * copt28 * copt8;
  Real copt485 = copt4 * copt428;
  Real copt486 = copt426 * copt8;
  Real copt487 = copt482 + copt483 + copt484 + copt485 + copt486;
  Real copt488 = copt2 * copt487;
  Real copt489 = copt468 + copt469 + copt470 + copt471 + copt472 + copt473 +
                 copt477 + copt478 + copt479 + copt480 + copt481 + copt488;
  Real copt490 = -(copt23 * copt489);
  Real copt491 = copt389 + copt390 + copt391 + copt392 + copt393 + copt394 +
                 copt395 + copt396 + copt398 + copt399 + copt400 + copt401 +
                 copt402 + copt403 + copt404 + copt405 + copt416 + copt417 +
                 copt418 + copt419 + copt420 + copt421 + copt422 + copt423 +
                 copt424 + copt425 + copt431 + copt432 + copt433 + copt434 +
                 copt435 + copt436 + copt437 + copt438 + copt439 + copt440 +
                 copt467 + copt490;
  Real copt492 = copt388 * copt491;
  Real copt493 = -2 * copt2 * copt8;
  Real copt494 = -2 * copt13 * copt18;
  Real copt495 = -2 * copt23 * copt28;
  Real copt496 = copt102 + copt121 + copt232 + copt238 + copt248 + copt493 +
                 copt494 + copt495 + copt71;
  Real copt497 = copt18 * copt25 * copt397;
  Real copt498 = -(copt15 * copt28 * copt397);
  Real copt499 = copt2 * copt25 * copt411;
  Real copt500 = -(copt25 * copt411 * copt8);
  Real copt501 = -(copt2 * copt28 * copt411);
  Real copt502 = copt28 * copt4 * copt411;
  Real copt503 = -(copt18 * copt397);
  Real copt504 = copt194 + copt397;
  Real copt505 = copt15 * copt504;
  Real copt506 = copt458 + copt459 + copt503 + copt505;
  Real copt507 = copt23 * copt506;
  Real copt508 = -(copt15 * copt2 * copt426);
  Real copt509 = copt15 * copt426 * copt8;
  Real copt510 = copt18 * copt2 * copt426;
  Real copt511 = -(copt18 * copt4 * copt426);
  Real copt512 = -(copt426 * copt8);
  Real copt513 = copt210 + copt426;
  Real copt514 = copt4 * copt513;
  Real copt515 = copt483 + copt484 + copt512 + copt514;
  Real copt516 = copt13 * copt515;
  Real copt517 = copt186 + copt187 + copt497 + copt498 + copt499 + copt500 +
                 copt501 + copt502 + copt507 + copt508 + copt509 + copt510 +
                 copt511 + copt516;
  Real copt518 = copt496 * copt517;
  Real copt519 = ArcTan(copt492, copt518);
  Real copt529 = t2(1);
  Real copt534 = Power(copt523, 2);
  Real copt537 = Power(copt526, 2);
  Real copt540 = Power(copt529, 2);
  out(0) = copt12 + copt22 + copt32;
  out(1) = copt11 * copt38 + copt21 * copt42 + copt31 * copt46;
  out(2) = copt49 + copt50 + copt51;
  out(3) = (copt53 * copt54 * copt55 * copt56 *
            (copt379 * copt519 * l0 * l1 + copt220 * copt376 * l0 * l2 +
             copt217 * copt58 * l1 * l2 - copt58 * l1 * l2 * thetarest0 -
             copt220 * l0 * l2 * thetarest1 - copt379 * l0 * l1 * thetarest2)) /
           2.;
  out(4) = (copt53 * copt54 * copt55 * copt56 *
            (copt378 * copt519 * copt529 * l0 * l1 +
             copt219 * copt376 * copt526 * l0 * l2 +
             copt217 * copt523 * copt57 * l1 * l2 -
             copt523 * copt57 * l1 * l2 * thetarest0 -
             copt219 * copt526 * l0 * l2 * thetarest1 -
             copt378 * copt529 * l0 * l1 * thetarest2)) /
           2.;
  out(5) = (copt53 * copt54 * copt55 * copt56 *
            (copt519 * copt540 * l0 * l1 + copt376 * copt537 * l0 * l2 +
             copt217 * copt534 * l1 * l2 - copt534 * l1 * l2 * thetarest0 -
             copt537 * l0 * l2 * thetarest1 - copt540 * l0 * l1 * thetarest2)) /
           2.;

  return out;
}
