#include "strain.hpp"

using namespace hylc;
using namespace hylc::mathematica;

std::tuple<Mat6x18, Vec6> hylc::mathematica::eklinear_valgrad(
    const Vec18 &xloc, const Mat2x2 &invDm, const Real &A,
    const Real &thetarest0, const Real &thetarest1, const Real &thetarest2,
    const Real &l0, const Real &l1, const Real &l2, const Vec2 &t0,
    const Vec2 &t1, const Vec2 &t2) {
  // define output
  Mat6x18 grad(0);
  Vec6 val(0);
  auto out1 = [&](int i) -> Real & { return val[i]; };
  auto out2 = [&](int i, int j) -> Real & { return grad(i, j); };

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
  Real copt16 = -copt13;
  Real copt17 = xloc(4);
  Real copt18 = copt16 + copt17;
  Real copt19 = copt1 * copt18;
  Real copt20 = xloc(7);
  Real copt21 = copt16 + copt20;
  Real copt22 = copt21 * copt7;
  Real copt23 = copt19 + copt22;
  Real copt24 = Power(copt23, 2);
  Real copt26 = xloc(2);
  Real copt28 = -copt26;
  Real copt29 = xloc(5);
  Real copt30 = copt28 + copt29;
  Real copt31 = copt1 * copt30;
  Real copt32 = xloc(8);
  Real copt34 = copt28 + copt32;
  Real copt35 = copt34 * copt7;
  Real copt36 = copt31 + copt35;
  Real copt38 = Power(copt36, 2);
  Real copt42 = invDm(0, 1);
  Real copt45 = invDm(1, 1);
  Real copt44 = copt42 * copt5;
  Real copt49 = copt45 * copt9;
  Real copt50 = copt44 + copt49;
  Real copt70 = Power(copt50, 2);
  Real copt54 = copt18 * copt42;
  Real copt55 = copt21 * copt45;
  Real copt57 = copt54 + copt55;
  Real copt72 = Power(copt57, 2);
  Real copt60 = copt30 * copt42;
  Real copt61 = copt34 * copt45;
  Real copt63 = copt60 + copt61;
  Real copt73 = Power(copt63, 2);
  Real copt85 = 1 / A;
  Real copt86 = 1 / l0;
  Real copt120 = 1 / l1;
  Real copt121 = 1 / l2;
  Real copt123 = t0(0);
  Real copt128 = Power(copt123, 2);
  Real copt134 = -copt4;
  Real copt138 = copt134 + copt2;
  Real copt139 = Power(copt138, 2);
  Real copt140 = -copt17;
  Real copt173 = copt13 + copt140;
  Real copt174 = Power(copt173, 2);
  Real copt175 = -copt29;
  Real copt177 = copt175 + copt26;
  Real copt178 = Power(copt177, 2);
  Real copt179 = copt139 + copt174 + copt178;
  Real copt181 = Sqrt(copt179);
  Real copt182 = Power(copt2, 2);
  Real copt183 = Power(copt17, 2);
  Real copt189 = Power(copt29, 2);
  Real copt197 = xloc(9);
  Real copt206 = xloc(10);
  Real copt213 = Power(copt4, 2);
  Real copt217 = Power(copt26, 2);
  Real copt228 = xloc(11);
  Real copt238 = Power(copt13, 2);
  Real copt219 = copt197 * copt8;
  Real copt220 = copt197 + copt8;
  Real copt222 = -(copt220 * copt4);
  Real copt240 = copt228 + copt32;
  Real copt254 = copt220 * copt29;
  Real copt224 = copt20 + copt206;
  Real copt260 = -2 * copt29;
  Real copt261 = copt228 + copt260 + copt32;
  Real copt278 = -(copt197 * copt20);
  Real copt253 = copt197 * copt32;
  Real copt259 = -(copt228 * copt8);
  Real copt336 = t1(0);
  Real copt337 = Power(copt336, 2);
  Real copt311 = -copt8;
  Real copt340 = copt311 + copt4;
  Real copt341 = Power(copt340, 2);
  Real copt342 = -copt20;
  Real copt343 = copt17 + copt342;
  Real copt344 = Power(copt343, 2);
  Real copt327 = -copt32;
  Real copt345 = copt29 + copt327;
  Real copt346 = Power(copt345, 2);
  Real copt347 = copt341 + copt344 + copt346;
  Real copt348 = Sqrt(copt347);
  Real copt350 = Power(copt8, 2);
  Real copt356 = Power(copt20, 2);
  Real copt366 = Power(copt32, 2);
  Real copt369 = xloc(12);
  Real copt384 = xloc(13);
  Real copt398 = xloc(14);
  Real copt416 = -copt369;
  Real copt417 = copt416 + copt8;
  Real copt425 = -copt384;
  Real copt426 = copt20 + copt425;
  Real copt433 = -copt398;
  Real copt434 = copt32 + copt433;
  Real copt303 = -(copt2 * copt20 * copt29);
  Real copt304 = copt17 * copt2 * copt32;
  Real copt451 = copt4 * copt426;
  Real copt499 = t2(0);
  Real copt500 = Power(copt499, 2);
  Real copt502 = copt2 + copt311;
  Real copt503 = Power(copt502, 2);
  Real copt504 = copt13 + copt342;
  Real copt505 = Power(copt504, 2);
  Real copt506 = copt26 + copt327;
  Real copt507 = Power(copt506, 2);
  Real copt508 = copt503 + copt505 + copt507;
  Real copt509 = Sqrt(copt508);
  Real copt519 = xloc(15);
  Real copt732 = xloc(16);
  Real copt557 = -copt350;
  Real copt561 = -copt519;
  Real copt562 = copt561 + copt8;
  Real copt566 = copt4 * copt562;
  Real copt573 = copt519 * copt8;
  Real copt1044 = xloc(17);
  Real copt735 = -copt732;
  Real copt736 = copt20 + copt735;
  Real copt1045 = -copt1044;
  Real copt1046 = copt1045 + copt32;
  Real copt1080 = copt4 * copt736;
  Real copt1081 = copt732 * copt8;
  Real copt1105 = copt29 * copt562;
  Real copt1106 = copt32 * copt519;
  Real copt188 = copt182 * copt183;
  Real copt190 = copt182 * copt189;
  Real copt191 = -(copt183 * copt2 * copt8);
  Real copt192 = -(copt189 * copt2 * copt8);
  Real copt193 = -(copt17 * copt182 * copt20);
  Real copt194 = copt17 * copt2 * copt20 * copt4;
  Real copt195 = -(copt182 * copt29 * copt32);
  Real copt196 = copt2 * copt29 * copt32 * copt4;
  Real copt198 = -(copt183 * copt197 * copt2);
  Real copt199 = -(copt189 * copt197 * copt2);
  Real copt200 = copt183 * copt197 * copt8;
  Real copt201 = copt189 * copt197 * copt8;
  Real copt202 = copt17 * copt197 * copt2 * copt20;
  Real copt203 = -(copt17 * copt197 * copt20 * copt4);
  Real copt204 = copt197 * copt2 * copt29 * copt32;
  Real copt205 = -(copt197 * copt29 * copt32 * copt4);
  Real copt207 = -(copt17 * copt182 * copt206);
  Real copt208 = copt17 * copt2 * copt206 * copt4;
  Real copt209 = copt17 * copt2 * copt206 * copt8;
  Real copt210 = -(copt17 * copt206 * copt4 * copt8);
  Real copt211 = copt182 * copt20 * copt206;
  Real copt212 = -2 * copt2 * copt20 * copt206 * copt4;
  Real copt214 = copt20 * copt206 * copt213;
  Real copt215 = copt189 * copt20 * copt206;
  Real copt216 = -(copt17 * copt206 * copt29 * copt32);
  Real copt223 = copt20 * copt206;
  Real copt225 = -(copt17 * copt224);
  Real copt226 = copt183 + copt213 + copt219 + copt222 + copt223 + copt225;
  Real copt227 = copt217 * copt226;
  Real copt229 = -(copt182 * copt228 * copt29);
  Real copt230 = copt2 * copt228 * copt29 * copt4;
  Real copt231 = copt2 * copt228 * copt29 * copt8;
  Real copt232 = -(copt228 * copt29 * copt4 * copt8);
  Real copt233 = -(copt17 * copt20 * copt228 * copt29);
  Real copt234 = copt182 * copt228 * copt32;
  Real copt235 = -2 * copt2 * copt228 * copt32 * copt4;
  Real copt236 = copt213 * copt228 * copt32;
  Real copt237 = copt183 * copt228 * copt32;
  Real copt239 = copt228 * copt32;
  Real copt241 = -(copt240 * copt29);
  Real copt242 = copt189 + copt213 + copt219 + copt222 + copt239 + copt241;
  Real copt243 = copt238 * copt242;
  Real copt244 = copt17 * copt20 * copt29;
  Real copt245 = -(copt183 * copt32);
  Real copt246 = -2 * copt197 * copt29 * copt8;
  Real copt247 = copt17 * copt206 * copt29;
  Real copt248 = -2 * copt20 * copt206 * copt29;
  Real copt249 = copt17 * copt206 * copt32;
  Real copt250 = -(copt183 * copt228);
  Real copt251 = copt17 * copt20 * copt228;
  Real copt252 = -(copt213 * copt240);
  Real copt255 = copt228 * copt8;
  Real copt256 = copt253 + copt254 + copt255;
  Real copt257 = copt256 * copt4;
  Real copt258 = -(copt197 * copt32);
  Real copt262 = copt261 * copt4;
  Real copt263 = copt254 + copt258 + copt259 + copt262;
  Real copt264 = copt2 * copt263;
  Real copt265 = copt244 + copt245 + copt246 + copt247 + copt248 + copt249 +
                 copt250 + copt251 + copt252 + copt257 + copt264;
  Real copt266 = copt26 * copt265;
  Real copt267 = copt17 * copt4 * copt8;
  Real copt268 = -(copt20 * copt213);
  Real copt269 = -(copt189 * copt20);
  Real copt270 = copt17 * copt29 * copt32;
  Real copt271 = copt17 * copt197 * copt4;
  Real copt272 = -2 * copt17 * copt197 * copt8;
  Real copt273 = copt197 * copt20 * copt4;
  Real copt274 = -(copt206 * copt213);
  Real copt275 = -(copt189 * copt206);
  Real copt276 = copt206 * copt4 * copt8;
  Real copt277 = copt206 * copt29 * copt32;
  Real copt279 = copt17 * copt220;
  Real copt280 = -(copt206 * copt8);
  Real copt281 = -2 * copt17;
  Real copt282 = copt20 + copt206 + copt281;
  Real copt283 = copt282 * copt4;
  Real copt284 = copt278 + copt279 + copt280 + copt283;
  Real copt285 = copt2 * copt284;
  Real copt286 = copt17 * copt228 * copt29;
  Real copt287 = copt20 * copt228 * copt29;
  Real copt288 = -2 * copt17 * copt228 * copt32;
  Real copt289 = -(copt206 * copt32);
  Real copt290 = copt224 * copt29;
  Real copt291 = -(copt20 * copt228);
  Real copt292 = copt17 * copt261;
  Real copt293 = copt289 + copt290 + copt291 + copt292;
  Real copt294 = copt26 * copt293;
  Real copt295 = copt267 + copt268 + copt269 + copt270 + copt271 + copt272 +
                 copt273 + copt274 + copt275 + copt276 + copt277 + copt285 +
                 copt286 + copt287 + copt288 + copt294;
  Real copt296 = copt13 * copt295;
  Real copt297 = copt188 + copt190 + copt191 + copt192 + copt193 + copt194 +
                 copt195 + copt196 + copt198 + copt199 + copt200 + copt201 +
                 copt202 + copt203 + copt204 + copt205 + copt207 + copt208 +
                 copt209 + copt210 + copt211 + copt212 + copt214 + copt215 +
                 copt216 + copt227 + copt229 + copt230 + copt231 + copt232 +
                 copt233 + copt234 + copt235 + copt236 + copt237 + copt243 +
                 copt266 + copt296;
  Real copt298 = -(copt181 * copt297);
  Real copt299 = -2 * copt2 * copt4;
  Real copt300 = -2 * copt13 * copt17;
  Real copt301 = -2 * copt26 * copt29;
  Real copt302 = copt182 + copt183 + copt189 + copt213 + copt217 + copt238 +
                 copt299 + copt300 + copt301;
  Real copt305 = copt197 * copt20 * copt29;
  Real copt306 = -(copt17 * copt197 * copt32);
  Real copt307 = copt2 * copt206 * copt29;
  Real copt308 = -(copt206 * copt29 * copt8);
  Real copt309 = -(copt2 * copt206 * copt32);
  Real copt310 = copt206 * copt32 * copt4;
  Real copt312 = copt197 + copt311;
  Real copt313 = copt17 * copt312;
  Real copt314 = -copt206;
  Real copt315 = copt20 + copt314;
  Real copt316 = copt315 * copt4;
  Real copt317 = copt206 * copt8;
  Real copt318 = copt278 + copt313 + copt316 + copt317;
  Real copt319 = copt26 * copt318;
  Real copt320 = -(copt17 * copt2 * copt228);
  Real copt321 = copt17 * copt228 * copt8;
  Real copt322 = copt2 * copt20 * copt228;
  Real copt323 = -(copt20 * copt228 * copt4);
  Real copt324 = -copt197;
  Real copt325 = copt324 + copt8;
  Real copt326 = copt29 * copt325;
  Real copt328 = copt228 + copt327;
  Real copt329 = copt328 * copt4;
  Real copt330 = copt253 + copt259 + copt326 + copt329;
  Real copt331 = copt13 * copt330;
  Real copt332 = copt303 + copt304 + copt305 + copt306 + copt307 + copt308 +
                 copt309 + copt310 + copt319 + copt320 + copt321 + copt322 +
                 copt323 + copt331;
  Real copt333 = copt302 * copt332;
  Real copt334 = ArcTan(copt298, copt333);
  Real copt1145 = t0(1);
  Real copt349 = -(copt26 * copt29 * copt4 * copt8);
  Real copt351 = -(copt183 * copt350);
  Real copt352 = copt26 * copt29 * copt350;
  Real copt353 = -(copt189 * copt350);
  Real copt354 = -(copt17 * copt20 * copt26 * copt29);
  Real copt355 = 2 * copt17 * copt20 * copt4 * copt8;
  Real copt357 = -(copt213 * copt356);
  Real copt358 = copt26 * copt29 * copt356;
  Real copt359 = -(copt189 * copt356);
  Real copt360 = copt213 * copt26 * copt32;
  Real copt361 = copt183 * copt26 * copt32;
  Real copt362 = -(copt26 * copt32 * copt4 * copt8);
  Real copt363 = 2 * copt29 * copt32 * copt4 * copt8;
  Real copt364 = -(copt17 * copt20 * copt26 * copt32);
  Real copt365 = 2 * copt17 * copt20 * copt29 * copt32;
  Real copt367 = -(copt213 * copt366);
  Real copt368 = -(copt183 * copt366);
  Real copt370 = copt26 * copt29 * copt369 * copt4;
  Real copt371 = copt183 * copt369 * copt8;
  Real copt372 = -(copt26 * copt29 * copt369 * copt8);
  Real copt373 = copt189 * copt369 * copt8;
  Real copt374 = -(copt17 * copt20 * copt369 * copt4);
  Real copt375 = -(copt17 * copt20 * copt369 * copt8);
  Real copt376 = copt356 * copt369 * copt4;
  Real copt378 = -(copt26 * copt32 * copt369 * copt4);
  Real copt379 = -(copt29 * copt32 * copt369 * copt4);
  Real copt381 = copt26 * copt32 * copt369 * copt8;
  Real copt382 = -(copt29 * copt32 * copt369 * copt8);
  Real copt383 = copt366 * copt369 * copt4;
  Real copt385 = copt17 * copt26 * copt29 * copt384;
  Real copt386 = -(copt17 * copt384 * copt4 * copt8);
  Real copt387 = copt17 * copt350 * copt384;
  Real copt388 = copt20 * copt213 * copt384;
  Real copt389 = -(copt20 * copt26 * copt29 * copt384);
  Real copt390 = copt189 * copt20 * copt384;
  Real copt391 = -(copt20 * copt384 * copt4 * copt8);
  Real copt392 = -(copt17 * copt26 * copt32 * copt384);
  Real copt393 = -(copt17 * copt29 * copt32 * copt384);
  Real copt394 = copt20 * copt26 * copt32 * copt384;
  Real copt396 = -(copt20 * copt29 * copt32 * copt384);
  Real copt397 = copt17 * copt366 * copt384;
  Real copt399 = -(copt213 * copt26 * copt398);
  Real copt400 = -(copt183 * copt26 * copt398);
  Real copt401 = 2 * copt26 * copt398 * copt4 * copt8;
  Real copt402 = -(copt29 * copt398 * copt4 * copt8);
  Real copt403 = -(copt26 * copt350 * copt398);
  Real copt404 = copt29 * copt350 * copt398;
  Real copt405 = 2 * copt17 * copt20 * copt26 * copt398;
  Real copt406 = -(copt17 * copt20 * copt29 * copt398);
  Real copt407 = -(copt26 * copt356 * copt398);
  Real copt408 = copt29 * copt356 * copt398;
  Real copt409 = copt213 * copt32 * copt398;
  Real copt410 = copt183 * copt32 * copt398;
  Real copt411 = -(copt32 * copt398 * copt4 * copt8);
  Real copt412 = -(copt17 * copt20 * copt32 * copt398);
  Real copt413 = copt189 * copt20;
  Real copt414 = -(copt20 * copt29 * copt32);
  Real copt415 = copt20 * copt369 * copt8;
  Real copt418 = copt17 * copt417;
  Real copt419 = copt20 * copt369;
  Real copt420 = -2 * copt384;
  Real copt421 = copt20 + copt420;
  Real copt422 = copt421 * copt8;
  Real copt423 = copt418 + copt419 + copt422;
  Real copt424 = -(copt4 * copt423);
  Real copt427 = copt213 * copt426;
  Real copt428 = -(copt189 * copt384);
  Real copt429 = -(copt350 * copt384);
  Real copt430 = 2 * copt29 * copt32 * copt384;
  Real copt431 = -(copt366 * copt384);
  Real copt432 = -(copt369 * copt8);
  Real copt435 = -(copt345 * copt434);
  Real copt436 = copt350 + copt432 + copt435;
  Real copt437 = copt17 * copt436;
  Real copt438 = -(copt20 * copt29 * copt398);
  Real copt439 = copt20 * copt32 * copt398;
  Real copt440 = copt413 + copt414 + copt415 + copt424 + copt427 + copt428 +
                 copt429 + copt430 + copt431 + copt437 + copt438 + copt439;
  Real copt441 = copt13 * copt440;
  Real copt442 = copt356 * copt4;
  Real copt443 = copt366 * copt4;
  Real copt444 = copt183 * copt417;
  Real copt445 = copt189 * copt417;
  Real copt446 = -(copt356 * copt369);
  Real copt447 = -(copt366 * copt369);
  Real copt448 = -(copt20 * copt384 * copt4);
  Real copt449 = copt20 * copt384 * copt8;
  Real copt450 = -2 * copt20 * copt369;
  Real copt452 = copt20 + copt384;
  Real copt453 = copt452 * copt8;
  Real copt454 = copt450 + copt451 + copt453;
  Real copt455 = -(copt17 * copt454);
  Real copt456 = -(copt32 * copt398 * copt4);
  Real copt457 = copt32 * copt398 * copt8;
  Real copt458 = -2 * copt32 * copt369;
  Real copt459 = copt4 * copt434;
  Real copt460 = copt32 + copt398;
  Real copt461 = copt460 * copt8;
  Real copt462 = copt458 + copt459 + copt461;
  Real copt463 = -(copt29 * copt462);
  Real copt464 = copt442 + copt443 + copt444 + copt445 + copt446 + copt447 +
                 copt448 + copt449 + copt455 + copt456 + copt457 + copt463;
  Real copt465 = copt2 * copt464;
  Real copt466 = copt349 + copt351 + copt352 + copt353 + copt354 + copt355 +
                 copt357 + copt358 + copt359 + copt360 + copt361 + copt362 +
                 copt363 + copt364 + copt365 + copt367 + copt368 + copt370 +
                 copt371 + copt372 + copt373 + copt374 + copt375 + copt376 +
                 copt378 + copt379 + copt381 + copt382 + copt383 + copt385 +
                 copt386 + copt387 + copt388 + copt389 + copt390 + copt391 +
                 copt392 + copt393 + copt394 + copt396 + copt397 + copt399 +
                 copt400 + copt401 + copt402 + copt403 + copt404 + copt405 +
                 copt406 + copt407 + copt408 + copt409 + copt410 + copt411 +
                 copt412 + copt441 + copt465;
  Real copt467 = copt348 * copt466;
  Real copt468 = -2 * copt4 * copt8;
  Real copt469 = -2 * copt17 * copt20;
  Real copt470 = -2 * copt29 * copt32;
  Real copt471 = copt183 + copt189 + copt213 + copt350 + copt356 + copt366 +
                 copt468 + copt469 + copt470;
  Real copt472 = copt20 * copt29 * copt369;
  Real copt473 = -(copt17 * copt32 * copt369);
  Real copt474 = copt2 * copt29 * copt384;
  Real copt475 = -(copt29 * copt384 * copt8);
  Real copt476 = -(copt2 * copt32 * copt384);
  Real copt477 = copt32 * copt384 * copt4;
  Real copt478 = -(copt20 * copt369);
  Real copt479 = copt311 + copt369;
  Real copt480 = copt17 * copt479;
  Real copt481 = copt384 * copt8;
  Real copt482 = copt451 + copt478 + copt480 + copt481;
  Real copt483 = copt26 * copt482;
  Real copt484 = -(copt17 * copt2 * copt398);
  Real copt485 = copt17 * copt398 * copt8;
  Real copt486 = copt2 * copt20 * copt398;
  Real copt487 = -(copt20 * copt398 * copt4);
  Real copt488 = copt29 * copt417;
  Real copt489 = copt32 * copt369;
  Real copt490 = -(copt398 * copt8);
  Real copt491 = copt327 + copt398;
  Real copt492 = copt4 * copt491;
  Real copt493 = copt488 + copt489 + copt490 + copt492;
  Real copt494 = copt13 * copt493;
  Real copt495 = copt303 + copt304 + copt472 + copt473 + copt474 + copt475 +
                 copt476 + copt477 + copt483 + copt484 + copt485 + copt486 +
                 copt487 + copt494;
  Real copt496 = copt471 * copt495;
  Real copt497 = ArcTan(copt467, copt496);
  Real copt1148 = t1(1);
  Real copt510 = copt17 * copt182 * copt20;
  Real copt512 = -(copt17 * copt2 * copt20 * copt8);
  Real copt513 = -(copt182 * copt356);
  Real copt514 = copt2 * copt356 * copt4;
  Real copt515 = copt182 * copt29 * copt32;
  Real copt516 = -(copt2 * copt29 * copt32 * copt8);
  Real copt517 = -(copt182 * copt366);
  Real copt518 = copt2 * copt366 * copt4;
  Real copt523 = -(copt17 * copt2 * copt20 * copt519);
  Real copt526 = copt17 * copt20 * copt519 * copt8;
  Real copt529 = copt2 * copt356 * copt519;
  Real copt534 = -(copt356 * copt4 * copt519);
  Real copt537 = -(copt2 * copt29 * copt32 * copt519);
  Real copt540 = copt29 * copt32 * copt519 * copt8;
  Real copt547 = copt2 * copt366 * copt519;
  Real copt553 = -(copt366 * copt4 * copt519);
  Real copt791 = copt343 * copt736;
  Real copt891 = copt557 + copt566 + copt573 + copt791;
  Real copt894 = copt217 * copt891;
  Real copt914 = -(copt17 * copt182 * copt732);
  Real copt968 = 2 * copt17 * copt2 * copt732 * copt8;
  Real copt1034 = -(copt17 * copt350 * copt732);
  Real copt1035 = copt182 * copt20 * copt732;
  Real copt1036 = -(copt2 * copt20 * copt4 * copt732);
  Real copt1038 = -(copt2 * copt20 * copt732 * copt8);
  Real copt1039 = copt20 * copt4 * copt732 * copt8;
  Real copt1041 = copt20 * copt29 * copt32 * copt732;
  Real copt1042 = -(copt17 * copt366 * copt732);
  Real copt1047 = copt1046 * copt345;
  Real copt1049 = copt1047 + copt557 + copt566 + copt573;
  Real copt1050 = copt1049 * copt238;
  Real copt1052 = -(copt1044 * copt182 * copt29);
  Real copt1053 = 2 * copt1044 * copt2 * copt29 * copt8;
  Real copt1055 = -(copt1044 * copt29 * copt350);
  Real copt1056 = -(copt1044 * copt29 * copt356);
  Real copt1057 = copt1044 * copt182 * copt32;
  Real copt1058 = -(copt1044 * copt2 * copt32 * copt4);
  Real copt1059 = -(copt1044 * copt2 * copt32 * copt8);
  Real copt1060 = copt1044 * copt32 * copt4 * copt8;
  Real copt1061 = copt1044 * copt17 * copt20 * copt32;
  Real copt1062 = -(copt17 * copt350);
  Real copt1063 = copt20 * copt26 * copt29;
  Real copt1064 = copt20 * copt4 * copt8;
  Real copt1065 = -2 * copt20 * copt26 * copt32;
  Real copt1066 = copt20 * copt29 * copt32;
  Real copt1067 = copt17 * copt519 * copt8;
  Real copt1069 = -2 * copt20 * copt4 * copt519;
  Real copt1070 = copt20 * copt519 * copt8;
  Real copt1071 = -(copt26 * copt29 * copt732);
  Real copt1072 = copt4 * copt732 * copt8;
  Real copt1073 = -(copt350 * copt732);
  Real copt1074 = copt26 * copt32 * copt732;
  Real copt1075 = copt29 * copt32 * copt732;
  Real copt1076 = -(copt366 * copt732);
  Real copt1077 = -2 * copt20 * copt8;
  Real copt1078 = copt17 * copt562;
  Real copt1079 = copt20 * copt519;
  Real copt1082 = copt1077 + copt1078 + copt1079 + copt1080 + copt1081;
  Real copt1083 = copt1082 * copt2;
  Real copt1084 = copt1046 * copt17 * copt506;
  Real copt1085 = copt1044 * copt20 * copt26;
  Real copt1086 = -2 * copt1044 * copt20 * copt29;
  Real copt1087 = copt1044 * copt20 * copt32;
  Real copt1088 = copt1062 + copt1063 + copt1064 + copt1065 + copt1066 +
                  copt1067 + copt1069 + copt1070 + copt1071 + copt1072 +
                  copt1073 + copt1074 + copt1075 + copt1076 + copt1083 +
                  copt1084 + copt1085 + copt1086 + copt1087;
  Real copt1089 = -(copt1088 * copt13);
  Real copt1090 = copt32 * copt4 * copt8;
  Real copt1091 = copt17 * copt20 * copt32;
  Real copt1092 = -2 * copt32 * copt4 * copt519;
  Real copt1093 = copt32 * copt519 * copt8;
  Real copt1094 = -2 * copt17 * copt32 * copt732;
  Real copt1095 = copt20 * copt32 * copt732;
  Real copt1096 = -(copt519 * copt8);
  Real copt1097 = -(copt20 * copt732);
  Real copt1098 = copt1096 + copt1097 + copt350 + copt356;
  Real copt1099 = -(copt1098 * copt29);
  Real copt1100 = copt1044 * copt4 * copt8;
  Real copt1101 = -(copt1044 * copt350);
  Real copt1102 = copt1044 * copt17 * copt20;
  Real copt1103 = -(copt1044 * copt356);
  Real copt1104 = -2 * copt32 * copt8;
  Real copt1107 = copt1046 * copt4;
  Real copt1108 = copt1044 * copt8;
  Real copt1109 = copt1104 + copt1105 + copt1106 + copt1107 + copt1108;
  Real copt1110 = copt1109 * copt2;
  Real copt1111 = copt1090 + copt1091 + copt1092 + copt1093 + copt1094 +
                  copt1095 + copt1099 + copt1100 + copt1101 + copt1102 +
                  copt1103 + copt1110;
  Real copt1112 = -(copt1111 * copt26);
  Real copt1113 = copt1034 + copt1035 + copt1036 + copt1038 + copt1039 +
                  copt1041 + copt1042 + copt1050 + copt1052 + copt1053 +
                  copt1055 + copt1056 + copt1057 + copt1058 + copt1059 +
                  copt1060 + copt1061 + copt1089 + copt1112 + copt510 +
                  copt512 + copt513 + copt514 + copt515 + copt516 + copt517 +
                  copt518 + copt523 + copt526 + copt529 + copt534 + copt537 +
                  copt540 + copt547 + copt553 + copt894 + copt914 + copt968;
  Real copt1114 = copt1113 * copt509;
  Real copt1115 = -2 * copt2 * copt8;
  Real copt1116 = -2 * copt13 * copt20;
  Real copt1117 = -2 * copt26 * copt32;
  Real copt1118 = copt1115 + copt1116 + copt1117 + copt182 + copt217 + copt238 +
                  copt350 + copt356 + copt366;
  Real copt1119 = copt20 * copt29 * copt519;
  Real copt1120 = -(copt17 * copt32 * copt519);
  Real copt1121 = copt2 * copt29 * copt732;
  Real copt1122 = -(copt29 * copt732 * copt8);
  Real copt1123 = -(copt2 * copt32 * copt732);
  Real copt1124 = copt32 * copt4 * copt732;
  Real copt1125 = -(copt20 * copt519);
  Real copt1126 = copt311 + copt519;
  Real copt1127 = copt1126 * copt17;
  Real copt1128 = copt1080 + copt1081 + copt1125 + copt1127;
  Real copt1129 = copt1128 * copt26;
  Real copt1130 = -(copt1044 * copt17 * copt2);
  Real copt1131 = copt1044 * copt17 * copt8;
  Real copt1132 = copt1044 * copt2 * copt20;
  Real copt1133 = -(copt1044 * copt20 * copt4);
  Real copt1134 = -(copt1044 * copt8);
  Real copt1135 = copt1044 + copt327;
  Real copt1136 = copt1135 * copt4;
  Real copt1137 = copt1105 + copt1106 + copt1134 + copt1136;
  Real copt1138 = copt1137 * copt13;
  Real copt1139 = copt1119 + copt1120 + copt1121 + copt1122 + copt1123 +
                  copt1124 + copt1129 + copt1130 + copt1131 + copt1132 +
                  copt1133 + copt1138 + copt303 + copt304;
  Real copt1140 = copt1118 * copt1139;
  Real copt1141 = ArcTan(copt1114, copt1140);
  Real copt1151 = t2(1);
  Real copt1156 = Power(copt1145, 2);
  Real copt1159 = Power(copt1148, 2);
  Real copt1162 = Power(copt1151, 2);
  Real copt1167 = copt1 + copt7;
  Real copt1168 = copt1 * copt138;
  Real copt1169 = copt502 * copt7;
  Real copt1170 = copt1168 + copt1169;
  Real copt1186 = copt42 + copt45;
  Real copt1172 = copt1 * copt173;
  Real copt1173 = copt504 * copt7;
  Real copt1174 = copt1172 + copt1173;
  Real copt1176 = copt1 * copt177;
  Real copt1177 = copt506 * copt7;
  Real copt1178 = copt1176 + copt1177;
  Real copt1188 = copt138 * copt42;
  Real copt1189 = copt45 * copt502;
  Real copt1190 = copt1188 + copt1189;
  Real copt1194 = copt173 * copt42;
  Real copt1195 = copt45 * copt504;
  Real copt1196 = copt1194 + copt1195;
  Real copt1203 = copt177 * copt42;
  Real copt1204 = copt45 * copt506;
  Real copt1205 = copt1203 + copt1204;
  Real copt1259 = Power(copt302, 2);
  Real copt1263 = Power(copt332, 2);
  Real copt1264 = copt1259 * copt1263;
  Real copt1265 = Power(copt297, 2);
  Real copt1266 = copt1265 * copt179;
  Real copt1268 = copt1264 + copt1266;
  Real copt1269 = 1 / copt1268;
  Real copt1271 = 1 / copt181;
  Real copt1306 = Power(copt471, 2);
  Real copt1307 = Power(copt495, 2);
  Real copt1308 = copt1306 * copt1307;
  Real copt1309 = Power(copt466, 2);
  Real copt1310 = copt1309 * copt347;
  Real copt1311 = copt1308 + copt1310;
  Real copt1312 = 1 / copt1311;
  Real copt1317 = -(copt20 * copt29);
  Real copt1318 = copt17 * copt32;
  Real copt1342 = Power(copt1118, 2);
  Real copt1343 = Power(copt1139, 2);
  Real copt1344 = copt1342 * copt1343;
  Real copt1345 = Power(copt1113, 2);
  Real copt1346 = copt1345 * copt508;
  Real copt1347 = copt1344 + copt1346;
  Real copt1348 = 1 / copt1347;
  Real copt1376 = 1 / copt509;
  Real copt1385 = 2 * copt13;
  Real copt1442 = 2 * copt26;
  Real copt1465 = -(copt32 * copt4 * copt8);
  Real copt1466 = -(copt17 * copt20 * copt32);
  Real copt1333 = copt1044 * copt20;
  Real copt1530 = 2 * copt4;
  Real copt1537 = copt1530 + copt311 + copt324;
  Real copt1338 = -2 * copt8;
  Real copt1618 = 1 / copt348;
  Real copt1508 = copt32 * copt732;
  Real copt1666 = 2 * copt17;
  Real copt1504 = copt20 * copt29;
  Real copt1522 = copt206 * copt32;
  Real copt1254 = copt20 * copt228;
  Real copt1660 = copt2 * copt32;
  Real copt1406 = -2 * copt20;
  Real copt1412 = -(copt20 * copt26 * copt29);
  Real copt1549 = copt197 * copt20;
  Real copt1700 = -2 * copt197 * copt8;
  Real copt1701 = -2 * copt4;
  Real copt1702 = copt1701 + copt197 + copt8;
  Real copt1703 = copt1702 * copt2;
  Real copt1799 = 2 * copt29;
  Real copt1248 = -copt228;
  Real copt1793 = -(copt2 * copt20);
  Real copt1487 = -2 * copt32;
  Real copt1469 = -(copt32 * copt369 * copt4);
  Real copt1473 = -(copt17 * copt32 * copt384);
  Real copt1572 = -(copt20 * copt398);
  Real copt1770 = copt2 * copt562;
  Real copt1774 = copt20 * copt32;
  Real copt1929 = copt134 + copt197;
  Real copt1542 = copt17 * copt2 * copt206;
  Real copt1555 = copt2 * copt228 * copt29;
  Real copt1942 = copt1248 + copt29;
  Real copt1834 = copt17 * copt228;
  Real copt1592 = -(copt17 * copt20 * copt369);
  Real copt1595 = -(copt29 * copt32 * copt369);
  Real copt1970 = 2 * copt8;
  Real copt1507 = -(copt29 * copt732);
  Real copt1529 = -2 * copt2;
  Real copt2030 = copt1338 + copt4 + copt519;
  Real copt1632 = -(copt2 * copt20 * copt732);
  Real copt1639 = -(copt1044 * copt2 * copt32);
  Real copt1952 = copt140 + copt206;
  Real copt1832 = -2 * copt206 * copt29;
  Real copt1699 = copt197 * copt4;
  Real copt1705 = copt228 * copt29;
  Real copt2088 = copt324 + copt4;
  Real copt2096 = -(copt2 * copt29);
  Real copt1738 = -(copt384 * copt4 * copt8);
  Real copt1741 = -(copt29 * copt32 * copt384);
  Real copt2002 = copt32 * copt398;
  Real copt1665 = -2 * copt13;
  Real copt2112 = 2 * copt20;
  Real copt2041 = copt1406 + copt17 + copt732;
  Real copt2018 = copt1044 * copt17;
  Real copt1907 = -2 * copt1044 * copt20;
  Real copt1694 = copt4 * copt8;
  Real copt1698 = copt29 * copt32;
  Real copt2086 = -copt213;
  Real copt2089 = copt2 * copt2088;
  Real copt1817 = copt17 * copt206;
  Real copt2083 = copt17 * copt29;
  Real copt1935 = copt17 + copt314;
  Real copt1690 = copt206 * copt29;
  Real copt1691 = -2 * copt17 * copt228;
  Real copt1937 = copt206 * copt4;
  Real copt2230 = copt17 * copt2;
  Real copt1853 = -(copt26 * copt4 * copt8);
  Real copt1856 = -(copt17 * copt20 * copt26);
  Real copt1468 = -(copt29 * copt369 * copt8);
  Real copt1472 = -(copt20 * copt29 * copt384);
  Real copt1871 = -(copt398 * copt4 * copt8);
  Real copt1873 = -(copt17 * copt20 * copt398);
  Real copt2134 = -2 * copt369;
  Real copt2136 = copt2134 + copt4 + copt8;
  Real copt2109 = -(copt398 * copt4);
  Real copt1720 = copt398 * copt8;
  Real copt2357 = 2 * copt32;
  Real copt1322 = copt20 * copt398;
  Real copt2039 = copt4 * copt732;
  Real copt1798 = -2 * copt26;
  Real copt1815 = copt17 * copt20;
  Real copt2191 = -2 * copt4 * copt519;
  Real copt2192 = copt2 * copt2030;
  Real copt1899 = copt20 * copt732;
  Real copt2056 = copt1044 + copt1487 + copt29;
  Real copt1330 = copt29 * copt732;
  Real copt1775 = -2 * copt32 * copt732;
  Real copt1925 = -(copt183 * copt2);
  Real copt1926 = -(copt189 * copt2);
  Real copt3656 = copt134 + copt8;
  Real copt1535 = copt17 * copt2 * copt20;
  Real copt1933 = copt17 * copt4;
  Real copt2036 = -2 * copt17 * copt8;
  Real copt2037 = copt20 * copt4;
  Real copt1536 = copt2 * copt29 * copt32;
  Real copt2072 = -(copt17 * copt182);
  Real copt2073 = copt17 * copt2 * copt4;
  Real copt1764 = copt182 * copt20;
  Real copt1828 = -2 * copt20 * copt29;
  Real copt2087 = -copt189;
  Real copt2048 = copt32 * copt4;
  Real copt2208 = -(copt182 * copt29);
  Real copt2209 = copt2 * copt29 * copt4;
  Real copt1459 = -(copt29 * copt4 * copt8);
  Real copt1461 = -(copt17 * copt20 * copt29);
  Real copt2213 = -copt183;
  Real copt3894 = copt2 * copt340;
  Real copt1895 = copt182 * copt32;
  Real copt1463 = copt213 * copt32;
  Real copt1464 = copt183 * copt32;
  Real copt3826 = copt175 + copt32;
  Real copt3815 = copt26 * copt343;
  Real copt1689 = -2 * copt17 * copt32;
  Real copt1546 = copt17 * copt8;
  Real copt3637 = copt183 * copt8;
  Real copt1584 = -(copt26 * copt29 * copt8);
  Real copt3646 = copt189 * copt8;
  Real copt3702 = -(copt17 * copt20 * copt4);
  Real copt1351 = -(copt17 * copt20 * copt8);
  Real copt1634 = copt20 * copt8;
  Real copt3875 = copt140 + copt20;
  Real copt1980 = -(copt26 * copt32 * copt4);
  Real copt3761 = -(copt29 * copt32 * copt4);
  Real copt1354 = -(copt29 * copt32 * copt8);
  Real copt1898 = -copt356;
  Real copt3825 = -(copt17 * copt32);
  Real copt3827 = copt13 * copt3826;
  Real copt3835 = copt1504 + copt3815 + copt3825 + copt3827;
  Real copt3855 = -(copt17 * copt4 * copt8);
  Real copt1411 = copt17 * copt350;
  Real copt3874 = copt20 * copt213;
  Real copt1413 = -(copt20 * copt4 * copt8);
  Real copt4041 = -(copt20 * copt4);
  Real copt2123 = -(copt17 * copt26 * copt32);
  Real copt3879 = -(copt17 * copt29 * copt32);
  Real copt4102 = 2 * copt29 * copt32;
  Real copt4103 = -copt366;
  Real copt3922 = copt2 * copt29;
  Real copt3936 = -(copt29 * copt8);
  Real copt3937 = copt26 * copt3656;
  Real copt3938 = -(copt2 * copt32);
  Real copt3943 = copt2048 + copt3922 + copt3936 + copt3937 + copt3938;
  Real copt1460 = copt29 * copt350;
  Real copt1462 = copt29 * copt356;
  Real copt1641 = copt32 * copt8;
  Real copt4019 = -(copt17 * copt2);
  Real copt4020 = copt13 * copt340;
  Real copt4026 = copt2 * copt20;
  Real copt4045 = copt1546 + copt4019 + copt4020 + copt4026 + copt4041;
  Real copt3675 = copt238 * copt3656;
  Real copt3677 = copt217 * copt3656;
  Real copt2025 = -(copt17 * copt2 * copt20);
  Real copt1625 = copt2 * copt356;
  Real copt1547 = -2 * copt20 * copt4;
  Real copt2026 = -(copt2 * copt29 * copt32);
  Real copt1626 = copt2 * copt366;
  Real copt1765 = -(copt2 * copt20 * copt8);
  Real copt3876 = copt217 * copt3875;
  Real copt4572 = copt2 * copt3656;
  Real copt1896 = -(copt2 * copt32 * copt8);
  Real copt3988 = copt238 * copt3826;
  Real copt1905 = copt20 * copt26;
  Real copt1246 = copt206 + copt342;
  Real copt1247 = copt1246 * copt29;
  Real copt1249 = copt1248 + copt32;
  Real copt1253 = copt1249 * copt17;
  Real copt1255 = copt1247 + copt1253 + copt1254 + copt289;
  Real copt1256 = copt1255 * copt302;
  Real copt1257 = 2 * copt138 * copt332;
  Real copt1258 = copt1256 + copt1257;
  Real copt1270 = -(copt1258 * copt1269 * copt181 * copt297);
  Real copt1272 = -2 * copt26 * copt29 * copt4;
  Real copt1273 = -(copt183 * copt8);
  Real copt1274 = copt26 * copt29 * copt8;
  Real copt1275 = -(copt189 * copt8);
  Real copt1276 = copt17 * copt20 * copt4;
  Real copt1278 = copt26 * copt32 * copt4;
  Real copt1279 = copt29 * copt32 * copt4;
  Real copt1280 = -(copt183 * copt197);
  Real copt1282 = copt197 * copt26 * copt29;
  Real copt1283 = -(copt189 * copt197);
  Real copt1284 = copt17 * copt197 * copt20;
  Real copt1285 = -(copt197 * copt26 * copt32);
  Real copt1286 = copt197 * copt29 * copt32;
  Real copt1287 = copt17 * copt206 * copt4;
  Real copt1288 = copt17 * copt206 * copt8;
  Real copt1289 = -2 * copt20 * copt206 * copt4;
  Real copt1290 = copt13 * copt284;
  Real copt1291 = copt228 * copt26 * copt4;
  Real copt1292 = copt228 * copt29 * copt4;
  Real copt1293 = -(copt228 * copt26 * copt8);
  Real copt1295 = copt228 * copt29 * copt8;
  Real copt1296 = -2 * copt228 * copt32 * copt4;
  Real copt1297 = copt183 + copt189 + copt223 + copt225 + copt239 + copt241;
  Real copt1298 = 2 * copt1297 * copt2;
  Real copt1299 = copt1272 + copt1273 + copt1274 + copt1275 + copt1276 +
                  copt1278 + copt1279 + copt1280 + copt1282 + copt1283 +
                  copt1284 + copt1285 + copt1286 + copt1287 + copt1288 +
                  copt1289 + copt1290 + copt1291 + copt1292 + copt1293 +
                  copt1295 + copt1296 + copt1298;
  Real copt1300 = -(copt1299 * copt179);
  Real copt1301 = -(copt138 * copt297);
  Real copt1302 = copt1300 + copt1301;
  Real copt1303 = -(copt1269 * copt1271 * copt1302 * copt302 * copt332);
  Real copt1304 = copt1270 + copt1303;
  Real copt1313 = -(copt1312 * copt348 * copt464 * copt471 * copt495);
  Real copt1319 = copt29 * copt384;
  Real copt1320 = -(copt32 * copt384);
  Real copt1321 = -(copt17 * copt398);
  Real copt1323 =
      copt1317 + copt1318 + copt1319 + copt1320 + copt1321 + copt1322;
  Real copt1324 = copt1312 * copt1323 * copt348 * copt466 * copt471;
  Real copt1328 = copt1313 + copt1324;
  Real copt1331 = -(copt32 * copt732);
  Real copt1332 = -(copt1044 * copt17);
  Real copt1335 =
      copt1317 + copt1318 + copt1330 + copt1331 + copt1332 + copt1333;
  Real copt1336 = copt1118 * copt1335;
  Real copt1337 = 2 * copt2;
  Real copt1339 = copt1337 + copt1338;
  Real copt1340 = copt1139 * copt1339;
  Real copt1341 = copt1336 + copt1340;
  Real copt1349 = copt1113 * copt1341 * copt1348 * copt509;
  Real copt1350 = 2 * copt17 * copt2 * copt20;
  Real copt1352 = -2 * copt2 * copt356;
  Real copt1353 = 2 * copt2 * copt29 * copt32;
  Real copt1355 = -2 * copt2 * copt366;
  Real copt1356 = -(copt17 * copt20 * copt519);
  Real copt1357 = copt356 * copt519;
  Real copt1358 = -(copt29 * copt32 * copt519);
  Real copt1359 = copt366 * copt519;
  Real copt1360 = -2 * copt17 * copt2 * copt732;
  Real copt1362 = 2 * copt17 * copt732 * copt8;
  Real copt1363 = 2 * copt2 * copt20 * copt732;
  Real copt1364 = -(copt20 * copt4 * copt732);
  Real copt1365 = -(copt20 * copt732 * copt8);
  Real copt1366 = -(copt1082 * copt13);
  Real copt1367 = -2 * copt1044 * copt2 * copt29;
  Real copt1368 = 2 * copt1044 * copt29 * copt8;
  Real copt1369 = 2 * copt1044 * copt2 * copt32;
  Real copt1370 = -(copt1044 * copt32 * copt4);
  Real copt1371 = -(copt1044 * copt32 * copt8);
  Real copt1373 = -(copt1109 * copt26);
  Real copt1374 = copt1350 + copt1351 + copt1352 + copt1353 + copt1354 +
                  copt1355 + copt1356 + copt1357 + copt1358 + copt1359 +
                  copt1360 + copt1362 + copt1363 + copt1364 + copt1365 +
                  copt1366 + copt1367 + copt1368 + copt1369 + copt1370 +
                  copt1371 + copt1373 + copt442 + copt443;
  Real copt1375 = copt1374 * copt509;
  Real copt1377 = copt1113 * copt1376 * copt502;
  Real copt1378 = copt1375 + copt1377;
  Real copt1379 = -(copt1118 * copt1139 * copt1348 * copt1378);
  Real copt1380 = copt1349 + copt1379;
  Real copt1384 = copt302 * copt330;
  Real copt1386 = copt1385 + copt281;
  Real copt1387 = copt1386 * copt332;
  Real copt1388 = copt1384 + copt1387;
  Real copt1389 = -(copt1269 * copt1388 * copt181 * copt297);
  Real copt1390 = 2 * copt13 * copt242;
  Real copt1391 = copt1390 + copt267 + copt268 + copt269 + copt270 + copt271 +
                  copt272 + copt273 + copt274 + copt275 + copt276 + copt277 +
                  copt285 + copt286 + copt287 + copt288 + copt294;
  Real copt1392 = -(copt1391 * copt181);
  Real copt1396 = -(copt1271 * copt173 * copt297);
  Real copt1397 = copt1392 + copt1396;
  Real copt1398 = -(copt1269 * copt1397 * copt302 * copt332);
  Real copt1399 = copt1389 + copt1398;
  Real copt1401 = -(copt1312 * copt348 * copt440 * copt471 * copt495);
  Real copt1402 = copt1312 * copt348 * copt466 * copt471 * copt493;
  Real copt1403 = copt1401 + copt1402;
  Real copt1405 = copt1118 * copt1137;
  Real copt1407 = copt1385 + copt1406;
  Real copt1408 = copt1139 * copt1407;
  Real copt1409 = copt1405 + copt1408;
  Real copt1410 = copt1113 * copt1348 * copt1409 * copt509;
  Real copt1414 = 2 * copt20 * copt26 * copt32;
  Real copt1415 = -(copt17 * copt519 * copt8);
  Real copt1416 = 2 * copt20 * copt4 * copt519;
  Real copt1417 = -(copt20 * copt519 * copt8);
  Real copt1418 = copt26 * copt29 * copt732;
  Real copt1419 = -(copt4 * copt732 * copt8);
  Real copt1420 = copt350 * copt732;
  Real copt1421 = -(copt26 * copt32 * copt732);
  Real copt1423 = -(copt29 * copt32 * copt732);
  Real copt1424 = copt366 * copt732;
  Real copt1425 = -(copt1082 * copt2);
  Real copt1426 = 2 * copt1049 * copt13;
  Real copt1427 = -(copt1046 * copt17 * copt506);
  Real copt1428 = -(copt1044 * copt20 * copt26);
  Real copt1429 = 2 * copt1044 * copt20 * copt29;
  Real copt1430 = -(copt1044 * copt20 * copt32);
  Real copt1431 = copt1411 + copt1412 + copt1413 + copt1414 + copt1415 +
                  copt1416 + copt1417 + copt1418 + copt1419 + copt1420 +
                  copt1421 + copt1423 + copt1424 + copt1425 + copt1426 +
                  copt1427 + copt1428 + copt1429 + copt1430 + copt414;
  Real copt1432 = copt1431 * copt509;
  Real copt1433 = copt1113 * copt1376 * copt504;
  Real copt1434 = copt1432 + copt1433;
  Real copt1435 = -(copt1118 * copt1139 * copt1348 * copt1434);
  Real copt1436 = copt1410 + copt1435;
  Real copt1440 = copt302 * copt318;
  Real copt1443 = copt1442 + copt260;
  Real copt1444 = copt1443 * copt332;
  Real copt1445 = copt1440 + copt1444;
  Real copt1446 = -(copt1269 * copt1445 * copt181 * copt297);
  Real copt1447 = 2 * copt226 * copt26;
  Real copt1448 = copt13 * copt293;
  Real copt1449 = copt1447 + copt1448 + copt244 + copt245 + copt246 + copt247 +
                  copt248 + copt249 + copt250 + copt251 + copt252 + copt257 +
                  copt264;
  Real copt1450 = -(copt1449 * copt181);
  Real copt1454 = -(copt1271 * copt177 * copt297);
  Real copt1455 = copt1450 + copt1454;
  Real copt1456 = -(copt1269 * copt1455 * copt302 * copt332);
  Real copt1457 = copt1446 + copt1456;
  Real copt1467 = copt29 * copt369 * copt4;
  Real copt1470 = copt32 * copt369 * copt8;
  Real copt1471 = copt17 * copt29 * copt384;
  Real copt1474 = copt20 * copt32 * copt384;
  Real copt1475 = -(copt213 * copt398);
  Real copt1476 = -(copt183 * copt398);
  Real copt1477 = 2 * copt398 * copt4 * copt8;
  Real copt1478 = -(copt350 * copt398);
  Real copt1479 = 2 * copt17 * copt20 * copt398;
  Real copt1480 = -(copt356 * copt398);
  Real copt1481 = copt1459 + copt1460 + copt1461 + copt1462 + copt1463 +
                  copt1464 + copt1465 + copt1466 + copt1467 + copt1468 +
                  copt1469 + copt1470 + copt1471 + copt1472 + copt1473 +
                  copt1474 + copt1475 + copt1476 + copt1477 + copt1478 +
                  copt1479 + copt1480;
  Real copt1482 = -(copt1312 * copt1481 * copt348 * copt471 * copt495);
  Real copt1483 = copt1312 * copt348 * copt466 * copt471 * copt482;
  Real copt1484 = copt1482 + copt1483;
  Real copt1486 = copt1118 * copt1128;
  Real copt1488 = copt1442 + copt1487;
  Real copt1489 = copt1139 * copt1488;
  Real copt1490 = copt1486 + copt1489;
  Real copt1491 = copt1113 * copt1348 * copt1490 * copt509;
  Real copt1492 = 2 * copt32 * copt4 * copt519;
  Real copt1493 = -(copt32 * copt519 * copt8);
  Real copt1495 = 2 * copt26 * copt891;
  Real copt1496 = 2 * copt17 * copt32 * copt732;
  Real copt1497 = -(copt20 * copt32 * copt732);
  Real copt1498 = copt1098 * copt29;
  Real copt1499 = -(copt1044 * copt4 * copt8);
  Real copt1500 = copt1044 * copt350;
  Real copt1501 = -(copt1044 * copt17 * copt20);
  Real copt1502 = copt1044 * copt356;
  Real copt1503 = -(copt1109 * copt2);
  Real copt1506 = -2 * copt20 * copt32;
  Real copt1509 = copt1046 * copt17;
  Real copt1510 =
      copt1333 + copt1504 + copt1506 + copt1507 + copt1508 + copt1509;
  Real copt1511 = -(copt13 * copt1510);
  Real copt1512 = copt1465 + copt1466 + copt1492 + copt1493 + copt1495 +
                  copt1496 + copt1497 + copt1498 + copt1499 + copt1500 +
                  copt1501 + copt1502 + copt1503 + copt1511;
  Real copt1513 = copt1512 * copt509;
  Real copt1514 = copt1113 * copt1376 * copt506;
  Real copt1515 = copt1513 + copt1514;
  Real copt1516 = -(copt1118 * copt1139 * copt1348 * copt1515);
  Real copt1517 = copt1491 + copt1516;
  Real copt1521 = copt26 * copt315;
  Real copt1523 = copt13 * copt328;
  Real copt1527 = copt1521 + copt1522 + copt1523 + copt291;
  Real copt1528 = copt1527 * copt302;
  Real copt1531 = copt1529 + copt1530;
  Real copt1532 = copt1531 * copt332;
  Real copt1533 = copt1528 + copt1532;
  Real copt1534 = -(copt1269 * copt1533 * copt181 * copt297);
  Real copt1538 = copt1537 * copt238;
  Real copt1539 = copt1537 * copt217;
  Real copt1540 = -(copt17 * copt197 * copt20);
  Real copt1541 = -(copt197 * copt29 * copt32);
  Real copt1543 = -(copt17 * copt206 * copt8);
  Real copt1544 = -2 * copt2 * copt20 * copt206;
  Real copt1545 = 2 * copt20 * copt206 * copt4;
  Real copt1548 = copt17 * copt197;
  Real copt1551 = -2 * copt206 * copt4;
  Real copt1552 = copt2 * copt282;
  Real copt1553 =
      copt1546 + copt1547 + copt1548 + copt1549 + copt1551 + copt1552 + copt317;
  Real copt1554 = copt13 * copt1553;
  Real copt1556 = -(copt228 * copt29 * copt8);
  Real copt1557 = -2 * copt2 * copt228 * copt32;
  Real copt1558 = 2 * copt228 * copt32 * copt4;
  Real copt1559 = -2 * copt240 * copt4;
  Real copt1560 = copt2 * copt261;
  Real copt1561 = copt1559 + copt1560 + copt253 + copt254 + copt255;
  Real copt1562 = copt1561 * copt26;
  Real copt1563 = copt1535 + copt1536 + copt1538 + copt1539 + copt1540 +
                  copt1541 + copt1542 + copt1543 + copt1544 + copt1545 +
                  copt1554 + copt1555 + copt1556 + copt1557 + copt1558 +
                  copt1562;
  Real copt1564 = -(copt1563 * copt181);
  Real copt1565 = copt1271 * copt138 * copt297;
  Real copt1566 = copt1564 + copt1565;
  Real copt1567 = -(copt1269 * copt1566 * copt302 * copt332);
  Real copt1568 = copt1534 + copt1567;
  Real copt1570 = copt26 * copt426;
  Real copt1571 = copt32 * copt384;
  Real copt1576 = copt13 * copt491;
  Real copt1577 = copt1570 + copt1571 + copt1572 + copt1576;
  Real copt1578 = copt1577 * copt471;
  Real copt1579 = copt1338 + copt1530;
  Real copt1580 = copt1579 * copt495;
  Real copt1581 = copt1578 + copt1580;
  Real copt1583 = copt1312 * copt1581 * copt348 * copt466;
  Real copt1585 = 2 * copt17 * copt20 * copt8;
  Real copt1586 = -2 * copt356 * copt4;
  Real copt1587 = 2 * copt26 * copt32 * copt4;
  Real copt1588 = -(copt26 * copt32 * copt8);
  Real copt1589 = 2 * copt29 * copt32 * copt8;
  Real copt1590 = -2 * copt366 * copt4;
  Real copt1591 = copt26 * copt29 * copt369;
  Real copt1593 = copt356 * copt369;
  Real copt1594 = -(copt26 * copt32 * copt369);
  Real copt1596 = copt366 * copt369;
  Real copt1597 = -(copt17 * copt417);
  Real copt1598 = -(copt421 * copt8);
  Real copt1599 = 2 * copt4 * copt426;
  Real copt1600 = copt1597 + copt1598 + copt1599 + copt478;
  Real copt1601 = copt13 * copt1600;
  Real copt1602 = -(copt17 * copt384 * copt8);
  Real copt1603 = 2 * copt20 * copt384 * copt4;
  Real copt1604 = -(copt20 * copt384 * copt8);
  Real copt1605 = -2 * copt26 * copt398 * copt4;
  Real copt1606 = 2 * copt26 * copt398 * copt8;
  Real copt1607 = -(copt29 * copt398 * copt8);
  Real copt1608 = 2 * copt32 * copt398 * copt4;
  Real copt1609 = -(copt32 * copt398 * copt8);
  Real copt1610 = -(copt17 * copt426);
  Real copt1611 = -(copt20 * copt384);
  Real copt1612 = -(copt29 * copt434);
  Real copt1613 = -(copt32 * copt398);
  Real copt1614 = copt1610 + copt1611 + copt1612 + copt1613 + copt356 + copt366;
  Real copt1615 = copt1614 * copt2;
  Real copt1616 = copt1584 + copt1585 + copt1586 + copt1587 + copt1588 +
                  copt1589 + copt1590 + copt1591 + copt1592 + copt1593 +
                  copt1594 + copt1595 + copt1596 + copt1601 + copt1602 +
                  copt1603 + copt1604 + copt1605 + copt1606 + copt1607 +
                  copt1608 + copt1609 + copt1615;
  Real copt1617 = copt1616 * copt348;
  Real copt1619 = copt1618 * copt340 * copt466;
  Real copt1620 = copt1617 + copt1619;
  Real copt1621 = -(copt1312 * copt1620 * copt471 * copt495);
  Real copt1623 = copt1583 + copt1621;
  Real copt1627 = copt238 * copt562;
  Real copt1629 = copt217 * copt562;
  Real copt1630 = -(copt356 * copt519);
  Real copt1631 = -(copt366 * copt519);
  Real copt1633 = copt20 * copt732 * copt8;
  Real copt1635 = -2 * copt20 * copt519;
  Real copt1636 = copt2 * copt736;
  Real copt1637 = copt1081 + copt1634 + copt1635 + copt1636;
  Real copt1638 = -(copt13 * copt1637);
  Real copt1640 = copt1044 * copt32 * copt8;
  Real copt1642 = -2 * copt32 * copt519;
  Real copt1643 = copt1046 * copt2;
  Real copt1644 = copt1108 + copt1641 + copt1642 + copt1643;
  Real copt1645 = -(copt1644 * copt26);
  Real copt1646 = copt1625 + copt1626 + copt1627 + copt1629 + copt1630 +
                  copt1631 + copt1632 + copt1633 + copt1638 + copt1639 +
                  copt1640 + copt1645;
  Real copt1647 = -(copt1118 * copt1139 * copt1348 * copt1646 * copt509);
  Real copt1648 = copt26 * copt736;
  Real copt1649 = -(copt1044 * copt20);
  Real copt1650 = copt1135 * copt13;
  Real copt1654 = copt1508 + copt1648 + copt1649 + copt1650;
  Real copt1655 = copt1113 * copt1118 * copt1348 * copt1654 * copt509;
  Real copt1656 = copt1647 + copt1655;
  Real copt1661 = copt26 * copt312;
  Real copt1662 = -(copt2 * copt228);
  Real copt1663 = copt1660 + copt1661 + copt1662 + copt255 + copt258;
  Real copt1664 = copt1663 * copt302;
  Real copt1667 = copt1665 + copt1666;
  Real copt1668 = copt1667 * copt332;
  Real copt1669 = copt1664 + copt1668;
  Real copt1670 = -(copt1269 * copt1669 * copt181 * copt297);
  Real copt1672 = 2 * copt17 * copt182;
  Real copt1673 = -2 * copt17 * copt2 * copt8;
  Real copt1674 = -(copt182 * copt20);
  Real copt1675 = copt2 * copt20 * copt4;
  Real copt1676 = -2 * copt17 * copt197 * copt2;
  Real copt1677 = 2 * copt17 * copt197 * copt8;
  Real copt1678 = copt197 * copt2 * copt20;
  Real copt1679 = -(copt197 * copt20 * copt4);
  Real copt1680 = copt1666 + copt314 + copt342;
  Real copt1681 = copt1680 * copt217;
  Real copt1682 = -(copt182 * copt206);
  Real copt1683 = copt2 * copt206 * copt4;
  Real copt1684 = copt2 * copt206 * copt8;
  Real copt1685 = -(copt206 * copt4 * copt8);
  Real copt1686 = -(copt206 * copt29 * copt32);
  Real copt1687 = -(copt20 * copt228 * copt29);
  Real copt1688 = 2 * copt17 * copt228 * copt32;
  Real copt1692 =
      copt1254 + copt1504 + copt1522 + copt1689 + copt1690 + copt1691;
  Real copt1693 = copt1692 * copt26;
  Real copt1706 = -2 * copt228 * copt32;
  Real copt1707 = copt26 * copt261;
  Real copt1708 = copt1694 + copt1698 + copt1699 + copt1700 + copt1703 +
                  copt1705 + copt1706 + copt1707;
  Real copt1709 = copt13 * copt1708;
  Real copt1710 = copt1672 + copt1673 + copt1674 + copt1675 + copt1676 +
                  copt1677 + copt1678 + copt1679 + copt1681 + copt1682 +
                  copt1683 + copt1684 + copt1685 + copt1686 + copt1687 +
                  copt1688 + copt1693 + copt1709;
  Real copt1711 = -(copt1710 * copt181);
  Real copt1712 = copt1271 * copt173 * copt297;
  Real copt1713 = copt1711 + copt1712;
  Real copt1714 = -(copt1269 * copt1713 * copt302 * copt332);
  Real copt1715 = copt1670 + copt1714;
  Real copt1717 = -(copt32 * copt369);
  Real copt1718 = copt26 * copt479;
  Real copt1719 = -(copt2 * copt398);
  Real copt1721 = copt1660 + copt1717 + copt1718 + copt1719 + copt1720;
  Real copt1722 = copt1721 * copt471;
  Real copt1723 = copt1406 + copt1666;
  Real copt1724 = copt1723 * copt495;
  Real copt1725 = copt1722 + copt1724;
  Real copt1726 = copt1312 * copt1725 * copt348 * copt466;
  Real copt1727 = -2 * copt17 * copt350;
  Real copt1728 = 2 * copt20 * copt4 * copt8;
  Real copt1729 = 2 * copt17 * copt26 * copt32;
  Real copt1731 = -(copt20 * copt26 * copt32);
  Real copt1732 = 2 * copt20 * copt29 * copt32;
  Real copt1733 = -2 * copt17 * copt366;
  Real copt1734 = 2 * copt17 * copt369 * copt8;
  Real copt1735 = -(copt20 * copt369 * copt4);
  Real copt1736 = -(copt20 * copt369 * copt8);
  Real copt1737 = copt26 * copt29 * copt384;
  Real copt1739 = copt350 * copt384;
  Real copt1740 = -(copt26 * copt32 * copt384);
  Real copt1742 = copt366 * copt384;
  Real copt1743 = 2 * copt17 * copt417;
  Real copt1744 = 2 * copt20 * copt369;
  Real copt1745 = -(copt4 * copt426);
  Real copt1746 = -(copt452 * copt8);
  Real copt1747 = copt1743 + copt1744 + copt1745 + copt1746;
  Real copt1748 = copt1747 * copt2;
  Real copt1749 = -(copt4 * copt417);
  Real copt1750 = copt1749 + copt350 + copt432 + copt435;
  Real copt1751 = copt13 * copt1750;
  Real copt1752 = -2 * copt17 * copt26 * copt398;
  Real copt1753 = 2 * copt20 * copt26 * copt398;
  Real copt1755 = 2 * copt17 * copt32 * copt398;
  Real copt1756 = -(copt20 * copt32 * copt398);
  Real copt1757 = copt1412 + copt1727 + copt1728 + copt1729 + copt1731 +
                  copt1732 + copt1733 + copt1734 + copt1735 + copt1736 +
                  copt1737 + copt1738 + copt1739 + copt1740 + copt1741 +
                  copt1742 + copt1748 + copt1751 + copt1752 + copt1753 +
                  copt1755 + copt1756 + copt438;
  Real copt1758 = copt1757 * copt348;
  Real copt1759 = copt1618 * copt343 * copt466;
  Real copt1760 = copt1758 + copt1759;
  Real copt1761 = -(copt1312 * copt1760 * copt471 * copt495);
  Real copt1762 = copt1726 + copt1761;
  Real copt1766 = -(copt2 * copt20 * copt519);
  Real copt1767 = copt217 * copt736;
  Real copt1768 = -(copt182 * copt732);
  Real copt1769 = 2 * copt2 * copt732 * copt8;
  Real copt1771 = copt1046 * copt506;
  Real copt1772 = copt1770 + copt1771 + copt557 + copt573;
  Real copt1773 = -(copt13 * copt1772);
  Real copt1776 = copt1333 + copt1774 + copt1775;
  Real copt1777 = -(copt1776 * copt26);
  Real copt1782 = copt1070 + copt1073 + copt1076 + copt1087 + copt1764 +
                  copt1765 + copt1766 + copt1767 + copt1768 + copt1769 +
                  copt1773 + copt1777;
  Real copt1783 = -(copt1118 * copt1139 * copt1348 * copt1782 * copt509);
  Real copt1784 = -(copt32 * copt519);
  Real copt1785 = copt1126 * copt26;
  Real copt1786 = -(copt1044 * copt2);
  Real copt1787 = copt1108 + copt1660 + copt1784 + copt1785 + copt1786;
  Real copt1788 = copt1113 * copt1118 * copt1348 * copt1787 * copt509;
  Real copt1789 = copt1783 + copt1788;
  Real copt1794 = copt13 * copt325;
  Real copt1795 = copt2 * copt206;
  Real copt1796 = copt1549 + copt1793 + copt1794 + copt1795 + copt280;
  Real copt1797 = copt1796 * copt302;
  Real copt1800 = copt1798 + copt1799;
  Real copt1801 = copt1800 * copt332;
  Real copt1802 = copt1797 + copt1801;
  Real copt1803 = -(copt1269 * copt1802 * copt181 * copt297);
  Real copt1804 = 2 * copt182 * copt29;
  Real copt1805 = -2 * copt2 * copt29 * copt8;
  Real copt1806 = -(copt182 * copt32);
  Real copt1807 = copt2 * copt32 * copt4;
  Real copt1808 = -2 * copt197 * copt2 * copt29;
  Real copt1809 = 2 * copt197 * copt29 * copt8;
  Real copt1810 = copt197 * copt2 * copt32;
  Real copt1811 = -(copt197 * copt32 * copt4);
  Real copt1812 = 2 * copt20 * copt206 * copt29;
  Real copt1813 = -(copt17 * copt206 * copt32);
  Real copt1816 = copt220 * copt4;
  Real copt1818 = -2 * copt20 * copt206;
  Real copt1819 =
      copt1700 + copt1703 + copt1815 + copt1816 + copt1817 + copt1818;
  Real copt1820 = copt1819 * copt26;
  Real copt1821 = copt1248 + copt1799 + copt327;
  Real copt1822 = copt1821 * copt238;
  Real copt1823 = -(copt182 * copt228);
  Real copt1824 = copt2 * copt228 * copt4;
  Real copt1825 = copt2 * copt228 * copt8;
  Real copt1826 = -(copt228 * copt4 * copt8);
  Real copt1827 = -(copt17 * copt20 * copt228);
  Real copt1833 = copt26 * copt282;
  Real copt1835 = copt1254 + copt1318 + copt1522 + copt1828 + copt1832 +
                  copt1833 + copt1834;
  Real copt1836 = copt13 * copt1835;
  Real copt1837 = copt1804 + copt1805 + copt1806 + copt1807 + copt1808 +
                  copt1809 + copt1810 + copt1811 + copt1812 + copt1813 +
                  copt1820 + copt1822 + copt1823 + copt1824 + copt1825 +
                  copt1826 + copt1827 + copt1836;
  Real copt1838 = -(copt181 * copt1837);
  Real copt1839 = copt1271 * copt177 * copt297;
  Real copt1840 = copt1838 + copt1839;
  Real copt1841 = -(copt1269 * copt1840 * copt302 * copt332);
  Real copt1842 = copt1803 + copt1841;
  Real copt1844 = copt13 * copt417;
  Real copt1845 = copt2 * copt384;
  Real copt1846 = -(copt384 * copt8);
  Real copt1847 = copt1793 + copt1844 + copt1845 + copt1846 + copt419;
  Real copt1848 = copt1847 * copt471;
  Real copt1849 = copt1487 + copt1799;
  Real copt1850 = copt1849 * copt495;
  Real copt1851 = copt1848 + copt1850;
  Real copt1852 = copt1312 * copt1851 * copt348 * copt466;
  Real copt1854 = copt26 * copt350;
  Real copt1855 = -2 * copt29 * copt350;
  Real copt1857 = copt26 * copt356;
  Real copt1858 = -2 * copt29 * copt356;
  Real copt1859 = 2 * copt32 * copt4 * copt8;
  Real copt1860 = 2 * copt17 * copt20 * copt32;
  Real copt1861 = copt26 * copt369 * copt4;
  Real copt1863 = -(copt26 * copt369 * copt8);
  Real copt1864 = 2 * copt29 * copt369 * copt8;
  Real copt1865 = -(copt32 * copt369 * copt8);
  Real copt1866 = copt17 * copt26 * copt384;
  Real copt1867 = -(copt20 * copt26 * copt384);
  Real copt1869 = 2 * copt20 * copt29 * copt384;
  Real copt1870 = -(copt20 * copt32 * copt384);
  Real copt1872 = copt350 * copt398;
  Real copt1874 = copt356 * copt398;
  Real copt1875 = 2 * copt20 * copt29;
  Real copt1876 = -(copt20 * copt32);
  Real copt1877 = -2 * copt29 * copt384;
  Real copt1878 = 2 * copt32 * copt384;
  Real copt1879 = -(copt17 * copt434);
  Real copt1880 =
      copt1572 + copt1875 + copt1876 + copt1877 + copt1878 + copt1879;
  Real copt1881 = copt13 * copt1880;
  Real copt1882 = 2 * copt29 * copt417;
  Real copt1883 = 2 * copt32 * copt369;
  Real copt1884 = -(copt4 * copt434);
  Real copt1885 = -(copt460 * copt8);
  Real copt1886 = copt1882 + copt1883 + copt1884 + copt1885;
  Real copt1887 = copt1886 * copt2;
  Real copt1888 = copt1469 + copt1473 + copt1853 + copt1854 + copt1855 +
                  copt1856 + copt1857 + copt1858 + copt1859 + copt1860 +
                  copt1861 + copt1863 + copt1864 + copt1865 + copt1866 +
                  copt1867 + copt1869 + copt1870 + copt1871 + copt1872 +
                  copt1873 + copt1874 + copt1881 + copt1887;
  Real copt1889 = copt1888 * copt348;
  Real copt1890 = copt1618 * copt345 * copt466;
  Real copt1891 = copt1889 + copt1890;
  Real copt1892 = -(copt1312 * copt1891 * copt471 * copt495);
  Real copt1893 = copt1852 + copt1892;
  Real copt1897 = -(copt2 * copt32 * copt519);
  Real copt1900 = copt1770 + copt1898 + copt1899 + copt557 + copt573;
  Real copt1901 = -(copt1900 * copt26);
  Real copt1902 = copt1046 * copt238;
  Real copt1903 = -(copt1044 * copt182);
  Real copt1904 = 2 * copt1044 * copt2 * copt8;
  Real copt1906 = -(copt26 * copt732);
  Real copt1908 = copt1508 + copt1774 + copt1905 + copt1906 + copt1907;
  Real copt1909 = -(copt13 * copt1908);
  Real copt1910 = copt1093 + copt1095 + copt1101 + copt1103 + copt1895 +
                  copt1896 + copt1897 + copt1901 + copt1902 + copt1903 +
                  copt1904 + copt1909;
  Real copt1911 = -(copt1118 * copt1139 * copt1348 * copt1910 * copt509);
  Real copt1912 = copt13 * copt562;
  Real copt1913 = copt2 * copt732;
  Real copt1915 = -(copt732 * copt8);
  Real copt1916 = copt1079 + copt1793 + copt1912 + copt1913 + copt1915;
  Real copt1917 = copt1113 * copt1118 * copt1348 * copt1916 * copt509;
  Real copt1921 = copt1911 + copt1917;
  Real copt1927 = copt183 * copt197;
  Real copt1928 = copt189 * copt197;
  Real copt1930 = copt1929 * copt238;
  Real copt1931 = copt1929 * copt217;
  Real copt1932 = -(copt17 * copt206 * copt4);
  Real copt1934 = -2 * copt17 * copt197;
  Real copt1936 = copt1935 * copt2;
  Real copt1938 = copt1933 + copt1934 + copt1936 + copt1937;
  Real copt1939 = copt13 * copt1938;
  Real copt1940 = -(copt228 * copt29 * copt4);
  Real copt1941 = -2 * copt197 * copt29;
  Real copt1943 = copt1942 * copt2;
  Real copt1944 = copt228 + copt29;
  Real copt1945 = copt1944 * copt4;
  Real copt1946 = copt1941 + copt1943 + copt1945;
  Real copt1947 = copt1946 * copt26;
  Real copt1948 = copt1542 + copt1555 + copt1925 + copt1926 + copt1927 +
                  copt1928 + copt1930 + copt1931 + copt1932 + copt1939 +
                  copt1940 + copt1947;
  Real copt1949 = copt1269 * copt181 * copt1948 * copt302 * copt332;
  Real copt1950 = -(copt206 * copt29);
  Real copt1953 = copt1952 * copt26;
  Real copt1954 = copt13 * copt1942;
  Real copt1958 = copt1834 + copt1950 + copt1953 + copt1954;
  Real copt1959 = -(copt1269 * copt181 * copt1958 * copt297 * copt302);
  Real copt1960 = copt1949 + copt1959;
  Real copt1962 = -(copt29 * copt384);
  Real copt1963 = copt140 + copt384;
  Real copt1964 = copt1963 * copt26;
  Real copt1965 = copt29 + copt433;
  Real copt1966 = copt13 * copt1965;
  Real copt1967 = copt17 * copt398;
  Real copt1968 = copt1962 + copt1964 + copt1966 + copt1967;
  Real copt1969 = copt1968 * copt471;
  Real copt1971 = copt1701 + copt1970;
  Real copt1972 = copt1971 * copt495;
  Real copt1973 = copt1969 + copt1972;
  Real copt1974 = copt1312 * copt1973 * copt348 * copt466;
  Real copt1975 = -(copt26 * copt29 * copt4);
  Real copt1976 = -2 * copt183 * copt8;
  Real copt1977 = 2 * copt26 * copt29 * copt8;
  Real copt1978 = -2 * copt189 * copt8;
  Real copt1979 = 2 * copt17 * copt20 * copt4;
  Real copt1981 = 2 * copt29 * copt32 * copt4;
  Real copt1982 = copt183 * copt369;
  Real copt1983 = -(copt26 * copt29 * copt369);
  Real copt1985 = copt189 * copt369;
  Real copt1986 = copt26 * copt32 * copt369;
  Real copt1987 = -(copt17 * copt384 * copt4);
  Real copt1988 = 2 * copt17 * copt384 * copt8;
  Real copt1989 = copt1970 + copt416;
  Real copt1990 = copt17 * copt1989;
  Real copt1991 = copt17 + copt20 + copt420;
  Real copt1992 = -(copt1991 * copt4);
  Real copt1993 = -2 * copt384 * copt8;
  Real copt1994 = copt1990 + copt1992 + copt1993 + copt419;
  Real copt1995 = copt13 * copt1994;
  Real copt1996 = 2 * copt26 * copt398 * copt4;
  Real copt1997 = -(copt29 * copt398 * copt4);
  Real copt1998 = -2 * copt26 * copt398 * copt8;
  Real copt1999 = 2 * copt29 * copt398 * copt8;
  Real copt2000 = copt20 * copt384;
  Real copt2001 = -(copt17 * copt452);
  Real copt2003 = -(copt29 * copt460);
  Real copt2004 = copt183 + copt189 + copt2000 + copt2001 + copt2002 + copt2003;
  Real copt2005 = copt2 * copt2004;
  Real copt2006 = copt1592 + copt1595 + copt1975 + copt1976 + copt1977 +
                  copt1978 + copt1979 + copt1980 + copt1981 + copt1982 +
                  copt1983 + copt1985 + copt1986 + copt1987 + copt1988 +
                  copt1995 + copt1996 + copt1997 + copt1998 + copt1999 +
                  copt2005 + copt448 + copt456;
  Real copt2008 = copt2006 * copt348;
  Real copt2009 = -(copt1618 * copt340 * copt466);
  Real copt2010 = copt2008 + copt2009;
  Real copt2011 = -(copt1312 * copt2010 * copt471 * copt495);
  Real copt2012 = copt1974 + copt2011;
  Real copt2014 = copt140 + copt732;
  Real copt2015 = copt2014 * copt26;
  Real copt2016 = copt1045 + copt29;
  Real copt2017 = copt13 * copt2016;
  Real copt2019 = copt1507 + copt2015 + copt2017 + copt2018;
  Real copt2020 = copt1118 * copt2019;
  Real copt2021 = copt1529 + copt1970;
  Real copt2022 = copt1139 * copt2021;
  Real copt2023 = copt2020 + copt2022;
  Real copt2024 = copt1113 * copt1348 * copt2023 * copt509;
  Real copt2027 = copt17 * copt20 * copt519;
  Real copt2028 = copt29 * copt32 * copt519;
  Real copt2031 = copt2030 * copt238;
  Real copt2032 = copt2030 * copt217;
  Real copt2033 = 2 * copt17 * copt2 * copt732;
  Real copt2034 = -2 * copt17 * copt732 * copt8;
  Real copt2035 = copt20 * copt4 * copt732;
  Real copt2038 = copt17 * copt519;
  Real copt2040 = -2 * copt732 * copt8;
  Real copt2042 = copt2 * copt2041;
  Real copt2043 = copt1079 + copt2036 + copt2037 + copt2038 + copt2039 +
                  copt2040 + copt2042;
  Real copt2044 = -(copt13 * copt2043);
  Real copt2045 = 2 * copt1044 * copt2 * copt29;
  Real copt2046 = -2 * copt1044 * copt29 * copt8;
  Real copt2047 = copt1044 * copt32 * copt4;
  Real copt2052 = copt1970 + copt561;
  Real copt2053 = -(copt2052 * copt29);
  Real copt2054 = copt1044 * copt4;
  Real copt2055 = -2 * copt1044 * copt8;
  Real copt2057 = copt2 * copt2056;
  Real copt2058 =
      copt1106 + copt2048 + copt2053 + copt2054 + copt2055 + copt2057;
  Real copt2059 = -(copt2058 * copt26);
  Real copt2060 = copt1632 + copt1639 + copt2025 + copt2026 + copt2027 +
                  copt2028 + copt2031 + copt2032 + copt2033 + copt2034 +
                  copt2035 + copt2044 + copt2045 + copt2046 + copt2047 +
                  copt2059;
  Real copt2061 = copt2060 * copt509;
  Real copt2062 = -(copt1113 * copt1376 * copt502);
  Real copt2063 = copt2061 + copt2062;
  Real copt2064 = -(copt1118 * copt1139 * copt1348 * copt2063);
  Real copt2065 = copt2024 + copt2064;
  Real copt2074 = copt17 * copt197 * copt2;
  Real copt2075 = -(copt17 * copt197 * copt4);
  Real copt2076 = copt182 * copt206;
  Real copt2077 = -2 * copt2 * copt206 * copt4;
  Real copt2078 = copt206 * copt213;
  Real copt2079 = copt189 * copt206;
  Real copt2080 = copt1952 * copt217;
  Real copt2082 = -(copt17 * copt228 * copt29);
  Real copt2084 = copt1832 + copt1834 + copt2083;
  Real copt2085 = copt2084 * copt26;
  Real copt2090 = copt1942 * copt26;
  Real copt2091 =
      copt1699 + copt1705 + copt2086 + copt2087 + copt2089 + copt2090;
  Real copt2092 = copt13 * copt2091;
  Real copt2093 = copt2072 + copt2073 + copt2074 + copt2075 + copt2076 +
                  copt2077 + copt2078 + copt2079 + copt2080 + copt2082 +
                  copt2085 + copt2092;
  Real copt2094 = copt1269 * copt181 * copt2093 * copt302 * copt332;
  Real copt2097 = copt2088 * copt26;
  Real copt2098 = copt197 * copt29;
  Real copt2099 = copt2 * copt228;
  Real copt2100 = -(copt228 * copt4);
  Real copt2101 = copt2096 + copt2097 + copt2098 + copt2099 + copt2100;
  Real copt2102 = -(copt1269 * copt181 * copt2101 * copt297 * copt302);
  Real copt2103 = copt2094 + copt2102;
  Real copt2105 = copt4 + copt416;
  Real copt2106 = copt2105 * copt26;
  Real copt2107 = copt29 * copt369;
  Real copt2108 = copt2 * copt398;
  Real copt2110 = copt2096 + copt2106 + copt2107 + copt2108 + copt2109;
  Real copt2111 = copt2110 * copt471;
  Real copt2113 = copt2112 + copt281;
  Real copt2114 = copt2113 * copt495;
  Real copt2115 = copt2111 + copt2114;
  Real copt2116 = copt1312 * copt2115 * copt348 * copt466;
  Real copt2117 = -(copt17 * copt26 * copt29);
  Real copt2119 = 2 * copt17 * copt4 * copt8;
  Real copt2120 = -2 * copt20 * copt213;
  Real copt2121 = 2 * copt20 * copt26 * copt29;
  Real copt2122 = -2 * copt189 * copt20;
  Real copt2124 = 2 * copt17 * copt29 * copt32;
  Real copt2125 = -(copt17 * copt369 * copt4);
  Real copt2126 = -(copt17 * copt369 * copt8);
  Real copt2127 = 2 * copt20 * copt369 * copt4;
  Real copt2129 = copt213 * copt384;
  Real copt2130 = -(copt26 * copt29 * copt384);
  Real copt2131 = copt189 * copt384;
  Real copt2132 = copt26 * copt32 * copt384;
  Real copt2133 = 2 * copt20 * copt4;
  Real copt2137 = -(copt17 * copt2136);
  Real copt2138 = -(copt384 * copt4);
  Real copt2139 = copt2133 + copt2137 + copt2138 + copt450 + copt481;
  Real copt2140 = copt2 * copt2139;
  Real copt2141 = 2 * copt17 * copt26 * copt398;
  Real copt2142 = -(copt17 * copt29 * copt398);
  Real copt2143 = -2 * copt20 * copt26 * copt398;
  Real copt2144 = 2 * copt20 * copt29 * copt398;
  Real copt2145 = -(copt17 * copt32 * copt398);
  Real copt2146 = -(copt29 * copt32);
  Real copt2147 = copt369 * copt8;
  Real copt2148 = copt369 + copt8;
  Real copt2149 = -(copt2148 * copt4);
  Real copt2150 = -(copt29 * copt398);
  Real copt2152 =
      copt189 + copt2002 + copt213 + copt2146 + copt2147 + copt2149 + copt2150;
  Real copt2153 = copt13 * copt2152;
  Real copt2154 = copt1738 + copt1741 + copt2117 + copt2119 + copt2120 +
                  copt2121 + copt2122 + copt2123 + copt2124 + copt2125 +
                  copt2126 + copt2127 + copt2129 + copt2130 + copt2131 +
                  copt2132 + copt2140 + copt2141 + copt2142 + copt2143 +
                  copt2144 + copt2145 + copt2153;
  Real copt2155 = copt2154 * copt348;
  Real copt2156 = -(copt1618 * copt343 * copt466);
  Real copt2157 = copt2155 + copt2156;
  Real copt2158 = -(copt1312 * copt2157 * copt471 * copt495);
  Real copt2159 = copt2116 + copt2158;
  Real copt2161 = copt4 + copt561;
  Real copt2162 = copt2161 * copt26;
  Real copt2163 = copt29 * copt519;
  Real copt2164 = copt1044 * copt2;
  Real copt2166 = -(copt1044 * copt4);
  Real copt2167 = copt2096 + copt2162 + copt2163 + copt2164 + copt2166;
  Real copt2168 = copt1118 * copt2167;
  Real copt2171 = copt1665 + copt2112;
  Real copt2172 = copt1139 * copt2171;
  Real copt2173 = copt2168 + copt2172;
  Real copt2174 = copt1113 * copt1348 * copt2173 * copt509;
  Real copt2175 = copt17 * copt182;
  Real copt2176 = -(copt17 * copt2 * copt8);
  Real copt2177 = -2 * copt182 * copt20;
  Real copt2178 = 2 * copt2 * copt20 * copt4;
  Real copt2179 = -(copt17 * copt2 * copt519);
  Real copt2180 = 2 * copt2 * copt20 * copt519;
  Real copt2181 = copt182 * copt732;
  Real copt2182 = -(copt2 * copt4 * copt732);
  Real copt2183 = -(copt2 * copt732 * copt8);
  Real copt2184 = copt2041 * copt217;
  Real copt2185 = copt1044 * copt17 * copt32;
  Real copt2186 = copt2112 + copt735;
  Real copt2187 = -(copt2186 * copt29);
  Real copt2188 = copt1318 + copt1508 + copt1907 + copt2018 + copt2187;
  Real copt2189 = -(copt2188 * copt26);
  Real copt2190 = copt26 * copt29;
  Real copt2193 = copt1044 * copt26;
  Real copt2194 = -2 * copt1044 * copt29;
  Real copt2195 = copt1044 * copt32;
  Real copt2196 = copt1117 + copt1694 + copt1698 + copt2190 + copt2191 +
                  copt2192 + copt2193 + copt2194 + copt2195 + copt573;
  Real copt2197 = -(copt13 * copt2196);
  Real copt2198 = copt1067 + copt1069 + copt1072 + copt1075 + copt1086 +
                  copt2175 + copt2176 + copt2177 + copt2178 + copt2179 +
                  copt2180 + copt2181 + copt2182 + copt2183 + copt2184 +
                  copt2185 + copt2189 + copt2197;
  Real copt2199 = copt2198 * copt509;
  Real copt2201 = -(copt1113 * copt1376 * copt504);
  Real copt2202 = copt2199 + copt2201;
  Real copt2203 = -(copt1118 * copt1139 * copt1348 * copt2202);
  Real copt2204 = copt2174 + copt2203;
  Real copt2210 = copt197 * copt2 * copt29;
  Real copt2211 = -(copt197 * copt29 * copt4);
  Real copt2212 = -(copt17 * copt206 * copt29);
  Real copt2214 = copt1699 + copt1817 + copt2086 + copt2089 + copt2213;
  Real copt2216 = copt2214 * copt26;
  Real copt2217 = copt182 * copt228;
  Real copt2218 = -2 * copt2 * copt228 * copt4;
  Real copt2219 = copt213 * copt228;
  Real copt2220 = copt183 * copt228;
  Real copt2221 = copt175 + copt228;
  Real copt2222 = copt2221 * copt238;
  Real copt2224 = copt1935 * copt26;
  Real copt2225 = copt1690 + copt1691 + copt2083 + copt2224;
  Real copt2226 = copt13 * copt2225;
  Real copt2227 = copt2208 + copt2209 + copt2210 + copt2211 + copt2212 +
                  copt2216 + copt2217 + copt2218 + copt2219 + copt2220 +
                  copt2222 + copt2226;
  Real copt2228 = copt1269 * copt181 * copt2227 * copt302 * copt332;
  Real copt2232 = -(copt17 * copt197);
  Real copt2234 = copt13 * copt1929;
  Real copt2255 = -(copt2 * copt206);
  Real copt2265 = copt1937 + copt2230 + copt2232 + copt2234 + copt2255;
  Real copt2287 = -(copt1269 * copt181 * copt2265 * copt297 * copt302);
  Real copt2289 = copt2228 + copt2287;
  Real copt2313 = -(copt17 * copt369);
  Real copt2314 = copt134 + copt369;
  Real copt2315 = copt13 * copt2314;
  Real copt2328 = -(copt2 * copt384);
  Real copt2329 = copt384 * copt4;
  Real copt2330 = copt2230 + copt2313 + copt2315 + copt2328 + copt2329;
  Real copt2349 = copt2330 * copt471;
  Real copt2367 = copt2357 + copt260;
  Real copt2383 = copt2367 * copt495;
  Real copt2413 = copt2349 + copt2383;
  Real copt2414 = copt1312 * copt2413 * copt348 * copt466;
  Real copt2429 = copt213 * copt26;
  Real copt2430 = copt183 * copt26;
  Real copt2431 = 2 * copt29 * copt4 * copt8;
  Real copt2445 = 2 * copt17 * copt20 * copt29;
  Real copt2477 = -2 * copt213 * copt32;
  Real copt2502 = -2 * copt183 * copt32;
  Real copt2521 = -(copt26 * copt369 * copt4);
  Real copt2522 = -(copt29 * copt369 * copt4);
  Real copt2523 = copt26 * copt369 * copt8;
  Real copt2561 = 2 * copt32 * copt369 * copt4;
  Real copt2592 = -(copt17 * copt26 * copt384);
  Real copt2605 = -(copt17 * copt29 * copt384);
  Real copt2606 = copt20 * copt26 * copt384;
  Real copt2607 = 2 * copt17 * copt32 * copt384;
  Real copt2648 = copt213 * copt398;
  Real copt2677 = copt183 * copt398;
  Real copt2691 = 2 * copt32 * copt4;
  Real copt2692 = -(copt2136 * copt29);
  Real copt2693 = copt1720 + copt2109 + copt2691 + copt2692 + copt458;
  Real copt2732 = copt2 * copt2693;
  Real copt2762 = 2 * copt29 * copt384;
  Real copt2775 = -2 * copt32 * copt384;
  Real copt2776 = copt175 + copt2357 + copt433;
  Real copt2777 = copt17 * copt2776;
  Real copt2814 = copt1317 + copt1322 + copt2762 + copt2775 + copt2777;
  Real copt2845 = copt13 * copt2814;
  Real copt2856 = copt1468 + copt1472 + copt1853 + copt1856 + copt1871 +
                  copt1873 + copt2429 + copt2430 + copt2431 + copt2445 +
                  copt2477 + copt2502 + copt2521 + copt2522 + copt2523 +
                  copt2561 + copt2592 + copt2605 + copt2606 + copt2607 +
                  copt2648 + copt2677 + copt2732 + copt2845;
  Real copt2858 = copt2856 * copt348;
  Real copt2904 = -(copt1618 * copt345 * copt466);
  Real copt2925 = copt2858 + copt2904;
  Real copt2949 = -(copt1312 * copt2925 * copt471 * copt495);
  Real copt2950 = copt2414 + copt2949;
  Real copt2982 = -(copt17 * copt519);
  Real copt2998 = copt134 + copt519;
  Real copt3026 = copt13 * copt2998;
  Real copt3027 = -(copt2 * copt732);
  Real copt3028 = copt2039 + copt2230 + copt2982 + copt3026 + copt3027;
  Real copt3057 = copt1118 * copt3028;
  Real copt3074 = copt1798 + copt2357;
  Real copt3126 = copt1139 * copt3074;
  Real copt3127 = copt3057 + copt3126;
  Real copt3128 = copt1113 * copt1348 * copt3127 * copt509;
  Real copt3130 = copt182 * copt29;
  Real copt3175 = -(copt2 * copt29 * copt8);
  Real copt3198 = -2 * copt182 * copt32;
  Real copt3199 = 2 * copt2 * copt32 * copt4;
  Real copt3200 = -(copt2 * copt29 * copt519);
  Real copt3237 = copt29 * copt519 * copt8;
  Real copt3246 = 2 * copt2 * copt32 * copt519;
  Real copt3301 = copt20 * copt29 * copt732;
  Real copt3302 = -2 * copt17 * copt732;
  Real copt3303 =
      copt1694 + copt1815 + copt1899 + copt2191 + copt2192 + copt3302 + copt573;
  Real copt3305 = -(copt26 * copt3303);
  Real copt3348 = copt1044 * copt182;
  Real copt3369 = -(copt1044 * copt2 * copt4);
  Real copt3370 = -(copt1044 * copt2 * copt8);
  Real copt3403 = copt2056 * copt238;
  Real copt3414 = -2 * copt20 * copt26;
  Real copt3463 = copt17 * copt506;
  Real copt3464 = copt26 * copt732;
  Real copt3465 = -(copt1046 * copt17);
  Real copt3467 = copt1330 + copt1333 + copt1504 + copt1775 + copt3414 +
                  copt3463 + copt3464 + copt3465;
  Real copt3504 = -(copt13 * copt3467);
  Real copt3526 = copt1092 + copt1094 + copt1100 + copt1102 + copt3130 +
                  copt3175 + copt3198 + copt3199 + copt3200 + copt3237 +
                  copt3246 + copt3301 + copt3305 + copt3348 + copt3369 +
                  copt3370 + copt3403 + copt3504;
  Real copt3527 = copt3526 * copt509;
  Real copt3528 = -(copt1113 * copt1376 * copt506);
  Real copt3545 = copt3527 + copt3528;
  Real copt3566 = -(copt1118 * copt1139 * copt1348 * copt3545);
  Real copt3579 = copt3128 + copt3566;
  Real copt3708 = copt2 * copt343;
  Real copt3730 = copt1933 + copt2036 + copt2037 + copt3708;
  Real copt3740 = copt13 * copt3730;
  Real copt3764 = -2 * copt29 * copt8;
  Real copt3770 = copt2 * copt345;
  Real copt3784 = copt29 + copt32;
  Real copt3785 = copt3784 * copt4;
  Real copt3786 = copt3764 + copt3770 + copt3785;
  Real copt3792 = copt26 * copt3786;
  Real copt3804 = copt1535 + copt1536 + copt1925 + copt1926 + copt3637 +
                  copt3646 + copt3675 + copt3677 + copt3702 + copt3740 +
                  copt3761 + copt3792;
  Real copt3808 = copt1269 * copt181 * copt302 * copt332 * copt3804;
  Real copt3836 = -(copt1269 * copt181 * copt297 * copt302 * copt3835);
  Real copt3837 = copt3808 + copt3836;
  Real copt3854 = copt17 * copt2 * copt8;
  Real copt3858 = -2 * copt2 * copt20 * copt4;
  Real copt3880 = copt1318 + copt1828 + copt2083;
  Real copt3881 = copt26 * copt3880;
  Real copt3900 = copt26 * copt345;
  Real copt3904 =
      copt1694 + copt1698 + copt2086 + copt2087 + copt3894 + copt3900;
  Real copt3916 = copt13 * copt3904;
  Real copt3920 = copt1764 + copt2072 + copt2073 + copt3854 + copt3855 +
                  copt3858 + copt3874 + copt3876 + copt3879 + copt3881 +
                  copt3916 + copt413;
  Real copt3921 = copt1269 * copt181 * copt302 * copt332 * copt3920;
  Real copt3955 = -(copt1269 * copt181 * copt297 * copt302 * copt3943);
  Real copt3959 = copt3921 + copt3955;
  Real copt3972 = copt2 * copt29 * copt8;
  Real copt3974 = copt1694 + copt1815 + copt2086 + copt2213 + copt3894;
  Real copt3986 = copt26 * copt3974;
  Real copt3987 = -2 * copt2 * copt32 * copt4;
  Real copt3993 = copt1504 + copt1689 + copt2083 + copt3815;
  Real copt4004 = copt13 * copt3993;
  Real copt4011 = copt1459 + copt1461 + copt1463 + copt1464 + copt1895 +
                  copt2208 + copt2209 + copt3972 + copt3986 + copt3987 +
                  copt3988 + copt4004;
  Real copt4018 = copt1269 * copt181 * copt302 * copt332 * copt4011;
  Real copt4046 = -(copt1269 * copt181 * copt297 * copt302 * copt4045);
  Real copt4047 = copt4018 + copt4046;
  Real copt4064 = copt26 * copt29 * copt4;
  Real copt4071 = -(copt17 * copt8);
  Real copt4072 = -(copt3875 * copt4);
  Real copt4073 = copt1634 + copt4071 + copt4072;
  Real copt4079 = copt13 * copt4073;
  Real copt4089 = copt26 * copt32 * copt8;
  Real copt4101 = 2 * copt17 * copt20;
  Real copt4104 =
      copt1898 + copt2087 + copt2213 + copt4101 + copt4102 + copt4103;
  Real copt4115 = copt2 * copt4104;
  Real copt4122 = copt1351 + copt1354 + copt1584 + copt1980 + copt3637 +
                  copt3646 + copt3702 + copt3761 + copt4064 + copt4079 +
                  copt4089 + copt4115 + copt442 + copt443;
  Real copt4124 = -(copt1312 * copt348 * copt4122 * copt471 * copt495);
  Real copt4128 = copt1312 * copt348 * copt3835 * copt466 * copt471;
  Real copt4135 = copt4124 + copt4128;
  Real copt4154 = copt17 * copt26 * copt29;
  Real copt4155 = -(copt17 * copt3656);
  Real copt4156 = copt1634 + copt4041 + copt4155;
  Real copt4170 = copt2 * copt4156;
  Real copt4171 = copt20 * copt26 * copt32;
  Real copt4172 = copt17 * copt366;
  Real copt4175 = 2 * copt4 * copt8;
  Real copt4187 =
      copt2086 + copt2087 + copt4102 + copt4103 + copt4175 + copt557;
  Real copt4188 = copt13 * copt4187;
  Real copt4200 = copt1411 + copt1412 + copt1413 + copt2123 + copt3855 +
                  copt3874 + copt3879 + copt413 + copt414 + copt4154 +
                  copt4170 + copt4171 + copt4172 + copt4188;
  Real copt4212 = -(copt1312 * copt348 * copt4200 * copt471 * copt495);
  Real copt4213 = copt1312 * copt348 * copt3943 * copt466 * copt471;
  Real copt4214 = copt4212 + copt4213;
  Real copt4224 = -(copt213 * copt26);
  Real copt4225 = -(copt183 * copt26);
  Real copt4233 = 2 * copt26 * copt4 * copt8;
  Real copt4241 = -(copt26 * copt350);
  Real copt4245 = 2 * copt17 * copt20 * copt26;
  Real copt4252 = -(copt26 * copt356);
  Real copt4262 = -(copt29 * copt3656);
  Real copt4270 = -(copt32 * copt4);
  Real copt4276 = copt1641 + copt4262 + copt4270;
  Real copt4278 = copt2 * copt4276;
  Real copt4286 = copt17 * copt345;
  Real copt4291 = copt1317 + copt1774 + copt4286;
  Real copt4303 = copt13 * copt4291;
  Real copt4316 = copt1459 + copt1460 + copt1461 + copt1462 + copt1463 +
                  copt1464 + copt1465 + copt1466 + copt4224 + copt4225 +
                  copt4233 + copt4241 + copt4245 + copt4252 + copt4278 +
                  copt4303;
  Real copt4320 = -(copt1312 * copt348 * copt4316 * copt471 * copt495);
  Real copt4332 = copt1312 * copt348 * copt4045 * copt466 * copt471;
  Real copt4355 = copt4320 + copt4332;
  Real copt4375 = copt17 * copt20 * copt8;
  Real copt4376 = -(copt356 * copt4);
  Real copt4377 = copt2 * copt3875;
  Real copt4385 = copt1546 + copt1547 + copt1634 + copt4377;
  Real copt4391 = -(copt13 * copt4385);
  Real copt4402 = copt29 * copt32 * copt8;
  Real copt4403 = -(copt366 * copt4);
  Real copt4404 = copt29 * copt8;
  Real copt4420 = -2 * copt32 * copt4;
  Real copt4424 = copt2 * copt3826;
  Real copt4440 = copt1641 + copt4404 + copt4420 + copt4424;
  Real copt4441 = -(copt26 * copt4440);
  Real copt4442 = copt1625 + copt1626 + copt2025 + copt2026 + copt3675 +
                  copt3677 + copt4375 + copt4376 + copt4391 + copt4402 +
                  copt4403 + copt4441;
  Real copt4463 = -(copt1118 * copt1139 * copt1348 * copt4442 * copt509);
  Real copt4478 = copt1113 * copt1118 * copt1348 * copt3835 * copt509;
  Real copt4495 = copt4463 + copt4478;
  Real copt4497 = 2 * copt17 * copt2 * copt8;
  Real copt4518 = -(copt2 * copt20 * copt4);
  Real copt4534 = -(copt17 * copt366);
  Real copt4551 = copt1504 + copt1689 + copt1774;
  Real copt4552 = -(copt26 * copt4551);
  Real copt4553 = -(copt26 * copt29);
  Real copt4586 = copt26 * copt32;
  Real copt4598 =
      copt1694 + copt1698 + copt4103 + copt4553 + copt4572 + copt4586 + copt557;
  Real copt4599 = -(copt13 * copt4598);
  Real copt4600 = copt1062 + copt1064 + copt1066 + copt1764 + copt1765 +
                  copt2072 + copt3876 + copt4497 + copt4518 + copt4534 +
                  copt4552 + copt4599;
  Real copt4616 = -(copt1118 * copt1139 * copt1348 * copt4600 * copt509);
  Real copt4633 = copt1113 * copt1118 * copt1348 * copt3943 * copt509;
  Real copt4655 = copt4616 + copt4633;
  Real copt4657 = 2 * copt2 * copt29 * copt8;
  Real copt4673 = -(copt29 * copt350);
  Real copt4689 = -(copt29 * copt356);
  Real copt4711 = copt1694 + copt1815 + copt1898 + copt4572 + copt557;
  Real copt4712 = -(copt26 * copt4711);
  Real copt4713 = -(copt2 * copt32 * copt4);
  Real copt4725 = -(copt17 * copt506);
  Real copt4739 = copt1774 + copt1828 + copt1905 + copt4725;
  Real copt4758 = -(copt13 * copt4739);
  Real copt4759 = copt1090 + copt1091 + copt1895 + copt1896 + copt2208 +
                  copt3988 + copt4657 + copt4673 + copt4689 + copt4712 +
                  copt4713 + copt4758;
  Real copt4760 = -(copt1118 * copt1139 * copt1348 * copt4759 * copt509);
  Real copt4776 = copt1113 * copt1118 * copt1348 * copt4045 * copt509;
  Real copt4793 = copt4760 + copt4776;
  out1(0) = copt12 + copt24 + copt38;
  out1(1) = copt11 * copt50 + copt23 * copt57 + copt36 * copt63;
  out1(2) = copt70 + copt72 + copt73;
  out1(3) =
      (copt120 * copt121 * copt85 * copt86 *
       (copt1141 * copt500 * l0 * l1 + copt337 * copt497 * l0 * l2 +
        copt128 * copt334 * l1 * l2 - copt128 * l1 * l2 * thetarest0 -
        copt337 * l0 * l2 * thetarest1 - copt500 * l0 * l1 * thetarest2)) /
      2.;
  out1(4) = (copt120 * copt121 * copt85 * copt86 *
             (copt1141 * copt1151 * copt499 * l0 * l1 +
              copt1148 * copt336 * copt497 * l0 * l2 +
              copt1145 * copt123 * copt334 * l1 * l2 -
              copt1145 * copt123 * l1 * l2 * thetarest0 -
              copt1148 * copt336 * l0 * l2 * thetarest1 -
              copt1151 * copt499 * l0 * l1 * thetarest2)) /
            2.;
  out1(5) =
      (copt120 * copt121 * copt85 * copt86 *
       (copt1141 * copt1162 * l0 * l1 + copt1159 * copt497 * l0 * l2 +
        copt1156 * copt334 * l1 * l2 - copt1156 * l1 * l2 * thetarest0 -
        copt1159 * l0 * l2 * thetarest1 - copt1162 * l0 * l1 * thetarest2)) /
      2.;
  out2(0, 0) = 2 * copt1167 * copt1170;
  out2(0, 1) = 2 * copt1167 * copt1174;
  out2(0, 2) = 2 * copt1167 * copt1178;
  out2(0, 3) = 2 * copt1 * copt11;
  out2(0, 4) = 2 * copt1 * copt23;
  out2(0, 5) = 2 * copt1 * copt36;
  out2(0, 6) = 2 * copt11 * copt7;
  out2(0, 7) = 2 * copt23 * copt7;
  out2(0, 8) = 2 * copt36 * copt7;
  out2(0, 9) = 0;
  out2(0, 10) = 0;
  out2(0, 11) = 0;
  out2(0, 12) = 0;
  out2(0, 13) = 0;
  out2(0, 14) = 0;
  out2(0, 15) = 0;
  out2(0, 16) = 0;
  out2(0, 17) = 0;
  out2(1, 0) = copt1170 * copt1186 + copt1167 * copt1190;
  out2(1, 1) = copt1174 * copt1186 + copt1167 * copt1196;
  out2(1, 2) = copt1178 * copt1186 + copt1167 * copt1205;
  out2(1, 3) = 2 * copt1 * copt42 * copt5 + copt1 * copt45 * copt9 +
               copt42 * copt7 * copt9;
  out2(1, 4) = 2 * copt1 * copt18 * copt42 + copt1 * copt21 * copt45 +
               copt21 * copt42 * copt7;
  out2(1, 5) = 2 * copt1 * copt30 * copt42 + copt1 * copt34 * copt45 +
               copt34 * copt42 * copt7;
  out2(1, 6) = copt42 * copt5 * copt7 + copt45 * (copt6 + 2 * copt7 * copt9);
  out2(1, 7) = copt18 * copt42 * copt7 + copt45 * (copt19 + 2 * copt21 * copt7);
  out2(1, 8) = copt30 * copt42 * copt7 + copt45 * (copt31 + 2 * copt34 * copt7);
  out2(1, 9) = 0;
  out2(1, 10) = 0;
  out2(1, 11) = 0;
  out2(1, 12) = 0;
  out2(1, 13) = 0;
  out2(1, 14) = 0;
  out2(1, 15) = 0;
  out2(1, 16) = 0;
  out2(1, 17) = 0;
  out2(2, 0) = 2 * copt1186 * copt1190;
  out2(2, 1) = 2 * copt1186 * copt1196;
  out2(2, 2) = 2 * copt1186 * copt1205;
  out2(2, 3) = 2 * copt42 * copt50;
  out2(2, 4) = 2 * copt42 * copt57;
  out2(2, 5) = 2 * copt42 * copt63;
  out2(2, 6) = 2 * copt45 * copt50;
  out2(2, 7) = 2 * copt45 * copt57;
  out2(2, 8) = 2 * copt45 * copt63;
  out2(2, 9) = 0;
  out2(2, 10) = 0;
  out2(2, 11) = 0;
  out2(2, 12) = 0;
  out2(2, 13) = 0;
  out2(2, 14) = 0;
  out2(2, 15) = 0;
  out2(2, 16) = 0;
  out2(2, 17) = 0;
  out2(3, 0) = (copt120 * copt121 * copt85 * copt86 *
                (copt1380 * copt500 * l0 * l1 + copt1328 * copt337 * l0 * l2 +
                 copt128 * copt1304 * l1 * l2)) /
               2.;
  out2(3, 1) = (copt120 * copt121 * copt85 * copt86 *
                (copt1436 * copt500 * l0 * l1 + copt1403 * copt337 * l0 * l2 +
                 copt128 * copt1399 * l1 * l2)) /
               2.;
  out2(3, 2) = (copt120 * copt121 * copt85 * copt86 *
                (copt1517 * copt500 * l0 * l1 + copt1484 * copt337 * l0 * l2 +
                 copt128 * copt1457 * l1 * l2)) /
               2.;
  out2(3, 3) = (copt120 * copt121 * copt85 * copt86 *
                (copt1656 * copt500 * l0 * l1 + copt1623 * copt337 * l0 * l2 +
                 copt128 * copt1568 * l1 * l2)) /
               2.;
  out2(3, 4) = (copt120 * copt121 * copt85 * copt86 *
                (copt1789 * copt500 * l0 * l1 + copt1762 * copt337 * l0 * l2 +
                 copt128 * copt1715 * l1 * l2)) /
               2.;
  out2(3, 5) = (copt120 * copt121 * copt85 * copt86 *
                (copt1921 * copt500 * l0 * l1 + copt1893 * copt337 * l0 * l2 +
                 copt128 * copt1842 * l1 * l2)) /
               2.;
  out2(3, 6) = (copt120 * copt121 * copt85 * copt86 *
                (copt2065 * copt500 * l0 * l1 + copt2012 * copt337 * l0 * l2 +
                 copt128 * copt1960 * l1 * l2)) /
               2.;
  out2(3, 7) = (copt120 * copt121 * copt85 * copt86 *
                (copt2204 * copt500 * l0 * l1 + copt2159 * copt337 * l0 * l2 +
                 copt128 * copt2103 * l1 * l2)) /
               2.;
  out2(3, 8) = (copt120 * copt121 * copt85 * copt86 *
                (copt3579 * copt500 * l0 * l1 + copt2950 * copt337 * l0 * l2 +
                 copt128 * copt2289 * l1 * l2)) /
               2.;
  out2(3, 9) = (copt128 * copt3837 * copt85 * copt86) / 2.;
  out2(3, 10) = (copt128 * copt3959 * copt85 * copt86) / 2.;
  out2(3, 11) = (copt128 * copt4047 * copt85 * copt86) / 2.;
  out2(3, 12) = (copt120 * copt337 * copt4135 * copt85) / 2.;
  out2(3, 13) = (copt120 * copt337 * copt4214 * copt85) / 2.;
  out2(3, 14) = (copt120 * copt337 * copt4355 * copt85) / 2.;
  out2(3, 15) = (copt121 * copt4495 * copt500 * copt85) / 2.;
  out2(3, 16) = (copt121 * copt4655 * copt500 * copt85) / 2.;
  out2(3, 17) = (copt121 * copt4793 * copt500 * copt85) / 2.;
  out2(4, 0) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt1380 * copt499 * l0 * l1 +
                 copt1148 * copt1328 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt1304 * l1 * l2)) /
               2.;
  out2(4, 1) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt1436 * copt499 * l0 * l1 +
                 copt1148 * copt1403 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt1399 * l1 * l2)) /
               2.;
  out2(4, 2) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt1517 * copt499 * l0 * l1 +
                 copt1148 * copt1484 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt1457 * l1 * l2)) /
               2.;
  out2(4, 3) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt1656 * copt499 * l0 * l1 +
                 copt1148 * copt1623 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt1568 * l1 * l2)) /
               2.;
  out2(4, 4) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt1789 * copt499 * l0 * l1 +
                 copt1148 * copt1762 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt1715 * l1 * l2)) /
               2.;
  out2(4, 5) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt1921 * copt499 * l0 * l1 +
                 copt1148 * copt1893 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt1842 * l1 * l2)) /
               2.;
  out2(4, 6) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt2065 * copt499 * l0 * l1 +
                 copt1148 * copt2012 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt1960 * l1 * l2)) /
               2.;
  out2(4, 7) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt2204 * copt499 * l0 * l1 +
                 copt1148 * copt2159 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt2103 * l1 * l2)) /
               2.;
  out2(4, 8) = (copt120 * copt121 * copt85 * copt86 *
                (copt1151 * copt3579 * copt499 * l0 * l1 +
                 copt1148 * copt2950 * copt336 * l0 * l2 +
                 copt1145 * copt123 * copt2289 * l1 * l2)) /
               2.;
  out2(4, 9) = (copt1145 * copt123 * copt3837 * copt85 * copt86) / 2.;
  out2(4, 10) = (copt1145 * copt123 * copt3959 * copt85 * copt86) / 2.;
  out2(4, 11) = (copt1145 * copt123 * copt4047 * copt85 * copt86) / 2.;
  out2(4, 12) = (copt1148 * copt120 * copt336 * copt4135 * copt85) / 2.;
  out2(4, 13) = (copt1148 * copt120 * copt336 * copt4214 * copt85) / 2.;
  out2(4, 14) = (copt1148 * copt120 * copt336 * copt4355 * copt85) / 2.;
  out2(4, 15) = (copt1151 * copt121 * copt4495 * copt499 * copt85) / 2.;
  out2(4, 16) = (copt1151 * copt121 * copt4655 * copt499 * copt85) / 2.;
  out2(4, 17) = (copt1151 * copt121 * copt4793 * copt499 * copt85) / 2.;
  out2(5, 0) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt1380 * l0 * l1 + copt1159 * copt1328 * l0 * l2 +
                 copt1156 * copt1304 * l1 * l2)) /
               2.;
  out2(5, 1) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt1436 * l0 * l1 + copt1159 * copt1403 * l0 * l2 +
                 copt1156 * copt1399 * l1 * l2)) /
               2.;
  out2(5, 2) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt1517 * l0 * l1 + copt1159 * copt1484 * l0 * l2 +
                 copt1156 * copt1457 * l1 * l2)) /
               2.;
  out2(5, 3) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt1656 * l0 * l1 + copt1159 * copt1623 * l0 * l2 +
                 copt1156 * copt1568 * l1 * l2)) /
               2.;
  out2(5, 4) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt1789 * l0 * l1 + copt1159 * copt1762 * l0 * l2 +
                 copt1156 * copt1715 * l1 * l2)) /
               2.;
  out2(5, 5) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt1921 * l0 * l1 + copt1159 * copt1893 * l0 * l2 +
                 copt1156 * copt1842 * l1 * l2)) /
               2.;
  out2(5, 6) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt2065 * l0 * l1 + copt1159 * copt2012 * l0 * l2 +
                 copt1156 * copt1960 * l1 * l2)) /
               2.;
  out2(5, 7) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt2204 * l0 * l1 + copt1159 * copt2159 * l0 * l2 +
                 copt1156 * copt2103 * l1 * l2)) /
               2.;
  out2(5, 8) = (copt120 * copt121 * copt85 * copt86 *
                (copt1162 * copt3579 * l0 * l1 + copt1159 * copt2950 * l0 * l2 +
                 copt1156 * copt2289 * l1 * l2)) /
               2.;
  out2(5, 9) = (copt1156 * copt3837 * copt85 * copt86) / 2.;
  out2(5, 10) = (copt1156 * copt3959 * copt85 * copt86) / 2.;
  out2(5, 11) = (copt1156 * copt4047 * copt85 * copt86) / 2.;
  out2(5, 12) = (copt1159 * copt120 * copt4135 * copt85) / 2.;
  out2(5, 13) = (copt1159 * copt120 * copt4214 * copt85) / 2.;
  out2(5, 14) = (copt1159 * copt120 * copt4355 * copt85) / 2.;
  out2(5, 15) = (copt1162 * copt121 * copt4495 * copt85) / 2.;
  out2(5, 16) = (copt1162 * copt121 * copt4655 * copt85) / 2.;
  out2(5, 17) = (copt1162 * copt121 * copt4793 * copt85) / 2.;
  return std::make_tuple(grad, val);
}
