#include "strain.hpp"
#ifndef hylc_strain_II

using namespace hylc;
using namespace hylc::mathematica;

std::tuple<Mat6x18, Vec6> hylc::mathematica::strain_valgrad(
    const Vec18 &xloc, const Mat2x2 &invDm, const Real &A,
    const Real &thetarest0, const Real &thetarest1, const Real &thetarest2,
    const Real &l0, const Real &l1, const Real &l2, const Vec2 &t0,
    const Vec2 &t1, const Vec2 &t2) {
  // define output
  Mat6x18 grad(0);
  Vec6 val(0);
  auto out1 = [&](int i) -> Real & { return val[i]; };
  auto out2 = [&](int i, int j) -> Real & { return grad(i, j); };

  Real copt41  = invDm(0, 0);
  Real copt45  = xloc(0);
  Real copt49  = -copt45;
  Real copt50  = xloc(3);
  Real copt57  = copt49 + copt50;
  Real copt64  = copt41 * copt57;
  Real copt162 = invDm(1, 0);
  Real copt165 = xloc(6);
  Real copt232 = copt165 + copt49;
  Real copt235 = copt162 * copt232;
  Real copt300 = copt235 + copt64;
  Real copt301 = Power(copt300, 2);
  Real copt302 = xloc(1);
  Real copt304 = -copt302;
  Real copt305 = xloc(4);
  Real copt307 = copt304 + copt305;
  Real copt308 = copt307 * copt41;
  Real copt310 = xloc(7);
  Real copt311 = copt304 + copt310;
  Real copt312 = copt162 * copt311;
  Real copt313 = copt308 + copt312;
  Real copt315 = Power(copt313, 2);
  Real copt316 = xloc(2);
  Real copt318 = -copt316;
  Real copt319 = xloc(5);
  Real copt321 = copt318 + copt319;
  Real copt322 = copt321 * copt41;
  Real copt323 = xloc(8);
  Real copt324 = copt318 + copt323;
  Real copt325 = copt162 * copt324;
  Real copt326 = copt322 + copt325;
  Real copt327 = Power(copt326, 2);
  Real copt328 = copt301 + copt315 + copt327;
  Real copt329 = Sqrt(copt328);
  Real copt330 = 1 / copt329;
  Real copt331 = invDm(0, 1);
  Real copt333 = invDm(1, 1);
  Real copt332 = copt331 * copt57;
  Real copt334 = copt232 * copt333;
  Real copt335 = copt332 + copt334;
  Real copt346 = Power(copt335, 2);
  Real copt337 = copt307 * copt331;
  Real copt338 = copt311 * copt333;
  Real copt339 = copt337 + copt338;
  Real copt347 = Power(copt339, 2);
  Real copt341 = copt321 * copt331;
  Real copt342 = copt324 * copt333;
  Real copt343 = copt341 + copt342;
  Real copt348 = Power(copt343, 2);
  Real copt349 = copt346 + copt347 + copt348;
  Real copt350 = Sqrt(copt349);
  Real copt351 = 1 / copt350;
  Real copt353 = 1 / A;
  Real copt354 = 1 / l0;
  Real copt355 = 1 / l1;
  Real copt356 = 1 / l2;
  Real copt357 = t0(0);
  Real copt358 = Power(copt357, 2);
  Real copt369 = xloc(9);
  Real copt374 = xloc(10);
  Real copt361 = -copt50;
  Real copt393 = xloc(11);
  Real copt376 = -copt305;
  Real copt386 = -copt319;
  Real copt418 = copt361 + copt45;
  Real copt419 = Power(copt418, 2);
  Real copt420 = copt302 + copt376;
  Real copt421 = Power(copt420, 2);
  Real copt422 = copt316 + copt386;
  Real copt423 = Power(copt422, 2);
  Real copt424 = copt419 + copt421 + copt423;
  Real copt425 = Sqrt(copt424);
  Real copt381 = -copt165;
  Real copt408 = -copt374;
  Real copt370 = -copt369;
  Real copt403 = -copt323;
  Real copt458 = t1(0);
  Real copt459 = Power(copt458, 2);
  Real copt382 = copt381 + copt50;
  Real copt401 = copt310 + copt376;
  Real copt466 = xloc(12);
  Real copt470 = xloc(13);
  Real copt383 = copt316 * copt382;
  Real copt384 = copt165 * copt319;
  Real copt385 = -(copt323 * copt50);
  Real copt387 = copt323 + copt386;
  Real copt388 = copt387 * copt45;
  Real copt389 = copt383 + copt384 + copt385 + copt388;
  Real copt468 = copt381 + copt466;
  Real copt479 = xloc(14);
  Real copt364 = -copt310;
  Real copt365 = copt305 + copt364;
  Real copt480 = -copt479;
  Real copt481 = copt323 + copt480;
  Real copt499 = Power(copt382, 2);
  Real copt500 = Power(copt365, 2);
  Real copt404 = copt319 + copt403;
  Real copt501 = Power(copt404, 2);
  Real copt502 = copt499 + copt500 + copt501;
  Real copt503 = Sqrt(copt502);
  Real copt426 = -(copt310 * copt319 * copt45);
  Real copt427 = copt305 * copt323 * copt45;
  Real copt467 = -(copt310 * copt466);
  Real copt469 = copt305 * copt468;
  Real copt471 = -copt470;
  Real copt472 = copt310 + copt471;
  Real copt473 = copt472 * copt50;
  Real copt474 = copt165 * copt470;
  Real copt475 = copt467 + copt469 + copt473 + copt474;
  Real copt528 = t2(0);
  Real copt529 = Power(copt528, 2);
  Real copt461 = copt302 * copt382;
  Real copt462 = copt165 * copt305;
  Real copt463 = -(copt310 * copt50);
  Real copt464 = copt401 * copt45;
  Real copt465 = copt461 + copt462 + copt463 + copt464;
  Real copt531 = xloc(15);
  Real copt536 = xloc(16);
  Real copt532 = -copt531;
  Real copt533 = copt165 + copt532;
  Real copt544 = xloc(17);
  Real copt486 = copt316 * copt365;
  Real copt487 = copt310 * copt319;
  Real copt488 = -(copt305 * copt323);
  Real copt489 = copt302 * copt387;
  Real copt490 = copt486 + copt487 + copt488 + copt489;
  Real copt546 = copt403 + copt544;
  Real copt559 = copt381 + copt45;
  Real copt560 = Power(copt559, 2);
  Real copt561 = copt302 + copt364;
  Real copt562 = Power(copt561, 2);
  Real copt563 = copt316 + copt403;
  Real copt564 = Power(copt563, 2);
  Real copt565 = copt560 + copt562 + copt564;
  Real copt566 = Sqrt(copt565);
  Real copt535 = copt310 * copt531;
  Real copt537 = -(copt165 * copt536);
  Real copt538 = copt364 + copt536;
  Real copt360 = -(copt165 * copt305);
  Real copt362 = copt165 + copt361;
  Real copt363 = copt302 * copt362;
  Real copt366 = copt365 * copt45;
  Real copt367 = copt310 * copt50;
  Real copt368 = copt360 + copt363 + copt366 + copt367;
  Real copt371 = copt370 + copt50;
  Real copt372 = copt302 * copt371;
  Real copt373 = copt305 * copt369;
  Real copt375 = -(copt374 * copt50);
  Real copt377 = copt374 + copt376;
  Real copt378 = copt377 * copt45;
  Real copt379 = copt372 + copt373 + copt375 + copt378;
  Real copt380 = copt368 * copt379;
  Real copt390 = -(copt319 * copt369);
  Real copt391 = copt361 + copt369;
  Real copt392 = copt316 * copt391;
  Real copt394 = -copt393;
  Real copt395 = copt319 + copt394;
  Real copt396 = copt395 * copt45;
  Real copt397 = copt393 * copt50;
  Real copt398 = copt390 + copt392 + copt396 + copt397;
  Real copt399 = copt389 * copt398;
  Real copt400 = -(copt310 * copt319);
  Real copt402 = copt316 * copt401;
  Real copt405 = copt302 * copt404;
  Real copt406 = copt305 * copt323;
  Real copt407 = copt400 + copt402 + copt405 + copt406;
  Real copt409 = copt305 + copt408;
  Real copt410 = copt316 * copt409;
  Real copt411 = copt319 * copt374;
  Real copt412 = -(copt305 * copt393);
  Real copt413 = copt386 + copt393;
  Real copt414 = copt302 * copt413;
  Real copt415 = copt410 + copt411 + copt412 + copt414;
  Real copt416 = copt407 * copt415;
  Real copt417 = copt380 + copt399 + copt416;
  Real copt428 = copt310 * copt319 * copt369;
  Real copt429 = -(copt305 * copt323 * copt369);
  Real copt430 = copt319 * copt374 * copt45;
  Real copt431 = -(copt165 * copt319 * copt374);
  Real copt432 = -(copt323 * copt374 * copt45);
  Real copt433 = copt323 * copt374 * copt50;
  Real copt434 = -(copt310 * copt369);
  Real copt435 = copt369 + copt381;
  Real copt436 = copt305 * copt435;
  Real copt437 = copt310 + copt408;
  Real copt438 = copt437 * copt50;
  Real copt439 = copt165 * copt374;
  Real copt440 = copt434 + copt436 + copt438 + copt439;
  Real copt441 = copt316 * copt440;
  Real copt442 = -(copt305 * copt393 * copt45);
  Real copt443 = copt165 * copt305 * copt393;
  Real copt444 = copt310 * copt393 * copt45;
  Real copt445 = -(copt310 * copt393 * copt50);
  Real copt446 = copt165 + copt370;
  Real copt447 = copt319 * copt446;
  Real copt448 = copt323 * copt369;
  Real copt449 = -(copt165 * copt393);
  Real copt450 = copt393 + copt403;
  Real copt451 = copt450 * copt50;
  Real copt452 = copt447 + copt448 + copt449 + copt451;
  Real copt453 = copt302 * copt452;
  Real copt454 = copt426 + copt427 + copt428 + copt429 + copt430 + copt431 +
                 copt432 + copt433 + copt441 + copt442 + copt443 + copt444 +
                 copt445 + copt453;
  Real copt455 = copt425 * copt454;
  Real copt456 = ArcTan(copt417, copt455);
  Real copt598 = t0(1);
  Real copt476 = copt465 * copt475;
  Real copt477 = -(copt323 * copt466);
  Real copt478 = copt319 * copt468;
  Real copt482 = copt481 * copt50;
  Real copt483 = copt165 * copt479;
  Real copt484 = copt477 + copt478 + copt482 + copt483;
  Real copt485 = copt389 * copt484;
  Real copt491 = -(copt323 * copt470);
  Real copt492 = copt364 + copt470;
  Real copt493 = copt319 * copt492;
  Real copt494 = copt305 * copt481;
  Real copt495 = copt310 * copt479;
  Real copt496 = copt491 + copt493 + copt494 + copt495;
  Real copt497 = copt490 * copt496;
  Real copt498 = copt476 + copt485 + copt497;
  Real copt504 = copt310 * copt319 * copt466;
  Real copt505 = -(copt305 * copt323 * copt466);
  Real copt506 = copt319 * copt45 * copt470;
  Real copt507 = -(copt165 * copt319 * copt470);
  Real copt508 = -(copt323 * copt45 * copt470);
  Real copt509 = copt323 * copt470 * copt50;
  Real copt510 = copt316 * copt475;
  Real copt511 = -(copt305 * copt45 * copt479);
  Real copt512 = copt165 * copt305 * copt479;
  Real copt513 = copt310 * copt45 * copt479;
  Real copt514 = -(copt310 * copt479 * copt50);
  Real copt515 = -copt466;
  Real copt516 = copt165 + copt515;
  Real copt517 = copt319 * copt516;
  Real copt518 = copt323 * copt466;
  Real copt519 = -(copt165 * copt479);
  Real copt520 = copt403 + copt479;
  Real copt521 = copt50 * copt520;
  Real copt522 = copt517 + copt518 + copt519 + copt521;
  Real copt523 = copt302 * copt522;
  Real copt524 = copt426 + copt427 + copt504 + copt505 + copt506 + copt507 +
                 copt508 + copt509 + copt510 + copt511 + copt512 + copt513 +
                 copt514 + copt523;
  Real copt525 = copt503 * copt524;
  Real copt526 = ArcTan(copt498, copt525);
  Real copt601 = t1(1);
  Real copt534 = copt302 * copt533;
  Real copt539 = copt45 * copt538;
  Real copt540 = copt534 + copt535 + copt537 + copt539;
  Real copt541 = copt465 * copt540;
  Real copt542 = copt316 * copt533;
  Real copt543 = copt323 * copt531;
  Real copt545 = -(copt165 * copt544);
  Real copt547 = copt45 * copt546;
  Real copt548 = copt542 + copt543 + copt545 + copt547;
  Real copt549 = copt389 * copt548;
  Real copt550 = -copt536;
  Real copt551 = copt310 + copt550;
  Real copt552 = copt316 * copt551;
  Real copt553 = copt323 * copt536;
  Real copt554 = -(copt310 * copt544);
  Real copt555 = copt302 * copt546;
  Real copt556 = copt552 + copt553 + copt554 + copt555;
  Real copt557 = copt490 * copt556;
  Real copt558 = copt541 + copt549 + copt557;
  Real copt567 = copt310 * copt319 * copt45;
  Real copt568 = -(copt305 * copt323 * copt45);
  Real copt569 = -(copt310 * copt319 * copt531);
  Real copt570 = copt305 * copt323 * copt531;
  Real copt571 = -(copt319 * copt45 * copt536);
  Real copt572 = copt165 * copt319 * copt536;
  Real copt573 = copt323 * copt45 * copt536;
  Real copt574 = -(copt323 * copt50 * copt536);
  Real copt575 = copt305 * copt533;
  Real copt576 = copt50 * copt538;
  Real copt577 = copt535 + copt537 + copt575 + copt576;
  Real copt578 = copt316 * copt577;
  Real copt579 = copt305 * copt45 * copt544;
  Real copt580 = -(copt165 * copt305 * copt544);
  Real copt581 = -(copt310 * copt45 * copt544);
  Real copt582 = copt310 * copt50 * copt544;
  Real copt583 = -(copt323 * copt531);
  Real copt584 = copt381 + copt531;
  Real copt585 = copt319 * copt584;
  Real copt586 = -copt544;
  Real copt587 = copt323 + copt586;
  Real copt588 = copt50 * copt587;
  Real copt589 = copt165 * copt544;
  Real copt590 = copt583 + copt585 + copt588 + copt589;
  Real copt591 = copt302 * copt590;
  Real copt592 = copt567 + copt568 + copt569 + copt570 + copt571 + copt572 +
                 copt573 + copt574 + copt578 + copt579 + copt580 + copt581 +
                 copt582 + copt591;
  Real copt593  = -(copt566 * copt592);
  Real copt594  = ArcTan(copt558, copt593);
  Real copt604  = t2(1);
  Real copt609  = Power(copt598, 2);
  Real copt612  = Power(copt601, 2);
  Real copt615  = Power(copt604, 2);
  Real copt620  = copt162 + copt41;
  Real copt639  = copt328 * copt329;
  Real copt640  = 1 / copt639;
  Real copt641  = copt349 * copt350;
  Real copt642  = 1 / copt641;
  Real copt336  = copt300 * copt335;
  Real copt340  = copt313 * copt339;
  Real copt344  = copt326 * copt343;
  Real copt345  = copt336 + copt340 + copt344;
  Real copt643  = copt331 + copt333;
  Real copt621  = copt41 * copt418;
  Real copt622  = copt162 * copt559;
  Real copt623  = copt621 + copt622;
  Real copt625  = copt41 * copt420;
  Real copt626  = copt162 * copt561;
  Real copt627  = copt625 + copt626;
  Real copt629  = copt41 * copt422;
  Real copt630  = copt162 * copt563;
  Real copt631  = copt629 + copt630;
  Real copt732  = -copt331;
  Real copt733  = -copt333;
  Real copt734  = copt732 + copt733;
  Real copt749  = Power(copt417, 2);
  Real copt750  = Power(copt454, 2);
  Real copt751  = copt424 * copt750;
  Real copt752  = copt749 + copt751;
  Real copt753  = 1 / copt752;
  Real copt759  = 1 / copt425;
  Real copt768  = Power(copt498, 2);
  Real copt769  = Power(copt524, 2);
  Real copt770  = copt502 * copt769;
  Real copt771  = copt768 + copt770;
  Real copt772  = 1 / copt771;
  Real copt785  = Power(copt592, 2);
  Real copt786  = copt565 * copt785;
  Real copt787  = Power(copt558, 2);
  Real copt788  = copt786 + copt787;
  Real copt789  = 1 / copt788;
  Real copt795  = 1 / copt566;
  Real copt900  = 1 / copt503;
  Real copt928  = copt323 * copt45;
  Real copt975  = -(copt310 * copt45);
  Real copt1031 = copt376 + copt470;
  Real copt1034 = copt386 + copt479;
  Real copt1074 = -(copt319 * copt45);
  Real copt1087 = copt50 + copt515;
  Real copt1054 = copt316 + copt586;
  Real copt1127 = copt305 * copt45;
  Real copt1105 = copt49 + copt531;
  Real copt1112 = copt319 * copt45;
  Real copt959  = -(copt323 * copt45);
  Real copt1166 = -(copt305 * copt45);
  Real copt1008 = copt310 * copt45;
  Real copt1188 = -(copt165 * copt319);
  Real copt1189 = copt316 * copt362;
  Real copt1190 = copt323 * copt50;
  Real copt1191 = copt1112 + copt1188 + copt1189 + copt1190 + copt959;
  Real copt1199 = copt1008 + copt1166 + copt461 + copt462 + copt463;
  Real copt744  = copt368 * copt377;
  Real copt745  = copt365 * copt379;
  Real copt746  = copt389 * copt395;
  Real copt747  = copt387 * copt398;
  Real copt748  = copt744 + copt745 + copt746 + copt747;
  Real copt754  = -(copt425 * copt454 * copt748 * copt753);
  Real copt755  = -(copt323 * copt374);
  Real copt756  = copt310 * copt393;
  Real copt757  = copt400 + copt406 + copt411 + copt412 + copt755 + copt756;
  Real copt758  = copt425 * copt757;
  Real copt760  = copt418 * copt454 * copt759;
  Real copt761  = copt758 + copt760;
  Real copt762  = copt417 * copt753 * copt761;
  Real copt763  = copt754 + copt762;
  Real copt765  = copt319 * copt470;
  Real copt766  = -(copt305 * copt479);
  Real copt767  = copt400 + copt406 + copt491 + copt495 + copt765 + copt766;
  Real copt773  = copt498 * copt503 * copt767 * copt772;
  Real copt774  = copt401 * copt475;
  Real copt775  = copt387 * copt484;
  Real copt776  = copt774 + copt775;
  Real copt777  = -(copt503 * copt524 * copt772 * copt776);
  Real copt778  = copt773 + copt777;
  Real copt780  = copt465 * copt538;
  Real copt781  = copt401 * copt540;
  Real copt782  = copt389 * copt546;
  Real copt783  = copt387 * copt548;
  Real copt784  = copt780 + copt781 + copt782 + copt783;
  Real copt790  = copt566 * copt592 * copt784 * copt789;
  Real copt791  = -(copt319 * copt536);
  Real copt792  = copt305 * copt544;
  Real copt793  = copt487 + copt488 + copt553 + copt554 + copt791 + copt792;
  Real copt794  = -(copt566 * copt793);
  Real copt796  = -(copt559 * copt592 * copt795);
  Real copt797  = copt794 + copt796;
  Real copt798  = copt558 * copt789 * copt797;
  Real copt799  = copt790 + copt798;
  Real copt803  = copt368 * copt371;
  Real copt804  = copt362 * copt379;
  Real copt805  = copt407 * copt413;
  Real copt806  = copt404 * copt415;
  Real copt807  = copt803 + copt804 + copt805 + copt806;
  Real copt808  = -(copt425 * copt454 * copt753 * copt807);
  Real copt809  = copt425 * copt452;
  Real copt810  = copt420 * copt454 * copt759;
  Real copt811  = copt809 + copt810;
  Real copt812  = copt417 * copt753 * copt811;
  Real copt813  = copt808 + copt812;
  Real copt815  = copt498 * copt503 * copt522 * copt772;
  Real copt816  = copt382 * copt475;
  Real copt817  = copt387 * copt496;
  Real copt818  = copt816 + copt817;
  Real copt819  = -(copt503 * copt524 * copt772 * copt818);
  Real copt820  = copt815 + copt819;
  Real copt822  = copt465 * copt533;
  Real copt823  = copt382 * copt540;
  Real copt824  = copt490 * copt546;
  Real copt825  = copt387 * copt556;
  Real copt826  = copt822 + copt823 + copt824 + copt825;
  Real copt827  = copt566 * copt592 * copt789 * copt826;
  Real copt828  = -(copt566 * copt590);
  Real copt829  = -(copt561 * copt592 * copt795);
  Real copt830  = copt828 + copt829;
  Real copt831  = copt558 * copt789 * copt830;
  Real copt832  = copt827 + copt831;
  Real copt836  = copt389 * copt391;
  Real copt837  = copt407 * copt409;
  Real copt838  = copt382 * copt398;
  Real copt839  = copt401 * copt415;
  Real copt840  = copt836 + copt837 + copt838 + copt839;
  Real copt841  = -(copt425 * copt454 * copt753 * copt840);
  Real copt842  = copt425 * copt440;
  Real copt843  = copt422 * copt454 * copt759;
  Real copt844  = copt842 + copt843;
  Real copt845  = copt417 * copt753 * copt844;
  Real copt846  = copt841 + copt845;
  Real copt848  = copt475 * copt498 * copt503 * copt772;
  Real copt849  = copt382 * copt484;
  Real copt850  = copt365 * copt496;
  Real copt851  = copt849 + copt850;
  Real copt852  = -(copt503 * copt524 * copt772 * copt851);
  Real copt853  = copt848 + copt852;
  Real copt855  = copt389 * copt533;
  Real copt856  = copt490 * copt551;
  Real copt857  = copt382 * copt548;
  Real copt858  = copt365 * copt556;
  Real copt859  = copt855 + copt856 + copt857 + copt858;
  Real copt860  = copt566 * copt592 * copt789 * copt859;
  Real copt861  = -(copt566 * copt577);
  Real copt862  = -(copt563 * copt592 * copt795);
  Real copt863  = copt861 + copt862;
  Real copt864  = copt558 * copt789 * copt863;
  Real copt865  = copt860 + copt864;
  Real copt869  = copt302 + copt408;
  Real copt870  = copt368 * copt869;
  Real copt871  = copt311 * copt379;
  Real copt872  = copt318 + copt393;
  Real copt873  = copt389 * copt872;
  Real copt874  = copt398 * copt563;
  Real copt875  = copt870 + copt871 + copt873 + copt874;
  Real copt876  = -(copt425 * copt454 * copt753 * copt875);
  Real copt877  = copt316 * copt437;
  Real copt878  = copt323 * copt374;
  Real copt879  = -(copt310 * copt393);
  Real copt880  = copt302 * copt450;
  Real copt881  = copt877 + copt878 + copt879 + copt880;
  Real copt882  = copt425 * copt881;
  Real copt883  = -(copt418 * copt454 * copt759);
  Real copt884  = copt882 + copt883;
  Real copt885  = copt417 * copt753 * copt884;
  Real copt886  = copt876 + copt885;
  Real copt888  = copt465 * copt472;
  Real copt889  = copt475 * copt561;
  Real copt890  = copt389 * copt481;
  Real copt891  = copt484 * copt563;
  Real copt892  = copt888 + copt889 + copt890 + copt891;
  Real copt893  = -(copt503 * copt524 * copt772 * copt892);
  Real copt894  = copt316 * copt472;
  Real copt895  = copt323 * copt470;
  Real copt896  = -(copt310 * copt479);
  Real copt897  = copt302 * copt520;
  Real copt898  = copt894 + copt895 + copt896 + copt897;
  Real copt899  = copt503 * copt898;
  Real copt901  = copt382 * copt524 * copt900;
  Real copt902  = copt899 + copt901;
  Real copt903  = copt498 * copt772 * copt902;
  Real copt904  = copt893 + copt903;
  Real copt906  = copt540 * copt561;
  Real copt907  = copt548 * copt563;
  Real copt908  = copt906 + copt907;
  Real copt909  = copt566 * copt592 * copt789 * copt908;
  Real copt910  = -(copt323 * copt536);
  Real copt911  = copt316 * copt538;
  Real copt912  = copt302 * copt587;
  Real copt913  = copt310 * copt544;
  Real copt914  = copt910 + copt911 + copt912 + copt913;
  Real copt915  = -(copt558 * copt566 * copt789 * copt914);
  Real copt916  = copt909 + copt915;
  Real copt920  = copt369 + copt49;
  Real copt921  = copt368 * copt920;
  Real copt922  = copt379 * copt559;
  Real copt923  = copt316 + copt394;
  Real copt924  = copt407 * copt923;
  Real copt925  = copt324 * copt415;
  Real copt926  = copt921 + copt922 + copt924 + copt925;
  Real copt927  = -(copt425 * copt454 * copt753 * copt926);
  Real copt929  = -(copt323 * copt369);
  Real copt930  = copt316 * copt435;
  Real copt931  = -(copt393 * copt45);
  Real copt932  = copt165 * copt393;
  Real copt933  = copt928 + copt929 + copt930 + copt931 + copt932;
  Real copt934  = copt425 * copt933;
  Real copt935  = -(copt420 * copt454 * copt759);
  Real copt936  = copt934 + copt935;
  Real copt937  = copt417 * copt753 * copt936;
  Real copt938  = copt927 + copt937;
  Real copt940  = copt465 * copt468;
  Real copt941  = copt232 * copt475;
  Real copt942  = copt481 * copt490;
  Real copt943  = copt496 * copt563;
  Real copt944  = copt940 + copt941 + copt942 + copt943;
  Real copt945  = -(copt503 * copt524 * copt772 * copt944);
  Real copt946  = copt316 * copt468;
  Real copt947  = -(copt45 * copt479);
  Real copt948  = copt477 + copt483 + copt928 + copt946 + copt947;
  Real copt949  = copt503 * copt948;
  Real copt950  = copt365 * copt524 * copt900;
  Real copt951  = copt949 + copt950;
  Real copt952  = copt498 * copt772 * copt951;
  Real copt953  = copt945 + copt952;
  Real copt955  = copt232 * copt540;
  Real copt956  = copt556 * copt563;
  Real copt957  = copt955 + copt956;
  Real copt958  = copt566 * copt592 * copt789 * copt957;
  Real copt960  = copt45 * copt544;
  Real copt961  = copt542 + copt543 + copt545 + copt959 + copt960;
  Real copt962  = -(copt558 * copt566 * copt789 * copt961);
  Real copt963  = copt958 + copt962;
  Real copt967  = copt370 + copt45;
  Real copt968  = copt389 * copt967;
  Real copt969  = copt304 + copt374;
  Real copt970  = copt407 * copt969;
  Real copt971  = copt232 * copt398;
  Real copt972  = copt415 * copt561;
  Real copt973  = copt968 + copt970 + copt971 + copt972;
  Real copt974  = -(copt425 * copt454 * copt753 * copt973);
  Real copt976  = copt302 * copt446;
  Real copt977  = copt310 * copt369;
  Real copt978  = copt374 * copt45;
  Real copt979  = -(copt165 * copt374);
  Real copt980  = copt975 + copt976 + copt977 + copt978 + copt979;
  Real copt981  = copt425 * copt980;
  Real copt982  = -(copt422 * copt454 * copt759);
  Real copt983  = copt981 + copt982;
  Real copt984  = copt417 * copt753 * copt983;
  Real copt985  = copt974 + copt984;
  Real copt987  = copt389 * copt468;
  Real copt988  = copt490 * copt492;
  Real copt989  = copt232 * copt484;
  Real copt990  = copt311 * copt496;
  Real copt991  = copt987 + copt988 + copt989 + copt990;
  Real copt992  = -(copt503 * copt524 * copt772 * copt991);
  Real copt993  = copt302 * copt516;
  Real copt994  = copt310 * copt466;
  Real copt995  = copt45 * copt470;
  Real copt996  = -(copt165 * copt470);
  Real copt997  = copt975 + copt993 + copt994 + copt995 + copt996;
  Real copt998  = copt503 * copt997;
  Real copt999  = copt404 * copt524 * copt900;
  Real copt1000 = copt998 + copt999;
  Real copt1001 = copt1000 * copt498 * copt772;
  Real copt1002 = copt1001 + copt992;
  Real copt1004 = copt232 * copt548;
  Real copt1005 = copt311 * copt556;
  Real copt1006 = copt1004 + copt1005;
  Real copt1007 = copt1006 * copt566 * copt592 * copt789;
  Real copt1009 = -(copt310 * copt531);
  Real copt1010 = copt302 * copt584;
  Real copt1011 = -(copt45 * copt536);
  Real copt1012 = copt165 * copt536;
  Real copt1013 = copt1008 + copt1009 + copt1010 + copt1011 + copt1012;
  Real copt1014 = -(copt1013 * copt558 * copt566 * copt789);
  Real copt1015 = copt1007 + copt1014;
  Real copt1019 = -(copt319 * copt374);
  Real copt1020 = copt316 * copt377;
  Real copt1021 = copt302 * copt395;
  Real copt1022 = copt305 * copt393;
  Real copt1023 = copt1019 + copt1020 + copt1021 + copt1022;
  Real copt1024 = copt1023 * copt417 * copt425 * copt753;
  Real copt1025 = copt379 * copt420;
  Real copt1026 = copt321 * copt398;
  Real copt1027 = copt1025 + copt1026;
  Real copt1028 = -(copt1027 * copt425 * copt454 * copt753);
  Real copt1029 = copt1024 + copt1028;
  Real copt1032 = copt1031 * copt465;
  Real copt1033 = copt307 * copt475;
  Real copt1035 = copt1034 * copt389;
  Real copt1036 = copt321 * copt484;
  Real copt1037 = copt1032 + copt1033 + copt1035 + copt1036;
  Real copt1038 = -(copt1037 * copt503 * copt524 * copt772);
  Real copt1039 = -(copt319 * copt470);
  Real copt1040 = copt1031 * copt316;
  Real copt1041 = copt319 + copt480;
  Real copt1042 = copt1041 * copt302;
  Real copt1043 = copt305 * copt479;
  Real copt1044 = copt1039 + copt1040 + copt1042 + copt1043;
  Real copt1045 = copt1044 * copt503;
  Real copt1046 = -(copt382 * copt524 * copt900);
  Real copt1047 = copt1045 + copt1046;
  Real copt1048 = copt1047 * copt498 * copt772;
  Real copt1049 = copt1038 + copt1048;
  Real copt1051 = copt302 + copt550;
  Real copt1052 = copt1051 * copt465;
  Real copt1053 = copt307 * copt540;
  Real copt1055 = copt1054 * copt389;
  Real copt1056 = copt321 * copt548;
  Real copt1057 = copt1052 + copt1053 + copt1055 + copt1056;
  Real copt1058 = copt1057 * copt566 * copt592 * copt789;
  Real copt1059 = copt305 + copt550;
  Real copt1060 = copt1059 * copt316;
  Real copt1061 = copt319 * copt536;
  Real copt1062 = -(copt305 * copt544);
  Real copt1063 = copt386 + copt544;
  Real copt1064 = copt1063 * copt302;
  Real copt1065 = copt1060 + copt1061 + copt1062 + copt1064;
  Real copt1066 = -(copt1065 * copt566);
  Real copt1067 = copt559 * copt592 * copt795;
  Real copt1068 = copt1066 + copt1067;
  Real copt1069 = copt1068 * copt558 * copt789;
  Real copt1070 = copt1058 + copt1069;
  Real copt1075 = copt316 * copt371;
  Real copt1076 = copt319 * copt369;
  Real copt1077 = copt393 * copt45;
  Real copt1078 = -(copt393 * copt50);
  Real copt1079 = copt1074 + copt1075 + copt1076 + copt1077 + copt1078;
  Real copt1080 = copt1079 * copt417 * copt425 * copt753;
  Real copt1081 = copt379 * copt57;
  Real copt1082 = copt415 * copt422;
  Real copt1083 = copt1081 + copt1082;
  Real copt1084 = -(copt1083 * copt425 * copt454 * copt753);
  Real copt1085 = copt1080 + copt1084;
  Real copt1088 = copt1087 * copt465;
  Real copt1089 = copt418 * copt475;
  Real copt1090 = copt1034 * copt490;
  Real copt1091 = copt321 * copt496;
  Real copt1092 = copt1088 + copt1089 + copt1090 + copt1091;
  Real copt1093 = -(copt1092 * copt503 * copt524 * copt772);
  Real copt1094 = copt1087 * copt316;
  Real copt1095 = copt319 * copt466;
  Real copt1096 = copt45 * copt479;
  Real copt1097 = -(copt479 * copt50);
  Real copt1098 = copt1074 + copt1094 + copt1095 + copt1096 + copt1097;
  Real copt1099 = copt1098 * copt503;
  Real copt1100 = -(copt365 * copt524 * copt900);
  Real copt1101 = copt1099 + copt1100;
  Real copt1102 = copt1101 * copt498 * copt772;
  Real copt1103 = copt1093 + copt1102;
  Real copt1106 = copt1105 * copt465;
  Real copt1107 = copt418 * copt540;
  Real copt1108 = copt1054 * copt490;
  Real copt1109 = copt321 * copt556;
  Real copt1110 = copt1106 + copt1107 + copt1108 + copt1109;
  Real copt1111 = copt1110 * copt566 * copt592 * copt789;
  Real copt1113 = -(copt319 * copt531);
  Real copt1114 = copt361 + copt531;
  Real copt1115 = copt1114 * copt316;
  Real copt1116 = -(copt45 * copt544);
  Real copt1117 = copt50 * copt544;
  Real copt1118 = copt1112 + copt1113 + copt1115 + copt1116 + copt1117;
  Real copt1119 = -(copt1118 * copt566);
  Real copt1120 = copt561 * copt592 * copt795;
  Real copt1121 = copt1119 + copt1120;
  Real copt1122 = copt1121 * copt558 * copt789;
  Real copt1123 = copt1111 + copt1122;
  Real copt1128 = -(copt305 * copt369);
  Real copt1129 = copt302 * copt391;
  Real copt1130 = -(copt374 * copt45);
  Real copt1131 = copt374 * copt50;
  Real copt1132 = copt1127 + copt1128 + copt1129 + copt1130 + copt1131;
  Real copt1133 = copt1132 * copt417 * copt425 * copt753;
  Real copt1134 = copt398 * copt418;
  Real copt1135 = copt307 * copt415;
  Real copt1136 = copt1134 + copt1135;
  Real copt1137 = -(copt1136 * copt425 * copt454 * copt753);
  Real copt1138 = copt1133 + copt1137;
  Real copt1140 = copt1087 * copt389;
  Real copt1141 = copt305 + copt471;
  Real copt1142 = copt1141 * copt490;
  Real copt1143 = copt418 * copt484;
  Real copt1144 = copt420 * copt496;
  Real copt1145 = copt1140 + copt1142 + copt1143 + copt1144;
  Real copt1146 = -(copt1145 * copt503 * copt524 * copt772);
  Real copt1147 = -(copt305 * copt466);
  Real copt1148 = copt361 + copt466;
  Real copt1149 = copt1148 * copt302;
  Real copt1150 = -(copt45 * copt470);
  Real copt1151 = copt470 * copt50;
  Real copt1152 = copt1127 + copt1147 + copt1149 + copt1150 + copt1151;
  Real copt1153 = copt1152 * copt503;
  Real copt1154 = -(copt404 * copt524 * copt900);
  Real copt1155 = copt1153 + copt1154;
  Real copt1156 = copt1155 * copt498 * copt772;
  Real copt1157 = copt1146 + copt1156;
  Real copt1159 = copt1105 * copt389;
  Real copt1160 = copt304 + copt536;
  Real copt1161 = copt1160 * copt490;
  Real copt1162 = copt418 * copt548;
  Real copt1163 = copt420 * copt556;
  Real copt1164 = copt1159 + copt1161 + copt1162 + copt1163;
  Real copt1165 = copt1164 * copt566 * copt592 * copt789;
  Real copt1167 = copt50 + copt532;
  Real copt1168 = copt1167 * copt302;
  Real copt1169 = copt305 * copt531;
  Real copt1170 = copt45 * copt536;
  Real copt1171 = -(copt50 * copt536);
  Real copt1172 = copt1166 + copt1168 + copt1169 + copt1170 + copt1171;
  Real copt1173 = -(copt1172 * copt566);
  Real copt1174 = copt563 * copt592 * copt795;
  Real copt1175 = copt1173 + copt1174;
  Real copt1176 = copt1175 * copt558 * copt789;
  Real copt1177 = copt1165 + copt1176;
  Real copt1181 = copt417 * copt425 * copt490 * copt753;
  Real copt1182 = copt307 * copt368;
  Real copt1183 = copt389 * copt422;
  Real copt1184 = copt1182 + copt1183;
  Real copt1185 = -(copt1184 * copt425 * copt454 * copt753);
  Real copt1186 = copt1181 + copt1185;
  Real copt1192 = copt1191 * copt417 * copt425 * copt753;
  Real copt1193 = copt368 * copt418;
  Real copt1194 = copt321 * copt407;
  Real copt1195 = copt1193 + copt1194;
  Real copt1196 = -(copt1195 * copt425 * copt454 * copt753);
  Real copt1197 = copt1192 + copt1196;
  Real copt1200 = copt1199 * copt417 * copt425 * copt753;
  Real copt1201 = copt407 * copt420;
  Real copt1202 = copt389 * copt57;
  Real copt1203 = copt1201 + copt1202;
  Real copt1204 = -(copt1203 * copt425 * copt454 * copt753);
  Real copt1205 = copt1200 + copt1204;
  Real copt1207 = copt490 * copt498 * copt503 * copt772;
  Real copt1208 = copt365 * copt465;
  Real copt1209 = copt389 * copt404;
  Real copt1210 = copt1208 + copt1209;
  Real copt1211 = -(copt1210 * copt503 * copt524 * copt772);
  Real copt1212 = copt1207 + copt1211;
  Real copt1214 = copt1191 * copt498 * copt503 * copt772;
  Real copt1215 = copt362 * copt465;
  Real copt1216 = copt404 * copt490;
  Real copt1217 = copt1215 + copt1216;
  Real copt1218 = -(copt1217 * copt503 * copt524 * copt772);
  Real copt1219 = copt1214 + copt1218;
  Real copt1221 = copt1199 * copt498 * copt503 * copt772;
  Real copt1222 = copt362 * copt389;
  Real copt1223 = copt401 * copt490;
  Real copt1224 = copt1222 + copt1223;
  Real copt1225 = -(copt1224 * copt503 * copt524 * copt772);
  Real copt1226 = copt1221 + copt1225;
  Real copt1228 = copt311 * copt465;
  Real copt1229 = copt324 * copt389;
  Real copt1230 = copt1228 + copt1229;
  Real copt1231 = copt1230 * copt566 * copt592 * copt789;
  Real copt1232 = -(copt407 * copt558 * copt566 * copt789);
  Real copt1233 = copt1231 + copt1232;
  Real copt1235 = copt465 * copt559;
  Real copt1236 = copt324 * copt490;
  Real copt1237 = copt1235 + copt1236;
  Real copt1238 = copt1237 * copt566 * copt592 * copt789;
  Real copt1239 = copt1074 + copt383 + copt384 + copt385 + copt928;
  Real copt1240 = -(copt1239 * copt558 * copt566 * copt789);
  Real copt1241 = copt1238 + copt1240;
  Real copt1243 = copt389 * copt559;
  Real copt1244 = copt490 * copt561;
  Real copt1245 = copt1243 + copt1244;
  Real copt1246 = copt1245 * copt566 * copt592 * copt789;
  Real copt1247 = copt1127 + copt360 + copt363 + copt367 + copt975;
  Real copt1248 = -(copt1247 * copt558 * copt566 * copt789);
  Real copt1249 = copt1246 + copt1248;
  out1(0)       = copt329;
  out1(1)       = copt330 * copt345 * copt351;
  out1(2)       = copt350;
  out1(3) =
      -(copt353 * copt354 * copt355 * copt356 *
        (copt529 * copt594 * l0 * l1 + copt459 * copt526 * l0 * l2 +
         copt358 * copt456 * l1 * l2 + copt358 * l1 * l2 * thetarest0 +
         copt459 * l0 * l2 * thetarest1 + copt529 * l0 * l1 * thetarest2)) /
      2.;
  out1(4) = -(copt353 * copt354 * copt355 * copt356 *
              (copt528 * copt594 * copt604 * l0 * l1 +
               copt458 * copt526 * copt601 * l0 * l2 +
               copt357 * copt456 * copt598 * l1 * l2 +
               copt357 * copt598 * l1 * l2 * thetarest0 +
               copt458 * copt601 * l0 * l2 * thetarest1 +
               copt528 * copt604 * l0 * l1 * thetarest2)) /
            2.;
  out1(5) =
      -(copt353 * copt354 * copt355 * copt356 *
        (copt594 * copt615 * l0 * l1 + copt526 * copt612 * l0 * l2 +
         copt456 * copt609 * l1 * l2 + copt609 * l1 * l2 * thetarest0 +
         copt612 * l0 * l2 * thetarest1 + copt615 * l0 * l1 * thetarest2)) /
      2.;
  out2(0, 0)  = copt330 * copt620 * copt623;
  out2(0, 1)  = copt330 * copt620 * copt627;
  out2(0, 2)  = copt330 * copt620 * copt631;
  out2(0, 3)  = copt300 * copt330 * copt41;
  out2(0, 4)  = copt313 * copt330 * copt41;
  out2(0, 5)  = copt326 * copt330 * copt41;
  out2(0, 6)  = copt162 * copt300 * copt330;
  out2(0, 7)  = copt162 * copt313 * copt330;
  out2(0, 8)  = copt162 * copt326 * copt330;
  out2(0, 9)  = 0;
  out2(0, 10) = 0;
  out2(0, 11) = 0;
  out2(0, 12) = 0;
  out2(0, 13) = 0;
  out2(0, 14) = 0;
  out2(0, 15) = 0;
  out2(0, 16) = 0;
  out2(0, 17) = 0;
  out2(1, 0)  = copt640 * copt642 *
               (copt300 * copt345 * copt349 * copt620 +
                copt328 * copt335 * copt345 * copt643 +
                copt328 * copt349 *
                    ((copt331 * copt418 + copt333 * copt559) * copt620 +
                     copt623 * copt643));
  out2(1, 1) = copt640 * copt642 *
               (copt313 * copt345 * copt349 * copt620 +
                copt328 * copt339 * copt345 * copt643 +
                copt328 * copt349 *
                    ((copt331 * copt420 + copt333 * copt561) * copt620 +
                     copt627 * copt643));
  out2(1, 2) = copt640 * copt642 *
               (copt326 * copt345 * copt349 * copt620 +
                copt328 * copt343 * copt345 * copt643 +
                copt328 * copt349 *
                    ((copt331 * copt422 + copt333 * copt563) * copt620 +
                     copt631 * copt643));
  out2(1, 3) = (-(copt328 * copt331 * copt335 * copt345) -
                copt300 * copt345 * copt349 * copt41 +
                copt328 * copt349 *
                    (copt162 * copt232 * copt331 +
                     copt41 * (copt334 - 2 * copt331 * copt418))) *
               copt640 * copt642;
  out2(1, 4) = (-(copt328 * copt331 * copt339 * copt345) -
                copt313 * copt345 * copt349 * copt41 +
                copt328 * copt349 *
                    (copt162 * copt311 * copt331 +
                     copt41 * (copt338 - 2 * copt331 * copt420))) *
               copt640 * copt642;
  out2(1, 5) = (-(copt328 * copt331 * copt343 * copt345) -
                copt326 * copt345 * copt349 * copt41 +
                copt328 * copt349 *
                    (copt162 * copt324 * copt331 +
                     copt41 * (copt342 - 2 * copt331 * copt422))) *
               copt640 * copt642;
  out2(1, 6) = (-(copt328 * copt333 * copt335 * copt345) -
                copt162 * copt300 * copt345 * copt349 +
                copt328 * copt349 *
                    (copt162 * copt331 * copt57 +
                     copt333 * (2 * copt162 * copt232 + copt64))) *
               copt640 * copt642;
  out2(1, 7) = (-(copt328 * copt333 * copt339 * copt345) +
                copt328 * (copt313 * copt333 + copt162 * copt339) * copt349 -
                copt162 * copt313 * copt345 * copt349) *
               copt640 * copt642;
  out2(1, 8) = copt330 * (copt326 * copt333 + copt162 * copt343) * copt351 -
               copt162 * copt326 * copt345 * copt351 * copt640 -
               copt330 * copt333 * copt343 * copt345 * copt642;
  out2(1, 9)  = 0;
  out2(1, 10) = 0;
  out2(1, 11) = 0;
  out2(1, 12) = 0;
  out2(1, 13) = 0;
  out2(1, 14) = 0;
  out2(1, 15) = 0;
  out2(1, 16) = 0;
  out2(1, 17) = 0;
  out2(2, 0)  = copt335 * copt351 * copt734;
  out2(2, 1)  = copt339 * copt351 * copt734;
  out2(2, 2)  = copt343 * copt351 * copt734;
  out2(2, 3)  = copt331 * copt335 * copt351;
  out2(2, 4)  = copt331 * copt339 * copt351;
  out2(2, 5)  = copt331 * copt343 * copt351;
  out2(2, 6)  = copt333 * copt335 * copt351;
  out2(2, 7)  = copt333 * copt339 * copt351;
  out2(2, 8)  = copt333 * copt343 * copt351;
  out2(2, 9)  = 0;
  out2(2, 10) = 0;
  out2(2, 11) = 0;
  out2(2, 12) = 0;
  out2(2, 13) = 0;
  out2(2, 14) = 0;
  out2(2, 15) = 0;
  out2(2, 16) = 0;
  out2(2, 17) = 0;
  out2(3, 0)  = -(copt353 * copt354 * copt355 * copt356 *
                 (copt529 * copt799 * l0 * l1 + copt459 * copt778 * l0 * l2 +
                  copt358 * copt763 * l1 * l2)) /
               2.;
  out2(3, 1) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt529 * copt832 * l0 * l1 + copt459 * copt820 * l0 * l2 +
                  copt358 * copt813 * l1 * l2)) /
               2.;
  out2(3, 2) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt529 * copt865 * l0 * l1 + copt459 * copt853 * l0 * l2 +
                  copt358 * copt846 * l1 * l2)) /
               2.;
  out2(3, 3) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt529 * copt916 * l0 * l1 + copt459 * copt904 * l0 * l2 +
                  copt358 * copt886 * l1 * l2)) /
               2.;
  out2(3, 4) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt529 * copt963 * l0 * l1 + copt459 * copt953 * l0 * l2 +
                  copt358 * copt938 * l1 * l2)) /
               2.;
  out2(3, 5) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1015 * copt529 * l0 * l1 + copt1002 * copt459 * l0 * l2 +
                  copt358 * copt985 * l1 * l2)) /
               2.;
  out2(3, 6) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1070 * copt529 * l0 * l1 + copt1049 * copt459 * l0 * l2 +
                  copt1029 * copt358 * l1 * l2)) /
               2.;
  out2(3, 7) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1123 * copt529 * l0 * l1 + copt1103 * copt459 * l0 * l2 +
                  copt1085 * copt358 * l1 * l2)) /
               2.;
  out2(3, 8) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1177 * copt529 * l0 * l1 + copt1157 * copt459 * l0 * l2 +
                  copt1138 * copt358 * l1 * l2)) /
               2.;
  out2(3, 9)  = -(copt1186 * copt353 * copt354 * copt358) / 2.;
  out2(3, 10) = -(copt1197 * copt353 * copt354 * copt358) / 2.;
  out2(3, 11) = -(copt1205 * copt353 * copt354 * copt358) / 2.;
  out2(3, 12) = -(copt1212 * copt353 * copt355 * copt459) / 2.;
  out2(3, 13) = -(copt1219 * copt353 * copt355 * copt459) / 2.;
  out2(3, 14) = -(copt1226 * copt353 * copt355 * copt459) / 2.;
  out2(3, 15) = -(copt1233 * copt353 * copt356 * copt529) / 2.;
  out2(3, 16) = -(copt1241 * copt353 * copt356 * copt529) / 2.;
  out2(3, 17) = -(copt1249 * copt353 * copt356 * copt529) / 2.;
  out2(4, 0)  = -(copt353 * copt354 * copt355 * copt356 *
                 (copt528 * copt604 * copt799 * l0 * l1 +
                  copt458 * copt601 * copt778 * l0 * l2 +
                  copt357 * copt598 * copt763 * l1 * l2)) /
               2.;
  out2(4, 1) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt528 * copt604 * copt832 * l0 * l1 +
                  copt458 * copt601 * copt820 * l0 * l2 +
                  copt357 * copt598 * copt813 * l1 * l2)) /
               2.;
  out2(4, 2) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt528 * copt604 * copt865 * l0 * l1 +
                  copt458 * copt601 * copt853 * l0 * l2 +
                  copt357 * copt598 * copt846 * l1 * l2)) /
               2.;
  out2(4, 3) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt528 * copt604 * copt916 * l0 * l1 +
                  copt458 * copt601 * copt904 * l0 * l2 +
                  copt357 * copt598 * copt886 * l1 * l2)) /
               2.;
  out2(4, 4) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt528 * copt604 * copt963 * l0 * l1 +
                  copt458 * copt601 * copt953 * l0 * l2 +
                  copt357 * copt598 * copt938 * l1 * l2)) /
               2.;
  out2(4, 5) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1015 * copt528 * copt604 * l0 * l1 +
                  copt1002 * copt458 * copt601 * l0 * l2 +
                  copt357 * copt598 * copt985 * l1 * l2)) /
               2.;
  out2(4, 6) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1070 * copt528 * copt604 * l0 * l1 +
                  copt1049 * copt458 * copt601 * l0 * l2 +
                  copt1029 * copt357 * copt598 * l1 * l2)) /
               2.;
  out2(4, 7) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1123 * copt528 * copt604 * l0 * l1 +
                  copt1103 * copt458 * copt601 * l0 * l2 +
                  copt1085 * copt357 * copt598 * l1 * l2)) /
               2.;
  out2(4, 8) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1177 * copt528 * copt604 * l0 * l1 +
                  copt1157 * copt458 * copt601 * l0 * l2 +
                  copt1138 * copt357 * copt598 * l1 * l2)) /
               2.;
  out2(4, 9)  = -(copt1186 * copt353 * copt354 * copt357 * copt598) / 2.;
  out2(4, 10) = -(copt1197 * copt353 * copt354 * copt357 * copt598) / 2.;
  out2(4, 11) = -(copt1205 * copt353 * copt354 * copt357 * copt598) / 2.;
  out2(4, 12) = -(copt1212 * copt353 * copt355 * copt458 * copt601) / 2.;
  out2(4, 13) = -(copt1219 * copt353 * copt355 * copt458 * copt601) / 2.;
  out2(4, 14) = -(copt1226 * copt353 * copt355 * copt458 * copt601) / 2.;
  out2(4, 15) = -(copt1233 * copt353 * copt356 * copt528 * copt604) / 2.;
  out2(4, 16) = -(copt1241 * copt353 * copt356 * copt528 * copt604) / 2.;
  out2(4, 17) = -(copt1249 * copt353 * copt356 * copt528 * copt604) / 2.;
  out2(5, 0)  = -(copt353 * copt354 * copt355 * copt356 *
                 (copt615 * copt799 * l0 * l1 + copt612 * copt778 * l0 * l2 +
                  copt609 * copt763 * l1 * l2)) /
               2.;
  out2(5, 1) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt615 * copt832 * l0 * l1 + copt612 * copt820 * l0 * l2 +
                  copt609 * copt813 * l1 * l2)) /
               2.;
  out2(5, 2) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt615 * copt865 * l0 * l1 + copt612 * copt853 * l0 * l2 +
                  copt609 * copt846 * l1 * l2)) /
               2.;
  out2(5, 3) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt615 * copt916 * l0 * l1 + copt612 * copt904 * l0 * l2 +
                  copt609 * copt886 * l1 * l2)) /
               2.;
  out2(5, 4) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt615 * copt963 * l0 * l1 + copt612 * copt953 * l0 * l2 +
                  copt609 * copt938 * l1 * l2)) /
               2.;
  out2(5, 5) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1015 * copt615 * l0 * l1 + copt1002 * copt612 * l0 * l2 +
                  copt609 * copt985 * l1 * l2)) /
               2.;
  out2(5, 6) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1070 * copt615 * l0 * l1 + copt1049 * copt612 * l0 * l2 +
                  copt1029 * copt609 * l1 * l2)) /
               2.;
  out2(5, 7) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1123 * copt615 * l0 * l1 + copt1103 * copt612 * l0 * l2 +
                  copt1085 * copt609 * l1 * l2)) /
               2.;
  out2(5, 8) = -(copt353 * copt354 * copt355 * copt356 *
                 (copt1177 * copt615 * l0 * l1 + copt1157 * copt612 * l0 * l2 +
                  copt1138 * copt609 * l1 * l2)) /
               2.;
  out2(5, 9)  = -(copt1186 * copt353 * copt354 * copt609) / 2.;
  out2(5, 10) = -(copt1197 * copt353 * copt354 * copt609) / 2.;
  out2(5, 11) = -(copt1205 * copt353 * copt354 * copt609) / 2.;
  out2(5, 12) = -(copt1212 * copt353 * copt355 * copt612) / 2.;
  out2(5, 13) = -(copt1219 * copt353 * copt355 * copt612) / 2.;
  out2(5, 14) = -(copt1226 * copt353 * copt355 * copt612) / 2.;
  out2(5, 15) = -(copt1233 * copt353 * copt356 * copt615) / 2.;
  out2(5, 16) = -(copt1241 * copt353 * copt356 * copt615) / 2.;
  out2(5, 17) = -(copt1249 * copt353 * copt356 * copt615) / 2.;
  return std::make_tuple(grad, val);
}

#endif  // hylc_strain_II
