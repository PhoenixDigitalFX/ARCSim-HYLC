#include "strain.hpp"
#ifdef hylc_strain_II

using namespace hylc;
using namespace hylc::mathematica;

Vec6 hylc::mathematica::strain(const Vec18 &xloc, const Mat2x2 &invDm,
                               const Real &A, const Real &thetarest0,
                               const Real &thetarest1, const Real &thetarest2,
                               const Real &l0, const Real &l1, const Real &l2,
                               const Vec2 &t0, const Vec2 &t1, const Vec2 &t2) {
  // define output
  Vec6 out;

  Real copt1   = invDm(0, 0);
  Real copt2   = xloc(0);
  Real copt3   = -copt2;
  Real copt4   = xloc(3);
  Real copt5   = copt3 + copt4;
  Real copt6   = copt1 * copt5;
  Real copt7   = invDm(1, 0);
  Real copt8   = xloc(6);
  Real copt9   = copt3 + copt8;
  Real copt10  = copt7 * copt9;
  Real copt11  = copt10 + copt6;
  Real copt12  = Power(copt11, 2);
  Real copt13  = xloc(1);
  Real copt14  = -copt13;
  Real copt15  = xloc(4);
  Real copt16  = copt14 + copt15;
  Real copt17  = copt1 * copt16;
  Real copt18  = xloc(7);
  Real copt19  = copt14 + copt18;
  Real copt20  = copt19 * copt7;
  Real copt21  = copt17 + copt20;
  Real copt22  = Power(copt21, 2);
  Real copt23  = xloc(2);
  Real copt24  = -copt23;
  Real copt25  = xloc(5);
  Real copt26  = copt24 + copt25;
  Real copt27  = copt1 * copt26;
  Real copt28  = xloc(8);
  Real copt29  = copt24 + copt28;
  Real copt30  = copt29 * copt7;
  Real copt31  = copt27 + copt30;
  Real copt32  = Power(copt31, 2);
  Real copt33  = copt12 + copt22 + copt32;
  Real copt34  = Sqrt(copt33);
  Real copt35  = 1 / copt34;
  Real copt36  = invDm(0, 1);
  Real copt38  = invDm(1, 1);
  Real copt37  = copt36 * copt5;
  Real copt39  = copt38 * copt9;
  Real copt40  = copt37 + copt39;
  Real copt51  = Power(copt40, 2);
  Real copt42  = copt16 * copt36;
  Real copt43  = copt19 * copt38;
  Real copt44  = copt42 + copt43;
  Real copt52  = Power(copt44, 2);
  Real copt46  = copt26 * copt36;
  Real copt47  = copt29 * copt38;
  Real copt48  = copt46 + copt47;
  Real copt53  = Power(copt48, 2);
  Real copt54  = copt51 + copt52 + copt53;
  Real copt55  = Sqrt(copt54);
  Real copt56  = 1 / copt55;
  Real copt58  = 1 / A;
  Real copt59  = 1 / l0;
  Real copt60  = 1 / l1;
  Real copt61  = 1 / l2;
  Real copt62  = -copt4;
  Real copt63  = copt2 + copt62;
  Real copt64  = Power(copt63, 2);
  Real copt65  = -copt15;
  Real copt66  = copt13 + copt65;
  Real copt67  = Power(copt66, 2);
  Real copt68  = -copt25;
  Real copt69  = copt23 + copt68;
  Real copt70  = Power(copt69, 2);
  Real copt71  = copt64 + copt67 + copt70;
  Real copt72  = Sqrt(copt71);
  Real copt73  = Power(copt2, 2);
  Real copt74  = Power(copt15, 2);
  Real copt76  = Power(copt25, 2);
  Real copt84  = xloc(9);
  Real copt93  = xloc(10);
  Real copt100 = Power(copt4, 2);
  Real copt104 = Power(copt23, 2);
  Real copt113 = xloc(11);
  Real copt123 = Power(copt13, 2);
  Real copt105 = copt8 * copt84;
  Real copt106 = copt8 + copt84;
  Real copt107 = -(copt106 * copt4);
  Real copt125 = copt113 + copt28;
  Real copt139 = copt106 * copt25;
  Real copt109 = copt18 + copt93;
  Real copt145 = -2 * copt25;
  Real copt146 = copt113 + copt145 + copt28;
  Real copt163 = -(copt18 * copt84);
  Real copt138 = copt28 * copt84;
  Real copt144 = -(copt113 * copt8);
  Real copt196 = -copt8;
  Real copt225 = copt196 + copt4;
  Real copt226 = Power(copt225, 2);
  Real copt227 = -copt18;
  Real copt228 = copt15 + copt227;
  Real copt229 = Power(copt228, 2);
  Real copt212 = -copt28;
  Real copt230 = copt212 + copt25;
  Real copt231 = Power(copt230, 2);
  Real copt232 = copt226 + copt229 + copt231;
  Real copt233 = Sqrt(copt232);
  Real copt235 = Power(copt8, 2);
  Real copt241 = Power(copt18, 2);
  Real copt251 = Power(copt28, 2);
  Real copt254 = xloc(12);
  Real copt267 = xloc(13);
  Real copt280 = xloc(14);
  Real copt298 = -copt254;
  Real copt299 = copt298 + copt8;
  Real copt307 = -copt267;
  Real copt308 = copt18 + copt307;
  Real copt315 = -copt280;
  Real copt316 = copt28 + copt315;
  Real copt188 = -(copt18 * copt2 * copt25);
  Real copt189 = copt15 * copt2 * copt28;
  Real copt333 = copt308 * copt4;
  Real copt385 = copt196 + copt2;
  Real copt386 = Power(copt385, 2);
  Real copt387 = copt13 + copt227;
  Real copt388 = Power(copt387, 2);
  Real copt389 = copt212 + copt23;
  Real copt390 = Power(copt389, 2);
  Real copt391 = copt386 + copt388 + copt390;
  Real copt392 = Sqrt(copt391);
  Real copt401 = xloc(15);
  Real copt415 = xloc(16);
  Real copt410 = -copt235;
  Real copt411 = -copt401;
  Real copt412 = copt411 + copt8;
  Real copt413 = copt4 * copt412;
  Real copt414 = copt401 * copt8;
  Real copt430 = xloc(17);
  Real copt416 = -copt415;
  Real copt417 = copt18 + copt416;
  Real copt431 = -copt430;
  Real copt432 = copt28 + copt431;
  Real copt462 = copt4 * copt417;
  Real copt463 = copt415 * copt8;
  Real copt487 = copt25 * copt412;
  Real copt488 = copt28 * copt401;
  Real copt41  = copt11 * copt40;
  Real copt45  = copt21 * copt44;
  Real copt49  = copt31 * copt48;
  Real copt50  = copt41 + copt45 + copt49;
  Real copt75  = copt73 * copt74;
  Real copt77  = copt73 * copt76;
  Real copt78  = -(copt2 * copt74 * copt8);
  Real copt79  = -(copt2 * copt76 * copt8);
  Real copt80  = -(copt15 * copt18 * copt73);
  Real copt81  = copt15 * copt18 * copt2 * copt4;
  Real copt82  = -(copt25 * copt28 * copt73);
  Real copt83  = copt2 * copt25 * copt28 * copt4;
  Real copt85  = -(copt2 * copt74 * copt84);
  Real copt86  = -(copt2 * copt76 * copt84);
  Real copt87  = copt74 * copt8 * copt84;
  Real copt88  = copt76 * copt8 * copt84;
  Real copt89  = copt15 * copt18 * copt2 * copt84;
  Real copt90  = -(copt15 * copt18 * copt4 * copt84);
  Real copt91  = copt2 * copt25 * copt28 * copt84;
  Real copt92  = -(copt25 * copt28 * copt4 * copt84);
  Real copt94  = -(copt15 * copt73 * copt93);
  Real copt95  = copt15 * copt2 * copt4 * copt93;
  Real copt96  = copt15 * copt2 * copt8 * copt93;
  Real copt97  = -(copt15 * copt4 * copt8 * copt93);
  Real copt98  = copt18 * copt73 * copt93;
  Real copt99  = -2 * copt18 * copt2 * copt4 * copt93;
  Real copt101 = copt100 * copt18 * copt93;
  Real copt102 = copt18 * copt76 * copt93;
  Real copt103 = -(copt15 * copt25 * copt28 * copt93);
  Real copt108 = copt18 * copt93;
  Real copt110 = -(copt109 * copt15);
  Real copt111 = copt100 + copt105 + copt107 + copt108 + copt110 + copt74;
  Real copt112 = copt104 * copt111;
  Real copt114 = -(copt113 * copt25 * copt73);
  Real copt115 = copt113 * copt2 * copt25 * copt4;
  Real copt116 = copt113 * copt2 * copt25 * copt8;
  Real copt117 = -(copt113 * copt25 * copt4 * copt8);
  Real copt118 = -(copt113 * copt15 * copt18 * copt25);
  Real copt119 = copt113 * copt28 * copt73;
  Real copt120 = -2 * copt113 * copt2 * copt28 * copt4;
  Real copt121 = copt100 * copt113 * copt28;
  Real copt122 = copt113 * copt28 * copt74;
  Real copt124 = copt113 * copt28;
  Real copt126 = -(copt125 * copt25);
  Real copt127 = copt100 + copt105 + copt107 + copt124 + copt126 + copt76;
  Real copt128 = copt123 * copt127;
  Real copt129 = copt15 * copt18 * copt25;
  Real copt130 = -(copt28 * copt74);
  Real copt131 = -2 * copt25 * copt8 * copt84;
  Real copt132 = copt15 * copt25 * copt93;
  Real copt133 = -2 * copt18 * copt25 * copt93;
  Real copt134 = copt15 * copt28 * copt93;
  Real copt135 = -(copt113 * copt74);
  Real copt136 = copt113 * copt15 * copt18;
  Real copt137 = -(copt100 * copt125);
  Real copt140 = copt113 * copt8;
  Real copt141 = copt138 + copt139 + copt140;
  Real copt142 = copt141 * copt4;
  Real copt143 = -(copt28 * copt84);
  Real copt147 = copt146 * copt4;
  Real copt148 = copt139 + copt143 + copt144 + copt147;
  Real copt149 = copt148 * copt2;
  Real copt150 = copt129 + copt130 + copt131 + copt132 + copt133 + copt134 +
                 copt135 + copt136 + copt137 + copt142 + copt149;
  Real copt151 = copt150 * copt23;
  Real copt152 = copt15 * copt4 * copt8;
  Real copt153 = -(copt100 * copt18);
  Real copt154 = -(copt18 * copt76);
  Real copt155 = copt15 * copt25 * copt28;
  Real copt156 = copt15 * copt4 * copt84;
  Real copt157 = -2 * copt15 * copt8 * copt84;
  Real copt158 = copt18 * copt4 * copt84;
  Real copt159 = -(copt100 * copt93);
  Real copt160 = -(copt76 * copt93);
  Real copt161 = copt4 * copt8 * copt93;
  Real copt162 = copt25 * copt28 * copt93;
  Real copt164 = copt106 * copt15;
  Real copt165 = -(copt8 * copt93);
  Real copt166 = -2 * copt15;
  Real copt167 = copt166 + copt18 + copt93;
  Real copt168 = copt167 * copt4;
  Real copt169 = copt163 + copt164 + copt165 + copt168;
  Real copt170 = copt169 * copt2;
  Real copt171 = copt113 * copt15 * copt25;
  Real copt172 = copt113 * copt18 * copt25;
  Real copt173 = -2 * copt113 * copt15 * copt28;
  Real copt174 = -(copt28 * copt93);
  Real copt175 = copt109 * copt25;
  Real copt176 = -(copt113 * copt18);
  Real copt177 = copt146 * copt15;
  Real copt178 = copt174 + copt175 + copt176 + copt177;
  Real copt179 = copt178 * copt23;
  Real copt180 = copt152 + copt153 + copt154 + copt155 + copt156 + copt157 +
                 copt158 + copt159 + copt160 + copt161 + copt162 + copt170 +
                 copt171 + copt172 + copt173 + copt179;
  Real copt181 = copt13 * copt180;
  Real copt182 = copt101 + copt102 + copt103 + copt112 + copt114 + copt115 +
                 copt116 + copt117 + copt118 + copt119 + copt120 + copt121 +
                 copt122 + copt128 + copt151 + copt181 + copt75 + copt77 +
                 copt78 + copt79 + copt80 + copt81 + copt82 + copt83 + copt85 +
                 copt86 + copt87 + copt88 + copt89 + copt90 + copt91 + copt92 +
                 copt94 + copt95 + copt96 + copt97 + copt98 + copt99;
  Real copt183 = -(copt182 * copt72);
  Real copt184 = -2 * copt2 * copt4;
  Real copt185 = -2 * copt13 * copt15;
  Real copt186 = -2 * copt23 * copt25;
  Real copt187 = copt100 + copt104 + copt123 + copt184 + copt185 + copt186 +
                 copt73 + copt74 + copt76;
  Real copt190 = copt18 * copt25 * copt84;
  Real copt191 = -(copt15 * copt28 * copt84);
  Real copt192 = copt2 * copt25 * copt93;
  Real copt193 = -(copt25 * copt8 * copt93);
  Real copt194 = -(copt2 * copt28 * copt93);
  Real copt195 = copt28 * copt4 * copt93;
  Real copt197 = copt196 + copt84;
  Real copt198 = copt15 * copt197;
  Real copt199 = -copt93;
  Real copt200 = copt18 + copt199;
  Real copt201 = copt200 * copt4;
  Real copt202 = copt8 * copt93;
  Real copt203 = copt163 + copt198 + copt201 + copt202;
  Real copt204 = copt203 * copt23;
  Real copt205 = -(copt113 * copt15 * copt2);
  Real copt206 = copt113 * copt15 * copt8;
  Real copt207 = copt113 * copt18 * copt2;
  Real copt208 = -(copt113 * copt18 * copt4);
  Real copt209 = -copt84;
  Real copt210 = copt209 + copt8;
  Real copt211 = copt210 * copt25;
  Real copt213 = copt113 + copt212;
  Real copt214 = copt213 * copt4;
  Real copt215 = copt138 + copt144 + copt211 + copt214;
  Real copt216 = copt13 * copt215;
  Real copt217 = copt188 + copt189 + copt190 + copt191 + copt192 + copt193 +
                 copt194 + copt195 + copt204 + copt205 + copt206 + copt207 +
                 copt208 + copt216;
  Real copt218 = copt187 * copt217;
  Real copt219 = ArcTan(copt183, copt218);
  Real copt220 = -copt219;
  Real copt221 = copt220 + thetarest0;
  Real copt222 = t0(0);
  Real copt531 = Power(copt222, 2);
  Real copt234 = -(copt23 * copt25 * copt4 * copt8);
  Real copt236 = -(copt235 * copt74);
  Real copt237 = copt23 * copt235 * copt25;
  Real copt238 = -(copt235 * copt76);
  Real copt239 = -(copt15 * copt18 * copt23 * copt25);
  Real copt240 = 2 * copt15 * copt18 * copt4 * copt8;
  Real copt242 = -(copt100 * copt241);
  Real copt243 = copt23 * copt241 * copt25;
  Real copt244 = -(copt241 * copt76);
  Real copt245 = copt100 * copt23 * copt28;
  Real copt246 = copt23 * copt28 * copt74;
  Real copt247 = -(copt23 * copt28 * copt4 * copt8);
  Real copt248 = 2 * copt25 * copt28 * copt4 * copt8;
  Real copt249 = -(copt15 * copt18 * copt23 * copt28);
  Real copt250 = 2 * copt15 * copt18 * copt25 * copt28;
  Real copt252 = -(copt100 * copt251);
  Real copt253 = -(copt251 * copt74);
  Real copt255 = copt23 * copt25 * copt254 * copt4;
  Real copt256 = copt254 * copt74 * copt8;
  Real copt257 = -(copt23 * copt25 * copt254 * copt8);
  Real copt258 = copt254 * copt76 * copt8;
  Real copt259 = -(copt15 * copt18 * copt254 * copt4);
  Real copt260 = -(copt15 * copt18 * copt254 * copt8);
  Real copt261 = copt241 * copt254 * copt4;
  Real copt262 = -(copt23 * copt254 * copt28 * copt4);
  Real copt263 = -(copt25 * copt254 * copt28 * copt4);
  Real copt264 = copt23 * copt254 * copt28 * copt8;
  Real copt265 = -(copt25 * copt254 * copt28 * copt8);
  Real copt266 = copt251 * copt254 * copt4;
  Real copt268 = copt15 * copt23 * copt25 * copt267;
  Real copt269 = -(copt15 * copt267 * copt4 * copt8);
  Real copt270 = copt15 * copt235 * copt267;
  Real copt271 = copt100 * copt18 * copt267;
  Real copt272 = -(copt18 * copt23 * copt25 * copt267);
  Real copt273 = copt18 * copt267 * copt76;
  Real copt274 = -(copt18 * copt267 * copt4 * copt8);
  Real copt275 = -(copt15 * copt23 * copt267 * copt28);
  Real copt276 = -(copt15 * copt25 * copt267 * copt28);
  Real copt277 = copt18 * copt23 * copt267 * copt28;
  Real copt278 = -(copt18 * copt25 * copt267 * copt28);
  Real copt279 = copt15 * copt251 * copt267;
  Real copt281 = -(copt100 * copt23 * copt280);
  Real copt282 = -(copt23 * copt280 * copt74);
  Real copt283 = 2 * copt23 * copt280 * copt4 * copt8;
  Real copt284 = -(copt25 * copt280 * copt4 * copt8);
  Real copt285 = -(copt23 * copt235 * copt280);
  Real copt286 = copt235 * copt25 * copt280;
  Real copt287 = 2 * copt15 * copt18 * copt23 * copt280;
  Real copt288 = -(copt15 * copt18 * copt25 * copt280);
  Real copt289 = -(copt23 * copt241 * copt280);
  Real copt290 = copt241 * copt25 * copt280;
  Real copt291 = copt100 * copt28 * copt280;
  Real copt292 = copt28 * copt280 * copt74;
  Real copt293 = -(copt28 * copt280 * copt4 * copt8);
  Real copt294 = -(copt15 * copt18 * copt28 * copt280);
  Real copt295 = copt18 * copt76;
  Real copt296 = -(copt18 * copt25 * copt28);
  Real copt297 = copt18 * copt254 * copt8;
  Real copt300 = copt15 * copt299;
  Real copt301 = copt18 * copt254;
  Real copt302 = -2 * copt267;
  Real copt303 = copt18 + copt302;
  Real copt304 = copt303 * copt8;
  Real copt305 = copt300 + copt301 + copt304;
  Real copt306 = -(copt305 * copt4);
  Real copt309 = copt100 * copt308;
  Real copt310 = -(copt267 * copt76);
  Real copt311 = -(copt235 * copt267);
  Real copt312 = 2 * copt25 * copt267 * copt28;
  Real copt313 = -(copt251 * copt267);
  Real copt314 = -(copt254 * copt8);
  Real copt317 = -(copt230 * copt316);
  Real copt318 = copt235 + copt314 + copt317;
  Real copt319 = copt15 * copt318;
  Real copt320 = -(copt18 * copt25 * copt280);
  Real copt321 = copt18 * copt28 * copt280;
  Real copt322 = copt295 + copt296 + copt297 + copt306 + copt309 + copt310 +
                 copt311 + copt312 + copt313 + copt319 + copt320 + copt321;
  Real copt323 = copt13 * copt322;
  Real copt324 = copt241 * copt4;
  Real copt325 = copt251 * copt4;
  Real copt326 = copt299 * copt74;
  Real copt327 = copt299 * copt76;
  Real copt328 = -(copt241 * copt254);
  Real copt329 = -(copt251 * copt254);
  Real copt330 = -(copt18 * copt267 * copt4);
  Real copt331 = copt18 * copt267 * copt8;
  Real copt332 = -2 * copt18 * copt254;
  Real copt334 = copt18 + copt267;
  Real copt335 = copt334 * copt8;
  Real copt336 = copt332 + copt333 + copt335;
  Real copt337 = -(copt15 * copt336);
  Real copt338 = -(copt28 * copt280 * copt4);
  Real copt339 = copt28 * copt280 * copt8;
  Real copt340 = -2 * copt254 * copt28;
  Real copt341 = copt316 * copt4;
  Real copt342 = copt28 + copt280;
  Real copt343 = copt342 * copt8;
  Real copt344 = copt340 + copt341 + copt343;
  Real copt345 = -(copt25 * copt344);
  Real copt346 = copt324 + copt325 + copt326 + copt327 + copt328 + copt329 +
                 copt330 + copt331 + copt337 + copt338 + copt339 + copt345;
  Real copt347 = copt2 * copt346;
  Real copt348 = copt234 + copt236 + copt237 + copt238 + copt239 + copt240 +
                 copt242 + copt243 + copt244 + copt245 + copt246 + copt247 +
                 copt248 + copt249 + copt250 + copt252 + copt253 + copt255 +
                 copt256 + copt257 + copt258 + copt259 + copt260 + copt261 +
                 copt262 + copt263 + copt264 + copt265 + copt266 + copt268 +
                 copt269 + copt270 + copt271 + copt272 + copt273 + copt274 +
                 copt275 + copt276 + copt277 + copt278 + copt279 + copt281 +
                 copt282 + copt283 + copt284 + copt285 + copt286 + copt287 +
                 copt288 + copt289 + copt290 + copt291 + copt292 + copt293 +
                 copt294 + copt323 + copt347;
  Real copt349 = copt233 * copt348;
  Real copt350 = -2 * copt4 * copt8;
  Real copt351 = -2 * copt15 * copt18;
  Real copt352 = -2 * copt25 * copt28;
  Real copt353 = copt100 + copt235 + copt241 + copt251 + copt350 + copt351 +
                 copt352 + copt74 + copt76;
  Real copt354 = copt18 * copt25 * copt254;
  Real copt355 = -(copt15 * copt254 * copt28);
  Real copt356 = copt2 * copt25 * copt267;
  Real copt357 = -(copt25 * copt267 * copt8);
  Real copt358 = -(copt2 * copt267 * copt28);
  Real copt359 = copt267 * copt28 * copt4;
  Real copt360 = -(copt18 * copt254);
  Real copt361 = copt196 + copt254;
  Real copt362 = copt15 * copt361;
  Real copt363 = copt267 * copt8;
  Real copt364 = copt333 + copt360 + copt362 + copt363;
  Real copt365 = copt23 * copt364;
  Real copt366 = -(copt15 * copt2 * copt280);
  Real copt367 = copt15 * copt280 * copt8;
  Real copt368 = copt18 * copt2 * copt280;
  Real copt369 = -(copt18 * copt280 * copt4);
  Real copt370 = copt25 * copt299;
  Real copt371 = copt254 * copt28;
  Real copt372 = -(copt280 * copt8);
  Real copt373 = copt212 + copt280;
  Real copt374 = copt373 * copt4;
  Real copt375 = copt370 + copt371 + copt372 + copt374;
  Real copt376 = copt13 * copt375;
  Real copt377 = copt188 + copt189 + copt354 + copt355 + copt356 + copt357 +
                 copt358 + copt359 + copt365 + copt366 + copt367 + copt368 +
                 copt369 + copt376;
  Real copt378 = copt353 * copt377;
  Real copt379 = ArcTan(copt349, copt378);
  Real copt380 = -copt379;
  Real copt381 = copt380 + thetarest1;
  Real copt382 = t1(0);
  Real copt533 = Power(copt382, 2);
  Real copt393 = copt15 * copt18 * copt73;
  Real copt394 = -(copt15 * copt18 * copt2 * copt8);
  Real copt395 = -(copt241 * copt73);
  Real copt396 = copt2 * copt241 * copt4;
  Real copt397 = copt25 * copt28 * copt73;
  Real copt398 = -(copt2 * copt25 * copt28 * copt8);
  Real copt399 = -(copt251 * copt73);
  Real copt400 = copt2 * copt251 * copt4;
  Real copt402 = -(copt15 * copt18 * copt2 * copt401);
  Real copt403 = copt15 * copt18 * copt401 * copt8;
  Real copt404 = copt2 * copt241 * copt401;
  Real copt405 = -(copt241 * copt4 * copt401);
  Real copt406 = -(copt2 * copt25 * copt28 * copt401);
  Real copt407 = copt25 * copt28 * copt401 * copt8;
  Real copt408 = copt2 * copt251 * copt401;
  Real copt409 = -(copt251 * copt4 * copt401);
  Real copt418 = copt228 * copt417;
  Real copt419 = copt410 + copt413 + copt414 + copt418;
  Real copt420 = copt104 * copt419;
  Real copt421 = -(copt15 * copt415 * copt73);
  Real copt422 = 2 * copt15 * copt2 * copt415 * copt8;
  Real copt423 = -(copt15 * copt235 * copt415);
  Real copt424 = copt18 * copt415 * copt73;
  Real copt425 = -(copt18 * copt2 * copt4 * copt415);
  Real copt426 = -(copt18 * copt2 * copt415 * copt8);
  Real copt427 = copt18 * copt4 * copt415 * copt8;
  Real copt428 = copt18 * copt25 * copt28 * copt415;
  Real copt429 = -(copt15 * copt251 * copt415);
  Real copt433 = copt230 * copt432;
  Real copt434 = copt410 + copt413 + copt414 + copt433;
  Real copt435 = copt123 * copt434;
  Real copt436 = -(copt25 * copt430 * copt73);
  Real copt437 = 2 * copt2 * copt25 * copt430 * copt8;
  Real copt438 = -(copt235 * copt25 * copt430);
  Real copt439 = -(copt241 * copt25 * copt430);
  Real copt440 = copt28 * copt430 * copt73;
  Real copt441 = -(copt2 * copt28 * copt4 * copt430);
  Real copt442 = -(copt2 * copt28 * copt430 * copt8);
  Real copt443 = copt28 * copt4 * copt430 * copt8;
  Real copt444 = copt15 * copt18 * copt28 * copt430;
  Real copt445 = -(copt15 * copt235);
  Real copt446 = copt18 * copt23 * copt25;
  Real copt447 = copt18 * copt4 * copt8;
  Real copt448 = -2 * copt18 * copt23 * copt28;
  Real copt449 = copt18 * copt25 * copt28;
  Real copt450 = copt15 * copt401 * copt8;
  Real copt451 = -2 * copt18 * copt4 * copt401;
  Real copt452 = copt18 * copt401 * copt8;
  Real copt453 = -(copt23 * copt25 * copt415);
  Real copt454 = copt4 * copt415 * copt8;
  Real copt455 = -(copt235 * copt415);
  Real copt456 = copt23 * copt28 * copt415;
  Real copt457 = copt25 * copt28 * copt415;
  Real copt458 = -(copt251 * copt415);
  Real copt459 = -2 * copt18 * copt8;
  Real copt460 = copt15 * copt412;
  Real copt461 = copt18 * copt401;
  Real copt464 = copt459 + copt460 + copt461 + copt462 + copt463;
  Real copt465 = copt2 * copt464;
  Real copt466 = copt15 * copt389 * copt432;
  Real copt467 = copt18 * copt23 * copt430;
  Real copt468 = -2 * copt18 * copt25 * copt430;
  Real copt469 = copt18 * copt28 * copt430;
  Real copt470 = copt445 + copt446 + copt447 + copt448 + copt449 + copt450 +
                 copt451 + copt452 + copt453 + copt454 + copt455 + copt456 +
                 copt457 + copt458 + copt465 + copt466 + copt467 + copt468 +
                 copt469;
  Real copt471 = -(copt13 * copt470);
  Real copt472 = copt28 * copt4 * copt8;
  Real copt473 = copt15 * copt18 * copt28;
  Real copt474 = -2 * copt28 * copt4 * copt401;
  Real copt475 = copt28 * copt401 * copt8;
  Real copt476 = -(copt401 * copt8);
  Real copt477 = copt18 * copt417;
  Real copt478 = copt235 + copt476 + copt477;
  Real copt479 = -(copt25 * copt478);
  Real copt480 = -2 * copt15 * copt28 * copt415;
  Real copt481 = copt18 * copt28 * copt415;
  Real copt482 = copt4 * copt430 * copt8;
  Real copt483 = -(copt235 * copt430);
  Real copt484 = copt15 * copt18 * copt430;
  Real copt485 = -(copt241 * copt430);
  Real copt486 = -2 * copt28 * copt8;
  Real copt489 = copt4 * copt432;
  Real copt490 = copt430 * copt8;
  Real copt491 = copt486 + copt487 + copt488 + copt489 + copt490;
  Real copt492 = copt2 * copt491;
  Real copt493 = copt472 + copt473 + copt474 + copt475 + copt479 + copt480 +
                 copt481 + copt482 + copt483 + copt484 + copt485 + copt492;
  Real copt494 = -(copt23 * copt493);
  Real copt495 = copt393 + copt394 + copt395 + copt396 + copt397 + copt398 +
                 copt399 + copt400 + copt402 + copt403 + copt404 + copt405 +
                 copt406 + copt407 + copt408 + copt409 + copt420 + copt421 +
                 copt422 + copt423 + copt424 + copt425 + copt426 + copt427 +
                 copt428 + copt429 + copt435 + copt436 + copt437 + copt438 +
                 copt439 + copt440 + copt441 + copt442 + copt443 + copt444 +
                 copt471 + copt494;
  Real copt496 = copt392 * copt495;
  Real copt497 = -2 * copt2 * copt8;
  Real copt498 = -2 * copt13 * copt18;
  Real copt499 = -2 * copt23 * copt28;
  Real copt500 = copt104 + copt123 + copt235 + copt241 + copt251 + copt497 +
                 copt498 + copt499 + copt73;
  Real copt501 = copt18 * copt25 * copt401;
  Real copt502 = -(copt15 * copt28 * copt401);
  Real copt503 = copt2 * copt25 * copt415;
  Real copt504 = -(copt25 * copt415 * copt8);
  Real copt505 = -(copt2 * copt28 * copt415);
  Real copt506 = copt28 * copt4 * copt415;
  Real copt507 = -(copt18 * copt401);
  Real copt508 = copt196 + copt401;
  Real copt509 = copt15 * copt508;
  Real copt510 = copt462 + copt463 + copt507 + copt509;
  Real copt511 = copt23 * copt510;
  Real copt512 = -(copt15 * copt2 * copt430);
  Real copt513 = copt15 * copt430 * copt8;
  Real copt514 = copt18 * copt2 * copt430;
  Real copt515 = -(copt18 * copt4 * copt430);
  Real copt516 = -(copt430 * copt8);
  Real copt517 = copt212 + copt430;
  Real copt518 = copt4 * copt517;
  Real copt519 = copt487 + copt488 + copt516 + copt518;
  Real copt520 = copt13 * copt519;
  Real copt521 = copt188 + copt189 + copt501 + copt502 + copt503 + copt504 +
                 copt505 + copt506 + copt511 + copt512 + copt513 + copt514 +
                 copt515 + copt520;
  Real copt522 = copt500 * copt521;
  Real copt523 = ArcTan(copt496, copt522);
  Real copt524 = -copt523;
  Real copt525 = copt524 + thetarest2;
  Real copt526 = t2(0);
  Real copt535 = Power(copt526, 2);
  Real copt540 = Power(copt50, 2);
  Real copt541 = -copt540;
  Real copt542 = copt33 * copt54;
  Real copt543 = copt541 + copt542;
  Real copt544 = 1 / copt543;
  Real copt546 = copt36 * copt7;
  Real copt547 = -(copt1 * copt38);
  Real copt548 = copt546 + copt547;
  Real copt549 = Power(copt548, 2);
  Real copt550 = 1 / copt549;
  Real copt223 = t0(1);
  Real copt551 = Power(copt223, 2);
  Real copt383 = t1(1);
  Real copt554 = Power(copt383, 2);
  Real copt527 = t2(1);
  Real copt557 = Power(copt527, 2);
  Real copt569 = -2 * copt2 * copt74 * copt8;
  Real copt570 = -2 * copt2 * copt76 * copt8;
  Real copt571 = copt235 * copt74;
  Real copt572 = copt235 * copt76;
  Real copt573 = -2 * copt15 * copt18 * copt73;
  Real copt574 = 2 * copt15 * copt18 * copt2 * copt4;
  Real copt575 = 2 * copt15 * copt18 * copt2 * copt8;
  Real copt576 = -2 * copt15 * copt18 * copt4 * copt8;
  Real copt577 = copt241 * copt73;
  Real copt578 = -2 * copt2 * copt241 * copt4;
  Real copt579 = copt100 * copt241;
  Real copt580 = copt241 * copt76;
  Real copt581 = copt100 + copt235 + copt241 + copt350 + copt351 + copt74;
  Real copt582 = copt104 * copt581;
  Real copt583 = -2 * copt25 * copt28 * copt73;
  Real copt584 = 2 * copt2 * copt25 * copt28 * copt4;
  Real copt585 = 2 * copt2 * copt25 * copt28 * copt8;
  Real copt586 = -2 * copt25 * copt28 * copt4 * copt8;
  Real copt587 = -2 * copt15 * copt18 * copt25 * copt28;
  Real copt588 = copt251 * copt73;
  Real copt589 = -2 * copt2 * copt251 * copt4;
  Real copt590 = copt100 * copt251;
  Real copt591 = copt251 * copt74;
  Real copt592 = copt100 + copt235 + copt251 + copt350 + copt352 + copt76;
  Real copt593 = copt123 * copt592;
  Real copt594 = -(copt15 * copt4 * copt8);
  Real copt595 = copt15 * copt235;
  Real copt596 = copt2 * copt225 * copt228;
  Real copt597 = copt100 * copt18;
  Real copt598 = -(copt18 * copt4 * copt8);
  Real copt599 = copt228 * copt23 * copt230;
  Real copt600 = -(copt15 * copt25 * copt28);
  Real copt601 = copt15 * copt251;
  Real copt602 = copt295 + copt296 + copt594 + copt595 + copt596 + copt597 +
                 copt598 + copt599 + copt600 + copt601;
  Real copt603 = -2 * copt13 * copt602;
  Real copt604 = copt235 * copt25;
  Real copt605 = -(copt15 * copt18 * copt25);
  Real copt606 = copt241 * copt25;
  Real copt607 = copt2 * copt225 * copt230;
  Real copt608 = copt100 * copt28;
  Real copt609 = copt28 * copt74;
  Real copt610 = -(copt15 * copt18 * copt28);
  Real copt611 = copt25 + copt28;
  Real copt612 = -(copt4 * copt611 * copt8);
  Real copt613 = copt604 + copt605 + copt606 + copt607 + copt608 + copt609 +
                 copt610 + copt612;
  Real copt614 = -2 * copt23 * copt613;
  Real copt615 = copt569 + copt570 + copt571 + copt572 + copt573 + copt574 +
                 copt575 + copt576 + copt577 + copt578 + copt579 + copt580 +
                 copt582 + copt583 + copt584 + copt585 + copt586 + copt587 +
                 copt588 + copt589 + copt590 + copt591 + copt593 + copt603 +
                 copt614 + copt75 + copt77;
  Real copt616 = 1 / copt615;
  Real copt618 = -(copt23 * copt25);
  Real copt619 = copt4 * copt8;
  Real copt620 = copt4 + copt8;
  Real copt621 = -(copt2 * copt620);
  Real copt622 = copt15 * copt18;
  Real copt623 = copt15 + copt18;
  Real copt624 = -(copt13 * copt623);
  Real copt625 = -(copt23 * copt28);
  Real copt626 = copt25 * copt28;
  Real copt627 = copt104 + copt123 + copt618 + copt619 + copt621 + copt622 +
                 copt624 + copt625 + copt626 + copt73;
  Real copt637 = -thetarest0;
  Real copt638 = copt219 + copt637;
  Real copt640 = -thetarest1;
  Real copt641 = copt379 + copt640;
  Real copt643 = -thetarest2;
  Real copt644 = copt523 + copt643;
  Real copt639 = (copt222 * copt223 * copt58 * copt59 * copt638) / 2.;
  Real copt642 = (copt382 * copt383 * copt58 * copt60 * copt641) / 2.;
  Real copt645 = (copt526 * copt527 * copt58 * copt61 * copt644) / 2.;
  Real copt646 = copt639 + copt642 + copt645;
  out(0)       = copt34;
  out(1)       = copt35 * copt50 * copt56;
  out(2)       = copt55;
  out(3) =
      (copt544 * copt58 * copt59 * copt60 * copt61 *
       (-(copt50 * (copt525 * copt526 * copt527 * l0 * l1 +
                    copt381 * copt382 * copt383 * l0 * l2 +
                    copt221 * copt222 * copt223 * l1 * l2)) +
        copt54 * (copt525 * copt535 * l0 * l1 + copt381 * copt533 * l0 * l2 +
                  copt221 * copt531 * l1 * l2))) /
      2.;
  out(4) =
      -(copt54 * copt544 * copt646) -
      (copt550 * copt58 * copt59 * copt60 * copt61 * copt616 *
       (copt1 * (copt187 * copt36 + copt38 * copt627) +
        (copt38 * copt500 + copt36 * copt627) * copt7) *
       (-(copt379 * copt554 * l0 * l2) - copt219 * copt551 * l1 * l2 +
        copt551 * l1 * l2 * thetarest0 + copt554 * l0 * l2 * thetarest1 +
        copt557 * l0 * l1 * thetarest2 -
        copt557 * l0 * l1 *
            ArcTan(copt392 * (copt393 + copt394 + copt395 + copt396 + copt397 +
                              copt398 + copt399 + copt400 + copt402 + copt403 +
                              copt404 + copt405 + copt406 + copt407 + copt408 +
                              copt409 + copt420 + copt421 + copt422 + copt423 +
                              copt424 + copt425 + copt426 + copt427 + copt428 +
                              copt429 + copt435 + copt436 + copt437 + copt438 +
                              copt439 + copt440 + copt441 + copt442 + copt443 +
                              copt444 + copt471 -
                              copt23 * (copt472 + copt473 + copt474 + copt475 -
                                        copt25 * (copt235 + copt241 -
                                                  copt18 * copt415 + copt476) +
                                        copt480 + copt481 + copt482 + copt483 +
                                        copt484 + copt485 + copt492)),
                   copt522))) /
          2.;
  out(5) = -(copt33 * copt544 *
             ((copt551 * copt58 * copt59 * copt638) / 2. +
              (copt554 * copt58 * copt60 * copt641) / 2. +
              (copt557 * copt58 * copt61 * copt644) / 2.)) -
           (-(copt11 * copt40) - copt21 * copt44 - copt31 * copt48) * copt544 *
               copt646;
  return out;
}

#endif  // hylc_strain_II
