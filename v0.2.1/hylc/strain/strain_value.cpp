#include "strain.hpp"

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
  Real copt62  = t0(0);
  Real copt63  = Power(copt62, 2);
  Real copt65  = -copt4;
  Real copt66  = copt2 + copt65;
  Real copt67  = Power(copt66, 2);
  Real copt68  = -copt15;
  Real copt69  = copt13 + copt68;
  Real copt70  = Power(copt69, 2);
  Real copt71  = -copt25;
  Real copt72  = copt23 + copt71;
  Real copt73  = Power(copt72, 2);
  Real copt74  = copt67 + copt70 + copt73;
  Real copt75  = Sqrt(copt74);
  Real copt76  = Power(copt2, 2);
  Real copt77  = Power(copt15, 2);
  Real copt79  = Power(copt25, 2);
  Real copt87  = xloc(9);
  Real copt96  = xloc(10);
  Real copt103 = Power(copt4, 2);
  Real copt107 = Power(copt23, 2);
  Real copt116 = xloc(11);
  Real copt126 = Power(copt13, 2);
  Real copt108 = copt8 * copt87;
  Real copt109 = copt8 + copt87;
  Real copt110 = -(copt109 * copt4);
  Real copt128 = copt116 + copt28;
  Real copt142 = copt109 * copt25;
  Real copt112 = copt18 + copt96;
  Real copt148 = -2 * copt25;
  Real copt149 = copt116 + copt148 + copt28;
  Real copt166 = -(copt18 * copt87);
  Real copt141 = copt28 * copt87;
  Real copt147 = -(copt116 * copt8);
  Real copt224 = t1(0);
  Real copt225 = Power(copt224, 2);
  Real copt199 = -copt8;
  Real copt227 = copt199 + copt4;
  Real copt228 = Power(copt227, 2);
  Real copt229 = -copt18;
  Real copt230 = copt15 + copt229;
  Real copt231 = Power(copt230, 2);
  Real copt215 = -copt28;
  Real copt232 = copt215 + copt25;
  Real copt233 = Power(copt232, 2);
  Real copt234 = copt228 + copt231 + copt233;
  Real copt235 = Sqrt(copt234);
  Real copt237 = Power(copt8, 2);
  Real copt243 = Power(copt18, 2);
  Real copt253 = Power(copt28, 2);
  Real copt256 = xloc(12);
  Real copt269 = xloc(13);
  Real copt282 = xloc(14);
  Real copt300 = -copt256;
  Real copt301 = copt300 + copt8;
  Real copt309 = -copt269;
  Real copt310 = copt18 + copt309;
  Real copt317 = -copt282;
  Real copt318 = copt28 + copt317;
  Real copt191 = -(copt18 * copt2 * copt25);
  Real copt192 = copt15 * copt2 * copt28;
  Real copt335 = copt310 * copt4;
  Real copt383 = t2(0);
  Real copt384 = Power(copt383, 2);
  Real copt386 = copt199 + copt2;
  Real copt387 = Power(copt386, 2);
  Real copt388 = copt13 + copt229;
  Real copt389 = Power(copt388, 2);
  Real copt390 = copt215 + copt23;
  Real copt391 = Power(copt390, 2);
  Real copt392 = copt387 + copt389 + copt391;
  Real copt393 = Sqrt(copt392);
  Real copt402 = xloc(15);
  Real copt416 = xloc(16);
  Real copt411 = -copt237;
  Real copt412 = -copt402;
  Real copt413 = copt412 + copt8;
  Real copt414 = copt4 * copt413;
  Real copt415 = copt402 * copt8;
  Real copt431 = xloc(17);
  Real copt417 = -copt416;
  Real copt418 = copt18 + copt417;
  Real copt432 = -copt431;
  Real copt433 = copt28 + copt432;
  Real copt463 = copt4 * copt418;
  Real copt464 = copt416 * copt8;
  Real copt488 = copt25 * copt413;
  Real copt489 = copt28 * copt402;
  Real copt78  = copt76 * copt77;
  Real copt80  = copt76 * copt79;
  Real copt81  = -(copt2 * copt77 * copt8);
  Real copt82  = -(copt2 * copt79 * copt8);
  Real copt83  = -(copt15 * copt18 * copt76);
  Real copt84  = copt15 * copt18 * copt2 * copt4;
  Real copt85  = -(copt25 * copt28 * copt76);
  Real copt86  = copt2 * copt25 * copt28 * copt4;
  Real copt88  = -(copt2 * copt77 * copt87);
  Real copt89  = -(copt2 * copt79 * copt87);
  Real copt90  = copt77 * copt8 * copt87;
  Real copt91  = copt79 * copt8 * copt87;
  Real copt92  = copt15 * copt18 * copt2 * copt87;
  Real copt93  = -(copt15 * copt18 * copt4 * copt87);
  Real copt94  = copt2 * copt25 * copt28 * copt87;
  Real copt95  = -(copt25 * copt28 * copt4 * copt87);
  Real copt97  = -(copt15 * copt76 * copt96);
  Real copt98  = copt15 * copt2 * copt4 * copt96;
  Real copt99  = copt15 * copt2 * copt8 * copt96;
  Real copt100 = -(copt15 * copt4 * copt8 * copt96);
  Real copt101 = copt18 * copt76 * copt96;
  Real copt102 = -2 * copt18 * copt2 * copt4 * copt96;
  Real copt104 = copt103 * copt18 * copt96;
  Real copt105 = copt18 * copt79 * copt96;
  Real copt106 = -(copt15 * copt25 * copt28 * copt96);
  Real copt111 = copt18 * copt96;
  Real copt113 = -(copt112 * copt15);
  Real copt114 = copt103 + copt108 + copt110 + copt111 + copt113 + copt77;
  Real copt115 = copt107 * copt114;
  Real copt117 = -(copt116 * copt25 * copt76);
  Real copt118 = copt116 * copt2 * copt25 * copt4;
  Real copt119 = copt116 * copt2 * copt25 * copt8;
  Real copt120 = -(copt116 * copt25 * copt4 * copt8);
  Real copt121 = -(copt116 * copt15 * copt18 * copt25);
  Real copt122 = copt116 * copt28 * copt76;
  Real copt123 = -2 * copt116 * copt2 * copt28 * copt4;
  Real copt124 = copt103 * copt116 * copt28;
  Real copt125 = copt116 * copt28 * copt77;
  Real copt127 = copt116 * copt28;
  Real copt129 = -(copt128 * copt25);
  Real copt130 = copt103 + copt108 + copt110 + copt127 + copt129 + copt79;
  Real copt131 = copt126 * copt130;
  Real copt132 = copt15 * copt18 * copt25;
  Real copt133 = -(copt28 * copt77);
  Real copt134 = -2 * copt25 * copt8 * copt87;
  Real copt135 = copt15 * copt25 * copt96;
  Real copt136 = -2 * copt18 * copt25 * copt96;
  Real copt137 = copt15 * copt28 * copt96;
  Real copt138 = -(copt116 * copt77);
  Real copt139 = copt116 * copt15 * copt18;
  Real copt140 = -(copt103 * copt128);
  Real copt143 = copt116 * copt8;
  Real copt144 = copt141 + copt142 + copt143;
  Real copt145 = copt144 * copt4;
  Real copt146 = -(copt28 * copt87);
  Real copt150 = copt149 * copt4;
  Real copt151 = copt142 + copt146 + copt147 + copt150;
  Real copt152 = copt151 * copt2;
  Real copt153 = copt132 + copt133 + copt134 + copt135 + copt136 + copt137 +
                 copt138 + copt139 + copt140 + copt145 + copt152;
  Real copt154 = copt153 * copt23;
  Real copt155 = copt15 * copt4 * copt8;
  Real copt156 = -(copt103 * copt18);
  Real copt157 = -(copt18 * copt79);
  Real copt158 = copt15 * copt25 * copt28;
  Real copt159 = copt15 * copt4 * copt87;
  Real copt160 = -2 * copt15 * copt8 * copt87;
  Real copt161 = copt18 * copt4 * copt87;
  Real copt162 = -(copt103 * copt96);
  Real copt163 = -(copt79 * copt96);
  Real copt164 = copt4 * copt8 * copt96;
  Real copt165 = copt25 * copt28 * copt96;
  Real copt167 = copt109 * copt15;
  Real copt168 = -(copt8 * copt96);
  Real copt169 = -2 * copt15;
  Real copt170 = copt169 + copt18 + copt96;
  Real copt171 = copt170 * copt4;
  Real copt172 = copt166 + copt167 + copt168 + copt171;
  Real copt173 = copt172 * copt2;
  Real copt174 = copt116 * copt15 * copt25;
  Real copt175 = copt116 * copt18 * copt25;
  Real copt176 = -2 * copt116 * copt15 * copt28;
  Real copt177 = -(copt28 * copt96);
  Real copt178 = copt112 * copt25;
  Real copt179 = -(copt116 * copt18);
  Real copt180 = copt149 * copt15;
  Real copt181 = copt177 + copt178 + copt179 + copt180;
  Real copt182 = copt181 * copt23;
  Real copt183 = copt155 + copt156 + copt157 + copt158 + copt159 + copt160 +
                 copt161 + copt162 + copt163 + copt164 + copt165 + copt173 +
                 copt174 + copt175 + copt176 + copt182;
  Real copt184 = copt13 * copt183;
  Real copt185 = copt100 + copt101 + copt102 + copt104 + copt105 + copt106 +
                 copt115 + copt117 + copt118 + copt119 + copt120 + copt121 +
                 copt122 + copt123 + copt124 + copt125 + copt131 + copt154 +
                 copt184 + copt78 + copt80 + copt81 + copt82 + copt83 + copt84 +
                 copt85 + copt86 + copt88 + copt89 + copt90 + copt91 + copt92 +
                 copt93 + copt94 + copt95 + copt97 + copt98 + copt99;
  Real copt186 = -(copt185 * copt75);
  Real copt187 = -2 * copt2 * copt4;
  Real copt188 = -2 * copt13 * copt15;
  Real copt189 = -2 * copt23 * copt25;
  Real copt190 = copt103 + copt107 + copt126 + copt187 + copt188 + copt189 +
                 copt76 + copt77 + copt79;
  Real copt193 = copt18 * copt25 * copt87;
  Real copt194 = -(copt15 * copt28 * copt87);
  Real copt195 = copt2 * copt25 * copt96;
  Real copt196 = -(copt25 * copt8 * copt96);
  Real copt197 = -(copt2 * copt28 * copt96);
  Real copt198 = copt28 * copt4 * copt96;
  Real copt200 = copt199 + copt87;
  Real copt201 = copt15 * copt200;
  Real copt202 = -copt96;
  Real copt203 = copt18 + copt202;
  Real copt204 = copt203 * copt4;
  Real copt205 = copt8 * copt96;
  Real copt206 = copt166 + copt201 + copt204 + copt205;
  Real copt207 = copt206 * copt23;
  Real copt208 = -(copt116 * copt15 * copt2);
  Real copt209 = copt116 * copt15 * copt8;
  Real copt210 = copt116 * copt18 * copt2;
  Real copt211 = -(copt116 * copt18 * copt4);
  Real copt212 = -copt87;
  Real copt213 = copt212 + copt8;
  Real copt214 = copt213 * copt25;
  Real copt216 = copt116 + copt215;
  Real copt217 = copt216 * copt4;
  Real copt218 = copt141 + copt147 + copt214 + copt217;
  Real copt219 = copt13 * copt218;
  Real copt220 = copt191 + copt192 + copt193 + copt194 + copt195 + copt196 +
                 copt197 + copt198 + copt207 + copt208 + copt209 + copt210 +
                 copt211 + copt219;
  Real copt221 = copt190 * copt220;
  Real copt222 = ArcTan(copt186, copt221);
  Real copt528 = t0(1);
  Real copt236 = -(copt23 * copt25 * copt4 * copt8);
  Real copt238 = -(copt237 * copt77);
  Real copt239 = copt23 * copt237 * copt25;
  Real copt240 = -(copt237 * copt79);
  Real copt241 = -(copt15 * copt18 * copt23 * copt25);
  Real copt242 = 2 * copt15 * copt18 * copt4 * copt8;
  Real copt244 = -(copt103 * copt243);
  Real copt245 = copt23 * copt243 * copt25;
  Real copt246 = -(copt243 * copt79);
  Real copt247 = copt103 * copt23 * copt28;
  Real copt248 = copt23 * copt28 * copt77;
  Real copt249 = -(copt23 * copt28 * copt4 * copt8);
  Real copt250 = 2 * copt25 * copt28 * copt4 * copt8;
  Real copt251 = -(copt15 * copt18 * copt23 * copt28);
  Real copt252 = 2 * copt15 * copt18 * copt25 * copt28;
  Real copt254 = -(copt103 * copt253);
  Real copt255 = -(copt253 * copt77);
  Real copt257 = copt23 * copt25 * copt256 * copt4;
  Real copt258 = copt256 * copt77 * copt8;
  Real copt259 = -(copt23 * copt25 * copt256 * copt8);
  Real copt260 = copt256 * copt79 * copt8;
  Real copt261 = -(copt15 * copt18 * copt256 * copt4);
  Real copt262 = -(copt15 * copt18 * copt256 * copt8);
  Real copt263 = copt243 * copt256 * copt4;
  Real copt264 = -(copt23 * copt256 * copt28 * copt4);
  Real copt265 = -(copt25 * copt256 * copt28 * copt4);
  Real copt266 = copt23 * copt256 * copt28 * copt8;
  Real copt267 = -(copt25 * copt256 * copt28 * copt8);
  Real copt268 = copt253 * copt256 * copt4;
  Real copt270 = copt15 * copt23 * copt25 * copt269;
  Real copt271 = -(copt15 * copt269 * copt4 * copt8);
  Real copt272 = copt15 * copt237 * copt269;
  Real copt273 = copt103 * copt18 * copt269;
  Real copt274 = -(copt18 * copt23 * copt25 * copt269);
  Real copt275 = copt18 * copt269 * copt79;
  Real copt276 = -(copt18 * copt269 * copt4 * copt8);
  Real copt277 = -(copt15 * copt23 * copt269 * copt28);
  Real copt278 = -(copt15 * copt25 * copt269 * copt28);
  Real copt279 = copt18 * copt23 * copt269 * copt28;
  Real copt280 = -(copt18 * copt25 * copt269 * copt28);
  Real copt281 = copt15 * copt253 * copt269;
  Real copt283 = -(copt103 * copt23 * copt282);
  Real copt284 = -(copt23 * copt282 * copt77);
  Real copt285 = 2 * copt23 * copt282 * copt4 * copt8;
  Real copt286 = -(copt25 * copt282 * copt4 * copt8);
  Real copt287 = -(copt23 * copt237 * copt282);
  Real copt288 = copt237 * copt25 * copt282;
  Real copt289 = 2 * copt15 * copt18 * copt23 * copt282;
  Real copt290 = -(copt15 * copt18 * copt25 * copt282);
  Real copt291 = -(copt23 * copt243 * copt282);
  Real copt292 = copt243 * copt25 * copt282;
  Real copt293 = copt103 * copt28 * copt282;
  Real copt294 = copt28 * copt282 * copt77;
  Real copt295 = -(copt28 * copt282 * copt4 * copt8);
  Real copt296 = -(copt15 * copt18 * copt28 * copt282);
  Real copt297 = copt18 * copt79;
  Real copt298 = -(copt18 * copt25 * copt28);
  Real copt299 = copt18 * copt256 * copt8;
  Real copt302 = copt15 * copt301;
  Real copt303 = copt18 * copt256;
  Real copt304 = -2 * copt269;
  Real copt305 = copt18 + copt304;
  Real copt306 = copt305 * copt8;
  Real copt307 = copt302 + copt303 + copt306;
  Real copt308 = -(copt307 * copt4);
  Real copt311 = copt103 * copt310;
  Real copt312 = -(copt269 * copt79);
  Real copt313 = -(copt237 * copt269);
  Real copt314 = 2 * copt25 * copt269 * copt28;
  Real copt315 = -(copt253 * copt269);
  Real copt316 = -(copt256 * copt8);
  Real copt319 = -(copt232 * copt318);
  Real copt320 = copt237 + copt316 + copt319;
  Real copt321 = copt15 * copt320;
  Real copt322 = -(copt18 * copt25 * copt282);
  Real copt323 = copt18 * copt28 * copt282;
  Real copt324 = copt297 + copt298 + copt299 + copt308 + copt311 + copt312 +
                 copt313 + copt314 + copt315 + copt321 + copt322 + copt323;
  Real copt325 = copt13 * copt324;
  Real copt326 = copt243 * copt4;
  Real copt327 = copt253 * copt4;
  Real copt328 = copt301 * copt77;
  Real copt329 = copt301 * copt79;
  Real copt330 = -(copt243 * copt256);
  Real copt331 = -(copt253 * copt256);
  Real copt332 = -(copt18 * copt269 * copt4);
  Real copt333 = copt18 * copt269 * copt8;
  Real copt334 = -2 * copt18 * copt256;
  Real copt336 = copt18 + copt269;
  Real copt337 = copt336 * copt8;
  Real copt338 = copt334 + copt335 + copt337;
  Real copt339 = -(copt15 * copt338);
  Real copt340 = -(copt28 * copt282 * copt4);
  Real copt341 = copt28 * copt282 * copt8;
  Real copt342 = -2 * copt256 * copt28;
  Real copt343 = copt318 * copt4;
  Real copt344 = copt28 + copt282;
  Real copt345 = copt344 * copt8;
  Real copt346 = copt342 + copt343 + copt345;
  Real copt347 = -(copt25 * copt346);
  Real copt348 = copt326 + copt327 + copt328 + copt329 + copt330 + copt331 +
                 copt332 + copt333 + copt339 + copt340 + copt341 + copt347;
  Real copt349 = copt2 * copt348;
  Real copt350 = copt236 + copt238 + copt239 + copt240 + copt241 + copt242 +
                 copt244 + copt245 + copt246 + copt247 + copt248 + copt249 +
                 copt250 + copt251 + copt252 + copt254 + copt255 + copt257 +
                 copt258 + copt259 + copt260 + copt261 + copt262 + copt263 +
                 copt264 + copt265 + copt266 + copt267 + copt268 + copt270 +
                 copt271 + copt272 + copt273 + copt274 + copt275 + copt276 +
                 copt277 + copt278 + copt279 + copt280 + copt281 + copt283 +
                 copt284 + copt285 + copt286 + copt287 + copt288 + copt289 +
                 copt290 + copt291 + copt292 + copt293 + copt294 + copt295 +
                 copt296 + copt325 + copt349;
  Real copt351 = copt235 * copt350;
  Real copt352 = -2 * copt4 * copt8;
  Real copt353 = -2 * copt15 * copt18;
  Real copt354 = -2 * copt25 * copt28;
  Real copt355 = copt103 + copt237 + copt243 + copt253 + copt352 + copt353 +
                 copt354 + copt77 + copt79;
  Real copt356 = copt18 * copt25 * copt256;
  Real copt357 = -(copt15 * copt256 * copt28);
  Real copt358 = copt2 * copt25 * copt269;
  Real copt359 = -(copt25 * copt269 * copt8);
  Real copt360 = -(copt2 * copt269 * copt28);
  Real copt361 = copt269 * copt28 * copt4;
  Real copt362 = -(copt18 * copt256);
  Real copt363 = copt199 + copt256;
  Real copt364 = copt15 * copt363;
  Real copt365 = copt269 * copt8;
  Real copt366 = copt335 + copt362 + copt364 + copt365;
  Real copt367 = copt23 * copt366;
  Real copt368 = -(copt15 * copt2 * copt282);
  Real copt369 = copt15 * copt282 * copt8;
  Real copt370 = copt18 * copt2 * copt282;
  Real copt371 = -(copt18 * copt282 * copt4);
  Real copt372 = copt25 * copt301;
  Real copt373 = copt256 * copt28;
  Real copt374 = -(copt282 * copt8);
  Real copt375 = copt215 + copt282;
  Real copt376 = copt375 * copt4;
  Real copt377 = copt372 + copt373 + copt374 + copt376;
  Real copt378 = copt13 * copt377;
  Real copt379 = copt191 + copt192 + copt356 + copt357 + copt358 + copt359 +
                 copt360 + copt361 + copt367 + copt368 + copt369 + copt370 +
                 copt371 + copt378;
  Real copt380 = copt355 * copt379;
  Real copt381 = ArcTan(copt351, copt380);
  Real copt531 = t1(1);
  Real copt394 = copt15 * copt18 * copt76;
  Real copt395 = -(copt15 * copt18 * copt2 * copt8);
  Real copt396 = -(copt243 * copt76);
  Real copt397 = copt2 * copt243 * copt4;
  Real copt398 = copt25 * copt28 * copt76;
  Real copt399 = -(copt2 * copt25 * copt28 * copt8);
  Real copt400 = -(copt253 * copt76);
  Real copt401 = copt2 * copt253 * copt4;
  Real copt403 = -(copt15 * copt18 * copt2 * copt402);
  Real copt404 = copt15 * copt18 * copt402 * copt8;
  Real copt405 = copt2 * copt243 * copt402;
  Real copt406 = -(copt243 * copt4 * copt402);
  Real copt407 = -(copt2 * copt25 * copt28 * copt402);
  Real copt408 = copt25 * copt28 * copt402 * copt8;
  Real copt409 = copt2 * copt253 * copt402;
  Real copt410 = -(copt253 * copt4 * copt402);
  Real copt419 = copt230 * copt418;
  Real copt420 = copt411 + copt414 + copt415 + copt419;
  Real copt421 = copt107 * copt420;
  Real copt422 = -(copt15 * copt416 * copt76);
  Real copt423 = 2 * copt15 * copt2 * copt416 * copt8;
  Real copt424 = -(copt15 * copt237 * copt416);
  Real copt425 = copt18 * copt416 * copt76;
  Real copt426 = -(copt18 * copt2 * copt4 * copt416);
  Real copt427 = -(copt18 * copt2 * copt416 * copt8);
  Real copt428 = copt18 * copt4 * copt416 * copt8;
  Real copt429 = copt18 * copt25 * copt28 * copt416;
  Real copt430 = -(copt15 * copt253 * copt416);
  Real copt434 = copt232 * copt433;
  Real copt435 = copt411 + copt414 + copt415 + copt434;
  Real copt436 = copt126 * copt435;
  Real copt437 = -(copt25 * copt431 * copt76);
  Real copt438 = 2 * copt2 * copt25 * copt431 * copt8;
  Real copt439 = -(copt237 * copt25 * copt431);
  Real copt440 = -(copt243 * copt25 * copt431);
  Real copt441 = copt28 * copt431 * copt76;
  Real copt442 = -(copt2 * copt28 * copt4 * copt431);
  Real copt443 = -(copt2 * copt28 * copt431 * copt8);
  Real copt444 = copt28 * copt4 * copt431 * copt8;
  Real copt445 = copt15 * copt18 * copt28 * copt431;
  Real copt446 = -(copt15 * copt237);
  Real copt447 = copt18 * copt23 * copt25;
  Real copt448 = copt18 * copt4 * copt8;
  Real copt449 = -2 * copt18 * copt23 * copt28;
  Real copt450 = copt18 * copt25 * copt28;
  Real copt451 = copt15 * copt402 * copt8;
  Real copt452 = -2 * copt18 * copt4 * copt402;
  Real copt453 = copt18 * copt402 * copt8;
  Real copt454 = -(copt23 * copt25 * copt416);
  Real copt455 = copt4 * copt416 * copt8;
  Real copt456 = -(copt237 * copt416);
  Real copt457 = copt23 * copt28 * copt416;
  Real copt458 = copt25 * copt28 * copt416;
  Real copt459 = -(copt253 * copt416);
  Real copt460 = -2 * copt18 * copt8;
  Real copt461 = copt15 * copt413;
  Real copt462 = copt18 * copt402;
  Real copt465 = copt460 + copt461 + copt462 + copt463 + copt464;
  Real copt466 = copt2 * copt465;
  Real copt467 = copt15 * copt390 * copt433;
  Real copt468 = copt18 * copt23 * copt431;
  Real copt469 = -2 * copt18 * copt25 * copt431;
  Real copt470 = copt18 * copt28 * copt431;
  Real copt471 = copt446 + copt447 + copt448 + copt449 + copt450 + copt451 +
                 copt452 + copt453 + copt454 + copt455 + copt456 + copt457 +
                 copt458 + copt459 + copt466 + copt467 + copt468 + copt469 +
                 copt470;
  Real copt472 = -(copt13 * copt471);
  Real copt473 = copt28 * copt4 * copt8;
  Real copt474 = copt15 * copt18 * copt28;
  Real copt475 = -2 * copt28 * copt4 * copt402;
  Real copt476 = copt28 * copt402 * copt8;
  Real copt477 = -2 * copt15 * copt28 * copt416;
  Real copt478 = copt18 * copt28 * copt416;
  Real copt479 = -(copt402 * copt8);
  Real copt480 = -(copt18 * copt416);
  Real copt481 = copt237 + copt243 + copt479 + copt480;
  Real copt482 = -(copt25 * copt481);
  Real copt483 = copt4 * copt431 * copt8;
  Real copt484 = -(copt237 * copt431);
  Real copt485 = copt15 * copt18 * copt431;
  Real copt486 = -(copt243 * copt431);
  Real copt487 = -2 * copt28 * copt8;
  Real copt490 = copt4 * copt433;
  Real copt491 = copt431 * copt8;
  Real copt492 = copt487 + copt488 + copt489 + copt490 + copt491;
  Real copt493 = copt2 * copt492;
  Real copt494 = copt473 + copt474 + copt475 + copt476 + copt477 + copt478 +
                 copt482 + copt483 + copt484 + copt485 + copt486 + copt493;
  Real copt495 = -(copt23 * copt494);
  Real copt496 = copt394 + copt395 + copt396 + copt397 + copt398 + copt399 +
                 copt400 + copt401 + copt403 + copt404 + copt405 + copt406 +
                 copt407 + copt408 + copt409 + copt410 + copt421 + copt422 +
                 copt423 + copt424 + copt425 + copt426 + copt427 + copt428 +
                 copt429 + copt430 + copt436 + copt437 + copt438 + copt439 +
                 copt440 + copt441 + copt442 + copt443 + copt444 + copt445 +
                 copt472 + copt495;
  Real copt497 = copt393 * copt496;
  Real copt498 = -2 * copt2 * copt8;
  Real copt499 = -2 * copt13 * copt18;
  Real copt500 = -2 * copt23 * copt28;
  Real copt501 = copt107 + copt126 + copt237 + copt243 + copt253 + copt498 +
                 copt499 + copt500 + copt76;
  Real copt502 = copt18 * copt25 * copt402;
  Real copt503 = -(copt15 * copt28 * copt402);
  Real copt504 = copt2 * copt25 * copt416;
  Real copt505 = -(copt25 * copt416 * copt8);
  Real copt506 = -(copt2 * copt28 * copt416);
  Real copt507 = copt28 * copt4 * copt416;
  Real copt508 = -(copt18 * copt402);
  Real copt509 = copt199 + copt402;
  Real copt510 = copt15 * copt509;
  Real copt511 = copt463 + copt464 + copt508 + copt510;
  Real copt512 = copt23 * copt511;
  Real copt513 = -(copt15 * copt2 * copt431);
  Real copt514 = copt15 * copt431 * copt8;
  Real copt515 = copt18 * copt2 * copt431;
  Real copt516 = -(copt18 * copt4 * copt431);
  Real copt517 = -(copt431 * copt8);
  Real copt518 = copt215 + copt431;
  Real copt519 = copt4 * copt518;
  Real copt520 = copt488 + copt489 + copt517 + copt519;
  Real copt521 = copt13 * copt520;
  Real copt522 = copt191 + copt192 + copt502 + copt503 + copt504 + copt505 +
                 copt506 + copt507 + copt512 + copt513 + copt514 + copt515 +
                 copt516 + copt521;
  Real copt523 = copt501 * copt522;
  Real copt524 = ArcTan(copt497, copt523);
  Real copt534 = t2(1);
  Real copt539 = Power(copt528, 2);
  Real copt542 = Power(copt531, 2);
  Real copt545 = Power(copt534, 2);
  out(0)       = copt34;
  out(1) =
      copt35 * (copt11 * copt40 + copt21 * copt44 + copt31 * copt48) * copt56;
  out(2) = copt55;
  out(3) =
      -(copt58 * copt59 * copt60 * copt61 *
        (copt384 * copt524 * l0 * l1 + copt225 * copt381 * l0 * l2 +
         copt222 * copt63 * l1 * l2 + copt63 * l1 * l2 * thetarest0 +
         copt225 * l0 * l2 * thetarest1 + copt384 * l0 * l1 * thetarest2)) /
      2.;
  out(4) = -(copt58 * copt59 * copt60 * copt61 *
             (copt383 * copt524 * copt534 * l0 * l1 +
              copt224 * copt381 * copt531 * l0 * l2 +
              copt222 * copt528 * copt62 * l1 * l2 +
              copt528 * copt62 * l1 * l2 * thetarest0 +
              copt224 * copt531 * l0 * l2 * thetarest1 +
              copt383 * copt534 * l0 * l1 * thetarest2)) /
           2.;
  out(5) =
      -(copt58 * copt59 * copt60 * copt61 *
        (copt524 * copt545 * l0 * l1 + copt381 * copt542 * l0 * l2 +
         copt222 * copt539 * l1 * l2 + copt539 * l1 * l2 * thetarest0 +
         copt542 * l0 * l2 * thetarest1 + copt545 * l0 * l1 * thetarest2)) /
      2.;
  return out;
}
