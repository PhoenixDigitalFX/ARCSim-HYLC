#include "strain.hpp"
#ifndef hylc_strain_II

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
  Real copt74  = xloc(9);
  Real copt79  = xloc(10);
  Real copt66  = -copt4;
  Real copt98  = xloc(11);
  Real copt81  = -copt15;
  Real copt91  = -copt25;
  Real copt123 = copt2 + copt66;
  Real copt124 = Power(copt123, 2);
  Real copt125 = copt13 + copt81;
  Real copt126 = Power(copt125, 2);
  Real copt127 = copt23 + copt91;
  Real copt128 = Power(copt127, 2);
  Real copt129 = copt124 + copt126 + copt128;
  Real copt130 = Sqrt(copt129);
  Real copt86  = -copt8;
  Real copt113 = -copt79;
  Real copt75  = -copt74;
  Real copt108 = -copt28;
  Real copt163 = t1(0);
  Real copt164 = Power(copt163, 2);
  Real copt87  = copt4 + copt86;
  Real copt106 = copt18 + copt81;
  Real copt171 = xloc(12);
  Real copt175 = xloc(13);
  Real copt88  = copt23 * copt87;
  Real copt89  = copt25 * copt8;
  Real copt90  = -(copt28 * copt4);
  Real copt92  = copt28 + copt91;
  Real copt93  = copt2 * copt92;
  Real copt94  = copt88 + copt89 + copt90 + copt93;
  Real copt173 = copt171 + copt86;
  Real copt184 = xloc(14);
  Real copt69  = -copt18;
  Real copt70  = copt15 + copt69;
  Real copt185 = -copt184;
  Real copt186 = copt185 + copt28;
  Real copt204 = Power(copt87, 2);
  Real copt205 = Power(copt70, 2);
  Real copt109 = copt108 + copt25;
  Real copt206 = Power(copt109, 2);
  Real copt207 = copt204 + copt205 + copt206;
  Real copt208 = Sqrt(copt207);
  Real copt131 = -(copt18 * copt2 * copt25);
  Real copt132 = copt15 * copt2 * copt28;
  Real copt172 = -(copt171 * copt18);
  Real copt174 = copt15 * copt173;
  Real copt176 = -copt175;
  Real copt177 = copt176 + copt18;
  Real copt178 = copt177 * copt4;
  Real copt179 = copt175 * copt8;
  Real copt180 = copt172 + copt174 + copt178 + copt179;
  Real copt233 = t2(0);
  Real copt234 = Power(copt233, 2);
  Real copt166 = copt13 * copt87;
  Real copt167 = copt15 * copt8;
  Real copt168 = -(copt18 * copt4);
  Real copt169 = copt106 * copt2;
  Real copt170 = copt166 + copt167 + copt168 + copt169;
  Real copt236 = xloc(15);
  Real copt241 = xloc(16);
  Real copt237 = -copt236;
  Real copt238 = copt237 + copt8;
  Real copt249 = xloc(17);
  Real copt191 = copt23 * copt70;
  Real copt192 = copt18 * copt25;
  Real copt193 = -(copt15 * copt28);
  Real copt194 = copt13 * copt92;
  Real copt195 = copt191 + copt192 + copt193 + copt194;
  Real copt251 = copt108 + copt249;
  Real copt264 = copt2 + copt86;
  Real copt265 = Power(copt264, 2);
  Real copt266 = copt13 + copt69;
  Real copt267 = Power(copt266, 2);
  Real copt268 = copt108 + copt23;
  Real copt269 = Power(copt268, 2);
  Real copt270 = copt265 + copt267 + copt269;
  Real copt271 = Sqrt(copt270);
  Real copt240 = copt18 * copt236;
  Real copt242 = -(copt241 * copt8);
  Real copt243 = copt241 + copt69;
  Real copt65  = -(copt15 * copt8);
  Real copt67  = copt66 + copt8;
  Real copt68  = copt13 * copt67;
  Real copt71  = copt2 * copt70;
  Real copt72  = copt18 * copt4;
  Real copt73  = copt65 + copt68 + copt71 + copt72;
  Real copt76  = copt4 + copt75;
  Real copt77  = copt13 * copt76;
  Real copt78  = copt15 * copt74;
  Real copt80  = -(copt4 * copt79);
  Real copt82  = copt79 + copt81;
  Real copt83  = copt2 * copt82;
  Real copt84  = copt77 + copt78 + copt80 + copt83;
  Real copt85  = copt73 * copt84;
  Real copt95  = -(copt25 * copt74);
  Real copt96  = copt66 + copt74;
  Real copt97  = copt23 * copt96;
  Real copt99  = -copt98;
  Real copt100 = copt25 + copt99;
  Real copt101 = copt100 * copt2;
  Real copt102 = copt4 * copt98;
  Real copt103 = copt101 + copt102 + copt95 + copt97;
  Real copt104 = copt103 * copt94;
  Real copt105 = -(copt18 * copt25);
  Real copt107 = copt106 * copt23;
  Real copt110 = copt109 * copt13;
  Real copt111 = copt15 * copt28;
  Real copt112 = copt105 + copt107 + copt110 + copt111;
  Real copt114 = copt113 + copt15;
  Real copt115 = copt114 * copt23;
  Real copt116 = copt25 * copt79;
  Real copt117 = -(copt15 * copt98);
  Real copt118 = copt91 + copt98;
  Real copt119 = copt118 * copt13;
  Real copt120 = copt115 + copt116 + copt117 + copt119;
  Real copt121 = copt112 * copt120;
  Real copt122 = copt104 + copt121 + copt85;
  Real copt133 = copt18 * copt25 * copt74;
  Real copt134 = -(copt15 * copt28 * copt74);
  Real copt135 = copt2 * copt25 * copt79;
  Real copt136 = -(copt25 * copt79 * copt8);
  Real copt137 = -(copt2 * copt28 * copt79);
  Real copt138 = copt28 * copt4 * copt79;
  Real copt139 = -(copt18 * copt74);
  Real copt140 = copt74 + copt86;
  Real copt141 = copt140 * copt15;
  Real copt142 = copt113 + copt18;
  Real copt143 = copt142 * copt4;
  Real copt144 = copt79 * copt8;
  Real copt145 = copt139 + copt141 + copt143 + copt144;
  Real copt146 = copt145 * copt23;
  Real copt147 = -(copt15 * copt2 * copt98);
  Real copt148 = copt15 * copt8 * copt98;
  Real copt149 = copt18 * copt2 * copt98;
  Real copt150 = -(copt18 * copt4 * copt98);
  Real copt151 = copt75 + copt8;
  Real copt152 = copt151 * copt25;
  Real copt153 = copt28 * copt74;
  Real copt154 = -(copt8 * copt98);
  Real copt155 = copt108 + copt98;
  Real copt156 = copt155 * copt4;
  Real copt157 = copt152 + copt153 + copt154 + copt156;
  Real copt158 = copt13 * copt157;
  Real copt159 = copt131 + copt132 + copt133 + copt134 + copt135 + copt136 +
                 copt137 + copt138 + copt146 + copt147 + copt148 + copt149 +
                 copt150 + copt158;
  Real copt160 = copt130 * copt159;
  Real copt161 = ArcTan(copt122, copt160);
  Real copt303 = t0(1);
  Real copt181 = copt170 * copt180;
  Real copt182 = -(copt171 * copt28);
  Real copt183 = copt173 * copt25;
  Real copt187 = copt186 * copt4;
  Real copt188 = copt184 * copt8;
  Real copt189 = copt182 + copt183 + copt187 + copt188;
  Real copt190 = copt189 * copt94;
  Real copt196 = -(copt175 * copt28);
  Real copt197 = copt175 + copt69;
  Real copt198 = copt197 * copt25;
  Real copt199 = copt15 * copt186;
  Real copt200 = copt18 * copt184;
  Real copt201 = copt196 + copt198 + copt199 + copt200;
  Real copt202 = copt195 * copt201;
  Real copt203 = copt181 + copt190 + copt202;
  Real copt209 = copt171 * copt18 * copt25;
  Real copt210 = -(copt15 * copt171 * copt28);
  Real copt211 = copt175 * copt2 * copt25;
  Real copt212 = -(copt175 * copt25 * copt8);
  Real copt213 = -(copt175 * copt2 * copt28);
  Real copt214 = copt175 * copt28 * copt4;
  Real copt215 = copt180 * copt23;
  Real copt216 = -(copt15 * copt184 * copt2);
  Real copt217 = copt15 * copt184 * copt8;
  Real copt218 = copt18 * copt184 * copt2;
  Real copt219 = -(copt18 * copt184 * copt4);
  Real copt220 = -copt171;
  Real copt221 = copt220 + copt8;
  Real copt222 = copt221 * copt25;
  Real copt223 = copt171 * copt28;
  Real copt224 = -(copt184 * copt8);
  Real copt225 = copt108 + copt184;
  Real copt226 = copt225 * copt4;
  Real copt227 = copt222 + copt223 + copt224 + copt226;
  Real copt228 = copt13 * copt227;
  Real copt229 = copt131 + copt132 + copt209 + copt210 + copt211 + copt212 +
                 copt213 + copt214 + copt215 + copt216 + copt217 + copt218 +
                 copt219 + copt228;
  Real copt230 = copt208 * copt229;
  Real copt231 = ArcTan(copt203, copt230);
  Real copt306 = t1(1);
  Real copt239 = copt13 * copt238;
  Real copt244 = copt2 * copt243;
  Real copt245 = copt239 + copt240 + copt242 + copt244;
  Real copt246 = copt170 * copt245;
  Real copt247 = copt23 * copt238;
  Real copt248 = copt236 * copt28;
  Real copt250 = -(copt249 * copt8);
  Real copt252 = copt2 * copt251;
  Real copt253 = copt247 + copt248 + copt250 + copt252;
  Real copt254 = copt253 * copt94;
  Real copt255 = -copt241;
  Real copt256 = copt18 + copt255;
  Real copt257 = copt23 * copt256;
  Real copt258 = copt241 * copt28;
  Real copt259 = -(copt18 * copt249);
  Real copt260 = copt13 * copt251;
  Real copt261 = copt257 + copt258 + copt259 + copt260;
  Real copt262 = copt195 * copt261;
  Real copt263 = copt246 + copt254 + copt262;
  Real copt272 = copt18 * copt2 * copt25;
  Real copt273 = -(copt15 * copt2 * copt28);
  Real copt274 = -(copt18 * copt236 * copt25);
  Real copt275 = copt15 * copt236 * copt28;
  Real copt276 = -(copt2 * copt241 * copt25);
  Real copt277 = copt241 * copt25 * copt8;
  Real copt278 = copt2 * copt241 * copt28;
  Real copt279 = -(copt241 * copt28 * copt4);
  Real copt280 = copt15 * copt238;
  Real copt281 = copt243 * copt4;
  Real copt282 = copt240 + copt242 + copt280 + copt281;
  Real copt283 = copt23 * copt282;
  Real copt284 = copt15 * copt2 * copt249;
  Real copt285 = -(copt15 * copt249 * copt8);
  Real copt286 = -(copt18 * copt2 * copt249);
  Real copt287 = copt18 * copt249 * copt4;
  Real copt288 = -(copt236 * copt28);
  Real copt289 = copt236 + copt86;
  Real copt290 = copt25 * copt289;
  Real copt291 = -copt249;
  Real copt292 = copt28 + copt291;
  Real copt293 = copt292 * copt4;
  Real copt294 = copt249 * copt8;
  Real copt295 = copt288 + copt290 + copt293 + copt294;
  Real copt296 = copt13 * copt295;
  Real copt297 = copt272 + copt273 + copt274 + copt275 + copt276 + copt277 +
                 copt278 + copt279 + copt283 + copt284 + copt285 + copt286 +
                 copt287 + copt296;
  Real copt298 = -(copt271 * copt297);
  Real copt299 = ArcTan(copt263, copt298);
  Real copt309 = t2(1);
  Real copt314 = Power(copt303, 2);
  Real copt317 = Power(copt306, 2);
  Real copt320 = Power(copt309, 2);
  out(0)       = copt34;
  out(1) =
      copt35 * (copt11 * copt40 + copt21 * copt44 + copt31 * copt48) * copt56;
  out(2) = copt55;
  out(3) =
      -(copt58 * copt59 * copt60 * copt61 *
        (copt234 * copt299 * l0 * l1 + copt164 * copt231 * l0 * l2 +
         copt161 * copt63 * l1 * l2 + copt63 * l1 * l2 * thetarest0 +
         copt164 * l0 * l2 * thetarest1 + copt234 * l0 * l1 * thetarest2)) /
      2.;
  out(4) = -(copt58 * copt59 * copt60 * copt61 *
             (copt233 * copt299 * copt309 * l0 * l1 +
              copt163 * copt231 * copt306 * l0 * l2 +
              copt161 * copt303 * copt62 * l1 * l2 +
              copt303 * copt62 * l1 * l2 * thetarest0 +
              copt163 * copt306 * l0 * l2 * thetarest1 +
              copt233 * copt309 * l0 * l1 * thetarest2)) /
           2.;
  out(5) =
      -(copt58 * copt59 * copt60 * copt61 *
        (copt299 * copt320 * l0 * l1 + copt231 * copt317 * l0 * l2 +
         copt161 * copt314 * l1 * l2 + copt314 * l1 * l2 * thetarest0 +
         copt317 * l0 * l2 * thetarest1 + copt320 * l0 * l1 * thetarest2)) /
      2.;
  return out;
}

#endif  // hylc_strain_II
