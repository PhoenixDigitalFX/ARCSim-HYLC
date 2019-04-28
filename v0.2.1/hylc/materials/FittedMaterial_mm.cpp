#include "FittedMaterial.hpp"

using namespace hylc;
using namespace hylc::mathematica;

typedef Vec6 strain_type;
typedef double value_type;
typedef Vec6 grad_type;
typedef std::pair<Mat6x6, Vec6> gradhess_type;

// 0th derivative of psi, grad, tuple<grad,hess>
value_type FittedMaterial::psi_taylor_0(const strain_type &straincpp) {
  Real copt1   = straincpp(0);
  Real copt2   = straincpp(2);
  Real copt3   = copt1 * copt2;
  Real copt4   = Ccompr1 + copt3;
  Real copt5   = 1 / copt4;
  Real copt7   = -1 + copt1;
  Real copt8   = Power(copt7, 2);
  Real copt9   = Power(copt8, 2);
  Real copt10  = strainscale(0);
  Real copt11  = Power(copt10, 2);
  Real copt12  = Power(copt11, 2);
  Real copt13  = 1 / copt12;
  Real copt15  = copt7 * copt8;
  Real copt16  = copt10 * copt11;
  Real copt17  = 1 / copt16;
  Real copt19  = 1 / copt11;
  Real copt21  = 1 / copt10;
  Real copt23  = straincpp(1);
  Real copt24  = Power(copt23, 2);
  Real copt25  = Power(copt24, 2);
  Real copt26  = strainscale(1);
  Real copt27  = Power(copt26, 2);
  Real copt28  = Power(copt27, 2);
  Real copt29  = 1 / copt28;
  Real copt31  = copt23 * copt24;
  Real copt32  = copt26 * copt27;
  Real copt33  = 1 / copt32;
  Real copt36  = 1 / copt27;
  Real copt40  = 1 / copt26;
  Real copt45  = -1 + copt2;
  Real copt46  = Power(copt45, 2);
  Real copt47  = Power(copt46, 2);
  Real copt48  = strainscale(2);
  Real copt49  = Power(copt48, 2);
  Real copt50  = Power(copt49, 2);
  Real copt51  = 1 / copt50;
  Real copt53  = copt45 * copt46;
  Real copt54  = copt48 * copt49;
  Real copt55  = 1 / copt54;
  Real copt59  = 1 / copt49;
  Real copt65  = 1 / copt48;
  Real copt73  = straincpp(3);
  Real copt74  = Power(copt73, 2);
  Real copt75  = Power(copt74, 2);
  Real copt76  = strainscale(3);
  Real copt77  = Power(copt76, 2);
  Real copt78  = Power(copt77, 2);
  Real copt79  = 1 / copt78;
  Real copt81  = copt73 * copt74;
  Real copt82  = copt76 * copt77;
  Real copt83  = 1 / copt82;
  Real copt88  = 1 / copt77;
  Real copt96  = 1 / copt76;
  Real copt107 = straincpp(4);
  Real copt108 = Power(copt107, 2);
  Real copt109 = Power(copt108, 2);
  Real copt110 = strainscale(4);
  Real copt111 = Power(copt110, 2);
  Real copt112 = Power(copt111, 2);
  Real copt113 = 1 / copt112;
  Real copt115 = copt107 * copt108;
  Real copt116 = copt110 * copt111;
  Real copt117 = 1 / copt116;
  Real copt123 = 1 / copt111;
  Real copt133 = 1 / copt110;
  Real copt147 = straincpp(5);
  Real copt148 = Power(copt147, 2);
  Real copt149 = Power(copt148, 2);
  Real copt150 = strainscale(5);
  Real copt151 = Power(copt150, 2);
  Real copt152 = Power(copt151, 2);
  Real copt153 = 1 / copt152;
  Real copt155 = copt147 * copt148;
  Real copt156 = copt150 * copt151;
  Real copt157 = 1 / copt156;
  Real copt164 = 1 / copt151;
  Real copt176 = 1 / copt150;
  return C0 + C44 * copt109 * copt113 + C43 * copt115 * copt117 +
         C42 * copt108 * copt123 + C41 * copt107 * copt133 +
         C54 * copt149 * copt153 + C53 * copt155 * copt157 +
         C4513 * copt107 * copt133 * copt155 * copt157 +
         C52 * copt148 * copt164 +
         C4522 * copt108 * copt123 * copt148 * copt164 +
         C4512 * copt107 * copt133 * copt148 * copt164 + C03 * copt15 * copt17 +
         C0431 * copt107 * copt133 * copt15 * copt17 + C51 * copt147 * copt176 +
         C4531 * copt115 * copt117 * copt147 * copt176 +
         C4521 * copt108 * copt123 * copt147 * copt176 +
         C4511 * copt107 * copt133 * copt147 * copt176 +
         C0531 * copt147 * copt15 * copt17 * copt176 + C14 * copt25 * copt29 +
         C13 * copt31 * copt33 + C1431 * copt107 * copt133 * copt31 * copt33 +
         C1531 * copt147 * copt176 * copt31 * copt33 + C12 * copt24 * copt36 +
         C1422 * copt108 * copt123 * copt24 * copt36 +
         C1421 * copt107 * copt133 * copt24 * copt36 +
         C1522 * copt148 * copt164 * copt24 * copt36 +
         C1521 * copt147 * copt176 * copt24 * copt36 + C11 * copt23 * copt40 +
         C1413 * copt115 * copt117 * copt23 * copt40 +
         C1412 * copt108 * copt123 * copt23 * copt40 +
         C1411 * copt107 * copt133 * copt23 * copt40 +
         C1513 * copt155 * copt157 * copt23 * copt40 +
         C1512 * copt148 * copt164 * copt23 * copt40 +
         C0131 * copt15 * copt17 * copt23 * copt40 +
         C1511 * copt147 * copt176 * copt23 * copt40 + Ccompr0 * copt5 +
         C24 * copt47 * copt51 + C23 * copt53 * copt55 +
         C2431 * copt107 * copt133 * copt53 * copt55 +
         C2531 * copt147 * copt176 * copt53 * copt55 +
         C1213 * copt23 * copt40 * copt53 * copt55 + C22 * copt46 * copt59 +
         C2422 * copt108 * copt123 * copt46 * copt59 +
         C2421 * copt107 * copt133 * copt46 * copt59 +
         C2522 * copt148 * copt164 * copt46 * copt59 +
         C2521 * copt147 * copt176 * copt46 * copt59 +
         C1222 * copt24 * copt36 * copt46 * copt59 +
         C1212 * copt23 * copt40 * copt46 * copt59 + C21 * copt45 * copt65 +
         C2413 * copt115 * copt117 * copt45 * copt65 +
         C2412 * copt108 * copt123 * copt45 * copt65 +
         C2411 * copt107 * copt133 * copt45 * copt65 +
         C2513 * copt155 * copt157 * copt45 * copt65 +
         C2512 * copt148 * copt164 * copt45 * copt65 +
         C0231 * copt15 * copt17 * copt45 * copt65 +
         C2511 * copt147 * copt176 * copt45 * copt65 +
         C1231 * copt31 * copt33 * copt45 * copt65 +
         C1221 * copt24 * copt36 * copt45 * copt65 +
         C1211 * copt23 * copt40 * copt45 * copt65 + C01 * copt21 * copt7 +
         C0413 * copt115 * copt117 * copt21 * copt7 +
         C0412 * copt108 * copt123 * copt21 * copt7 +
         C0411 * copt107 * copt133 * copt21 * copt7 +
         C0513 * copt155 * copt157 * copt21 * copt7 +
         C0512 * copt148 * copt164 * copt21 * copt7 +
         C0511 * copt147 * copt176 * copt21 * copt7 +
         C0113 * copt21 * copt31 * copt33 * copt7 +
         C0112 * copt21 * copt24 * copt36 * copt7 +
         C0111 * copt21 * copt23 * copt40 * copt7 +
         C0213 * copt21 * copt53 * copt55 * copt7 +
         C0212 * copt21 * copt46 * copt59 * copt7 +
         C0211 * copt21 * copt45 * copt65 * copt7 + C34 * copt75 * copt79 +
         C02 * copt19 * copt8 + C0422 * copt108 * copt123 * copt19 * copt8 +
         C0421 * copt107 * copt133 * copt19 * copt8 +
         C0522 * copt148 * copt164 * copt19 * copt8 +
         C0521 * copt147 * copt176 * copt19 * copt8 +
         C0122 * copt19 * copt24 * copt36 * copt8 +
         C0121 * copt19 * copt23 * copt40 * copt8 +
         C0222 * copt19 * copt46 * copt59 * copt8 +
         C0221 * copt19 * copt45 * copt65 * copt8 + C33 * copt81 * copt83 +
         C3431 * copt107 * copt133 * copt81 * copt83 +
         C3531 * copt147 * copt176 * copt81 * copt83 +
         C1313 * copt23 * copt40 * copt81 * copt83 +
         C2313 * copt45 * copt65 * copt81 * copt83 +
         C0313 * copt21 * copt7 * copt81 * copt83 + C32 * copt74 * copt88 +
         C3422 * copt108 * copt123 * copt74 * copt88 +
         C3421 * copt107 * copt133 * copt74 * copt88 +
         C3522 * copt148 * copt164 * copt74 * copt88 +
         C3521 * copt147 * copt176 * copt74 * copt88 +
         C1322 * copt24 * copt36 * copt74 * copt88 +
         C1312 * copt23 * copt40 * copt74 * copt88 +
         C2322 * copt46 * copt59 * copt74 * copt88 +
         C2312 * copt45 * copt65 * copt74 * copt88 +
         C0312 * copt21 * copt7 * copt74 * copt88 +
         C0322 * copt19 * copt74 * copt8 * copt88 + C04 * copt13 * copt9 +
         C31 * copt73 * copt96 + C3413 * copt115 * copt117 * copt73 * copt96 +
         C3412 * copt108 * copt123 * copt73 * copt96 +
         C3411 * copt107 * copt133 * copt73 * copt96 +
         C3513 * copt155 * copt157 * copt73 * copt96 +
         C3512 * copt148 * copt164 * copt73 * copt96 +
         C0331 * copt15 * copt17 * copt73 * copt96 +
         C3511 * copt147 * copt176 * copt73 * copt96 +
         C1331 * copt31 * copt33 * copt73 * copt96 +
         C1321 * copt24 * copt36 * copt73 * copt96 +
         C1311 * copt23 * copt40 * copt73 * copt96 +
         C2331 * copt53 * copt55 * copt73 * copt96 +
         C2321 * copt46 * copt59 * copt73 * copt96 +
         C2311 * copt45 * copt65 * copt73 * copt96 +
         C0311 * copt21 * copt7 * copt73 * copt96 +
         C0321 * copt19 * copt73 * copt8 * copt96;
}

grad_type FittedMaterial::grad_taylor_0(const strain_type &straincpp) {
  Vec6 out(0);
  Real copt2   = straincpp(0);
  Real copt1   = straincpp(2);
  Real copt3   = copt1 * copt2;
  Real copt4   = Ccompr1 + copt3;
  Real copt5   = Power(copt4, 2);
  Real copt6   = 1 / copt5;
  Real copt8   = -1 + copt2;
  Real copt9   = Power(copt8, 2);
  Real copt10  = copt8 * copt9;
  Real copt11  = strainscale(0);
  Real copt12  = Power(copt11, 2);
  Real copt13  = Power(copt12, 2);
  Real copt14  = 1 / copt13;
  Real copt16  = copt11 * copt12;
  Real copt17  = 1 / copt16;
  Real copt19  = 1 / copt12;
  Real copt21  = 1 / copt11;
  Real copt23  = straincpp(1);
  Real copt24  = Power(copt23, 2);
  Real copt25  = copt23 * copt24;
  Real copt26  = strainscale(1);
  Real copt27  = Power(copt26, 2);
  Real copt28  = copt26 * copt27;
  Real copt29  = 1 / copt28;
  Real copt31  = 1 / copt27;
  Real copt34  = 1 / copt26;
  Real copt38  = -1 + copt1;
  Real copt39  = Power(copt38, 2);
  Real copt40  = copt38 * copt39;
  Real copt41  = strainscale(2);
  Real copt42  = Power(copt41, 2);
  Real copt43  = copt41 * copt42;
  Real copt44  = 1 / copt43;
  Real copt46  = 1 / copt42;
  Real copt49  = 1 / copt41;
  Real copt53  = straincpp(3);
  Real copt54  = Power(copt53, 2);
  Real copt55  = copt53 * copt54;
  Real copt56  = strainscale(3);
  Real copt57  = Power(copt56, 2);
  Real copt58  = copt56 * copt57;
  Real copt59  = 1 / copt58;
  Real copt61  = 1 / copt57;
  Real copt64  = 1 / copt56;
  Real copt68  = straincpp(4);
  Real copt69  = Power(copt68, 2);
  Real copt70  = copt68 * copt69;
  Real copt71  = strainscale(4);
  Real copt72  = Power(copt71, 2);
  Real copt73  = copt71 * copt72;
  Real copt74  = 1 / copt73;
  Real copt76  = 1 / copt72;
  Real copt79  = 1 / copt71;
  Real copt83  = straincpp(5);
  Real copt84  = Power(copt83, 2);
  Real copt85  = copt83 * copt84;
  Real copt86  = strainscale(5);
  Real copt87  = Power(copt86, 2);
  Real copt88  = copt86 * copt87;
  Real copt89  = 1 / copt88;
  Real copt91  = 1 / copt87;
  Real copt94  = 1 / copt86;
  Real copt99  = Power(copt27, 2);
  Real copt100 = 1 / copt99;
  Real copt138 = Power(copt42, 2);
  Real copt139 = 1 / copt138;
  Real copt175 = Power(copt57, 2);
  Real copt176 = 1 / copt175;
  Real copt213 = Power(copt72, 2);
  Real copt214 = 1 / copt213;
  Real copt251 = Power(copt87, 2);
  Real copt252 = 1 / copt251;
  out(0)       = 4 * C04 * copt10 * copt14 + C01 * copt21 +
           C0113 * copt21 * copt25 * copt29 + C0112 * copt21 * copt24 * copt31 +
           C0111 * copt21 * copt23 * copt34 + C0213 * copt21 * copt40 * copt44 +
           C0212 * copt21 * copt39 * copt46 + C0211 * copt21 * copt38 * copt49 +
           C0313 * copt21 * copt55 * copt59 - Ccompr0 * copt1 * copt6 +
           C0312 * copt21 * copt54 * copt61 + C0311 * copt21 * copt53 * copt64 +
           C0413 * copt21 * copt70 * copt74 + C0412 * copt21 * copt69 * copt76 +
           C0411 * copt21 * copt68 * copt79 + 2 * C02 * copt19 * copt8 +
           2 * C0122 * copt19 * copt24 * copt31 * copt8 +
           2 * C0121 * copt19 * copt23 * copt34 * copt8 +
           2 * C0222 * copt19 * copt39 * copt46 * copt8 +
           2 * C0221 * copt19 * copt38 * copt49 * copt8 +
           2 * C0322 * copt19 * copt54 * copt61 * copt8 +
           2 * C0321 * copt19 * copt53 * copt64 * copt8 +
           2 * C0422 * copt19 * copt69 * copt76 * copt8 +
           2 * C0421 * copt19 * copt68 * copt79 * copt8 +
           C0513 * copt21 * copt85 * copt89 + 3 * C03 * copt17 * copt9 +
           3 * C0131 * copt17 * copt23 * copt34 * copt9 +
           3 * C0231 * copt17 * copt38 * copt49 * copt9 +
           3 * C0331 * copt17 * copt53 * copt64 * copt9 +
           3 * C0431 * copt17 * copt68 * copt79 * copt9 +
           C0512 * copt21 * copt84 * copt91 +
           2 * C0522 * copt19 * copt8 * copt84 * copt91 +
           C0511 * copt21 * copt83 * copt94 +
           2 * C0521 * copt19 * copt8 * copt83 * copt94 +
           3 * C0531 * copt17 * copt83 * copt9 * copt94;
  out(1) =
      copt100 *
      (4 * C14 * copt25 + 3 * C13 * copt24 * copt26 +
       2 * C12 * copt23 * copt27 + C11 * copt28 +
       C0131 * copt10 * copt17 * copt28 + C1213 * copt28 * copt40 * copt44 +
       2 * C1222 * copt23 * copt27 * copt39 * copt46 +
       C1212 * copt28 * copt39 * copt46 +
       3 * C1231 * copt24 * copt26 * copt38 * copt49 +
       2 * C1221 * copt23 * copt27 * copt38 * copt49 +
       C1211 * copt28 * copt38 * copt49 + C1313 * copt28 * copt55 * copt59 +
       2 * C1322 * copt23 * copt27 * copt54 * copt61 +
       C1312 * copt28 * copt54 * copt61 +
       3 * C1331 * copt24 * copt26 * copt53 * copt64 +
       2 * C1321 * copt23 * copt27 * copt53 * copt64 +
       C1311 * copt28 * copt53 * copt64 + C1413 * copt28 * copt70 * copt74 +
       2 * C1422 * copt23 * copt27 * copt69 * copt76 +
       C1412 * copt28 * copt69 * copt76 +
       3 * C1431 * copt24 * copt26 * copt68 * copt79 +
       2 * C1421 * copt23 * copt27 * copt68 * copt79 +
       C1411 * copt28 * copt68 * copt79 +
       3 * C0113 * copt21 * copt24 * copt26 * copt8 +
       2 * C0112 * copt21 * copt23 * copt27 * copt8 +
       C0111 * copt21 * copt28 * copt8 + C1513 * copt28 * copt85 * copt89 +
       2 * C0122 * copt19 * copt23 * copt27 * copt9 +
       C0121 * copt19 * copt28 * copt9 +
       2 * C1522 * copt23 * copt27 * copt84 * copt91 +
       C1512 * copt28 * copt84 * copt91 +
       3 * C1531 * copt24 * copt26 * copt83 * copt94 +
       2 * C1521 * copt23 * copt27 * copt83 * copt94 +
       C1511 * copt28 * copt83 * copt94);
  out(2) = 4 * C24 * copt139 * copt40 + 3 * C23 * copt39 * copt44 +
           3 * C1213 * copt23 * copt34 * copt39 * copt44 +
           2 * C22 * copt38 * copt46 +
           2 * C1222 * copt24 * copt31 * copt38 * copt46 +
           2 * C1212 * copt23 * copt34 * copt38 * copt46 + C21 * copt49 +
           C0231 * copt10 * copt17 * copt49 + C1231 * copt25 * copt29 * copt49 +
           C1221 * copt24 * copt31 * copt49 + C1211 * copt23 * copt34 * copt49 +
           C2313 * copt49 * copt55 * copt59 - Ccompr0 * copt2 * copt6 +
           2 * C2322 * copt38 * copt46 * copt54 * copt61 +
           C2312 * copt49 * copt54 * copt61 +
           3 * C2331 * copt39 * copt44 * copt53 * copt64 +
           2 * C2321 * copt38 * copt46 * copt53 * copt64 +
           C2311 * copt49 * copt53 * copt64 + C2413 * copt49 * copt70 * copt74 +
           2 * C2422 * copt38 * copt46 * copt69 * copt76 +
           C2412 * copt49 * copt69 * copt76 +
           3 * C2431 * copt39 * copt44 * copt68 * copt79 +
           2 * C2421 * copt38 * copt46 * copt68 * copt79 +
           C2411 * copt49 * copt68 * copt79 +
           3 * C0213 * copt21 * copt39 * copt44 * copt8 +
           2 * C0212 * copt21 * copt38 * copt46 * copt8 +
           C0211 * copt21 * copt49 * copt8 + C2513 * copt49 * copt85 * copt89 +
           2 * C0222 * copt19 * copt38 * copt46 * copt9 +
           C0221 * copt19 * copt49 * copt9 +
           2 * C2522 * copt38 * copt46 * copt84 * copt91 +
           C2512 * copt49 * copt84 * copt91 +
           3 * C2531 * copt39 * copt44 * copt83 * copt94 +
           2 * C2521 * copt38 * copt46 * copt83 * copt94 +
           C2511 * copt49 * copt83 * copt94;
  out(3) =
      copt176 *
      (4 * C34 * copt55 + 3 * C33 * copt54 * copt56 +
       3 * C1313 * copt23 * copt34 * copt54 * copt56 +
       3 * C2313 * copt38 * copt49 * copt54 * copt56 +
       2 * C32 * copt53 * copt57 +
       2 * C1322 * copt24 * copt31 * copt53 * copt57 +
       2 * C1312 * copt23 * copt34 * copt53 * copt57 +
       2 * C2322 * copt39 * copt46 * copt53 * copt57 +
       2 * C2312 * copt38 * copt49 * copt53 * copt57 + C31 * copt58 +
       C0331 * copt10 * copt17 * copt58 + C1331 * copt25 * copt29 * copt58 +
       C1321 * copt24 * copt31 * copt58 + C1311 * copt23 * copt34 * copt58 +
       C2331 * copt40 * copt44 * copt58 + C2321 * copt39 * copt46 * copt58 +
       C2311 * copt38 * copt49 * copt58 + C3413 * copt58 * copt70 * copt74 +
       2 * C3422 * copt53 * copt57 * copt69 * copt76 +
       C3412 * copt58 * copt69 * copt76 +
       3 * C3431 * copt54 * copt56 * copt68 * copt79 +
       2 * C3421 * copt53 * copt57 * copt68 * copt79 +
       C3411 * copt58 * copt68 * copt79 +
       3 * C0313 * copt21 * copt54 * copt56 * copt8 +
       2 * C0312 * copt21 * copt53 * copt57 * copt8 +
       C0311 * copt21 * copt58 * copt8 + C3513 * copt58 * copt85 * copt89 +
       2 * C0322 * copt19 * copt53 * copt57 * copt9 +
       C0321 * copt19 * copt58 * copt9 +
       2 * C3522 * copt53 * copt57 * copt84 * copt91 +
       C3512 * copt58 * copt84 * copt91 +
       3 * C3531 * copt54 * copt56 * copt83 * copt94 +
       2 * C3521 * copt53 * copt57 * copt83 * copt94 +
       C3511 * copt58 * copt83 * copt94);
  out(4) =
      copt214 *
      (4 * C44 * copt70 + 3 * C43 * copt69 * copt71 +
       3 * C1413 * copt23 * copt34 * copt69 * copt71 +
       3 * C2413 * copt38 * copt49 * copt69 * copt71 +
       3 * C3413 * copt53 * copt64 * copt69 * copt71 +
       2 * C42 * copt68 * copt72 +
       2 * C1422 * copt24 * copt31 * copt68 * copt72 +
       2 * C1412 * copt23 * copt34 * copt68 * copt72 +
       2 * C2422 * copt39 * copt46 * copt68 * copt72 +
       2 * C2412 * copt38 * copt49 * copt68 * copt72 +
       2 * C3422 * copt54 * copt61 * copt68 * copt72 +
       2 * C3412 * copt53 * copt64 * copt68 * copt72 + C41 * copt73 +
       C0431 * copt10 * copt17 * copt73 + C1431 * copt25 * copt29 * copt73 +
       C1421 * copt24 * copt31 * copt73 + C1411 * copt23 * copt34 * copt73 +
       C2431 * copt40 * copt44 * copt73 + C2421 * copt39 * copt46 * copt73 +
       C2411 * copt38 * copt49 * copt73 + C3431 * copt55 * copt59 * copt73 +
       C3421 * copt54 * copt61 * copt73 + C3411 * copt53 * copt64 * copt73 +
       3 * C0413 * copt21 * copt69 * copt71 * copt8 +
       2 * C0412 * copt21 * copt68 * copt72 * copt8 +
       C0411 * copt21 * copt73 * copt8 + C4513 * copt73 * copt85 * copt89 +
       2 * C0422 * copt19 * copt68 * copt72 * copt9 +
       C0421 * copt19 * copt73 * copt9 +
       2 * C4522 * copt68 * copt72 * copt84 * copt91 +
       C4512 * copt73 * copt84 * copt91 +
       3 * C4531 * copt69 * copt71 * copt83 * copt94 +
       2 * C4521 * copt68 * copt72 * copt83 * copt94 +
       C4511 * copt73 * copt83 * copt94);
  out(5) =
      copt252 *
      (4 * C54 * copt85 + 3 * C53 * copt84 * copt86 +
       3 * C1513 * copt23 * copt34 * copt84 * copt86 +
       3 * C2513 * copt38 * copt49 * copt84 * copt86 +
       3 * C3513 * copt53 * copt64 * copt84 * copt86 +
       3 * C4513 * copt68 * copt79 * copt84 * copt86 +
       3 * C0513 * copt21 * copt8 * copt84 * copt86 +
       2 * C52 * copt83 * copt87 +
       2 * C1522 * copt24 * copt31 * copt83 * copt87 +
       2 * C1512 * copt23 * copt34 * copt83 * copt87 +
       2 * C2522 * copt39 * copt46 * copt83 * copt87 +
       2 * C2512 * copt38 * copt49 * copt83 * copt87 +
       2 * C3522 * copt54 * copt61 * copt83 * copt87 +
       2 * C3512 * copt53 * copt64 * copt83 * copt87 +
       2 * C4522 * copt69 * copt76 * copt83 * copt87 +
       2 * C4512 * copt68 * copt79 * copt83 * copt87 +
       2 * C0512 * copt21 * copt8 * copt83 * copt87 + C51 * copt88 +
       C0531 * copt10 * copt17 * copt88 + C1531 * copt25 * copt29 * copt88 +
       C1521 * copt24 * copt31 * copt88 + C1511 * copt23 * copt34 * copt88 +
       C2531 * copt40 * copt44 * copt88 + C2521 * copt39 * copt46 * copt88 +
       C2511 * copt38 * copt49 * copt88 + C3531 * copt55 * copt59 * copt88 +
       C3521 * copt54 * copt61 * copt88 + C3511 * copt53 * copt64 * copt88 +
       C4531 * copt70 * copt74 * copt88 + C4521 * copt69 * copt76 * copt88 +
       C4511 * copt68 * copt79 * copt88 + C0511 * copt21 * copt8 * copt88 +
       2 * C0522 * copt19 * copt83 * copt87 * copt9 +
       C0521 * copt19 * copt88 * copt9);

  return out;
}

gradhess_type FittedMaterial::gradhess_taylor_0(const strain_type &straincpp) {
  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real copt2   = straincpp(0);
  Real copt1   = straincpp(2);
  Real copt3   = copt1 * copt2;
  Real copt4   = Ccompr1 + copt3;
  Real copt5   = Power(copt4, 2);
  Real copt6   = 1 / copt5;
  Real copt8   = -1 + copt2;
  Real copt9   = Power(copt8, 2);
  Real copt10  = copt8 * copt9;
  Real copt11  = strainscale(0);
  Real copt12  = Power(copt11, 2);
  Real copt13  = Power(copt12, 2);
  Real copt14  = 1 / copt13;
  Real copt16  = copt11 * copt12;
  Real copt17  = 1 / copt16;
  Real copt19  = 1 / copt12;
  Real copt21  = 1 / copt11;
  Real copt23  = straincpp(1);
  Real copt24  = Power(copt23, 2);
  Real copt25  = copt23 * copt24;
  Real copt26  = strainscale(1);
  Real copt27  = Power(copt26, 2);
  Real copt28  = copt26 * copt27;
  Real copt29  = 1 / copt28;
  Real copt31  = 1 / copt27;
  Real copt34  = 1 / copt26;
  Real copt38  = -1 + copt1;
  Real copt39  = Power(copt38, 2);
  Real copt40  = copt38 * copt39;
  Real copt41  = strainscale(2);
  Real copt42  = Power(copt41, 2);
  Real copt43  = copt41 * copt42;
  Real copt44  = 1 / copt43;
  Real copt46  = 1 / copt42;
  Real copt49  = 1 / copt41;
  Real copt53  = straincpp(3);
  Real copt54  = Power(copt53, 2);
  Real copt55  = copt53 * copt54;
  Real copt56  = strainscale(3);
  Real copt57  = Power(copt56, 2);
  Real copt58  = copt56 * copt57;
  Real copt59  = 1 / copt58;
  Real copt61  = 1 / copt57;
  Real copt64  = 1 / copt56;
  Real copt68  = straincpp(4);
  Real copt69  = Power(copt68, 2);
  Real copt70  = copt68 * copt69;
  Real copt71  = strainscale(4);
  Real copt72  = Power(copt71, 2);
  Real copt73  = copt71 * copt72;
  Real copt74  = 1 / copt73;
  Real copt76  = 1 / copt72;
  Real copt79  = 1 / copt71;
  Real copt83  = straincpp(5);
  Real copt84  = Power(copt83, 2);
  Real copt85  = copt83 * copt84;
  Real copt86  = strainscale(5);
  Real copt87  = Power(copt86, 2);
  Real copt88  = copt86 * copt87;
  Real copt89  = 1 / copt88;
  Real copt91  = 1 / copt87;
  Real copt94  = 1 / copt86;
  Real copt99  = Power(copt27, 2);
  Real copt100 = 1 / copt99;
  Real copt138 = Power(copt42, 2);
  Real copt139 = 1 / copt138;
  Real copt175 = Power(copt57, 2);
  Real copt176 = 1 / copt175;
  Real copt213 = Power(copt72, 2);
  Real copt214 = 1 / copt213;
  Real copt251 = Power(copt87, 2);
  Real copt252 = 1 / copt251;
  Real copt289 = Power(copt1, 2);
  Real copt290 = copt4 * copt5;
  Real copt291 = 1 / copt290;
  Real copt313 = 3 * C0113 * copt12 * copt24;
  Real copt314 = 4 * C0122 * copt11 * copt23 * copt8;
  Real copt315 = 2 * C0112 * copt12 * copt23;
  Real copt316 = 3 * C0131 * copt9;
  Real copt317 = 2 * C0121 * copt8;
  Real copt318 = C0111 * copt11;
  Real copt319 = copt317 + copt318;
  Real copt320 = copt11 * copt319;
  Real copt321 = copt316 + copt320;
  Real copt322 = copt26 * copt321;
  Real copt323 = copt314 + copt315 + copt322;
  Real copt324 = copt26 * copt323;
  Real copt325 = copt313 + copt324;
  Real copt326 = copt17 * copt29 * copt325;
  Real copt327 = -Ccompr1;
  Real copt328 = copt3 + copt327;
  Real copt329 = Ccompr0 * copt291 * copt328;
  Real copt330 = 3 * C0213 * copt12 * copt39;
  Real copt331 = 4 * C0222 * copt11 * copt38 * copt8;
  Real copt332 = 2 * C0212 * copt12 * copt38;
  Real copt333 = 3 * C0231 * copt9;
  Real copt334 = 2 * C0221 * copt8;
  Real copt335 = C0211 * copt11;
  Real copt336 = copt334 + copt335;
  Real copt337 = copt11 * copt336;
  Real copt338 = copt333 + copt337;
  Real copt339 = copt338 * copt41;
  Real copt340 = copt331 + copt332 + copt339;
  Real copt341 = copt340 * copt41;
  Real copt342 = copt330 + copt341;
  Real copt343 = copt17 * copt342 * copt44;
  Real copt344 = copt329 + copt343;
  Real copt407 = 3 * C1213 * copt27 * copt39;
  Real copt408 = 4 * C1222 * copt23 * copt26 * copt38;
  Real copt409 = 2 * C1212 * copt27 * copt38;
  Real copt410 = 3 * C1231 * copt24;
  Real copt411 = 2 * C1221 * copt23 * copt26;
  Real copt412 = C1211 * copt27;
  Real copt413 = copt410 + copt411 + copt412;
  Real copt414 = copt41 * copt413;
  Real copt415 = copt408 + copt409 + copt414;
  Real copt416 = copt41 * copt415;
  Real copt417 = copt407 + copt416;
  Real copt418 = copt29 * copt417 * copt44;
  Real copt455 = Power(copt2, 2);
  Real copt345 = 3 * C0313 * copt12 * copt54;
  Real copt346 = 4 * C0322 * copt11 * copt53 * copt8;
  Real copt347 = 2 * C0312 * copt12 * copt53;
  Real copt348 = 3 * C0331 * copt9;
  Real copt349 = 2 * C0321 * copt8;
  Real copt350 = C0311 * copt11;
  Real copt351 = copt349 + copt350;
  Real copt352 = copt11 * copt351;
  Real copt353 = copt348 + copt352;
  Real copt354 = copt353 * copt56;
  Real copt355 = copt346 + copt347 + copt354;
  Real copt356 = copt355 * copt56;
  Real copt357 = copt345 + copt356;
  Real copt358 = copt17 * copt357 * copt59;
  Real copt419 = 3 * C1313 * copt27 * copt54;
  Real copt420 = 4 * C1322 * copt23 * copt26 * copt53;
  Real copt421 = 2 * C1312 * copt27 * copt53;
  Real copt422 = 3 * C1331 * copt24;
  Real copt423 = 2 * C1321 * copt23 * copt26;
  Real copt424 = C1311 * copt27;
  Real copt425 = copt422 + copt423 + copt424;
  Real copt426 = copt425 * copt56;
  Real copt427 = copt420 + copt421 + copt426;
  Real copt428 = copt427 * copt56;
  Real copt429 = copt419 + copt428;
  Real copt430 = copt29 * copt429 * copt59;
  Real copt477 = 3 * C2313 * copt42 * copt54;
  Real copt478 = 4 * C2322 * copt38 * copt41 * copt53;
  Real copt479 = 2 * C2312 * copt42 * copt53;
  Real copt480 = 3 * C2331 * copt39;
  Real copt481 = 2 * C2321 * copt38;
  Real copt482 = C2311 * copt41;
  Real copt483 = copt481 + copt482;
  Real copt484 = copt41 * copt483;
  Real copt485 = copt480 + copt484;
  Real copt486 = copt485 * copt56;
  Real copt487 = copt478 + copt479 + copt486;
  Real copt488 = copt487 * copt56;
  Real copt489 = copt477 + copt488;
  Real copt490 = copt44 * copt489 * copt59;
  Real copt359 = 3 * C0413 * copt12 * copt69;
  Real copt360 = 4 * C0422 * copt11 * copt68 * copt8;
  Real copt361 = 2 * C0412 * copt12 * copt68;
  Real copt362 = 3 * C0431 * copt9;
  Real copt363 = 2 * C0421 * copt8;
  Real copt364 = C0411 * copt11;
  Real copt365 = copt363 + copt364;
  Real copt366 = copt11 * copt365;
  Real copt367 = copt362 + copt366;
  Real copt368 = copt367 * copt71;
  Real copt369 = copt360 + copt361 + copt368;
  Real copt370 = copt369 * copt71;
  Real copt371 = copt359 + copt370;
  Real copt372 = copt17 * copt371 * copt74;
  Real copt431 = 3 * C1413 * copt27 * copt69;
  Real copt432 = 4 * C1422 * copt23 * copt26 * copt68;
  Real copt433 = 2 * C1412 * copt27 * copt68;
  Real copt434 = 3 * C1431 * copt24;
  Real copt435 = 2 * C1421 * copt23 * copt26;
  Real copt436 = C1411 * copt27;
  Real copt437 = copt434 + copt435 + copt436;
  Real copt438 = copt437 * copt71;
  Real copt439 = copt432 + copt433 + copt438;
  Real copt440 = copt439 * copt71;
  Real copt441 = copt431 + copt440;
  Real copt442 = copt29 * copt441 * copt74;
  Real copt491 = 3 * C2413 * copt42 * copt69;
  Real copt492 = 4 * C2422 * copt38 * copt41 * copt68;
  Real copt493 = 2 * C2412 * copt42 * copt68;
  Real copt494 = 3 * C2431 * copt39;
  Real copt495 = 2 * C2421 * copt38;
  Real copt496 = C2411 * copt41;
  Real copt497 = copt495 + copt496;
  Real copt498 = copt41 * copt497;
  Real copt499 = copt494 + copt498;
  Real copt500 = copt499 * copt71;
  Real copt501 = copt492 + copt493 + copt500;
  Real copt502 = copt501 * copt71;
  Real copt503 = copt491 + copt502;
  Real copt504 = copt44 * copt503 * copt74;
  Real copt539 = 3 * C3413 * copt57 * copt69;
  Real copt540 = 4 * C3422 * copt53 * copt56 * copt68;
  Real copt541 = 2 * C3412 * copt57 * copt68;
  Real copt542 = 3 * C3431 * copt54;
  Real copt543 = 2 * C3421 * copt53 * copt56;
  Real copt544 = C3411 * copt57;
  Real copt545 = copt542 + copt543 + copt544;
  Real copt546 = copt545 * copt71;
  Real copt547 = copt540 + copt541 + copt546;
  Real copt548 = copt547 * copt71;
  Real copt549 = copt539 + copt548;
  Real copt550 = copt549 * copt59 * copt74;
  Real copt373 = 3 * C0513 * copt12 * copt84;
  Real copt374 = 4 * C0522 * copt11 * copt8 * copt83;
  Real copt375 = 2 * C0512 * copt12 * copt83;
  Real copt376 = 3 * C0531 * copt9;
  Real copt377 = 2 * C0521 * copt8;
  Real copt378 = C0511 * copt11;
  Real copt379 = copt377 + copt378;
  Real copt380 = copt11 * copt379;
  Real copt381 = copt376 + copt380;
  Real copt382 = copt381 * copt86;
  Real copt383 = copt374 + copt375 + copt382;
  Real copt384 = copt383 * copt86;
  Real copt385 = copt373 + copt384;
  Real copt386 = copt17 * copt385 * copt89;
  Real copt443 = 3 * C1513 * copt27 * copt84;
  Real copt444 = 4 * C1522 * copt23 * copt26 * copt83;
  Real copt445 = 2 * C1512 * copt27 * copt83;
  Real copt446 = 3 * C1531 * copt24;
  Real copt447 = 2 * C1521 * copt23 * copt26;
  Real copt448 = C1511 * copt27;
  Real copt449 = copt446 + copt447 + copt448;
  Real copt450 = copt449 * copt86;
  Real copt451 = copt444 + copt445 + copt450;
  Real copt452 = copt451 * copt86;
  Real copt453 = copt443 + copt452;
  Real copt454 = copt29 * copt453 * copt89;
  Real copt505 = 3 * C2513 * copt42 * copt84;
  Real copt506 = 4 * C2522 * copt38 * copt41 * copt83;
  Real copt507 = 2 * C2512 * copt42 * copt83;
  Real copt508 = 3 * C2531 * copt39;
  Real copt509 = 2 * C2521 * copt38;
  Real copt510 = C2511 * copt41;
  Real copt511 = copt509 + copt510;
  Real copt512 = copt41 * copt511;
  Real copt513 = copt508 + copt512;
  Real copt514 = copt513 * copt86;
  Real copt515 = copt506 + copt507 + copt514;
  Real copt516 = copt515 * copt86;
  Real copt517 = copt505 + copt516;
  Real copt518 = copt44 * copt517 * copt89;
  Real copt551 = 3 * C3513 * copt57 * copt84;
  Real copt552 = 4 * C3522 * copt53 * copt56 * copt83;
  Real copt553 = 2 * C3512 * copt57 * copt83;
  Real copt554 = 3 * C3531 * copt54;
  Real copt555 = 2 * C3521 * copt53 * copt56;
  Real copt556 = C3511 * copt57;
  Real copt557 = copt554 + copt555 + copt556;
  Real copt558 = copt557 * copt86;
  Real copt559 = copt552 + copt553 + copt558;
  Real copt560 = copt559 * copt86;
  Real copt561 = copt551 + copt560;
  Real copt562 = copt561 * copt59 * copt89;
  Real copt583 = 3 * C4513 * copt72 * copt84;
  Real copt584 = 4 * C4522 * copt68 * copt71 * copt83;
  Real copt585 = 2 * C4512 * copt72 * copt83;
  Real copt586 = 3 * C4531 * copt69;
  Real copt587 = 2 * C4521 * copt68 * copt71;
  Real copt588 = C4511 * copt72;
  Real copt589 = copt586 + copt587 + copt588;
  Real copt590 = copt589 * copt86;
  Real copt591 = copt584 + copt585 + copt590;
  Real copt592 = copt591 * copt86;
  Real copt593 = copt583 + copt592;
  Real copt594 = copt593 * copt74 * copt89;
  out1(0) =
      4 * C04 * copt10 * copt14 + C01 * copt21 +
      C0113 * copt21 * copt25 * copt29 + C0112 * copt21 * copt24 * copt31 +
      C0111 * copt21 * copt23 * copt34 + C0213 * copt21 * copt40 * copt44 +
      C0212 * copt21 * copt39 * copt46 + C0211 * copt21 * copt38 * copt49 +
      C0313 * copt21 * copt55 * copt59 - Ccompr0 * copt1 * copt6 +
      C0312 * copt21 * copt54 * copt61 + C0311 * copt21 * copt53 * copt64 +
      C0413 * copt21 * copt70 * copt74 + C0412 * copt21 * copt69 * copt76 +
      C0411 * copt21 * copt68 * copt79 + 2 * C02 * copt19 * copt8 +
      2 * C0122 * copt19 * copt24 * copt31 * copt8 +
      2 * C0121 * copt19 * copt23 * copt34 * copt8 +
      2 * C0222 * copt19 * copt39 * copt46 * copt8 +
      2 * C0221 * copt19 * copt38 * copt49 * copt8 +
      2 * C0322 * copt19 * copt54 * copt61 * copt8 +
      2 * C0321 * copt19 * copt53 * copt64 * copt8 +
      2 * C0422 * copt19 * copt69 * copt76 * copt8 +
      2 * C0421 * copt19 * copt68 * copt79 * copt8 +
      C0513 * copt21 * copt85 * copt89 + 3 * C03 * copt17 * copt9 +
      3 * C0131 * copt17 * copt23 * copt34 * copt9 +
      3 * C0231 * copt17 * copt38 * copt49 * copt9 +
      3 * C0331 * copt17 * copt53 * copt64 * copt9 +
      3 * C0431 * copt17 * copt68 * copt79 * copt9 +
      C0512 * copt21 * copt84 * copt91 +
      2 * C0522 * copt19 * copt8 * copt84 * copt91 +
      C0511 * copt21 * copt83 * copt94 +
      2 * C0521 * copt19 * copt8 * copt83 * copt94 +
      3 * C0531 * copt17 * copt83 * copt9 * copt94;
  out1(1) =
      copt100 *
      (4 * C14 * copt25 + 3 * C13 * copt24 * copt26 +
       2 * C12 * copt23 * copt27 + C11 * copt28 +
       C0131 * copt10 * copt17 * copt28 + C1213 * copt28 * copt40 * copt44 +
       2 * C1222 * copt23 * copt27 * copt39 * copt46 +
       C1212 * copt28 * copt39 * copt46 +
       3 * C1231 * copt24 * copt26 * copt38 * copt49 +
       2 * C1221 * copt23 * copt27 * copt38 * copt49 +
       C1211 * copt28 * copt38 * copt49 + C1313 * copt28 * copt55 * copt59 +
       2 * C1322 * copt23 * copt27 * copt54 * copt61 +
       C1312 * copt28 * copt54 * copt61 +
       3 * C1331 * copt24 * copt26 * copt53 * copt64 +
       2 * C1321 * copt23 * copt27 * copt53 * copt64 +
       C1311 * copt28 * copt53 * copt64 + C1413 * copt28 * copt70 * copt74 +
       2 * C1422 * copt23 * copt27 * copt69 * copt76 +
       C1412 * copt28 * copt69 * copt76 +
       3 * C1431 * copt24 * copt26 * copt68 * copt79 +
       2 * C1421 * copt23 * copt27 * copt68 * copt79 +
       C1411 * copt28 * copt68 * copt79 +
       3 * C0113 * copt21 * copt24 * copt26 * copt8 +
       2 * C0112 * copt21 * copt23 * copt27 * copt8 +
       C0111 * copt21 * copt28 * copt8 + C1513 * copt28 * copt85 * copt89 +
       2 * C0122 * copt19 * copt23 * copt27 * copt9 +
       C0121 * copt19 * copt28 * copt9 +
       2 * C1522 * copt23 * copt27 * copt84 * copt91 +
       C1512 * copt28 * copt84 * copt91 +
       3 * C1531 * copt24 * copt26 * copt83 * copt94 +
       2 * C1521 * copt23 * copt27 * copt83 * copt94 +
       C1511 * copt28 * copt83 * copt94);
  out1(2) =
      4 * C24 * copt139 * copt40 + 3 * C23 * copt39 * copt44 +
      3 * C1213 * copt23 * copt34 * copt39 * copt44 +
      2 * C22 * copt38 * copt46 +
      2 * C1222 * copt24 * copt31 * copt38 * copt46 +
      2 * C1212 * copt23 * copt34 * copt38 * copt46 + C21 * copt49 +
      C0231 * copt10 * copt17 * copt49 + C1231 * copt25 * copt29 * copt49 +
      C1221 * copt24 * copt31 * copt49 + C1211 * copt23 * copt34 * copt49 +
      C2313 * copt49 * copt55 * copt59 - Ccompr0 * copt2 * copt6 +
      2 * C2322 * copt38 * copt46 * copt54 * copt61 +
      C2312 * copt49 * copt54 * copt61 +
      3 * C2331 * copt39 * copt44 * copt53 * copt64 +
      2 * C2321 * copt38 * copt46 * copt53 * copt64 +
      C2311 * copt49 * copt53 * copt64 + C2413 * copt49 * copt70 * copt74 +
      2 * C2422 * copt38 * copt46 * copt69 * copt76 +
      C2412 * copt49 * copt69 * copt76 +
      3 * C2431 * copt39 * copt44 * copt68 * copt79 +
      2 * C2421 * copt38 * copt46 * copt68 * copt79 +
      C2411 * copt49 * copt68 * copt79 +
      3 * C0213 * copt21 * copt39 * copt44 * copt8 +
      2 * C0212 * copt21 * copt38 * copt46 * copt8 +
      C0211 * copt21 * copt49 * copt8 + C2513 * copt49 * copt85 * copt89 +
      2 * C0222 * copt19 * copt38 * copt46 * copt9 +
      C0221 * copt19 * copt49 * copt9 +
      2 * C2522 * copt38 * copt46 * copt84 * copt91 +
      C2512 * copt49 * copt84 * copt91 +
      3 * C2531 * copt39 * copt44 * copt83 * copt94 +
      2 * C2521 * copt38 * copt46 * copt83 * copt94 +
      C2511 * copt49 * copt83 * copt94;
  out1(3) =
      copt176 *
      (4 * C34 * copt55 + 3 * C33 * copt54 * copt56 +
       3 * C1313 * copt23 * copt34 * copt54 * copt56 +
       3 * C2313 * copt38 * copt49 * copt54 * copt56 +
       2 * C32 * copt53 * copt57 +
       2 * C1322 * copt24 * copt31 * copt53 * copt57 +
       2 * C1312 * copt23 * copt34 * copt53 * copt57 +
       2 * C2322 * copt39 * copt46 * copt53 * copt57 +
       2 * C2312 * copt38 * copt49 * copt53 * copt57 + C31 * copt58 +
       C0331 * copt10 * copt17 * copt58 + C1331 * copt25 * copt29 * copt58 +
       C1321 * copt24 * copt31 * copt58 + C1311 * copt23 * copt34 * copt58 +
       C2331 * copt40 * copt44 * copt58 + C2321 * copt39 * copt46 * copt58 +
       C2311 * copt38 * copt49 * copt58 + C3413 * copt58 * copt70 * copt74 +
       2 * C3422 * copt53 * copt57 * copt69 * copt76 +
       C3412 * copt58 * copt69 * copt76 +
       3 * C3431 * copt54 * copt56 * copt68 * copt79 +
       2 * C3421 * copt53 * copt57 * copt68 * copt79 +
       C3411 * copt58 * copt68 * copt79 +
       3 * C0313 * copt21 * copt54 * copt56 * copt8 +
       2 * C0312 * copt21 * copt53 * copt57 * copt8 +
       C0311 * copt21 * copt58 * copt8 + C3513 * copt58 * copt85 * copt89 +
       2 * C0322 * copt19 * copt53 * copt57 * copt9 +
       C0321 * copt19 * copt58 * copt9 +
       2 * C3522 * copt53 * copt57 * copt84 * copt91 +
       C3512 * copt58 * copt84 * copt91 +
       3 * C3531 * copt54 * copt56 * copt83 * copt94 +
       2 * C3521 * copt53 * copt57 * copt83 * copt94 +
       C3511 * copt58 * copt83 * copt94);
  out1(4) =
      copt214 *
      (4 * C44 * copt70 + 3 * C43 * copt69 * copt71 +
       3 * C1413 * copt23 * copt34 * copt69 * copt71 +
       3 * C2413 * copt38 * copt49 * copt69 * copt71 +
       3 * C3413 * copt53 * copt64 * copt69 * copt71 +
       2 * C42 * copt68 * copt72 +
       2 * C1422 * copt24 * copt31 * copt68 * copt72 +
       2 * C1412 * copt23 * copt34 * copt68 * copt72 +
       2 * C2422 * copt39 * copt46 * copt68 * copt72 +
       2 * C2412 * copt38 * copt49 * copt68 * copt72 +
       2 * C3422 * copt54 * copt61 * copt68 * copt72 +
       2 * C3412 * copt53 * copt64 * copt68 * copt72 + C41 * copt73 +
       C0431 * copt10 * copt17 * copt73 + C1431 * copt25 * copt29 * copt73 +
       C1421 * copt24 * copt31 * copt73 + C1411 * copt23 * copt34 * copt73 +
       C2431 * copt40 * copt44 * copt73 + C2421 * copt39 * copt46 * copt73 +
       C2411 * copt38 * copt49 * copt73 + C3431 * copt55 * copt59 * copt73 +
       C3421 * copt54 * copt61 * copt73 + C3411 * copt53 * copt64 * copt73 +
       3 * C0413 * copt21 * copt69 * copt71 * copt8 +
       2 * C0412 * copt21 * copt68 * copt72 * copt8 +
       C0411 * copt21 * copt73 * copt8 + C4513 * copt73 * copt85 * copt89 +
       2 * C0422 * copt19 * copt68 * copt72 * copt9 +
       C0421 * copt19 * copt73 * copt9 +
       2 * C4522 * copt68 * copt72 * copt84 * copt91 +
       C4512 * copt73 * copt84 * copt91 +
       3 * C4531 * copt69 * copt71 * copt83 * copt94 +
       2 * C4521 * copt68 * copt72 * copt83 * copt94 +
       C4511 * copt73 * copt83 * copt94);
  out1(5) =
      copt252 *
      (4 * C54 * copt85 + 3 * C53 * copt84 * copt86 +
       3 * C1513 * copt23 * copt34 * copt84 * copt86 +
       3 * C2513 * copt38 * copt49 * copt84 * copt86 +
       3 * C3513 * copt53 * copt64 * copt84 * copt86 +
       3 * C4513 * copt68 * copt79 * copt84 * copt86 +
       3 * C0513 * copt21 * copt8 * copt84 * copt86 +
       2 * C52 * copt83 * copt87 +
       2 * C1522 * copt24 * copt31 * copt83 * copt87 +
       2 * C1512 * copt23 * copt34 * copt83 * copt87 +
       2 * C2522 * copt39 * copt46 * copt83 * copt87 +
       2 * C2512 * copt38 * copt49 * copt83 * copt87 +
       2 * C3522 * copt54 * copt61 * copt83 * copt87 +
       2 * C3512 * copt53 * copt64 * copt83 * copt87 +
       2 * C4522 * copt69 * copt76 * copt83 * copt87 +
       2 * C4512 * copt68 * copt79 * copt83 * copt87 +
       2 * C0512 * copt21 * copt8 * copt83 * copt87 + C51 * copt88 +
       C0531 * copt10 * copt17 * copt88 + C1531 * copt25 * copt29 * copt88 +
       C1521 * copt24 * copt31 * copt88 + C1511 * copt23 * copt34 * copt88 +
       C2531 * copt40 * copt44 * copt88 + C2521 * copt39 * copt46 * copt88 +
       C2511 * copt38 * copt49 * copt88 + C3531 * copt55 * copt59 * copt88 +
       C3521 * copt54 * copt61 * copt88 + C3511 * copt53 * copt64 * copt88 +
       C4531 * copt70 * copt74 * copt88 + C4521 * copt69 * copt76 * copt88 +
       C4511 * copt68 * copt79 * copt88 + C0511 * copt21 * copt8 * copt88 +
       2 * C0522 * copt19 * copt83 * copt87 * copt9 +
       C0521 * copt19 * copt88 * copt9);
  out2(0, 0) =
      2 *
      (C02 * copt19 + Ccompr0 * copt289 * copt291 +
       C0122 * copt19 * copt24 * copt31 + C0121 * copt19 * copt23 * copt34 +
       C0222 * copt19 * copt39 * copt46 + C0221 * copt19 * copt38 * copt49 +
       C0322 * copt19 * copt54 * copt61 + C0321 * copt19 * copt53 * copt64 +
       C0422 * copt19 * copt69 * copt76 + C0421 * copt19 * copt68 * copt79 +
       3 * C03 * copt17 * copt8 + 3 * C0131 * copt17 * copt23 * copt34 * copt8 +
       3 * C0231 * copt17 * copt38 * copt49 * copt8 +
       3 * C0331 * copt17 * copt53 * copt64 * copt8 +
       3 * C0431 * copt17 * copt68 * copt79 * copt8 + 6 * C04 * copt14 * copt9 +
       C0522 * copt19 * copt84 * copt91 + C0521 * copt19 * copt83 * copt94 +
       3 * C0531 * copt17 * copt8 * copt83 * copt94);
  out2(0, 1) = copt326;
  out2(0, 2) = copt344;
  out2(0, 3) = copt358;
  out2(0, 4) = copt372;
  out2(0, 5) = copt386;
  out2(1, 0) = copt326;
  out2(1, 1) =
      2 * copt100 *
      (6 * C14 * copt24 + 3 * C13 * copt23 * copt26 + C12 * copt27 +
       C1222 * copt27 * copt39 * copt46 +
       3 * C1231 * copt23 * copt26 * copt38 * copt49 +
       C1221 * copt27 * copt38 * copt49 + C1322 * copt27 * copt54 * copt61 +
       3 * C1331 * copt23 * copt26 * copt53 * copt64 +
       C1321 * copt27 * copt53 * copt64 + C1422 * copt27 * copt69 * copt76 +
       3 * C1431 * copt23 * copt26 * copt68 * copt79 +
       C1421 * copt27 * copt68 * copt79 +
       3 * C0113 * copt21 * copt23 * copt26 * copt8 +
       C0112 * copt21 * copt27 * copt8 + C0122 * copt19 * copt27 * copt9 +
       C1522 * copt27 * copt84 * copt91 +
       3 * C1531 * copt23 * copt26 * copt83 * copt94 +
       C1521 * copt27 * copt83 * copt94);
  out2(1, 2) = copt418;
  out2(1, 3) = copt430;
  out2(1, 4) = copt442;
  out2(1, 5) = copt454;
  out2(2, 0) = copt344;
  out2(2, 1) = copt418;
  out2(2, 2) =
      2 * (6 * C24 * copt139 * copt39 + 3 * C23 * copt38 * copt44 +
           3 * C1213 * copt23 * copt34 * copt38 * copt44 +
           Ccompr0 * copt291 * copt455 + C22 * copt46 +
           C1222 * copt24 * copt31 * copt46 + C1212 * copt23 * copt34 * copt46 +
           C2322 * copt46 * copt54 * copt61 +
           3 * C2331 * copt38 * copt44 * copt53 * copt64 +
           C2321 * copt46 * copt53 * copt64 + C2422 * copt46 * copt69 * copt76 +
           3 * C2431 * copt38 * copt44 * copt68 * copt79 +
           C2421 * copt46 * copt68 * copt79 +
           3 * C0213 * copt21 * copt38 * copt44 * copt8 +
           C0212 * copt21 * copt46 * copt8 + C0222 * copt19 * copt46 * copt9 +
           C2522 * copt46 * copt84 * copt91 +
           3 * C2531 * copt38 * copt44 * copt83 * copt94 +
           C2521 * copt46 * copt83 * copt94);
  out2(2, 3) = copt490;
  out2(2, 4) = copt504;
  out2(2, 5) = copt518;
  out2(3, 0) = copt358;
  out2(3, 1) = copt430;
  out2(3, 2) = copt490;
  out2(3, 3) =
      2 * copt176 *
      (6 * C34 * copt54 + 3 * C33 * copt53 * copt56 +
       3 * C1313 * copt23 * copt34 * copt53 * copt56 +
       3 * C2313 * copt38 * copt49 * copt53 * copt56 + C32 * copt57 +
       C1322 * copt24 * copt31 * copt57 + C1312 * copt23 * copt34 * copt57 +
       C2322 * copt39 * copt46 * copt57 + C2312 * copt38 * copt49 * copt57 +
       C3422 * copt57 * copt69 * copt76 +
       3 * C3431 * copt53 * copt56 * copt68 * copt79 +
       C3421 * copt57 * copt68 * copt79 +
       3 * C0313 * copt21 * copt53 * copt56 * copt8 +
       C0312 * copt21 * copt57 * copt8 + C0322 * copt19 * copt57 * copt9 +
       C3522 * copt57 * copt84 * copt91 +
       3 * C3531 * copt53 * copt56 * copt83 * copt94 +
       C3521 * copt57 * copt83 * copt94);
  out2(3, 4) = copt550;
  out2(3, 5) = copt562;
  out2(4, 0) = copt372;
  out2(4, 1) = copt442;
  out2(4, 2) = copt504;
  out2(4, 3) = copt550;
  out2(4, 4) =
      2 * copt214 *
      (6 * C44 * copt69 + 3 * C43 * copt68 * copt71 +
       3 * C1413 * copt23 * copt34 * copt68 * copt71 +
       3 * C2413 * copt38 * copt49 * copt68 * copt71 +
       3 * C3413 * copt53 * copt64 * copt68 * copt71 + C42 * copt72 +
       C1422 * copt24 * copt31 * copt72 + C1412 * copt23 * copt34 * copt72 +
       C2422 * copt39 * copt46 * copt72 + C2412 * copt38 * copt49 * copt72 +
       C3422 * copt54 * copt61 * copt72 + C3412 * copt53 * copt64 * copt72 +
       3 * C0413 * copt21 * copt68 * copt71 * copt8 +
       C0412 * copt21 * copt72 * copt8 + C0422 * copt19 * copt72 * copt9 +
       C4522 * copt72 * copt84 * copt91 +
       3 * C4531 * copt68 * copt71 * copt83 * copt94 +
       C4521 * copt72 * copt83 * copt94);
  out2(4, 5) = copt594;
  out2(5, 0) = copt386;
  out2(5, 1) = copt454;
  out2(5, 2) = copt518;
  out2(5, 3) = copt562;
  out2(5, 4) = copt594;
  out2(5, 5) =
      2 * copt252 *
      (6 * C54 * copt84 + 3 * C53 * copt83 * copt86 +
       3 * C1513 * copt23 * copt34 * copt83 * copt86 +
       3 * C2513 * copt38 * copt49 * copt83 * copt86 +
       3 * C3513 * copt53 * copt64 * copt83 * copt86 +
       3 * C4513 * copt68 * copt79 * copt83 * copt86 +
       3 * C0513 * copt21 * copt8 * copt83 * copt86 + C52 * copt87 +
       C1522 * copt24 * copt31 * copt87 + C1512 * copt23 * copt34 * copt87 +
       C2522 * copt39 * copt46 * copt87 + C2512 * copt38 * copt49 * copt87 +
       C3522 * copt54 * copt61 * copt87 + C3512 * copt53 * copt64 * copt87 +
       C4522 * copt69 * copt76 * copt87 + C4512 * copt68 * copt79 * copt87 +
       C0512 * copt21 * copt8 * copt87 + C0522 * copt19 * copt87 * copt9);

  return std::make_pair(hess, grad);
}

#define BARRIERPOWER 6  // should be even integer!
value_type FittedMaterial::psi_barrier(const strain_type &straincpp,
                                       const strain_type &strainclamped) {
  Real copt22  = strainscale(0);
  Real copt30  = 1 / copt22;
  Real copt47  = strainscale(1);
  Real copt48  = 1 / copt47;
  Real copt65  = strainscale(2);
  Real copt66  = 1 / copt65;
  Real copt82  = strainscale(3);
  Real copt90  = 1 / copt82;
  Real copt101 = strainscale(4);
  Real copt102 = 1 / copt101;
  Real copt109 = strainscale(5);
  Real copt110 = 1 / copt109;
  return bscale * (Power(bspeed * copt30 * (-strainclamped(0) + straincpp(0)),
                         BARRIERPOWER) +
                   Power(bspeed * copt48 * (-strainclamped(1) + straincpp(1)),
                         BARRIERPOWER) +
                   Power(bspeed * copt66 * (-strainclamped(2) + straincpp(2)),
                         BARRIERPOWER) +
                   Power(bspeed * copt90 * (-strainclamped(3) + straincpp(3)),
                         BARRIERPOWER) +
                   Power(bspeed * copt102 * (-strainclamped(4) + straincpp(4)),
                         BARRIERPOWER) +
                   Power(bspeed * copt110 * (-strainclamped(5) + straincpp(5)),
                         BARRIERPOWER));
}

grad_type FittedMaterial::grad_barrier(const strain_type &straincpp,
                                       const strain_type &strainclamped) {
  Vec6 out(0);
  Real copt22  = strainscale(0);
  Real copt30  = 1 / copt22;
  Real copt50  = strainscale(1);
  Real copt51  = 1 / copt50;
  Real copt33  = -1 + BARRIERPOWER;
  Real copt75  = strainscale(2);
  Real copt77  = 1 / copt75;
  Real copt95  = strainscale(3);
  Real copt96  = 1 / copt95;
  Real copt106 = strainscale(4);
  Real copt107 = 1 / copt106;
  Real copt115 = strainscale(5);
  Real copt116 = 1 / copt115;
  out(0)       = BARRIERPOWER * bscale * bspeed * copt30 *
           Power(bspeed * copt30 * (-strainclamped(0) + straincpp(0)), copt33);
  out(1) = BARRIERPOWER * bscale * bspeed * copt51 *
           Power(bspeed * copt51 * (-strainclamped(1) + straincpp(1)), copt33);
  out(2) = BARRIERPOWER * bscale * bspeed * copt77 *
           Power(bspeed * copt77 * (-strainclamped(2) + straincpp(2)), copt33);
  out(3) = BARRIERPOWER * bscale * bspeed * copt96 *
           Power(bspeed * copt96 * (-strainclamped(3) + straincpp(3)), copt33);
  out(4) = BARRIERPOWER * bscale * bspeed * copt107 *
           Power(bspeed * copt107 * (-strainclamped(4) + straincpp(4)), copt33);
  out(5) = BARRIERPOWER * bscale * bspeed * copt116 *
           Power(bspeed * copt116 * (-strainclamped(5) + straincpp(5)), copt33);

  return out;
}
gradhess_type FittedMaterial::gradhess_barrier(
    const strain_type &straincpp, const strain_type &strainclamped) {
  // define output
  Mat6x6 hess(0);
  Vec6 grad(0);
  auto out1 = [&](int i) -> Real & { return grad[i]; };
  auto out2 = [&](int i, int j) -> Real & { return hess(i, j); };

  Real copt22  = strainscale(0);
  Real copt30  = 1 / copt22;
  Real copt50  = strainscale(1);
  Real copt51  = 1 / copt50;
  Real copt33  = -1 + BARRIERPOWER;
  Real copt75  = strainscale(2);
  Real copt77  = 1 / copt75;
  Real copt95  = strainscale(3);
  Real copt96  = 1 / copt95;
  Real copt106 = strainscale(4);
  Real copt107 = 1 / copt106;
  Real copt115 = strainscale(5);
  Real copt116 = 1 / copt115;
  Real copt120 = Power(bspeed, 2);
  Real copt7   = strainclamped(0);
  Real copt15  = -copt7;
  Real copt18  = straincpp(0);
  Real copt20  = copt15 + copt18;
  Real copt32  = bspeed * copt20 * copt30;
  Real copt123 = Power(copt22, 2);
  Real copt124 = 1 / copt123;
  Real copt37  = strainclamped(1);
  Real copt45  = -copt37;
  Real copt47  = straincpp(1);
  Real copt48  = copt45 + copt47;
  Real copt52  = bspeed * copt48 * copt51;
  Real copt121 = -2 + BARRIERPOWER;
  Real copt127 = Power(copt50, 2);
  Real copt128 = 1 / copt127;
  Real copt63  = strainclamped(2);
  Real copt65  = -copt63;
  Real copt66  = straincpp(2);
  Real copt67  = copt65 + copt66;
  Real copt78  = bspeed * copt67 * copt77;
  Real copt131 = Power(copt75, 2);
  Real copt132 = 1 / copt131;
  Real copt82  = strainclamped(3);
  Real copt90  = -copt82;
  Real copt92  = straincpp(3);
  Real copt93  = copt90 + copt92;
  Real copt97  = bspeed * copt93 * copt96;
  Real copt135 = Power(copt95, 2);
  Real copt136 = 1 / copt135;
  Real copt102 = strainclamped(4);
  Real copt103 = -copt102;
  Real copt104 = straincpp(4);
  Real copt105 = copt103 + copt104;
  Real copt108 = bspeed * copt105 * copt107;
  Real copt141 = Power(copt106, 2);
  Real copt142 = 1 / copt141;
  Real copt111 = strainclamped(5);
  Real copt112 = -copt111;
  Real copt113 = straincpp(5);
  Real copt114 = copt112 + copt113;
  Real copt117 = bspeed * copt114 * copt116;
  Real copt145 = Power(copt115, 2);
  Real copt146 = 1 / copt145;
  out1(0) = BARRIERPOWER * bscale * bspeed * copt30 * Power(copt32, copt33);
  out1(1) = BARRIERPOWER * bscale * bspeed * copt51 * Power(copt52, copt33);
  out1(2) = BARRIERPOWER * bscale * bspeed * copt77 * Power(copt78, copt33);
  out1(3) = BARRIERPOWER * bscale * bspeed * copt96 * Power(copt97, copt33);
  out1(4) = BARRIERPOWER * bscale * bspeed * copt107 * Power(copt108, copt33);
  out1(5) = BARRIERPOWER * bscale * bspeed * copt116 * Power(copt117, copt33);
  out2(0, 0) = BARRIERPOWER * bscale * copt120 * copt124 *
               Power(copt32, copt121) * copt33;
  out2(0, 1) = 0;
  out2(0, 2) = 0;
  out2(0, 3) = 0;
  out2(0, 4) = 0;
  out2(0, 5) = 0;
  out2(1, 0) = 0;
  out2(1, 1) = BARRIERPOWER * bscale * copt120 * copt128 * copt33 *
               Power(copt52, copt121);
  out2(1, 2) = 0;
  out2(1, 3) = 0;
  out2(1, 4) = 0;
  out2(1, 5) = 0;
  out2(2, 0) = 0;
  out2(2, 1) = 0;
  out2(2, 2) = BARRIERPOWER * bscale * copt120 * copt132 * copt33 *
               Power(copt78, copt121);
  out2(2, 3) = 0;
  out2(2, 4) = 0;
  out2(2, 5) = 0;
  out2(3, 0) = 0;
  out2(3, 1) = 0;
  out2(3, 2) = 0;
  out2(3, 3) = BARRIERPOWER * bscale * copt120 * copt136 * copt33 *
               Power(copt97, copt121);
  out2(3, 4) = 0;
  out2(3, 5) = 0;
  out2(4, 0) = 0;
  out2(4, 1) = 0;
  out2(4, 2) = 0;
  out2(4, 3) = 0;
  out2(4, 4) = BARRIERPOWER * bscale * Power(copt108, copt121) * copt120 *
               copt142 * copt33;
  out2(4, 5) = 0;
  out2(5, 0) = 0;
  out2(5, 1) = 0;
  out2(5, 2) = 0;
  out2(5, 3) = 0;
  out2(5, 4) = 0;
  out2(5, 5) = BARRIERPOWER * bscale * Power(copt117, copt121) * copt120 *
               copt146 * copt33;

  return std::make_pair(hess, grad);
}
