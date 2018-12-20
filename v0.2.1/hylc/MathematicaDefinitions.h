#ifndef _MATHEMATICADEFINITIONS_H_
#define _MATHEMATICADEFINITIONS_H_

// #include "EigenDefinitions.h"
#include <math.h>

namespace hylc {
namespace mathematica {

typedef double Real;
inline Real Power(const Real &A, const Real &B) { return std::pow(A, B); }
// inline Real Power(const Real &A, const int &B) { return std::pow(A, B); }
inline Real Power(const Real &A, const int &B) {
  assert(abs(B) < 8 && "require higher pow, implement using std.");
  Real out = 1.0;
  if (B >= 0) {
    for (int i = 0; i < B; ++i)
      out *= A;
  } else {
    Real invA = 1.0 / A;
    for (int i = 0; i < -B; ++i)
      out *= invA;
  }
  return out;
}
inline Real Abs(const Real &A) { return std::abs(A); }
inline Real Sqrt(const Real &A) { return std::sqrt(A); }
inline Real Exp(const Real &A) { return std::exp(A); }
inline Real Sin(const Real &A) { return std::sin(A); }
inline Real Cos(const Real &A) { return std::cos(A); }
inline Real ArcCos(const Real &A) { return std::acos(A); }
inline Real ArcTan(const Real &X, const Real &Y) { return std::atan2(Y, X); }

} // namespace Mathematica
} // namespace HYLC

/*
OptimizeExpressionToC[expr_, tuple_: False, sum_: False,
   pow2_: True] :=
  Module[{optimizedExpr, mainExpr, n, m, defs, output},
   optimizedExpr = Experimental`OptimizeExpression[expr];
   If[ToString@optimizedExpr[[1, 0]] ==
     "Block", {n = Length[optimizedExpr[[1, 1]]];
     mainExpr = optimizedExpr[[1, 2, n + 1]];}, {n = 0;
     mainExpr = Flatten@{optimizedExpr[[1]]};}];
   m = Length[mainExpr];
   defs =
    Table["Real " <> ToString@CForm@optimizedExpr[[1, 2, i, 1]] <>
      " = " <> ToString@CForm@optimizedExpr[[1, 2, i, 2]] <>
      ";\n", {i, 1, n}];
   output = If[tuple,
     Flatten[
      Table[Flatten@
        MapIndexed[
         "out" <> ToString[i] <> "(" <>
           StringJoin@Riffle[ToString /@ (#2 - 1), ","] <> ") = " <>
           ToString@CForm@#1 <> ";\n" &,
         mainExpr[[i]], {ArrayDepth[mainExpr[[i]]]}], {i, 1,
        Length[mainExpr]}]],
     Flatten@
      MapIndexed[
       "out(" <> StringJoin@Riffle[ToString /@ (#2 - 1), ","] <>
         ") = " <> ToString@CForm@#1 <> ";\n" &,
       mainExpr, {ArrayDepth[mainExpr]}]
     ];
   output = Join[defs, output];
   output = StringReplace[output, "Power(E," -> "Exp("]; (*
   Rename exponent function *)

   output = StringReplace[output, "Compile_$" -> "c"]; (*
   Rename temporary variables *)

   output =
    If[pow2,
     StringReplace[output,
      RegularExpression["Power\\(([^\\)\\(,\\s]*),2\\)"] -> "$1*$1"],
     output]; (* Replace Power(v,2) with multiplication v*v *)

   output =
    If[sum, StringReplace[output,
      RegularExpression["(out\\([^\\)]*\\)) (=)"] -> "$1 +="],
     output]; (* out(..) += instead of out(..) = .. *)

   StringJoin[output]];
*/

#endif /* end of include guard: _MATHEMATICADEFINITIONS_H_ */
