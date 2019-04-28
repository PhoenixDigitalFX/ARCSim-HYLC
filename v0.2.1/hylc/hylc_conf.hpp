#ifndef _HYLC_CONF_HPP_
#define _HYLC_CONF_HPP_

#include <json/json.h>

namespace hylc {

// default values will be used if not replaced through config file
struct Config {
  bool enabled = false; // false .. use arcsim material, true .. use hylc
  int material_type = 0; // 0 .. analytic material, 1 .. fitted data
  // bool eklinear = true; // toggle between linear theta vs. 2 tan(theta/2)
  // NOTE: ek is my shorthand notitation for epsilonkappa
  // which i use to describe the combined vector of the 6 dof of the
  // stretching tensor epsilon and the shape operator/curvature tensor kappa
};

extern Config config; // global struct, defined in cpp

void parse_hylc(const Json::Value &json);
} // namespace hylc

#endif /* end of include guard: _HYLC_CONF_HPP_ */
