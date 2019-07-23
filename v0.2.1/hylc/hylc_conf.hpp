#ifndef _HYLC_CONF_HPP_
#define _HYLC_CONF_HPP_

#include <json/json.h>
#include <memory>
#include "materials/BaseMaterial.hpp"

namespace hylc {

// default values will be used if not replaced through config file
struct Config {
  bool enabled = false; // false .. use arcsim material, true .. use hylc
  std::shared_ptr<BaseMaterial> material = nullptr;
  double seam_stiffness = 0.0;
};

extern Config config; // global struct, defined in cpp, meh.. should try and use arcsim stuff instead.
// like actually make my material a replacement for arcsim material


struct Debug {
  int debug_color = 0; // cycle through colorcoded strain rendering (0 = none) 
};
extern Debug debug;

void parse_hylc(const Json::Value &json);
} // namespace hylc

#endif /* end of include guard: _HYLC_CONF_HPP_ */
