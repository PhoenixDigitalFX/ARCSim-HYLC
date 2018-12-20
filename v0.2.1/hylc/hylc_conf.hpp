#ifndef _HYLC_CONF_HPP_
#define _HYLC_CONF_HPP_

namespace hylc {

// default values will be used if not replaced through config file
struct Material {
  bool enabled = false;
  double a0 = 5e-2;
  double a1 = 5e-2;
  double b0 = 1e-2;
  double b1 = 1e-2;
};

extern Material material;
}

#endif /* end of include guard: _HYLC_CONF_HPP_ */
