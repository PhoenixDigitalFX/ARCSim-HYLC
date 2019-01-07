#include "hylc_conf.hpp"
using namespace hylc;

hylc::Config hylc::config{};

extern void parse(double &x, const Json::Value &json);
extern void parse(int &x, const Json::Value &json);
extern void parse(bool &x, const Json::Value &json);

template <typename T> void parse_optional(T &x, const Json::Value &json) {
  if (json.isNull())
    return; // leave x as its default value
  else
    parse(x, json);
}

void parse(hylc::Config &params, const Json::Value &json) {
  parse_optional(params.enabled, json["enabled"]);
  parse_optional(params.a0, json["a0"]);
  parse_optional(params.a1, json["a1"]);
  parse_optional(params.b0, json["b0"]);
  parse_optional(params.b1, json["b1"]);
  parse_optional(params.material_type, json["type"]);
}

void hylc::parse_hylc(const Json::Value &json) { parse(hylc::config, json); }
