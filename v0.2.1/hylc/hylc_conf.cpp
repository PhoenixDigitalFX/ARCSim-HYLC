#include "hylc_conf.hpp"
#include "materials/SplineMaterial.hpp"
#include <fstream>

using namespace hylc;

hylc::Config hylc::config{};

extern void complain(const Json::Value &json, const std::string &expected);
extern void parse(double &x, const Json::Value &json);
extern void parse(int &x, const Json::Value &json);
extern void parse(bool &x, const Json::Value &json);
extern void parse(std::string &x, const Json::Value &json);


void parse(SplineMaterial::Spline1D &spline1d, const Json::Value &json);
void parse(SplineMaterial::Spline2D &spline2d, const Json::Value &json);

template <int n> void parse(Vec<n> &v, const Json::Value &json) {
  if (!json.isArray())
    complain(json, "array");
  assert(json.size() == n);
  for (int i = 0; i < n; i++)
    parse(v[i], json[i]);
}
template <typename T> void parse(std::vector<T> &v, const Json::Value &json) {
  if (!json.isArray())
    complain(json, "array");
  v.resize(json.size());
  for (int i = 0; i < json.size(); i++)
    parse(v[i], json[i]);
}

template <typename T> void parse_optional(T &x, const Json::Value &json) {
  if (json.isNull())
    return; // leave x as its default value
  else
    parse(x, json);
}


void parse(SplineMaterial::Spline1D &spline1d, const Json::Value &json) {
  std::vector<double> tx, c;
  int degree = 3;

  parse(spline1d.k, json["k"]);
  parse(tx, json["tx"]);
  parse(c, json["c"]);

  spline1d.spline = std::make_shared<fitpackpp::BSplineCurve>(tx,c,degree);


  // printf("k %d\n",spline1d.k);
  // for(double a = -1; a <= 1; a+= 0.333333)
  // {
  //   printf("%f %f\n",a, spline1d.spline->eval(a));
  // }
  
}

void parse(SplineMaterial::Spline2D &spline2d, const Json::Value &json) {
  std::vector<double> tx, ty, c;
  int degree = 3;

  parse(spline2d.k0, json["k0"]);
  parse(spline2d.k1, json["k1"]);
  parse(tx, json["tx"]);
  parse(ty, json["ty"]);
  parse(c, json["c"]);

  spline2d.spline = std::make_shared<fitpackpp::BSplineSurface>(tx,ty,c,degree);
}

void parse(Poly2D &poly, const Json::Value &json) {
  parse(poly.k0, json["k0"]);
  parse(poly.k1, json["k1"]);
  parse(poly.c, json["C"]);
}


std::shared_ptr<SplineMaterial> load_material(const std::string &filename) {
  Json::Value json;
  Json::Reader reader;
  std::ifstream file(filename.c_str());
  bool parsingSuccessful = reader.parse(file, json);
  if (!parsingSuccessful) {
    fprintf(stderr, "Error reading file: %s\n", filename.c_str());
    fprintf(stderr, "%s", reader.getFormatedErrorMessages().c_str());
    abort();
  }
  file.close();

  auto material = std::make_shared<SplineMaterial>();

  parse(material->density, json["density"]);
  // TODO damping (mass,stretch for expl and stiffness for impl)

  parse(material->strainscale, json["strain scale"]);
  parse(material->strainshift, json["strain shift"]);
  parse(material->strain_min, json["strain min"]);
  parse(material->strain_max, json["strain max"]);

  // for(int i = 0; i < 6; i++)
  //   printf("%.2f  ", material->strainscale(i));
  // printf("\n");
  // for(int i = 0; i < 6; i++)
  //   printf("%.2f  ", material->strainshift(i));
  // printf("\n");
  // for(int i = 0; i < 6; i++)
  //   printf("%.2f  ", material->strain_min(i));
  // printf("\n");
  // for(int i = 0; i < 6; i++)
  //   printf("%.2f  ", material->strain_max(i));
  // printf("\n");

  Json::Value jsoncoeff = json["coeffs"];
  parse(material->C0, jsoncoeff["const"]);
  parse(material->splines_1d, jsoncoeff["1D"]);
  // parse(material->splines_2d, jsoncoeff["2D"]); // deprecated
  parse(material->polys_2d, jsoncoeff["2D"]);

  // for(double a = 0.3; a <= 2.0; a+= 0.1888888888888)
  // {
  //   Vec6 s(0);
  //   s(0) = 1.0;
  //   s(2) = 1.0;
  //   s(0) = a;
  //   // printf("%f %f\n",a, material->psi(s));
  //   // printf("%f %f\n",a, material->psi_grad(s)(0));
  //   printf("%f %f\n",a, material->TST(s));
  // }



  // for(double a = -60; a <= 60.0; a+= 13.3333333333333)
  // {
  //   Vec6 s(0);
  //   s(0) = 1.0;
  //   s(2) = 1.0;
  //   s(3) = a;
  //   // printf("%f %f\n",a, material->psi(s));
  //   // printf("%f %f\n",a, material->psi_grad(s)(3));
  //   printf("%f %f\n",a, material->TST(s));
  // }


  material->initialized = true;
  return material;
}


void parse(hylc::Config &params, const Json::Value &json) {
  parse_optional(params.enabled, json["enabled"]);

  std::string filename;
  parse(filename, json["material"]);

  params.material = load_material(filename);
}

void hylc::parse_hylc(const Json::Value &json) { parse(hylc::config, json); }
