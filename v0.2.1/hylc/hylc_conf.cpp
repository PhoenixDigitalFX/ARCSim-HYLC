#include "hylc_conf.hpp"
#include <fstream>
#include "materials/SplineMaterial.hpp"

using namespace hylc;

hylc::Config hylc::config{};
hylc::Debug hylc::debug{};

extern void complain(const Json::Value &json, const std::string &expected);
extern void parse(double &x, const Json::Value &json);
extern void parse(int &x, const Json::Value &json);
extern void parse(bool &x, const Json::Value &json);
extern void parse(std::string &x, const Json::Value &json);

template <typename T>
void parse(T &x, const Json::Value &json, const T &x0) {
  if (json.isNull())
    x = x0;
  else
    parse(x, json);
}

void parse(SplineMaterial::Spline1D &spline1d, const Json::Value &json);
void parse(SplineMaterial::Spline2D &spline2d, const Json::Value &json);
void parse(SplineMaterial::HSpline2Das1D &, const Json::Value &);

template <int n>
void parse(Vec<n> &v, const Json::Value &json) {
  if (!json.isArray())
    complain(json, "array");
  assert(json.size() == n);
  for (int i = 0; i < n; i++) parse(v[i], json[i]);
}
template <typename T>
void parse(std::vector<T> &v, const Json::Value &json) {
  if (!json.isArray())
    complain(json, "array");
  v.resize(json.size());
  for (int i = 0; i < json.size(); i++) parse(v[i], json[i]);
}

template <typename T>
void parse_optional(T &x, const Json::Value &json) {
  if (json.isNull())
    return;  // leave x as its default value
  else
    parse(x, json);
}

void parse(SplineMaterial::Spline1D &spline1d, const Json::Value &json) {
  std::vector<double> tx, c;
  int degree = 3;

  parse(spline1d.k, json["k"]);
  parse(tx, json["tx"]);
  parse(c, json["c"]);

  spline1d.spline = std::make_shared<fitpackpp::BSplineCurve>(tx, c, degree);

  parse(spline1d.clampx, json["clampx"]);
  if (spline1d.clampx) {
    parse(spline1d.xmin, json["xmin"]);
    parse(spline1d.xmax, json["xmax"]);
  }

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

  // TODO clamp

  spline2d.spline =
      std::make_shared<fitpackpp::BSplineSurface>(tx, ty, c, degree);
}

void parse(Poly1D &poly, const Json::Value &json) {
  parse(poly.c, json["k"]);

  parse(poly.compr, json["compr"]);
  parse(poly.clampx, json["clampx"]);
  if (poly.clampx) {
    parse(poly.xmin, json["xmin"]);
    parse(poly.xmax, json["xmax"]);
  }
}

void parse(HermiteSpline1D &spline, const Json::Value &json) {
  parse(spline.k, json["k"]);
  parse(spline.t, json["t"]);
  parse(spline.p, json["p"]);
  parse(spline.m, json["m"]);
  parse(spline.ext, json["ext"], 0);
}
void parse(HermiteSpline2D &spline, const Json::Value &json) {
  parse(spline.k0, json["k0"]);
  parse(spline.k1, json["k1"]);
  parse(spline.tu, json["tu"]);
  parse(spline.tv, json["tv"]);
  parse(spline.p, json["p"]);
  parse(spline.mu, json["mu"]);
  parse(spline.mv, json["mv"]);
  if (!json["muv"].isNull())
    parse(spline.muv, json["muv"]);
  else
    spline.muv.resize(spline.mu.size(), 0);
  parse(spline.ext, json["ext"], 0);
}

void parse(SplineMaterial::HSpline2Das1D &spline, const Json::Value &json) {
  parse(spline.k0, json["k0"]);
  parse(spline.k1, json["k1"]);
  parse(spline.fun.t, json["t"]);
  parse(spline.fun.p, json["p"]);
  parse(spline.fun.m, json["m"]);
  parse(spline.fun.ext, json["ext"], 0);

  // parse(spline.fun.c, json["c"]); // DEBUG POLY1D INSTEAD
  // spline.fun.c[0] *= 0.5;
  // spline.fun.compr = false;
}

void parse(Poly2D &poly, const Json::Value &json) {
  parse(poly.k0, json["k0"]);
  parse(poly.k1, json["k1"]);
  parse(poly.c, json["C"]);

  parse(poly.clampx, json["clampx"]);
  if (poly.clampx) {
    parse(poly.xmin, json["xmin"]);
    parse(poly.xmax, json["xmax"]);
  }
  parse(poly.clampy, json["clampy"]);
  if (poly.clampy) {
    parse(poly.ymin, json["ymin"]);
    parse(poly.ymax, json["ymax"]);
  }
}

#include <random>
void testmat(std::shared_ptr<SplineMaterial> mat) {
  // TEST PLOT EXTRAPOLATION
  Vec6 straintst(0);
  straintst(0) = 1.0;
  straintst(2) = 1.0;

  // // std::cout << "\n\n\n\n";
  // // for (int i = 0; i < 100; ++i) {
  // //   Vec6 straincopy = straintst;
  // //   double a = i * 1.0 / 100;
  // //   straincopy(0) = (1 - a) * 0.1 + a * 3.5;
  // //   double psi = global_material->psi(straincopy);
  // //   std::cout << straincopy(0) << ", " << psi << ", ";
  // // }
  // // std::cout << "\n\n\n\n";
  // // for (int i = 0; i < 100; ++i) {
  // //   Vec6 straincopy = straintst;
  // //   double a = i * 1.0 / 100;
  // //   straincopy(1) = (1 - a) * -10 + a * 10;
  // //   double psi = global_material->psi(straincopy);
  // //   std::cout << straincopy(1) << ", " << psi << ", ";
  // // }
  // // std::cout << "\n\n\n\n";
  // // for (int i = 0; i < 100; ++i) {
  // //   Vec6 straincopy = straintst;
  // //   double a = i * 1.0 / 100;
  // //   straincopy(5) = (1 - a) * -150 + a * 150;
  // //   double psi = global_material->psi(straincopy);
  // //   std::cout << straincopy(5) << ", " << psi << ", ";
  // // }
  // // std::cout << "\n\n\n\n";
  // int k0=0,k1=2,n=25;
  // for (int i = 0; i < n; ++i) {
  //   double a = i * 1.0 / n;
  //   for (int j = 0; j < n; ++j) {
  //     Vec6 straincopy = straintst;
  //     double b = j * 1.0 / n;
  //     straincopy(k0) = (1 - a) * 0.3 + a * 1.2;
  //     straincopy(k1) = (1 - b) * -0.2 + b * 0.0;
  //     straincopy(k0) += 1;
  //     straincopy(k1) += 1;
  //     // straincopy(k0) = (1 - a) * -160 + a * 150;
  //     // straincopy(k1) = (1 - b) * -200 + b * 200;
  //     double psi = material->psi(straincopy);
  //     printf("%.10e, %.10e, %.10e, ", straincopy(k0),straincopy(k1), psi);
  //   }
  // }
  // std::cout << "\n\n\n\n";

  // TEST GRADIENTS
  // std::mt19937 rng(1991);
  // std::uniform_real_distribution<> rnd(-0.5, 0.5);

  // std::cout << "\n\n\n\n<clip>";
  // for (float p = -9; p < -1.01; p += 0.1f) {
  //   Vec6 x;
  //   for (int i = 0; i < 6; i++) x(i) = rnd(rng) * (i < 3 ? 1 : 200);
  //   x(0) += 1.0;
  //   x(2) += 1.0;

  //   Vec6 dx;
  //   for (int i = 0; i < 6; i++)
  //     // dx(i) = rnd(rng) * (i < 3 ? 1 : 200);
  //     dx(i) = rnd(rng);

  //   // for (int i = 3; i < 6; i++) {
  //   //   dx(i) = 0;
  //   //   x(i)  = 0;
  //   // }
  //   dx = dx * 1.0 / norm(dx);

  //   double eps = std::pow(10, p);

  //   // psi (s + eps ds) - [psi(s) + eps ds dpsids(s)] ~~ O(eps^2)
  //   double psi0 = mat->psi(x);
  //   double psi1 = mat->psi(x + eps * dx);

  //   std::pair<Mat6x6, Vec6> drv = mat->psi_drv(x);
  //   double epsdsdpsids          = eps * dot(dx, drv.second);
  //   double err                  = std::abs(psi0 + epsdsdpsids - psi1);
  //   printf("%.2e, %.10e, ", eps, err);
  // }
  // std::cout << "<clip>\n\n\n\n";
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
  // material->density *= 2.0;
  // TODO damping (mass,stretch for expl and stiffness for impl)

  parse(material->strainscale, json["strain scale"]);
  parse(material->strainshift, json["strain shift"]);
  // parse(material->strain_min, json["strain min"]);
  // parse(material->strain_max, json["strain max"]);

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
  // parse(material->polys_1d, jsoncoeff["1D"]);
  parse(material->hsplines_1d, jsoncoeff["1D"]);
  // parse(material->splines_1d, jsoncoeff["1D"]);
  // parse(material->splines_2d, jsoncoeff["2D"]); // deprecated
  // parse(material->polys_2d, jsoncoeff["2D"]);
  parse(material->hsplines_2d, jsoncoeff["2D"]);
  // parse(material->hsplines_2d1d, jsoncoeff["2D"]);

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

  testmat(material);

  return material;
}

void parse(hylc::Config &params, const Json::Value &json) {
  parse_optional(params.enabled, json["enabled"]);

  std::string filename;
  parse(filename, json["material"]);

  params.material = load_material(filename);
  parse(config.seam_stiffness, json["seam_stiffness"], 0.0);
}

void hylc::parse_hylc(const Json::Value &json) { parse(hylc::config, json); }
