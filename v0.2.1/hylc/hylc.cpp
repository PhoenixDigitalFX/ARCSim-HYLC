

#include "hylc.hpp"
#include "../blockvectors.hpp"

namespace hylc {

bool hylc_enabled() { return config.enabled; }

// bad design? store somewhere locally, but rest of arccsim is so global
std::shared_ptr<BaseMaterial> global_material = nullptr;
std::shared_ptr<BaseMaterial> get_material() {
  if (global_material == nullptr) {
    // if (config.material_type < 0)
    //   global_material = std::make_shared<AnalyticMaterial>(
    //       config.a0, config.a1, config.b0, config.b1);
    // else
    
    global_material = std::make_shared<FittedMaterial>(config.material_type);

    // // TEST PLOT EXTRAPOLATION
    // Vec6 straintst(0);
    // straintst(0) = 1.0;
    // straintst(2) = 1.0;

    // std::cout << "\n\n\n\n";
    // for (int i = 0; i < 100; ++i) {
    //   Vec6 straincopy = straintst;
    //   double a = i * 1.0 / 100;
    //   straincopy(0) = (1 - a) * 0.1 + a * 3.5;
    //   double psi = global_material->psi(straincopy);
    //   std::cout << straincopy(0) << ", " << psi << ", ";
    // }
    // std::cout << "\n\n\n\n";
    // for (int i = 0; i < 100; ++i) {
    //   Vec6 straincopy = straintst;
    //   double a = i * 1.0 / 100;
    //   straincopy(1) = (1 - a) * -10 + a * 10;
    //   double psi = global_material->psi(straincopy);
    //   std::cout << straincopy(1) << ", " << psi << ", ";
    // }
    // std::cout << "\n\n\n\n";
    // for (int i = 0; i < 100; ++i) {
    //   Vec6 straincopy = straintst;
    //   double a = i * 1.0 / 100;
    //   straincopy(5) = (1 - a) * -150 + a * 150;
    //   double psi = global_material->psi(straincopy);
    //   std::cout << straincopy(5) << ", " << psi << ", ";
    // }
    // std::cout << "\n\n\n\n";
    // for (int i = 0; i < 15; ++i) {
    //   double a = i * 1.0 / 15;
    //   for (int j = 0; j < 15; ++j) {
    //     Vec6 straincopy = straintst;
    //     double b = j * 1.0 / 15;
    //     straincopy(0) = (1 - a) * -0.1 + a * 0.5;
    //     straincopy(2) = (1 - b) * -0.1 + b * 0.5;
    //     // straincopy(5) = (1 - b) * -150 + b * 150;
    //     double psi = global_material->psi(straincopy);
    //     printf("%.10e, %.10e, %.10e, ", straincopy(0),straincopy(2), psi);
    //   }
    // }
    // std::cout << "\n\n\n\n";
  }
  return global_material;
}

double get_density() { return get_material()->density; }

Node *find_across(const Edge *edge, const Node *oppv1) {
  Node *found = nullptr;
  Vert *v     = edge_opp_vert(edge, 0);
  if (v)
    if (v->node == oppv1)  // wrong side
      v = edge_opp_vert(edge, 1);
  if (v)
    found = v->node;
  return found;
}

template <Space s>
void compute_local_information(const Face *face, Vec18 &xlocal, double &l0,
                               double &l1, double &l2, double &theta_rest0,
                               double &theta_rest1, double &theta_rest2,
                               bool &nn0_exists, bool &nn1_exists,
                               bool &nn2_exists) {
  xlocal         = Vec18(0);
  const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
             *n2 = face->v[2]->node;
  Vec3 x0 = pos<s>(n0), x1 = pos<s>(n1), x2 = pos<s>(n2);
  for (int j = 0; j < 3; ++j) {
    xlocal[0 + j] = x0[j];
    xlocal[3 + j] = x1[j];
    xlocal[6 + j] = x2[j];
  }

  // potential neighbors
  Node *nn0 = nullptr, *nn1 = nullptr, *nn2 = nullptr;
  l0 = 1, l1 = 1, l2 = 1;
  theta_rest0 = 0, theta_rest1 = 0, theta_rest2 = 0;
  {
    Edge *e = get_edge(n0, n1);
    nn0     = find_across(e, n2);  // vtx opposite of n2 across edge n0n1
    if (nn0) {
      Vec3 x = pos<s>(nn0);
      for (int j = 0; j < 3; ++j) xlocal[9 + j] = x[j];
    }
    l0          = e->l;
    theta_rest0 = e->theta_ideal;
  }
  {
    Edge *e = get_edge(n1, n2);
    nn1     = find_across(e, n0);
    if (nn1) {
      Vec3 x = pos<s>(nn1);
      for (int j = 0; j < 3; ++j) xlocal[12 + j] = x[j];
    }
    l1          = e->l;
    theta_rest1 = e->theta_ideal;
  }
  {
    Edge *e = get_edge(n2, n0);
    nn2     = find_across(e, n1);
    if (nn2) {
      Vec3 x = pos<s>(nn2);
      for (int j = 0; j < 3; ++j) xlocal[15 + j] = x[j];
    }
    l2          = e->l;
    theta_rest2 = e->theta_ideal;
  }

  nn0_exists = nn0 ? true : false;
  nn1_exists = nn1 ? true : false;
  nn2_exists = nn2 ? true : false;
}

template <Space s>
double hylc_local_energy(const Face *face) {
  namespace mm = hylc::mathematica;
  // 0. compute local primitive values
  Vec18 xlocal;
  double l0, l1, l2;
  double theta_rest0, theta_rest1, theta_rest2;
  bool nn0_exists, nn1_exists, nn2_exists;
  compute_local_information<s>(face, xlocal, l0, l1, l2, theta_rest0,
                               theta_rest1, theta_rest2, nn0_exists, nn1_exists,
                               nn2_exists);

  double A     = face->a;      // material space / reference config area
  Mat2x2 invDm = face->invDm;  // shape matrix
  Vec2 t0 = face->t0, t1 = face->t1, t2 = face->t2;

  t0 *= (double)nn0_exists;
  t1 *= (double)nn1_exists;
  t2 *= (double)nn2_exists;

  // 1. mmcpp compute epsilon(x),kappa(x) as vec6
  Vec6 strain;
  strain = mm::strain(xlocal, invDm, A, theta_rest0, theta_rest1, theta_rest2,
                      l0, l1, l2, t0, t1, t2);

  // std::cout<<"strain\n"<<strain<<"\n\n";
  // std::cout<<"t\n"<<t0<<"\n"<<t1<<"\n"<<t2<<"\n\n";
  // std::cout<<"l\n"<<l0<<"\n"<<l1<<"\n"<<l2<<"\n\n";
  // std::cout<<"A\n"<<A<<"\n";
  // 2. mmcpp compute psi(epskappa)
  double E = get_material()->psi(strain);
  // psi is energy density, constant over triangle
  E *= A;
  return E;

  // NOTE: probably not optimal recomputing normals within mmcpp of strains
  // same goes for dihedral angle theta
}

bool debug_print_strain_range =
    true;  // DEBUG print out running min max strains
double strain0a = 10, strain0b = 0, strain1a = 10, strain1b = -10,
       strain2a = 10, strain2b = 0;
double strain3a = 1e5, strain3b = -1e5, strain4a = 1e5, strain4b = -1e5,
       strain5a = 1e5, strain5b = -1e5;
template <Space s>
std::pair<Mat18x18, Vec18> hylc_local_forces(const Face *face) {
  namespace mm = hylc::mathematica;
  using namespace hylc::mathematica;
  // 0. compute local primitive values
  Vec18 xlocal;
  double l0, l1, l2;
  double theta_rest0, theta_rest1, theta_rest2;
  bool nn0_exists, nn1_exists, nn2_exists;
  compute_local_information<s>(face, xlocal, l0, l1, l2, theta_rest0,
                               theta_rest1, theta_rest2, nn0_exists, nn1_exists,
                               nn2_exists);

  double A     = face->a;      // material space / reference config area
  Mat2x2 invDm = face->invDm;  // shape matrix
  Vec2 t0 = face->t0, t1 = face->t1, t2 = face->t2;
  t0 *= (double)nn0_exists;
  t1 *= (double)nn1_exists;
  t2 *= (double)nn2_exists;

  double Acpy         = A;
  bool debug_localize = false;  // TODO also other methods
  Mat2x2 tile;
  tile(0, 0) = 0.032;
  tile(1, 0) = 0.0;
  tile(0, 1) = 0.0;
  tile(1, 1) = 0.02;
  if (debug_localize) {
    if (nn0_exists) {
      t0 = normalize(t0);
      l0 = norm(tile * t0);
    }
    if (nn1_exists) {
      t1 = normalize(t1);
      l1 = norm(tile * t1);
    }
    if (nn2_exists) {
      t2 = normalize(t2);
      l2 = norm(tile * t2);
    }
    Acpy = 1.0;  // sum
    // Acpy = 1.0/3.0; // avg
  }

  // std::cout<<"TRI INFO\n"<<l0<<" "<<norm(t0)<<", "<<l1<<" "<<norm(t1)<<",
  // "<<l2<<" "<<norm(t2)<<"; "<<A<<"\n\n";

  // 1. mmcpp compute epsilon, kappa as vec6 strain, and simulatenous grad and
  // hess
  std::tuple<std::vector<Mat18x18>, Mat6x18, Vec6> strain_hgv;
  strain_hgv = mm::strain_valdrv(xlocal, invDm, Acpy, theta_rest0, theta_rest1,
                                 theta_rest2, l0, l1, l2, t0, t1, t2);

  Vec6 &strain                       = std::get<2>(strain_hgv);
  Mat6x18 &strain_grad               = std::get<1>(strain_hgv);
  std::vector<Mat18x18> &strain_hess = std::get<0>(strain_hgv);

  // DEBUG not threadsafe but good enough for testing
  if (debug_print_strain_range) {
    strain0a = std::min(strain0a, strain[0]);
    strain0b = std::max(strain0b, strain[0]);
    strain1a = std::min(strain1a, strain[1]);
    strain1b = std::max(strain1b, strain[1]);
    strain2a = std::min(strain2a, strain[2]);
    strain2b = std::max(strain2b, strain[2]);
    strain3a = std::min(strain3a, strain[3]);
    strain3b = std::max(strain3b, strain[3]);
    strain4a = std::min(strain4a, strain[4]);
    strain4b = std::max(strain4b, strain[4]);
    strain5a = std::min(strain5a, strain[5]);
    strain5b = std::max(strain5b, strain[5]);
  }

  // 2. mmcpp compute dpsidek(ek),d2psidekdek(ek) at the same time
  std::pair<Mat6x6, Vec6> psidrv = get_material()->psi_drv(strain);
  // 4. assemble local grad: dekdx.T dpsidek
  Vec18 g = transpose(strain_grad) * psidrv.second;
  // 5. assemble local hess: dpsidek.T d2ekdxdx + dekdx.T d2psidekdek dekdx
  // second part 18x6 * 6x6 * 6x18
  Mat18x18 H = transpose(strain_grad) * psidrv.first * strain_grad;
  // first part 1x6 * 6x18x18
  for (int i = 0; i < 6; ++i) H += psidrv.second[i] * strain_hess[i];

  // DEBUG SPDify 18x18 matrix
  // TODO maybe deal with the boundary triangles with < 18 dof
  // make Eigen map to H.. cant, need to copy...
  typedef Eigen::Matrix<double, 18, 18> emat18;
  typedef Eigen::Matrix<double, 18, 1> evec18;
  emat18 eH;
  for (int i = 0; i < 18; ++i) {
    for (int j = 0; j < 18; ++j) {
      eH(i, j) = H(i, j);
    }
  }

  // construct SVD
  // compare time eigendecomp
  // BDCSVD "one of the fastest SVD algorithms"
  // Eigen::BDCSVD	<emat18> svd(eH, Eigen::ComputeThinU |
  // Eigen::ComputeThinV); JacobiSVD "slow but fast for small matrices" "<16"
  // Eigen::JacobiSVD<emat18, Eigen::NoQRPreconditioner> svd(
  //     eH, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // SelfAdjointEigenSolver is much faster apparently, robustness idk
  Eigen::SelfAdjointEigenSolver<emat18> slv;
  slv.compute(eH);  // 2x2 and 3x3 use computeDirect()

  // compute V and D
  // double emax = slv.eigenvalues().cwiseAbs().maxCoeff();
  // double emin = slv.eigenvalues().cwiseAbs().minCoeff();
  // printf("%.2e %.2e\n", emin, emax);
  double mineig = 1e-10;
  if (slv.eigenvalues()(0) < mineig) {
    evec18 D = slv.eigenvalues().cwiseMax(mineig);

    eH = slv.eigenvectors() * D.asDiagonal() * slv.eigenvectors().transpose();

    for (int i = 0; i < 18; ++i) {
      for (int j = 0; j < 18; ++j) {
        H(i, j) = eH(i, j);
      }
    }
  }
  // END DEBUG

  return std::make_pair(A * H, A * g);
}

template <Space s>
Vec18 hylc_local_forces_nojac(const Face *face) {
  namespace mm = hylc::mathematica;
  using namespace hylc::mathematica;
  // 0. compute local primitive values
  Vec18 xlocal;
  double l0, l1, l2;
  double theta_rest0, theta_rest1, theta_rest2;
  bool nn0_exists, nn1_exists, nn2_exists;
  compute_local_information<s>(face, xlocal, l0, l1, l2, theta_rest0,
                               theta_rest1, theta_rest2, nn0_exists, nn1_exists,
                               nn2_exists);

  double A     = face->a;      // material space / reference config area
  Mat2x2 invDm = face->invDm;  // shape matrix
  Vec2 t0 = face->t0, t1 = face->t1, t2 = face->t2;
  t0 *= (double)nn0_exists;
  t1 *= (double)nn1_exists;
  t2 *= (double)nn2_exists;

  // 1. mmcpp compute epsilon, kappa as vec6 ek, and simulatenous grad and hess
  std::tuple<Mat6x18, Vec6> strain_gv;
  strain_gv = mm::strain_valgrad(xlocal, invDm, A, theta_rest0, theta_rest1,
                                 theta_rest2, l0, l1, l2, t0, t1, t2);

  Vec6 &strain         = std::get<1>(strain_gv);
  Mat6x18 &strain_grad = std::get<0>(strain_gv);

  // DEBUG not threadsafe but good enough for testing
  if (debug_print_strain_range) {
    strain0a = std::min(strain0a, strain[0]);
    strain0b = std::max(strain0b, strain[0]);
    strain1a = std::min(strain1a, strain[1]);
    strain1b = std::max(strain1b, strain[1]);
    strain2a = std::min(strain2a, strain[2]);
    strain2b = std::max(strain2b, strain[2]);
    strain3a = std::min(strain3a, strain[3]);
    strain3b = std::max(strain3b, strain[3]);
    strain4a = std::min(strain4a, strain[4]);
    strain4b = std::max(strain4b, strain[4]);
    strain5a = std::min(strain5a, strain[5]);
    strain5b = std::max(strain5b, strain[5]);
  }

  Vec6 psi_grad = get_material()->psi_grad(strain);
  Vec18 g       = transpose(strain_grad) * psi_grad;

  return A * g;
}

template <int m, int n>
Mat<3, 3> submat3(const Mat<m, n> &A, int i, int j) {
  Mat3x3 Asub;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++) Asub(k, l) = A(i * 3 + k, j * 3 + l);
  return Asub;
}

template <int n>
Vec<3> subvec3(const Vec<n> &b, int i) {
  Vec3 bsub;
  for (int k = 0; k < 3; k++) bsub[k] = b[i * 3 + k];
  return bsub;
}

template <int m>
void add_submat(const Mat<m * 3, m * 3> &Asub, const Vec<m, int> &ix,
                SpMat<Mat3x3> &A) {
  for (int i = 0; i < m; i++)
    if (ix[i] >= 0)
      for (int j = 0; j < m; j++)
        if (ix[j] >= 0)
          A(ix[i], ix[j]) += submat3(Asub, i, j);
}

template <int m>
void add_subvec(const Vec<m * 3> &bsub, const Vec<m, int> &ix,
                std::vector<Vec3> &b) {
  for (int i = 0; i < m; i++)
    if (ix[i] >= 0)
      b[ix[i]] += subvec3(bsub, i);
}

template <int n>
Mat<1, n> rowmat(const Vec<n> &v) {
  Mat<1, n> A;
  for (int i = 0; i < n; i++) A(0, i) = v[i];
  return A;
}

template <int m, int n, int p, int q>
Mat<m * p, n * q> kronecker(const Mat<m, n> &A, const Mat<p, q> &B) {
  Mat<m * p, n * q> C;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
        for (int l = 0; l < q; l++) C(i * p + k, j * q + l) = A(i, j) * B(k, l);
  return C;
}

Vec<6, int> indices(const Node *n0, const Node *n1, const Node *n2,
                    const Node *nn0, const Node *nn1, const Node *nn2) {
  Vec<6, int> ix;
  ix[0] = n0->index;
  ix[1] = n1->index;
  ix[2] = n2->index;
  ix[3] = (nn0) ? nn0->index : -1;
  ix[4] = (nn1) ? nn1->index : -1;
  ix[5] = (nn2) ? nn2->index : -1;
  return ix;
}

Vec<3, int> indices(const Node *n0, const Node *n1, const Node *n2) {
  Vec<3, int> ix;
  ix[0] = n0->index;
  ix[1] = n1->index;
  ix[2] = n2->index;
  return ix;
}

}  // namespace hylc

template <Space s>
double hylc::hylc_internal_energy(const Cloth &cloth) {
  // local, parallel per triangle
  const Mesh &mesh = cloth.mesh;
  int n_triangles  = mesh.faces.size();
  std::vector<double> local_E(n_triangles);
#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    local_E[i] = hylc_local_energy<s>(mesh.faces[i]);
  }
  // global sequential / omp reduction
  double E = 0;
#pragma omp parallel for reduction(+ : E)
  for (size_t i = 0; i < local_E.size(); i++) {
    E += local_E[i];
  }

  return E;
}

template <Space s>
void hylc::hylc_add_internal_forces(const Cloth &cloth, SpMat<Mat3x3> &A,
                                    std::vector<Vec3> &b, double dt) {
  // local, parallel per triangle
  const Mesh &mesh = cloth.mesh;
  int n_triangles  = mesh.faces.size();
  std::vector<std::pair<Mat18x18, Vec18>> local_Hg(n_triangles);

  // std::chrono::high_resolution_clock::time_point t0 =
  //     std::chrono::high_resolution_clock::now();

#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    local_Hg[i] = hylc_local_forces<s>(mesh.faces[i]);
  }

  if (debug_print_strain_range) {
    printf(
        "strain: [%.2f, %.2f] [%.2f, %.2f] [%.2f, %.2f] [%.2f, %.2f] [%.2f, "
        "%.2f] [%.2f, %.2f]\n",
        strain0a, strain0b, strain1a, strain1b, strain2a, strain2b, strain3a,
        strain3b, strain4a, strain4b, strain5a, strain5b);
  }
  // int nt = std::chrono::duration_cast<std::chrono::milliseconds>(
  //              std::chrono::high_resolution_clock::now() - t0)
  //              .count();
  // printf("Time (ms): %d\n", nt);

  // global sequential assembly
  // sum force and hess with correct indexing
  for (int i = 0; i < n_triangles; ++i) {
    auto Hg = local_Hg[i];

    Mat18x18 &H = Hg.first;
    Vec18 &g    = Hg.second;

    const Face *face = mesh.faces[i];
    const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
               *n2 = face->v[2]->node;

    // try finding neighbors
    const Node *nn0 = find_across(get_edge(n0, n1), n2),
               *nn1 = find_across(get_edge(n1, n2), n0),
               *nn2 = find_across(get_edge(n2, n0), n1);

    // get velocities (missing neighbors set to 0)
    Vec18 vs;
    for (int j = 0; j < 3; ++j) {
      vs[3 * 0 + j] = n0->v[j];
      vs[3 * 1 + j] = n1->v[j];
      vs[3 * 2 + j] = n2->v[j];
      vs[3 * 3 + j] = (nn0) ? nn0->v[j] : 0;
      vs[3 * 4 + j] = (nn1) ? nn1->v[j] : 0;
      vs[3 * 5 + j] = (nn2) ? nn2->v[j] : 0;
    }
    // get indices (missing neighbors set to -1)
    auto ixs = hylc::indices(n0, n1, n2, nn0, nn1, nn2);
    // add to global matrix (missing neighbors with ix -1 ignored)
    if (dt == 0) {
      hylc::add_submat(H, ixs, A);   // -J
      hylc::add_subvec(-g, ixs, b);  // F
    } else {
      double damping = cloth.materials[face->label]->damping;  // or EDGEbased ?
      //-dt * (dt + damping) * J
      hylc::add_submat(dt * (dt + damping) * H, ixs, A);
      // dt * (F + (dt + damping) * J * vs)
      hylc::add_subvec(-dt * (g + (dt + damping) * H * vs), ixs, b);
    }
  }
}

template <Space s>
void hylc::hylc_add_internal_forces(const Cloth &cloth, std::vector<Vec3> &b,
                                    double dt, double stretchdamping) {
  typedef Mat<9, 9> Mat9x9;
  typedef Vec<9> Vec9;

  // local, parallel per triangle
  const Mesh &mesh = cloth.mesh;
  int n_triangles  = mesh.faces.size();
  std::vector<Vec18> local_g(n_triangles);
  std::vector<Vec9> local_Jfake_vs(n_triangles);

#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    local_g[i] = hylc_local_forces_nojac<s>(mesh.faces[i]);
  }

  if (debug_print_strain_range) {
    printf(
        "strain: [%.2f, %.2f] [%.2f, %.2f] [%.2f, %.2f] [%.2f, %.2f] [%.2f, "
        "%.2f] [%.2f, %.2f]\n",
        strain0a, strain0b, strain1a, strain1b, strain2a, strain2b, strain3a,
        strain3b, strain4a, strain4b, strain5a, strain5b);
  }

#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    auto &face = mesh.faces[i];
    Mat3x2 F   = derivative(pos<s>(face->v[0]->node), pos<s>(face->v[1]->node),
                          pos<s>(face->v[2]->node), face);
    Mat2x2 G   = (F.t() * F - Mat2x2(1)) / 2.;
    // NOTE dependency on non-hylc material params removed!
    // Vec4 k = stretching_stiffness(G,
    // (*::materials)[face->label]->stretching); double weakening =
    // (*::materials)[face->label]->weakening; k *= 1 / (1 + weakening *
    // face->damage);
    Mat2x3 D = derivative(face);
    Vec3 du = D.row(0), dv = D.row(1);
    Mat<3, 9> Du   = kronecker(rowmat(du), Mat3x3(1)),
              Dv   = kronecker(rowmat(dv), Mat3x3(1));
    const Vec3 &xu = F.col(0), &xv = F.col(1);  // should equal Du*mat_to_vec(X)
    Vec9 fuu = Du.t() * xu, fvv = Dv.t() * xv,
         fuv      = (Du.t() * xv + Dv.t() * xu) / 2.;
    Mat9x9 hess_e = (outer(fuu, fuu) + std::max(G(0, 0), 0.) * Du.t() * Du) +
                    (outer(fvv, fvv) + std::max(G(1, 1), 0.) * Dv.t() * Dv) +
                    (outer(fuu, fvv) + std::max(G(0, 0), 0.) * Dv.t() * Dv +
                     outer(fvv, fuu) + std::max(G(1, 1), 0.) * Du.t() * Du) +
                    2. * (outer(fuv, fuv));

    const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
               *n2    = face->v[2]->node;
    Vec9 vs           = mat_to_vec(Mat3x3(n0->v, n1->v, n2->v));
    local_Jfake_vs[i] = -face->a * hess_e * vs;
  }

  // if (debug_print_strain_range) {
  //   printf("strain: [%.2f, %.2f] [%.2f, %.2f] [%.2f, %.2f] [%.2f, %.2f]
  //   [%.2f, "
  //          "%.2f] [%.2f, %.2f]\n",
  //          strain0a, strain0b, strain1a, strain1b, strain2a, strain2b,
  //          strain3a, strain3b, strain4a, strain4b, strain5a, strain5b);
  // }

  // global sequential assembly
  // sum force and hess with correct indexing
  for (int i = 0; i < n_triangles; ++i) {
    Vec18 &g = local_g[i];

    const Face *face = mesh.faces[i];
    const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
               *n2 = face->v[2]->node;

    // try finding neighbors
    const Node *nn0 = find_across(get_edge(n0, n1), n2),
               *nn1 = find_across(get_edge(n1, n2), n0),
               *nn2 = find_across(get_edge(n2, n0), n1);

    // get indices (missing neighbors set to -1)
    auto ixs = hylc::indices(n0, n1, n2, nn0, nn1, nn2);
    // add to global matrix (missing neighbors with ix -1 ignored)
    if (dt == 0) {
      hylc::add_subvec(-g, ixs, b);  // F
    } else {
      hylc::add_subvec(-dt * g, ixs, b);  // dt * F
      // dt*(dt+damp)*J*vs, only with "fake" stretching J (from default arcsim)
      hylc::add_subvec(dt * ((dt + stretchdamping) * local_Jfake_vs[i]),
                       hylc::indices(n0, n1, n2), b);
    }
  }
}

template double hylc::hylc_internal_energy<WS>(const Cloth &);
template double hylc::hylc_internal_energy<PS>(const Cloth &);
template void hylc::hylc_add_internal_forces<WS>(const Cloth &, SpMat<Mat3x3> &,
                                                 std::vector<Vec3> &, double);
template void hylc::hylc_add_internal_forces<PS>(const Cloth &, SpMat<Mat3x3> &,
                                                 std::vector<Vec3> &, double);
template void hylc::hylc_add_internal_forces<WS>(const Cloth &,
                                                 std::vector<Vec3> &, double,
                                                 double);
template void hylc::hylc_add_internal_forces<PS>(const Cloth &,
                                                 std::vector<Vec3> &, double,
                                                 double);
