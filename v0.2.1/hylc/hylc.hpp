#ifndef _HYLC_H_
#define _HYLC_H_

#include "../cloth.hpp"
#include <utility>
#include <vector>

#include "../sparse.hpp"

// to turn off remeshing
// -> "disable": ["remeshing"] in config files

#include "hylc_conf.hpp"
#include "hylc_mm.hpp"

#include <iostream>

#define USE_SPARSE3

namespace hylc {

bool hylc_enabled() {
  return material.enabled;
}

Node *find_across(const Edge *edge, const Node *oppv1) {
  Node *found = nullptr;
  Vert *v = edge_opp_vert(edge, 0);
  if (v)
    if (v->node == oppv1) // wrong side
      v = edge_opp_vert(edge, 1);
  if (v)
    found = v->node;
  return found;
}

template <Space s>
void compute_local_information(const Face *face, Vec18 &xlocal, double &l0,
                               double &l1, double &l2, double &theta_rest0,
                               double &theta_rest1, double &theta_rest2, bool &nn0_exists, bool& nn1_exists, bool&nn2_exists) {
  xlocal = Vec18(0);
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
    nn0 = find_across(e, n2); // vtx opposite of n2 across edge n0n1
    if (nn0) {
      Vec3 x = pos<s>(nn0);
      for (int j = 0; j < 3; ++j)
        xlocal[9 + j] = x[j];
    }
    l0 = e->l;
    theta_rest0 = e->theta_ideal;
  }
  {
    Edge *e = get_edge(n1, n2);
    nn1 = find_across(e, n0);
    if (nn1) {
      Vec3 x = pos<s>(nn1);
      for (int j = 0; j < 3; ++j)
        xlocal[12 + j] = x[j];
    }
    l1 = e->l;
    theta_rest1 = e->theta_ideal;
  }
  {
    Edge *e = get_edge(n2, n0);
    nn2 = find_across(e, n1);
    if (nn2) {
      Vec3 x = pos<s>(nn2);
      for (int j = 0; j < 3; ++j)
        xlocal[15 + j] = x[j];
    }
    l2 = e->l;
    theta_rest2 = e->theta_ideal;
  }

  nn0_exists = nn0 ? true:false;
  nn1_exists = nn1 ? true:false;
  nn2_exists = nn2 ? true:false;
}

template <Space s> double hylc_local_energy(const Face *face) {
  namespace mm = hylc::mathematica;
  // 0. compute local primitive values
  Vec18 xlocal;
  double l0, l1, l2;
  double theta_rest0, theta_rest1, theta_rest2;
  bool nn0,nn1,nn2;
  compute_local_information<s>(face, xlocal, l0, l1, l2, theta_rest0,
                               theta_rest1, theta_rest2, nn0,nn1,nn2);

  double A = face->a;         // material space / reference config area
  Mat2x2 invDm = face->invDm; // shape matrix
  Vec2 t0 = face->t0, t1 = face->t1, t2 = face->t2;

  t0 *= (double) nn0;
  t1 *= (double) nn1;
  t2 *= (double) nn2;

  // 1. mmcpp compute epsilon(x),kappa(x) as vec6
  // NOTE not using theta_ideal, TODO
  Vec6 ek = mm::ek(xlocal, invDm, A, l0, l1, l2, t0, t1, t2);
  // std::cout<<"ek\n"<<ek<<"\n\n";
  // std::cout<<"t\n"<<t0<<"\n"<<t1<<"\n"<<t2<<"\n\n";
  // std::cout<<"l\n"<<l0<<"\n"<<l1<<"\n"<<l2<<"\n\n";
  // std::cout<<"A\n"<<A<<"\n";
  // 2. mmcpp compute psi(epskappa)
  double E = mm::psi(ek, material.a0, material.a1, material.b0, material.b1);
  return E;

  // NOTE: not optimal recomputing normals within mmcpp of strains
  // same goes for dihedral angle theta
}

template <Space s> double hylc_internal_energy(const Cloth &cloth) {
  // local, parallel per triangle
  const Mesh &mesh = cloth.mesh;
  int n_triangles = mesh.faces.size();
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
std::pair<Mat18x18, Vec18> hylc_local_forces(const Face *face) {
  namespace mm = hylc::mathematica;
  using namespace hylc::mathematica;
  // 0. compute local primitive values
  Vec18 xlocal;
  double l0, l1, l2;
  double theta_rest0, theta_rest1, theta_rest2;
  bool nn0,nn1,nn2;
  compute_local_information<s>(face, xlocal, l0, l1, l2, theta_rest0,
                               theta_rest1, theta_rest2,  nn0,nn1,nn2);

  double A = face->a;         // material space / reference config area
  Mat2x2 invDm = face->invDm; // shape matrix
  Vec2 t0 = face->t0, t1 = face->t1, t2 = face->t2;
  t0 *= (double) nn0;
  t1 *= (double) nn1;
  t2 *= (double) nn2;

  // 1. mmcpp compute epsilon, kappa as vec6 ek
  // TODO could also compute at the same time as derivatives and avoid some
  // double computation
  Vec6 ek = mm::ek(xlocal, invDm, A, l0, l1, l2, t0, t1, t2);
  // 2. mmcpp compute dpsidek(ek),d2psidekdek(ek) at the same time
  std::pair<Mat6x6, Vec6> psidrv =
      mm::psi_drv(ek, material.a0, material.a1, material.b0, material.b1);
  // 3. mmcpp compute dekdx and d2ekdxdx at the same time
  std::pair<std::vector<Mat18x18>, Mat6x18> ekdrv =
      mm::ek_drv(xlocal, invDm, A, l0, l1, l2, t0, t1, t2);
  // 4. assemble local grad: dekdx.T dpsidek
  Vec18 g = transpose(ekdrv.second) * psidrv.second;
  // 5. assemble local hess: dpsidek.T d2ekdxdx + dekdx.T d2psidekdek dekdx
  // second part 18x6 * 6x6 * 6x18
  Mat18x18 H = transpose(ekdrv.second) * psidrv.first * ekdrv.second;
  // first part 1x6 * 6x18x18
  for (int i = 0; i < 6; ++i)
    H += psidrv.second[i] * ekdrv.first[i];

  return std::make_pair(H, g);
}

// Vec<3, int> indices(const Node *n0, const Node *n1, const Node *n2) {
//   Vec<3, int> ix;
//   ix[0] = n0->index;
//   ix[1] = n1->index;
//   ix[2] = n2->index;
//   return ix;
// }

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

template <Space s>
void hylc_add_internal_forces(const Cloth &cloth, SpMat<Mat3x3> &A,
                              std::vector<Vec3> &b, double dt) {
  // local, parallel per triangle
  const Mesh &mesh = cloth.mesh;
  int n_triangles = mesh.faces.size();
  std::vector<std::pair<Mat18x18, Vec18>> local_Hg(n_triangles);
#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    local_Hg[i] = hylc_local_forces<s>(mesh.faces[i]);
  }
  // global sequential assembly
  // sum force and hess with correct indexing
  for (int i = 0; i < n_triangles; ++i) {
    auto Hg = local_Hg[i];

    Mat18x18 &H = Hg.first;
    Vec18 &g = Hg.second;

    // H = H * 0.0;
    // g = g * 0.0;

    // std::cout<<"H\n"<<H<<"\n\n";
    // std::cout<<"g\n"<<g<<"\n\n";
    // TODO ASK RAHUL
    // which damping to use, face or some adjf stuff
    // can i assume that e0 goes from x0 to x1 etc. ? <- give him localprimcode
    // should i ignore dihedral and normal already in arcsim, or take it out,
    // bc mm recomputes it
    // why physics static materials pointer..? do i need that?
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
      hylc::add_submat(H, ixs, A);  // -J
      hylc::add_subvec(-g, ixs, b); // F
    } else {
      double damping = cloth.materials[face->label]->damping; // or EDGEbased ?
      //-dt * (dt + damping) * J
      hylc::add_submat(dt * (dt + damping) * H, ixs, A);
      // dt * (F + (dt + damping) * J * vs)
      hylc::add_subvec(-dt * (g + (dt + damping) * H * vs), ixs, b);
    }
  }
}
} // namespace hylc
#endif /* end of include guard: _HYLC_H_ */
