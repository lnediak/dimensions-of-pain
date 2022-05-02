#ifndef HYPERVOXEL_PHYSICS_CONVEX_HPP_
#define HYPERVOXEL_PHYSICS_CONVEX_HPP_

#include <vector>

#include "util.hpp"
#include "vector.hpp"

namespace hypervoxel {

template <std::size_t N, class Attr> struct HalfSpace {
  v::DVec<N> n;
  double t;
  Attr attr;

  bool contains(const v::DVec<N> &p) const { return v::dot(p, n) <= t; }
};

/// represents a convex polytope in N dimensions
template <std::size_t N, class Attr> struct Polytope {
  std::vector<HalfSpace<N, Attr>> halfSpaces;

  bool contains(const v::DVec<N> &p) const {
    for (const HalfSpace &h : halfSpaces) {
      if (!h.contains(p)) {
        return false;
      }
    }
    return true;
  }

private:
  static bool isSmol(double d) { return -1e-8 < d && d < 1e-8; }

  static void combine2d(const std::vector<v::DVec<3>> &l0,
                        const std::vector<v::DVec<3>> &r0,
                        const std::vector<v::DVec<3>> &l1,
                        const std::vector<v::DVec<3>> &r1,
                        const std::vector<v::DVec<3>> &l2,
                        const std::vector<v::DVec<3>> &r2) {
    // TODO: WRITE
  }

  std::vector<v::DVec<4>> tmp;

public:
  template <class Fun> writeTriangles(const SliceDirs &sd, Fun &&fun) const {
    tmp.clear();
    tmp.resize(halfSpaces.size());
    for (auto &half : halfSpaces) {
      v::DVec<3> n3 = {v::dot(sd.r, half.n), v::dot(sd.u, half.n),
                       v::dot(sd.f, half.n)};
      double t3 = half.t - v::dot(sd.c, half.n);
      // I don't actually need to normalize it, but I don't like smol numbers
      double norm1 = v::norm1(n3);
      if (norm1 < 1e-8) {
        return;
      }
      n3 /= norm1;
      t3 /= norm1;
      tmp.push_back({n3[0], n3[1], n3[2], t3});
    }
    for (auto it = tmp.begin(), ite = tmp.end(); it != ite; ++it) {
      v::DVec<3> n3 = {(*it)[0], (*it)[1], (*it)[2]};
      double t3 = (*it)[3];
      v::DVec<3> a = {1, 0, 0};
      if (isSmol(n3[1]) && isSmol(n3[2])) {
        a[1] = 1;
      }
      double nor2 = v::norm2(n3);
      a -= (v::dot(a, n3) / nor2) * n3;
      v::DVec<3> b = cross3(a, n3);
      v::DVec<3> c = (t3 / nor2) * n3;
      // process da 2d face here
    }
  }
};

} // namespace hypervoxel

#endif // HYPERVOXEL_PHYSICS_CONVEX_HPP_

