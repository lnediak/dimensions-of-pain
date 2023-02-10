#ifndef CONVEX_HPP_
#define CONVEX_HPP_

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <vector>

#include "util.hpp"
#include "vector.hpp"

template <std::size_t N, class Attr> struct HalfSpace {
  v::DVec<N> n;
  double t;
  Attr attr;

  bool contains(const v::DVec<N> &p) const { return v::dot(p, n) <= t; }
};

bool isSmol(double d, double err) { return -err < d && d < err; }

namespace geom {

/// my vec3s are just [a, b, t] where ax+by<=t is the half-space represented
bool getCorner(const v::DVec<3> &e0, const v::DVec<3> &e1, v::DVec<2> &out) {
  double det = e0[0] * e1[1] - e0[1] * e1[0];
  if (isSmol(det, 1e-5)) {
    return false;
  }
  out[0] = e1[1] * e0[2] - e0[1] * e1[2];
  out[1] = e0[0] * e1[2] - e1[0] * e0[2];
  out /= det;
  return true;
}
/// for sorting
double getAngleThing(v::DVec<2> n) {
  double n0 = n[0] / v::norm1(n);
  if (n[1] >= 0) {
    return 1 - n0;
  }
  return 3 + n0;
}
/// 0 - keep all, 1 - infeasible, 2 - destroy b
/// precondition: a -> b -> c is clockwise, each turning <180deg
int checkStatus(const v::DVec<3> &a, const v::DVec<3> &b, const v::DVec<3> &c) {
  double detac = a[0] * c[1] - a[1] * c[0];
  // compare with getCorner's result dot product with {b[0], b[1]}
  double det3 =
      a[2] * (b[0] * c[1] - b[1] * c[0]) + c[2] * (a[0] * b[1] - a[1] * b[0]);
  det3 -= b[2] * detac;
  if (isSmol(det3, 1e-12) && isSmol(detac, 1e-10) &&
      a[0] * c[0] + a[1] * c[1] > 0 && a[0] * b[0] + a[1] * b[1] > 0) {
    double anorm = v::norm2(v::DVec<2>{a[0], a[1]});
    double bnorm = v::norm2(v::DVec<2>{b[0], b[1]});
    double cnorm = v::norm2(v::DVec<2>{c[0], c[1]});
    return 2 * (b[2] * anorm >= a[2] * bnorm && b[2] * cnorm >= c[2] * bnorm);
  }
  if (detac < 0) {
    return 2 * (det3 >= -1e-12);
  }
  return det3 > 0;
}
/// helper class
struct HalfSpace2D {
  v::DVec<2> n;
  double t;
  double tmp;
  HalfSpace2D *nextp, *prevp;
  HalfSpace2D() : n{}, t{}, tmp{} {}
  HalfSpace2D(const v::DVec<2> &n, double t) : n(n), t(t) {}
  HalfSpace2D(const v::DVec<3> &abt) : n{abt[0], abt[1]}, t(abt[2]) {}
  bool operator<(const HalfSpace2D &o) const { return tmp < o.tmp; }
  v::DVec<3> as3() const { return {n[0], n[1], t}; }
};
#define DO_PRINT_HALFS 0
/// for debugging, obviously
void printHalfs(HalfSpace2D *ptr) {
  (void)ptr;
#if DO_PRINT_HALFS
  std::cout << "halfs: " << std::endl;
  HalfSpace2D *p = ptr;
  int count = 0;
  do {
    std::cout << p->as3() << std::endl;
    p = p->nextp;
    count++;
    if (count > 1000) {
      std::cerr << "halfs is broken" << std::endl;
      throw std::exception();
    }
  } while (p != ptr);
  std::cout << std::endl;
#endif
}
/// ptr points to `b` as passed to checkStatus
template <HalfSpace2D *HalfSpace2D::*prevp = &HalfSpace2D::prevp>
HalfSpace2D *clearCheckStat3(HalfSpace2D *ptr, HalfSpace2D *begp) {
  while (ptr != begp) {
    int status = checkStatus(ptr->nextp->as3(), ptr->as3(), ptr->prevp->as3());
    switch (status) {
    case 0:
      return ptr;
    case 1:
      return nullptr;
    default: // case 2
      ptr->prevp->nextp = ptr->nextp;
      ptr->nextp->prevp = ptr->prevp;
      ptr = ptr->*prevp;
    }
  }
  return ptr;
}
/// 0 - success, 1 - infeasible, 2 - unbounded
/// out is a vector of indices in halfs (as it is after the function)
int evaluateFace(std::vector<HalfSpace2D> &halfs, std::vector<int> &out) {
  if (halfs.size() < 3) {
    return 2;
  }
  for (HalfSpace2D &half : halfs) {
    half.tmp = getAngleThing(half.n);
  }
  std::sort(halfs.begin(), halfs.end());
  for (std::size_t i = 1, sz = halfs.size(); i < sz; i++) {
    if (halfs[i].tmp - halfs[i - 1].tmp >= 2 - 1e-8) {
      return 2;
    }
  }
  if (halfs[0].tmp - halfs.back().tmp >= -2 - 1e-8) {
    return 2;
  }
  // i love circular doubly linked lists
  for (std::size_t ni = 0, pi = halfs.size(); pi--; ni = pi) {
    halfs[ni].prevp = &halfs[pi];
    halfs[pi].nextp = &halfs[ni];
  }
  HalfSpace2D *begp = &halfs[0];
  printHalfs(begp);
  HalfSpace2D *bptr = begp->nextp;
  while (bptr != begp) {
    if (begp == begp->nextp->nextp) {
      return 1;
    }
    bptr = clearCheckStat3(bptr, begp);
    if (!bptr) {
      return 1;
    }
    bptr = bptr->nextp;
  }
  HalfSpace2D *tmp = clearCheckStat3(begp, nullptr);
  if (!tmp) {
    return 1;
  }
  // now we know that it is feasible yay
  bool isClear = tmp == begp;
  bool isCtmp;
  HalfSpace2D *ptr;
  if (!isClear) {
    while (true) {
      ptr = tmp->nextp;
      begp = clearCheckStat3<&HalfSpace2D::nextp>(ptr, nullptr);
      if ((isCtmp = ptr == begp) && isClear) {
        break;
      }
      ptr = begp->prevp;
      tmp = clearCheckStat3(ptr, nullptr);
      if ((isClear = ptr == tmp) && isCtmp) {
        break;
      }
    }
  }
  out.clear();
  printHalfs(begp);
  ptr = begp;
  do {
    out.push_back(ptr - &halfs[0]);
    ptr = ptr->nextp;
  } while (ptr != begp);
  return 0;
}

} // namespace geom

/// represents a convex polytope in N dimensions
template <std::size_t N, class Attr> struct Polytope {
  std::vector<HalfSpace<N, Attr>> halfSpaces;

  Polytope() {}
  Polytope(std::initializer_list<HalfSpace<N, Attr>> lst) : halfSpaces(lst) {}

  bool contains(const v::DVec<N> &p) const {
    for (const HalfSpace<N, Attr> &h : halfSpaces) {
      if (!h.contains(p)) {
        return false;
      }
    }
    return true;
  }

private:
  struct HalfSpace3D {
    v::DVec<3> n;
    double t;
    Attr *attr;
  };
  std::vector<HalfSpace3D> tmp;
  std::vector<geom::HalfSpace2D> tmp2;
  std::vector<int> tmp3;
  std::vector<v::DVec<3>> tmp4;

public:
  /// calls fun(v::DVec<3>, v::DVec<3>, v::DVec<3>, Attr &) for each triangle
  /// triangles produced are ccw
  template <class Fun> void writeTriangles(const SliceDirs<N> &sd, Fun &&fun) {
    tmp.clear();
    tmp.reserve(halfSpaces.size());
    tmp2.reserve(halfSpaces.size());
    for (auto &half : halfSpaces) {
      v::DVec<3> n3 = {v::dot(sd.r, half.n), v::dot(sd.u, half.n),
                       v::dot(sd.f, half.n)};
      double t3 = half.t - v::dot(sd.c, half.n);
      double norm1 = v::norm1(n3);
      if (norm1 < 1e-8) {
        return;
      }
      n3 /= norm1;
      t3 /= norm1;
      tmp.push_back({n3, t3, &half.attr});
    }
    for (int i = 0, sz = tmp.size(); i < sz; i++) {
      v::DVec<3> n3 = tmp[i].n;
      double t3 = tmp[i].t;
      Attr &attr = *tmp[i].attr;
      v::DVec<3> a = {1, 0, 0};
      if (isSmol(n3[1], 1e-5) && isSmol(n3[2], 1e-5)) {
        a[1] = 1;
      }
      double nor2 = v::norm2(n3);
      a -= (v::dot(a, n3) / nor2) * n3;
      v::DVec<3> b = cross3(a, n3); // since we are in left-handed 3D
      v::DVec<3> c = (t3 / nor2) * n3;
      tmp2.clear();
      // FIXME: use vertex graph in 3D instead of following slow
      for (int j = 0; j < sz; j++) {
        if (j == i) {
          continue;
        }
        v::DVec<2> n2 = {v::dot(a, tmp[j].n), v::dot(b, tmp[j].n)};
        double t2 = tmp[j].t - v::dot(c, tmp[j].n);
        double norm2 = v::norm2(n2);
        if (isSmol(norm2, 1e-10)) {
          continue;
        }
        double norm = std::sqrt(norm2);
        geom::HalfSpace2D toadd;
        toadd.n = n2 / norm;
        toadd.t = t2 / norm;
        tmp2.push_back(toadd);
      }
      if (tmp2.size() < 3) {
        continue;
      }
      // infeasable or unbounded
      if (evaluateFace(tmp2, tmp3)) {
        continue;
      }
      tmp4.clear();
      for (std::size_t i = 0, j = tmp3.size() - 1, sz = tmp3.size(); i < sz;
           j = i++) {
        v::DVec<2> corner;
        if (!geom::getCorner(tmp2[tmp3[i]].as3(), tmp2[tmp3[j]].as3(),
                             corner)) {
          continue;
        }
        tmp4.push_back(corner[0] * a + corner[1] * b + c);
      }
      if (tmp4.size() < 3) {
        continue;
      }
      for (std::size_t i = 2, sz = tmp4.size(); i < sz; i++) {
        fun(tmp4[0], tmp4[i - 1], tmp4[i], attr);
      }
    }
  }
};

#endif // CONVEX_HPP_
