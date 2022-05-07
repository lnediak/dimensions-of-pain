#ifndef CONVEX_HPP_
#define CONVEX_HPP_

#include <algorithm>
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

namespace geom {

bool isSmol(double d) { return -1e-8 < d && d < 1e-8; }

/// my vec3s are just [a, b, t] where ax+by<=t is the half-space represented
bool getCorner(const v::DVec<3> &e0, const v::DVec<3> &e1, v::DVec<2> &out) {
  double det = e0[0] * e1[1] - e0[1] * e1[0];
  if (isSmol(det)) {
    return false;
  }
  out[0] = e1[1] * e0[2] - e0[1] * e1[2];
  out[1] = e0[0] * e1[2] - e1[0] * e0[2];
  out /= det;
  return true;
}
/// for sorting
double getAngleThing(v::DVec<2> n) {
  n /= v::norm1(n);
  if (n[1] >= 0) {
    return 1 - n[0];
  }
  return 3 + n[0];
}
/// 0 - normal, 1 - infeasible, 2 - destroy a, 3 - destroy b
int checkStatus(const v::DVec<3> &a, const v::DVec<3> &b) {
  if (isSmol(a[0] * b[1] - a[1] * b[0])) {
    double dt = a[0] * b[0] + a[1] * b[1];
    double da = a[0] * a[0] + a[1] * a[1];
    if (dt < 0) {
      return b[2] * da < a[2] * dt;
    }
    return 2 + (b[2] * da >= a[2] * dt);
  }
  return 0;
}
/// 0 - keep all, 1 - infeasible, 2 - destroy b
/// assumes !checkStatus(a, b) and !checkStatus(b, c)
int checkStatus(const v::DVec<3> &a, const v::DVec<3> &b, const v::DVec<3> &c) {
  int acstat = checkStatus(a, c);
  if (acstat) {
    // because the face is bounded
    return 2;
  }
  v::DVec<2> d;
  if (!getCorner(a, c, d)) {
    // e.g. a square
    return 0;
  }
  double dt = b[0] * d[0] + b[1] * d[1];
  if (a[0] * c[1] - a[1] * c[0] < 0) {
    return 2 * (dt <= b[2]);
  }
  // because the face is bounded
  return dt > b[2];
}
/// helper class
struct HalfSpace2D {
  v::DVec<2> n;
  double t;
  double tmp;
  HalfSpace2D *nextp, *prevp;
  HalfSpace2D(const v::DVec<2> &n, double t) : n(n), t(t) {}
  HalfSpace2D(const v::DVec<3> &abt) : n{abt[0], abt[1]}, t(abt[2]) {}
  bool operator<(const HalfSpace2D &o) const { return tmp < o.tmp; }
  v::DVec<3> as3() const { return {n[0], n[1], t}; }
};
/// assumes ptr has nextp and prevp set up, links to bounded face
template <HalfSpace2D *HalfSpace2D::*nextp = &HalfSpace2D::nextp,
          HalfSpace2D *HalfSpace2D::*prevp = &HalfSpace2D::prevp>
HalfSpace2D *clearCheckStat2(HalfSpace2D *ptr) {
  while (true) {
    int status = checkStatus((ptr->*nextp)->as3(), ptr->as3());
    switch (status) {
    case 0:
      return ptr;
    case 1:
      return nullptr;
    case 2:
      ptr->*nextp = ptr->*nextp->*nextp;
      ptr->*nextp->*prevp = ptr;
      break;
    default: // case 3
      ptr->*prevp->*nextp = ptr->*nextp;
      ptr->*nextp->*prevp = ptr->*prevp;
      ptr = ptr->*nextp;
    }
  }
  return ptr;
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
/// for debugging, obviously
void printHalfs(HalfSpace2D *ptr) {
  std::cout << "halfs: " << std::endl;
  HalfSpace2D *p = ptr;
  do {
    std::cout << p->as3() << std::endl;
    p = p->nextp;
  } while (p != ptr);
  std::cout << std::endl;
}
/// precondition: halfs describes a bounded face or is infeasible
/// returns true on success, false on infeasibility
bool evaluateFace(std::vector<HalfSpace2D> &halfs,
                  std::vector<v::DVec<3>> &out) {
  for (HalfSpace2D &half : halfs) {
    half.tmp = getAngleThing(half.n);
  }
  std::sort(halfs.begin(), halfs.end());
  // i love circular doubly linked lists
  for (std::size_t ni = 0, pi = halfs.size(); pi--; ni = pi) {
    halfs[ni].prevp = &halfs[pi];
    halfs[pi].nextp = &halfs[ni];
  }
  HalfSpace2D *begp = &halfs[0];
  printHalfs(begp);
  begp = clearCheckStat2(begp);
  if (!begp) {
    return false;
  }
  // this is only necessary for a pathological case
  begp = clearCheckStat2<&HalfSpace2D::prevp, &HalfSpace2D::nextp>(begp);
  if (!begp) {
    return false;
  }
  printHalfs(begp);
  HalfSpace2D *bptr = begp->nextp;
  while (bptr != begp) {
    printHalfs(begp);
    bptr = clearCheckStat2(bptr);
    if (!bptr) {
      return false;
    }
    bptr = bptr->nextp;
  }
  printHalfs(begp);
  std::cout << "entering clearCheckStat3 loop" << std::endl;
  bptr = bptr->nextp;
  while (bptr != begp) {
    printHalfs(begp);
    bptr = clearCheckStat3(bptr, begp);
    if (!bptr) {
      return false;
    }
    bptr = bptr->nextp;
  }
  std::cout << "left loop" << std::endl;
  printHalfs(begp);
  begp = clearCheckStat3(begp, nullptr);
  if (!begp) {
    return false;
  }
  begp = clearCheckStat3<&HalfSpace2D::nextp>(begp->nextp, nullptr);
  if (!begp) {
    return false;
  }
  printHalfs(begp);
  out.clear();
  HalfSpace2D *ptr = begp;
  do {
    out.push_back(ptr->as3());
    ptr = ptr->nextp;
  } while (ptr != begp);
  return true;
}

} // namespace geom

/// represents a convex polytope in N dimensions
template <std::size_t N, class Attr> struct Polytope {
  std::vector<HalfSpace<N, Attr>> halfSpaces;

  bool contains(const v::DVec<N> &p) const {
    for (const HalfSpace<N, Attr> &h : halfSpaces) {
      if (!h.contains(p)) {
        return false;
      }
    }
    return true;
  }

  std::vector<v::DVec<4>> tmp;
  template <class Fun> void writeTriangles(const SliceDirs<N> &sd, Fun &&fun) {
    tmp.clear();
    tmp.reserve(halfSpaces.size());
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
      if (geom::isSmol(n3[1]) && geom::isSmol(n3[2])) {
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

#endif // CONVEX_HPP_

