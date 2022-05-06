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
      std::cout << "strange stats" << std::endl;
      return b[2] * da < a[2] * dt;
    }
    std::cout << "strange stats + 2" << std::endl;
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
    std::cout << "strange stats checkStatus 3" << std::endl;
    return 0;
  }
  double dt = b[0] * d[0] + b[1] * d[1];
  if (a[0] * c[1] - a[1] * c[0] < 0) {
    return 2 * (dt <= b[2]);
  }
  // because the face is bounded
  return dt > b[2];
}
struct HalfSpace2D {
  v::DVec<2> n;
  double t;
  double tmp;
  bool operator<(const HalfSpace2D &o) const { return tmp < o.tmp; }
  v::DVec<3> as3() const { return {n[0], n[1], t}; }
};
/// precondition: halfs describes a bounded face
/// returns true on success, false on infeasibility
bool evaluateFace(std::vector<HalfSpace2D> &halfs,
                  std::vector<v::DVec<3>> &out) {
  out.clear();
  out.reserve(halfs.size());
  for (HalfSpace2D &half : halfs) {
    half.tmp = getAngleThing(half.n);
  }
  std::sort(halfs.begin(), halfs.end());
  // halfs has at least 3 planes cuz face is bounded
  std::size_t ci = 0, bi = 1, sz = halfs.size();
  int status;
  while ((status = checkStatus(halfs[bi].as3(), halfs[ci].as3()))) {
    switch (status) {
    case 1:
      return false;
    case 3:
      ci = bi;
    }
    // i won't go out of range cuz face is bounded
    bi++;
  }
  std::list<std::size_t> outis = {ci, bi};
  for (std::size_t ai = bi + 1; ai < sz; ai++) {
    v::DVec<3> ha3 = halfs[ai].as3();
    v::DVec<3> hb3 = halfs[bi].as3();
    status = checkStatus(ha3, hb3);
    int statu;
    switch (status) {
    case 0:
      statu = checkStatus(ha3, hb3, halfs[ci].as3());
      switch (statu) {
      case 0:
        outis.push_back(ai);
        ci = bi;
        bi = ai;
        break;
      case 1:
        return false;
      default: { // case 2
        auto itb = --outis.end();
        auto itc = itb;
        --itc;
        outis.erase(itb);
        itb = itc;
        auto ite = outis.begin();
        while (itb != ite) {
          --itc;
          status = checkStatus(ha3, halfs[*itb].as3(), halfs[*itc].as3());
          if (!status) {
            break;
          }
          if (status == 1) {
            return false;
          }
          outis.erase(itb);
          itb = itc;
        }
        outis.push_back(ai);
      } // default
      } // switch
      break;
    case 1:
      return false;
    case 3:
      outis.back() = ai;
    }
  }
  // TODO: CONTINUE AS ABOVE, OR MODIFY?
  status = checkStatus(a, b);
  int statu;
  switch (status) {
  case 0:
    statu = checkStatus(a, b, c);
    switch (statu) {
    case 0:
      out.push_back(b);
      break;
    case 1:
      return false;
    }
    out.push_back(a);
    return true;
  case 1:
    return false;
  case 2:
    out.push_back(b);
    return true;
  }
  out.push_back(a);
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

