#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "convex.hpp"
#include "union_find.hpp"

#define DEBUG_PRINT(...) std::cout << __VA_ARGS__ << std::endl

struct OneIterConf {
  int startIter, endIter, iterPrintTick, seed;
};

void testGetCorner(OneIterConf iterC) {
  std::cout << "entering testGetCorner" << std::endl;
  std::uniform_real_distribution<double> distr(-100, 100);
  for (int spam = iterC.startIter; spam < iterC.endIter; spam++) {
    std::mt19937 mtrand(iterC.seed + spam);
    if (spam % iterC.iterPrintTick == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    v::DVec<3> e0 = {distr(mtrand), distr(mtrand), distr(mtrand)};
    while (isSmol(e0[0] * e0[0] + e0[1] * e0[1], 1e-5)) {
      e0 = {distr(mtrand), distr(mtrand), distr(mtrand)};
    }
    v::DVec<3> e1 = {distr(mtrand), distr(mtrand), distr(mtrand)};
    while (isSmol(e1[0] * e1[0] + e1[1] * e1[1], 1e-5)) {
      e1 = {distr(mtrand), distr(mtrand), distr(mtrand)};
    }
    v::DVec<2> out;
    if (!geom::getCorner(e0, e1, out)) {
      if (!isSmol(e0[0] * e1[1] - e0[1] * e1[0], 1e-5)) {
        std::cerr << "determinant not small???" << std::endl;
        goto fail;
      }
    } else {
      if (!isSmol(out[0] * e0[0] + out[1] * e0[1] - e0[2], 1e-3)) {
        std::cerr << "result not lying on edge of e0" << std::endl;
        goto fail;
      }
      if (!isSmol(out[0] * e1[0] + out[1] * e1[1] - e1[2], 1e-3)) {
        std::cerr << "result not lying on edge of e1" << std::endl;
        goto fail;
      }
    }
    continue;
  fail:
    std::cerr << "Failed on iteration #" << spam << std::endl;
    throw std::exception();
  }
}

/// 0 - bounded, 1 - unbounded, 2 - borderline
int getBoundedness(const std::vector<geom::HalfSpace2D> &halfs,
                   double low = 3.1, double high = 3.2) {
  std::size_t nh = halfs.size();
  std::vector<double> angles(nh);
  for (std::size_t i = nh; i--;) {
    angles[i] = std::atan2(halfs[i].n[1], halfs[i].n[0]);
  }
  std::sort(angles.begin(), angles.end());
  for (std::size_t i = nh; --i;) {
    if (angles[i] - angles[i - 1] >= low) {
      if (angles[i] - angles[i - 1] >= high) {
        return 1;
      }
      return 2;
    }
  }
  if (angles[0] - angles.back() >= -high) {
    if (angles[0] - angles.back() >= -low) {
      return 1;
    }
    return 2;
  }
  return 0;
}

// absolutely cursed, but I don't know of an easier way for this to be stable
double binsearchBorder(const v::DVec<2> &base, const v::DVec<2> &dir,
                       const v::DVec<2> &n, double t, bool &isIncreasing) {
  double low = -1e9;
  double high = 1e9;
  double det = v::dot(dir, n);
  isIncreasing = det >= 0;
  while (high - low > 1e-8) {
    double mid = (low + high) / 2;
    double crit = v::dot(base + mid * dir, n);
    if ((crit <= t && isIncreasing) || (crit > t && !isIncreasing)) {
      low = mid;
    } else if ((crit <= t && !isIncreasing) || (crit > t && isIncreasing)) {
      high = mid;
    }
  }
  return isIncreasing ? low : high;
}
/*
  indices are the indices in halfs we consider, index is of working edge
  0 - found proper interior point, 1 - found borderline point, 2 - no point
*/
int getRelevantPoint(const std::vector<geom::HalfSpace2D> &halfs,
                     const std::vector<int> &indices, int index,
                     v::DVec<2> &out) {
  v::DVec<2> n = halfs[index].n;
  double t = halfs[index].t;
  v::DVec<2> base;
  if (isSmol(n[0], 1e-2)) {
    base = {0, t / n[1]};
  } else {
    base = {t / n[0], 0};
  }
  v::DVec<2> dir = {-n[1], n[0]};
  double lowmax = -1e9;
  double lowmin = -1e9;
  double highmax = 1e9;
  double highmin = 1e9;
  for (int i : indices) {
    if (i == index) {
      continue;
    }
    bool isIncreasing, isIncr2;
    double t1 =
        binsearchBorder(base, dir, halfs[i].n, halfs[i].t - 1e-5, isIncreasing);
    double t2 =
        binsearchBorder(base, dir, -halfs[i].n, -halfs[i].t - 1e-5, isIncr2);
    if (isIncreasing) {
      if (t1 < highmin) {
        highmin = t1;
      }
    } else {
      if (t1 > lowmax) {
        lowmax = t1;
      }
    }
    if (isIncr2) {
      if (t2 > lowmin) {
        lowmin = t2;
      }
    } else {
      if (t2 < highmax) {
        highmax = t2;
      }
    }
  }
  if (lowmin > highmax) {
    return 2;
  }
  if (lowmax > highmin) {
    out = base + lowmin * dir;
    return 1;
  }
  out = base + lowmax * dir;
  return 0;
}
struct RandomStuffs {
  std::mt19937 mt;
  std::uniform_real_distribution<double> unif;
  RandomStuffs(int seed) : mt(seed), unif(0, 1) {}
};
v::DVec<2> randVec(RandomStuffs &rands) {
  v::DVec<2> n;
  do {
    n[0] = rands.unif(rands.mt) * 14 - 7;
    n[1] = rands.unif(rands.mt) * 14 - 7;
  } while (v::norm2(n) < 1e-2);
  return n;
}
double randSmol(RandomStuffs &rands, int n = 2) {
  if (rands.mt() % n) {
    return 0;
  }
  double tmp = 1 + 1e-10 - rands.unif(rands.mt);
  return std::exp(-1 / (tmp * tmp)) * (2 * (rands.mt() % 2) - 1);
}
/// replaces out with face containing the origin
void generateFeasibleFace(RandomStuffs &rands, int sz,
                          std::vector<geom::HalfSpace2D> &out) {
  out.clear();
  for (int i = 0; i < sz; i++) {
    v::DVec<2> n;
    if (i >= 2 && !(rands.mt() % 3)) {
      int j = rands.mt() % i;
      int k = rands.mt() % (i - 1);
      k += (k >= j);
      do {
        double w = rands.unif(rands.mt) * 5 - 3;
        if (w < 0) {
          w = w > -1 ? 0 : w + 1;
        }
        if (w > 1) {
          w = w < 2 ? 1 : w - 1;
        }
        n = out[j].n + w * (out[k].n - out[j].n);
        n[0] += randSmol(rands);
        n[1] += randSmol(rands);
      } while (v::norm2(n) < 1e-2 || v::norm2(n) >= 100);
    } else {
      n = randVec(rands);
    }
    out.emplace_back(n, 1);
  }
  for (int i = 0; i < sz; i++) {
    double norm = v::norm2(out[i].n);
    out[i].n /= norm;
    out[i].t /= norm;
  }
}
/// adds 3 edges to out that are infeasible
void addInfeasible(RandomStuffs &rands, std::vector<geom::HalfSpace2D> &out) {
  v::DVec<2> a[3];
  while (true) {
    a[0] = randVec(rands);
    a[0] /= v::norm2(a[0]);
    a[1] = randVec(rands);
    a[1] /= v::norm2(a[1]);
    a[2] = randVec(rands);
    a[2] /= v::norm2(a[2]);
    double x = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    double y = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    double z = a[2][0] * a[0][1] - a[2][1] * a[0][0];
    if ((x > 5e-2 && y > 5e-2 && z > 5e-2) ||
        (x < -5e-2 && y < -5e-2 && z < -5e-2)) {
      break;
    }
  }
  out.emplace_back(a[0], -rands.unif(rands.mt) * 9.99 - 1e-2);
  out.emplace_back(a[1], -rands.unif(rands.mt) * 9.99 - 1e-2);
  out.emplace_back(a[2], -rands.unif(rands.mt) * 9.99 - 1e-2);
}
/// adds n nearly parallel edges passing very close to the origin
void addDegenerate(RandomStuffs &rands, int n,
                   std::vector<geom::HalfSpace2D> &out) {
  v::DVec<2> vec = randVec(rands);
  vec /= v::norm2(vec);
  while (n--) {
    out.emplace_back(vec, randSmol(rands, 6));
    vec[0] += randSmol(rands, 3);
    vec[1] += randSmol(rands, 3);
  }
}
v::DVec<2> offsetFace(RandomStuffs &rands,
                      std::vector<geom::HalfSpace2D> &face) {
  double x = rands.unif(rands.mt) * 20 - 10;
  double y = rands.unif(rands.mt) * 20 - 10;
  v::DVec<2> offset = {x, y};
  for (std::size_t i = face.size(); i--;) {
    face[i].t += v::dot(face[i].n, offset);
  }
  return offset;
}

int subBoundedness(const std::vector<geom::HalfSpace2D> &halfs,
                   const std::vector<int> &indices, double low, double high) {
  std::vector<geom::HalfSpace2D> sub(indices.size());
  for (std::size_t i = indices.size(); i--;) {
    sub[i] = halfs[indices[i]];
  }
  return getBoundedness(sub, low, high);
}
bool isValidSolution(const std::vector<geom::HalfSpace2D> &halfs,
                     const std::vector<int> &indices) {
  std::vector<int> allI(halfs.size());
  for (std::size_t i = halfs.size(); i--;) {
    allI[i] = i;
  }
  std::vector<int> results(halfs.size());
  for (std::size_t i = halfs.size(); i--;) {
    v::DVec<2> tmp;
    results[i] = getRelevantPoint(halfs, allI, i, tmp);
  }
  for (int i : indices) {
    if (results[i] < 0) {
      // is a duplicated index
      return false;
    }
    results[i] = ~results[i];
  }
  for (std::size_t i = halfs.size(); i--;) {
    int r = results[i];
    if (!r || r == ~2) {
      return false;
    }
  }
  for (int i : indices) {
    v::DVec<2> point;
    if (!getRelevantPoint(halfs, indices, i, point)) {
      for (const geom::HalfSpace2D &edge : halfs) {
        if (v::dot(point, edge.n) > edge.t + 1e-5) {
          return false;
        }
      }
    }
  }
  return true;
}

void printFail(int iteri) {
  std::cerr << "Failed on iteration #" << iteri << std::endl;
  throw std::exception();
}
void testEvaluateFace(OneIterConf iterC) {
  std::cout << "entering testEvaluateFace" << std::endl;
  std::vector<geom::HalfSpace2D> halfs;
  std::vector<int> indices;
  for (int iteri = iterC.startIter; iteri < iterC.endIter; iteri++) {
    RandomStuffs rands(iterC.seed + iteri);
    if (iteri % iterC.iterPrintTick == 0) {
      std::cout << "Iteration #" << iteri << std::endl;
    }
    int nh = rands.mt() % 20;
    generateFeasibleFace(rands, nh, halfs);
    int boundedness = getBoundedness(halfs);
    int res = geom::evaluateFace(halfs, indices);
    if (!res) {
      if (subBoundedness(halfs, indices, 3.14159, 3.141595)) {
        printFail(iteri);
      }
    }
    if (boundedness == 1 && res != 2) {
      printFail(iteri);
    }
    if (!boundedness) {
      if (res) {
        printFail(iteri);
      }
      if (!isValidSolution(halfs, indices)) {
        printFail(iteri);
      }
    }

    addDegenerate(rands, 3, halfs);
    boundedness = getBoundedness(halfs);
    res = geom::evaluateFace(halfs, indices);
    if (!res) {
      if (subBoundedness(halfs, indices, 3.141595, 3.14159)) {
        printFail(iteri);
      }
    }
    if (boundedness == 1 && res != 2) {
      printFail(iteri);
    }
    if (!boundedness && res == 2) {
      printFail(iteri);
    }

    addInfeasible(rands, halfs);
    boundedness = getBoundedness(halfs);
    res = geom::evaluateFace(halfs, indices);
    if (!res) {
      printFail(iteri);
    }
    if (!boundedness && res == 2) {
      printFail(iteri);
    }
  }
}

int main() {
  std::cout << std::setprecision(15);
  // testGetCorner({0, 10000000, 1000000, 1});

  testEvaluateFace({0, 10000000, 100000, 1});
}
