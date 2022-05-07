#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

#include "convex.hpp"

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
    while (geom::isSmol(e0[0] * e0[0] + e0[1] * e0[1])) {
      e0 = {distr(mtrand), distr(mtrand), distr(mtrand)};
    }
    v::DVec<3> e1 = {distr(mtrand), distr(mtrand), distr(mtrand)};
    while (geom::isSmol(e1[0] * e1[0] + e1[1] * e1[1])) {
      e1 = {distr(mtrand), distr(mtrand), distr(mtrand)};
    }
    v::DVec<2> out;
    if (!geom::getCorner(e0, e1, out)) {
      if (!geom::isSmol(e0[0] * e1[1] - e0[1] * e1[0])) {
        std::cerr << "determinant not small???" << std::endl;
        goto fail;
      }
    } else {
      if (!geom::isSmol((out[0] * e0[0] + out[1] * e0[1] - e0[2]) / 100)) {
        std::cerr << "result not lying on edge of e0" << std::endl;
        goto fail;
      }
      if (!geom::isSmol((out[0] * e1[0] + out[1] * e1[1] - e1[2]) / 100)) {
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

/// 0 - succcess, 1 - infeasible, 2 - unbounded, 3 - questionable feasibility
/// out[k][2] are all irrelevant btw
int slowEvalFacePoints(const std::vector<geom::HalfSpace2D> &halfs,
                       std::vector<v::DVec<3>> &out) {
  out.clear();
  std::size_t nh = halfs.size();
  std::cout << "nh: " << nh << std::endl;
  v::DVec<2> mid{0, 0};
  for (std::size_t i = nh; i--;) {
    for (std::size_t j = i; j--;) {
      v::DVec<2> toadd;
      if (!geom::getCorner(halfs[j].as3(), halfs[i].as3(), toadd)) {
        continue;
      }
      bool shouldAdd = true;
      for (std::size_t k = nh; k--;) {
        if (v::dot(toadd, halfs[k].n) > halfs[k].t + 1e-6) {
          shouldAdd = false;
          break;
        }
      }
      if (shouldAdd) {
        out.push_back({toadd[0], toadd[1], 0});
        mid += toadd;
      }
    }
  }
  if (!out.size()) {
    std::cout << "empty out?" << std::endl;
    for (std::size_t i = nh; i--;) {
      v::DVec<2> toadd;
      if (geom::isSmol(halfs[i].n[0])) {
        toadd = {0, halfs[i].t / halfs[i].n[1]};
      } else {
        toadd = {halfs[i].t / halfs[i].n[0], 0};
      }
      bool shouldAdd = true;
      for (std::size_t k = nh; k--;) {
        if (v::dot(toadd, halfs[k].n) > halfs[k].t + 1e-6) {
          shouldAdd = false;
          break;
        }
      }
      if (shouldAdd) {
        return 2;
      }
    }
    std::cout << "returned 1" << std::endl;
    return 1;
  }
  mid /= out.size();
  for (std::size_t k = nh; k--;) {
    double df = v::dot(mid, halfs[k].n) - halfs[k].t;
    if (df > 1e-6) {
      std::cout << "mid nonfunctional?" << std::endl;
      return 1;
    } else if (df > -1e-5) {
      return 3;
    }
  }
  for (v::DVec<3> &e : out) {
    e[2] = std::atan2(e[1] - mid[1], e[0] - mid[0]);
  }
  std::sort(out.begin(), out.end(),
            [](const v::DVec<3> &a, const v::DVec<3> &b) -> bool {
              return a[2] < b[2];
            });
  for (std::size_t i0 = out.size() - 1, i1 = 0, sz = out.size(); i1 < sz;
       i0 = i1++) {
    v::DVec<2> dir = {out[i1][1] - out[i0][1], out[i0][0] - out[i1][0]};
    bool isUnb = true;
    for (std::size_t k = nh; k--;) {
      if (v::dot(dir, halfs[k].n) > 1e6) {
        isUnb = false;
        break;
      }
    }
    if (isUnb) {
      return 2;
    }
  }
  return 0;
}

void testEvaluateFace(OneIterConf iterC) {
  std::cout << "entering testEvaluateFace" << std::endl;
  std::uniform_int_distribution<int> disti(3, 100);
  std::uniform_real_distribution<double> distr(-100, 100);
  std::vector<geom::HalfSpace2D> halfs;
  std::vector<v::DVec<3>> out0, out1;
  for (int spam = iterC.startIter; spam < iterC.endIter; spam++) {
    std::mt19937 mtrand(iterC.seed + spam);
    if (spam % iterC.iterPrintTick == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    int nh = disti(mtrand);
    halfs.clear();
    for (int i = nh; i--;) {
      v::DVec<3> randv = {distr(mtrand), distr(mtrand), distr(mtrand)};
      while (geom::isSmol(randv[0] * randv[0] + randv[1] * randv[1])) {
        randv = {distr(mtrand), distr(mtrand), distr(mtrand)};
      }
      halfs.push_back({{randv[0], randv[1]}, randv[2]});
    }
    int res = slowEvalFacePoints(halfs, out0);
    if (res >= 2) {
      continue;
    }
    std::size_t out0i = 0;
    if (geom::evaluateFace(halfs, out1) == res) {
      std::cerr << res << std::endl;
      std::cerr << "incorrect evaluation of feasibility" << std::endl;
      goto fail;
    }
    for (std::size_t i0 = out1.size() - 1, i1 = 0, sz = out1.size(); i1 < sz;
         i0 = i1++) {
      v::DVec<2> corner;
      if (!geom::getCorner(out1[i0], out1[i1], corner)) {
        std::cerr << "geom::evaluateFace has produced an abomination"
                  << std::endl;
        goto fail;
      }
      bool isGood = false;
      for (std::size_t i2 = out0i, sz0 = out0.size(); i2 < sz0; i2++) {
        if (v::norm2(corner - v::DVec<2>{out0[i2][0], out0[i2][1]}) < 1e-8) {
          out0i = i2;
          isGood = true;
          break;
        }
      }
      if (isGood) {
        continue;
      }
      for (std::size_t i2 = 0; i2 < out0i; i2++) {
        if (v::norm2(corner - v::DVec<2>{out0[i2][0], out0[i2][1]}) < 1e-8) {
          out0i = i2;
          isGood = true;
          break;
        }
      }
      if (!isGood) {
        std::cerr << "a corner from evaluateFace is funny" << std::endl;
        goto fail;
      }
    }
    continue;
  fail:
    std::cerr << "Failed on iteration #" << spam << std::endl;
    throw std::exception();
  }
}

int main() {
  // testGetCorner({0, 10000000, 1000000, 1});
  testEvaluateFace({0, 1000000, 100000, 1});
}

