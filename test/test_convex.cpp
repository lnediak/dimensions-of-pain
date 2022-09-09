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

bool isUnbounded(const std::vector<geom::HalfSpace2D> &halfs) {
  std::size_t nh = halfs.size();
  std::vector<double> angles(nh);
  for (std::size_t i = nh; i--;) {
    angles[i] = std::atan2(halfs[i].n[1], halfs[i].n[0]);
  }
  std::sort(angles.begin(), angles.end());
  for (std::size_t i = nh; --i;) {
    if (angles[i] - angles[i - 1] >= 3.14) {
      return true;
    }
  }
  if (angles[0] - angles.back() >= -3.143) {
    return true;
  }
  return false;
}

/// 0 - succcess, 1 - infeasible, 2 - questionable feasibility
/// at least one pair (two inds in halfs) is required in each ele of out
int slowEvalFaceHalfs(
    const std::vector<geom::HalfSpace2D> &halfs,
    std::vector<std::unordered_set<v::IVec<2>, v::IVecHash<2>,
                                   v::EqualFunctor<v::IVec<2>, v::IVec<2>>>>
        &out) {
  std::size_t nh = halfs.size();
  struct Entry {
    int i, j;
    v::DVec<2> p;
  };
  std::vector<Entry> verts;
  std::vector<v::IVec<2>> ufh(nh);
  auto hGetP = [&ufh](int i) -> int & { return ufh[i][0]; };
  auto hGetSz = [&ufh](int i) -> int & { return ufh[i][1]; };
  for (int i = ufh.size(); i--;) {
    ufh[i] = {i, 1};
  }
  v::DVec<2> mid{0, 0};
  int numVerts = 0;
  for (int i = nh; i--;) {
    for (int j = i; j--;) {
      if (geom::isSmol(1e-2 * (halfs[i].n[0] * halfs[j].n[1] -
                               halfs[i].n[1] * halfs[j].n[0]))) {
        double vali = halfs[i].t / v::norm1(halfs[i].n);
        double valj = halfs[j].t / v::norm1(halfs[i].n);
        if (geom::isSmol(1e-2 * (vali - valj))) {
          uf::ufUnion(i, j, hGetP, hGetSz);
        }
      }
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
        mid += toadd;
        numVerts++;
        verts.push_back({i, j, toadd});
      }
    }
  }
  if (!verts.size()) {
    // infeasible
    return 1;
  }
  mid /= numVerts;
  for (std::size_t k = nh; k--;) {
    double df = v::dot(mid, halfs[k].n) - halfs[k].t;
    if (df > 1e-2) {
      return 1;
    } else if (df > -1e-2) {
      // questionable feasibility
      return 2;
    }
  }
  std::vector<v::IVec<2>> ufv(verts.size());
  auto vGetP = [&ufv](int i) -> int & { return ufv[i][0]; };
  auto vGetSz = [&ufv](int i) -> int & { return ufv[i][1]; };
  for (int i = ufv.size(); i--;) {
    ufv[i] = {i, 1};
  }
  for (std::size_t i = verts.size(); i--;) {
    for (std::size_t j = i; j--;) {
      if (geom::isSmol(1e-2 * v::norm1(verts[i].p - verts[j].p))) {
        uf::ufUnion(i, j, vGetP, vGetSz);
      }
    }
  }
  std::unordered_map<int, std::unordered_set<int>> hparts;
  for (std::size_t i = nh; i--;) {
    hparts[uf::ufFind(i, hGetP)].insert(i);
  }
  std::vector<std::unordered_set<v::IVec<2>, v::IVecHash<2>,
                                 v::EqualFunctor<v::IVec<2>, v::IVec<2>>>>
      tmp(verts.size());
  for (std::size_t i = verts.size(); i--;) {
    std::size_t vi = uf::ufFind(i, vGetP);
    std::size_t hi = uf::ufFind(verts[i].i, hGetP);
    std::size_t hj = uf::ufFind(verts[i].j, hGetP);
    if (hi == hj) {
      for (int i0 : hparts[hi]) {
        tmp[vi].insert({i0, i0});
      }
    } else {
      for (int i0 : hparts[hi]) {
        for (int i1 : hparts[hj]) {
          if (i0 >= i1) {
            tmp[vi].insert({i0, i1});
          } else {
            tmp[vi].insert({i1, i0});
          }
        }
      }
    }
  }
  out.clear();
  for (auto s : tmp) {
    if (s.size()) {
      out.push_back(s);
    }
  }
  return 0;
}

void testEvaluateFace(OneIterConf iterC) {
  std::cout << "entering testEvaluateFace" << std::endl;
  // std::uniform_int_distribution<int> disti(3, 20);
  std::uniform_int_distribution<int> disti(4, 4);
  std::uniform_real_distribution<double> distr(-100, 100);
  std::uniform_int_distribution<int> distc(0, 23);
  std::uniform_real_distribution<double> distsmol(1, 2);
  std::vector<geom::HalfSpace2D> halfs;
  std::vector<std::unordered_set<v::IVec<2>, v::IVecHash<2>,
                                 v::EqualFunctor<v::IVec<2>, v::IVec<2>>>>
      expected;
  std::vector<int> result;
  std::unordered_set<int> unoresult;
  std::vector<int> allowed;
  for (int spam = iterC.startIter; spam < iterC.endIter; spam++) {
    std::mt19937 mtrand(iterC.seed + spam);
    if (spam % iterC.iterPrintTick == 0) {
      std::cout << "Iteration #" << spam << std::endl;
    }
    int nh = disti(mtrand);
    halfs.clear();
    for (int i = nh; i--;) {
      v::DVec<3> randv;
      int ident = distc(mtrand);
      if (halfs.size() && ident <= 3) {
        randv = (2 * (ident % 2) - 1.) * distsmol(mtrand) *
                halfs[mtrand() % halfs.size()].as3();
        randv[1] += distr(mtrand) * 1e-6;
        if (ident < 2) {
          randv[2] = distr(mtrand);
        } else {
          randv[2] += distr(mtrand) * 1e-6;
        }
      } else if (ident <= 5) {
        // my beloved pathological half-planes
        randv = {1., distr(mtrand) * 1e-8, distr(mtrand)};
      } else {
        randv = {distr(mtrand), distr(mtrand), distr(mtrand)};
        while (randv[0] * randv[0] + randv[1] * randv[1] < 1e-2) {
          randv = {distr(mtrand), distr(mtrand), distr(mtrand)};
        }
      }
      halfs.emplace_back(randv);
    }
    if (isUnbounded(halfs)) {
      continue;
    }
    bool res = geom::evaluateFace(halfs, result); // false for infeasibility
    // this is after geom::evaluateFace cuz geom sorts halfs
    int status = slowEvalFaceHalfs(halfs, expected);
    // improvement: WRITE A TEST FOR THIS CASE
    if (status == 2) {
      continue;
    }
    if (res == status) {
      std::cerr << status << std::endl;
      std::cerr << "incorrect evaluation of feasibility" << std::endl;
      goto fail;
    }
    if (status) {
      // infeasible, rip
      continue;
    }
    for (std::size_t i0 = result.size(), i1 = 0; i0--; i1 = i0) {
      v::DVec<2> corner;
      if (!geom::getCorner(halfs[result[i0]].as3(), halfs[result[i1]].as3(),
                           corner)) {
        std::cerr << "geom::evaluateFace has produced an abomination"
                  << std::endl;
        goto fail;
      }
    }
    unoresult.clear();
    for (int i : result) {
      auto p = unoresult.insert(i);
      if (!p.second) {
        // this means there's a duplicate in result
        std::cerr << "geom::evaluateFace has a duplicate" << std::endl;
        goto fail;
      }
    }
    allowed.resize(nh);
    std::memset(&allowed[0], 0, sizeof(int) * nh);
    for (auto s : expected) {
      bool isGood = false;
      for (v::IVec<2> e : s) {
        allowed[e[0]] = allowed[e[1]] = 1;
        if (unoresult.count(e[0]) && unoresult.count(e[1])) {
          isGood = true;
        }
      }
      if (!isGood) {
        std::cerr << "geom::evaluateFace is missing something" << std::endl;
        goto fail;
      }
    }
    for (int i : result) {
      if (!allowed[i]) {
        std::cerr << "geom::evaluateFace has something extra" << std::endl;
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
  std::cout << std::setprecision(15);
  // testGetCorner({0, 10000000, 1000000, 1});

  // testEvaluateFace({0, 1, 1, 1});
  testEvaluateFace({52947, 10000000, 1, 1});
  // testEvaluateFace({0, 10000000, 100000, 1});
}
