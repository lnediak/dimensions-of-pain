#include <iostream>
#include <random>

#include "convex.hpp"

struct OneIterConf {
  int startIter, endIter, iterPrintTick;
};

void testGetCorner(OneIterConf iterC, int seed) {
  std::uniform_real_distribution<double> distr(-100, 100);
  for (int spam = iterC.startIter; spam < iterC.endIter; spam++) {
    std::mt19937 mtrand(seed + spam);
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
        std::cout << "determinant not small???" << std::endl;
        goto fail;
      }
    } else {
      if (!geom::isSmol((out[0] * e0[0] + out[1] * e0[1] - e0[2]) / 100)) {
        std::cout << "result not lying on edge of e0" << std::endl;
        goto fail;
      }
      if (!geom::isSmol((out[0] * e1[0] + out[1] * e1[1] - e1[2]) / 100)) {
        std::cout << "result not lying on edge of e1" << std::endl;
        goto fail;
      }
    }
    continue;
  fail:
    std::cerr << "Failed on iteration #" << spam << std::endl;
    throw std::exception();
  }
}

bool ridiculousEvalFace(const std::vector<HalfSpace2D> &halfs, std::vector<v::DVec<3>> &out) {
  out.clear();
  std::size_t nh = half.size();
  for (std::size_t i = nh; i--;) {
    for (std::size_t j = i; j--;) {
      // TODO: CONTINUE
    }
  }
}

int main() { testGetCorner({0, 10000000, 1000000}, 1); }

