#ifndef UTIL_HPP_
#define UTIL_HPP_

#include "vector.hpp"

/// this is for easily passing into a gl shader
union Color {
  unsigned char c[4];
  float f;
};

template <std::size_t N> struct SliceDirs {
  v::DVec<N> c;       /// camera
  v::DVec<N> r, u, f; /// right, up, forward
  double rm, um;      /// r and u multipliers (at f dist. 1) (aspect ratio+fov)
  double fm;          /// forward multiplier (render distance)
};

v::DVec<3> cross3(const v::DVec<3> &a, const v::DVec<3> &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

// ---------------------- ACTUAL UTILITY FUNCTIONS BELOW -----------------------

// Source: https://www.lomont.org/papers/2003/InvSqrt.pdf
double fastInvSqrt(double a) {
  double a2 = 0.5 * a;
  long i = *(long *)&a;
  i = 0x5fe6ec85e7de30da - (i >> 1);
  double b = *(double *)&i;
  // don't remove any of these. I need the precision for the quaternions
  b = b * (1.5 - a2 * b * b);
  b = b * (1.5 - a2 * b * b);
  b = b * (1.5 - a2 * b * b);
  return b;
}

#define GETFLOOR_ENDIANESS_INDEX 0

/// do not use values out of the range of a short
int fastFloor(float val) {
  val -= 0.5;
  val += 8388608 * 1.5;
  return ((short *)&val)[GETFLOOR_ENDIANESS_INDEX];
}

int fastFloorD(double val) {
  val -= 0.5;
  val += 6755399441055744;
  return ((int *)&val)[GETFLOOR_ENDIANESS_INDEX];
}

/// This is basically floor(val - small_number)
int harshFloor(float val) {
  val -= 0.50001;
  val += 8388608 * 1.5;
  return ((short *)&val)[GETFLOOR_ENDIANESS_INDEX];
}

/// do not use values out of the range of a short
int fastRound(float val) {
  val += 8388608 * 1.5;
  return ((short *)&val)[GETFLOOR_ENDIANESS_INDEX];
}

#endif // UTIL_HPP_

