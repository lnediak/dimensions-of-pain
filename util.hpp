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

#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wuninitialized"

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

#pragma GCC diagnostic pop

/*
// ---------- NOW FOR SOME STUFF TO SUPPLEMENT VVEC ----------

template <class T> struct FloorU {
  static int apply(T a) {
    int i = a;
    return i - (a < i);
  }
};
struct FastFloorU {
  static int apply(float a) { return fastFloor(a); }
};
struct FastFloorDU {
  static int apply(double a) { return fastFloorD(a); }
};
struct FastRoundU {
  static int apply(float a) { return fastRound(a); }
};
struct Int2FloatU {
  static float apply(int a) { return a; }
};
struct Int2DoubleU {
  static double apply(int a) { return a; }
};

template <class T> struct AbsU {
  static T apply(T a) { return a > 0 ? a : -a; }
};
template <class T> struct L1NormU {
  static T apply(T a, T b) { return a + AbsU<T>::apply(b); }
};

template <class T> struct EI {
  T a;
  std::size_t b;

  EI<T> get(const EI<T> &c) const {
    if (c.a < a) {
      return c;
    }
    return *this;
  }
};
template <class T, std::size_t N, class A> struct MinI {

  typedef MinI<T, N - 1, A> smaller;
  const A &a;
  EI<T> evaluate() const {
    return smaller{a}.evaluate().get({a[N - 1], N - 1});
  }
};
template <class T, class A> struct MinI<T, 1, A> {
  const A &a;
  EI<T> evaluate() const { return {a[0], 0}; }
};

template <class T, std::size_t N, class A>
using Floor = v::UnaryOp<int, N, FloorU<T>, A>;
template <std::size_t N, class A>
using FastFloor = v::UnaryOp<int, N, FastFloorU, A>;
template <std::size_t N, class A>
using FastFloorD = v::UnaryOp<int, N, FastFloorDU, A>;
template <std::size_t N, class A>
using FastRound = v::UnaryOp<int, N, FastRoundU, A>;
template <std::size_t N, class A>
using Int2Float = v::UnaryOp<float, N, Int2FloatU, A>;
template <std::size_t N, class A>
using Int2Double = v::UnaryOp<double, N, Int2DoubleU, A>;

template <class A, class = typename rem_cvr<A>::thisisavvec>
Floor<typename A::value_type, A::size, A> vfloor(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, float>::value>::type>
FastFloor<A::size, A> vfastFloor(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, double>::value>::type>
FastFloorD<A::size, A> vfastFloorD(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, float>::value>::type>
FastRound<A::size, A> vfastRound(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, int>::value>::type>
Int2Float<A::size, A> vint2float(const A &a) {
  return {a};
}
template <class A, class = typename rem_cvr<A>::thisisavvec,
          class = typename std::enable_if<
              std::is_same<typename A::value_type, int>::value>::type>
Int2Double<A::size, A> vint2double(const A &a) {
  return {a};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
EI<typename A::value_type> vminI(const A &a) {
  return MinI<typename A::value_type, A::size, A>{a}.evaluate();
}
*/

#endif // UTIL_HPP_
