#ifndef V_VECTOR_HPP_
#define V_VECTOR_HPP_

#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <type_traits>

template <class T>
using rem_cvr =
    typename std::remove_cv<typename std::remove_reference<T>::type>::type;

namespace v {

// ------------------------BASIC OPERATIONS-------------

template <class T, std::size_t N, class U, class A, class B> struct BinaryOp {

  const A &a;
  const B &b;

  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  BinaryOp(const A &a, const B &b) : a(a), b(b) {}

  T operator[](std::size_t i) const { return U::apply(a[i], b[i]); }
};

template <class T, std::size_t N, class U, class A> struct ScalarOp {

  const A &a;
  T s;

  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  ScalarOp(const A &a, T s) : a(a), s(s) {}

  T operator[](std::size_t i) const { return U::apply(a[i], s); }
};

template <class T, std::size_t N, class U, class A> struct ReverseScalarOp {

  const A &a;
  T s;

  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  ReverseScalarOp(const A &a, T s) : a(a), s(s) {}

  T operator[](std::size_t i) const { return U::apply(s, a[i]); }
};

template <class T, std::size_t N, class U1, class U2, class A> struct Combine {

  const A &a;

  Combine(const A &a) : a(a) {}

  T evaluate() const {
    T toret = U2::apply(a[0]);
    for (std::size_t i = 1; i < N; i++) {
      toret = U1::apply(toret, a[i]);
    }
    return toret;
  }
};

template <std::size_t N, class A, class B> struct IsEqual {
  bool operator()(const A &a, const B &b) const noexcept {
    for (std::size_t i = N; i--;) {
      if (a[i] != b[i]) {
        return false;
      }
    }
    return true;
  }
};

template <class A, class B> struct EqualFunctor {
  bool operator()(const A &a, const B &b) const noexcept {
    return v::IsEqual<rem_cvr<A>::size, A, B>{}(a, b);
  }
};

template <class T, std::size_t N, class U, class A> struct UnaryOp {

  const A &a;

  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  UnaryOp(const A &a) : a(a) {}

  T operator[](std::size_t i) const { return U::apply(a[i]); }
};

// ------------------------TYPEDEFS------------------

template <class T> struct IdU {
  static T apply(T a) { return a; }
};

template <class T> struct AddU {
  static T apply(T a, T b) { return a + b; }
};
template <class T, std::size_t N, class A, class B>
using Add = BinaryOp<T, N, AddU<T>, A, B>;
template <class T, std::size_t N, class A>
using SAdd = ScalarOp<T, N, AddU<T>, A>;
template <class T, std::size_t N, class A>
using Sum = Combine<T, N, AddU<T>, IdU<T>, A>;

template <class T> struct SubU {
  static T apply(T a, T b) { return a - b; }
};
template <class T, std::size_t N, class A, class B>
using Sub = BinaryOp<T, N, SubU<T>, A, B>;
template <class T, std::size_t N, class A>
using SSub = ScalarOp<T, N, SubU<T>, A>;
template <class T, std::size_t N, class A>
using ReverseSSub = ReverseScalarOp<T, N, SubU<T>, A>;

template <class T> struct MultU {
  static T apply(T a, T b) { return a * b; }
};
template <class T, std::size_t N, class A, class B>
using Mult = BinaryOp<T, N, MultU<T>, A, B>;
template <class T, std::size_t N, class A>
using SMult = ScalarOp<T, N, MultU<T>, A>;
template <class T, std::size_t N, class A>
using Product = Combine<T, N, MultU<T>, IdU<T>, A>;

template <class T> struct DivU {
  static T apply(T a, T b) { return a / b; }
};
template <class T, std::size_t N, class A, class B>
using Div = BinaryOp<T, N, DivU<T>, A, B>;
template <class T, std::size_t N, class A>
using SDiv = ScalarOp<T, N, DivU<T>, A>;
template <class T, std::size_t N, class A>
using ReverseSDiv = ReverseScalarOp<T, N, DivU<T>, A>;

template <class T> struct NegU {
  static T apply(T a) { return -a; }
};
template <class T, std::size_t N, class A>
using UnaryNeg = UnaryOp<T, N, NegU<T>, A>;
template <class T> struct AbsU {
  static T apply(T a) { return a > 0 ? a : -a; }
};
template <class T, std::size_t N, class A>
using Abs = UnaryOp<T, N, AbsU<T>, A>;

template <class T> struct MinU {
  static T apply(T a, T b) { return a > b ? b : a; }
};
template <class T, std::size_t N, class A, class B>
using ElementwiseMin = BinaryOp<T, N, MinU<T>, A, B>;
template <class T, std::size_t N, class A>
using Min = Combine<T, N, MinU<T>, IdU<T>, A>;

template <class T> struct MaxU {
  static T apply(T a, T b) { return a > b ? a : b; }
};
template <class T, std::size_t N, class A, class B>
using ElementwiseMax = BinaryOp<T, N, MaxU<T>, A, B>;
template <class T, std::size_t N, class A>
using Max = Combine<T, N, MaxU<T>, IdU<T>, A>;

// -----------------------A THING?-----------------------

template <class U, class T, std::size_t N>
using UVVecSpec =
    std::pair<typename rem_cvr<U>::thisisavvec,
              typename std::enable_if<
                  std::is_same<T, typename rem_cvr<U>::value_type>::value &&
                  N == rem_cvr<U>::size>::type>;

template <class A, class B>
using ABVVecSame =
    std::pair<typename rem_cvr<A>::thisisavvec,
              UVVecSpec<B, typename rem_cvr<A>::value_type, rem_cvr<A>::size>>;

// -----------------------FUNCTIONS----------------------

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename rem_cvr<A>::value_type sum(A &&a) {
  return Sum<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>(a)
      .evaluate();
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename rem_cvr<A>::value_type product(A &&a) {
  return Product<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>(a)
      .evaluate();
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename rem_cvr<A>::value_type min(A &&a) {
  return Min<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>(a)
      .evaluate();
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename rem_cvr<A>::value_type max(A &&a) {
  return Max<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>(a)
      .evaluate();
}

template <class A, class B, class = ABVVecSame<A, B>>
ElementwiseMin<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A, B>
elementwiseMin(A &&a, B &&b) {
  return {a, b};
}

template <class A, class B, class = ABVVecSame<A, B>>
ElementwiseMax<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A, B>
elementwiseMax(A &&a, B &&b) {
  return {a, b};
}

template <class A, class B, class = ABVVecSame<A, B>>
typename rem_cvr<A>::value_type dot(A &&a, B &&b) {
  return sum(a * b);
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename rem_cvr<A>::value_type norm2(A &&a) {
  return sum(a * a);
}

template <class A, class B, class = ABVVecSame<A, B>>
typename rem_cvr<A>::value_type dist2(A &&a, B &&b) {
  return norm2(a - b);
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
Abs<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A> vabs(A &&a) {
  return {a};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
typename rem_cvr<A>::value_type norm1(A &&a) {
  return sum(vabs(a));
}

// ----------------------VEC DEFINITION-----------------

template <class T, std::size_t N> struct Vec {

  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  T data[N];

  Vec() {}

  Vec(const std::initializer_list<T> &inil) {
    if (inil.size() != N) {
      return;
    }
    std::size_t i = 0;
    for (const T &val : inil) {
      data[i++] = val;
    }
  }

  template <class U, class = UVVecSpec<U, T, N>> Vec(const U &other) {
    copyFrom<U>(other);
  }

  explicit Vec(const T *other) { copyFrom(other); }

  Vec<T, N> &operator=(const std::initializer_list<T> &inil) {
    if (inil.size() != N) {
      return *this;
    }
    std::size_t i = 0;
    for (const T &val : inil) {
      data[i++] = val;
    }
    return *this;
  }

  template <class U, class = UVVecSpec<U, T, N>>
  Vec<T, N> &operator=(const U &other) {
    copyFrom<U>(other);
    return *this;
  }

  T &operator[](std::size_t i) { return data[i]; }
  const T &operator[](std::size_t i) const { return data[i]; }

  void copyFrom(const T *other) {
    for (std::size_t i = N; i--;) {
      data[i] = other[i];
    }
  }

  template <class U, class = UVVecSpec<U, T, N>> void copyFrom(const U &other) {
    for (std::size_t i = N; i--;) {
      data[i] = other[i];
    }
  }

  template <class U, class = UVVecSpec<U, T, N>>
  Vec<T, N> &operator+=(const U &other) {
    return operator=(*this + other);
  }
  template <class U, class = UVVecSpec<U, T, N>>
  Vec<T, N> &operator-=(const U &other) {
    return operator=(*this - other);
  }
  template <class U, class = UVVecSpec<U, T, N>>
  Vec<T, N> &operator*=(const U &other) {
    return operator=(*this *other);
  }
  template <class U, class = UVVecSpec<U, T, N>>
  Vec<T, N> &operator/=(const U &other) {
    return operator=(*this / other);
  }

  Vec<T, N> &operator+=(T val) { return operator=(*this + val); }
  Vec<T, N> &operator-=(T val) { return operator=(*this - val); }
  Vec<T, N> &operator*=(T val) { return operator=(*this *val); }
  Vec<T, N> &operator/=(T val) { return operator=(*this / val); }
};

template <class T, T value, std::size_t N> struct ConstantVec {

  typedef void thisisavvec;
  typedef T value_type;
  static const std::size_t size = N;

  constexpr T operator[](std::size_t) const { return value; }
};

// -----------------------VEC HELPERS-------------------

template <std::size_t N> using IVec = Vec<std::int32_t, N>;
template <std::size_t N> using FVec = Vec<float, N>;
template <std::size_t N> using DVec = Vec<double, N>;

// Source: https://en.wikipedia.org/wiki/MurmurHash
template <std::size_t N> struct IVecHash {

protected:
  std::uint32_t murmurscram(std::uint32_t val) const noexcept {
    val *= 0xcc9e2d51;
    val = (val << 15) | (val >> 17);
    val *= 0x1b873593;
    return val;
  }

  std::uint32_t murmurstep(std::uint32_t h, std::uint32_t val) const noexcept {
    h ^= murmurscram(val);
    h = (h << 13) | (h >> 19);
    return h * 5 + 0xe6546b64;
  }

public:
  std::uint32_t impl(const std::int32_t *argp) const noexcept {
    return murmurstep(IVecHash<N - 1>().impl(argp), argp[N - 1]);
  }

  std::size_t operator()(const IVec<N> &arg) const noexcept {
    return IVecHash<N - 1>().impl(&arg[0]) ^ murmurscram(arg[N - 1]);
  }
};
template <> struct IVecHash<1> : protected IVecHash<2> {

  std::uint32_t impl(const std::int32_t *argp) const noexcept {
    return murmurstep(0 /* seed */, argp[0]);
  }

  std::size_t operator()(const IVec<1> &arg) const noexcept {
    return murmurscram(arg[0]);
  }
};

} // namespace v

// -----------------------GENERIC OPERATORS-------------

template <class A, class B, class = v::ABVVecSame<A, B>>
v::Add<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A, B>
operator+(A &&a, B &&b) noexcept {
  return {a, b};
}

template <class A, class B, class = v::ABVVecSame<A, B>>
v::Sub<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A, B>
operator-(A &&a, B &&b) noexcept {
  return {a, b};
}

template <class A, class B, class = v::ABVVecSame<A, B>>
v::Mult<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A, B>
operator*(A &&a, B &&b) noexcept {
  return {a, b};
}

template <class A, class B, class = v::ABVVecSame<A, B>>
v::Div<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A, B>
operator/(A &&a, B &&b) noexcept {
  return {a, b};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::UnaryNeg<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator-(A &&a) noexcept {
  return {a};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SAdd<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator+(A &&a, typename rem_cvr<A>::value_type s) noexcept {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SAdd<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator+(typename rem_cvr<A>::value_type s, A &&a) noexcept {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SSub<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator-(A &&a, typename rem_cvr<A>::value_type s) noexcept {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::ReverseSSub<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator-(typename rem_cvr<A>::value_type s, A &&a) noexcept {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SMult<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator*(A &&a, typename rem_cvr<A>::value_type s) noexcept {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SMult<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator*(typename rem_cvr<A>::value_type s, A &&a) noexcept {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::SDiv<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator/(A &&a, typename rem_cvr<A>::value_type s) noexcept {
  return {a, s};
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
v::ReverseSDiv<typename rem_cvr<A>::value_type, rem_cvr<A>::size, A>
operator/(typename rem_cvr<A>::value_type s, A &&a) noexcept {
  return {a, s};
}

template <class A, class B, class = v::ABVVecSame<A, B>>
bool operator==(A &&a, B &&b) noexcept {
  return v::EqualFunctor<A, B>{}(a, b);
}

template <class A, class B, class = v::ABVVecSame<A, B>>
bool operator!=(A &&a, B &&b) noexcept {
  return !(a == b);
}

template <class A, class = typename rem_cvr<A>::thisisavvec>
std::ostream &operator<<(std::ostream &os, A &&a) {
  os << "v[";
  for (std::size_t i = 0; i < rem_cvr<A>::size - 1; i++) {
    os << a[i] << " ";
  }
  os << a[rem_cvr<A>::size - 1] << "] ";
  return os;
}

#endif // V_VECTOR_HPP_
