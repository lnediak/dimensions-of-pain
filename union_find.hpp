#ifndef UNION_FIND_HPP_
#define UNION_FIND_HPP_

namespace uf {

/*
  GetP is functor, int &getP(int i)
  GetSz is functor, int &getSz(int i)
*/

/// find representative index of the set i is in
template <class GetP> int ufFind(int i, GetP &&getP) {
  // path halving
  while (getP(i) != i) {
    i = getP(i) = getP(getP(i));
  }
  return i;
}
/// returns true if a union was performed
template <class GetP, class GetSz>
bool ufUnion(int i, int j, GetP &&getP, GetSz &&getSz) {
  int pi = ufFind(i, getP);
  int pj = ufFind(j, getP);
  if (pi == pj) {
    return false;
  }
  if (getSz(pi) < getSz(pj)) {
    int tmp = pi;
    pi = pj;
    pj = tmp;
  }
  getP(pj) = pi;
  getSz(pi) += getSz(pj);
  return true;
}

} // namespace uf

#endif // UNION_FIND_HPP_
