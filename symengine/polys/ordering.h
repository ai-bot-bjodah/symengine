#ifndef SYMENGINE_POLYS_ORDERING_H
#define SYMENGINE_POLYS_ORDERING_H

#include <symengine/monomials.h>

namespace SymEngine
{

enum class MonomialOrder
{
    Lex,
    GrLex,
    GRevLex
};

// Returns -1 if a < b, 0 if a == b, +1 if a > b under the selected order.
inline int monomial_compare(const vec_int &a, const vec_int &b,
                            MonomialOrder order)
{
    SYMENGINE_ASSERT(a.size() == b.size());

    if (order == MonomialOrder::Lex) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) {
                return a[i] > b[i] ? 1 : -1;
            }
        }
        return 0;
    }

    const int degree_a = monomial_total_degree(a);
    const int degree_b = monomial_total_degree(b);
    if (degree_a != degree_b) {
        return degree_a > degree_b ? 1 : -1;
    }

    if (order == MonomialOrder::GrLex) {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) {
                return a[i] > b[i] ? 1 : -1;
            }
        }
        return 0;
    }

    for (size_t i = a.size(); i-- > 0;) {
        if (a[i] != b[i]) {
            return a[i] < b[i] ? 1 : -1;
        }
    }
    return 0;
}

} // namespace SymEngine

#endif
