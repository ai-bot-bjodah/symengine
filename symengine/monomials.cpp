#include <symengine/monomials.h>

namespace SymEngine
{

// This is the fastest implementation:
void monomial_mul(const vec_int &A, const vec_int &B, vec_int &C)
{
    SYMENGINE_ASSERT(A.size() == B.size());
    SYMENGINE_ASSERT(C.size() == A.size());
    size_t n = A.size();
    for (size_t i = 0; i < n; ++i) {
        C[i] = A[i] + B[i];
    }
}

void monomial_lcm(const vec_int &A, const vec_int &B, vec_int &C)
{
    SYMENGINE_ASSERT(A.size() == B.size());
    SYMENGINE_ASSERT(C.size() == A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        C[i] = std::max(A[i], B[i]);
    }
}

void monomial_gcd(const vec_int &A, const vec_int &B, vec_int &C)
{
    SYMENGINE_ASSERT(A.size() == B.size());
    SYMENGINE_ASSERT(C.size() == A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        C[i] = std::min(A[i], B[i]);
    }
}

bool monomial_divides(const vec_int &A, const vec_int &B)
{
    SYMENGINE_ASSERT(A.size() == B.size());
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i] > B[i]) {
            return false;
        }
    }
    return true;
}

bool monomial_div(const vec_int &B, const vec_int &A, vec_int &C)
{
    SYMENGINE_ASSERT(A.size() == B.size());
    SYMENGINE_ASSERT(C.size() == A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i] > B[i]) {
            return false;
        }
        C[i] = B[i] - A[i];
    }
    return true;
}

int monomial_total_degree(const vec_int &A)
{
    int degree = 0;
    for (const auto exponent : A) {
        degree += exponent;
    }
    return degree;
}

/*
// Other implementation of monomial_mul() are below. Those are slightly slower,
// so they are commented out.

// This is slightly slower than monomial_mul
void monomial_mul2(const vec_int &A, const vec_int &B, vec_int &C)
{
    std::transform(A.begin(), A.end(), B.begin(), C.begin(), std::plus<int>());
}

// The same as monomial_mul2
void monomial_mul3(const vec_int &A, const vec_int &B, vec_int &C)
{
    std::transform(A.begin(), A.end(), B.begin(), C.begin(),
        [] (int a, int b) { return a + b; });
}
*/

} // namespace SymEngine
