/**
 *  \file monomials.h
 *  Monomial Multiplication
 *
 **/

#ifndef SYMENGINE_MONOMIALS_H
#define SYMENGINE_MONOMIALS_H

#include <symengine/basic.h>

namespace SymEngine
{
//! Monomial multiplication
void monomial_mul(const vec_int &A, const vec_int &B, vec_int &C);

//! Least common multiple of two monomials
void monomial_lcm(const vec_int &A, const vec_int &B, vec_int &C);

//! Greatest common divisor of two monomials
void monomial_gcd(const vec_int &A, const vec_int &B, vec_int &C);

//! True iff A divides B, i.e. A[i] <= B[i] for all i
bool monomial_divides(const vec_int &A, const vec_int &B);

//! If A divides B, set C = B / A and return true, else return false
bool monomial_div(const vec_int &B, const vec_int &A, vec_int &C);

//! Total degree of a monomial
int monomial_total_degree(const vec_int &A);

} // namespace SymEngine

#endif
