#ifndef SYMENGINE_POLYS_GROEBNER_DETAIL_H
#define SYMENGINE_POLYS_GROEBNER_DETAIL_H

#include <vector>

#include <symengine/polys/groebner.h>

namespace SymEngine
{

namespace detail
{

struct Term
{
    vec_int monom;
    Expression coeff;
};

using ExpressionVector = std::vector<Expression>;
using ExpressionMatrix = std::vector<ExpressionVector>;

Expression normalize_coeff(const Expression &c);
void insert_term(GPoly &poly, const vec_int &monom, const Expression &coeff);
Term leading_term(const GPoly &poly);
GPoly spoly(const GPoly &f, const GPoly &g);
std::vector<GPoly> reduce_generators(std::vector<GPoly> polys);
std::vector<GPoly> red_groebner(const std::vector<GPoly> &basis);

std::vector<vec_int> fglm_basis_monomials(const std::vector<GPoly> &basis);
std::vector<ExpressionMatrix>
fglm_representing_matrices(const std::vector<vec_int> &basis,
                           const std::vector<GPoly> &groebner_basis);

} // namespace detail

} // namespace SymEngine

#endif
