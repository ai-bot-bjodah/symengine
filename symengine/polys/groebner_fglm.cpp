#include <symengine/polys/groebner.h>

#include <algorithm>
#include <map>
#include <utility>

#include <symengine/polys/groebner_detail.h>
#include <symengine/symengine_exception.h>

namespace SymEngine
{

namespace
{

using detail::ExpressionMatrix;
using detail::ExpressionVector;

vec_int incr_k(const vec_int &monom, unsigned int variable)
{
    vec_int result = monom;
    ++result[variable];
    return result;
}

ExpressionMatrix identity_matrix(size_t size)
{
    ExpressionMatrix result(size, ExpressionVector(size, Expression(0)));
    for (size_t i = 0; i < size; ++i) {
        result[i][i] = Expression(1);
    }
    return result;
}

ExpressionVector matrix_mul(const ExpressionMatrix &matrix,
                            const ExpressionVector &vector)
{
    ExpressionVector result(matrix.size(), Expression(0));
    for (size_t row = 0; row < matrix.size(); ++row) {
        Expression accum(0);
        for (size_t col = 0; col < vector.size(); ++col) {
            accum += detail::normalize_coeff(matrix[row][col] * vector[col]);
        }
        result[row] = detail::normalize_coeff(accum);
    }
    return result;
}

ExpressionMatrix update(size_t s, const ExpressionVector &lambda,
                        ExpressionMatrix P)
{
    size_t pivot = lambda.size();
    for (size_t i = s; i < lambda.size(); ++i) {
        if (not is_zero_coeff(lambda[i])) {
            pivot = i;
            break;
        }
    }
    if (pivot == lambda.size()) {
        throw SymEngineException("fglm expected a non-zero pivot");
    }

    const Expression pivot_value = lambda[pivot];
    for (size_t row = 0; row < P.size(); ++row) {
        if (row == pivot) {
            continue;
        }
        for (size_t col = 0; col < P[row].size(); ++col) {
            P[row][col] = detail::normalize_coeff(
                P[row][col] - (P[pivot][col] * lambda[row]) / pivot_value);
        }
    }

    for (auto &entry : P[pivot]) {
        entry = detail::normalize_coeff(entry / pivot_value);
    }
    std::swap(P[pivot], P[s]);
    return P;
}

bool is_pure_power(const vec_int &monom, unsigned int variable)
{
    if (monom[variable] == 0) {
        return false;
    }
    for (unsigned int i = 0; i < monom.size(); ++i) {
        if (i != variable and monom[i] != 0) {
            return false;
        }
    }
    return true;
}

void require_zero_dimensional(const std::vector<GPoly> &basis)
{
    if (basis.empty()) {
        throw SymEngineException("fglm requires a zero-dimensional ideal");
    }

    const unsigned int nvars = basis.front().nvars;
    for (unsigned int variable = 0; variable < nvars; ++variable) {
        bool found = false;
        for (const auto &poly : basis) {
            if (is_pure_power(poly.LM(), variable)) {
                found = true;
                break;
            }
        }
        if (not found) {
            throw SymEngineException("fglm requires a zero-dimensional ideal");
        }
    }
}

void sort_candidates(std::vector<std::pair<unsigned int, size_t>> &candidates,
                     const std::vector<vec_int> &basis, MonomialOrder order)
{
    std::sort(candidates.begin(), candidates.end(),
              [&basis, order](const std::pair<unsigned int, size_t> &a,
                              const std::pair<unsigned int, size_t> &b) {
                  return monomial_compare(incr_k(basis[a.second], a.first),
                                          incr_k(basis[b.second], b.first), order)
                         > 0;
              });
}

void dedup_candidates(std::vector<std::pair<unsigned int, size_t>> &candidates)
{
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
                     candidates.end());
}

void prune_candidates(std::vector<std::pair<unsigned int, size_t>> &candidates,
                      const std::vector<vec_int> &basis,
                      const std::vector<GPoly> &groebner_basis)
{
    candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                    [&basis, &groebner_basis](
                                        const std::pair<unsigned int, size_t> &t) {
                                        const vec_int monom
                                            = incr_k(basis[t.second], t.first);
                                        for (const auto &poly : groebner_basis) {
                                            if (monomial_divides(poly.LM(), monom)) {
                                                return true;
                                            }
                                        }
                                        return false;
                                    }),
                     candidates.end());
}

GPoly dependency_polynomial(const vec_int &leading,
                            const std::vector<vec_int> &basis,
                            const ExpressionVector &lambda, size_t s,
                            unsigned int nvars, MonomialOrder order)
{
    GPoly poly(nvars, order);
    detail::insert_term(poly, leading, Expression(1));
    for (size_t i = 0; i < s; ++i) {
        if (not is_zero_coeff(lambda[i])) {
            detail::insert_term(poly, basis[i],
                                detail::normalize_coeff(-lambda[i]));
        }
    }
    return poly;
}

} // namespace

GroebnerBasis groebner_fglm(const GroebnerBasis &source, MonomialOrder target)
{
    if (source.basis.empty()) {
        throw SymEngineException("fglm requires a zero-dimensional ideal");
    }
    if (source.basis.size() == 1 and source.basis.front().LM()
                                      == vec_int(source.basis.front().nvars, 0)) {
        return GroebnerBasis{{gpoly_reorder(source.basis.front(), target)},
                             source.vars, target};
    }

    require_zero_dimensional(source.basis);

    const unsigned int nvars = source.basis.front().nvars;
    const std::vector<vec_int> staircase
        = detail::fglm_basis_monomials(source.basis);
    const auto matrices = detail::fglm_representing_matrices(staircase, source.basis);
    const size_t dim = staircase.size();

    std::vector<vec_int> S = {vec_int(nvars, 0)};
    std::vector<ExpressionVector> V = {ExpressionVector(dim, Expression(0))};
    V.front().front() = Expression(1);
    std::vector<GPoly> G;

    std::vector<std::pair<unsigned int, size_t>> L;
    for (unsigned int i = 0; i < nvars; ++i) {
        L.emplace_back(i, 0);
    }
    sort_candidates(L, S, target);

    ExpressionMatrix P = identity_matrix(dim);

    while (true) {
        const std::pair<unsigned int, size_t> t = L.back();
        L.pop_back();

        const size_t s = S.size();
        const ExpressionVector v = matrix_mul(matrices[t.first], V[t.second]);
        const ExpressionVector lambda = matrix_mul(P, v);

        bool dependent = true;
        for (size_t i = s; i < lambda.size(); ++i) {
            if (not is_zero_coeff(lambda[i])) {
                dependent = false;
                break;
            }
        }

        const vec_int candidate = incr_k(S[t.second], t.first);
        if (dependent) {
            GPoly poly
                = dependency_polynomial(candidate, S, lambda, s, nvars, target);
            if (not poly.is_zero()) {
                G.push_back(gpoly_monic(poly));
            }
        } else {
            P = update(s, lambda, std::move(P));
            S.push_back(candidate);
            V.push_back(v);
            for (unsigned int i = 0; i < nvars; ++i) {
                L.emplace_back(i, s);
            }
            dedup_candidates(L);
        }

        prune_candidates(L, S, G);
        if (L.empty()) {
            return GroebnerBasis{detail::red_groebner(G), source.vars, target};
        }
        sort_candidates(L, S, target);
    }
}

} // namespace SymEngine
