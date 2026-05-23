#ifndef SYMENGINE_POLYS_GROEBNER_H
#define SYMENGINE_POLYS_GROEBNER_H

#include <map>
#include <vector>

#include <symengine/expression.h>
#include <symengine/polys/msymenginepoly.h>
#include <symengine/polys/ordering.h>

namespace SymEngine
{

class GPoly
{
public:
    struct Cmp
    {
        MonomialOrder order;

        bool operator()(const vec_int &x, const vec_int &y) const
        {
            return monomial_compare(x, y, order) > 0;
        }
    };

    MonomialOrder order;
    unsigned int nvars;
    std::map<vec_int, Expression, Cmp> terms;

    explicit GPoly(unsigned int nvars_, MonomialOrder order_)
        : order(order_), nvars(nvars_), terms(Cmp{order_})
    {
    }

    bool is_zero() const
    {
        return terms.empty();
    }

    const vec_int &LM() const
    {
        SYMENGINE_ASSERT(not is_zero());
        return terms.begin()->first;
    }

    const Expression &LC() const
    {
        SYMENGINE_ASSERT(not is_zero());
        return terms.begin()->second;
    }

    bool operator==(const GPoly &other) const
    {
        return order == other.order and nvars == other.nvars
               and terms == other.terms;
    }

    bool operator!=(const GPoly &other) const
    {
        return not(*this == other);
    }
};

//! Coefficients are treated as zero iff expand(c) == 0 structurally.
bool is_zero_coeff(const Expression &c);

GPoly gpoly_from_mexprpoly(const RCP<const MExprPoly> &p, const vec_basic &vars,
                           MonomialOrder order);
RCP<const Basic> gpoly_to_basic(const GPoly &p, const vec_basic &vars);
GPoly gpoly_reorder(const GPoly &p, MonomialOrder order);

GPoly gpoly_add(const GPoly &a, const GPoly &b);
GPoly gpoly_sub(const GPoly &a, const GPoly &b);
GPoly gpoly_mul_term(const GPoly &a, const vec_int &m, const Expression &c);
GPoly gpoly_monic(const GPoly &a);
GPoly gpoly_rem(const GPoly &f, const std::vector<GPoly> &G);

std::vector<GPoly> buchberger(std::vector<GPoly> polys);
std::vector<GPoly> f5b(std::vector<GPoly> polys);
std::vector<GPoly> groebner_f4(std::vector<GPoly> polys);

struct GroebnerBasis
{
    std::vector<GPoly> basis;
    vec_basic vars;
    MonomialOrder order;
};

enum class GroebnerAlgorithm {
    Buchberger,
    F5B,
    F4,
    FGLM,
};

struct GroebnerOptions
{
    MonomialOrder order = MonomialOrder::GRevLex;
    GroebnerAlgorithm algorithm = GroebnerAlgorithm::F5B;
    MonomialOrder start_order = MonomialOrder::GRevLex;
};

// Parameters are assumed generic: everything not listed in vars is treated as a
// coefficient in the rational function field K(params). Specializing
// coefficients later can change the basis. F4 and FGLM use structural zero
// tests on Expression coefficients, so rational/numeric coefficients are the
// reliable path and symbolic coefficients are best-effort.
GroebnerBasis groebner_fglm(const GroebnerBasis &source, MonomialOrder target);

GroebnerBasis groebner_basis(const std::vector<GPoly> &polys,
                             const vec_basic &vars,
                             const GroebnerOptions &options);
GroebnerBasis groebner_basis(const std::vector<GPoly> &polys,
                             const vec_basic &vars,
                             MonomialOrder order = MonomialOrder::GRevLex);

// Parameters are assumed generic: everything not listed in vars is treated as a
// coefficient in the rational function field K(params). Specializing
// coefficients later can change the basis. For Buchberger/F5B/F4 the
// computation is performed in `options.order`; the dispatcher rebuilds the
// internal GPolys under that order before invoking the engine. For FGLM the
// input is built under `options.start_order`, reduced there, and converted to
// `options.order`; FGLM is restricted to reduced zero-dimensional ideals.
std::vector<RCP<const Basic>>
groebner(const std::vector<RCP<const Basic>> &exprs, const vec_basic &vars,
         const GroebnerOptions &options);
std::vector<RCP<const Basic>>
groebner(const std::vector<RCP<const Basic>> &exprs, const vec_basic &vars,
         MonomialOrder order = MonomialOrder::GRevLex);

} // namespace SymEngine

#endif
