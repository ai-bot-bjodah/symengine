#include <symengine/polys/groebner.h>

#include <algorithm>
#include <cstddef>
#include <map>
#include <utility>

#include <symengine/polys/basic_conversions.h>
#include <symengine/symengine_exception.h>

namespace SymEngine
{

namespace
{

struct Term
{
    vec_int monom;
    Expression coeff;
};

struct Signature
{
    vec_int monom;
    unsigned int index;

    bool operator==(const Signature &other) const
    {
        return index == other.index and monom == other.monom;
    }
};

struct LabeledPoly
{
    Signature sig;
    GPoly poly;
    unsigned int num;

    bool operator==(const LabeledPoly &other) const
    {
        return sig == other.sig and poly == other.poly and num == other.num;
    }
};

struct CriticalPair
{
    Signature first_sig;
    Term first_mul;
    LabeledPoly first_poly;
    Signature second_sig;
    Term second_mul;
    LabeledPoly second_poly;
};

void assert_compatible(const GPoly &a, const GPoly &b)
{
    SYMENGINE_ASSERT(a.nvars == b.nvars);
    SYMENGINE_ASSERT(a.order == b.order);
}

Expression expanded_coeff(const Expression &c)
{
    return expand(c);
}

void insert_term(GPoly &poly, const vec_int &monom, const Expression &coeff)
{
    SYMENGINE_ASSERT(monom.size() == poly.nvars);

    Expression expanded = expanded_coeff(coeff);
    if (is_zero_coeff(expanded)) {
        return;
    }

    auto it = poly.terms.find(monom);
    if (it == poly.terms.end()) {
        poly.terms.insert(std::make_pair(monom, expanded));
    } else {
        it->second = expanded_coeff(it->second + expanded);
        if (is_zero_coeff(it->second)) {
            poly.terms.erase(it);
        }
    }
}

Term leading_term(const GPoly &poly)
{
    return Term{poly.LM(), poly.LC()};
}

GPoly leading_term_poly(const GPoly &poly)
{
    GPoly result(poly.nvars, poly.order);
    insert_term(result, poly.LM(), poly.LC());
    return result;
}

bool gpoly_lm_less(const GPoly &a, const GPoly &b)
{
    return monomial_compare(a.LM(), b.LM(), a.order) < 0;
}

void sort_by_lm_desc(std::vector<GPoly> &polys)
{
    std::sort(polys.begin(), polys.end(),
              [](const GPoly &a, const GPoly &b) { return gpoly_lm_less(b, a); });
}

void sort_labeled_by_lm_desc(std::vector<LabeledPoly> &basis)
{
    std::sort(basis.begin(), basis.end(),
              [](const LabeledPoly &a, const LabeledPoly &b) {
                  return gpoly_lm_less(b.poly, a.poly);
              });
}

GPoly negate_gpoly(const GPoly &poly)
{
    GPoly result(poly.nvars, poly.order);
    for (const auto &term : poly.terms) {
        insert_term(result, term.first, -term.second);
    }
    return result;
}

Term term_div(const Term &numerator, const Term &denominator)
{
    Term quotient = {vec_int(numerator.monom.size(), 0), Expression(0)};
    const bool divisible
        = monomial_div(numerator.monom, denominator.monom, quotient.monom);
    SYMENGINE_ASSERT(divisible);
    quotient.coeff = expanded_coeff(numerator.coeff / denominator.coeff);
    return quotient;
}

GPoly spoly(const GPoly &f, const GPoly &g)
{
    assert_compatible(f, g);
    vec_int lcm(f.nvars, 0);
    monomial_lcm(f.LM(), g.LM(), lcm);

    vec_int mf(f.nvars, 0), mg(f.nvars, 0);
    const bool f_divides = monomial_div(lcm, f.LM(), mf);
    const bool g_divides = monomial_div(lcm, g.LM(), mg);
    SYMENGINE_ASSERT(f_divides and g_divides);

    return gpoly_sub(gpoly_mul_term(f, mf, Expression(1) / f.LC()),
                     gpoly_mul_term(g, mg, Expression(1) / g.LC()));
}

std::vector<GPoly> reduce_generators(std::vector<GPoly> polys)
{
    if (polys.empty()) {
        return polys;
    }

    std::vector<GPoly> current;
    current.reserve(polys.size());
    for (const auto &poly : polys) {
        if (not poly.is_zero()) {
            current.push_back(gpoly_monic(poly));
        }
    }

    while (true) {
        std::vector<GPoly> next;
        next.reserve(current.size());

        for (size_t i = 0; i < current.size(); ++i) {
            std::vector<GPoly> reducers(current.begin(),
                                        current.begin()
                                            + static_cast<ptrdiff_t>(i));
            GPoly remainder = gpoly_rem(current[i], reducers);
            if (not remainder.is_zero()) {
                next.push_back(gpoly_monic(remainder));
            }
        }

        if (next == current) {
            return next;
        }
        current.swap(next);
    }
}

std::vector<GPoly> red_groebner(const std::vector<GPoly> &basis)
{
    std::vector<GPoly> F;
    for (const auto &poly : basis) {
        if (not poly.is_zero()) {
            F.push_back(gpoly_monic(poly));
        }
    }

    std::vector<GPoly> H;
    while (not F.empty()) {
        GPoly f0 = F.back();
        F.pop_back();

        bool redundant = false;
        for (const auto &poly : F) {
            if (monomial_divides(poly.LM(), f0.LM())) {
                redundant = true;
                break;
            }
        }
        if (not redundant) {
            for (const auto &poly : H) {
                if (monomial_divides(poly.LM(), f0.LM())) {
                    redundant = true;
                    break;
                }
            }
        }
        if (not redundant) {
            H.push_back(f0);
        }
    }

    std::vector<GPoly> reduced;
    for (size_t i = 0; i < H.size(); ++i) {
        std::vector<GPoly> others;
        others.reserve(H.size() - 1);
        for (size_t j = 0; j < H.size(); ++j) {
            if (i != j) {
                others.push_back(H[j]);
            }
        }
        GPoly remainder = gpoly_rem(H[i], others);
        if (not remainder.is_zero()) {
            reduced.push_back(gpoly_monic(remainder));
        }
    }

    sort_by_lm_desc(reduced);
    return reduced;
}

Signature sig_mult(const Signature &sig, const vec_int &monom)
{
    vec_int product(sig.monom.size(), 0);
    monomial_mul(sig.monom, monom, product);
    return Signature{product, sig.index};
}

int sig_compare(const Signature &u, const Signature &v, MonomialOrder order)
{
    if (u.index != v.index) {
        return u.index > v.index ? -1 : 1;
    }

    const int cmp = monomial_compare(u.monom, v.monom, order);
    if (cmp < 0) {
        return -1;
    }
    if (cmp > 0) {
        return 1;
    }
    return 0;
}

int lbp_compare(const LabeledPoly &f, const LabeledPoly &g)
{
    const int sig_cmp = sig_compare(f.sig, g.sig, f.poly.order);
    if (sig_cmp != 0) {
        return sig_cmp;
    }
    if (f.num > g.num) {
        return -1;
    }
    if (f.num < g.num) {
        return 1;
    }
    return 0;
}

LabeledPoly lbp_mul_term(const LabeledPoly &f, const Term &term)
{
    return LabeledPoly{
        sig_mult(f.sig, term.monom),
        gpoly_mul_term(f.poly, term.monom, term.coeff),
        f.num,
    };
}

LabeledPoly lbp_sub(const LabeledPoly &f, const LabeledPoly &g)
{
    const LabeledPoly &max_poly = lbp_compare(f, g) < 0 ? g : f;
    return LabeledPoly{max_poly.sig, gpoly_sub(f.poly, g.poly), max_poly.num};
}

CriticalPair critical_pair(const LabeledPoly &f, const LabeledPoly &g)
{
    vec_int lcm(f.poly.nvars, 0);
    monomial_lcm(f.poly.LM(), g.poly.LM(), lcm);

    const Term lt = {lcm, Expression(1)};
    const Term uf = term_div(lt, leading_term(f.poly));
    const Term vg = term_div(lt, leading_term(g.poly));

    const LabeledPoly fr = lbp_mul_term(
        LabeledPoly{f.sig, leading_term_poly(f.poly), f.num}, uf);
    const LabeledPoly gr = lbp_mul_term(
        LabeledPoly{g.sig, leading_term_poly(g.poly), g.num}, vg);

    if (lbp_compare(fr, gr) < 0) {
        return CriticalPair{gr.sig, vg, g, fr.sig, uf, f};
    }
    return CriticalPair{fr.sig, uf, f, gr.sig, vg, g};
}

int cp_compare(const CriticalPair &c, const CriticalPair &d)
{
    const LabeledPoly c0 = {c.first_sig, GPoly(c.first_poly.poly.nvars,
                                               c.first_poly.poly.order),
                            c.first_poly.num};
    const LabeledPoly d0 = {d.first_sig, GPoly(d.first_poly.poly.nvars,
                                               d.first_poly.poly.order),
                            d.first_poly.num};

    const int first = lbp_compare(c0, d0);
    if (first != 0) {
        return first;
    }

    const LabeledPoly c1 = {c.second_sig, GPoly(c.second_poly.poly.nvars,
                                                c.second_poly.poly.order),
                            c.second_poly.num};
    const LabeledPoly d1 = {d.second_sig, GPoly(d.second_poly.poly.nvars,
                                                d.second_poly.poly.order),
                            d.second_poly.num};
    return lbp_compare(c1, d1);
}

void sort_critical_pairs_desc(std::vector<CriticalPair> &pairs)
{
    std::sort(pairs.begin(), pairs.end(), [](const CriticalPair &a,
                                             const CriticalPair &b) {
        return cp_compare(a, b) > 0;
    });
}

LabeledPoly s_poly(const CriticalPair &cp)
{
    return lbp_sub(
        lbp_mul_term(cp.first_poly, cp.first_mul),
        lbp_mul_term(cp.second_poly, cp.second_mul));
}

bool is_rewritable_or_comparable(const Signature &sign, unsigned int num,
                                 const std::vector<LabeledPoly> &basis)
{
    for (const auto &poly : basis) {
        if (sign.index < poly.sig.index
            and monomial_divides(poly.poly.LM(), sign.monom)) {
            return true;
        }
        if (sign.index == poly.sig.index and num < poly.num
            and monomial_divides(poly.sig.monom, sign.monom)) {
            return true;
        }
    }
    return false;
}

LabeledPoly f5_reduce_impl(LabeledPoly poly, const std::vector<LabeledPoly> &basis)
{
    if (poly.poly.is_zero()) {
        return poly;
    }

    while (true) {
        bool changed = false;

        for (const auto &candidate : basis) {
            if (candidate.poly.is_zero()) {
                continue;
            }
            if (not monomial_divides(candidate.poly.LM(), poly.poly.LM())) {
                continue;
            }

            const Term multiplier
                = term_div(leading_term(poly.poly), leading_term(candidate.poly));
            if (sig_compare(sig_mult(candidate.sig, multiplier.monom), poly.sig,
                            poly.poly.order)
                < 0) {
                poly = lbp_sub(poly, lbp_mul_term(candidate, multiplier));
                changed = true;
                break;
            }
        }

        if (not changed or poly.poly.is_zero()) {
            return poly;
        }
    }
}

} // namespace

bool is_zero_coeff(const Expression &c)
{
    return expand(c) == Expression(0);
}

GPoly gpoly_from_mexprpoly(const RCP<const MExprPoly> &p, const vec_basic &vars,
                           MonomialOrder order)
{
    GPoly result(numeric_cast<unsigned int>(vars.size()), order);

    std::map<RCP<const Basic>, unsigned int, RCPBasicKeyLess> user_positions;
    for (unsigned int i = 0; i < vars.size(); ++i) {
        auto inserted = user_positions.insert(std::make_pair(vars[i], i));
        if (not inserted.second) {
            throw SymEngineException("groebner variables must be distinct");
        }
    }

    std::vector<unsigned int> positions;
    positions.reserve(p->get_vars().size());
    for (const auto &var : p->get_vars()) {
        auto it = user_positions.find(var);
        if (it == user_positions.end()) {
            throw SymEngineException(
                "polynomial variable missing from requested Groebner variable order");
        }
        positions.push_back(it->second);
    }

    for (const auto &term : p->get_poly().get_dict()) {
        vec_int exponents(vars.size(), 0);
        SYMENGINE_ASSERT(term.first.size() == positions.size());
        for (size_t i = 0; i < positions.size(); ++i) {
            exponents[positions[i]] = term.first[i];
        }
        insert_term(result, exponents, term.second);
    }

    return result;
}

RCP<const Basic> gpoly_to_basic(const GPoly &p, const vec_basic &vars)
{
    SYMENGINE_ASSERT(vars.size() == p.nvars);

    if (p.is_zero()) {
        return zero;
    }

    vec_basic terms;
    for (const auto &term : p.terms) {
        RCP<const Basic> expr = term.second.get_basic();
        for (size_t i = 0; i < vars.size(); ++i) {
            if (term.first[i] != 0) {
                expr = mul(expr, pow(vars[i], integer(term.first[i])));
            }
        }
        terms.push_back(expr);
    }
    return add(terms);
}

GPoly gpoly_add(const GPoly &a, const GPoly &b)
{
    assert_compatible(a, b);
    GPoly result = a;
    for (const auto &term : b.terms) {
        insert_term(result, term.first, term.second);
    }
    return result;
}

GPoly gpoly_sub(const GPoly &a, const GPoly &b)
{
    return gpoly_add(a, negate_gpoly(b));
}

GPoly gpoly_mul_term(const GPoly &a, const vec_int &m, const Expression &c)
{
    SYMENGINE_ASSERT(m.size() == a.nvars);
    GPoly result(a.nvars, a.order);
    for (const auto &term : a.terms) {
        vec_int monomial(a.nvars, 0);
        monomial_mul(term.first, m, monomial);
        insert_term(result, monomial, term.second * c);
    }
    return result;
}

GPoly gpoly_monic(const GPoly &a)
{
    SYMENGINE_ASSERT(not a.is_zero());
    GPoly result(a.nvars, a.order);
    const Expression leading = a.LC();
    for (const auto &term : a.terms) {
        insert_term(result, term.first, term.second / leading);
    }
    return result;
}

GPoly gpoly_rem(const GPoly &f, const std::vector<GPoly> &G)
{
    GPoly remainder(f.nvars, f.order);
    GPoly current = f;

    while (not current.is_zero()) {
        bool divided = false;

        for (const auto &divisor : G) {
            if (divisor.is_zero()) {
                continue;
            }
            assert_compatible(current, divisor);
            if (monomial_divides(divisor.LM(), current.LM())) {
                vec_int quotient(current.nvars, 0);
                const bool ok = monomial_div(current.LM(), divisor.LM(), quotient);
                SYMENGINE_ASSERT(ok);
                const Expression coeff
                    = expanded_coeff(current.LC() / divisor.LC());
                current = gpoly_sub(current,
                                    gpoly_mul_term(divisor, quotient, coeff));
                divided = true;
                break;
            }
        }

        if (not divided) {
            auto leading = current.terms.begin();
            insert_term(remainder, leading->first, leading->second);
            current.terms.erase(leading);
        }
    }

    return remainder;
}

std::vector<GPoly> buchberger(std::vector<GPoly> polys)
{
    polys = reduce_generators(std::move(polys));
    if (polys.empty()) {
        return polys;
    }

    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t i = 0; i < polys.size(); ++i) {
        for (size_t j = i + 1; j < polys.size(); ++j) {
            pairs.push_back(std::make_pair(i, j));
        }
    }

    while (not pairs.empty()) {
        const auto pair = pairs.back();
        pairs.pop_back();

        GPoly remainder
            = gpoly_rem(spoly(polys[pair.first], polys[pair.second]), polys);
        if (remainder.is_zero()) {
            continue;
        }

        remainder = gpoly_monic(remainder);
        const size_t new_index = polys.size();
        for (size_t i = 0; i < new_index; ++i) {
            pairs.push_back(std::make_pair(i, new_index));
        }
        polys.push_back(remainder);
    }

    return red_groebner(polys);
}

std::vector<GPoly> f5b(std::vector<GPoly> polys)
{
    std::vector<GPoly> reduced = reduce_generators(std::move(polys));
    if (reduced.empty()) {
        return reduced;
    }

    const unsigned int nvars = reduced.front().nvars;
    const MonomialOrder order = reduced.front().order;
    const vec_int zero_monom(nvars, 0);

    std::vector<LabeledPoly> basis;
    basis.reserve(reduced.size());
    for (unsigned int i = 0; i < reduced.size(); ++i) {
        basis.push_back(
            LabeledPoly{Signature{zero_monom, i + 1}, reduced[i], i + 1});
    }
    sort_labeled_by_lm_desc(basis);

    std::vector<CriticalPair> pairs;
    for (size_t i = 0; i < basis.size(); ++i) {
        for (size_t j = i + 1; j < basis.size(); ++j) {
            pairs.push_back(critical_pair(basis[i], basis[j]));
        }
    }
    sort_critical_pairs_desc(pairs);

    unsigned int k = numeric_cast<unsigned int>(basis.size());
    while (not pairs.empty()) {
        CriticalPair cp = pairs.back();
        pairs.pop_back();

        if (is_rewritable_or_comparable(cp.first_sig, cp.first_poly.num, basis)
            or is_rewritable_or_comparable(cp.second_sig, cp.second_poly.num,
                                           basis)) {
            continue;
        }

        LabeledPoly poly = f5_reduce_impl(s_poly(cp), basis);
        if (poly.poly.is_zero()) {
            continue;
        }

        poly.poly = gpoly_monic(poly.poly);
        poly.num = k + 1;

        std::vector<size_t> redundant_indices;
        for (size_t i = 0; i < pairs.size(); ++i) {
            if (is_rewritable_or_comparable(pairs[i].first_sig,
                                            pairs[i].first_poly.num,
                                            std::vector<LabeledPoly>{poly})
                or is_rewritable_or_comparable(
                    pairs[i].second_sig, pairs[i].second_poly.num,
                    std::vector<LabeledPoly>{poly})) {
                redundant_indices.push_back(i);
            }
        }
        for (size_t i = redundant_indices.size(); i-- > 0;) {
            pairs.erase(pairs.begin()
                        + static_cast<ptrdiff_t>(redundant_indices[i]));
        }

        for (const auto &candidate : basis) {
            if (candidate.poly.is_zero()) {
                continue;
            }
            CriticalPair new_pair = critical_pair(poly, candidate);
            if (is_rewritable_or_comparable(new_pair.first_sig,
                                            new_pair.first_poly.num,
                                            std::vector<LabeledPoly>{poly})
                or is_rewritable_or_comparable(new_pair.second_sig,
                                               new_pair.second_poly.num,
                                               std::vector<LabeledPoly>{poly})) {
                continue;
            }
            pairs.push_back(new_pair);
        }
        sort_critical_pairs_desc(pairs);

        if (basis.empty()
            or monomial_compare(poly.poly.LM(), basis.back().poly.LM(), order)
                   <= 0) {
            basis.push_back(poly);
        } else {
            auto it = basis.begin();
            while (it != basis.end()
                   and monomial_compare(poly.poly.LM(), it->poly.LM(), order)
                           <= 0) {
                ++it;
            }
            basis.insert(it, poly);
        }
        ++k;
    }

    std::vector<GPoly> raw_basis;
    for (const auto &poly : basis) {
        if (not poly.poly.is_zero()) {
            raw_basis.push_back(gpoly_monic(poly.poly));
        }
    }
    return red_groebner(raw_basis);
}

GroebnerBasis groebner_basis(const std::vector<GPoly> &polys,
                             const vec_basic &vars, MonomialOrder order)
{
    GroebnerBasis result;
    result.vars = vars;
    result.order = order;
    result.basis = f5b(polys);
    return result;
}

std::vector<RCP<const Basic>>
groebner(const std::vector<RCP<const Basic>> &exprs, const vec_basic &vars,
         MonomialOrder order)
{
    set_basic gens(vars.begin(), vars.end());
    if (gens.size() != vars.size()) {
        throw SymEngineException("groebner variables must be distinct");
    }

    std::vector<GPoly> polys;
    polys.reserve(exprs.size());
    for (const auto &expr : exprs) {
        auto poly = gpoly_from_mexprpoly(
            from_basic<MExprPoly>(expr, gens, true), vars, order);
        if (not poly.is_zero()) {
            polys.push_back(poly);
        }
    }

    const auto basis = groebner_basis(polys, vars, order);
    std::vector<RCP<const Basic>> result;
    result.reserve(basis.basis.size());
    for (const auto &poly : basis.basis) {
        result.push_back(gpoly_to_basic(poly, vars));
    }
    return result;
}

} // namespace SymEngine
