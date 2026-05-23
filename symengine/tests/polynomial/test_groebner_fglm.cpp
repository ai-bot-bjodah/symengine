#include "catch.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include <symengine/polys/basic_conversions.h>
#include <symengine/polys/groebner.h>
#include <symengine/polys/groebner_detail.h>
#include <symengine/printers/strprinter.h>

using SymEngine::Basic;
using SymEngine::eq;
using SymEngine::expand;
using SymEngine::Expression;
using SymEngine::from_basic;
using SymEngine::GPoly;
using SymEngine::groebner;
using SymEngine::groebner_basis;
using SymEngine::GroebnerAlgorithm;
using SymEngine::GroebnerOptions;
using SymEngine::gpoly_from_mexprpoly;
using SymEngine::gpoly_to_basic;
using SymEngine::integer;
using SymEngine::MExprPoly;
using SymEngine::MonomialOrder;
using SymEngine::RCP;
using SymEngine::set_basic;
using SymEngine::str;
using SymEngine::symbol;
using SymEngine::vec_basic;
using SymEngine::vec_int;

namespace
{

RCP<const Basic> normalized(const RCP<const Basic> &expr)
{
    return expand(expr);
}

bool same_poly(const RCP<const Basic> &lhs, const RCP<const Basic> &rhs)
{
    return eq(*normalized(lhs), *normalized(rhs));
}

bool same_basis(const std::vector<RCP<const Basic>> &lhs,
                const std::vector<RCP<const Basic>> &rhs)
{
    if (lhs.size() != rhs.size()) {
        return false;
    }

    std::vector<bool> matched(rhs.size(), false);
    for (const auto &poly : lhs) {
        bool found = false;
        for (size_t i = 0; i < rhs.size(); ++i) {
            if (not matched[i] and same_poly(poly, rhs[i])) {
                matched[i] = true;
                found = true;
                break;
            }
        }
        if (not found) {
            return false;
        }
    }
    return true;
}

std::vector<std::string> basis_strings(
    const std::vector<RCP<const Basic>> &basis)
{
    std::vector<std::string> result;
    result.reserve(basis.size());
    for (const auto &poly : basis) {
        result.push_back(str(*normalized(poly)));
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<RCP<const Basic>> basis_to_basic(const std::vector<GPoly> &basis,
                                             const vec_basic &vars)
{
    std::vector<RCP<const Basic>> result;
    result.reserve(basis.size());
    for (const auto &poly : basis) {
        result.push_back(gpoly_to_basic(poly, vars));
    }
    return result;
}

GPoly as_gpoly(const RCP<const Basic> &expr, const vec_basic &vars,
               MonomialOrder order)
{
    set_basic gens(vars.begin(), vars.end());
    return gpoly_from_mexprpoly(from_basic<MExprPoly>(expr, gens, true), vars,
                                order);
}

void require_basis_eq(const std::vector<RCP<const Basic>> &actual,
                      const std::vector<RCP<const Basic>> &expected)
{
    CAPTURE(basis_strings(actual));
    CAPTURE(basis_strings(expected));
    REQUIRE(same_basis(actual, expected));
}

void require_matrix_eq(
    const SymEngine::detail::ExpressionMatrix &actual,
    const std::vector<std::vector<Expression>> &expected)
{
    REQUIRE(actual.size() == expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        REQUIRE(actual[i].size() == expected[i].size());
        for (size_t j = 0; j < expected[i].size(); ++j) {
            REQUIRE(eq(*expand(actual[i][j].get_basic()),
                       *expand(expected[i][j].get_basic())));
        }
    }
}

} // namespace

TEST_CASE("FGLM matches direct F5B on SymPy reference systems", "[groebner][fglm]")
{
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto d = symbol("d");
    auto t = symbol("t");
    auto x = symbol("x");
    auto y = symbol("y");

    struct Case
    {
        std::vector<RCP<const Basic>> exprs;
        vec_basic vars;
        MonomialOrder start_order;
        MonomialOrder target_order;
    };

    const std::vector<Case> cases = {
        {{SymEngine::add({a, b, c, d}),
          SymEngine::add({SymEngine::mul(a, b), SymEngine::mul(a, d),
                          SymEngine::mul(b, c), SymEngine::mul(b, d)}),
          SymEngine::add({SymEngine::mul(SymEngine::mul(a, b), c),
                          SymEngine::mul(SymEngine::mul(a, b), d),
                          SymEngine::mul(SymEngine::mul(a, c), d),
                          SymEngine::mul(SymEngine::mul(b, c), d)}),
          SymEngine::sub(SymEngine::mul(SymEngine::mul(SymEngine::mul(a, b), c), d),
                         integer(1))},
         {a, b, c, d},
         MonomialOrder::GrLex,
         MonomialOrder::Lex},
        {{SymEngine::add(
              {SymEngine::mul(integer(9), SymEngine::pow(x, integer(8))),
               SymEngine::mul(integer(36), SymEngine::pow(x, integer(7))),
               SymEngine::mul(integer(-32), SymEngine::pow(x, integer(6))),
               SymEngine::mul(integer(-252), SymEngine::pow(x, integer(5))),
               SymEngine::mul(integer(-78), SymEngine::pow(x, integer(4))),
               SymEngine::mul(integer(468), SymEngine::pow(x, integer(3))),
               SymEngine::mul(integer(288), SymEngine::pow(x, integer(2))),
               SymEngine::mul(integer(-108), x), integer(9)}),
          SymEngine::add(
              {SymEngine::mul(integer(-72),
                              SymEngine::mul(t, SymEngine::pow(x, integer(7)))),
               SymEngine::mul(integer(-252),
                              SymEngine::mul(t, SymEngine::pow(x, integer(6)))),
               SymEngine::mul(integer(192),
                              SymEngine::mul(t, SymEngine::pow(x, integer(5)))),
               SymEngine::mul(integer(1260),
                              SymEngine::mul(t, SymEngine::pow(x, integer(4)))),
               SymEngine::mul(integer(312),
                              SymEngine::mul(t, SymEngine::pow(x, integer(3)))),
               SymEngine::mul(integer(-404),
                              SymEngine::mul(t, SymEngine::pow(x, integer(2)))),
               SymEngine::mul(integer(-576), SymEngine::mul(t, x)),
               SymEngine::mul(integer(108), t),
               SymEngine::mul(integer(-72), SymEngine::pow(x, integer(7))),
               SymEngine::mul(integer(-256), SymEngine::pow(x, integer(6))),
               SymEngine::mul(integer(192), SymEngine::pow(x, integer(5))),
               SymEngine::mul(integer(1280), SymEngine::pow(x, integer(4))),
               SymEngine::mul(integer(312), SymEngine::pow(x, integer(3))),
               SymEngine::mul(integer(-576), x), integer(96)})},
         {t, x},
         MonomialOrder::GrLex,
         MonomialOrder::Lex},
        {{SymEngine::add({SymEngine::pow(x, integer(2)),
                          SymEngine::mul(integer(-1), x),
                          SymEngine::mul(integer(-3), y), integer(1)}),
          SymEngine::add({SymEngine::mul(integer(-2), x),
                          SymEngine::pow(y, integer(2)), y,
                          integer(-1)})},
         {x, y},
         MonomialOrder::Lex,
         MonomialOrder::GrLex},
    };

    for (const auto &test_case : cases) {
        const auto direct
            = groebner(test_case.exprs, test_case.vars, test_case.target_order);
        const auto converted = groebner(
            test_case.exprs, test_case.vars,
            GroebnerOptions{test_case.target_order, GroebnerAlgorithm::FGLM,
                            test_case.start_order});
        require_basis_eq(converted, direct);
    }
}

TEST_CASE("FGLM representing matrices follow the SymPy reference example",
          "[groebner][fglm]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    vec_basic vars = {x, y};

    std::vector<GPoly> polys = {
        as_gpoly(SymEngine::add({SymEngine::pow(x, integer(2)),
                                 SymEngine::mul(integer(-1), x),
                                 SymEngine::mul(integer(-3), y), integer(1)}),
                 vars, MonomialOrder::GrLex),
        as_gpoly(SymEngine::add({SymEngine::mul(integer(-2), x),
                                 SymEngine::pow(y, integer(2)), y,
                                 integer(-1)}),
                 vars, MonomialOrder::GrLex),
    };

    const auto basis = groebner_basis(polys, vars, MonomialOrder::GrLex);
    const auto staircase = SymEngine::detail::fglm_basis_monomials(basis.basis);
    REQUIRE(staircase
            == std::vector<vec_int>{{0, 0}, {0, 1}, {1, 0}, {1, 1}});

    const auto matrices
        = SymEngine::detail::fglm_representing_matrices(staircase, basis.basis);
    REQUIRE(matrices.size() == 2);

    require_matrix_eq(
        matrices[0],
        {{Expression(0), Expression(0), Expression(-1), Expression(3)},
         {Expression(0), Expression(0), Expression(3), Expression(-4)},
         {Expression(1), Expression(0), Expression(1), Expression(6)},
         {Expression(0), Expression(1), Expression(0), Expression(1)}});
    require_matrix_eq(
        matrices[1],
        {{Expression(0), Expression(1), Expression(0), Expression(-2)},
         {Expression(1), Expression(-1), Expression(0), Expression(6)},
         {Expression(0), Expression(2), Expression(0), Expression(3)},
         {Expression(0), Expression(0), Expression(1), Expression(-1)}});
}

TEST_CASE("FGLM handles failure paths and edge cases", "[groebner][fglm]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    vec_basic xy = {x, y};
    vec_basic xyz = {x, y, z};

    REQUIRE_THROWS(groebner({}, xy, GroebnerOptions{MonomialOrder::Lex,
                                                    GroebnerAlgorithm::FGLM,
                                                    MonomialOrder::GRevLex}));

    REQUIRE_THROWS(groebner({SymEngine::sub(x, SymEngine::pow(z, integer(2))),
                             SymEngine::sub(y, SymEngine::pow(z, integer(3)))},
                            xyz,
                            GroebnerOptions{MonomialOrder::Lex,
                                            GroebnerAlgorithm::FGLM,
                                            MonomialOrder::GRevLex}));

    require_basis_eq(
        groebner({integer(1)}, xy,
                 GroebnerOptions{MonomialOrder::Lex, GroebnerAlgorithm::FGLM,
                                 MonomialOrder::GRevLex}),
        {integer(1)});

    const auto direct = groebner(
        {SymEngine::add({SymEngine::pow(x, integer(2)),
                         SymEngine::mul(integer(-1), x),
                         SymEngine::mul(integer(-3), y), integer(1)}),
         SymEngine::add({SymEngine::mul(integer(-2), x),
                         SymEngine::pow(y, integer(2)), y, integer(-1)})},
        xy, MonomialOrder::GrLex);
    const auto converted = groebner(
        {SymEngine::add({SymEngine::pow(x, integer(2)),
                         SymEngine::mul(integer(-1), x),
                         SymEngine::mul(integer(-3), y), integer(1)}),
         SymEngine::add({SymEngine::mul(integer(-2), x),
                         SymEngine::pow(y, integer(2)), y, integer(-1)})},
        xy, GroebnerOptions{MonomialOrder::GrLex, GroebnerAlgorithm::FGLM,
                            MonomialOrder::GrLex});
    require_basis_eq(converted, direct);
}

TEST_CASE("FGLM works on expand-canonical symbolic coefficients",
          "[groebner][fglm]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto c1 = symbol("c1");
    auto c2 = symbol("c2");
    vec_basic vars = {x, y};

    const std::vector<RCP<const Basic>> exprs = {
        SymEngine::add({SymEngine::mul(c1, x), SymEngine::mul(c2, y)}),
        SymEngine::add({x, y}),
    };

    const auto direct = groebner(exprs, vars, MonomialOrder::Lex);
    const auto converted
        = groebner(exprs, vars,
                   GroebnerOptions{MonomialOrder::Lex, GroebnerAlgorithm::FGLM,
                                   MonomialOrder::GRevLex});
    require_basis_eq(converted, direct);
}
