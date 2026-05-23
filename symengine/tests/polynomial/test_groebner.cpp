#include "catch.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include <symengine/monomials.h>
#include <symengine/polys/basic_conversions.h>
#include <symengine/polys/groebner.h>
#include <symengine/printers/strprinter.h>

using SymEngine::Basic;
using SymEngine::buchberger;
using SymEngine::eq;
using SymEngine::expand;
using SymEngine::Expression;
using SymEngine::from_basic;
using SymEngine::GPoly;
using SymEngine::gpoly_add;
using SymEngine::gpoly_from_mexprpoly;
using SymEngine::gpoly_monic;
using SymEngine::gpoly_mul_term;
using SymEngine::gpoly_rem;
using SymEngine::gpoly_to_basic;
using SymEngine::groebner;
using SymEngine::integer;
using SymEngine::MExprPoly;
using SymEngine::MonomialOrder;
using SymEngine::monomial_compare;
using SymEngine::monomial_div;
using SymEngine::monomial_divides;
using SymEngine::monomial_gcd;
using SymEngine::monomial_lcm;
using SymEngine::monomial_total_degree;
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

} // namespace

TEST_CASE("Monomial helpers", "[groebner]")
{
    vec_int a = {1, 2, 0};
    vec_int b = {0, 1, 3};
    vec_int out(3, 0);

    monomial_lcm(a, b, out);
    REQUIRE(out == vec_int({1, 2, 3}));

    monomial_gcd(a, b, out);
    REQUIRE(out == vec_int({0, 1, 0}));

    REQUIRE(monomial_divides(vec_int({1, 1, 0}), vec_int({2, 3, 1})));
    REQUIRE(not monomial_divides(vec_int({1, 2, 0}), vec_int({1, 1, 0})));

    REQUIRE(monomial_div(vec_int({2, 3, 1}), vec_int({1, 1, 0}), out));
    REQUIRE(out == vec_int({1, 2, 1}));
    REQUIRE(not monomial_div(vec_int({1, 1, 0}), vec_int({1, 2, 0}), out));

    REQUIRE(monomial_total_degree(vec_int({2, 0, 3})) == 5);
}

TEST_CASE("Monomial orderings", "[groebner]")
{
    REQUIRE(monomial_compare(vec_int({1, 0, 0}), vec_int({0, 2, 0}),
                             MonomialOrder::Lex)
            > 0);
    REQUIRE(monomial_compare(vec_int({0, 2, 0}), vec_int({1, 0, 0}),
                             MonomialOrder::GrLex)
            > 0);
    REQUIRE(monomial_compare(vec_int({1, 2, 0}), vec_int({2, 0, 1}),
                             MonomialOrder::GRevLex)
            > 0);
}

TEST_CASE("GPoly conversion and arithmetic", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto c1 = symbol("c1");
    auto c2 = symbol("c2");
    vec_basic vars = {x, y};

    auto expr = SymEngine::add({
        SymEngine::mul(c1, x),
        y,
    });
    auto poly = as_gpoly(expr, vars, MonomialOrder::Lex);

    REQUIRE(poly.LM() == vec_int({1, 0}));
    REQUIRE(poly.LC() == Expression(c1));
    REQUIRE(same_poly(gpoly_to_basic(poly, vars), expr));

    auto cancel = as_gpoly(SymEngine::mul(integer(-1), SymEngine::mul(c1, x)),
                           vars, MonomialOrder::Lex);
    auto sum = gpoly_add(poly, cancel);
    REQUIRE(same_poly(gpoly_to_basic(sum, vars), y));

    auto multiplied = gpoly_mul_term(poly, vec_int({1, 0}), Expression(c2));
    REQUIRE(same_poly(gpoly_to_basic(multiplied, vars),
                      SymEngine::expand(
                          SymEngine::mul(c2, SymEngine::mul(x, expr)))));

    auto monic = gpoly_monic(as_gpoly(
        SymEngine::add({SymEngine::mul(c1, SymEngine::pow(x, integer(2))),
                        SymEngine::mul(c2, x)}),
        vars, MonomialOrder::Lex));
    REQUIRE(same_poly(gpoly_to_basic(monic, vars),
                      SymEngine::add({SymEngine::pow(x, integer(2)),
                                      SymEngine::mul(SymEngine::div(c2, c1), x)})));
}

TEST_CASE("Multivariate division remainder", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto c1 = symbol("c1");
    auto c2 = symbol("c2");
    vec_basic vars = {x, y};

    auto f = as_gpoly(SymEngine::add(
                          {SymEngine::mul(SymEngine::pow(x, integer(2)), y),
                           SymEngine::mul(x, SymEngine::pow(y, integer(2))),
                           SymEngine::pow(y, integer(2))}),
                      vars, MonomialOrder::Lex);
    std::vector<GPoly> divisors = {
        as_gpoly(SymEngine::sub(SymEngine::mul(x, y), integer(1)), vars,
                 MonomialOrder::Lex),
        as_gpoly(SymEngine::sub(SymEngine::pow(y, integer(2)), integer(1)), vars,
                 MonomialOrder::Lex),
    };

    REQUIRE(same_poly(gpoly_to_basic(gpoly_rem(f, divisors), vars),
                      SymEngine::add({x, y, integer(1)})));

    auto fx = as_gpoly(SymEngine::mul(c1, SymEngine::pow(x, integer(2))),
                       vec_basic{x}, MonomialOrder::Lex);
    std::vector<GPoly> linear = {as_gpoly(SymEngine::sub(x, c2), vec_basic{x},
                                          MonomialOrder::Lex)};
    REQUIRE(same_poly(gpoly_to_basic(gpoly_rem(fx, linear), vec_basic{x}),
                      SymEngine::mul(c1, SymEngine::pow(c2, integer(2)))));
}

TEST_CASE("Buchberger and F5B agree on small systems", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto c1 = symbol("c1");
    auto c2 = symbol("c2");
    auto c3 = symbol("c3");

    const std::vector<std::vector<RCP<const Basic>>> systems = {
        {SymEngine::add({SymEngine::mul(c1, x), SymEngine::mul(c2, y)}),
         SymEngine::add({x, y})},
        {SymEngine::add({SymEngine::pow(x, integer(2)), SymEngine::pow(y, integer(2)),
                         SymEngine::pow(z, integer(2)), integer(-1)}),
         SymEngine::sub(SymEngine::add({SymEngine::pow(x, integer(2)),
                                        SymEngine::pow(z, integer(2))}),
                        y),
         SymEngine::sub(x, z)},
        {SymEngine::add({SymEngine::mul(c1, SymEngine::mul(x, y)),
                         SymEngine::mul(c2,
                                        SymEngine::mul(SymEngine::pow(y, integer(2)),
                                                       z)),
                         SymEngine::mul(c3, SymEngine::pow(z, integer(2)))}),
         SymEngine::add({x, y, z})},
    };

    const std::vector<vec_basic> variables = {
        {x, y},
        {x, y, z},
        {x, y, z},
    };

    for (size_t i = 0; i < systems.size(); ++i) {
        std::vector<GPoly> input;
        for (const auto &expr : systems[i]) {
            input.push_back(as_gpoly(expr, variables[i], MonomialOrder::GRevLex));
        }

        auto buch = basis_to_basic(buchberger(input), variables[i]);
        auto fast = basis_to_basic(SymEngine::f5b(input), variables[i]);
        REQUIRE(same_basis(buch, fast));
    }
}

TEST_CASE("Public groebner API handles symbolic coefficients", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto c1 = symbol("c1");
    auto c2 = symbol("c2");
    auto c3 = symbol("c3");
    vec_basic vars = {x, y};

    auto linear = groebner(
        {SymEngine::add({SymEngine::mul(c1, x), SymEngine::mul(c2, y)})}, vars,
        MonomialOrder::GRevLex);
    require_basis_eq(
        linear, {SymEngine::add({x, SymEngine::mul(SymEngine::div(c2, c1), y)})});

    auto quadratic = groebner(
        {SymEngine::add({SymEngine::mul(c1, SymEngine::pow(x, integer(2))),
                         SymEngine::mul(c2, SymEngine::mul(x, y)),
                         SymEngine::mul(c3, SymEngine::pow(y, integer(2)))}),
         SymEngine::add({x, y})},
        vars, MonomialOrder::Lex);
    require_basis_eq(
        quadratic,
        {SymEngine::add({x, y}), SymEngine::pow(y, integer(2))});

    const vec_basic xyz = {x, y, z};
    auto fractional_field_poly = groebner(
        {SymEngine::add(
             {SymEngine::mul(c1, SymEngine::mul(x, y)),
              SymEngine::mul(c2,
                             SymEngine::mul(SymEngine::pow(y, integer(2)), z)),
              SymEngine::mul(c3, SymEngine::pow(z, integer(2))),
              SymEngine::mul(integer(42), SymEngine::pow(x, integer(2)))})},
        xyz, MonomialOrder::Lex);
    require_basis_eq(
        fractional_field_poly,
        {SymEngine::add(
            {SymEngine::pow(x, integer(2)),
             SymEngine::mul(SymEngine::div(c1, integer(42)),
                            SymEngine::mul(x, y)),
             SymEngine::mul(SymEngine::div(c2, integer(42)),
                            SymEngine::mul(SymEngine::pow(y, integer(2)), z)),
             SymEngine::mul(SymEngine::div(c3, integer(42)),
                            SymEngine::pow(z, integer(2)))})});
}

TEST_CASE("Ordering changes the leading monomial", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto poly = SymEngine::add({x, SymEngine::pow(y, integer(2))});
    vec_basic vars = {x, y};

    auto lex = as_gpoly(poly, vars, MonomialOrder::Lex);
    auto grlex = as_gpoly(poly, vars, MonomialOrder::GrLex);

    REQUIRE(lex.LM() == vec_int({1, 0}));
    REQUIRE(grlex.LM() == vec_int({0, 2}));
}

TEST_CASE("Degenerate groebner inputs", "[groebner]")
{
    auto x = symbol("x");
    auto y = symbol("y");

    REQUIRE(groebner({}, vec_basic{x, y}).empty());

    auto constant = groebner({Expression(5).get_basic()}, vec_basic{x, y});
    require_basis_eq(constant, {integer(1)});

    auto zeros = groebner({integer(0), integer(0)}, vec_basic{x, y});
    REQUIRE(zeros.empty());
}
