#include "catch.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include <symengine/polys/groebner.h>
#include <symengine/printers/strprinter.h>

using SymEngine::Basic;
using SymEngine::eq;
using SymEngine::expand;
using SymEngine::GroebnerAlgorithm;
using SymEngine::GroebnerOptions;
using SymEngine::groebner;
using SymEngine::integer;
using SymEngine::MonomialOrder;
using SymEngine::RCP;
using SymEngine::str;
using SymEngine::symbol;
using SymEngine::vec_basic;

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

void require_basis_eq(const std::vector<RCP<const Basic>> &actual,
                      const std::vector<RCP<const Basic>> &expected)
{
    CAPTURE(basis_strings(actual));
    CAPTURE(basis_strings(expected));
    REQUIRE(same_basis(actual, expected));
}

} // namespace

TEST_CASE("F4 matches F5B on numeric benchmark systems", "[groebner][f4]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto u0 = symbol("u0");
    auto u1 = symbol("u1");
    auto u2 = symbol("u2");
    auto u3 = symbol("u3");

    struct Case
    {
        std::vector<RCP<const Basic>> exprs;
        vec_basic vars;
        MonomialOrder order;
    };

    const std::vector<Case> cases = {
        {{SymEngine::add({SymEngine::pow(x, integer(2)),
                          SymEngine::mul(integer(2),
                                         SymEngine::mul(x, SymEngine::pow(y, integer(2))))}),
          SymEngine::add({SymEngine::mul(x, y),
                          SymEngine::mul(integer(2), SymEngine::pow(y, integer(3))),
                          integer(-1)})},
         {x, y},
         MonomialOrder::Lex},
        {{SymEngine::add({SymEngine::pow(x, integer(2)), SymEngine::pow(y, integer(2)),
                          SymEngine::pow(z, integer(2)), integer(-1)}),
          SymEngine::sub(SymEngine::add({SymEngine::pow(x, integer(2)),
                                         SymEngine::pow(z, integer(2))}),
                         y),
          SymEngine::sub(x, z)},
         {x, y, z},
         MonomialOrder::GRevLex},
        {{SymEngine::add({u0, SymEngine::mul(integer(2), u1),
                          SymEngine::mul(integer(2), u2),
                          SymEngine::mul(integer(2), u3), integer(-1)}),
          SymEngine::add({SymEngine::pow(u1, integer(2)),
                          SymEngine::mul(integer(2), SymEngine::mul(u0, u1)),
                          SymEngine::mul(integer(2), SymEngine::mul(u2, u3)),
                          SymEngine::mul(integer(-1), u1)}),
          SymEngine::add({SymEngine::pow(u2, integer(2)),
                          SymEngine::mul(integer(2), SymEngine::mul(u0, u2)),
                          SymEngine::mul(integer(2), SymEngine::mul(u1, u3)),
                          SymEngine::mul(integer(-1), u2)}),
          SymEngine::add({SymEngine::pow(u3, integer(2)),
                          SymEngine::mul(integer(2), SymEngine::mul(u0, u3)),
                          SymEngine::mul(integer(-1), u3)})},
         {u0, u1, u2, u3},
         MonomialOrder::GRevLex},
    };

    for (const auto &test_case : cases) {
        const auto direct = groebner(test_case.exprs, test_case.vars, test_case.order);
        const auto f4 = groebner(
            test_case.exprs, test_case.vars,
            GroebnerOptions{test_case.order, GroebnerAlgorithm::F4,
                            test_case.order});
        require_basis_eq(f4, direct);
    }
}

TEST_CASE("F4 handles best-effort symbolic coefficients and degenerate inputs",
          "[groebner][f4]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto c1 = symbol("c1");
    auto c2 = symbol("c2");
    vec_basic vars = {x, y};

    const std::vector<RCP<const Basic>> symbolic = {
        SymEngine::add({SymEngine::mul(c1, x), SymEngine::mul(c2, y)}),
        SymEngine::add({x, y}),
    };

    require_basis_eq(
        groebner(symbolic, vars,
                 GroebnerOptions{MonomialOrder::GRevLex,
                                 GroebnerAlgorithm::F4,
                                 MonomialOrder::GRevLex}),
        groebner(symbolic, vars, MonomialOrder::GRevLex));

    REQUIRE(groebner({}, vars,
                     GroebnerOptions{MonomialOrder::Lex,
                                     GroebnerAlgorithm::F4,
                                     MonomialOrder::Lex})
                .empty());

    require_basis_eq(groebner({integer(1)}, vars,
                              GroebnerOptions{MonomialOrder::Lex,
                                              GroebnerAlgorithm::F4,
                                              MonomialOrder::Lex}),
                     {integer(1)});
}
