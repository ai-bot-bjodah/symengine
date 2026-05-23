#include "catch.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include <symengine/polys/groebner.h>
#include <symengine/printers/strprinter.h>

using SymEngine::Basic;
using SymEngine::eq;
using SymEngine::expand;
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

TEST_CASE("Reference Groebner systems from SymPy stay stable", "[groebner][refs]")
{
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");

    require_basis_eq(
        groebner({SymEngine::add({SymEngine::mul(integer(-1), SymEngine::pow(x, integer(2))),
                                  y}),
                  SymEngine::add({SymEngine::mul(integer(-1), SymEngine::pow(x, integer(3))),
                                  z})},
                 vec_basic{x, y, z}, MonomialOrder::Lex),
        {SymEngine::sub(SymEngine::pow(x, integer(2)), y),
         SymEngine::sub(SymEngine::mul(x, y), z),
         SymEngine::sub(SymEngine::mul(x, z), SymEngine::pow(y, integer(2))),
         SymEngine::sub(SymEngine::pow(y, integer(3)), SymEngine::pow(z, integer(2)))});

    require_basis_eq(
        groebner({SymEngine::add({SymEngine::mul(integer(-1), SymEngine::pow(x, integer(2))),
                                  y}),
                  SymEngine::add({SymEngine::mul(integer(-1), SymEngine::pow(x, integer(3))),
                                  z})},
                 vec_basic{x, y, z}, MonomialOrder::GrLex),
        {SymEngine::sub(SymEngine::pow(y, integer(3)), SymEngine::pow(z, integer(2))),
         SymEngine::sub(SymEngine::pow(x, integer(2)), y),
         SymEngine::sub(SymEngine::mul(x, y), z),
         SymEngine::sub(SymEngine::mul(x, z), SymEngine::pow(y, integer(2)))});

    require_basis_eq(
        groebner({SymEngine::sub(x, SymEngine::pow(z, integer(2))),
                  SymEngine::sub(y, SymEngine::pow(z, integer(3)))},
                 vec_basic{x, y, z}, MonomialOrder::Lex),
        {SymEngine::sub(x, SymEngine::pow(z, integer(2))),
         SymEngine::sub(y, SymEngine::pow(z, integer(3)))});

    require_basis_eq(
        groebner({SymEngine::sub(x, SymEngine::pow(y, integer(2))),
                  SymEngine::add({SymEngine::mul(integer(-1), SymEngine::pow(y, integer(3))),
                                  z})},
                 vec_basic{x, y, z}, MonomialOrder::GrLex),
        {SymEngine::sub(SymEngine::pow(x, integer(2)), SymEngine::mul(y, z)),
         SymEngine::sub(SymEngine::mul(x, y), z),
         SymEngine::sub(SymEngine::pow(y, integer(2)), x)});
}
