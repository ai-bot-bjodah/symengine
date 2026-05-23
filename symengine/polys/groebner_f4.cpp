#include <symengine/polys/groebner.h>

#include <algorithm>
#include <map>
#include <set>
#include <tuple>
#include <utility>

#include <symengine/matrix.h>
#include <symengine/polys/groebner_detail.h>

namespace SymEngine
{

namespace
{

enum class MonomialState {
    Pivot,
    Done,
    Todo,
};

struct F4Pair
{
    size_t first;
    size_t second;
    vec_int lcm;
    unsigned int degree;
};

struct PreprocessResult
{
    std::vector<GPoly> rows;
    std::set<vec_int> known_pivots;
};

bool monomial_less(const vec_int &a, const vec_int &b, MonomialOrder order)
{
    return monomial_compare(a, b, order) < 0;
}

bool monomial_greater(const vec_int &a, const vec_int &b, MonomialOrder order)
{
    return monomial_compare(a, b, order) > 0;
}

std::vector<F4Pair> initial_pairs(const std::vector<GPoly> &basis)
{
    std::vector<F4Pair> pairs;
    for (size_t i = 0; i < basis.size(); ++i) {
        for (size_t j = i + 1; j < basis.size(); ++j) {
            vec_int lcm(basis[i].nvars, 0);
            monomial_lcm(basis[i].LM(), basis[j].LM(), lcm);
            pairs.push_back(F4Pair{
                i, j, lcm, static_cast<unsigned int>(monomial_total_degree(lcm))});
        }
    }
    return pairs;
}

std::pair<std::vector<F4Pair>, std::vector<F4Pair>>
select_minimal_degree_pairs(std::vector<F4Pair> pairs, MonomialOrder order)
{
    if (pairs.empty()) {
        return {{}, {}};
    }

    unsigned int min_degree = pairs.front().degree;
    for (const auto &pair : pairs) {
        min_degree = std::min(min_degree, pair.degree);
    }

    std::vector<F4Pair> selected;
    std::vector<F4Pair> remaining;
    for (auto &pair : pairs) {
        if (pair.degree == min_degree) {
            selected.push_back(std::move(pair));
        } else {
            remaining.push_back(std::move(pair));
        }
    }

    std::sort(selected.begin(), selected.end(),
              [order](const F4Pair &a, const F4Pair &b) {
                  const int cmp = monomial_compare(a.lcm, b.lcm, order);
                  if (cmp != 0) {
                      return cmp > 0;
                  }
                  if (a.first != b.first) {
                      return a.first < b.first;
                  }
                  return a.second < b.second;
              });
    return {selected, remaining};
}

GPoly multiplied_row(const GPoly &poly, const vec_int &leading_lcm)
{
    vec_int quotient(poly.nvars, 0);
    const bool divisible = monomial_div(leading_lcm, poly.LM(), quotient);
    SYMENGINE_ASSERT(divisible);
    return gpoly_mul_term(poly, quotient, Expression(1) / poly.LC());
}

std::vector<vec_int> collect_row_tail(const GPoly &poly)
{
    std::vector<vec_int> monomials;
    bool first = true;
    for (const auto &term : poly.terms) {
        if (first) {
            first = false;
            continue;
        }
        monomials.push_back(term.first);
    }
    return monomials;
}

vec_int next_todo(const std::map<vec_int, MonomialState> &states,
                  MonomialOrder order)
{
    bool found = false;
    vec_int result;
    for (const auto &entry : states) {
        if (entry.second != MonomialState::Todo) {
            continue;
        }
        if (not found or monomial_greater(entry.first, result, order)) {
            result = entry.first;
            found = true;
        }
    }
    SYMENGINE_ASSERT(found);
    return result;
}

PreprocessResult symbolic_preprocess(const std::vector<F4Pair> &pairs,
                                     const std::vector<GPoly> &basis)
{
    PreprocessResult result;
    if (pairs.empty()) {
        return result;
    }

    const MonomialOrder order = basis.front().order;
    std::map<vec_int, std::set<size_t>> grouped;
    for (const auto &pair : pairs) {
        grouped[pair.lcm].insert(pair.first);
        grouped[pair.lcm].insert(pair.second);
    }

    std::map<vec_int, MonomialState> states;
    for (const auto &entry : grouped) {
        result.known_pivots.insert(entry.first);
        states[entry.first] = MonomialState::Pivot;
        bool first = true;
        for (size_t index : entry.second) {
            GPoly row = multiplied_row(basis[index], entry.first);
            result.rows.push_back(row);
            if (first) {
                first = false;
            }
            for (const auto &monom : collect_row_tail(row)) {
                states.emplace(monom, MonomialState::Todo);
            }
        }
    }

    while (true) {
        bool has_todo = false;
        for (const auto &entry : states) {
            if (entry.second == MonomialState::Todo) {
                has_todo = true;
                break;
            }
        }
        if (not has_todo) {
            return result;
        }

        vec_int monom = next_todo(states, order);
        bool found_reducer = false;
        for (const auto &poly : basis) {
            if (monomial_divides(poly.LM(), monom)) {
                GPoly row = multiplied_row(poly, monom);
                result.rows.push_back(row);
                states[monom] = MonomialState::Pivot;
                result.known_pivots.insert(monom);
                for (const auto &tail_monom : collect_row_tail(row)) {
                    states.emplace(tail_monom, MonomialState::Todo);
                }
                found_reducer = true;
                break;
            }
        }
        if (not found_reducer) {
            states[monom] = MonomialState::Done;
        }
    }
}

std::vector<vec_int> sorted_columns(const std::vector<GPoly> &rows,
                                    MonomialOrder order)
{
    std::set<vec_int> unique;
    for (const auto &row : rows) {
        for (const auto &term : row.terms) {
            unique.insert(term.first);
        }
    }

    std::vector<vec_int> columns(unique.begin(), unique.end());
    std::sort(columns.begin(), columns.end(),
              [order](const vec_int &a, const vec_int &b) {
                  return monomial_compare(a, b, order) > 0;
              });
    return columns;
}

DenseMatrix build_matrix(const std::vector<GPoly> &rows,
                         const std::vector<vec_int> &columns)
{
    vec_basic entries(rows.size() * columns.size(), zero);
    DenseMatrix matrix(numeric_cast<unsigned int>(rows.size()),
                       numeric_cast<unsigned int>(columns.size()), entries);

    std::map<vec_int, size_t> column_index;
    for (size_t i = 0; i < columns.size(); ++i) {
        column_index.emplace(columns[i], i);
    }

    for (size_t row = 0; row < rows.size(); ++row) {
        for (const auto &term : rows[row].terms) {
            const size_t col = column_index.at(term.first);
            matrix.set(numeric_cast<unsigned int>(row), numeric_cast<unsigned int>(col),
                       detail::normalize_coeff(term.second).get_basic());
        }
    }
    return matrix;
}

GPoly row_to_gpoly(const DenseMatrix &matrix, unsigned int row,
                   const std::vector<vec_int> &columns, unsigned int nvars,
                   MonomialOrder order)
{
    GPoly poly(nvars, order);
    for (unsigned int col = 0; col < matrix.ncols(); ++col) {
        Expression coeff(matrix.get(row, col));
        if (not is_zero_coeff(coeff)) {
            detail::insert_term(poly, columns[col], detail::normalize_coeff(coeff));
        }
    }
    return poly;
}

std::vector<GPoly> new_basis_rows(const DenseMatrix &matrix,
                                  const std::vector<vec_int> &columns,
                                  const std::set<vec_int> &known_pivots,
                                  const std::vector<GPoly> &basis)
{
    std::set<vec_int> existing_lms;
    for (const auto &poly : basis) {
        existing_lms.insert(poly.LM());
    }

    std::vector<GPoly> result;
    std::set<vec_int> new_lms;
    const unsigned int nvars = basis.front().nvars;
    const MonomialOrder order = basis.front().order;
    for (unsigned int row = 0; row < matrix.nrows(); ++row) {
        GPoly poly = row_to_gpoly(matrix, row, columns, nvars, order);
        if (poly.is_zero()) {
            continue;
        }
        poly = gpoly_monic(poly);
        if (known_pivots.count(poly.LM()) != 0) {
            continue;
        }
        if (existing_lms.count(poly.LM()) != 0
            or new_lms.insert(poly.LM()).second == false) {
            continue;
        }
        result.push_back(poly);
    }

    std::sort(result.begin(), result.end(),
              [order](const GPoly &a, const GPoly &b) {
                  return monomial_compare(a.LM(), b.LM(), order) > 0;
              });
    return result;
}

void append_new_pairs(std::vector<F4Pair> &pairs, size_t old_size, size_t new_size,
                      const std::vector<GPoly> &basis)
{
    for (size_t i = old_size; i < old_size + new_size; ++i) {
        for (size_t j = 0; j < i; ++j) {
            vec_int lcm(basis[i].nvars, 0);
            monomial_lcm(basis[i].LM(), basis[j].LM(), lcm);
            pairs.push_back(F4Pair{
                i, j, lcm, static_cast<unsigned int>(monomial_total_degree(lcm))});
        }
    }
}

} // namespace

std::vector<GPoly> groebner_f4(std::vector<GPoly> polys)
{
    std::vector<GPoly> basis = detail::reduce_generators(std::move(polys));
    if (basis.empty()) {
        return basis;
    }

    std::vector<F4Pair> pairs = initial_pairs(basis);
    while (not pairs.empty()) {
        auto selected_and_remaining
            = select_minimal_degree_pairs(std::move(pairs), basis.front().order);
        std::vector<F4Pair> selected = std::move(selected_and_remaining.first);
        pairs = std::move(selected_and_remaining.second);

        PreprocessResult prep = symbolic_preprocess(selected, basis);
        if (prep.rows.empty()) {
            continue;
        }

        const std::vector<vec_int> columns
            = sorted_columns(prep.rows, basis.front().order);
        DenseMatrix matrix = build_matrix(prep.rows, columns);
        DenseMatrix reduced(matrix.nrows(), matrix.ncols());
        vec_uint pivot_cols;
        // RREF pivoting uses SymEngine's structural zero checks. That is exact for
        // numeric/rational entries and only best-effort for symbolic fractions.
        reduced_row_echelon_form(matrix, reduced, pivot_cols);

        std::vector<GPoly> new_rows
            = new_basis_rows(reduced, columns, prep.known_pivots, basis);
        if (new_rows.empty()) {
            continue;
        }

        const size_t old_size = basis.size();
        for (auto &poly : new_rows) {
            basis.push_back(std::move(poly));
        }
        append_new_pairs(pairs, old_size, new_rows.size(), basis);
    }

    return detail::red_groebner(basis);
}

} // namespace SymEngine
