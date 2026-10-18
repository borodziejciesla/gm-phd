#ifndef GM_PHD_INCLUDE_OSPA_HPP_
#define GM_PHD_INCLUDE_OSPA_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace mot {
using OspaPosition = std::array<double, 2u>;

inline double OspaDistance(const std::vector<OspaPosition>& estimates,
                           const std::vector<OspaPosition>& truth, const double cutoff,
                           const double order) {
  if (!std::isfinite(cutoff) || cutoff <= 0.0 || !std::isfinite(order) || order <= 0.0) {
    throw std::invalid_argument("OSPA cutoff and order must be finite and positive");
  }
  for (const auto& point : estimates) {
    if (!std::isfinite(point[0]) || !std::isfinite(point[1])) {
      throw std::invalid_argument("OSPA positions must be finite");
    }
  }
  for (const auto& point : truth) {
    if (!std::isfinite(point[0]) || !std::isfinite(point[1])) {
      throw std::invalid_argument("OSPA positions must be finite");
    }
  }

  const auto max_size = std::max(estimates.size(), truth.size());
  if (max_size == 0u) {
    return 0.0;
  }
  if (std::min(estimates.size(), truth.size()) == 0u) {
    return cutoff;
  }

  const auto* rows = &estimates;
  const auto* columns = &truth;
  if (rows->size() > columns->size()) {
    std::swap(rows, columns);
  }
  const auto row_count = rows->size();
  const auto column_count = columns->size();
  const auto cutoff_cost = std::pow(cutoff, order);
  const auto assignment_cost = [rows, columns, cutoff, order](const size_t row,
                                                              const size_t column) {
    const auto dx = (*rows)[row][0] - (*columns)[column][0];
    const auto dy = (*rows)[row][1] - (*columns)[column][1];
    return std::pow(std::min(cutoff, std::hypot(dx, dy)), order);
  };

  std::vector<double> row_potential(row_count + 1u, 0.0);
  std::vector<double> column_potential(column_count + 1u, 0.0);
  std::vector<size_t> matched_row(column_count + 1u, 0u);
  std::vector<size_t> previous_column(column_count + 1u, 0u);

  for (size_t row = 1u; row <= row_count; ++row) {
    matched_row[0] = row;
    size_t current_column = 0u;
    std::vector<double> min_reduced_cost(column_count + 1u,
                                         std::numeric_limits<double>::infinity());
    std::vector<bool> used(column_count + 1u, false);

    do {
      used[current_column] = true;
      const auto current_row = matched_row[current_column];
      auto delta = std::numeric_limits<double>::infinity();
      size_t next_column = 0u;

      for (size_t column = 1u; column <= column_count; ++column) {
        if (used[column]) {
          continue;
        }
        const auto reduced_cost =
            assignment_cost(current_row - 1u, column - 1u) - row_potential[current_row] -
            column_potential[column];
        if (reduced_cost < min_reduced_cost[column]) {
          min_reduced_cost[column] = reduced_cost;
          previous_column[column] = current_column;
        }
        if (min_reduced_cost[column] < delta) {
          delta = min_reduced_cost[column];
          next_column = column;
        }
      }

      for (size_t column = 0u; column <= column_count; ++column) {
        if (used[column]) {
          row_potential[matched_row[column]] += delta;
          column_potential[column] -= delta;
        } else {
          min_reduced_cost[column] -= delta;
        }
      }
      current_column = next_column;
    } while (matched_row[current_column] != 0u);

    do {
      const auto next_column = previous_column[current_column];
      matched_row[current_column] = matched_row[next_column];
      current_column = next_column;
    } while (current_column != 0u);
  }

  auto total_cost = cutoff_cost * (max_size - row_count);
  for (size_t column = 1u; column <= column_count; ++column) {
    if (matched_row[column] != 0u) {
      total_cost += assignment_cost(matched_row[column] - 1u, column - 1u);
    }
  }

  return std::pow(total_cost / static_cast<double>(max_size), 1.0 / order);
}
}  // namespace mot

#endif  // GM_PHD_INCLUDE_OSPA_HPP_
