#include "differences.h"

#include <cassert>
#include <iostream>
#include <vector>

constexpr uint32_t P = 54;
constexpr uint32_t Q = 625;
constexpr uint32_t M = P * Q;

std::vector<uint128_t> compute_quartic_powers(uint32_t n, bool mod_n = true) {
  std::vector<uint128_t> quartic_powers(n);
  for (uint128_t i = 0; i < n; i++) {
    quartic_powers[i] = i * i * i * i;
    if (mod_n) {
      quartic_powers[i] %= n;
    }
  }
  return quartic_powers;
}

std::vector<bool> sums_of_quartic_residues_mod_m(uint32_t m) {
  auto quartic_powers_mod_m = compute_quartic_powers(m);
  std::vector<bool> is_sum_of_quartic_residue_mod_m(m, false);
  for (uint32_t i = 0; i < m; i++) {
    for (uint32_t j = 0; j < m; j++) {
      auto s = (quartic_powers_mod_m[i] + quartic_powers_mod_m[j]) % m;
      is_sum_of_quartic_residue_mod_m[s] = true;
    }
  }
  return is_sum_of_quartic_residue_mod_m;
}

auto find_good_pairs_mod_M() {
  auto pow4_mod_Q = compute_quartic_powers(Q);
  auto pow4_mod_P = compute_quartic_powers(P);
  auto is_sum_of_pow4_mod_P = sums_of_quartic_residues_mod_m(P);

  std::vector<std::pair<uint32_t, uint32_t>> good_pairs;
  // Let i = d % M and j = c % M.
  for (uint32_t i = 0; i < M; i++) {
    // MODULO 4
    // Filter out even numbers, since if d is even in a^4 + b^4 + c^4 = d^4,
    // then two of a, b, c must be odd, in which case LHS is 2 != 0 mod 4.
    if (i % 2 == 0) {
      continue;
    }
    // MODULO 5
    // Since x^4 % 5 is either 0 or 1, in any minimal solution, one has d % 5 !=
    // 0 and two of a, b, c are divisible by 5. We assume without loss of
    // generality that a and b are divisble by 5.
    if (i % 5 == 0) {
      continue;
    }

    for (uint32_t j = 0; j < M; j++) {
      // MODULO 5
      // If a % 5 == b % 5 == 0, then (d^4 - c^4) % 625 == 0.
      if (pow4_mod_Q[i % Q] != pow4_mod_Q[j % Q]) {
        continue;
      }
      // MODULO P
      // Ensure (d^4 - c^4) % P is a sum of two fourth powers modulo P.
      auto diff = (pow4_mod_P[i % P] - pow4_mod_P[j % P] + P) % P;
      if (!is_sum_of_pow4_mod_P[diff]) {
        continue;
      }

      good_pairs.push_back({i, j});
    }
  }

  return good_pairs;
}

std::vector<CandidateDifference> compute_differences(uint32_t max_d) {
  assert(Q == 625);                  // Q is 625
  assert(P % 2 == 0 && P % 5 != 0);  // P is even and coprime to Q

  auto good_pairs = find_good_pairs_mod_M();
  std::cout << "Found " << good_pairs.size() << " good pairs ("
            << 100 * float(good_pairs.size()) / M / M << "%)" << std::endl;

  // Precompute sums of quartic residues mod 2^8, 3^6, 13^2, 29^2,
  auto is_sum_of_quartic_residues_mod_256 = sums_of_quartic_residues_mod_m(256);
  auto is_sum_of_quartic_residues_mod_729 = sums_of_quartic_residues_mod_m(729);
  auto is_sum_of_quartic_residues_mod_169 = sums_of_quartic_residues_mod_m(169);
  auto is_sum_of_quartic_residues_mod_841 = sums_of_quartic_residues_mod_m(841);

  // Compute quartic powers
  auto quartic_powers = compute_quartic_powers(max_d + 1, false);

  // Compute differences
  std::vector<CandidateDifference> differences;
  for (int i = 0; M * i <= max_d; i++) {
    for (int j = 0; j <= i; j++) {
      for (auto [k, l] : good_pairs) {
        uint32_t d = M * i + k;
        uint32_t c = M * j + l;
        // Ensure 0 < c < d <= max_d
        if (c == 0 || c >= d || d > max_d) {
          continue;
        }

        auto diff = quartic_powers[d] - quartic_powers[c];

        // Ensure diff is a sum of two fourth powers mod 2^8, 3^6, 13^2, 29^2.
        if (!is_sum_of_quartic_residues_mod_256[diff % 256] ||
            !is_sum_of_quartic_residues_mod_729[diff % 729] ||
            !is_sum_of_quartic_residues_mod_169[diff % 169] ||
            !is_sum_of_quartic_residues_mod_841[diff % 841]) {
          continue;
        }

        differences.push_back({diff, c, d});
      }
    }
  }

  std::cout << "Found " << differences.size() << " candidate differences ("
            << 100 * float(differences.size()) / max_d / max_d << "%)"
            << std::endl;

  return differences;
}