#include "differences.h"

#include <cassert>
#include <iostream>
#include <vector>

constexpr uint32_t Q = 625;
constexpr uint32_t M = 24 * Q;

std::vector<uint32_t> compute_quartic_powers_mod_m(uint32_t m) {
  std::vector<uint32_t> quartic_powers(m);
  for (uint128_t i = 0; i < m; i++) {
    quartic_powers[i] = i * i * i * i % m;
  }
  return quartic_powers;
}

std::vector<uint32_t> sums_of_quartic_residues_mod_m(uint32_t m) {
  auto quartic_powers_mod_m = compute_quartic_powers_mod_m(m);
  std::vector<uint32_t> is_sum_of_quartic_residue_mod_m(m, 0);
  for (uint32_t i = 0; i < m; i++) {
    for (uint32_t j = 0; j < m; j++) {
      auto s = (quartic_powers_mod_m[i] + quartic_powers_mod_m[j]) % m;
      is_sum_of_quartic_residue_mod_m[s] = 1;
    }
  }
  return is_sum_of_quartic_residue_mod_m;
}

auto find_good_pairs_mod_M() {
  auto pow4_mod_Q = compute_quartic_powers_mod_m(Q);
  std::vector<std::pair<uint32_t, uint32_t>> good_pairs;
  // Let i = d % M and j = c % M.
  for (uint32_t i = 0; i < M; i++) {
    // MODULO 8
    // Filter out even numbers, since if d is even in a^4 + b^4 + c^4 = d^4,
    // then two of a, b, c must be odd, in which case LHS is 2 != 0 mod 4.
    // In fact, Ward (Duke Math. J., 1948) shows a fortiori that d = 1 mod 8.
    if (i % 8 != 1) {
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
      // MODULO 9
      // No sum of two fourth powers is 0 mod 9 unless both are 0 mod 3.
      if (i % 3 == 0 && j % 3 == 0) {
        continue;
      }
      // MODULO 8
      // Ward (Duke Math. J., 1948) shows that (a % 8, b % 8, c % 8) is either a
      // permutation of (0, 0, 1) or (0, 0, 7). Thus, j % 8 is in {0, 1, 7}.
      if (j % 8 != 0 && j % 8 != 1 && j % 8 != 7) {
        continue;
      }

      good_pairs.push_back({i, j});
    }
  }

  return good_pairs;
}

std::vector<CandidateDifference> compute_differences(uint32_t max_d) {
  assert(Q == 625);           // Q is 625
  assert(M % (24 * Q) == 0);  // 24 * Q divides M

  auto good_pairs = find_good_pairs_mod_M();
  std::cout << "Found " << good_pairs.size() << " good pairs ("
            << 100 * float(good_pairs.size()) / M / M << "%)" << std::endl;

  // Precompute (sums of) quartic residues mod 2^8, 3^6, 7^3, 11^2, 13^2, 29^2
  auto pow4_mod_256 = compute_quartic_powers_mod_m(256);
  auto pow4_mod_729 = compute_quartic_powers_mod_m(729);
  auto pow4_mod_343 = compute_quartic_powers_mod_m(343);
  auto pow4_mod_121 = compute_quartic_powers_mod_m(121);
  auto pow4_mod_169 = compute_quartic_powers_mod_m(169);
  auto pow4_mod_841 = compute_quartic_powers_mod_m(841);
  auto is_sum_of_quartic_residues_mod_256 = sums_of_quartic_residues_mod_m(256);
  auto is_sum_of_quartic_residues_mod_729 = sums_of_quartic_residues_mod_m(729);
  auto is_sum_of_quartic_residues_mod_343 = sums_of_quartic_residues_mod_m(343);
  auto is_sum_of_quartic_residues_mod_121 = sums_of_quartic_residues_mod_m(121);
  auto is_sum_of_quartic_residues_mod_169 = sums_of_quartic_residues_mod_m(169);
  auto is_sum_of_quartic_residues_mod_841 = sums_of_quartic_residues_mod_m(841);

  // Compute quartic powers
  auto quartic_powers = std::vector<uint128_t>(max_d + 1);
  for (uint128_t i = 0; i <= max_d; i++) {
    quartic_powers[i] = i * i * i * i;
  }

  // MODULO P^2 for P odd, P != 1 mod 8
  // Suppose v_P(d^4 - c^4) != 0 mod 4. Then, let a = P^m a' and b = P^m b' be
  // such that, wlog, a' != 0 mod P. Then, v_P(a'^4 + b'^4) = v_P(d^4 - c^4) -
  // 4m > 0, so a'^4 + b'^4 = 0 mod P. Then, b'/a' is a primitive 8th root of
  // unity mod P, which can only happen if P = 1 mod 8.
  //
  // By the lifting the exponent lemma, if P does not divide either d or c, then
  // v_P(d - c) > 0 implies v_P(d^4 - c^4) = v_P(d - c), and v_P(d + c) > 0
  // implies v_P(d^4 - (-c)^4) = v_P(d + c). If P divides both d and c, then
  // P divides a^4 + b^4, so P = 1 mod 8. So if P divides d - c (or d + c) not a
  // multiple of 4 times, then P = 1 mod 8.
  std::vector<bool> is_prime(2 * max_d + 1, true);
  std::vector<bool> is_bad(2 * max_d + 1, false);
  is_prime[0] = false;
  is_prime[1] = false;
  for (uint32_t i = 2; i <= 2 * max_d; i++) {
    if (!is_prime[i]) {
      continue;
    }
    for (uint32_t j = i + i; j <= 2 * max_d; j += i) {
      is_prime[j] = false;
    }
    if (i % 8 != 1 && i != 2) {
      for (uint32_t j = i; j <= 2 * max_d; j += i) {
        if ((j / i) % i != 0 || (j / i / i) % i != 0 ||
            (j / i / i / i) % i != 0) {
          is_bad[j] = true;
        }
      }
    }
  }

  // Compute differences
  std::vector<CandidateDifference> differences;
  for (auto [k, l] : good_pairs) {
    for (int i = 0; M * i <= max_d; i++) {
      for (int j = 0; j <= i; j++) {
        uint32_t d = M * i + k;
        uint32_t c = M * j + l;
        // Ensure 0 < c < d <= max_d
        if (c == 0 || c >= d || d > max_d) {
          continue;
        }

        // MODULO 4096
        // Morgan (Duke Math. J., 1948) shows that, if c is odd, d^4 - c^4 = 0
        // mod 4096.
        auto pow2_d_mod_4096 = (d % 4096) * (d % 4096) % 4096;
        auto pow2_c_mod_4096 = (c % 4096) * (c % 4096) % 4096;
        auto pow4_d_mod_4096 = pow2_d_mod_4096 * pow2_d_mod_4096 % 4096;
        auto pow4_c_mod_4096 = pow2_c_mod_4096 * pow2_c_mod_4096 % 4096;
        if (c % 2 == 1 && pow4_c_mod_4096 != pow4_d_mod_4096) {
          continue;
        }
        // v_P(d \pm c)
        // If v_P(d - c) != 0 mod 4 for an odd prime P != 1 mod 8 then skip.
        // If v_P(d + c) != 0 mod 4 for an odd prime P != 1 mod 8 then skip.
        if (is_bad[d - c] || is_bad[d + c]) {
          continue;
        }
        // Ensure diff is a sum of two fourth powers mod 2^8, 3^6, 7^3, 11^2,
        // 13^2, 29^2.
        if (!is_sum_of_quartic_residues_mod_256
                [(pow4_mod_256[d % 256] - pow4_mod_256[c % 256] + 256) % 256] ||
            !is_sum_of_quartic_residues_mod_729
                [(pow4_mod_729[d % 729] - pow4_mod_729[c % 729] + 729) % 729] ||
            !is_sum_of_quartic_residues_mod_343
                [(pow4_mod_343[d % 343] - pow4_mod_343[c % 343] + 343) % 343] ||
            !is_sum_of_quartic_residues_mod_121
                [(pow4_mod_121[d % 121] - pow4_mod_121[c % 121] + 121) % 121] ||
            !is_sum_of_quartic_residues_mod_169
                [(pow4_mod_169[d % 169] - pow4_mod_169[c % 169] + 169) % 169] ||
            !is_sum_of_quartic_residues_mod_841
                [(pow4_mod_841[d % 841] - pow4_mod_841[c % 841] + 841) % 841]) {
          continue;
        }

        auto diff = quartic_powers[d] - quartic_powers[c];
        // Check v_5(d^4 - c^4) = 0 mod 4 (noting that 625 | d^4 - c^4)
        if (diff % 3125 == 0 && diff % 390625 != 0) {
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