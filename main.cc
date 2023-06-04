#include <omp.h>

#include <bitset>
#include <cassert>
#include <iostream>

#include "differences.h"
#include "timer.h"

typedef unsigned __int128 uint128_t;

constexpr uint32_t PREFIX_LUT_K = 32;
constexpr uint64_t PREFIX_LUT_M = 1ULL << PREFIX_LUT_K;
std::bitset<PREFIX_LUT_M> prefix_lut;

constexpr uint32_t HASH_MAP_K = 27;
constexpr uint32_t HASH_MAP_M = 1 << HASH_MAP_K;
uint32_t hash_map[HASH_MAP_M + 16];

uint32_t inline hash_1(uint64_t x) { return x >> 32; }
uint32_t inline hash_2(uint64_t x) { return (x >> 24) ^ (x >> 16); }

void insert(uint64_t x) {
  assert(x != 0);
  prefix_lut[hash_1(x) >> (32 - PREFIX_LUT_K)] = true;

  auto i = hash_2(x) >> (32 - HASH_MAP_K);
  while (hash_map[i] != 0) {
    i++;
  }
  // Ensure a sentinel 0 value is left the end of the hash map. This assert will
  // fail the first time we try to write the last slot.
  assert(i + 1 < HASH_MAP_M + 16);
  hash_map[i] = x;
}

bool contains(uint64_t x) {
  if (prefix_lut[hash_1(x) >> (32 - PREFIX_LUT_K)]) {
    auto x_int32 = uint32_t(x);
    auto i = hash_2(x) >> (32 - HASH_MAP_K);
    while (hash_map[i] != 0 && hash_map[i] != x_int32) {
      i++;
    }
    return hash_map[i] != 0;
  }
  return false;
}

constexpr uint32_t MAX_D = 10000000;
uint64_t pow4_uint64[MAX_D + 1];

void verify_ab(std::vector<CandidateDifference> candidate_differences,
               uint32_t a, uint32_t b) {
  auto sum = uint128_t(1) * a * a * a * a + uint128_t(1) * b * b * b * b;
  for (auto [diff, c, d] : candidate_differences) {
    if (diff == sum) {
      std::cout << std::endl;
      std::cout << "Solution found: " << a << "^4 + " << b << "^4 + " << c
                << "^4 = " << d << "^4" << std::endl;
    }
  }
}

int main() {
  Timer t = Timer();
  std::cout << "Searching up to D = " << MAX_D << " with "
            << omp_get_max_threads() << " threads" << std::endl
            << std::endl;

  // Compute differences
  auto differences = compute_differences(MAX_D);
  t.log_task("Compute differences");

  // Populate hash map
  for (auto [diff, _c, _d] : differences) {
    insert(diff / 625);
  }
  t.log_task("Populate filter and hash map");

  // Precompute quartic powers
  for (uint128_t i = 0; i <= MAX_D; i++) {
    pow4_uint64[i] = i * i * i * i;
  }

  // Check all pairwise sums, assuming wlog that a = 5 * i and b = 5 * j
  // Ward (Duke Math. J., 1948) shows that (a % 8, b % 8, c % 8) is a
  // permutation of either (0, 0, 1) or (0, 0, 7). Therefore, (i % 8, j % 8) =
  // (0, 0), (5, 0), (0, 5), (3, 0), or (0, 3).
#pragma omp parallel for
  for (uint32_t i = 8; i <= MAX_D / 5; i += 8) {
    for (uint32_t j = 8; j <= i; j += 8) {
      if (contains(pow4_uint64[i] + pow4_uint64[j])) [[unlikely]] {
        verify_ab(differences, 5 * i, 5 * j);
      }
      if (contains(pow4_uint64[i + 5] + pow4_uint64[j])) [[unlikely]] {
        verify_ab(differences, 5 * (i + 5), 5 * j);
      }
      if (contains(pow4_uint64[i] + pow4_uint64[j + 5])) [[unlikely]] {
        verify_ab(differences, 5 * i, 5 * (j + 5));
      }
      if (contains(pow4_uint64[i + 3] + pow4_uint64[j])) [[unlikely]] {
        verify_ab(differences, 5 * (i + 3), 5 * j);
      }
      if (contains(pow4_uint64[i] + pow4_uint64[j + 3])) [[unlikely]] {
        verify_ab(differences, 5 * i, 5 * (j + 3));
      }
    }
  }
  t.log_task("Check pairwise sums");
}