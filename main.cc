#include <bitset>
#include <cassert>
#include <iostream>

#include "differences.h"
#include "timer.h"

typedef unsigned __int128 uint128_t;

constexpr uint32_t BLOOM_FILTER_M = 1 << 28;
std::bitset<BLOOM_FILTER_M> bloom_filter[2];

constexpr uint32_t HASH_MAP_M = 1 << 25;
uint32_t hash_map[HASH_MAP_M + 16];

uint32_t inline hash_1(uint64_t x) { return x >> 32; }
uint32_t inline hash_2(uint64_t x) { return x >> 8; }

void insert(uint64_t x) {
  assert(x != 0);
  bloom_filter[0][hash_1(x) >> 4] = true;
  bloom_filter[1][hash_2(x) >> 4] = true;

  auto h = hash_1(x ^ (x << 24)) & (HASH_MAP_M - 1);
  while (hash_map[h] != 0) {
    h++;
  }
  // Ensure a sentinel 0 value is left the end of the hash map. This assert will
  // fail the first time we try to write the last slot.
  assert(h + 1 < HASH_MAP_M + 16);
  hash_map[h] = x;
}

bool contains(uint64_t x) {
  if (bloom_filter[0][hash_1(x) >> 4] && bloom_filter[1][hash_2(x) >> 4]) {
    auto x_int32 = uint32_t(x);
    auto h = hash_1(x ^ (x << 24)) & (HASH_MAP_M - 1);
    while (hash_map[h] != 0 && hash_map[h] != x_int32) {
      h++;
    }
    return hash_map[h] != 0;
  }
  return false;
}

constexpr uint32_t MAX_D = 500000;
uint64_t quartic_powers_uint64[MAX_D + 1];

int main() {
  Timer t = Timer();
  std::cout << "Searching up to D = " << MAX_D << std::endl << std::endl;

  // Compute differences
  auto differences = compute_differences(MAX_D);
  t.log_task("Compute differences");

  // Populate hash map
  for (auto [diff, _c, _d] : differences) {
    insert(diff / 625);
  }
  t.log_task("Populate Bloom filter and hash map");

  // Precompute quartic powers
  for (uint128_t i = 0; i <= MAX_D; i++) {
    quartic_powers_uint64[i] = i * i * i * i;
  }

  // Check all pairwise sums, assuming wlog that a = 5 * i and b = 5 * j
  for (uint32_t i = 1; i <= MAX_D / 5; i++) {
    for (uint32_t j = 1; j <= i; j++) {
      // MODULO 4
      // If a and b are both odd, then a^4 + b^4 = 2 mod 4. However, no
      // difference of fourth powers is 2 mod 4.
      if (i % 2 == 1 && j % 2 == 1) {
        continue;
      }
      auto sum_uint64 = quartic_powers_uint64[i] + quartic_powers_uint64[j];
      if (contains(sum_uint64)) [[unlikely]] {
        uint32_t a = 5 * i;
        uint32_t b = 5 * j;
        auto sum = uint128_t(1) * a * a * a * a + uint128_t(1) * b * b * b * b;
        for (auto [diff, c, d] : differences) {
          if (diff == sum) {
            t.log_task("Check pairwise sums");
            std::cout << std::endl;
            std::cout << "Solution found: " << a << "^4 + " << b << "^4 + " << c
                      << "^4 = " << d << "^4" << std::endl;
            return 0;
          }
        }
      }
    }
  }
}