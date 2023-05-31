#include <bitset>
#include <iostream>

#include "differences.h"
#include "timer.h"

typedef unsigned __int128 uint128_t;

constexpr uint32_t BLOOM_FILTER_M = 1 << 29;
std::bitset<BLOOM_FILTER_M> bloom_filter;

constexpr uint32_t HASH_MAP_M = 1 << 25;
uint64_t hash_map[HASH_MAP_M];

std::pair<uint32_t, uint32_t> inline hash(uint64_t x) {
  return {x ^ (x >> 32), x >> 32};
}

void insert(uint64_t x) {
  assert(x != 0);
  auto [h1, h2] = hash(x);
  bloom_filter.set(h1 % BLOOM_FILTER_M);
  bloom_filter.set(h2 % BLOOM_FILTER_M);

  auto h = h1 % HASH_MAP_M;
  while (hash_map[h] != 0) {
    h = (h + 1) % HASH_MAP_M;
  }
  hash_map[h] = x;
}

bool contains(uint64_t x) {
  auto [h1, h2] = hash(x);
  if (bloom_filter.test(h1 % BLOOM_FILTER_M) &&
      bloom_filter.test(h2 % BLOOM_FILTER_M)) {
    auto h = h1 % HASH_MAP_M;
    while (hash_map[h] != 0 && hash_map[h] != x) {
      h = (h + 1) % HASH_MAP_M;
    }
    return hash_map[h] != 0;
  }
  return false;
}

constexpr uint32_t MAX_D = 500000;

int main() {
  Timer t = Timer();
  std::cout << "Searching up to D = " << MAX_D << std::endl << std::endl;

  // Compute differences
  auto differences = compute_differences(MAX_D);
  t.log_task("Compute differences");

  // Populate hash map
  for (auto [diff, _c, _d] : differences) {
    insert(diff);
  }
  t.log_task("Populate Bloom filter and hash map");

  // Precompute quartic powers
  std::vector<uint128_t> quartic_powers(MAX_D + 1);
  std::vector<uint64_t> quartic_powers_uint64(MAX_D + 1);
  for (uint128_t i = 0; i <= MAX_D; i++) {
    quartic_powers[i] = i * i * i * i;
    quartic_powers_uint64[i] = quartic_powers[i];
  }

  // Check all pairwise sums
  for (uint32_t a = 5; a <= MAX_D; a += 5) {
    for (uint32_t b = 5; b <= a; b += 5) {
      // MODULO 4
      // If a and b are both odd, then a^4 + b^4 = 2 mod 4. However, no
      // difference of fourth powers is 2 mod 4.
      if (a % 2 == 1 && b % 2 == 1) {
        continue;
      }
      auto sum_uint64 = quartic_powers_uint64[a] + quartic_powers_uint64[b];
      if (contains(sum_uint64)) {
        auto sum = quartic_powers[a] + quartic_powers[b];
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