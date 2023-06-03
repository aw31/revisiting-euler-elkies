#pragma once

#include <vector>

typedef unsigned __int128 uint128_t;

struct CandidateDifference {
  uint128_t diff;
  uint32_t c;
  uint32_t d;
};

std::vector<CandidateDifference> compute_differences(uint32_t max_d);