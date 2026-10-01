#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../board.hpp"

namespace minclicks::detail {

using Clock = std::chrono::steady_clock;

inline double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

inline void unique_sort(std::vector<int> &v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

inline bool subset(const std::vector<int> &a, const std::vector<int> &b) {
    return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

struct DSU {
    std::vector<int> p;

    explicit DSU(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }

    int root(int a) {
        while (p[a] != a) {
            p[a] = p[p[a]];
            a = p[a];
        }
        return a;
    }

    void join(int a, int b) { p[root(a)] = root(b); }
};

inline void flip(std::string &s, int i) {
    s[i / 8] ^= char(1u << (i % 8));
}

}  // namespace minclicks::detail