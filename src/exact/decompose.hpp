#pragma once

#include "model.hpp"

namespace minclicks::detail {

inline std::vector<std::vector<int>> independentParts(const Model &m) {
    int n = int(m.cells.size()), total = int(m.edges.size());
    auto graph = m.edges;

    for (const auto &factor : m.factors)
        for (std::size_t j = 1; j < factor.vars.size(); ++j) {
            int a = factor.vars[0], b = factor.vars[j];
            graph[a].push_back(b);
            graph[b].push_back(a);
        }

    for (auto &adjacent : graph) detail::unique_sort(adjacent);

    std::vector<int> discovered(total, -1), low(total);
    std::vector<bool> cut(total);
    int time = 0;

    auto dfs = [&](auto &&self, int v, int parent) -> void {
        discovered[v] = low[v] = time++;
        int children = 0;
        for (int u : graph[v])
            if (u != parent) {
                if (discovered[u] < 0) {
                    ++children;
                    self(self, u, v);
                    low[v] = std::min(low[v], low[u]);
                    if (v >= n && parent >= 0 && low[u] >= discovered[v]) cut[v] = true;
                } else low[v] = std::min(low[v], discovered[u]);
            }
        if (v >= n && parent < 0 && children > 1) cut[v] = true;
    };

    for (int v = 0; v < total; ++v)
        if (discovered[v] < 0) dfs(dfs, v, -1);

    std::vector<bool> seen(total);
    std::vector<std::vector<int>> parts;

    for (int v = 0; v < n; ++v)
        if (!seen[v]) {
            std::vector<int> part, stack{v};
            seen[v] = true;
            while (!stack.empty()) {
                int u = stack.back();
                stack.pop_back();
                if (u < n) part.push_back(u);
                for (int other : graph[u])
                    if (!cut[other] && !seen[other]) {
                        seen[other] = true;
                        stack.push_back(other);
                    }
            }
            std::sort(part.begin(), part.end());
            parts.push_back(std::move(part));
        }
        
    return parts;
}

}  // namespace minclicks::detail
