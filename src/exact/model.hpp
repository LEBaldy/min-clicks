#pragma once

#include "common.hpp"

namespace minclicks::detail {

struct Factor {
    std::vector<int> vars;
    bool mine;
    int weight = 1;
};

struct Model {
    const Board &b;
    std::vector<int> cells;
    std::vector<std::vector<int>> edges;
    std::vector<Factor> factors;
    int base = 0, removed = 0, isolatedZeros = 0;

    explicit Model(const Board &board) : b(board) {
        std::vector<std::vector<int>> safe(b.w * b.h), mines(b.w * b.h);
        std::vector<int> candidates;

        for (int p = 0; p < b.w * b.h; ++p)
            if (b.number[p] > 0) {
                safe[p].push_back(p);
                for (int q : b.adj[p]) (b.mine[q] ? mines[p] : safe[p]).push_back(q);
                unique_sort(safe[p]);
                unique_sort(mines[p]);
                candidates.push_back(p);
            }

        std::vector<std::vector<int>> flood = b.zeros;

        for (std::size_t z = 0; z < flood.size(); ++z) {
            flood[z].insert(flood[z].end(), b.borders[z].begin(), b.borders[z].end());
            unique_sort(flood[z]);
        }

        std::vector<bool> redundant(b.w * b.h);

        for (int p : candidates) {
            redundant[p] = safe[p].size() <= 2;

            for (const auto &region : flood)
                if (subset(safe[p], region)) redundant[p] = true;

            for (int q : b.adj[p])
                if (b.island[q] >= 0) {
                    const auto &region = flood[b.island[q]];
                    safe[p].insert(safe[p].end(), region.begin(), region.end());
                }

            unique_sort(safe[p]);
        }

        for (int a : candidates) {
            bool dominated = redundant[a];

            if (!dominated)
                for (int other : candidates)
                    if (a != other && subset(safe[a], safe[other]) &&
                        subset(mines[other], mines[a]) &&
                        (safe[a] != safe[other] || mines[a] != mines[other] ||
                         other < a)) {
                        dominated = true;
                        break;
                    }
            if (!dominated) cells.push_back(a);
            else ++removed;
        }

        int n = (int)cells.size(), z = (int)b.zeros.size();
        edges.resize(n + z);
        std::vector<int> id(b.w * b.h, -1);

        for (int i = 0; i < n; ++i) id[cells[i]] = i;

        for (int i = 0; i < n; ++i) {
            for (int q : b.adj[cells[i]]) {
                if (id[q] >= 0) edges[i].push_back(id[q]);
                if (b.island[q] >= 0) edges[i].push_back(n + b.island[q]);
            }
            unique_sort(edges[i]);
            for (int j : edges[i])
                if (j >= n) edges[j].push_back(i);
        }

        {
            auto original = edges;
            for (int i = 0; i < n; ++i) {
                edges[i].clear();
                for (int j : original[i]) {
                    bool redundant = false;
                    if (j < n)
                        for (int z : original[i])
                            if (z >= n && binary_search(original[j].begin(),
                                                        original[j].end(), z)) {
                                redundant = true;
                                break;
                            }
                    if (!redundant) edges[i].push_back(j);
                }
            }
        }

        for (int p = 0; p < b.w * b.h; ++p)
            if (b.mine[p]) {
                Factor f{{}, true};
                for (int q : b.adj[p])
                    if (id[q] >= 0) f.vars.push_back(id[q]);
                if (!f.vars.empty()) factors.push_back(std::move(f));
            }

        for (int p : b.targets) {
            Factor f{{}, false};
            if (id[p] >= 0) f.vars.push_back(id[p]);
            for (int q : b.adj[p])
                if (id[q] >= 0) f.vars.push_back(id[q]);
            unique_sort(f.vars);
            if (f.vars.empty()) ++base;
            else factors.push_back(std::move(f));
        }

        for (int i = n; i < n + z; ++i)
            if (edges[i].empty()) {
                ++base;
                ++isolatedZeros;
            }
        {
            for (auto &f : factors) unique_sort(f.vars);
            std::sort(factors.begin(), factors.end(),
                      [](const Factor &a, const Factor &c) { return a.vars < c.vars; });
            std::vector<Factor> merged;

            for (std::size_t i = 0; i < factors.size();) {
                std::size_t end = i;
                int positive = 0, negative = 0;
                while (end < factors.size() && factors[end].vars == factors[i].vars) {
                    (factors[end].mine ? positive : negative) += factors[end].weight;
                    ++end;
                }
                base += std::min(positive, negative);
                if (positive != negative)
                    merged.push_back({factors[i].vars, positive > negative,
                                      std::abs(positive - negative)});
                i = end;
            }
            factors = std::move(merged);
        }
    }

    Model(const Model &parent, const std::vector<int> &variables) : b(parent.b) {
        int oldN = int(parent.cells.size()), n = int(variables.size());
        std::vector<int> vertices = variables, id(parent.edges.size(), -1);

        for (int i = 0; i < n; ++i) {
            id[variables[i]] = i;
            cells.push_back(parent.cells[variables[i]]);
        }

        for (int z = oldN; z < int(parent.edges.size()); ++z) {
            bool incident = false;
            for (int v : parent.edges[z])
                if (id[v] >= 0) {
                    incident = true;
                    break;
                }
            if (incident) {
                id[z] = int(vertices.size());
                vertices.push_back(z);
            }
        }

        edges.resize(vertices.size());
        for (int i = 0; i < int(vertices.size()); ++i)
            for (int v : parent.edges[vertices[i]])
                if (id[v] >= 0) edges[i].push_back(id[v]);

        for (const auto &f : parent.factors)
            if (id[f.vars[0]] >= 0) {
                Factor local{{}, f.mine, f.weight};
                for (int v : f.vars) {
                    if (id[v] < 0)
                        throw std::logic_error("Factor crosses independent pieces");
                    local.vars.push_back(id[v]);
                }
                factors.push_back(std::move(local));
            }
    }

    int evaluate(const std::string &chosen,
                 std::array<int, 4> *counts = nullptr) const {
        int n = (int)cells.size();
        DSU d((int)edges.size());
        std::array<int, 4> c{};

        auto yes = [&](int i) {
            return i >= n || ((uint8_t)chosen[i / 8] >> (i % 8) & 1);
        };

        for (int i = 0; i < (int)edges.size(); ++i)
            if (yes(i)) {
                if (i < n) ++c[1];
                for (int j : edges[i])
                    if (yes(j)) d.join(i, j);
            }

        for (int i = 0; i < (int)edges.size(); ++i)
            if (yes(i) && d.root(i) == i) ++c[2];

        if (!counts) {
            int cost = base - isolatedZeros + c[1] + c[2];
            for (const auto &factor : factors) {
                bool active = false;
                for (int v : factor.vars)
                    if (yes(v)) {
                        active = true;
                        break;
                    }
                if (active == factor.mine) cost += factor.weight;
            }
            return cost;
        }

        std::vector<bool> flags(b.w * b.h), covered(b.w * b.h);
        for (int i = 0; i < n; ++i)
            if (yes(i)) {
                covered[cells[i]] = true;
                for (int q : b.adj[cells[i]])
                    if (b.mine[q]) flags[q] = true;
                    else covered[q] = true;
            }

        c[0] = (int)std::count(flags.begin(), flags.end(), true);

        for (int p : b.targets)
            if (!covered[p]) ++c[3];

        if (counts) *counts = c;
        
        return std::accumulate(c.begin(), c.end(), 0);
    }
};

}  // namespace minclicks::detail
