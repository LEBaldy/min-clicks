#include "schedule.hpp"

namespace minclicks::detail {

// Optimize the permutation inside each spatial strip.
void optimizeStrips(const Model &m, std::vector<int> &order, bool rows) {
    int n = (int)order.size();
    std::vector<int> pos(n);

    for (int t = 0; t < n; ++t) pos[order[t]] = t;

    auto strip = [&](int v) { return rows ? m.cells[v] / m.b.w : m.cells[v] % m.b.w; };

    struct Block {
        int begin, end;
        std::vector<int> width;
    };
    std::vector<Block> blocks;

    for (int begin = 0; begin < n;) {
        int end = begin + 1;
        while (end < n && strip(order[end]) == strip(order[begin])) ++end;
        if (end - begin > 18) return;
        blocks.push_back({begin, end, {}});
        begin = end;
    }

    int globalPeak = 0;

    for (auto &block : blocks) {
        int begin = block.begin, end = block.end, k = end - begin, size = 1 << k,
            full = size - 1;
        std::vector<int> direct(size), complement(size);
        int constant = 0;

        auto classify = [&](const std::vector<int> &scope) {
            int mask = 0;
            bool before = false, after = false;
            for (int v : scope)
                if (v < n) {
                    int p = pos[v];
                    if (p < begin) before = true;
                    else if (p >= end) after = true;
                    else mask |= 1 << (p - begin);
                }
            return std::tuple{mask, before, after};
        };

        auto hyperedge = [&](const std::vector<int> &scope, int weight) {
            auto [mask, before, after] = classify(scope);
            if (!mask && !before && !after) return;
            // [some processed AND some future] = 1 - [none processed]
            //                                      - [none future].
            constant += weight;
            if (!before) complement[mask] -= weight;
            if (!after) direct[mask] -= weight;
        };

        for (int v = 0; v < n; ++v)
            if (pos[v] < end) {
                auto [mask, before, after] = classify(m.edges[v]);
                (void)before;
                if (pos[v] < begin) {
                    constant += 3;
                    if (!after) direct[mask] -= 3;
                } else {
                    int self = 1 << (pos[v] - begin);
                    direct[self] += 3;
                    if (!after) direct[mask | self] -= 3;
                }
            }

        for (int z = n; z < (int)m.edges.size(); ++z) hyperedge(m.edges[z], 3);

        for (const auto &f : m.factors) hyperedge(f.vars, 2);

        // O(k 2^k) Subset zeta transform
        for (int bit = 1; bit < size; bit <<= 1)
            for (int mask = 0; mask < size; ++mask)
                if (mask & bit) {
                    direct[mask] += direct[mask ^ bit];
                    complement[mask] += complement[mask ^ bit];
                }

        block.width.resize(size);
        std::vector<int> peak(size, std::numeric_limits<int>::max());

        for (int mask = 0; mask < size; ++mask)
            block.width[mask] = constant + direct[mask] + complement[full ^ mask];

        peak[0] = block.width[0];

        for (int mask = 1; mask < size; ++mask) {
            for (unsigned bits = unsigned(mask); bits; bits &= bits - 1) {
                unsigned bit = bits & -bits;
                peak[mask] =
                    std::min(peak[mask], std::max(peak[mask ^ bit], block.width[mask]));
            }
        }

        globalPeak = std::max(globalPeak, peak.back());
    }

    std::vector<int> optimized;
    optimized.reserve(n);

    for (const auto &block : blocks) {
        int size = (int)block.width.size();
        std::vector<double> cost(size, std::numeric_limits<double>::infinity());
        std::vector<unsigned> parent(size);
        cost[0] = 0;

        for (int mask = 1; mask < size; ++mask)
            if (block.width[mask] <= globalPeak) {
                double increment = std::exp2(block.width[mask] / 8.0);
                for (unsigned bits = unsigned(mask); bits; bits &= bits - 1) {
                    unsigned bit = bits & -bits;
                    double value = cost[mask ^ bit] + increment;
                    if (value <= cost[mask]) {
                        cost[mask] = value;
                        parent[mask] = bit;
                    }
                }
            }

        std::vector<int> reverse;

        for (unsigned mask = unsigned(size - 1); mask;) {
            unsigned bit = parent[mask];
            if (!bit)
                throw std::logic_error("Strip planner failed to reconstruct an order");
            reverse.push_back(order[block.begin + std::countr_zero(bit)]);
            mask ^= bit;
        }
        optimized.insert(optimized.end(), reverse.rbegin(), reverse.rend());
    }
    order = std::move(optimized);
}

OrderPlan planOrder(const Model &m, std::string name) {
    int n = (int)m.cells.size(), total = (int)m.edges.size();
    std::vector<int> order(n), pos(n);
    iota(order.begin(), order.end(), 0);
    bool rows = name.find("rows") != std::string::npos,
         reverse = name.find("reverse") != std::string::npos;

    sort(order.begin(), order.end(), [&](int a, int b) {
        int p = m.cells[a], q = m.cells[b];
        int ka = rows ? p : (p % m.b.w) * m.b.h + p / m.b.w;
        int kb = rows ? q : (q % m.b.w) * m.b.h + q / m.b.w;
        return reverse ? ka > kb : ka < kb;
    });

    if (name.find("-dp") != std::string::npos) optimizeStrips(m, order, rows);

    for (int i = 0; i < n; ++i) pos[order[i]] = i;

    std::vector<int> first(total, n), last(total, -1);
    for (int i = 0; i < total; ++i) {
        if (i < n) first[i] = last[i] = pos[i];
        for (int j : m.edges[i])
            if (j < n) {
                first[i] = std::min(first[i], pos[j]);
                last[i] = std::max(last[i], pos[j]);
            }
        if (i < n) first[i] = pos[i];
    }

    std::vector<int> ff(m.factors.size(), n), fl(m.factors.size(), -1);
    for (int f = 0; f < (int)m.factors.size(); ++f)
        for (int i : m.factors[f].vars) {
            ff[f] = std::min(ff[f], pos[i]);
            fl[f] = std::max(fl[f], pos[i]);
        }

    std::vector<int> connDelta(n + 1), factorDelta(n + 1);
    for (int i = 0; i < total; ++i)
        if (first[i] < last[i]) {
            ++connDelta[first[i]];
            --connDelta[last[i]];
        }

    for (int f = 0; f < (int)m.factors.size(); ++f)
        if (ff[f] < fl[f]) {
            ++factorDelta[ff[f]];
            --factorDelta[fl[f]];
        }

    OrderPlan plan{std::move(order), std::move(name)};
    int connections = 0, factors = 0;

    for (int t = 0; t < n; ++t) {
        connections += connDelta[t];
        factors += factorDelta[t];
        double width = 1.5 * connections + factors;
        plan.score = std::max(plan.score, width);
        plan.totalScore += std::exp2(width / 4);
    }

    return plan;
}

Schedule schedule(const Model &m, std::string name) {
    return schedule(m, planOrder(m, std::move(name)));
}

Schedule schedule(const Model &m, const OrderPlan &plan) {
    int n = (int)m.cells.size(), total = (int)m.edges.size();
    const auto &order = plan.order;
    std::vector<int> pos(n);

    for (int i = 0; i < n; ++i) pos[order[i]] = i;

    std::vector<int> first(total, n), last(total, -1);
    for (int i = 0; i < total; ++i) {
        if (i < n) first[i] = last[i] = pos[i];
        for (int j : m.edges[i])
            if (j < n) {
                first[i] = std::min(first[i], pos[j]);
                last[i] = std::max(last[i], pos[j]);
            }
        if (i < n) first[i] = pos[i];
    }

    std::vector<int> ff(m.factors.size(), n), fl(m.factors.size(), -1);
    for (int f = 0; f < (int)m.factors.size(); ++f)
        for (int i : m.factors[f].vars) {
            ff[f] = std::min(ff[f], pos[i]);
            fl[f] = std::max(fl[f], pos[i]);
        }

    Schedule s;
    s.name = plan.name;
    std::vector<int> active, af;
    for (int t = 0; t < n; ++t) {
        Step st;
        st.var = order[t];
        st.before = active;
        st.expanded = active;
        st.oldFactors = af;

        for (int z = n; z < total; ++z)
            if (first[z] == t) {
                st.expanded.push_back(z);
                ++st.freshZeros;
            }

        st.expanded.push_back(st.var);
        for (int j = 0; j < (int)st.expanded.size() - 1; ++j)
            if (find(m.edges[st.var].begin(), m.edges[st.var].end(), st.expanded[j]) !=
                m.edges[st.var].end())
                st.neighbors.push_back(j);

        for (int j = 0; j < (int)st.expanded.size(); ++j)
            if (last[st.expanded[j]] > t) {
                st.keep.push_back(j);
                st.after.push_back(st.expanded[j]);
            }

        std::vector<int> all = af;
        for (int f = 0; f < (int)m.factors.size(); ++f)
            if (ff[f] == t) all.push_back(f);

        for (int f : all) {
            auto it = find(af.begin(), af.end(), f);
            st.factorOld.push_back(it == af.end() ? -1 : (int)(it - af.begin()));
            st.factorTouch.push_back(find(m.factors[f].vars.begin(),
                                          m.factors[f].vars.end(),
                                          st.var) != m.factors[f].vars.end());
            st.factorClose.push_back(fl[f] == t ? (m.factors[f].mine ? 1 : 0) : -1);
            st.factorWeight.push_back(m.factors[f].weight);
            if (fl[f] > t) {
                st.factorKeep.push_back((int)st.newFactors.size());
                st.newFactors.push_back(f);
            } else st.factorKeep.push_back(-1);
        }

        active = st.after;
        af = st.newFactors;
        if (st.after.size() <= 24)
            for (int v : st.after) {
                std::vector<uint32_t> contacts;
                for (int u : m.edges[v])
                    if (u < n && pos[u] > t) {
                        uint32_t mask = 0;
                        for (int j = 0; j < (int)st.after.size(); ++j)
                            if (find(m.edges[u].begin(), m.edges[u].end(),
                                     st.after[j]) != m.edges[u].end())
                                mask |= uint32_t(1) << j;
                        contacts.push_back(mask);
                    }
                st.contacts.push_back(std::move(contacts));
            }

        s.maxConn = std::max(s.maxConn, (int)active.size());
        s.maxFactors = std::max(s.maxFactors, (int)af.size());
        s.score = std::max(s.score, 1.5 * active.size() + af.size());
        s.totalScore += std::exp2((1.5 * active.size() + af.size()) / 4);
        s.steps.push_back(std::move(st));
    }

    std::vector<int> slots(m.factors.size(), -1);
    uint64_t busy = 0;
    int stepIndex = 0;
    for (auto &st : s.steps) {
        if (st.expanded.size() > 31 || st.after.size() > 24) s.packed = false;
        std::vector<int> live = st.oldFactors;
        for (int f : st.newFactors)
            if (find(live.begin(), live.end(), f) == live.end()) live.push_back(f);

        // Include factors that both start and finish at this candidate.
        for (int f = 0; f < (int)m.factors.size(); ++f)
            if (find(m.factors[f].vars.begin(), m.factors[f].vars.end(), st.var) !=
                    m.factors[f].vars.end() &&
                find(live.begin(), live.end(), f) == live.end())
                live.push_back(f);

        for (int f : live)
            if (slots[f] < 0) {
                if (busy == ~uint64_t(0)) {
                    s.packed = false;
                    return s;
                }
                slots[f] = std::countr_zero(~busy);
                busy |= uint64_t(1) << slots[f];
            }

        for (int f : live) {
            uint64_t bit = uint64_t(1) << slots[f];
            if (m.factors[f].mine) st.mineMask |= bit;
            if (m.factors[f].mine && m.factors[f].weight > 1)
                st.extraMineWeights.push_back({bit, m.factors[f].weight - 1});
            if (std::find(m.factors[f].vars.begin(), m.factors[f].vars.end(), st.var) !=
                m.factors[f].vars.end())
                st.touchMask |= bit;
            if (std::find(st.newFactors.begin(), st.newFactors.end(), f) ==
                st.newFactors.end()) {
                (m.factors[f].mine ? st.closeMines : st.closeTargets) |= bit;
                if (!m.factors[f].mine && m.factors[f].weight > 1)
                    st.extraTargetWeights.push_back({bit, m.factors[f].weight - 1});
                busy &= ~bit;
                slots[f] = -1;
            }
        }

        // Feasible dual
        for (int f : st.newFactors)
            st.liveWeights.push_back({uint64_t(1) << slots[f], m.factors[f].weight});

        std::vector<std::pair<int, int>> remaining;
        for (int f = 0; f < (int)m.factors.size(); ++f)
            if (!m.factors[f].mine && fl[f] > stepIndex) {
                int degree = 0;
                for (int v : m.factors[f].vars) degree += pos[v] > stepIndex;
                remaining.push_back({degree, f});
            }

        sort(remaining.begin(), remaining.end());
        std::vector<int> capacity(n, 12);
        std::array<uint64_t, 13> masks{};
        for (auto [degree, f] : remaining) {
            (void)degree;
            int weight = 12;
            for (int v : m.factors[f].vars)
                if (pos[v] > stepIndex) weight = std::min(weight, capacity[v]);
            if (!weight) continue;
            for (int v : m.factors[f].vars)
                if (pos[v] > stepIndex) capacity[v] -= weight;
            if (ff[f] > stepIndex) st.boundConstant += weight;
            else masks[weight] |= uint64_t(1) << slots[f];
        }

        for (int w = 1; w <= 12; ++w)
            if (masks[w]) st.boundMasks.push_back({masks[w], w});

        {
            std::map<std::pair<std::vector<int>, int>, FactorFold> groups;
            for (int f : st.newFactors) {
                std::vector<int> scope;
                for (int v : m.factors[f].vars)
                    if (pos[v] > stepIndex) scope.push_back(v);
                std::sort(scope.begin(), scope.end());
                auto &group = groups[{scope, m.factors[f].weight}];
                group.weight = m.factors[f].weight;
                (m.factors[f].mine ? group.mines : group.targets) |= uint64_t(1)
                                                                     << slots[f];
            }
            
            for (auto &[scope, group] : groups) {
                (void)scope;
                int mines = std::popcount(group.mines),
                    targets = std::popcount(group.targets);
                if (mines + targets < 2) continue;
                group.targetCount = targets;
                for (int balance = -targets; balance <= mines; ++balance) {
                    uint64_t remaining = balance < 0 ? group.targets : group.mines,
                             zeros = 0;
                    for (int k = 0; k < std::abs(balance); ++k) {
                        zeros |= remaining & -remaining;
                        remaining &= remaining - 1;
                    }
                    group.zeroMasks.push_back(zeros);
                }
                st.folds.push_back(std::move(group));
            }
        }

        for (int z = n; z < total; ++z)
            if (!m.edges[z].empty() && first[z] > stepIndex) st.futureZero = true;
        ++stepIndex;
    }
    return s;
}

}  // namespace minclicks::detail
