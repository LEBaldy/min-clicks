#include "../exact.hpp"
#include "decompose.hpp"
#include "state.hpp"

namespace minclicks::detail {

using Options = ExactOptions;

struct Entry {
    int cost;
    std::string chosen;
};

struct Result {
    bool exact = false;
    int best = 0;
    std::string chosen, reason;
    std::size_t peak = 1, transitions = 0;
    std::size_t dominated = 0;
    int completed = 0;
    double elapsed = 0;
};

static Result initial(const Model &m, const Options &opt) {
    Result result;
    int n = (int)m.cells.size();
    result.chosen = std::string((n + 7) / 8, '\0');
    result.best = m.evaluate(result.chosen);
    opt.report("heuristic", 0, n, 0, 0, result.best, 0, 0, true);

    // Deterministic improving one-flip search gives an initial feasible bound.
    for (;;) {
        int best = result.best, arg = -1;
        for (int i = 0; i < n; ++i) {
            flip(result.chosen, i);
            int c = m.evaluate(result.chosen);
            flip(result.chosen, i);
            if (c < best) {
                best = c;
                arg = i;
            }
            if ((i & 63) == 0) opt.report("heuristic", 0, n, 0, 0, result.best, 0, 0);
        }
        if (arg < 0) break;
        flip(result.chosen, arg);
        result.best = best;
    }

    std::mt19937 rng(0x5eed);
    std::vector<int> order(n);
    iota(order.begin(), order.end(), 0);

    for (int start = 0; start < 4; ++start) {
        std::string chosen((n + 7) / 8, '\0');
        for (int i = 0; i < n; ++i)
            if (start == 0 || rng() % 4 < unsigned(start)) flip(chosen, i);
        int cost = m.evaluate(chosen);
        bool improved = true;
        for (int pass = 0; pass < 12 && improved; ++pass) {
            improved = false;
            for (int j = n - 1; j > 0; --j)
                std::swap(order[j], order[rng() % unsigned(j + 1)]);
            for (int i : order) {
                bool selected = ((uint8_t)chosen[i / 8] >> (i % 8) & 1);
                flip(chosen, i);
                int value = m.evaluate(chosen);
                if (value < cost || (value == cost && selected)) {
                    improved |= value < cost;
                    cost = value;
                } else flip(chosen, i);
            }
        }
        if (cost < result.best) {
            result.best = cost;
            result.chosen = std::move(chosen);
        }
        opt.report("heuristic", 0, n, 0, 0, result.best, 0, 0);
    }

    return result;
}

static Result solve_reference(const Model &m, const Schedule &sched,
                              const Options &opt) {
    auto start = Clock::now();
    Result result = initial(m, opt);
    int n = (int)m.cells.size();

    if (opt.progress) std::cerr << "Initial upper bound: " << result.best << '\n';

    std::unordered_map<std::string, Entry> states, next;
    states.emplace("", Entry{m.base, std::string((n + 7) / 8, '\0')});
    double lastLog = seconds(start);

    for (const Step &st : sched.steps) {
        if (opt.timeLimit > 0 && seconds(start) >= opt.timeLimit) {
            result.reason = "time limit";
            break;
        }

        // Byte labels permit up to 254 simultaneous connectivity entries.
        if (st.expanded.size() > 254) {
            result.reason = "connectivity frontier exceeds 254 entries";
            break;
        }
        next.clear();
        next.reserve(opt.maxStates ? std::min(opt.maxStates, states.size() * 2)
                                   : states.size() * 2);
        bool stopped = false;

        for (const auto &[key, entry] : states) {
            for (int take = 0; take <= 1; ++take) {
                ++result.transitions;
                std::array<uint8_t, 256> labels{}, rename{}, present{}, survive{};
                uint8_t highest = 0;
                int old = (int)st.before.size(), count = (int)st.expanded.size();

                for (int j = 0; j < old; ++j) {
                    labels[j] = (uint8_t)key[j];
                    highest = std::max(highest, labels[j]);
                }

                for (int j = old; j < count - 1; ++j) labels[j] = ++highest;

                if (take) {
                    uint8_t component = ++highest;
                    labels[count - 1] = component;
                    for (int j : st.neighbors)
                        if (labels[j]) {
                            uint8_t merged = labels[j];
                            for (int k = 0; k < count; ++k)
                                if (labels[k] == merged) labels[k] = component;
                        }
                }

                for (int j = 0; j < count; ++j) present[labels[j]] = 1;

                for (int j : st.keep) survive[labels[j]] = 1;

                int cost = entry.cost + take;
                for (int j = 1; j <= highest; ++j) cost += present[j] && !survive[j];

                int nc = (int)st.after.size();
                std::string out(nc + (st.newFactors.size() + 7) / 8, '\0');
                uint8_t fresh = 0;
                for (int j = 0; j < nc; ++j) {
                    uint8_t label = labels[st.keep[j]];
                    if (label && !rename[label]) rename[label] = ++fresh;
                    out[j] = char(rename[label]);
                }

                for (int f = 0; f < (int)st.factorOld.size(); ++f) {
                    int prev = st.factorOld[f];
                    bool used =
                        (take && st.factorTouch[f]) ||
                        (prev >= 0 && ((uint8_t)key[old + prev / 8] >> (prev % 8) & 1));
                    if (st.factorClose[f] >= 0)
                        cost += (used == bool(st.factorClose[f])) * st.factorWeight[f];
                    else if (used) {
                        int j = st.factorKeep[f];
                        out[nc + j / 8] |= char(1u << (j % 8));
                    }
                }

                if (cost >= result.best) continue;

                auto it = next.find(out);
                if (it == next.end()) {
                    std::string chosen = entry.chosen;
                    if (take) flip(chosen, st.var);
                    next.emplace(std::move(out), Entry{cost, std::move(chosen)});
                    if (opt.maxStates && next.size() > opt.maxStates) {
                        result.reason = "state limit";
                        stopped = true;
                        break;
                    }
                } else if (cost < it->second.cost) {
                    it->second.cost = cost;
                    it->second.chosen = entry.chosen;
                    if (take) flip(it->second.chosen, st.var);
                }
            }

            if (stopped) break;

            if (opt.onProgress && (result.transitions & 16383) == 0)
                opt.report("search", result.completed, n, next.size(),
                           std::max(result.peak, next.size()), result.best,
                           seconds(start), result.transitions);

            if ((result.transitions & 16383) == 0 && opt.timeLimit > 0 &&
                seconds(start) >= opt.timeLimit) {
                result.reason = "time limit";
                stopped = true;
                break;
            }
        }

        result.peak = std::max(result.peak, next.size());
        if (stopped) break;

        states.swap(next);
        ++result.completed;

        opt.report("search", result.completed, n, states.size(), result.peak,
                   result.best, seconds(start), result.transitions);

        if (opt.progress && seconds(start) - lastLog >= 2) {
            std::cerr << "DP " << result.completed << '/' << n << ": " << states.size()
                      << " states, " << seconds(start) << " s\n";
            lastLog = seconds(start);
        }

        if (states.empty()) {
            result.exact = true;
            break;
        }
    }

    if (result.completed == n) {
        result.exact = true;
        if (!states.empty()) {
            const auto &entry = states.at("");
            if (entry.cost < result.best) {
                result.best = entry.cost;
                result.chosen = entry.chosen;
            }
        }
    }

    result.elapsed = seconds(start);
    return result;
}

static Result solve_packed(const Model &m, const Schedule &sched, const Options &opt) {
    auto start = Clock::now();
    Result result = initial(m, opt);
    int n = (int)m.cells.size(), words = (n + 63) / 64;

    if (opt.progress) std::cerr << "Initial upper bound: " << result.best << '\n';

    FlatTable states(words), next(words);
    std::vector<uint64_t> empty(words);
    // Direct mapped: collisions only cause recomputation, never approximation.
    std::vector<ConnectionCacheEntry> connectionCache(64);
    std::vector<uint64_t> future(words, ~uint64_t(0));
    states.insert({}, m.base, empty.data());
    double lastLog = seconds(start);

    for (const Step &st : sched.steps) {
        if (opt.timeLimit > 0 && seconds(start) >= opt.timeLimit) {
            result.reason = "time limit";
            break;
        }

        std::size_t cacheSize =
            states.nodes.size() >= 16384
                ? 32768
                : std::bit_ceil(std::max(std::size_t(64), states.nodes.size() * 2));
        if (connectionCache.size() < cacheSize) connectionCache.resize(cacheSize);
        next.clear();
        bool stopped = false;

        for (std::size_t index = 0; index < states.nodes.size(); ++index) {
            const auto &entry = states.nodes[index];
            std::array<ConnectionResult, 2> connections;
            {
                std::size_t slot = mix(entry.key.lo ^ std::rotl(entry.key.hi, 27)) &
                                   (connectionCache.size() - 1);
                auto &cached = connectionCache[slot];
                if (cached.generation != result.completed + 1 ||
                    cached.lo != entry.key.lo || cached.hi != entry.key.hi) {
                    cached.lo = entry.key.lo;
                    cached.hi = entry.key.hi;
                    cached.generation = result.completed + 1;
                    cached.next = connectionTransitions(entry.key, st);
                }
                connections = cached.next;
            }

            for (int take = 0; take <= 1; ++take) {
                ++result.transitions;
                const auto &connection = connections[take];
                PackedKey out{connection.lo, connection.hi, 0};
                uint64_t factors = entry.key.factors | (take ? st.touchMask : 0);

                // Charge a flag on its first use, while retaining its OR bit
                // until last use so later chords share it for free.
                int cost = entry.cost + take + connection.seeds +
                           std::popcount((factors & ~entry.key.factors) & st.mineMask) +
                           std::popcount(~factors & st.closeTargets);

                for (auto [mask, weight] : st.extraMineWeights)
                    if ((factors & ~entry.key.factors) & mask) cost += weight;

                for (auto [mask, weight] : st.extraTargetWeights)
                    if (~factors & mask) cost += weight;

                if (cost >= result.best) continue;

                out.factors = factors & ~(st.closeMines | st.closeTargets);

                for (const auto &group : st.folds) {
                    int mines = std::popcount(group.mines & ~out.factors),
                        targets = std::popcount(group.targets & ~out.factors);
                    cost += std::min(mines, targets) * group.weight;
                    out.factors = (out.factors | group.mines | group.targets) &
                                  ~group.zeroMasks[mines - targets + group.targetCount];
                }
                {
                    int lower = st.boundConstant;
                    for (auto [mask, weight] : st.boundMasks)
                        lower += std::popcount(~out.factors & mask) * weight;

                    lower = (lower + 11) / 12 + int(st.futureZero || connection.live);
                    if (cost + lower >= result.best) continue;
                }
                
                next.insert(out, cost, states.choices.data() + index * words,
                            take ? st.var : -1);

                if (opt.maxStates && next.nodes.size() > opt.maxStates) {
                    result.reason = "state limit";
                    stopped = true;
                    break;
                }
            }

            if (stopped) break;

            if (opt.onProgress && (result.transitions & 16383) == 0)
                opt.report("search", result.completed, n, next.nodes.size(),
                           std::max(result.peak, next.nodes.size()), result.best,
                           seconds(start), result.transitions);

            if ((result.transitions & 16383) == 0 && opt.timeLimit > 0 &&
                seconds(start) >= opt.timeLimit) {
                result.reason = "time limit";
                stopped = true;
                break;
            }
        }
        result.peak = std::max(result.peak, next.nodes.size());
        
        if (!stopped) result.dominated += next.pruneFactors(st);

        future[st.var / 64] &= ~(uint64_t(1) << (st.var % 64));

        if (!next.nodes.empty() &&
            (stopped || result.completed % 8 == 7 || result.completed == n - 1)) {
            std::size_t minimum = 0;
            for (std::size_t i = 1; i < next.nodes.size(); ++i)
                if (next.nodes[i].cost < next.nodes[minimum].cost) minimum = i;
            std::size_t stride = std::max(std::size_t(1), next.nodes.size() / 128);
            auto complete = [&](std::size_t index) {
                std::string chosen = result.chosen;
                for (int i = 0; i < n; ++i) {
                    if (!(future[i / 64] >> (i % 64) & 1)) {
                        bool take =
                            next.choices[index * words + i / 64] >> (i % 64) & 1;
                        bool had = (uint8_t)chosen[i / 8] >> (i % 8) & 1;
                        if (take != had) flip(chosen, i);
                    }
                }
                int cost = m.evaluate(chosen);
                if (cost < result.best) {
                    result.best = cost;
                    result.chosen = std::move(chosen);
                }
            };
            complete(minimum);
            for (std::size_t i = 0; i < next.nodes.size(); i += stride) complete(i);
        }

        if (stopped) break;
        std::swap(states, next);
        ++result.completed;

        opt.report("search", result.completed, n, states.nodes.size(), result.peak,
                   result.best, seconds(start), result.transitions);

        if (opt.progress && seconds(start) - lastLog >= 2) {
            std::cerr << "DP " << result.completed << '/' << n << ": "
                      << states.nodes.size() << " states, " << seconds(start) << " s\n";
            lastLog = seconds(start);
        }

        if (states.nodes.empty()) {
            result.exact = true;
            break;
        }
    }

    if (result.completed == n) {
        result.exact = true;
        if (!states.nodes.empty() && states.nodes[0].cost < result.best) {
            result.best = states.nodes[0].cost;
            fill(result.chosen.begin(), result.chosen.end(), '\0');
            for (int i = 0; i < n; ++i)
                if (states.choices[i / 64] >> (i % 64) & 1) flip(result.chosen, i);
        }
    }

    result.elapsed = seconds(start);
    return result;
}
static std::vector<Action> witness(const Model &m, const std::string &chosen) {
    const Board &b = m.b;
    int n = (int)m.cells.size();
    std::vector<bool> selected(b.w * b.h), flags(b.w * b.h), open(b.w * b.h),
        done(b.w * b.h);
    std::vector<Action> actions;

    for (int i = 0; i < n; ++i)
        if ((uint8_t)chosen[i / 8] >> (i % 8) & 1) {
            selected[m.cells[i]] = true;
            for (int q : b.adj[m.cells[i]])
                if (b.mine[q]) flags[q] = true;
        }

    for (int p = 0; p < b.w * b.h; ++p)
        if (flags[p]) actions.push_back({'F', p});

    auto reveal = [&](int p) {
        if (b.mine[p]) throw std::logic_error("Witness opens a mine");
        open[p] = true;
        if (b.island[p] >= 0) {
            int z = b.island[p];
            for (int q : b.zeros[z]) open[q] = true;
            for (int q : b.borders[z]) open[q] = true;
        }
    };

    auto propagate = [&]() {
        bool changed = true;
        while (changed) {
            changed = false;
            for (int p = 0; p < b.w * b.h; ++p)
                if (selected[p] && open[p] && !done[p]) {
                    actions.push_back({'C', p});
                    done[p] = true;
                    changed = true;
                    for (int q : b.adj[p])
                        if (!b.mine[q]) reveal(q);
                }
        }
    };

    // Prefer an island as seed whenever the component contains one.
    for (const auto &z : b.zeros)
        if (!open[z[0]]) {
            actions.push_back({'O', z[0]});
            reveal(z[0]);
            propagate();
        }

    for (int p = 0; p < b.w * b.h; ++p)
        if (selected[p] && !done[p]) {
            actions.push_back({'O', p});
            reveal(p);
            propagate();
        }

    for (int p = 0; p < b.w * b.h; ++p)
        if (!b.mine[p] && !open[p]) {
            actions.push_back({'O', p});
            reveal(p);
        }

    return actions;
}

static Result searchModel(const Model &m, const Options &opt, Schedule &sched) {
    int emptyCost = m.evaluate(std::string((m.cells.size() + 7) / 8, '\0'));
    OrderPlan plan = planOrder(m, "columns");
    opt.report("ordering", 0, int(m.cells.size()), 0, 0, emptyCost, 0, 0, true);
    std::vector<std::string> orders = {"rows", "columns-reverse", "rows-reverse"};

    if (m.b.h <= m.b.w)
        orders.insert(orders.end(), {"columns-dp", "columns-dp-reverse"});

    if (m.b.w <= m.b.h) orders.insert(orders.end(), {"rows-dp", "rows-dp-reverse"});

    for (const auto &name : orders) {
        opt.report("ordering", 0, int(m.cells.size()), 0, 0, emptyCost, 0, 0);
        auto candidate = planOrder(m, name);
        if (candidate.score < plan.score ||
            (candidate.score == plan.score && candidate.totalScore < plan.totalScore))
            plan = std::move(candidate);
    }

    sched = schedule(m, plan);
    return sched.packed ? solve_packed(m, sched, opt) : solve_reference(m, sched, opt);
}

}  // namespace minclicks::detail

namespace minclicks {

ExactResult solveExact(const Board &b, const ExactOptions &opt) {
    using namespace detail;
    auto started = Clock::now();
    opt.report("preparing", 0, 0, 0, 0, -1, 0, 0, true);
    Model m(b);
    auto groups = independentParts(m);
    ExactResult out;
    std::string chosen((m.cells.size() + 7) / 8, '\0');

    if (groups.size() <= 1) {
        Schedule sched;
        auto r = searchModel(m, opt, sched);
        out.exact = r.exact;
        out.clicks = r.best;
        out.reason = r.reason;
        chosen = std::move(r.chosen);
        out.peak = r.peak;
        out.transitions = r.transitions;
        out.dominated = r.dominated;
        out.completed = r.completed;
        out.elapsed = r.elapsed;
        out.order = sched.name;
        out.engine = sched.packed ? "packed" : "reference";
        out.maxConn = sched.maxConn;
        out.maxFactors = sched.maxFactors;
    } else {
        std::vector<Model> pieces;
        pieces.reserve(groups.size());
        std::vector<int> costs;
        int zeroCopies = 0;

        for (const auto &vars : groups) {
            pieces.emplace_back(m, vars);
            const auto &piece = pieces.back();
            zeroCopies += int(piece.edges.size() - piece.cells.size());
            costs.push_back(piece.evaluate(std::string((vars.size() + 7) / 8, '\0')));
        }

        int nonisolatedZeros = int(m.edges.size() - m.cells.size()) - m.isolatedZeros;
        int constant = m.base - (zeroCopies - nonisolatedZeros);
        int best = constant + accumulate(costs.begin(), costs.end(), 0);
        out.exact = true;
        out.engine = "packed";
        out.order = "components(" + std::to_string(pieces.size()) + ")";

        for (std::size_t i = 0; i < pieces.size(); ++i) {
            ExactOptions local = opt;

            if (opt.timeLimit > 0)
                local.timeLimit = std::max(1e-9, opt.timeLimit - out.elapsed);

            if (opt.onProgress)
                local.onProgress = [&](const Progress &p) {
                    bool hasBound = std::string(p.phase) == "heuristic" ||
                                    std::string(p.phase) == "search";
                    int upper =
                        hasBound && p.best >= 0 ? best - costs[i] + p.best : best;
                    opt.report(p.phase, out.completed + p.completed,
                               int(m.cells.size()), p.states,
                               std::max(out.peak, p.peak), upper, seconds(started),
                               out.transitions + p.transitions,
                               std::string(p.phase) == "heuristic" && p.completed == 0);
                };
            Schedule sched;
            auto r = searchModel(pieces[i], local, sched);
            best += r.best - costs[i];
            costs[i] = r.best;

            for (int v = 0; v < int(groups[i].size()); ++v)
                if ((uint8_t)r.chosen[v / 8] >> (v % 8) & 1) flip(chosen, groups[i][v]);

            out.peak = std::max(out.peak, r.peak);
            out.transitions += r.transitions;
            out.dominated += r.dominated;
            out.completed += r.completed;
            out.elapsed += r.elapsed;
            out.maxConn = std::max(out.maxConn, sched.maxConn);
            out.maxFactors = std::max(out.maxFactors, sched.maxFactors);

            if (!sched.packed) out.engine = "reference";

            if (!r.exact) {
                out.exact = false;
                out.reason = r.reason;
                break;
            }
        }
        out.clicks = best;
    }

    if (m.evaluate(chosen, &out.breakdown) != out.clicks)
        throw std::logic_error("Internal objective mismatch");

    if (opt.includeActions) {
        out.actions = witness(m, chosen);
        if (int(out.actions.size()) != out.clicks)
            throw std::logic_error("Internal witness cost mismatch");
    }

    out.candidates = int(m.cells.size());
    out.removed = m.removed;
    opt.report(out.exact ? "complete" : "limited", out.completed, out.candidates, 0,
               out.peak, out.clicks, out.elapsed, out.transitions, true);
    return out;
}

}  // namespace minclicks
