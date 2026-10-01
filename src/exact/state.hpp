#pragma once

#include "schedule.hpp"

namespace minclicks::detail {

struct PackedKey {
    uint64_t lo = 0, hi = 0, factors = 0;
    bool operator==(const PackedKey &) const = default;

    unsigned label(int i) const {
        return unsigned((i < 12 ? lo >> (5 * i) : hi >> (5 * (i - 12))) & 31);
    }

    void put(int i, unsigned label) {
        (i < 12 ? lo : hi) |= uint64_t(label) << (5 * (i < 12 ? i : i - 12));
    }
};

inline uint64_t mix(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline uint64_t hashKey(const PackedKey &k) {
    return mix(k.lo + 0x9e3779b97f4a7c15ULL) ^ std::rotl(mix(k.hi), 21) ^
           std::rotl(mix(k.factors), 42);
}

struct ConnectionResult {
    uint64_t lo = 0, hi = 0;
    int seeds = 0;
    bool live = false;
};

struct ConnectionCacheEntry {
    uint64_t lo = 0, hi = 0;
    int generation = 0;
    std::array<ConnectionResult, 2> next;
};

// A union of boundary contact sets is uniquely described by all boundary
// positions whose contacts it contains.
inline uint32_t contactClosure(uint32_t members, const Step &st) {
    uint32_t closure = members;
    for (int v = 0; v < int(st.contacts.size()); ++v) {
        if (members >> v & 1) continue;
        bool contained = true;
        for (uint32_t contact : st.contacts[v])
            if (!(contact & members)) {
                contained = false;
                break;
            }
        if (contained) closure |= uint32_t(1) << v;
    }
    return closure;
}

struct ContactPartition {
    std::vector<uint32_t> components;
    uint32_t united = 0;

    ContactPartition(const PackedKey &key, const Step &st) {
        std::array<uint32_t, 32> members{};
        uint32_t all = 0;

        for (int v = 0; v < int(st.contacts.size()); ++v)
            if (unsigned label = key.label(v))
                all |= members[label] |= uint32_t(1) << v;

        united = contactClosure(all, st);
        for (auto mask : members)
            if (mask) components.push_back(contactClosure(mask, st));

        std::sort(components.begin(), components.end());
    }
};

inline bool contactCoarsens(const ContactPartition &coarse,
                            const ContactPartition &fine) {
    if (coarse.united != fine.united ||
        coarse.components.size() > fine.components.size())
        return false;

    if (coarse.components.empty()) return fine.components.empty();

    if (coarse.components.size() == 1) return true;

    std::array<uint32_t, 24> eligible{};
    uint32_t covered = 0;

    for (std::size_t c = 0; c < coarse.components.size(); ++c)
        for (std::size_t f = 0; f < fine.components.size(); ++f)
            if (!(fine.components[f] & ~coarse.components[c]))
                covered |= eligible[c] |= uint32_t(1) << f;

    if (covered != (uint32_t(1) << fine.components.size()) - 1) return false;

    std::array<int, 24> match;
    match.fill(-1);

    auto augment = [&](auto &&self, int c, uint32_t &seen) -> bool {
        uint32_t available = eligible[c] & ~seen;
        while (available) {
            int f = std::countr_zero(available);
            available &= available - 1;
            seen |= uint32_t(1) << f;
            if (match[f] < 0 || self(self, match[f], seen)) {
                match[f] = c;
                return true;
            }
        }
        return false;
    };

    for (int c = 0; c < int(coarse.components.size()); ++c) {
        uint32_t seen = 0;
        if (!augment(augment, c, seen)) return false;
    }

    return true;
}

inline std::array<ConnectionResult, 2> connectionTransitions(const PackedKey &key,
                                                             const Step &st) {
    int old = (int)st.before.size(), count = (int)st.expanded.size();
    std::array<uint8_t, 32> original{};
    unsigned highest = 0;

    for (int j = 0; j < old; ++j) {
        original[j] = uint8_t(key.label(j));
        highest = std::max(highest, unsigned(original[j]));
    }

    for (int j = old; j < count - 1; ++j) original[j] = uint8_t(++highest);

    std::array<ConnectionResult, 2> result;
    for (int take = 0; take <= 1; ++take) {
        auto labels = original;
        if (take) {
            unsigned component = highest + 1;
            labels[count - 1] = uint8_t(component);
            uint32_t merged = 0;
            for (int j : st.neighbors) merged |= uint32_t(1) << labels[j];
            merged &= ~uint32_t(1);
            for (int j = 0; j < count - 1; ++j)
                if (merged >> labels[j] & 1) labels[j] = uint8_t(component);
        }
        uint32_t present = 0, surviving = 0;

        for (int j = 0; j < count; ++j) present |= uint32_t(1) << labels[j];
        PackedKey out;
        std::array<uint8_t, 32> rename{};
        unsigned fresh = 0;

        for (int j = 0; j < (int)st.keep.size(); ++j) {
            unsigned label = labels[st.keep[j]];
            surviving |= uint32_t(1) << label;
            if (label && !rename[label]) rename[label] = uint8_t(++fresh);
            out.put(j, rename[label]);
        }
        {
            std::array<uint32_t, 32> members{};
            for (int j = 0; j < (int)st.keep.size(); ++j)
                if (out.label(j)) members[out.label(j)] |= uint32_t(1) << j;
            for (int j = 0; j < (int)st.keep.size(); ++j)
                if (!out.label(j)) {
                    for (unsigned label = 1; label <= fresh; ++label) {
                        bool covered = true;
                        for (uint32_t contact : st.contacts[j])
                            if (!(members[label] & contact)) {
                                covered = false;
                                break;
                            }
                        if (covered) {
                            out.put(j, label);
                            members[label] |= uint32_t(1) << j;
                            break;
                        }
                    }
                }
            for (int j = (int)st.keep.size() - 1; j >= 0; --j)
                if (unsigned label = out.label(j)) {
                    uint32_t others = members[label] & ~(uint32_t(1) << j);
                    bool redundant = others != 0;
                    for (uint32_t contact : st.contacts[j])
                        if (!(others & contact)) {
                            redundant = false;
                            break;
                        }
                    if (redundant) members[label] = others;
                }
            PackedKey reduced;
            rename.fill(0);
            fresh = 0;
            for (int j = 0; j < (int)st.keep.size(); ++j) {
                unsigned label = out.label(j);
                if (!(members[label] >> j & 1)) continue;
                if (!rename[label]) rename[label] = uint8_t(++fresh);
                reduced.put(j, rename[label]);
            }
            out = reduced;
        }
        result[take] = {out.lo, out.hi,
                        std::popcount(present & ~surviving & ~uint32_t(1)),
                        bool(surviving & ~uint32_t(1))};
    }
    return result;
}

// Dense state storage with an open-addressed index.
struct FlatTable {
    struct Node {
        PackedKey key;
        int cost;
    };

    std::vector<Node> nodes;
    std::vector<uint64_t> choices;
    std::vector<uint32_t> buckets;
    int words;
    explicit FlatTable(int nw) : words(nw) { buckets.resize(16); }

    void clear() {
        nodes.clear();
        choices.clear();
        std::fill(buckets.begin(), buckets.end(), 0);
    }

    void rehash(std::size_t size) {
        buckets.assign(size, 0);
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            std::size_t slot = hashKey(nodes[i].key) & (size - 1);
            while (buckets[slot]) slot = (slot + 1) & (size - 1);
            buckets[slot] = uint32_t(i + 1);
        }
    }

    void insert(const PackedKey &key, int cost, const uint64_t *source, int take = -1) {
        if ((nodes.size() + 1) * 4 > buckets.size() * 3) rehash(buckets.size() * 2);
        std::size_t slot = hashKey(key) & (buckets.size() - 1);

        while (buckets[slot] && !(nodes[buckets[slot] - 1].key == key))
            slot = (slot + 1) & (buckets.size() - 1);
        std::size_t index;

        if (buckets[slot]) {
            index = buckets[slot] - 1;
            if (nodes[index].cost <= cost) return;
            nodes[index].cost = cost;
        } else {
            index = nodes.size();
            nodes.push_back({key, cost});
            buckets[slot] = uint32_t(index + 1);
            choices.resize(choices.size() + words);
        }

        for (int j = 0; j < words; ++j) choices[index * words + j] = source[j];

        if (take >= 0) choices[index * words + take / 64] |= uint64_t(1) << (take % 64);
    }

    std::size_t pruneFactors(const Step &st) {
        constexpr std::size_t crossThreshold = 1024;
        if (nodes.size() < 2 ||
            (st.liveWeights.empty() &&
             (st.contacts.empty() || nodes.size() < crossThreshold)))
            return 0;

        // For equal connectivity, c(A)+weight(P(B) minus P(A)) <= c(B)
        // proves that A is no worse for every continuation.
        std::vector<uint32_t> order(nodes.size());
        iota(order.begin(), order.end(), 0);
        auto cheaper = [&](uint32_t a, uint32_t b) {
            if (nodes[a].cost != nodes[b].cost) return nodes[a].cost < nodes[b].cost;
            return nodes[a].key.factors > nodes[b].key.factors;
        };

        if (order.size() < 256) {
            std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
                const auto &x = nodes[a];
                const auto &y = nodes[b];
                if (x.key.lo != y.key.lo) return x.key.lo < y.key.lo;
                if (x.key.hi != y.key.hi) return x.key.hi < y.key.hi;
                return cheaper(a, b);
            });
        } else {
            // Group equal connectivity with byte radix passes, then only sort
            // costs/masks within each group. Skip bytes constant in this layer.
            uint64_t differentLo = 0, differentHi = 0;
            for (const auto &node : nodes) {
                differentLo |= node.key.lo ^ nodes[0].key.lo;
                differentHi |= node.key.hi ^ nodes[0].key.hi;
            }

            std::vector<uint32_t> scratch(order.size());
            for (int high = 1; high >= 0; --high)
                for (int shift = 0; shift < 64; shift += 8) {
                    if (!(((high ? differentHi : differentLo) >> shift) & 255))
                        continue;
                    std::array<std::size_t, 256> positions{};
                    auto byte = [&](uint32_t index) {
                        const auto &key = nodes[index].key;
                        return unsigned(((high ? key.hi : key.lo) >> shift) & 255);
                    };
                    for (uint32_t index : order) ++positions[byte(index)];
                    std::size_t sum = 0;
                    for (auto &count : positions) {
                        std::size_t next = sum + count;
                        count = sum;
                        sum = next;
                    }
                    for (uint32_t index : order)
                        scratch[positions[byte(index)]++] = index;
                    order.swap(scratch);
                }

            for (std::size_t first = 0; first < order.size();) {
                std::size_t last = first + 1;
                const auto &key = nodes[order[first]].key;
                while (last < order.size() && nodes[order[last]].key.lo == key.lo &&
                       nodes[order[last]].key.hi == key.hi)
                    ++last;
                std::sort(order.begin() + first, order.begin() + last, cheaper);
                first = last;
            }
        }

        std::vector<std::pair<uint64_t, int>> extra;
        for (auto [bit, w] : st.liveWeights)
            if (w > 1) extra.push_back({bit, w - 1});

        std::vector<uint8_t> removed(nodes.size());
        std::vector<uint32_t> survivors;
        uint64_t lastLo = ~uint64_t(0), lastHi = ~uint64_t(0);
        for (uint32_t i : order) {
            const auto &node = nodes[i];
            if (lastLo != node.key.lo || lastHi != node.key.hi) {
                survivors.clear();
                lastLo = node.key.lo;
                lastHi = node.key.hi;
            }

            for (uint32_t j : survivors) {
                uint64_t missing = node.key.factors & ~nodes[j].key.factors;
                if (!missing) {
                    removed[i] = 1;
                    break;
                }
                if (nodes[j].cost == node.cost) continue;
                int penalty = std::popcount(missing);
                for (auto [bit, w] : extra)
                    if (missing & bit) penalty += w;
                if (nodes[j].cost + penalty <= node.cost) {
                    removed[i] = 1;
                    break;
                }
            }
            if (!removed[i]) survivors.push_back(i);
        }

        if (nodes.size() >= crossThreshold) {
            struct Group {
                ContactPartition contacts;
                std::vector<uint32_t> states;
                int population = 0;
            };

            std::vector<Group> groups;
            for (std::size_t first = 0; first < order.size();) {
                std::size_t last = first + 1;
                const auto &key = nodes[order[first]].key;
                while (last < order.size() && nodes[order[last]].key.lo == key.lo &&
                       nodes[order[last]].key.hi == key.hi)
                    ++last;
                Group group{ContactPartition(key, st), {}};
                for (uint32_t component : group.contacts.components)
                    group.population += std::popcount(component);
                for (std::size_t j = first; j < last; ++j)
                    if (!removed[order[j]]) group.states.push_back(order[j]);
                if (!group.states.empty()) groups.push_back(std::move(group));
                first = last;
            }

            std::sort(groups.begin(), groups.end(), [](const Group &a, const Group &b) {
                if (a.contacts.united != b.contacts.united)
                    return a.contacts.united < b.contacts.united;
                if (a.contacts.components.size() != b.contacts.components.size())
                    return a.contacts.components.size() < b.contacts.components.size();
                if (a.population != b.population) return a.population > b.population;
                return a.contacts.components < b.contacts.components;
            });

            for (std::size_t start = 0; start < groups.size();) {
                std::size_t end = start + 1;
                while (end < groups.size() &&
                       groups[end].contacts.united == groups[start].contacts.united)
                    ++end;
                for (std::size_t f = start + 1; f < end; ++f) {
                    auto &fine = groups[f];
                    for (std::size_t c = start; c < f; ++c) {
                        const auto &coarse = groups[c];
                        if (coarse.contacts.components.size() >
                            fine.contacts.components.size())
                            break;
                        if (nodes[coarse.states.front()].cost >
                                nodes[fine.states.back()].cost ||
                            !contactCoarsens(coarse.contacts, fine.contacts))
                            continue;
                        for (uint32_t target : fine.states) {
                            if (removed[target]) continue;
                            for (uint32_t source : coarse.states) {
                                int allowance = nodes[target].cost - nodes[source].cost;
                                if (allowance < 0) break;
                                uint64_t missing = nodes[target].key.factors &
                                                   ~nodes[source].key.factors;
                                int penalty = std::popcount(missing);
                                for (auto [bit, weight] : extra)
                                    if (missing & bit) penalty += weight;
                                if (penalty <= allowance) {
                                    removed[target] = 1;
                                    break;
                                }
                            }
                        }
                    }
                }
                start = end;
            }
        }

        std::size_t kept = 0, old = nodes.size();
        for (std::size_t i = 0; i < old; ++i)
            if (!removed[i]) {
                if (kept != i) {
                    nodes[kept] = nodes[i];
                    for (int j = 0; j < words; ++j)
                        choices[kept * words + j] = choices[i * words + j];
                }
                ++kept;
            }
        nodes.resize(kept);
        choices.resize(kept * words);

        return old - kept;
    }
};

}  // namespace minclicks::detail
