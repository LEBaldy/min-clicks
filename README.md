# Minesweeper minimum clicks

This project contains an exact minesweeper minimum-clicks solver, as well as 3 community heuristics: HZiNi, LZini, and 8-Way Zini.

The solver is a standalone C++20 program with no requirements, that can optionally compile to WebAssembly for use in web applications.

### Acknowledgements

The core DP idea was first introduced by **qqwref**, but was derived from scratch for this project. 8-way, LZiNi and HZiNi are already existing heuristics: LZiNi and HZiNi are implemented as a C++ translation of Llamasweeper's implementation, while 8-way was independently optimized for this project.

AI played a significant role in late-game optimizations, drafting this document, and writing the [min-clicks website](https://min-clicks.netlify.app/) which hosts an interactive version of the solver, allowing this project to be released in days rather than months.

## Building

Requires a C++20 compiler. Can be built like any standard C++ project, for example:

```sh
g++ -std=c++20 -O3 -DNDEBUG -march=native optimal.cpp -o optimal
```

### Compiling to WebAssembly

Compiling to WebAssembly requires Emscripten. The default compilation command is:

```sh
em++ -std=c++20 -O3 -DNDEBUG -fexceptions \
             -sMODULARIZE=1 -sEXPORT_ES6=1 -sINVOKE_RUN=0 \
             -sEXPORTED_RUNTIME_METHODS=callMain,UTF8ToString -sENVIRONMENT=worker \
             -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=67108864 \
             -sMAXIMUM_MEMORY=4294967296 -sSTACK_SIZE=5242880 \
             -sFILESYSTEM=0 -sDYNAMIC_EXECUTION=0 -sEXIT_RUNTIME=0
```

## Usage

The easiest way to run the solver is to provide a board from its Llamasweeper URL/encoded string:

```sh
./optimal [URL/ENCODED_STRING] [options]
```

Or specify board dimensions for a random board using the `--random` option: 

```sh
./optimal --random WIDTH HEIGHT MINES SEED [--64b] [options]
```

The CLI accepts any URL with compatible `b`/`m` parameters, or `b=...&m=...`.
Random dimensions are 1–99 and mine counts 0–width×height.
Seeds are unsigned 32-bit integers by default. Add `--64b` to `--random` to use `mt19937_64` and unsigned 64-bit seeds.

There are **no default time or state limits**. Explicit `--time-limit SECONDS` and `--max-states N` are available; zero means unlimited.

`--quiet` suppresses stderr progress.
`--witness` includes actions with zero-based x/y coordinates: `F` flag, `O` open, `C` chord. 

Exit codes are 0 success, 1 error, and 2 an explicitly limited exact search.

## Machine-readable comparison

`--all` calculates all four algorithms and writes exactly one JSON object to
stdout, with no progress output. Its fields include:

| Field | Meaning |
| --- | --- |
| `url`, `width`, `height`, `mines`, `three_bv` | Board description |
| `status` | `optimal`, or `resource limit` if an explicit limit was reached |
| `minclicks` | Certified minimum, or JSON `null` when limited |
| `upper_bound` | Exact solver's best feasible result |
| `8way`, `lzini`, `hzini` | Community heuristic click counts |
| `peak_states`, `solve_seconds` | Exact-search complexity metrics |

`--all --witness` also returns an `actions` object keyed by `optimal`, `8way`, `lzini`, and `hzini`, each containing `[kind,x,y]` triples.
`--json` alone gives machine-readable output for the selected algorithm.
A heuristic returns `status: "heuristic"`, `algorithm`, `clicks`, and optionally an actions array.
LZini always uses the better of raw LZiNi and HZiNi, including its action sequence. Thus `minclicks <= lzini <= hzini` and `minclicks <= 8way`.

# Algorithm explanation

## Problem description

A board is fully known ahead of time. The goal of the Minesweeper minimum-clicks problem is to determine the smallest number of clicks required to solve it, that is, reveal all safe tiles. A click can either open a cell, flag a cell, or chord (open all surrounding cells of a revealed number if all surrounding flags are placed). The exact solver guarantees the minimum number of clicks.

Call a connected region of zeros a **zero island**. Its **flood region** consists of the zeros and their border numbers. Call a numbered tile with no adjacent zero a **target**: automatic reveals from zeros cannot reveal it, so it needs a direct opening or help from a chord.

Without chords, we need one click per zero island and one per target. This number is called **3BV**, and serves as a baseline to compare against. For any min-clicks algorithm, exact or heuristic, we have `clicks <= 3BV`.

We make heavy use of the **chord set** concept, which states that a solution to the min-clicks problem is fully determined by the set of chorded tiles. Indeed, for such a set, we can charge:

| Type | What we count |
| --- | --- |
| Chords | One per chosen tile |
| Flags | One per mine adjacent to any chosen tile |
| Chain starts | One per connected component of chords |
| Remaining 3BV | One per target neither chosen as a center nor adjacent to a chosen center ("cleanup clicks")|

Let `C` be the set of chorded tiles. Write:

- `M(C)` for the distinct mines requiring flags;
- `k(C)` for the number of components;
- `U(C)` for the remaining uncovered targets;

The solver minimizes this expression over all choices C:

```text
cost(C) = |C| +   |M(C)| + k(C) +   |U(C)|
          chords  flags    starts   remaining openings
```

Choosing no chords gives `cost({}) = 3BV`.

This makes the problem similar to **minimum set cover**, with extra complexity, and structure induced by its minesweeper origins. As such, an exact polynomial-time algorithm is unlikely (albeit not formally proven impossible). Naively enumerating all chord sets is exponential in the number of tiles, which even on beginner board can greatly exceed computational limits. Polynomial-time heuristics do exist, improving greedy approaches slightly. Our algorithm uses row/column-sweep DP to quickly and exactly solve the min-clicks problem on all 3 standard board sizes in `O(max(w, h) * 2^min(w, h))`.

## Dynamic Programming (DP)

Consider processing candidates from left to right. After deciding about the left side, most details far behind us no longer affect the right side. What matters is the information that still crosses the boundary. This allows for a great reduction of the search space.

**Dynamic programming**, abbreviated **DP**, is the systematic use of this idea: solve a sequence of partial problems, store a summary of each distinct situation, and retain the cheapest way of reaching that situation. To implement this, we'll first defined a few key terms.

### Factors and Factor Bits

A **factor** is one small cost rule depending on a particular set of candidate decisions, called its **scope**.

A mine factor says: “Charge one flag if at least one of the neighboring candidates is in `C`.” 
A target factor says: “Charge one click if none of the candidates covering this target is in `C`.” 

A target's covering candidates are its own center and adjacent number tiles. Logically, "at least one" translates to **OR**. Thus, a mine factor contributes `OR(scope)` to the cost, while a target factor contributes `NOR(scope)` to the cost. This is efficiently implemented by a **factor bit**, which remembers the answer to "has any
candidate in this scope been selected so far?". Once the last candidate in the scope has been processed, we can settle the cost and forget the bit.

### Connectivity, Partitions

We also need to remember which partly processed components might connect to future candidates. The basic DP records which retained boundary chords belong to the same component. This grouping is called a **partition**. We give each item a label:

```text
boundary positions: A B C D
labels:             1 2 1 0
```

A and C are already connected, B belongs to a separate component, and D is a skipped chord candidate. The numbers are just group names: `7 4 7 0` describes the same grouping. Always renaming labels by first appearance gives `1 2 1 0` and prevents storing the same situation under different names. This is called **canonicalization**.

Remembering only the number of components would not work: if a future chord touches A and C, it joins one existing component in `1 2 1 0`, but joins two in `1 1 2 0`, which may change the cost.

Only when the last member of a connected component has been processed can we charge one click to start the chord chain, since that component can no longer merge with anything later.

### States, Keys

A **state** is one retained partial situation. In the basic algorithm it has:

- A connectivity partition;
- The current factor bits;
- The clicks count so far;
- Enough information to recover the chosen candidates later.

Partition and factor bits together form a **key**, a unique identifier for each boundary situation. This key allows the dynamic programming algorithm to efficiently look up and update the best known partial solutions.

## A step-by-step description of the DP

For now, use all candidates, no pruning, and charge factors when they finish. Choose any fixed order of candidates.

Start with one empty state. Zero islands without any candidate neighbors each contribute a fixed opening cost; targets without covering candidates do too. Other zero islands enter the frontier just before their first neighboring candidate and leave after their last.

For each candidate in the layer, extend the current state in two ways:

1. **Skip it**, or **select it**, charging one chord in the latter case.
2. If selected, join its component to each selected neighbor or island already represented on the boundary. Set the factor bits whose scopes include it.
3. Settle each factor whose last candidate has now been processed: a mine charges if its bit is 1, a target charges if its bit is 0. Forget those bits.
4. Forget vertices with no future neighbors. Charge one chain start for each component that loses its last representative. Canonicalize.
5. In the next layer, keep only the cheapest entry for each resulting key.

After the final candidate, no boundary information remains, and there is one state left, whose smallest cost is the final answer.

The optimality proof is a one-liner with the following invariant: **after every layer, each stored key has the smallest cost among all partial choices producing that key.**

This approach alone is sufficient to prove boards, but to achieve practicality and efficient integration within other tools, it must be significantly optimized. The following is a list of optimizations that allow the algorithm to break the 1-second barrier on the vast majority of expert boards.

## Optimizations

### Reducing the problem size

Several aspects of the problem can be simplified exactly. These reductions mostly live in [`src/exact/model.hpp`](src/exact/model.hpp).

**Useless chords**: A tile with at most one safe neighbor cannot save clicks. This includes 7s and 8s, as well as edge 4s and 5s and corner 2s and 3s.

**Candidate dominance**: If a candidate is strictly worse than another candidate in terms of coverage and cost, it can be safely ignored.

**Redundant edges**: Two adjacent tiles touching the same zero island are already connected through it. Removing the base edge does not affect connectivity and may allow connectivity items to expire earlier.

**Repeated factors**: Factors that share a scope are combined into a single weighted factor. Opposite costs with identical scopes cancel:

```text
X + (1 - X) = 1
```

For example, three mine factors and two target factors on the same scope become a constant 2 plus one mine factor.

### Splitting independent components

The solver also looks for independent components that may be solved separately. A **separator** is a vertex (currently, only zero islands) whose removal disconnects the graph. When separators exist, we can solve each component independently and recover the combined cost:

```text
combined cost = sum of component costs - corrections for duplicated islands
                + fixed costs
```

The correction is purely geometric, and independent of `C`. In practice, boards are rarely composed of multiple independent components, but this can help on particularly sparse boards or low-density large boards. Implementation lives in [`src/exact/decompose.hpp`](src/exact/decompose.hpp).

### Choosing a small boundary

While order does not affect the final answer, it can significantly impact the size of the state space. The solver compares row and column sweeps in both directions, using a crude `3 * connectivity + 2 * factor` heuristic to choose the best sweep.

Within a single layer, the solver plans the order in which candidates are considered in an inner small DP. One pass minimizes the peak width, a second minimizes the above heuristic. This section is responsible for a large part of the speed-up from the base DP, and lives in [`src/exact/schedule.cpp`](src/exact/schedule.cpp).

### More succint descriptions of boundaries

Suppose chords A and B belong to one component. A touches future candidates `{p, q}`, while B touches only `{q}`. B carries no extra future information and its label can be removed. More generally, by labelling components rather than chords, we increase the amount of equivalences we manage to catch, reducing the state space.

### Pruning states that are always worse

This improvement compares states with the same connectivity but different factor bits. First, the solver changes when mines are charged: it charges a mine factor when first used, retaining its bit to prevent charging it again. This makes a set bits a "credit" for both kinds of factor: a mine already paid for, or a target already covered.

Let `b(A)` denote state A's set bits, and `c(A)` its charged cost. Relative to state B, A can lose at most the total weight of the credits that B has and A lacks:

```text
loss(A, B) = sum of weights of bits in b(B) but not in b(A)

if c(A) + loss(A, B) <= c(B), discard B
```

For example, suppose A costs 20 and B costs 23. If B's only extra "credits" are 2 mine factors (meaning 2 extra flagged mines), A can pay at most two extra clicks later to compensate. A is still no worse under every continuation, so B is unnecessary.

The implementation, in [`src/exact/state.hpp`](src/exact/state.hpp), uses radix grouping to determine the states with equal connectivity on large tables.

A more complex comparison considers different connectivity partitions whose future contact set is equal. It is only used for large layers (> 1024 states).

### Making each state smaller

When possible, state keys are packed into three 64-bit words. This allows up to 24 retained connectivity positions, 31 temporary positions during a transition, and 64 simultaneous factor bits. Larger states use the general engine, which can handle up to 254 connectivity positions and an arbitrary number of factor bits. Many small improvements work together to reduce overhead:

- **Stable factor slots:** Factor bits keep their position during their lifetime.
- **Cached connectivity transitions:** Select/skip transitions are cached by partition, so states with equal connectivity can reuse previous calculations.
- **Small board-data overhead:** each tile's at-most-eight neighbors are stored inline, avoiding many tiny allocations during board construction.

# Results and limits

The solver can handle every beginner, intermediate and expert board tested. Approximate single-core throughput is:

| Board size | Approximate throughput (optimal) | Approximate throughput (8-way) | Approximate throughput (LZiNi/HZiNi) |
| --- | --- | --- | --- |
| Beginner | ~1000/s | ~50000/s | ~15000/s |
| Intermediate | ~25/s | ~12500/s | ~4000/s |
| Expert | ~3/s | ~5000/s | ~1250/s |

