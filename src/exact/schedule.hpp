#pragma once

#include "model.hpp"

namespace minclicks::detail {

struct FactorFold {
    uint64_t mines = 0, targets = 0;
    int weight = 1, targetCount = 0;
    std::vector<uint64_t> zeroMasks;
};

struct Step {
    int var;
    std::vector<int> before, expanded, after, keep, neighbors;
    std::vector<int> oldFactors, newFactors, factorOld, factorTouch, factorClose,
        factorKeep, factorWeight;
    int freshZeros = 0;
    uint64_t touchMask = 0, closeMines = 0, closeTargets = 0;
    uint64_t mineMask = 0;
    std::vector<std::pair<uint64_t, int>> extraMineWeights, extraTargetWeights;
    std::vector<std::pair<uint64_t, int>> boundMasks;
    int boundConstant = 0;
    bool futureZero = false;
    std::vector<FactorFold> folds;
    std::vector<std::vector<uint32_t>> contacts;
    std::vector<std::pair<uint64_t, int>> liveWeights;
};

struct Schedule {
    std::vector<Step> steps;
    int maxConn = 0, maxFactors = 0;
    double score = 0;
    double totalScore = 0;
    std::string name;
    bool packed = true;
};

struct OrderPlan {
    std::vector<int> order;
    std::string name;
    double score = 0, totalScore = 0;
};

OrderPlan planOrder(const Model &model, std::string name);

Schedule schedule(const Model &model, const OrderPlan &plan);

void optimizeStrips(const Model &model, std::vector<int> &order, bool rows);

Schedule schedule(const Model &model, std::string name);

} // namespace minclicks::detail
