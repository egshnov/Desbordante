#include "tane.h"

#include "config/error/option.h"
#include "config/error_measure/option.h"
#include "enums.h"
#include "fd/pli_based_fd_algorithm.h"
#include "model/table/column_data.h"

namespace algos {
using Cluster = model::PositionListIndex::Cluster;

Tane::Tane(std::optional<ColumnLayoutRelationDataManager> relation_manager)
    : tane::TaneCommon(relation_manager) {
    RegisterOption(config::kErrorMeasureOpt(&error_measure_));
}

void Tane::MakeExecuteOptsAvailableFDInternal() {
    MakeOptionsAvailable({config::kErrorOpt.GetName(), config::kErrorMeasureOpt.GetName()});
}

config::ErrorType Tane::CalculateZeroAryFdError(ColumnData const* rhs) {
    if (error_measure_ == +ErrorMeasure::g1) return CalculateZeroAryG1(rhs);
    return 1;
}

config::ErrorType Tane::CalculateFdError(model::PositionListIndex const* lhs_pli,
                                         model::PositionListIndex const* rhs_pli,
                                         model::PositionListIndex const* joint_pli) {
    switch (error_measure_) {
        case +ErrorMeasure::pdep:
            return 1 - CalculatePdepMeasure(lhs_pli, joint_pli);
        case +ErrorMeasure::tau:
            return 1 - CalculateTauMeasure(lhs_pli, rhs_pli, joint_pli);
        case +ErrorMeasure::mu_plus:
            return 1 - CalculateMuPlusMeasure(lhs_pli, rhs_pli, joint_pli);
        case +ErrorMeasure::rho:
            return 1 - CalculateRhoMeasure(lhs_pli, joint_pli);
    }
    return CalculateG1Error(lhs_pli, joint_pli);
}

config::ErrorType Tane::CalculateZeroAryG1(ColumnData const* rhs) {
    return 1 - rhs->GetPositionListIndex()->GetNepAsLong() /
                       static_cast<double>(relation_.get()->GetNumTuplePairs());
}

config::ErrorType Tane::CalculateG1Error(model::PositionListIndex const* lhs_pli,
                                         model::PositionListIndex const* joint_pli) {
    return (double)(lhs_pli->GetNepAsLong() - joint_pli->GetNepAsLong()) /
           static_cast<double>(relation_.get()->GetNumTuplePairs());
}

config::ErrorType Tane::PdepSelf(model::PositionListIndex const* x_pli) {
    // model::PositionListIndex const* x_pli = rhs->GetPositionListIndex();
    size_t N = x_pli->GetRelationSize();
    config::ErrorType sum = 0;
    std::size_t cluster_rows_count = 0;
    std::deque<Cluster> const& x_index = x_pli->GetIndex();
    for (Cluster const& x_cluster : x_index) {
        cluster_rows_count += x_cluster.size();
        sum += x_cluster.size() * x_cluster.size();
    }
    std::size_t unique_rows = x_pli->GetRelationSize() - cluster_rows_count;
    sum += unique_rows;
    return static_cast<double>(sum / (N * N));
}

config::ErrorType Tane::CalculatePdepMeasure(model::PositionListIndex const* x_pli,
                                             model::PositionListIndex const* xa_pli) {
    std::deque<Cluster> xa_index = xa_pli->GetIndex();
    std::deque<Cluster> x_index = x_pli->GetIndex();
    size_t N = x_pli->GetRelationSize();

    config::ErrorType sum = 0;

    std::unordered_map<int, unsigned> x_frequencies;

    int x_value_id = 1;
    for (Cluster const& x_cluster : x_index) {
        x_frequencies[x_value_id++] = x_cluster.size();
    }
    x_frequencies[0] = 1;

    auto x_prob = x_pli->CalculateAndGetProbingTable();

    auto get_x_freq_by_tuple_ind{[&x_prob, &x_frequencies](int tuple_ind) {
        int value_id = x_prob->at(tuple_ind);
        return static_cast<config::ErrorType>(x_frequencies[value_id]);
    }};

    for (Cluster const& xa_cluster : xa_index) {
        config::ErrorType num = xa_cluster.size() * xa_cluster.size();
        config::ErrorType denum = get_x_freq_by_tuple_ind(xa_cluster.front());
        sum += num / denum;
    }

    auto xa_prob = xa_pli->CalculateAndGetProbingTable();
    for (int i = 0; i < xa_prob->size(); i++) {
        if (xa_prob->at(i) == 0) {
            sum += 1 / get_x_freq_by_tuple_ind(i);
        }
    }
    return (sum / static_cast<config::ErrorType>(N));
}

config::ErrorType Tane::CalculateTauMeasure(model::PositionListIndex const* x_pli,
                                            model::PositionListIndex const* a_pli,
                                            model::PositionListIndex const* xa_pli) {
    config::ErrorType pdepY = Tane::PdepSelf(a_pli);
    if (pdepY == 1) return 1;
    
    config::ErrorType pdepXY = Tane::CalculatePdepMeasure(x_pli, xa_pli);
    
    return ((pdepXY - pdepY) / (1 - pdepY));
}

config::ErrorType Tane::CalculateMuPlusMeasure(model::PositionListIndex const* x_pli,
                                               model::PositionListIndex const* a_pli,
                                               model::PositionListIndex const* xa_pli) {
    config::ErrorType pdepY = Tane::PdepSelf(a_pli);
    if (pdepY == 1) return 1;
    
    config::ErrorType pdepXY = Tane::CalculatePdepMeasure(x_pli, xa_pli);
    
    size_t N = x_pli->GetRelationSize();
    std::size_t cluster_rows_count = 0;
    std::deque<Cluster> const& x_index = x_pli->GetIndex();
    int K = x_index.size();
    
    for (Cluster const& x_cluster : x_index) {
        cluster_rows_count += x_cluster.size();
    }
    
    std::size_t unique_rows = x_pli->GetRelationSize() - cluster_rows_count;
    K += unique_rows;
    
    if (K == N) return 1;
    
    config::ErrorType mu = 1 - (1 - pdepXY) / (1 - pdepY) * (N - 1) / (N - K);
    config::ErrorType mu_plus = std::max(0., mu);
    return mu_plus;
}

config::ErrorType Tane::CalculateRhoMeasure(model::PositionListIndex const* x_pli,
                                            model::PositionListIndex const* xa_pli) {
    auto CalculateDom = [](model::PositionListIndex const* pli) {
        auto index = pli->GetIndex();
        int dom = index.size();

        std::size_t cluster_rows_count = 0;
        for (Cluster const& cluster : index) {
            cluster_rows_count += cluster.size();
        }
        
        std::size_t unique_rows = pli->GetRelationSize() - cluster_rows_count;
        dom += unique_rows;
        return static_cast<config::ErrorType>(dom);
    };
    config::ErrorType domX = CalculateDom(x_pli);
    config::ErrorType domXA = CalculateDom(xa_pli);
    return domX / domXA;
}

}  // namespace algos
