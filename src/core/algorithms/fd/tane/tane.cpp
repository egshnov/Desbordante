#include "tane.h"

#include <chrono>
#include <iomanip>
#include <list>
#include <memory>

#include <easylogging++.h>

#include "config/error/option.h"
#include "config/error_measure/option.h"
#include "enums.h"
#include "fd/pli_based_fd_algorithm.h"
#include "fd/tane/lattice_level.h"
#include "fd/tane/lattice_vertex.h"
#include "model/table/column_data.h"
#include "model/table/column_layout_relation_data.h"
#include "model/table/relational_schema.h"

namespace algos {
using boost::dynamic_bitset;
using Cluster = model::PositionListIndex::Cluster;

void tane::TaneCommon::ResetStateFd() {}

tane::TaneCommon::TaneCommon(std::optional<ColumnLayoutRelationDataManager> relation_manager)
    : PliBasedFDAlgorithm({kDefaultPhaseName}, relation_manager) {}

double tane::TaneCommon::CalculateUccError(model::PositionListIndex const* pli,
                                           ColumnLayoutRelationData const* relation_data) {
    return pli->GetNepAsLong() / static_cast<double>(relation_data->GetNumTuplePairs());
}

void tane::TaneCommon::RegisterAndCountFd(Vertical const& lhs, Column const* rhs,
                                          [[maybe_unused]] config::ErrorType error,
                                          [[maybe_unused]] RelationalSchema const* schema) {
    dynamic_bitset<> lhs_bitset = lhs.GetColumnIndices();
    PliBasedFDAlgorithm::RegisterFd(lhs, *rhs);
}

void tane::TaneCommon::Prune(model::LatticeLevel* level) {
    RelationalSchema const* schema = relation_->GetSchema();
    std::list<model::LatticeVertex*> key_vertices;
    for (auto& [map_key, vertex] : level->GetVertices()) {
        Vertical columns = vertex->GetVertical();  // Originally it's a ColumnCombination

        if (vertex->GetIsKeyCandidate()) {
            double ucc_error = CalculateUccError(vertex->GetPositionListIndex(), relation_.get());
            if (ucc_error == 0) {  // if a (super) key
                vertex->SetKeyCandidate(false);

                for (std::size_t rhs_index = vertex->GetRhsCandidates().find_first();
                     rhs_index != boost::dynamic_bitset<>::npos;
                     rhs_index = vertex->GetRhsCandidates().find_next(rhs_index)) {
                    Vertical rhs = static_cast<Vertical>(*schema->GetColumn((int)rhs_index));
                    if (!columns.Contains(rhs)) {
                        bool is_rhs_candidate = true;
                        for (auto const& column : columns.GetColumns()) {
                            Vertical sibling =
                                    columns.Without(static_cast<Vertical>(*column)).Union(rhs);
                            auto sibling_vertex =
                                    level->GetLatticeVertex(sibling.GetColumnIndices());
                            if (sibling_vertex == nullptr ||
                                !sibling_vertex->GetConstRhsCandidates()[rhs.GetColumnIndices()
                                                                                 .find_first()]) {
                                is_rhs_candidate = false;
                                break;
                            }
                            // for each outer rhs: if there is a sibling s.t. it doesn't
                            // have this rhs, there is no FD: vertex->rhs
                        }
                        // Found fd: vertex->rhs => register it
                        if (is_rhs_candidate) {
                            RegisterAndCountFd(columns, schema->GetColumn(rhs_index), 0, schema);
                        }
                    }
                }
                key_vertices.push_back(vertex.get());
            }
        }
        for (auto key_vertex : key_vertices) {
            key_vertex->GetRhsCandidates() &= key_vertex->GetVertical().GetColumnIndices();
            key_vertex->SetInvalid(true);
        }
    }
}

void tane::TaneCommon::ComputeDependencies(model::LatticeLevel* level) {
    RelationalSchema const* schema = relation_->GetSchema();
    for (auto& [key_map, xa_vertex] : level->GetVertices()) {
        if (xa_vertex->GetIsInvalid()) {
            continue;
        }
        Vertical xa = xa_vertex->GetVertical();
        // Calculate XA PLI
        if (xa_vertex->GetPositionListIndex() == nullptr) {
            auto parent_pli_1 = xa_vertex->GetParents()[0]->GetPositionListIndex();
            auto parent_pli_2 = xa_vertex->GetParents()[1]->GetPositionListIndex();
            xa_vertex->AcquirePositionListIndex(parent_pli_1->Intersect(parent_pli_2));
        }

        dynamic_bitset<> xa_indices = xa.GetColumnIndices();
        dynamic_bitset<> a_candidates = xa_vertex->GetRhsCandidates();
        auto xa_pli = xa_vertex->GetPositionListIndex();
        for (auto const& x_vertex : xa_vertex->GetParents()) {
            Vertical const& lhs = x_vertex->GetVertical();

            // Find index of A in XA. If a is not a candidate, continue. TODO: possible to do it
            // easier??
            // like "a_index = xa_indices - x_indices;"
            int a_index = xa_indices.find_first();
            dynamic_bitset<> x_indices = lhs.GetColumnIndices();
            while (a_index >= 0 && x_indices[a_index]) {
                a_index = xa_indices.find_next(a_index);
            }
            if (!a_candidates[a_index]) {
                continue;
            }
            auto x_pli = x_vertex->GetPositionListIndex();
            auto a_pli = relation_->GetColumnData(a_index).GetPositionListIndex();
            // Check X -> A
            config::ErrorType error = CalculateFdError(x_pli, a_pli, xa_pli);
            if (error <= max_fd_error_) {
                Column const* rhs = schema->GetColumns()[a_index].get();

                RegisterAndCountFd(lhs, rhs, error, schema);
                xa_vertex->GetRhsCandidates().set(rhs->GetIndex(), false);
                if (error == 0) {
                    xa_vertex->GetRhsCandidates() &= lhs.GetColumnIndices();
                }
            }
        }
    }
}

unsigned long long tane::TaneCommon::ExecuteInternal() {
    long apriori_millis = 0;
    RelationalSchema const* schema = relation_->GetSchema();

    LOG(DEBUG) << schema->GetName() << " has " << relation_->GetNumColumns() << " columns, "
               << relation_->GetNumRows() << " rows, and a maximum NIP of " << std::setw(2)
               << relation_->GetMaximumNip() << ".";

    for (auto& column : schema->GetColumns()) {
        double avg_partners = relation_->GetColumnData(column->GetIndex())
                                      .GetPositionListIndex()
                                      ->GetNepAsLong() *
                              2.0 / relation_->GetNumRows();
        LOG(DEBUG) << "* " << column->ToString() << ": every tuple has " << std::setw(2)
                   << avg_partners << " partners on average.";
    }
    auto start_time = std::chrono::system_clock::now();
    double progress_step = 100.0 / (schema->GetNumColumns() + 1);

    // Initialize level 0
    std::vector<std::unique_ptr<model::LatticeLevel>> levels;
    auto level0 = std::make_unique<model::LatticeLevel>(0);
    // TODO: через указатели кажется надо переделать
    level0->Add(std::make_unique<model::LatticeVertex>(*(schema->empty_vertical_)));
    model::LatticeVertex const* empty_vertex = level0->GetVertices().begin()->second.get();
    levels.push_back(std::move(level0));
    AddProgress(progress_step);

    // Initialize level1
    dynamic_bitset<> zeroary_fd_rhs(schema->GetNumColumns());
    auto level1 = std::make_unique<model::LatticeLevel>(1);
    for (auto& column : schema->GetColumns()) {
        // for each attribute set vertex
        ColumnData const& column_data = relation_->GetColumnData(column->GetIndex());
        auto vertex = std::make_unique<model::LatticeVertex>(static_cast<Vertical>(*column));

        vertex->AddRhsCandidates(schema->GetColumns());
        vertex->GetParents().push_back(empty_vertex);
        vertex->SetKeyCandidate(true);
        vertex->SetPositionListIndex(column_data.GetPositionListIndex());

        // check FDs: 0->A
        double fd_error = CalculateZeroAryFdError(&column_data);
        if (fd_error <= max_fd_error_) {  // TODO: max_error
            zeroary_fd_rhs.set(column->GetIndex());
            RegisterAndCountFd(*schema->empty_vertical_, column.get(), fd_error, schema);

            vertex->GetRhsCandidates().set(column->GetIndex(), false);
            if (fd_error == 0) {
                vertex->GetRhsCandidates().reset();
            }
        }

        level1->Add(std::move(vertex));
    }

    for (auto& [key_map, vertex] : level1->GetVertices()) {
        Vertical column = vertex->GetVertical();
        vertex->GetRhsCandidates() &=
                ~zeroary_fd_rhs;  //~ returns flipped copy <- removed already discovered zeroary FDs

        // вот тут костыль, чтобы вытянуть индекс колонки из вершины, в которой только один индекс
        ColumnData const& column_data =
                relation_->GetColumnData(column.GetColumnIndices().find_first());
        double ucc_error = CalculateUccError(column_data.GetPositionListIndex(), relation_.get());

        if (ucc_error == 0 && max_lhs_ != 0) {  // if a (super) key
            vertex->SetKeyCandidate(false);
            for (unsigned long rhs_index = vertex->GetRhsCandidates().find_first();
                 rhs_index < vertex->GetRhsCandidates().size();
                 rhs_index = vertex->GetRhsCandidates().find_next(rhs_index)) {
                if (rhs_index != column.GetColumnIndices().find_first()) {
                    RegisterAndCountFd(column, schema->GetColumn(rhs_index), 0, schema);
                }
            }
            vertex->GetRhsCandidates() &= column.GetColumnIndices();
            vertex->SetInvalid(true);
        }
    }
    levels.push_back(std::move(level1));
    AddProgress(progress_step);

    unsigned int max_arity =
            max_lhs_ == std::numeric_limits<unsigned int>::max() ? max_lhs_ : max_lhs_ + 1;
    for (unsigned int arity = 2; arity <= max_arity; arity++) {
        model::LatticeLevel::ClearLevelsBelow(levels, arity - 1);
        model::LatticeLevel::GenerateNextLevel(levels);

        model::LatticeLevel* level = levels[arity].get();
        LOG(TRACE) << "Checking " << level->GetVertices().size() << " " << arity
                   << "-ary lattice vertices.";
        if (level->GetVertices().empty()) {
            break;
        }

        ComputeDependencies(level);

        if (arity == max_arity) {
            break;
        }

        Prune(level);
        // TODO: printProfilingData
        AddProgress(progress_step);
    }

    SetProgress(100);
    std::chrono::milliseconds elapsed_milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() -
                                                                  start_time);
    apriori_millis += elapsed_milliseconds.count();

    LOG(DEBUG) << "Time: " << apriori_millis << " milliseconds";
    LOG(DEBUG) << "Intersection time: " << model::PositionListIndex::micros_ / 1000 << "ms";
    LOG(DEBUG) << "Total intersections: " << model::PositionListIndex::intersection_count_
               << std::endl;
    LOG(DEBUG) << "Total FD count: " << fd_collection_.Size();
    LOG(DEBUG) << "HASH: " << Fletcher16();
    return apriori_millis;
}

Tane::Tane(std::optional<ColumnLayoutRelationDataManager> relation_manager)
    : tane::TaneCommon(relation_manager) {
    RegisterOptions();
}

void Tane::RegisterOptions() {
    RegisterOption(config::kErrorOpt(&max_fd_error_));
    RegisterOption(config::kErrorMeasureOpt(&error_measure_));
}

void Tane::MakeExecuteOptsAvailableFDInternal() {
    MakeOptionsAvailable({config::kErrorOpt.GetName(), config::kErrorMeasureOpt.GetName()});
}

config::ErrorType Tane::CalculateZeroAryG1(ColumnData const* rhs) {
    return 1 - rhs->GetPositionListIndex()->GetNepAsLong() /
                       static_cast<double>(relation_.get()->GetNumTuplePairs());
}

config::ErrorType Tane::CalculateZeroAryFdError(ColumnData const* rhs) {
    if (error_measure_ == +ErrorMeasure::g1) return CalculateZeroAryG1(rhs);
    return 1;
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

    // NOTE: Сломается если изменить механизм обхода в CalculateAndGetProbingTable
    // но должно работать быстрее чем вычислять частоты на основе PT
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
    auto CalculateDom{[](model::PositionListIndex const* pli) {
        auto index = pli->GetIndex();
        int dom = index.size();
        std::size_t cluster_rows_count = 0;
        for (Cluster const& cluster : index) {
            cluster_rows_count += cluster.size();
        }
        std::size_t unique_rows = pli->GetRelationSize() - cluster_rows_count;
        dom += unique_rows;
        return static_cast<config::ErrorType>(dom);
    }};
    config::ErrorType domX = CalculateDom(x_pli);
    config::ErrorType domXA = CalculateDom(xa_pli);
    return domX / domXA;
}

config::ErrorType Tane::CalculateFdError(model::PositionListIndex const* lhs_pli,
                                         model::PositionListIndex const* rhs_pli,
                                         model::PositionListIndex const* joint_pli) {
    switch (error_measure_) {
        case +ErrorMeasure::g1:
            return CalculateG1Error(lhs_pli, joint_pli);
        case +ErrorMeasure::pdep:
            return 1 - CalculatePdepMeasure(lhs_pli, joint_pli);
        case +ErrorMeasure::tau:
            return 1 - CalculateTauMeasure(lhs_pli, rhs_pli, joint_pli);
        case +ErrorMeasure::mu_plus:
            return 1 - CalculateMuPlusMeasure(lhs_pli, rhs_pli, joint_pli);
        case +ErrorMeasure::rho:
            return 1 - CalculateRhoMeasure(lhs_pli, joint_pli);
    }
}

PFDTane::PFDTane(std::optional<ColumnLayoutRelationDataManager> relation_manager)
    : tane::TaneCommon(relation_manager) {
    RegisterOptions();
}

config::ErrorType PFDTane::CalculateZeroAryPFDError(ColumnData const* rhs,
                                                    ColumnLayoutRelationData const*) {
    std::size_t max = 1;
    model::PositionListIndex const* x_pli = rhs->GetPositionListIndex();
    for (Cluster const& x_cluster : x_pli->GetIndex()) {
        max = std::max(max, x_cluster.size());
    }
    return 1.0 - static_cast<double>(max) / x_pli->GetRelationSize();
}

config::ErrorType PFDTane::CalculatePFDError(model::PositionListIndex const* x_pli,
                                             model::PositionListIndex const* xa_pli,
                                             ErrorMeasure measure,
                                             ColumnLayoutRelationData const*) {
    std::deque<Cluster> xa_index = xa_pli->GetIndex();
    std::shared_ptr<Cluster const> probing_table = x_pli->CalculateAndGetProbingTable();
    std::sort(xa_index.begin(), xa_index.end(),
              [&probing_table](Cluster const& a, Cluster const& b) {
                  return probing_table->at(a.front()) < probing_table->at(b.front());
              });
    double sum = 0.0;
    std::size_t cluster_rows_count = 0;
    std::deque<Cluster> const& x_index = x_pli->GetIndex();
    auto xa_cluster_it = xa_index.begin();

    for (Cluster const& x_cluster : x_index) {
        std::size_t max = 1;
        for (int x_row : x_cluster) {
            if (xa_cluster_it == xa_index.end()) {
                break;
            }
            if (x_row == xa_cluster_it->at(0)) {
                max = std::max(max, xa_cluster_it->size());
                xa_cluster_it++;
            }
        }
        sum += measure == +ErrorMeasure::per_tuple ? static_cast<double>(max)
                                                   : static_cast<double>(max) / x_cluster.size();
        cluster_rows_count += x_cluster.size();
    }
    unsigned int unique_rows =
            static_cast<unsigned int>(x_pli->GetRelationSize() - cluster_rows_count);
    double probability = static_cast<double>(sum + unique_rows) /
                         (measure == +ErrorMeasure::per_tuple ? x_pli->GetRelationSize()
                                                              : x_index.size() + unique_rows);
    return 1.0 - probability;
}

config::ErrorType PFDTane::CalculateZeroAryFdError(ColumnData const* rhs) {
    return CalculateZeroAryPFDError(rhs, relation_.get());
}

config::ErrorType PFDTane::CalculateFdError(model::PositionListIndex const* lhs_pli,
                                            model::PositionListIndex const* rhs_pli,
                                            model::PositionListIndex const* joint_pli) {
    return CalculatePFDError(lhs_pli, joint_pli, error_measure_, relation_.get());
}

void PFDTane::RegisterOptions() {
    RegisterOption(config::kErrorOpt(&max_fd_error_));
    RegisterOption(config::kErrorMeasureOpt(&error_measure_));
}

void PFDTane::MakeExecuteOptsAvailableFDInternal() {
    MakeOptionsAvailable({config::kErrorOpt.GetName(), config::kErrorMeasureOpt.GetName()});
}

}  // namespace algos
