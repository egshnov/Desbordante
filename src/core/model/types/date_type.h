#pragma once
#include <boost/date_time/gregorian/gregorian.hpp>

#include "type.h"

namespace model {
using Date = boost::gregorian::date;
using DateDelta = boost::gregorian::date_duration;
class DateType : public Type {
public:
    DateType() noexcept : Type(TypeId::kDate) {}

    [[nodiscard]] std::string ValueToString(std::byte const *value) const override {
        return boost::gregorian::to_iso_extended_string(GetValue<Date>(value));
    }

    [[nodiscard]] size_t Hash(std::byte const *value) const override {
        return GetValue<Date>(value).julian_day();
    }

    void ValueFromStr(std::byte *dest, std::string s) const override {
        new (dest) Date(boost::gregorian::from_simple_string(s));
    }

    [[nodiscard]] std::unique_ptr<Type> CloneType() const override {
        return std::make_unique<DateType>();
    }

    [[nodiscard]] CompareResult Compare(std::byte const *l, std::byte const *r) const override {
        auto const &l_val = GetValue<Date>(l);
        auto const &r_val = GetValue<Date>(r);
        return Compare(l_val, r_val);
    }

    static CompareResult Compare(Date const &l_val, Date const &r_val) {
        if (l_val == r_val) {
            return CompareResult::kEqual;
        }
        if (l_val < r_val) {
            return CompareResult::kLess;
        }
        return CompareResult::kGreater;
    }

    [[nodiscard]] size_t GetSize() const override {
        return sizeof(Date);
    }

    [[nodiscard]] std::byte *MakeValue(Date const date) const {
        auto *buf = new std::byte[GetSize()];
        GetValue<Date>(buf) = date;
        return buf;
    }

    std::byte *AddDelta(std::byte const *date, std::byte const *delta, std::byte *res) {
        GetValue<Date>(res) = GetValue<Date>(date) + GetValue<DateDelta>(delta);
        return res;
    }

    std::byte *SubDelta(std::byte const *date, std::byte const *delta, std::byte *res) {
        GetValue<Date>(res) = GetValue<Date>(date) - GetValue<DateDelta>(delta);
        return res;
    }

    std::byte *SubDate(std::byte const *left_date, std::byte const *right_date, std::byte *res) {
        GetValue<DateDelta>(res) = GetValue<Date>(left_date) - GetValue<Date>(right_date);
        return res;
    }
};
}  // namespace model