#pragma once
#include "type.h"
#include "date.h"
namespace model {
class DateType : public Type {
public:
    DateType() noexcept : Type(TypeId::kDate) {}
    std::byte *MakeValueOfOrdinal(TimeUnit const ordinal) const {
        auto *buf = new std::byte[GetSize()];
        GetValue<Date>(buf) = Date::DateFromOrdinal(ordinal);
        return buf;
    }
    [[nodiscard]] std::string ValueToString(std::byte const *value) const override {
        return GetValue<Date>(value).ToIsoFormatString();
    }
    [[nodiscard]] size_t Hash(std::byte const *value) const override {
        return std::hash<Date>{}(GetValue<Date>(value));
    }
    void ValueFromStr(std::byte *dest, std::string s) const override {
        new (dest) Date(s);
    }
    [[nodiscard]] std::unique_ptr<Type> CloneType() const override {
        return std::make_unique<DateType>();
    }
    [[nodiscard]] CompareResult Compare(std::byte const *l, std::byte const *r) const override {
        Date l_val = GetValue<Date>(l);
        Date r_val = GetValue<Date>(r);
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

    static std::byte *AddDelta(std::byte const *date, std::byte const *delta, std::byte *res) {
        GetValue<Date>(res) = GetValue<Date>(date) + GetValue<DateDelta>(delta);
        return res;
    }
    static std::byte *SubDelta(std::byte const *date, std::byte const *delta, std::byte *res) {
        GetValue<Date>(res) = GetValue<Date>(date) - GetValue<DateDelta>(delta);
        return res;
    }
    static std::byte *SubDate(std::byte const *left_date, std::byte const *right_date,
                              std::byte *res) {
        GetValue<DateDelta>(res) = GetValue<Date>(left_date) - GetValue<Date>(right_date);
        return res;
    }
};
}  // namespace model

// void OrdinalConversionTest() {
//     for (size_t y = Date::MINYEAR; y < Date::MAXYEAR; y++) {
//         for (size_t m = 1; m <= 12; m++) {
//             for (size_t d = 1; d <= GetDaysInMonth(m, y); d++) {
//                 Date date(y, m, d);
//                 TimeUnit tmp1 = date.ToOrdinal();
//                 Date tmp2;
//                 bool caught = false;
//                 try {
//                     tmp2 = Date::DateFromOrdinal(tmp1);
//                 } catch (...) {
//                     std::cout << "incorrect Date Format on " << date.ToIsoFormatString() << '\n'
//                     << "ordinal is: "
//                               << tmp1 << '\n';
//                     caught = true;
//                 }
//                 if (date != tmp2 && !caught) {
//                     std::cout << "incorrect answer on " << date.ToIsoFormatString() << '\n' <<
//                     "from ordinal got: "
//                               << tmp2.ToIsoFormatString()
//                               << '\n';
//                 }
//             }
//
//         }
//
//     }
// }
//
// void StringConversionTest() {
//     for (size_t y = Date::MINYEAR; y < Date::MAXYEAR; y++) {
//         for (size_t m = 1; m <= 12; m++) {
//             for (size_t d = 1; d <= GetDaysInMonth(m, y); d++) {
//                 Date date(y, m, d);
//                 std::string tmp1 = date.ToIsoFormatString();
//                 bool caught = false;
//                 Date tmp2;
//                 try {
//                     tmp2 = Date::DateFromString(tmp1);
//                 } catch (...) {
//                     std::cout << "incorrect string format " << date.ToIsoFormatString() << '\n';
//                     caught = true;
//                 }
//                 if (date != tmp2 && !caught) {
//                     std::cout << "incorrect answer on " << date.ToIsoFormatString() << '\n';
//                 }
//             }
//
//         }
//
//     }
// }