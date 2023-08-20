#pragma once
#include <iostream>
#include <string>

#include <boost/format.hpp>
#include <boost/functional/hash.hpp>
namespace model {
// ISO 8601 - YYYY-MM-DD
//  BOOST_STRONG_TYPEDEF или отдельный класс для TimeUnit и DateDelta?
using DateDelta = int;  // in days
using TimeUnit = int;
static const int daysInMonth[] = {-1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const int daysBeforeMonth[] = {-1, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
static const int MAXYEAR = 9999;
static const int MINYEAR = 1;
constexpr TimeUnit DI4Y = 4 * 365 + 1;       // 1461 days in 4 years
constexpr TimeUnit DI100Y = 25 * DI4Y - 1;   // 36524 days in 100 years
constexpr TimeUnit DI400Y = 4 * DI100Y + 1;  // 146097 days in 400 years

inline bool IsLeap(const TimeUnit &y) {
    if (y >= MINYEAR && y <= MAXYEAR) return y % 4 == 0 && (y % 100 != 0 || y % 400 == 0);
    throw std::invalid_argument("incorrect year");
}

inline TimeUnit GetDaysBeforeYear(const TimeUnit &year) {
    if (year >= MINYEAR && year <= MAXYEAR) {
        TimeUnit y = year - 1;
        return y * 365 + y / 4 - y / 100 + y / 400;
    }
    throw std::invalid_argument("incorrect year");
}

inline TimeUnit GetDaysBeforeMonth(const TimeUnit &year, const TimeUnit &month) {
    if (year >= MINYEAR && year <= MAXYEAR && month >= 1 && month <= 12)
        return daysBeforeMonth[month] + (month > 2 && IsLeap(year));
    throw std::invalid_argument("incorrect date");
}

// проверить что поданы корректные данные рандомный комментарий

inline TimeUnit GetDaysInMonth(TimeUnit m, TimeUnit y = MINYEAR) {
    if (m >= 1 && m <= 12 && y >= MINYEAR && y <= MAXYEAR) {
        if (IsLeap(y) && m == 2) {
            return 29;
        }
        return daysInMonth[m];
    }
    throw std::invalid_argument("incorrect date");
}
inline bool IsCorrect(TimeUnit y, TimeUnit m = 1, TimeUnit d = 1) {
    return y <= MAXYEAR && y >= MINYEAR && m >= 1 && m <= 12 && d >= 1 && d <= GetDaysInMonth(m, y);
}

class Date {
private:
    TimeUnit day;
    TimeUnit month;
    TimeUnit year;

public:
    Date() noexcept {
        this->day = 0;
        this->month = 0;
        this->year = 0;
    }
    explicit Date(const TimeUnit &year, const TimeUnit &month = 0, const TimeUnit &day = 0) {
        if (IsCorrect(year, month, day) || (year == 0 && month == 0 && day == 0)) {
            this->day = day;
            this->month = month;
            this->year = year;
        } else {
            throw std::invalid_argument("Incorrect date format");
        }
    }
    explicit Date(const std::string &str) {  // only iso format
        if (str.length() != 7 && str.length() != 8 && str.length() != 10) {
            throw std::invalid_argument("string doesn't satisfy date format");
        }
        year = std::stoi(str.substr(0, 4));
        bool has_sep = str[4] == '-';
        size_t pos = 4 + has_sep;
        month = std::stoi(str.substr(pos, 2));
        pos += 2;
        if ((str.substr(pos, 1) == "-") != has_sep) {
            throw std::invalid_argument(
                    "string doesn't satisfy date format: inconsistent use of dash separator");
        }
        pos += has_sep;
        day = std::stoi(str.substr(pos, 2));
    }

    static Date DateFromOrdinal(const TimeUnit &ord) {
        if (ord <= 0) throw std::invalid_argument("ordinal must be more than zero ");
        TimeUnit ordinal = ord - 1;
        TimeUnit year = 1;
        TimeUnit div[] = {DI400Y, DI100Y, DI4Y, 365};
        TimeUnit fact[] = {400, 100, 4, 1};
        TimeUnit res[4];
        for (size_t i = 0; i < 4; i++) {
            res[i] = (ordinal / div[i]);
            year += res[i] * fact[i];
            ordinal %= div[i];
        }
        if (res[3] == 4 || res[1] == 4) return Date(year - 1, 12, 31);
        bool leapyear = IsLeap(year);
        TimeUnit month = (ordinal + 50) >> 5;
        TimeUnit preceding = GetDaysBeforeMonth(year, month);
        if (preceding > ordinal) {
            month--;
            preceding -= daysInMonth[month] + (month == 2 && leapyear);
        }
        ordinal -= preceding;
        return Date(year, month, ordinal + 1);
    }
    [[nodiscard]] bool IsNull() const {
        return year == 0 && month == 0 && day == 0;
    }
    [[nodiscard]] TimeUnit ToOrdinal() const {  // nodiscard?
        return GetDaysBeforeYear(year) + GetDaysBeforeMonth(year, month) + day;
    }

    [[nodiscard]] std::string ToIsoFormatString() const {
        return (boost::format("%04d-%02d-%02d") % year % month % day).str();
    }

    [[nodiscard]] TimeUnit GetDay() const {
        return day;
    }

    [[nodiscard]] TimeUnit GetMonth() const {
        return month;
    }

    [[nodiscard]] TimeUnit GetYear() const {
        return year;
    }

    [[nodiscard]] TimeUnit GetWeekday() const {
        return (this->ToOrdinal() + 6) % 7;
    }

    bool operator==(Date const &right) const {
        return this->day == right.day && this->month == right.month && this->year == right.year;
    }

    bool operator!=(Date const &right) const {
        return !(*this == right);
    }

    bool operator<(Date const &right) const {
        return this->year < right.year || (this->year == right.year && this->month < right.month) ||
               (this->year == right.year && this->month == right.month && this->day < right.day);
    }

    bool operator<=(Date const &right) const {
        return (*this) < right || (*this) == right;
    }

    bool operator>(Date const &right) const {
        return !(*this <= right);
    }

    bool operator>=(Date const &right) const {
        return !(*this < right);
    }

    Date &operator-=(DateDelta const &right) {
        (*this) = DateFromOrdinal(this->ToOrdinal() - right);
        return *this;
    }

    Date &operator+=(DateDelta const &right) {
        (*this) = DateFromOrdinal(this->ToOrdinal() + right);
        return *this;
    }
};

DateDelta operator-(Date const &left, Date const &right);
Date operator+(Date const &left, DateDelta const &right);

Date operator-(Date const &left, DateDelta const &right);
}  // namespace model
namespace std {
template <>
struct hash<model::Date> {
    std::size_t operator()(model::Date const &date) {
        std::size_t seed = 0;
        boost::hash_combine(seed, date.GetDay());
        boost::hash_combine(seed, date.GetMonth());
        boost::hash_combine(seed, date.GetYear());
        return seed;
    }
};
}  // namespace std