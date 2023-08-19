#include "date.h"
namespace model {
DateDelta operator-(Date const &left, Date const &right) {
    return left.ToOrdinal() - right.ToOrdinal();
}

Date operator+(Date const &left, DateDelta const &right) {
    return Date::DateFromOrdinal(left.ToOrdinal() + right);
}

Date operator-(Date const &left, DateDelta const &right) {
    return Date::DateFromOrdinal(left.ToOrdinal() - right);
}
}  // namespace model