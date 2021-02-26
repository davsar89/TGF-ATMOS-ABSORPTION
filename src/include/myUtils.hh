#pragma once

#include <chrono>
#include <CoordTopocentric.h>
#include <CoordGeodetic.h>
#include <Observer.h>
#include <SGP4.h>
#include <DateTime.h>
#include <vector>
#include <string>
#include <G4ios.hh>
#include <iostream>
#include <fstream>
#include <thread>
#include "Settings.hh"

#include "Randomize.hh"

#include "sys/types.h"
#include "sys/sysinfo.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <uuid/uuid.h>

namespace myUtils
{

long generate_a_unique_ID();

std::vector<double> logspace(const double a, const double b, const int n);

double get_wall_time();

double get_ISS_altitude_at_time(const int year,
                                const int month,
                                const int day,
                                const int hour, // UTC
                                const int minute,
                                const int second,
                                const int microsecond);

///
namespace datetools
{
namespace details
{
constexpr unsigned int days_to_month[2][12] =
    {
        // non-leap year
        {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
        // leap year
        {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335},
};
}

constexpr bool is_leap(int const year) noexcept
{
    return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

constexpr unsigned int day_of_year(int const year, unsigned int const month, unsigned int const day)
{
    return details::days_to_month[is_leap(year)][month - 1] + day;
}
} // namespace datetools
} // namespace myUtils