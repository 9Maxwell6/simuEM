#pragma once
#include <string>
#include <vector>
#include <regex>


namespace util
{


/**
 * @brief Checks if string `a` contains all strings in at least one group.
 * 
 * Each group is a vector of strings that must ALL be found in `a`.
 * Returns true if ANY group is fully satisfied.
 * 
 * Example:
 *   groups = {{"aaa", "bbb"}, {"ccc"}}
 *   - Returns true if `a` contains both "aaa" AND "bbb" 
 *               OR if `a` contains "ccc"
 * 
 * Special cases:
 *   - {{"foo"}, {"bar"}}        -> equivalent to OR  (any of)
 *   - {{"foo", "bar"}}          -> equivalent to AND (all of)
 *   - {{}}                      -> returns true (empty group is trivially satisfied)
 *   - {}                        -> returns false (no groups to satisfy)
 * 
 * @param a The string to search in.
 * @param b A vector of groups, where each group is a vector of strings.
 * @return true  If all strings in at least one group are found in `a`.
 */
inline bool a_contains_b(const std::string& a, const std::vector<std::vector<std::string>>& b) 
{
    for (const auto& group : b) {
        bool all_found = true;
        for (const auto& s : group) {
            if (a.find(s) == std::string::npos) {
                all_found = false;
                break;
            }
        }
        if (all_found)
            return true;
    }
    return false;
}


/**
 * @brief Extracts the last integer token from an string.
 *
 * @param a  input string.
 *
 * @return   the last integer found in the string, if not integer found, return 0.
 *
 * @note     the integer must be separated with underscore "_".
 */
inline int extract_last_int(const std::string& a)
{
    static const std::regex pattern(R"((?:^|_)(-?\d+)(?:_|$))");

    std::sregex_iterator it(a.begin(), a.end(), pattern);
    std::sregex_iterator end;

    if (it == end)
        return 0;

    std::smatch last_int;
    while (it != end) 
    {
        last_int = *it;
        ++it;
    }

    return std::stoi(last_int[1].str());
}


}