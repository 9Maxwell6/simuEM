// terminal_utils.hpp
#pragma once

#include <string>
#include <iomanip>

namespace util {

// Terminal control codes
namespace terminal {
    // Cursor control
    inline const std::string CLEAR_LINE         = "\033[K";           // Clear from cursor to end of line
    inline const std::string CARRIAGE_RETURN    = "\r";          // Move cursor to beginning of line
    inline const std::string CLEAR_SCREEN       = "\033[2J";        // Clear entire screen
    inline const std::string CURSOR_HOME        = "\033[H";          // Move cursor to home position
    
    // Reset
    inline const std::string RESET              = "\033[0m";
    
    // Styles
    inline const std::string BOLD               = "\033[1m";
    inline const std::string DIM                = "\033[2m";
    inline const std::string ITALIC             = "\033[3m";
    inline const std::string UNDERLINE          = "\033[4m";
    inline const std::string BLINK              = "\033[5m";
    inline const std::string REVERSE            = "\033[7m";
    inline const std::string HIDDEN             = "\033[8m";
    inline const std::string STRIKETHROUGH      = "\033[9m";
    
    // Regular colors
    inline const std::string BLACK              = "\033[30m";
    inline const std::string RED                = "\033[31m";
    inline const std::string GREEN              = "\033[32m";
    inline const std::string YELLOW             = "\033[33m";
    inline const std::string BLUE               = "\033[34m";
    inline const std::string MAGENTA            = "\033[35m";
    inline const std::string CYAN               = "\033[36m";
    inline const std::string WHITE              = "\033[37m";
    
    // Bright colors
    inline const std::string BRIGHT_BLACK       = "\033[90m";
    inline const std::string BRIGHT_RED         = "\033[91m";
    inline const std::string BRIGHT_GREEN       = "\033[92m";
    inline const std::string BRIGHT_YELLOW      = "\033[93m";
    inline const std::string BRIGHT_BLUE        = "\033[94m";
    inline const std::string BRIGHT_MAGENTA     = "\033[95m";
    inline const std::string BRIGHT_CYAN        = "\033[96m";
    inline const std::string BRIGHT_WHITE       = "\033[97m";
    
    // Background colors
    inline const std::string BG_BLACK           = "\033[40m";
    inline const std::string BG_RED             = "\033[41m";
    inline const std::string BG_GREEN           = "\033[42m";
    inline const std::string BG_YELLOW          = "\033[43m";
    inline const std::string BG_BLUE            = "\033[44m";
    inline const std::string BG_MAGENTA         = "\033[45m";
    inline const std::string BG_CYAN            = "\033[46m";
    inline const std::string BG_WHITE           = "\033[47m";
    
    // Combined styles (commonly used)
    inline const std::string BOLD_RED           = "\033[1;31m";
    inline const std::string BOLD_GREEN         = "\033[1;32m";
    inline const std::string BOLD_YELLOW        = "\033[1;33m";
    inline const std::string BOLD_BLUE          = "\033[1;34m";
    inline const std::string BOLD_MAGENTA       = "\033[1;35m";
    inline const std::string BOLD_CYAN          = "\033[1;36m";
    inline const std::string BOLD_WHITE         = "\033[1;37m";
    
    // Helper function to clear and update line
    inline std::string update_line(const std::string& text) {
        return CARRIAGE_RETURN + CLEAR_LINE + text;
    }
    
    // Helper function to colorize text
    inline std::string colorize(const std::string& text, const std::string& color) {
        return color + text + RESET;
    }
    
    // Helper function for bold text
    inline std::string bold(const std::string& text) {
        return BOLD + text + RESET;
    }
    
    // Helper function for bold + color
    inline std::string bold_color(const std::string& text, const std::string& color) {
        return BOLD + color + text + RESET;
    }
    
} // namespace terminal

} // namespace utils