#include "utils/logger.h"



std::stringstream Logger::buffer;

void Logger::log(Level level, const std::string& message)
{
    switch (level)
    {
        case Level::INFO:
            std::cout << util::terminal::BOLD_CYAN << "[INFO] " << util::terminal::RESET;
            buffer << "[INFO] ";
            break;
        case Level::SUCCESS:
            std::cout << util::terminal::BOLD_GREEN << "[SUCCESS] " << util::terminal::RESET;
            buffer << "[SUCCESS] ";
            break;
        case Level::ERROR:
            std::cout << util::terminal::BOLD_RED << "[ERROR] " << util::terminal::RESET;
            buffer << "[ERROR] ";
            break;
        case Level::DEBUG:
            std::cout << util::terminal::BOLD_BLUE << "[DEBUG] " << util::terminal::RESET;
            buffer << "[DEBUG] ";
            break;
        case Level::WARNING:
            std::cout << util::terminal::BOLD_YELLOW << "[WARNING] " << util::terminal::RESET;
            buffer << "[WARNING] ";
            break;
        default: break;
    }
    std::cout << message << std::endl;

    buffer << message << "\n";
}


void Logger::mesh_entity(int dim, int tag, int key_id, int n_elements, int n_types, std::array<size_t, 4> e_size, const std::string& name)
{
    std::string key_str = "{" + std::to_string(dim) + ", " + std::to_string(key_id) + "}";
    
    std::cout << util::terminal::GREEN << "Dim: "       << util::terminal::RESET << std::left << std::setw(4)  << dim         << " | "
              << util::terminal::GREEN << "gmsh id: "   << util::terminal::RESET << std::left << std::setw(3)  << tag         << " | "
              << util::terminal::GREEN << "Key: "       << util::terminal::RESET << std::left << std::setw(6) << key_str     << " | "
              << util::terminal::GREEN << "#elements: " << util::terminal::RESET << std::left << std::setw(4)  << n_elements  << " | "
              << util::terminal::GREEN << "#e-types: "  << util::terminal::RESET << std::left << std::setw(3)  << n_types     << " | "
              << util::terminal::GREEN << "#[node/edge/face/volume]: "   
                    << util::terminal::RESET << std::left << std::setw(16) << (std::to_string(e_size[0]) + "/" 
                                                                             + std::to_string(e_size[1]) + "/" 
                                                                             + std::to_string(e_size[2]) + "/" 
                                                                             + std::to_string(e_size[3]))                     << " | "
              << util::terminal::GREEN << "Name: "      << util::terminal::RESET << "\"" << name << "\"" 
              << std::endl;
    
    buffer << "     GROUP | " 
           << "Dim: "  << std::left << std::setw(4) << dim      << " | "
           << "ID: "   << std::left << std::setw(3) << tag      << " | "
           << "Key: "  << std::left << std::setw(6) << key_str << " | "
           << "#elements: " << std::left << std::setw(4)  << n_elements  << " | "
           << "#e-types: "  << std::left << std::setw(3)  << n_types     << " | "
           << "#[node/edge/face/volume]: " << std::left << std::setw(16) << (std::to_string(e_size[0]) + "/" 
                                                                           + std::to_string(e_size[1]) + "/" 
                                                                           + std::to_string(e_size[2]) + "/" 
                                                                           + std::to_string(e_size[3]))  << " | "
           << "Name: " << name << "\n";
}


void Logger::block_info(size_t id, size_t row_offset, size_t col_offset, size_t row_size, size_t col_size)
{
    std::cout << util::terminal::GREEN << "block id: "   << util::terminal::RESET << std::left << std::setw(3)  << id          << " | "
              << util::terminal::GREEN << "row_offset: " << util::terminal::RESET << std::left << std::setw(3)  << row_offset  << " | "
              << util::terminal::GREEN << "col_offset: " << util::terminal::RESET << std::left << std::setw(3)  << col_offset  << " | "
              << util::terminal::GREEN << "#row: "       << util::terminal::RESET << std::left << std::setw(3)  << row_size  << " | "
              << util::terminal::GREEN << "#col: "       << util::terminal::RESET << std::left << std::setw(3)  << col_size  << " | "
              << std::endl;
    
    buffer << "     BLOCK | " 
           << "block id: "     << std::left << std::setw(3)  << id          << " | "
           << "row_offset: "   << std::left << std::setw(3)  << row_offset  << " | "
           << "col_offset: "   << std::left << std::setw(3)  << col_offset  << " | "
           << "#row: "         << std::left << std::setw(3)  << row_size    << " | "
           << "#col: "         << std::left << std::setw(3)  << col_size    << " | " << "\n";
}

void Logger::start_timer(const std::string& label)
{
    timers_[label] = std::chrono::high_resolution_clock::now();
}

void Logger::stop_timer(const std::string& label)
{
    auto it = timers_.find(label);
    if (it == timers_.end()) {
        Logger::error("Logger::stop_timer - Timer ["+label+"] not found.\n");
        return;
    }
    auto elapsed = std::chrono::high_resolution_clock::now() - it->second;
    auto ms = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    Logger::info("["+label+"] elapsed: " + std::to_string(ms / 1'000'000.0) + "s.");
    timers_.erase(it);
}


void Logger::export_to_file(const std::string& filename)
{
    std::ofstream file(filename);
    if (file.is_open()) {
        file << buffer.str();
        file.close();
        Logger::info("Log successfully exported to " + filename);
    }
}