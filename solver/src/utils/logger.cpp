#include "utils/logger.h"



std::stringstream Logger::buffer;

void Logger::log(Level level, const std::string& message) {
    switch (level) {
        case Level::INFO:
            std::cout << utils::terminal::BOLD_CYAN << "[INFO] " << utils::terminal::RESET;
            buffer << "[INFO] ";
            break;
        case Level::SUCCESS:
            std::cout << utils::terminal::BOLD_GREEN << "[SUCCESS] " << utils::terminal::RESET;
            buffer << "[SUCCESS] ";
            break;
        case Level::ERROR:
            std::cout << utils::terminal::BOLD_RED << "[ERROR] " << utils::terminal::RESET;
            buffer << "[ERROR] ";
            break;
        case Level::DEBUG:
            std::cout << utils::terminal::BOLD_BLUE << "[DEBUG] " << utils::terminal::RESET;
            buffer << "[DEBUG] ";
            break;
        case Level::WARNING:
            std::cout << utils::terminal::BOLD_YELLOW << "[WARNING] " << utils::terminal::RESET;
            buffer << "[WARNING] ";
            break;
        default: break;
    }
    std::cout << message << std::endl;

    buffer << message << "\n";
}


void Logger::mesh_entity(int dim, int tag, int key_id, const std::string& name) {
    std::string key_str = "{" + std::to_string(dim) + ", " + std::to_string(key_id) + "}";
    
    std::cout << utils::terminal::GREEN << "Dim: "      << utils::terminal::RESET << std::left << std::setw(4) << dim 
              << utils::terminal::GREEN << "gmsh id: "  << utils::terminal::RESET << std::left << std::setw(6) << tag 
              << utils::terminal::GREEN << "Key: "      << utils::terminal::RESET << std::left << std::setw(12) << key_str
              << utils::terminal::GREEN << "Name: "     << utils::terminal::RESET << "\"" << name << "\"" 
              << std::endl;
    
    buffer << "ENTITY | " 
           << "Dim: "  << std::left << std::setw(4) << dim     << " | "
           << "ID: "   << std::left << std::setw(6) << tag     << " | "
           << "Key: "  << std::left << std::setw(12) << key_str << " | "
           << "Name: " << name << "\n";
}


void Logger::export_to_file(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << buffer.str();
        file.close();
        Logger::info("Log successfully exported to " + filename);
    }
}