#pragma once

#include "terminal.h"
#include "config.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

class Logger 
{
private:
    static std::stringstream buffer;

public:
    enum class Level 
    { 
        INFO, 
        SUCCESS, 
        ERROR, 
        DEBUG,
        WARNING
    };

    // Standard log messages with [INFO] etc. prefixes
    static void log(Level level, const std::string& message);
    
    // Convenience wrappers
    static void info(const std::string& msg)      { log(Level::INFO,    msg); }
    static void success(const std::string& msg)   { log(Level::SUCCESS, msg); }
    static void warning(const std::string& msg)   { log(Level::WARNING, msg); }
    static void error(const std::string& msg)     
    { 
        log(Level::ERROR,   msg); 
        export_to_file(std::string(PROJECT_NAME)+".log");
    }

    // Specialized table row for your Mesh entities
    static void mesh_entity(int dim, int tag, int key_id, int n_elements, const std::string& name);

    // Save all collected logs to a file at once
    static void export_to_file(const std::string& filename);

    // Clear the history if needed
    static void clear() { buffer.str(""); buffer.clear(); }

};