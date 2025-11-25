#ifndef DISPLAY_H
#define DISPLAY_H

// GLFW will include its own definitions and automatically load the Vulkan header with it.
#include <vulkan/vulkan.h>

#define GLFW_INCLUDE_VULKAN

#include <GLFW/glfw3.h>

#include <iostream>
#include <stdexcept>
#include <cstdlib>

#include <vector>
#include <algorithm>

class APP_Window {

    public:
        void run();


    private:

        VkInstance instance;

        GLFWwindow* window;

        const uint32_t WIDTH = 800;
        const uint32_t HEIGHT = 600;

        void initWindow();
        void initVulkan();
        void mainLoop();
        void cleanup();
        void createInstance();
        void checkExtension(const char** glfwExtensions, int glfwExtensionCount);



};

#endif 