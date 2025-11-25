#include "visualizer/app_window.h"

void APP_Window::run(){
    initWindow();
    initVulkan();
    mainLoop();
    cleanup();
}

void APP_Window::initWindow() {
    // initialize the GLFW library
    glfwInit();

    // GLFW was originally designed to create an OpenGL context,
    // tell it to not create an openGL context with a subsequent call:
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);

    // disable resized window for now with another window hint call:
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

    // create the actual window:
    //      the 4th parameter allows you to optionally specify a monitor to open the window on
    //      the last parameter is only relevant to OpenGL.
    window = glfwCreateWindow(WIDTH, HEIGHT, "simu-PHY", nullptr, nullptr);
}

void APP_Window::initVulkan() {
    // The instance is the connection between your application and the Vulkan library, 
    // creating it involves specifying some details about your application to the driver.
    createInstance();
}

void APP_Window::mainLoop() {
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
    }
}

void APP_Window::cleanup() {
    // the allocation and deallocation functions in Vulkan have an optional allocator callback,
    // for now we'll ignore by passing nullptr.
    vkDestroyInstance(instance, nullptr);
    glfwDestroyWindow(window);
    glfwTerminate();
}

void APP_Window::createInstance() {
    VkApplicationInfo appInfo{};
    appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    appInfo.pApplicationName = "Triangle test";
    appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.pEngineName = "No Engine";
    appInfo.engineVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.apiVersion = VK_API_VERSION_1_0;

    // tells the Vulkan driver which global extensions and validation layers we want to use.
    VkInstanceCreateInfo createInfo{};
    // structs in Vulkan require you to explicitly specify the type in the sType member.
    createInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    createInfo.pApplicationInfo = &appInfo;

    uint32_t glfwExtensionCount = 0;
    const char** glfwExtensions;

    glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);

    checkExtension(glfwExtensions, glfwExtensionCount);

    createInfo.enabledExtensionCount = glfwExtensionCount;
    createInfo.ppEnabledExtensionNames = glfwExtensions;

    createInfo.enabledLayerCount = 0;

    // the general pattern that object creation function parameters in Vulkan follow is:
    //      - Pointer to struct with creation info
    //      - Pointer to custom allocator callbacks
    //      - Pointer to the variable that stores the handle to the new object
    if(vkCreateInstance(&createInfo, nullptr, &instance)!= VK_SUCCESS) {
        throw std::runtime_error("failed to create instance!");
    }

}

void APP_Window::checkExtension(const char** glfwExtensions, int glfwExtensionCount) {
    // To allocate an array to hold the extension details we first need to know how many there are.
    // one can request just the number of extensions by leaving the latter parameter empty:
    uint32_t extensionCount = 0;
    vkEnumerateInstanceExtensionProperties(nullptr, &extensionCount, nullptr);

    // Now allocate an array to hold the extension details
    std::vector<VkExtensionProperties> extensions(extensionCount);
    vkEnumerateInstanceExtensionProperties(nullptr, &extensionCount, extensions.data());

    std::cout << "required extensions:\n";
    for(int i=0; i<glfwExtensionCount; ++i){
        for(const auto& extension: extensions){
            if(std::string(*(glfwExtensions+i))==std::string(extension.extensionName)){
                goto found;
            }
        }
        std::cout << '\t' << *(glfwExtensions+i) << " - note found!" << '\n';
        continue;

        found:
        std::cout << '\t' << *(glfwExtensions+i) << " - found" << '\n';
    }

    std::cout << "available extensions:\n";

    for(const auto& extension : extensions){
        std::cout << '\t' << extension.extensionName << '\n';
    }
}

