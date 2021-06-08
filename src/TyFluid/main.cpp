#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <omp.h>
#include "TyFluid.h"
#include <cstdio>

#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
#include <GL/gl3w.h> // Initialize with gl3wInit()
#pragma message("C Preprocessor got here!")
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
#include <GL/glew.h> // Initialize with glewInit()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
#include <glad/glad.h> // Initialize with gladLoadGL()
#else
#include IMGUI_IMPL_OPENGL_LOADER_CUSTOM
#endif

#include <GLFW/glfw3.h>
#include <glad/glad.h>

#include <Eigen/Dense>
#include <iostream>

TyF::TyFluid *fem_main_p = NULL;

static void glfw_error_callback(int error, const char *description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

void framebuffer_size_callback(GLFWwindow *window, int width, int height){
    if (fem_main_p != NULL)
        fem_main_p->framebuffer_size_callback(window, width, height);
}

double lastX, lastY;

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods){
    if(action != GLFW_PRESS) return ;
    if(fem_main_p->erase_mode){
        double nowX, nowY;
        glfwGetCursorPos(window, &nowX, &nowY);
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        for(int i = 0; i < fem_main_p->SourceX.size(); i++){
            double X = 1.0 * fem_main_p->SourceX[i] / SCALE * display_w;
            double Y = 1.0 * (SCALE - fem_main_p->SourceY[i]) / SCALE * display_h;
            if(fabs(nowX - X) <= 15 && fabs(nowY - Y) <= 15){
                for(int j = i; j < fem_main_p->SourceX.size() - 1; j++){
                    fem_main_p->SourceX[j] = fem_main_p->SourceX[j + 1];
                    fem_main_p->SourceY[j] = fem_main_p->SourceY[j + 1];
                    fem_main_p->SourceVx[j] = fem_main_p->SourceVx[j + 1];
                    fem_main_p->SourceVy[j] = fem_main_p->SourceVy[j + 1];
                }
                fem_main_p->SourceX.pop_back();
                fem_main_p->SourceY.pop_back();
                fem_main_p->SourceVx.pop_back();
                fem_main_p->SourceVy.pop_back();
                i--;
            }
        }
    }
    if(!fem_main_p->click_mode) return;
    if(!fem_main_p->click){
        fem_main_p->click = true;
        glfwGetCursorPos(window, &lastX, &lastY);
    }
    else {
        double nowX, nowY;
        glfwGetCursorPos(window, &nowX, &nowY);
        float diffX = nowX - lastX;
        float diffY = lastY - nowY;
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        int X = lastX / display_w * SCALE;
        int Y = SCALE - lastY / display_h * SCALE;
        cout << X << " " << Y << endl;
        fem_main_p->click = false;
        fem_main_p->SourceX.push_back(X);
        fem_main_p->SourceY.push_back(Y);
        fem_main_p->SourceVx.push_back(diffX / SCALE);
        fem_main_p->SourceVy.push_back(diffY / SCALE);
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
    if(action != GLFW_PRESS) return ;
    switch(key){
        case GLFW_KEY_P:
            for(int i = 0; i < fem_main_p->f.num[fem_main_p->f.K] * 2; i++)
                if(fabs(fem_main_p->f._velocity(i)) > 1e-1)
                    cout << i << " " << fem_main_p->f._velocity(i) << endl;
            break;
        case GLFW_KEY_E:
            fem_main_p->erase_mode ^= 1;
            if(fem_main_p->erase_mode){
                fem_main_p->click_mode = false;
                fem_main_p->click = false;
            }
            break;

        case GLFW_KEY_C:
            fem_main_p->click_mode ^= 1;
            fem_main_p->click = false;
            if(fem_main_p->click_mode)
                fem_main_p->erase_mode = false;
            break;
        
    }
}

int main(int argc, char *argv[]){
    printf("flag\n");
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(
        4); // Use 4 threads for all consecutive parallel regions
    
    // Setup window
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    const char *glsl_version = "#version 330";
    const GLFWvidmode *mode = glfwGetVideoMode(glfwGetPrimaryMonitor());

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_SAMPLES, 16);

    glfwWindowHint(GLFW_RED_BITS, mode->redBits);
    glfwWindowHint(GLFW_GREEN_BITS, mode->greenBits);
    glfwWindowHint(GLFW_BLUE_BITS, mode->blueBits);
    glfwWindowHint(GLFW_REFRESH_RATE, mode->refreshRate);
    int window_w = mode->width;
    int window_h = mode->height;
    window_w = 1024;
    window_h = 1024;
    // Create window with graphics context
    GLFWwindow *window = glfwCreateWindow(
        window_w, window_h, "New TyFluid", NULL, NULL);
    if (window == NULL)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetKeyCallback(window, key_callback);

    glfwSwapInterval(1); // Enable vsync

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // Initialize OpenGL loader
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
    bool err = gl3wInit() != 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
    bool err = glewInit() != GLEW_OK;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
    bool err = gladLoadGL() == 0;
#else
    bool err =
        false; // If you use IMGUI_IMPL_OPENGL_LOADER_CUSTOM, your loader is
               // likely to requires some form of initialization.
#endif
    if (err)
    {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        return 1;
    }

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    ImGuiStyle &style = ImGui::GetStyle();
    style.Colors[2].w = 0.8;
    // ImGui::StyleColorsClassic();

    // Setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    bool show_demo_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
    //clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);

    TyF::TyFluid fem_main;
    //fem_main.makeContext();
    //fem_main.setCurrentContext();
    fem_main.initScreenBuffer(window); // init screen buffer objects
    
    /*if(xml_file_path!=NULL)
    fem_main.loadXMLfile(xml_file_path);*/


    fem_main_p = &fem_main;
    //fem_main.initialGL();

    //important
    // Main loop
    while (!glfwWindowShouldClose(window))
    {
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to
        // tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data
        // to your main application.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input
        // data to your main application. Generally you may always pass all
        // inputs to dear imgui, and hide them from your application based on
        // those two flags.
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();

        fem_main.windowLoop(window);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());


        glfwSwapInterval(fem_main.enablevsync); // Enable vsync

        glfwSwapBuffers(window);
    }

    //fem_main.deleteGL();

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
