#pragma once
#include "fluid/fluid.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <glad/glad.h>
#include "stb_image_write_o.h"

namespace TyF{

class TyFluid{
public:
    TyFluid();
    ~TyFluid();
    virtual void initScreenBuffer(GLFWwindow *window);
    
    virtual void windowLoop(GLFWwindow *window);
    virtual void drawImGui();
    //virtual void paintGL(/*TyFluidShader *shader*/);
    //virtual void deleteGL();
    //virtual void makeContext();
    //virtual void setCurrentContext();
    void framebuffer_size_callback(GLFWwindow* window, int width, int height);

    bool enablevsync;
    
    vector<int> SourceX, SourceY;
    vector<float> SourceVx, SourceVy;

    bool click;
    Fluid f;
    bool click_mode, erase_mode;

protected:

    //void showMainWindow(ImGui::FileBrowser *fileDialog, bool *p_open = NULL);
    //void saveCurScreenImagePNG(const char* imagename);
    //void saveCurScreenImagePNGAnimation(const char* imagebase);

    //Lobo::Camera camera;

    void drawRect(float x, float y, float h, float w);
    void change_density();

    ImVec4 clear_color;

    std::string config_file_path;

    unsigned int quadVAO;
    unsigned int quadVBO;
    unsigned int framebuffer;
    unsigned int textureColorbuffer;
    unsigned int textureColorBufferMultiSampled;
    unsigned int intermediateFBO;
    unsigned int screenTexture;
    unsigned int rbo;

    float ds;
    bool stop;

    bool use_screen_buffer;
    bool export_screen_buffer;

    int multisamples;
};

}
