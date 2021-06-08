#include <iomanip>
#include "TyFluid.h"

const int color[4][3] = {{1, 1, 1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

TyF::TyFluid::TyFluid(){
    //dynamic_scene = NULL;
    //scene = NULL;
    use_screen_buffer = false;
    export_screen_buffer = false;
    multisamples = 16;
    enablevsync = true;
    f.K = 2;
    f._dt = 0.01;
    ds = 50.0f;
    stop = false;
    click = false;
    click_mode = false;
    erase_mode = false;
    ((float*)&clear_color)[0] = 1.0;
    ((float*)&clear_color)[1] = 1.0;
    ((float*)&clear_color)[2] = 1.0;
    vector<int> a;
    a.push_back(SCALE);
    a.push_back(SCALE);
    f.initialize(&a);
}

TyF::TyFluid::~TyFluid(){
    //delete dynamic_scene;
    //delete scene;
}

void TyF::TyFluid::drawRect(float x, float y, float w, float h){
    //cout << x << " " << y << " " << h << " " << w << endl;
    glBegin(GL_QUADS);
    glVertex2f(x, y);
    glVertex2f(x + w, y);
    glVertex2f(x + w, y + h);
    glVertex2f(x, y + h);
    glEnd();
}


void TyF::TyFluid::change_density(){
    //cout << SourceX.size() << endl;
    for(int i = 0; i < SourceX.size(); i++){
        int X = SourceX[i];
        int Y = SourceY[i];
        float Vx = SourceVx[i];
        float Vy = SourceVy[i];
        VectorXd c(2);
        for(c(0) = X; c(0) <= X; c(0) ++)
            for(c(1) = Y; c(1) <= Y; c(1) ++){
                int id = f.get_index(c);
                f._density(id) += DENSITY_AMOUNT * ds;
                f._heat(id) += HEAT_AMOUNT;
                f._velocity(id * f.K) += VELOCITY_AMOUNT * Vx;
                f._velocity(id * f.K + 1) += VELOCITY_AMOUNT * Vy;
            }
                
        //cout << f._velocity(id * f.K) << " " << f._velocity(id * f.K + 1) << endl;
    }
}

void TyF::TyFluid::initScreenBuffer(GLFWwindow *window){
    
}

void TyF::TyFluid::windowLoop(GLFWwindow *window){
    //glClear(GL_)
    //printf("flag\n");
    ImGui::NewFrame();
    {
        

        ImGui::Begin("Control");                          // Create a window called "Hello, world!" and append into it.

        ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
        //ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
        //ImGui::Checkbox("Another Window", &show_another_window);

        ImGui::SliderFloat("density", &ds, 0.0f, 100.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
        ImGui::ColorEdit3("smoke color", (float*)&clear_color); // Edit 3 floats representing a color
        if(stop){
            if(ImGui::Button("start"))
                stop ^= true;
        }
        else {
            if(ImGui::Button("stop"))
                stop ^= true;
        }
        if(click_mode)
            ImGui::Text("Click to set a source");
        
        else 
            ImGui::Text("Type C to enter click mode.(set source)");

        if(erase_mode)
            ImGui::Text("click to remove a source");
        else 
            ImGui::Text("Type E to enter erase mode.(erase source)");
        //ImGui::SameLine();
        //ImGui::Text("counter = %d", counter);
        ImGui::SliderFloat("vorticity", &f._vorticityEPS, 0.0f, 10.0f);

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        ImGui::End();
    }
    ImGui::Render();
    int display_w, display_h;
    //glfwMakeContextCurrent(window);
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
    

    if(!stop){
        change_density();
        f.step();
        
    }

    for(int id = 0; id < f.num[f.K]; id++){
        VectorXd c = f.get_cell(id);
        glColor3f(f._density[id] / 32 * ((float*)&clear_color)[0], f._density[id] / 32 * ((float*)&clear_color)[1], f._density[id] / 32 * ((float*)&clear_color)[2]);
        //cout << 1.0 * (c(0) + 1) / display_w * 128 - 0.5 << "______" << endl;
        drawRect(1.0 * (c(0) + 1) / SCALE * 2 - 1, 1.0 * (c(1) + 1) / SCALE * 2 - 1, 2.0 / SCALE, 2.0 / SCALE);
    }
    //cout << f._density << endl << "-------------" << endl;
    
}

void TyF::TyFluid::drawImGui(){
    bool show_demo_window = false;

    //showMainWindow(&fileDialog);


}

void TyF::TyFluid::framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    /*glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, textureColorBufferMultiSampled);
    glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, multisamples, GL_RGB, width, height, GL_TRUE);
    glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo);
    glRenderbufferStorageMultisample(GL_RENDERBUFFER, multisamples, GL_DEPTH24_STENCIL8, width, height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, screenTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glBindTexture(GL_TEXTURE_2D, 0);*/
}