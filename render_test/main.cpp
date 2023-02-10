#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <chrono>
#include <iostream>

#define GLAD_GL_IMPLEMENTATION
#include "gl_program.hpp"

#include "convex.hpp"
#include "util.hpp"

namespace {

void errCallback(int, const char *err) { std::cerr << err << std::endl; }

} // namespace

int main() {
  glfwSetErrorCallback(&errCallback);
  if (!glfwInit()) {
    std::cout << "GLFW Initialization Failed." << std::endl;
    return 1;
  }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  const std::size_t width = 1000, height = 1000;
  GLFWwindow *window =
      glfwCreateWindow(width, height, "test lel", nullptr, nullptr);
  if (!window) {
    std::cout << "GLFW window creation failed." << std::endl;
    return 1;
  }
  glfwMakeContextCurrent(window);
  if (!gladLoadGL((GLADloadfunc)glfwGetProcAddress)) {
    std::cout << "GLAD loading failed." << std::endl;
  }
  glClearColor(0.f, 0.f, 0.f, 1.f);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glFrontFace(GL_CCW);

  SliceDirs<3> sd = {{3, 0, 3}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, 1.5, 1.5, 128};
  v::DVec<3> e1 = {1, 0, 0}, e2 = {0, 1, 0}, e3 = {0, 0, 1};
  v::DVec<3> e123 = e1 + e2 + e3;
  e123 /= std::sqrt(v::norm2(e123));
  Color red = {{255, 0, 0, 255}};
  Color green = {{0, 255, 0, 255}};
  Color blue = {{0, 0, 255, 255}};
  Color white = {{255, 255, 255, 255}};
  Polytope<3, Color> poly{
      {e123, 4, white}, {-e1, -2, red}, {-e2, -2, green}, {-e3, -2, blue}};
  std::vector<float> triangles;

  GLProgram prog;
  prog.compileProgram();
  GLTrianglesRecorder fun(prog);
  fun.setDefaultProjMat(sd.rm, sd.um, sd.fm);
  GLMesh mesh;

  std::size_t fpsCount = 0;
  auto beg = std::chrono::high_resolution_clock::now();
  while (!glfwWindowShouldClose(window)) {
    glfwSwapBuffers(window);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glfwPollEvents();

    // std::cout << "NEW FRAME\n\n\n\n\n";

    triangles.clear();
    poly.writeTriangles(sd,
                        [&triangles](const v::DVec<3> &a, const v::DVec<3> &b,
                                     const v::DVec<3> &c, Color col) -> void {
                          for (const v::DVec<3> &vec : {a, b, c}) {
                            triangles.push_back(vec[0]);
                            triangles.push_back(vec[1]);
                            triangles.push_back(vec[2]);
                            triangles.push_back(col.f);
                            // std::cout << vec << std::endl;
                          }
                        });
    fun(mesh, triangles.size(), &triangles[0]);
    fun.renderAll();
    fun.clear();

    GLint err = glGetError();
    if (err) {
      std::cout << "OpenGL Error: " << err << std::endl;
    }

    sd.c += 0.005;
    // sd.c[2] -= 0.01;

    auto dur = std::chrono::high_resolution_clock::now() - beg;
    double secs =
        std::chrono::duration_cast<std::chrono::duration<double>>(dur).count();
    fpsCount++;
    if (secs >= 1) {
      std::cout << fpsCount << "fps" << std::endl;
      fpsCount = 0;
      beg = std::chrono::high_resolution_clock::now();
    }
  }
}
