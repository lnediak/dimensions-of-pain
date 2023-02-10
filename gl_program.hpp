#ifndef GL_PROGRAM_HPP_
#define GL_PROGRAM_HPP_

#include "glad/gl.h"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct GLProgram {

  GLuint prog;
  GLuint posLoc;
  GLuint colLoc;
  GLuint projMatLoc;

  GLProgram() : prog(0) {}
  GLProgram(const GLProgram &) = delete;
  GLProgram(GLProgram &&o)
      : prog(o.prog), posLoc(o.posLoc), colLoc(o.colLoc),
        projMatLoc(o.projMatLoc) {
    o.prog = 0;
  }
  GLProgram &operator=(const GLProgram &) = delete;
  GLProgram &operator=(GLProgram &&other) {
    if (prog != other.prog) {
      this->~GLProgram();
      prog = other.prog;
      posLoc = other.posLoc;
      colLoc = other.colLoc;
      projMatLoc = other.projMatLoc;
    }
    other.prog = 0;
    return *this;
  }
  ~GLProgram() {
    if (prog) {
      glDeleteProgram(prog);
    }
  }

  void compileProgram() {
    std::string vsrc = "#version 130\n"
                       "in vec3 pos;"
                       "in vec4 col;"
                       "out vec4 outCol;"
                       "uniform mat4 projMat;"
                       "void main() {"
                       "  gl_Position = projMat * vec4(pos, 1.0);"
                       "  outCol = col;"
                       "}";
    std::string fsrc = "#version 130\n"
                       "in vec4 outCol;"
                       "out vec4 fragCol;"
                       "void main() {"
                       "  fragCol = outCol;"
                       "}";
    GLuint vs = createShader(vsrc, GL_VERTEX_SHADER);
    GLuint fs = createShader(fsrc, GL_FRAGMENT_SHADER);

    GLuint pro = glCreateProgram();
    if (!pro) {
      throw std::runtime_error("Failed to create OpenGL program.");
    }
    glAttachShader(pro, vs);
    glAttachShader(pro, fs);
    glLinkProgram(pro);
    glDetachShader(pro, vs);
    glDetachShader(pro, fs);
    glDeleteShader(vs);
    glDeleteShader(fs);

    GLint llen = 0;
    glGetProgramiv(pro, GL_INFO_LOG_LENGTH, &llen);
    if (llen) {
      std::string infoLog(llen, ' ');
      glGetProgramInfoLog(pro, llen, NULL, &infoLog[0]);
      std::cerr << infoLog << std::endl;
    }

    GLint success = 0;
    glGetProgramiv(pro, GL_LINK_STATUS, &success);
    if (success == GL_FALSE) {
      glDeleteProgram(pro);
      throw std::runtime_error("Failed to link program.");
    }
    if (prog) {
      glDeleteProgram(prog);
    }
    prog = pro;
    glUseProgram(prog);

    posLoc = glGetAttribLocation(prog, "pos");
    colLoc = glGetAttribLocation(prog, "col");
    projMatLoc = glGetUniformLocation(prog, "projMat");
  }

  /// row-major please
  void setProjMat(const double *p) {
    GLfloat mat[16];
    for (int i = 0; i < 16; i++) {
      mat[i] = p[i];
    }
    glUniformMatrix4fv(projMatLoc, 1, GL_TRUE, mat);
  }

private:
  GLuint createShader(const std::string &src, GLenum stype) const {
    const char *csrc = src.c_str();
    const GLint slen = src.length();

    GLuint shad = glCreateShader(stype);
    if (!shad) {
      throw std::runtime_error("Failed to create OpenGL shader.");
    }
    glShaderSource(shad, 1, &csrc, &slen);
    glCompileShader(shad);

    GLint llen = 0;
    glGetShaderiv(shad, GL_INFO_LOG_LENGTH, &llen);
    if (llen) {
      std::string infoLog(llen, ' ');
      glGetShaderInfoLog(shad, llen, NULL, &infoLog[0]);
      std::cerr << infoLog << std::endl;
    }

    GLint success = 0;
    glGetShaderiv(shad, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
      glDeleteShader(shad);
      throw std::runtime_error(std::string() +
                               "Failed to compile shader. Source:\n" + src);
    }
    return shad;
  }
};

#define BYTES_PER_VERT (4 * sizeof(float))

struct GLMesh {

  GLuint vao, buf;
  std::size_t sz, nverts;

  void initData(std::size_t nsz, const void *data) {
    sz = nsz;
    nverts = nsz / BYTES_PER_VERT;
    glGenBuffers(1, &buf);
    glBindBuffer(GL_ARRAY_BUFFER, buf);
    glBufferData(GL_ARRAY_BUFFER, nsz, data, GL_DYNAMIC_DRAW);
  }

  GLMesh() : vao(0), buf(0), sz(0), nverts(0) {}
  GLMesh(std::size_t sz, const void *data) {
    glGenVertexArrays(1, &vao);
    initData(sz, data);
  }
  GLMesh(const GLMesh &) = delete;
  GLMesh(GLMesh &&o) : vao(o.vao), buf(o.buf), sz(o.sz), nverts(o.nverts) {
    o.vao = o.buf = o.sz = o.nverts = 0;
  }
  GLMesh &operator=(const GLMesh &) = delete;
  GLMesh &operator=(GLMesh &&o) {
    if (buf != o.buf) {
      this->~GLMesh();
      vao = o.vao;
      buf = o.buf;
      sz = o.sz;
      nverts = o.nverts;
    }
    o.vao = o.buf = o.sz = o.nverts = 0;
    return *this;
  }

  operator bool() const { return !!buf; }

  void updateData(std::size_t nsz, const void *data) {
    if (!nsz) {
      if (buf) {
        glDeleteBuffers(1, &buf);
        glDeleteVertexArrays(1, &vao);
        buf = vao = sz = nverts = 0;
      }
      return;
    }
    if (!buf) {
      glGenVertexArrays(1, &vao);
      initData(nsz, data);
      return;
    }
    nverts = nsz / BYTES_PER_VERT;
    glBindBuffer(GL_ARRAY_BUFFER, buf);
    if (nsz > sz) {
      glBufferData(GL_ARRAY_BUFFER, nsz, data, GL_DYNAMIC_DRAW);
    } else {
      glBufferSubData(GL_ARRAY_BUFFER, 0, nsz, data);
    }
  }

  ~GLMesh() {
    if (buf) {
      glDeleteBuffers(1, &buf);
      glDeleteVertexArrays(1, &vao);
      buf = vao = sz = nverts = 0;
    }
  }
};

struct GLTrianglesRecorder {

  std::vector<GLMesh *> meshes;
  GLProgram &prog;

  GLTrianglesRecorder(GLProgram &prog) : prog(prog) {}

  void clear() { meshes.clear(); }

  void operator()(GLMesh &tag, std::size_t sz, const float *triangles) {
    if (!sz) {
      tag.~GLMesh();
      new (&tag) GLMesh();
      return;
    }
    if (!tag) {
      tag.~GLMesh();
      new (&tag) GLMesh(sz * sizeof(float), triangles);
    } else {
      tag.updateData(sz * sizeof(float), triangles);
    }
    meshes.push_back(&tag);
    glBindVertexArray(tag.vao);
    glBindBuffer(GL_ARRAY_BUFFER, tag.buf);
    glVertexAttribPointer(prog.posLoc, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float),
                          0);
    glVertexAttribPointer(prog.colLoc, 4, GL_UNSIGNED_BYTE, GL_TRUE,
                          4 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(prog.posLoc);
    glEnableVertexAttribArray(prog.colLoc);
  }
  void operator()(GLMesh &tag) {
    if (tag) {
      meshes.push_back(&tag);
    }
  }

  void setProjMat(const double *mat) { prog.setProjMat(mat); }
  void setDefaultProjMat(double rm, double um, double fm) {
    double mat[][4] = {{1. / rm, 0, 0, 0},
                       {0, 1. / um, 0, 0},
                       {0, 0, (fm + 1) / (fm - 1), -2 * fm / (fm - 1)},
                       {0, 0, 1, 0}};
    setProjMat(&mat[0][0]);
  }

  void renderAll() const {
    for (GLMesh *mesh : meshes) {
      glBindVertexArray(mesh->vao);
      glDrawArrays(GL_TRIANGLES, 0, mesh->nverts);
    }
  }
};

#endif // GL_PROGRAM_HPP_
