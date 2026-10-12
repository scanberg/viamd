// Optional standalone EGL regression test; see docs/density-volume.md.
#define GL_GLEXT_PROTOTYPES
#include <EGL/egl.h>
#include <EGL/eglext.h>
#include <GL/glcorearb.h>
#include <core/md_str.h>
#include "../../src/gfx/volume_overlay_shader.inl"
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

static void check(bool ok, const char* message) {
    if (!ok) { std::fprintf(stderr, "%s\n", message); std::exit(1); }
}
static GLuint shader(GLenum type, const char* text, int length) {
    GLuint id = glCreateShader(type);
    glShaderSource(id, 1, &text, &length);
    glCompileShader(id);
    GLint ok = 0; glGetShaderiv(id, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[8192]; glGetShaderInfoLog(id, sizeof(log), nullptr, log);
        std::fprintf(stderr, "%s\n", log);
    }
    check(ok, "Shader compilation failed");
    return id;
}
int main() {
    EGLDisplay display = eglGetPlatformDisplay(EGL_PLATFORM_SURFACELESS_MESA, EGL_DEFAULT_DISPLAY, nullptr);
    check(eglInitialize(display, nullptr, nullptr), "No surfaceless EGL display");
    check(eglBindAPI(EGL_OPENGL_API), "No EGL OpenGL API");
    EGLint config_attribs[] = {EGL_SURFACE_TYPE,EGL_PBUFFER_BIT,EGL_RENDERABLE_TYPE,EGL_OPENGL_BIT,EGL_NONE};
    EGLConfig config; EGLint count;
    check(eglChooseConfig(display, config_attribs, &config, 1, &count) && count, "No EGL configuration");
    EGLint context_attribs[] = {EGL_CONTEXT_MAJOR_VERSION,4,EGL_CONTEXT_MINOR_VERSION,1,
        EGL_CONTEXT_OPENGL_PROFILE_MASK,EGL_CONTEXT_OPENGL_CORE_PROFILE_BIT,EGL_NONE};
    EGLContext context = eglCreateContext(display, config, EGL_NO_CONTEXT, context_attribs);
    check(context != EGL_NO_CONTEXT && eglMakeCurrent(display,EGL_NO_SURFACE,EGL_NO_SURFACE,context), "No GL 4.1 context");
    const char* vertex = "#version 410 core\nvoid main(){uint i=uint(gl_VertexID);gl_Position=vec4(float(i&1u)*4.-1.,float((i>>1u)&1u)*4.-1.,0,1);}";
    GLuint vs=shader(GL_VERTEX_SHADER, vertex, (int)std::strlen(vertex));
    GLuint fs=shader(GL_FRAGMENT_SHADER, overlay_fragment_source.ptr, (int)overlay_fragment_source.len);
    GLuint program=glCreateProgram();
    glAttachShader(program,vs); glAttachShader(program,fs); glLinkProgram(program);
    GLint linked; glGetProgramiv(program,GL_LINK_STATUS,&linked);
    if (!linked) { char log[8192]; glGetProgramInfoLog(program,sizeof(log),nullptr,log); std::fprintf(stderr,"%s\n",log); }
    check(linked, "Shader link failed");
    glUseProgram(program);
    auto loc=[&](const char* name){ return glGetUniformLocation(program,name); };
    glUniform1i(loc("u_count"),2);
    glUniform2f(loc("u_inv_res"),1,1);
    glUniform3f(loc("u_env"),0.1f,0.1f,0.1f);
    float identity[16]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    float cube[16]={.5f,0,0,0,0,.5f,0,0,0,0,.5f,0,.5f,.5f,.5f,1};
    float normal[9]={1,0,0,0,1,0,0,0,1};
    glUniformMatrix4fv(loc("u_clip_to_view"),1,GL_FALSE,identity);
    GLuint textures[13]; glGenTextures(13,textures);
    float depth=1;
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D,textures[0]);
    glTexImage2D(GL_TEXTURE_2D,0,GL_R32F,1,1,0,GL_RED,GL_FLOAT,&depth);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glUniform1i(loc("u_depth"),0);
    for (int i=0;i<6;++i) {
        char name[64];
        auto at=[&](const char* base) { std::snprintf(name,sizeof(name),"%s[%d]",base,i); return loc(name); };
        glUniform1i(at("u_volume"),1+2*i); glUniform1i(at("u_tf"),2+2*i);
        glUniformMatrix4fv(at("u_clip_to_model"),1,GL_FALSE,cube);
        glUniformMatrix3fv(at("u_normal_to_view"),1,GL_FALSE,normal);
        glUniform3f(at("u_min"),0,0,0); glUniform3f(at("u_max"),1,1,1);
        glUniform2f(at("u_range"),0,1);
        glUniform1i(at("u_dvr"),1); glUniform1i(at("u_iso_count"),0);
        float grid[512];
        for (int z=0;z<8;++z) for(int y=0;y<8;++y) for(int x=0;x<8;++x) grid[z*64+y*8+x]=(z+.5f)/8.f;
        glActiveTexture(GL_TEXTURE0+1+2*i); glBindTexture(GL_TEXTURE_3D,textures[1+2*i]);
        glTexImage3D(GL_TEXTURE_3D,0,GL_R32F,8,8,8,0,GL_RED,GL_FLOAT,grid);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);
        float tf[4]={i==0 ? 1.f : 0.f,0,i==1 ? 1.f : 0.f,.01f};
        glActiveTexture(GL_TEXTURE0+2+2*i); glBindTexture(GL_TEXTURE_2D,textures[2+2*i]);
        glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F,1,1,0,GL_RGBA,GL_FLOAT,tf);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    }
    GLuint output,fbo,vao;
    glGenTextures(1,&output); glBindTexture(GL_TEXTURE_2D,output);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F,1,1,0,GL_RGBA,GL_FLOAT,nullptr);
    glGenFramebuffers(1,&fbo); glBindFramebuffer(GL_FRAMEBUFFER,fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,output,0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    check(glCheckFramebufferStatus(GL_FRAMEBUFFER)==GL_FRAMEBUFFER_COMPLETE,"Incomplete framebuffer");
    glGenVertexArrays(1,&vao); glBindVertexArray(vao); glViewport(0,0,1,1);
    auto draw=[&](float color[4]) {
        glDrawArrays(GL_TRIANGLES,0,3);
        glReadPixels(0,0,1,1,GL_RGBA,GL_FLOAT,color);
        check(glGetError()==GL_NO_ERROR,"GL error in overlay draw");
        for(int i=0;i<4;++i) check(std::isfinite(color[i]),"Nonfinite color");
    };
    auto reorder=[&](bool reversed) {
        glUniform1i(loc("u_tf[0]"),reversed ? 4 : 2);
        glUniform1i(loc("u_tf[1]"),reversed ? 2 : 4);
    };
    float before[4],after[4]; draw(before); reorder(true); draw(after);
    for(int i=0;i<4;++i) check(std::fabs(before[i]-after[i])<1e-5,"DVR depends on layer ordering");
    check(std::fabs(before[0]-before[2])<1e-5 && before[3]>.9f,"Overlapping DVR extinction is incorrect");
    depth=.5f;
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D,textures[0]);
    glTexSubImage2D(GL_TEXTURE_2D,0,0,0,1,1,GL_RED,GL_FLOAT,&depth);
    draw(after);
    check(after[3]<before[3] && after[3]>.5f,"Opaque depth does not truncate the shared ray");
    depth=1;
    glTexSubImage2D(GL_TEXTURE_2D,0,0,0,1,1,GL_RED,GL_FLOAT,&depth);
    float red[4]={1,0,0,1},blue[4]={0,0,1,1};
    glUniform1i(loc("u_dvr[0]"),0); glUniform1i(loc("u_dvr[1]"),0);
    glUniform1i(loc("u_iso_count[0]"),1); glUniform1i(loc("u_iso_count[1]"),1);
    glUniform1f(loc("u_values[0]"),.25f); glUniform4fv(loc("u_colors[0]"),1,red);
    glUniform1f(loc("u_values[8]"),.75f); glUniform4fv(loc("u_colors[8]"),1,blue);
    draw(before);
    glUniform1f(loc("u_values[0]"),.75f); glUniform4fv(loc("u_colors[0]"),1,blue);
    glUniform1f(loc("u_values[8]"),.25f); glUniform4fv(loc("u_colors[8]"),1,red);
    draw(after);
    for(int i=0;i<4;++i) check(std::fabs(before[i]-after[i])<1e-5,"Isosurfaces depend on layer ordering");
    check(before[0]>before[2] && before[3]>.99f,"Far surface hides the near surface");
    glUniform1i(loc("u_dvr[0]"),1); glUniform1i(loc("u_iso_count[0]"),0);
    glUniform1i(loc("u_tf[0]"),4);
    draw(before);
    glUniform1i(loc("u_dvr[0]"),0); glUniform1i(loc("u_iso_count[0]"),1);
    glUniform1f(loc("u_values[0]"),.25f); glUniform4fv(loc("u_colors[0]"),1,red);
    glUniform1i(loc("u_dvr[1]"),1); glUniform1i(loc("u_iso_count[1]"),0);
    glUniform1i(loc("u_tf[1]"),4);
    draw(after);
    for(int i=0;i<4;++i) check(std::fabs(before[i]-after[i])<1e-5,"Mixed DVR/isosurface depends on layer ordering");
    // Shading: an isosurface of a sphere is lit, so its brightness follows its normal
    {
        const int size = 64;
        std::vector<float> depths(size*size, 1.0f);
        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D,textures[0]);
        glTexImage2D(GL_TEXTURE_2D,0,GL_R32F,size,size,0,GL_RED,GL_FLOAT,depths.data());
        glActiveTexture(GL_TEXTURE0+13); glBindTexture(GL_TEXTURE_2D,output);
        glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F,size,size,0,GL_RGBA,GL_FLOAT,nullptr);
        glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,output,0);
        check(glCheckFramebufferStatus(GL_FRAMEBUFFER)==GL_FRAMEBUFFER_COMPLETE,"Incomplete framebuffer after resize");
        glViewport(0,0,size,size);
        glUniform2f(loc("u_inv_res"),1.0f/size,1.0f/size);
        // An orthographic view looking down -z: near 1, far 3
        float clip_to_view[16]={1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,-2,1};
        glUniformMatrix4fv(loc("u_clip_to_view"),1,GL_FALSE,clip_to_view);
        glUniform1i(loc("u_count"),1);
        const int n=32;
        std::vector<float> sphere(n*n*n);
        for (int z=0;z<n;++z) for (int y=0;y<n;++y) for (int x=0;x<n;++x) {
            const float dx=(x+.5f)/n-.5f, dy=(y+.5f)/n-.5f, dz=(z+.5f)/n-.5f;
            sphere[(z*n+y)*n+x]=1.0f-2.0f*std::sqrt(dx*dx+dy*dy+dz*dz);
        }
        glActiveTexture(GL_TEXTURE0+1); glBindTexture(GL_TEXTURE_3D,textures[1]);
        glTexImage3D(GL_TEXTURE_3D,0,GL_R32F,n,n,n,0,GL_RED,GL_FLOAT,sphere.data());
        float white[4]={1,1,1,1};
        glUniform1i(loc("u_dvr[0]"),0); glUniform1i(loc("u_iso_count[0]"),1);
        glUniform1f(loc("u_values[0]"),.5f); glUniform4fv(loc("u_colors[0]"),1,white);
        glClearColor(0,0,0,0); glClear(GL_COLOR_BUFFER_BIT);
        glDrawArrays(GL_TRIANGLES,0,3);
        check(glGetError()==GL_NO_ERROR,"GL error drawing the sphere");
        std::vector<float> image(size*size*4);
        glReadPixels(0,0,size,size,GL_RGBA,GL_FLOAT,image.data());
        float lo=1e30f, hi=0; int covered=0;
        for (int i=0;i<size*size;++i) {
            if (image[i*4+3]<.5f) continue;
            const float luma=(image[i*4]+image[i*4+1]+image[i*4+2])/3.0f;
            lo=std::min(lo,luma); hi=std::max(hi,luma); ++covered;
        }
        std::printf("Sphere: %d pixels, luminance %.4f to %.4f\n",covered,lo,hi);
        check(covered>100,"The sphere is not drawn");
        check(hi>1.5f*lo+1e-3f,"The isosurface is not shaded: its brightness does not follow its normal");
    }
    std::puts("Overlay shader: compilation, overlapping DVR, depth truncation, surface ordering and mixed rendering passed");
    glDeleteProgram(program); glDeleteShader(vs); glDeleteShader(fs);
    glDeleteTextures(13,textures); glDeleteTextures(1,&output);
    glDeleteFramebuffers(1,&fbo); glDeleteVertexArrays(1,&vao);
    eglMakeCurrent(display,EGL_NO_SURFACE,EGL_NO_SURFACE,EGL_NO_CONTEXT);
    eglDestroyContext(display,context); eglTerminate(display);
}
