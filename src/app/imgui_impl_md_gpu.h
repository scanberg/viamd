#pragma once

// Dear ImGui's renderer on md_gpu: the draw data drawn by md_post_ui into a texture (a surface's image, or one of the
// caller's), the textures ImGui makes and updates (its font atlas: ImGuiBackendFlags_RendererHasTextures) as md_gpu
// textures. An ImTextureID is an md_gpu texture's sampled handle (ImGui_ImplMdGpu_TextureID), so what the application
// draws with md_gpu goes into the UI as it is.
//
// Drawing state, per command list, set by a callback and back to the default at ImDrawCallback_ResetRenderState:
//     ImGui_ImplMdGpu_PushState(dl, ImGui_ImplMdGpu_Premultiplied);
//     dl->AddImage(...);
//     ImGui_ImplMdGpu_PopState(dl);
// Default: straight alpha, sampled linearly (between mips too, where the texture has them), clamped to the edge.

#include <imgui.h>
#include <core/md_gpu.h>

#include <stdint.h>

enum {
    ImGui_ImplMdGpu_Premultiplied = 1u << 0,    // the colour times the sample is premultiplied by its alpha
    ImGui_ImplMdGpu_Nearest       = 1u << 1,    // sampled at the nearest texel
};

// The stream draws and uploads: a GRAPHICS stream, the one the surface is acquired on
bool ImGui_ImplMdGpu_Init(md_gpu_device_t device, md_gpu_stream_t stream);
void ImGui_ImplMdGpu_Shutdown();

// The draw data into `target` (RENDER_TARGET, RGBA8_UNORM or BGRA8_UNORM, of the draw data's size in pixels), over
// what is there, or cleared to `clear` first (RGBA, given as a non-null pointer). Can be called again with the same
// draw data for another target (a screenshot).
bool ImGui_ImplMdGpu_RenderDrawData(ImDrawData* draw_data, md_gpu_texture_t target, const float* clear = nullptr);

// An md_gpu texture as ImGui draws it: 2D, SAMPLED, of a float format
static inline ImTextureID ImGui_ImplMdGpu_TextureID(md_gpu_texture_t tex) {
    return tex ? (ImTextureID)md_gpu_texture_sampled(tex).handle : ImTextureID_Invalid;
}

// The callback that sets the state (cmd->UserCallbackData: ImGui_ImplMdGpu_* flags)
void ImGui_ImplMdGpu_SetState(const ImDrawList* list, const ImDrawCmd* cmd);

static inline void ImGui_ImplMdGpu_PushState(ImDrawList* dl, uint32_t state) {
    dl->AddCallback(ImGui_ImplMdGpu_SetState, (void*)(uintptr_t)state);
}
static inline void ImGui_ImplMdGpu_PopState(ImDrawList* dl) {
    dl->AddCallback(ImDrawCallback_ResetRenderState, nullptr);
}
