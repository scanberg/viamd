#include <app/imgui_impl_md_gpu.h>

#include <md_post.h>
#include <core/md_log.h>

#include <string.h>

struct ImGui_ImplMdGpu_Data {
    md_gpu_device_t  device = nullptr;
    md_gpu_stream_t  stream = nullptr;
    md_post_ctx_t    post = nullptr;
    md_gpu_sampler_t linear = {};       // the default: between mips too
    md_gpu_sampler_t nearest = {};
    ImVector<md_post_ui_list_t> lists;
    ImVector<md_post_ui_cmd_t>  cmds;
    ImVector<unsigned char>     pixels; // an update's rows, packed
};

static ImGui_ImplMdGpu_Data* backend_data() {
    return ImGui::GetCurrentContext() ? (ImGui_ImplMdGpu_Data*)ImGui::GetIO().BackendRendererUserData : nullptr;
}

void ImGui_ImplMdGpu_SetState(const ImDrawList*, const ImDrawCmd*) {
    // Only marks the commands after it: RenderDrawData reads its data
}

static void destroy_texture(ImTextureData* tex) {
    if (tex->BackendUserData) md_gpu_texture_destroy((md_gpu_texture_t)tex->BackendUserData);   // deferred past the GPU's use
    tex->BackendUserData = nullptr;
    tex->SetTexID(ImTextureID_Invalid);
    tex->SetStatus(ImTextureStatus_Destroyed);
}

// The texture's region (x, y, w, h) from ImGui's pixels, its rows packed
static bool upload_region(ImGui_ImplMdGpu_Data* bd, ImTextureData* tex, int x, int y, int w, int h) {
    const int row = w * tex->BytesPerPixel;
    bd->pixels.resize(row * h);
    for (int r = 0; r < h; ++r) memcpy(bd->pixels.Data + (size_t)r * row, tex->GetPixelsAt(x, y + r), (size_t)row);
    md_gpu_tex_region_t rg = {};
    rg.offset[0] = (uint32_t)x;
    rg.offset[1] = (uint32_t)y;
    rg.extent[0] = (uint32_t)w;
    rg.extent[1] = (uint32_t)h;
    rg.extent[2] = 1;
    return md_gpu_upload_texture(bd->stream, (md_gpu_texture_t)tex->BackendUserData, &rg, bd->pixels.Data, (size_t)row * h);
}

static void update_texture(ImGui_ImplMdGpu_Data* bd, ImTextureData* tex) {
    if (tex->Status == ImTextureStatus_WantCreate) {
        IM_ASSERT(tex->TexID == ImTextureID_Invalid && tex->BackendUserData == nullptr);
        IM_ASSERT(tex->Format == ImTextureFormat_RGBA32);
        md_gpu_texture_desc_t d = {};
        d.type   = MD_GPU_TEX_2D;
        d.format = MD_GPU_FORMAT_RGBA8_UNORM;
        d.usage  = MD_GPU_TEX_SAMPLED;
        d.width  = (uint32_t)tex->Width;
        d.height = (uint32_t)tex->Height;
        d.label  = "imgui";
        md_gpu_texture_t t = md_gpu_texture_create(bd->stream, &d);
        if (!t) {
            MD_LOG_ERROR("ImGui's texture (%d x %d) could not be made: %s", tex->Width, tex->Height, md_gpu_last_error());
            return;
        }
        tex->BackendUserData = t;
        if (!upload_region(bd, tex, 0, 0, tex->Width, tex->Height)) {
            MD_LOG_ERROR("ImGui's texture could not be uploaded: %s", md_gpu_last_error());
        }
        tex->SetTexID(ImGui_ImplMdGpu_TextureID(t));
        tex->SetStatus(ImTextureStatus_OK);
    } else if (tex->Status == ImTextureStatus_WantUpdates) {
        for (const ImTextureRect& r : tex->Updates) {
            if (!upload_region(bd, tex, r.x, r.y, r.w, r.h)) {
                MD_LOG_ERROR("ImGui's texture could not be updated: %s", md_gpu_last_error());
                break;
            }
        }
        tex->SetStatus(ImTextureStatus_OK);
    } else if (tex->Status == ImTextureStatus_WantDestroy && tex->UnusedFrames > 0) {
        destroy_texture(tex);
    }
}

bool ImGui_ImplMdGpu_Init(md_gpu_device_t device, md_gpu_stream_t stream) {
    ImGuiIO& io = ImGui::GetIO();
    IM_ASSERT(io.BackendRendererUserData == nullptr && "Already initialized a renderer backend!");
    ImGui_ImplMdGpu_Data* bd = IM_NEW(ImGui_ImplMdGpu_Data)();
    bd->device = device;
    bd->stream = stream;
    bd->post   = md_post_create(device);
    md_gpu_sampler_desc_t sd = {};
    sd.min_filter = sd.mag_filter = sd.mip_filter = MD_GPU_FILTER_LINEAR;
    bd->linear = md_gpu_sampler(device, &sd);
    sd.min_filter = sd.mag_filter = sd.mip_filter = MD_GPU_FILTER_NEAREST;
    bd->nearest = md_gpu_sampler(device, &sd);
    if (!bd->post || !bd->linear.handle || !bd->nearest.handle) {
        MD_LOG_ERROR("ImGui's renderer could not be made: %s", md_gpu_last_error());
        if (bd->post) md_post_destroy(bd->post);
        IM_DELETE(bd);
        return false;
    }
    io.BackendRendererUserData = bd;
    io.BackendRendererName = "imgui_impl_md_gpu";
    io.BackendFlags |= ImGuiBackendFlags_RendererHasVtxOffset | ImGuiBackendFlags_RendererHasTextures;
    return true;
}

void ImGui_ImplMdGpu_Shutdown() {
    ImGui_ImplMdGpu_Data* bd = backend_data();
    if (!bd) return;
    for (ImTextureData* tex : ImGui::GetPlatformIO().Textures) {
        if (tex->RefCount == 1) destroy_texture(tex);
    }
    md_post_destroy(bd->post);
    ImGuiIO& io = ImGui::GetIO();
    io.BackendRendererName = nullptr;
    io.BackendRendererUserData = nullptr;
    io.BackendFlags &= ~(ImGuiBackendFlags_RendererHasVtxOffset | ImGuiBackendFlags_RendererHasTextures);
    IM_DELETE(bd);
}

bool ImGui_ImplMdGpu_RenderDrawData(ImDrawData* draw_data, md_gpu_texture_t target, const float* clear) {
    ImGui_ImplMdGpu_Data* bd = backend_data();
    if (!bd || !draw_data || !target) return false;
    if (draw_data->Textures) {
        for (ImTextureData* tex : *draw_data->Textures) {
            if (tex->Status != ImTextureStatus_OK) update_texture(bd, tex);
        }
    }

    // The lists' commands, each with the state its callbacks left
    bd->lists.resize(0);
    bd->cmds.resize(0);
    for (const ImDrawList* dl : draw_data->CmdLists) {
        uint32_t state = 0;
        for (const ImDrawCmd& c : dl->CmdBuffer) {
            if (c.UserCallback) {
                if (c.UserCallback == ImDrawCallback_ResetRenderState) state = 0;
                else if (c.UserCallback == ImGui_ImplMdGpu_SetState) state = (uint32_t)(uintptr_t)c.UserCallbackData;
                else c.UserCallback(dl, &c);
                continue;
            }
            md_post_ui_cmd_t m = {};
            m.texture.handle = (uint64_t)c.GetTexID();
            m.sampler = (state & ImGui_ImplMdGpu_Nearest) ? bd->nearest : bd->linear;
            m.clip[0] = c.ClipRect.x;
            m.clip[1] = c.ClipRect.y;
            m.clip[2] = c.ClipRect.z;
            m.clip[3] = c.ClipRect.w;
            m.first_index   = c.IdxOffset;
            m.index_count   = c.ElemCount;
            m.vertex_offset = c.VtxOffset;
            m.flags = (state & ImGui_ImplMdGpu_Premultiplied) ? (uint32_t)MD_POST_UI_PREMULTIPLIED : 0u;
            bd->cmds.push_back(m);
        }
        md_post_ui_list_t l = {};
        l.vertices     = (const md_post_ui_vertex_t*)dl->VtxBuffer.Data;
        l.vertex_count = (uint32_t)dl->VtxBuffer.Size;
        l.indices      = dl->IdxBuffer.Data;
        l.index_count  = (uint32_t)dl->IdxBuffer.Size;
        l.cmd_count    = (uint32_t)bd->cmds.Size;   // the end, for now: made counts below
        bd->lists.push_back(l);
    }
    // The commands' pointers, once the array has stopped growing
    uint32_t beg = 0;
    for (md_post_ui_list_t& l : bd->lists) {
        const uint32_t end = l.cmd_count;
        l.cmds = bd->cmds.Data + beg;
        l.cmd_count = end - beg;
        beg = end;
    }

    md_post_ui_t ui = {};
    ui.target = target;
    if (clear) {
        ui.clear = true;
        for (int k = 0; k < 4; ++k) ui.clear_color[k] = clear[k];
    }
    ui.origin[0]  = draw_data->DisplayPos.x;
    ui.origin[1]  = draw_data->DisplayPos.y;
    ui.scale[0]   = draw_data->FramebufferScale.x;
    ui.scale[1]   = draw_data->FramebufferScale.y;
    ui.index_type = sizeof(ImDrawIdx) == 2 ? MD_GPU_INDEX_U16 : MD_GPU_INDEX_U32;
    ui.lists      = bd->lists.Data;
    ui.list_count = (uint32_t)bd->lists.Size;
    return md_post_ui(bd->stream, bd->post, &ui);
}

static_assert(sizeof(ImDrawVert) == sizeof(md_post_ui_vertex_t), "ImDrawVert as md_post_ui's vertex");
static_assert(sizeof(ImTextureID) == sizeof(uint64_t), "an ImTextureID holds a sampled handle");
