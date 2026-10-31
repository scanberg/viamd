// The script interface (script_engine.h) for v1: mdlib's md_script and md_filter, call for call. What v1 hands out is
// converted to the interface's types where they are viamd's own, and passed through where only the name differs.

#include "script_engine_impl.h"

#include <md_script.h>
#include <md_filter.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_array.h>

#include <string.h>

namespace script {
namespace v1 {

static inline Range range_of(md_script_range_marker_t r) { return Range{r.beg, r.end}; }
static inline Ref ref_of(md_script_vis_ref_t r) { return Ref{r.ir_id, r.node_idx}; }
static inline md_script_vis_ref_t v1_ref(Ref r) { md_script_vis_ref_t v = {}; v.ir_id = r.program; v.node_idx = r.node; return v; }
static inline const md_script_ir_t* v1_ir(const Program* p) { return p ? p->v1 : nullptr; }

// ### COMPILE ###

Program* compile(Language lang, str_t src, const CompileDesc& desc, md_allocator_i* alloc) {
    Program* p = (Program*)md_alloc(alloc, sizeof(Program));
    memset(p, 0, sizeof(*p));
    p->lang = lang;
    p->alloc = alloc;
    p->v1 = md_script_ir_create(alloc);
    for (size_t i = 0; i < desc.num_sets; ++i) {
        md_script_ir_add_identifier_bitset(p->v1, desc.sets[i].name, desc.sets[i].atoms);
    }
    // Empty source is not compiled (v1 would log that it is empty): the program holds the sets alone
    if (!str_empty(src)) md_script_ir_compile_from_source(p->v1, src, desc.sys, nullptr);

    const size_t ne = md_script_ir_num_errors(p->v1);
    const md_log_token_t* e = md_script_ir_errors(p->v1);
    for (size_t i = 0; i < ne; ++i) md_array_push(p->errs, (Diagnostic{range_of(e[i].range), e[i].text, e[i].context}), alloc);
    const size_t nw = md_script_ir_num_warnings(p->v1);
    const md_log_token_t* w = md_script_ir_warnings(p->v1);
    for (size_t i = 0; i < nw; ++i) md_array_push(p->warns, (Diagnostic{range_of(w[i].range), w[i].text, w[i].context}), alloc);
    const size_t nt = md_script_ir_num_vis_tokens(p->v1);
    const md_script_vis_token_t* t = md_script_ir_vis_tokens(p->v1);
    for (size_t i = 0; i < nt; ++i) md_array_push(p->marks, (Mark{range_of(t[i].range), t[i].depth, t[i].text, ref_of(t[i].ref)}), alloc);
    return p;
}

void program_free(Program* p) {
    if (!p) return;
    md_allocator_i* alloc = p->alloc;
    md_array_free(p->errs, alloc);
    md_array_free(p->warns, alloc);
    md_array_free(p->marks, alloc);
    if (p->v1) md_script_ir_free(p->v1);
    md_free(alloc, p, sizeof(Program));
}

Language language(const Program* p) { return p ? p->lang : Language::V1; }
bool     valid(const Program* p) { return p && md_script_ir_valid(p->v1); }
uint64_t fingerprint(const Program* p) { return p ? md_script_ir_fingerprint(p->v1) : 0; }

size_t            num_errors(const Program* p) { return p ? md_array_size(p->errs) : 0; }
const Diagnostic* errors(const Program* p) { return p ? p->errs : nullptr; }
size_t            num_warnings(const Program* p) { return p ? md_array_size(p->warns) : 0; }
const Diagnostic* warnings(const Program* p) { return p ? p->warns : nullptr; }
size_t            num_marks(const Program* p) { return p ? md_array_size(p->marks) : 0; }
const Mark*       marks(const Program* p) { return p ? p->marks : nullptr; }

void write_diagnostics(md_writer_i* w, const Program* p) { md_script_ir_write_diagnostics(w, v1_ir(p)); }

size_t       num_identifiers(const Program* p) { return p ? md_script_ir_num_identifiers(p->v1) : 0; }
const str_t* identifiers(const Program* p) { return p ? md_script_ir_identifiers(p->v1) : nullptr; }

// ### OUTPUTS ###

size_t       num_outputs(const Program* p) { return md_script_ir_property_count(v1_ir(p)); }
const str_t* output_names(const Program* p) { return md_script_ir_property_names(v1_ir(p)); }
uint32_t     output_flags(const Program* p, str_t name) { return (uint32_t)md_script_ir_property_flags(v1_ir(p), name); }
Range        output_range(const Program* p, str_t name) { return range_of(md_script_ir_property_range(v1_ir(p), name)); }
uint64_t     output_fingerprint(const Program* p, str_t name) { return md_script_ir_property_fingerprint(v1_ir(p), name); }
Ref          output_ref(const Program* p, str_t name) { return ref_of(md_script_ir_property_vis_ref(v1_ir(p), name)); }

bool  ref_valid(const Program* p, Ref ref) { return md_script_vis_ref_valid(v1_ir(p), v1_ref(ref)); }
str_t ref_name(const Program* p, Ref ref) { return md_script_vis_ref_ident(v1_ir(p), v1_ref(ref)); }
int   ref_dim(const Program* p, Ref ref) { return md_script_vis_ref_dim(v1_ir(p), v1_ref(ref)); }

// ### EVALUATE ###

Eval* eval_create(const Eval* prev, const Program* prog, const EvalDesc& desc, md_allocator_i* alloc) {
    md_script_eval_desc_t d = {};
    d.num_frames = desc.num_frames;
    md_script_eval_t* v = md_script_eval_create_from(prev ? prev->v1 : nullptr, v1_ir(prog), &d, alloc);
    if (!v) return nullptr;
    Eval* e = (Eval*)md_alloc(alloc, sizeof(Eval));
    memset(e, 0, sizeof(*e));
    e->lang = Language::V1;
    e->alloc = alloc;
    e->v1 = v;
    return e;
}

void eval_free(Eval* e) {
    if (!e) return;
    md_script_eval_free(e->v1);
    md_free(e->alloc, e, sizeof(Eval));
}

uint64_t eval_fingerprint(const Eval* e) { return e ? md_script_eval_ir_fingerprint(e->v1) : 0; }
size_t   eval_num_outputs(const Eval* e) { return e ? md_script_eval_property_count(e->v1) : 0; }
uint64_t eval_output_fingerprint(const Eval* e, str_t name) { return e ? md_script_eval_property_fingerprint(e->v1, name) : 0; }
size_t   eval_num_frames(const Eval* e) { return e ? md_script_eval_frame_count(e->v1) : 0; }

md_bitset_t eval_output_frames(const Eval* e, str_t name, md_allocator_i* alloc) {
    return e ? md_script_eval_property_frames(e->v1, name, alloc) : md_bitset_empty_set();
}

md_bitset_t eval_pending_frames(Eval* e, md_allocator_i* alloc) {
    return e ? md_script_eval_pending_frames(e->v1, alloc) : md_bitset_empty_set();
}

bool eval_frames(Eval* e, const Program* prog, const md_system_t* sys, str_t run, const uint32_t* frames, size_t num_frames) {
    return e && md_script_eval_frames(e->v1, v1_ir(prog), sys, run, frames, num_frames);
}

bool eval_whole(const Eval*) { return false; }

// v1 logs what fails while evaluating; it keeps no diagnostics of it
EvalErrors eval_errors(const Eval*) { return EvalErrors{}; }

float eval_progress(const Eval* e) {
    const size_t n = e ? md_script_eval_frame_count(e->v1) : 0;
    if (!n) return 1.0f;
    md_allocator_i* temp = md_get_heap_allocator();
    md_bitset_t pending = md_script_eval_pending_frames(e->v1, temp);
    const size_t left = md_bitset_count(pending);
    md_bitset_free(&pending, temp);
    return 1.0f - (float)left / (float)n;
}

void eval_interrupt(Eval* e) { if (e) md_script_eval_interrupt(e->v1); }
void eval_reset_interrupt(Eval* e) { if (e) md_script_eval_reset_interrupt(e->v1); }

const md_attributes_t* attributes(const Eval* e) { return e ? md_script_eval_attributes(e->v1) : nullptr; }

// ### WHAT AN EDIT CHANGES ###

bool plan(Plan* out, size_t cap, md_bitset_t* out_pending, const Eval* eval, const Program* prog, size_t num_frames, md_allocator_i* alloc) {
    const size_t n = md_script_ir_property_count(v1_ir(prog));
    md_script_eval_plan_t* tmp = (md_script_eval_plan_t*)md_alloc(alloc, (n ? n : 1) * sizeof(md_script_eval_plan_t));
    const bool ok = md_script_eval_plan(tmp, cap, out_pending, eval ? eval->v1 : nullptr, v1_ir(prog), num_frames, alloc);
    if (ok) {
        for (size_t i = 0; i < n; ++i) {
            out[i].status      = tmp[i].status == MD_SCRIPT_EVAL_STATUS_HELD    ? Status::Held :
                                 tmp[i].status == MD_SCRIPT_EVAL_STATUS_CHANGED ? Status::Changed : Status::New;
            out[i].from        = tmp[i].from;
            out[i].num_frames  = tmp[i].num_frames;
            out[i].num_pending = tmp[i].num_pending;
        }
    }
    md_free(alloc, tmp, (n ? n : 1) * sizeof(md_script_eval_plan_t));
    return ok;
}

// ### AGGREGATES ###

static inline AggregateInfo info_of(const md_script_eval_slot_info_t& s) {
    AggregateInfo a;
    a.version = s.version;
    a.num_frames = s.num_frames;
    a.num_pending = s.num_pending;
    a.num_values = s.num_values;
    return a;
}

bool aggregate_request(Eval* e, uint64_t key, str_t name, const md_bitset_t* frames) {
    return e && md_script_eval_slot_request(e->v1, key, name, frames);
}

void aggregate_release(Eval* e, uint64_t key, str_t name) {
    if (e) md_script_eval_slot_release(e->v1, key, name);
}

bool aggregate_info(AggregateInfo* out, const Eval* e, uint64_t key, str_t name) {
    md_script_eval_slot_info_t s = {};
    if (!e || !md_script_eval_slot_info(&s, e->v1, key, name)) return false;
    *out = info_of(s);
    return true;
}

bool aggregate_read(AggregateInfo* out, float* values, float* weights, size_t num_values, const Eval* e, uint64_t key, str_t name) {
    md_script_eval_slot_info_t s = {};
    if (!e) return false;
    const bool ok = md_script_eval_slot_extract(&s, values, weights, num_values, e->v1, key, name);
    if (out) *out = info_of(s);
    return ok;
}

// ### DRAWING ###

static_assert(sizeof(Vertex) == sizeof(md_script_vis_vertex_t), "script::Vertex mirrors md_script_vis_vertex_t");

void drawing_init(Drawing* d, md_allocator_i* alloc) {
    *d = Drawing{};
    d->alloc = alloc;
}

void drawing_clear(Drawing* d) {
    md_allocator_i* alloc = d->alloc;
    md_array_free(d->points, alloc);
    md_array_free(d->lines, alloc);
    md_array_free(d->triangles, alloc);
    md_array_free(d->text, alloc);
    for (size_t i = 0; i < md_array_size(d->sdf.structures); ++i) md_bitset_free(&d->sdf.structures[i], alloc);
    md_array_free(d->sdf.structures, alloc);
    md_array_free(d->sdf.matrices, alloc);
    for (size_t i = 0; i < md_array_size(d->structures); ++i) md_bitset_free(&d->structures[i], alloc);
    md_array_free(d->structures, alloc);
    for (size_t i = 0; i < md_array_size(d->groups); ++i) md_bitset_free(&d->groups[i], alloc);
    md_array_free(d->groups, alloc);
    md_bitset_free(&d->atom_mask, alloc);
    v1::drawing_init(d, alloc);
}

static md_bitset_t copy_set(const md_bitset_t& s, md_allocator_i* alloc) {
    return md_bitset_copy(s, alloc);
}

// v1's drawing into d's arrays, appended
static void take(Drawing* d, const md_script_vis_t& v) {
    md_allocator_i* alloc = d->alloc;
    for (size_t i = 0; i < md_array_size(v.points); ++i) {
        Vertex x; memcpy(&x, &v.points[i], sizeof(x)); md_array_push(d->points, x, alloc);
    }
    for (size_t i = 0; i < md_array_size(v.lines); ++i) {
        Vertex x; memcpy(&x, &v.lines[i], sizeof(x)); md_array_push(d->lines, x, alloc);
    }
    for (size_t i = 0; i < md_array_size(v.triangles); ++i) {
        Vertex x; memcpy(&x, &v.triangles[i], sizeof(x)); md_array_push(d->triangles, x, alloc);
    }
    for (size_t i = 0; i < md_array_size(v.text); ++i) {
        Text t;
        t.pos = v.text[i].pos;
        t.str = str_copy(v.text[i].str, alloc);
        t.value = v.text[i].value;
        t.unit = v.text[i].unit;
        md_array_push(d->text, t, alloc);
    }
    for (size_t i = 0; i < md_array_size(v.sdf.matrices); ++i) md_array_push(d->sdf.matrices, v.sdf.matrices[i], alloc);
    for (size_t i = 0; i < md_array_size(v.sdf.structures); ++i) md_array_push(d->sdf.structures, copy_set(v.sdf.structures[i], alloc), alloc);
    d->sdf.extent = v.sdf.extent;
    for (size_t i = 0; i < md_array_size(v.structures); ++i) md_array_push(d->structures, copy_set(v.structures[i], alloc), alloc);
    for (size_t i = 0; i < md_array_size(v.groups); ++i) md_array_push(d->groups, copy_set(v.groups[i], alloc), alloc);
    md_bitset_t mask = md_bitset_or(d->atom_mask, v.atom_mask, alloc);
    md_bitset_free(&d->atom_mask, alloc);
    d->atom_mask = mask;
}

static md_script_vis_ctx_t v1_ctx(const DrawCtx& ctx) {
    md_script_vis_ctx_t c = {};
    c.ir = v1_ir(ctx.prog);
    c.sys = ctx.sys;
    c.state = ctx.state;
    return c;
}

bool draw_ref(Drawing* out, Ref ref, int subidx, const DrawCtx& ctx, uint32_t flags) {
    md_script_vis_t v = {};
    md_script_vis_init(&v, out->alloc);
    const md_script_vis_ctx_t c = v1_ctx(ctx);
    const bool ok = md_script_vis_eval_ref(&v, v1_ref(ref), subidx, &c, flags);
    if (ok) take(out, v);
    md_script_vis_free(&v);
    return ok;
}

bool draw_text(Drawing* out, str_t expr, const DrawCtx& ctx, uint32_t flags) {
    md_script_vis_t v = {};
    md_script_vis_init(&v, out->alloc);
    const md_script_vis_ctx_t c = v1_ctx(ctx);
    const bool ok = md_script_vis_eval_string(&v, expr, &c, flags);
    if (ok) take(out, v);
    md_script_vis_free(&v);
    return ok;
}

size_t value_text(char* buf, size_t cap, Ref ref, int subidx, const DrawCtx& ctx) {
    const md_script_vis_ctx_t c = v1_ctx(ctx);
    return md_script_vis_value_text(buf, cap, v1_ref(ref), subidx, &c);
}

// ### EDITOR ###

static_assert(sizeof(Completion) == sizeof(md_script_completion_t), "script::Completion mirrors md_script_completion_t");

static Completions completions_of(const md_script_completions_t& c, md_allocator_i* alloc) {
    Completions r;
    r.range = range_of(c.range);
    r.prefix = c.prefix;
    r.count = c.count;
    r.items = c.count ? (Completion*)md_alloc(alloc, c.count * sizeof(Completion)) : nullptr;
    for (size_t i = 0; i < c.count; ++i) {
        r.items[i].label = c.items[i].label;
        r.items[i].text  = c.items[i].text;
        r.items[i].kind  = (CompletionKind)c.items[i].kind;
        r.items[i].fits  = c.items[i].fits;
    }
    return r;
}

size_t       num_keywords(Language) { return md_script_num_keywords(); }
const str_t* keywords(Language) { return md_script_keywords(); }
size_t       builtin_names(Language, str_t* out, size_t cap) { return md_script_builtin_identifiers(out, cap); }

Completions complete(Language, str_t src, int cursor, const md_system_t* sys, md_allocator_i* alloc) {
    return completions_of(md_script_complete(src, cursor, sys, alloc), alloc);
}

str_t attr_path(Language, const md_attributes_t* table, const md_attribute_t* attr) { return md_script_attr_path(table, attr); }

bool name_valid(Language, str_t name) { return md_script_identifier_name_valid(name); }

// ### FILTERS ###

static md_filter_context_t v1_filter_ctx(const FilterCtx& ctx) {
    md_filter_context_t c = {};
    c.sys = ctx.sys;
    c.state = ctx.state;
    c.ir = v1_ir(ctx.prog);
    c.atoms = ctx.atoms;
    return c;
}

static FilterResult result_of(const md_filter_result_t& r) {
    FilterResult f;
    f.success = r.success;
    f.is_dynamic = r.is_dynamic;
    f.error_range = range_of(r.error_range);
    memcpy(f.error, r.error, sizeof(f.error) < sizeof(r.error) ? sizeof(f.error) : sizeof(r.error));
    f.error[sizeof(f.error) - 1] = '\0';
    return f;
}

FilterResult filter(Language, md_bitset_t* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc) {
    const md_filter_context_t c = v1_filter_ctx(ctx);
    return result_of(md_filter(dst, expr, &c, alloc));
}

FilterResult filter_sets(Language, md_array(md_bitset_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc) {
    const md_filter_context_t c = v1_filter_ctx(ctx);
    return result_of(md_filter_sets(dst, expr, &c, alloc));
}

FilterResult filter_bonds(Language, md_array(md_atom_pair_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc) {
    const md_filter_context_t c = v1_filter_ctx(ctx);
    return result_of(md_filter_bonds(dst, expr, &c, alloc));
}

Completions filter_complete(Language, str_t expr, int cursor, const FilterCtx& ctx, md_allocator_i* alloc) {
    const md_filter_context_t c = v1_filter_ctx(ctx);
    return completions_of(md_filter_complete(expr, cursor, &c, alloc), alloc);
}

}  // namespace v1
}  // namespace script
