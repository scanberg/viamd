// The script interface (script_engine.h) for v2: VIAL, mdlib's md_vial.

#include "script_engine_impl.h"

#include <md_vial.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_common.h>
#include <core/md_writer.h>

#include <md_attributes.h>

#include <string.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

namespace script {

// A program: VIAL's, and what the interface hands out of it, made when compiling
struct ProgramV2 {
    md_vial_program_t*          prog;       // NULL for empty source
    uint64_t                    id;         // Ref.program: this compilation, never 0
    const md_system_t*          sys;
    str_t                       run;
    md_array(md_vial_named_set_t) sets;     // the stored selections, as VIAL takes them
    md_array(str_t)             idents;     // every top-level name it binds, and the stored selections'
    md_array(str_t)             outputs;    // the published values plots read: not selections
    md_array(uint32_t)          output_idx; // into md_vial_outputs, by outputs
};

namespace v2 {


static inline ProgramV2* pv(const Program* p) { return p && p->lang == Language::V2 ? p->v2 : nullptr; }
static inline const md_vial_program_t* vp(const Program* p) { const ProgramV2* x = pv(p); return x ? x->prog : nullptr; }

static uint64_t next_program_id() {
    static volatile uint64_t counter = 0;
    return __atomic_add_fetch(&counter, 1, __ATOMIC_RELAXED) | (1ull << 62);  // apart from v1's ids
}

static inline Range range_of(md_vial_span_t s) { return Range{(int)s.beg, (int)s.end}; }

// The output of a name, among md_vial_outputs; nullptr if none
static const md_vial_output_t* output_of(const Program* p, str_t name) {
    const ProgramV2* x = pv(p);
    if (!x || !x->prog) return nullptr;
    const md_vial_output_t* o;
    const size_t n = md_vial_outputs(x->prog, &o);
    for (size_t i = 0; i < n; ++i) if (str_eq(o[i].name, name)) return &o[i];
    return nullptr;
}

static FilterResult not_yet() {
    FilterResult r;
    snprintf(r.error, sizeof(r.error), "not yet in VIAL");
    return r;
}

Program* compile(Language lang, str_t src, const CompileDesc& desc, md_allocator_i* alloc) {
    Program* p = (Program*)md_alloc(alloc, sizeof(Program));
    memset(p, 0, sizeof(*p));
    p->lang = lang;
    p->alloc = alloc;
    ProgramV2* x = (ProgramV2*)md_alloc(alloc, sizeof(ProgramV2));
    memset(x, 0, sizeof(*x));
    p->v2 = x;
    x->id = next_program_id();
    x->sys = desc.sys;
    x->run = str_copy(desc.run, alloc);
    // The stored selections, as lists of atoms
    for (size_t i = 0; i < desc.num_sets; ++i) {
        const md_bitset_t& b = desc.sets[i].atoms;
        const size_t n = md_bitset_count(b);
        uint32_t* atoms = (uint32_t*)md_alloc(alloc, (n ? n : 1) * sizeof(uint32_t));
        size_t k = 0;
        md_bitset_iter_t it = md_bitset_iter(b);
        while (k < n && md_bitset_iter_next(&it)) atoms[k++] = (uint32_t)it.idx;
        md_vial_named_set_t ns = {str_copy(desc.sets[i].name, alloc), atoms, k};
        md_array_push(x->sets, ns, alloc);
        md_array_push(x->idents, ns.name, alloc);
    }
    // Empty source is not compiled: the program holds the sets alone
    if (str_empty(src)) return p;
    md_vial_compile_desc_t cd = {};
    cd.sys = desc.sys;
    cd.run = x->run;
    cd.sets = x->sets;
    cd.num_sets = md_array_size(x->sets);
    x->prog = md_vial_compile(src, &cd, alloc);
    if (!x->prog) return p;

    const md_vial_diagnostic_t* dg;
    const size_t nd = md_vial_diagnostics(x->prog, &dg);
    for (size_t i = 0; i < nd; ++i) {
        // The code first, as the reference lists it
        char buf[512];
        const int len = snprintf(buf, sizeof(buf), "%.5s: %.*s", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
        Diagnostic d = {range_of(dg[i].span), str_copy(str_t{buf, (size_t)(len > 0 ? MIN(len, (int)sizeof(buf) - 1) : 0)}, alloc), nullptr};
        if (dg[i].severity == MD_VIAL_SEVERITY_ERROR) md_array_push(p->errs, d, alloc);
        else md_array_push(p->warns, d, alloc);
    }
    // Marks: every expression with a type, its type and, where known when compiling, its value
    const md_vial_mark_t* mk;
    const size_t nm = md_vial_marks(x->prog, &mk);
    for (size_t i = 0; i < nm; ++i) {
        md_vial_hover_t h;
        if (!md_vial_expr_info(x->prog, mk[i].expr, &h)) continue;
        char buf[256];
        const int len = h.value[0] ? snprintf(buf, sizeof(buf), "%s -> %s", h.type, h.value) : snprintf(buf, sizeof(buf), "%s", h.type);
        Mark m = {range_of(mk[i].span), (int)mk[i].depth, str_copy(str_t{buf, (size_t)(len > 0 ? MIN(len, (int)sizeof(buf) - 1) : 0)}, alloc), Ref{x->id, mk[i].expr}};
        md_array_push(p->marks, m, alloc);
    }
    const md_vial_output_t* o;
    const size_t no = md_vial_outputs(x->prog, &o);
    for (size_t i = 0; i < no; ++i) {
        md_array_push(x->idents, o[i].name, alloc);
        if (o[i].selection) continue;
        md_array_push(x->outputs, o[i].name, alloc);
        md_array_push(x->output_idx, (uint32_t)i, alloc);
    }
    return p;
}

void program_free(Program* p) {
    if (!p) return;
    md_allocator_i* alloc = p->alloc;
    if (ProgramV2* x = p->v2) {
        if (x->prog) md_vial_program_free(x->prog);
        md_array_free(x->sets, alloc);
        md_array_free(x->idents, alloc);
        md_array_free(x->outputs, alloc);
        md_array_free(x->output_idx, alloc);
        md_free(alloc, x, sizeof(ProgramV2));
    }
    md_array_free(p->errs, alloc);
    md_array_free(p->warns, alloc);
    md_array_free(p->marks, alloc);
    md_free(alloc, p, sizeof(Program));
}

Language language(const Program*) { return Language::V2; }

bool valid(const Program* p) { return vp(p) && md_array_size(p->errs) == 0; }

uint64_t fingerprint(const Program* p) {
    // What every output is computed from, and what it is called
    const md_vial_program_t* v = vp(p);
    if (!v) return 0;
    const md_vial_output_t* o;
    const size_t n = md_vial_outputs(v, &o);
    uint64_t h = 0xCBF29CE484222325ull;
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < o[i].name.len; ++k) h = (h ^ (uint8_t)o[i].name.ptr[k]) * 0x100000001B3ull;
        h = (h ^ o[i].fingerprint) * 0x100000001B3ull;
    }
    return h ? h : 1;
}

size_t num_errors(const Program* p) { return p ? md_array_size(p->errs) : 0; }

const Diagnostic* errors(const Program* p) { return p ? p->errs : nullptr; }

size_t num_warnings(const Program* p) { return p ? md_array_size(p->warns) : 0; }

const Diagnostic* warnings(const Program* p) { return p ? p->warns : nullptr; }

size_t num_marks(const Program* p) { return p ? md_array_size(p->marks) : 0; }

const Mark* marks(const Program* p) { return p ? p->marks : nullptr; }

void write_diagnostics(md_writer_i* w, const Program* p) {
    if (!w || !p) return;
    char buf[640];
    for (int pass = 0; pass < 2; ++pass) {
        const Diagnostic* d = pass == 0 ? p->errs : p->warns;
        for (size_t i = 0; i < md_array_size(d); ++i) {
            const int len = snprintf(buf, sizeof(buf), "%s [%d, %d): %.*s\n", pass == 0 ? "error" : "warning", d[i].range.beg, d[i].range.end, (int)d[i].text.len, d[i].text.ptr);
            if (len > 0) md_write(w, buf, (size_t)MIN(len, (int)sizeof(buf) - 1));
        }
    }
}

size_t num_identifiers(const Program* p) { return pv(p) ? md_array_size(pv(p)->idents) : 0; }

const str_t* identifiers(const Program* p) { return pv(p) ? pv(p)->idents : nullptr; }

size_t num_outputs(const Program* p) { return pv(p) ? md_array_size(pv(p)->outputs) : 0; }

const str_t* output_names(const Program* p) { return pv(p) ? pv(p)->outputs : nullptr; }

uint32_t output_flags(const Program* p, str_t name) {
    const md_vial_output_t* o = output_of(p, name);
    if (!o || o->selection) return 0;
    uint32_t f = 0;
    if (o->flags & MD_VIAL_OUTPUT_TEMPORAL)     f |= OUTPUT_TEMPORAL;
    if (o->flags & MD_VIAL_OUTPUT_DISTRIBUTION) f |= OUTPUT_DISTRIBUTION;
    if (o->flags & MD_VIAL_OUTPUT_PERIODIC)     f |= OUTPUT_PERIODIC;
    return f;
}

Range output_range(const Program* p, str_t name) { const md_vial_output_t* o = output_of(p, name); return o ? range_of(o->span) : Range{}; }

uint64_t output_fingerprint(const Program* p, str_t name) { const md_vial_output_t* o = output_of(p, name); return o ? o->fingerprint : 0; }

Ref output_ref(const Program* p, str_t name) {
    const md_vial_output_t* o = output_of(p, name);
    if (!o) return Ref{};
    const uint32_t e = md_vial_expr_at(vp(p), o->span.beg);
    return e ? Ref{pv(p)->id, e} : Ref{};
}

bool ref_valid(const Program* p, Ref ref) { return pv(p) && ref.program == pv(p)->id && ref.node != 0; }

str_t ref_name(const Program* p, Ref ref) {
    // The output it is the name of, as plots know it; empty for any other expression
    if (!v2::ref_valid(p, ref)) return str_t{};
    const md_vial_output_t* o;
    const size_t n = md_vial_outputs(vp(p), &o);
    for (size_t i = 0; i < n; ++i) {
        if (!o[i].selection && md_vial_expr_at(vp(p), o[i].span.beg) == ref.node) return o[i].name;
    }
    return str_t{};
}

int ref_dim(const Program* p, Ref ref) {
    // The outer extent of its type, as hover prints it: float[433][3] is 433, a scalar 1, an extent known only when
    // evaluating -1
    if (!v2::ref_valid(p, ref)) return 0;
    md_vial_hover_t h;
    if (!md_vial_expr_info(vp(p), ref.node, &h)) return 0;
    const char* b = strchr(h.type, '[');
    if (!b) return 1;
    char* end = nullptr;
    const long n = strtol(b + 1, &end, 10);
    return end && *end == ']' && n >= 0 ? (int)n : -1;
}

Eval* eval_create(const Eval* prev,  const Program* prog,  const EvalDesc& desc,  md_allocator_i* alloc) {
    (void)prev; (void)prog; (void)desc; (void)alloc;
    return nullptr;
}

void eval_free(Eval* eval) {
    (void)eval;
    
}

uint64_t eval_fingerprint(const Eval* eval) {
    (void)eval;
    return 0;
}

size_t eval_num_outputs(const Eval* eval) {
    (void)eval;
    return 0;
}

uint64_t eval_output_fingerprint(const Eval* eval,  str_t name) {
    (void)eval; (void)name;
    return 0;
}

size_t eval_num_frames(const Eval* eval) {
    (void)eval;
    return 0;
}

md_bitset_t eval_output_frames(const Eval* eval,  str_t name,  md_allocator_i* alloc) {
    (void)eval; (void)name; (void)alloc;
    return md_bitset_empty_set();
}

md_bitset_t eval_pending_frames(Eval* eval,  md_allocator_i* alloc) {
    (void)eval; (void)alloc;
    return md_bitset_empty_set();
}

bool eval_frames(Eval* eval,  const Program* prog,  const md_system_t* sys,  str_t run,  const uint32_t* frames,  size_t num_frames) {
    (void)eval; (void)prog; (void)sys; (void)run; (void)frames; (void)num_frames;
    return false;
}

bool eval_whole(const Eval* eval) {
    (void)eval;
    return false;
}

float eval_progress(const Eval* eval) {
    (void)eval;
    return 1.0f;
}

void eval_interrupt(Eval* eval) {
    (void)eval;
    
}

void eval_reset_interrupt(Eval* eval) {
    (void)eval;
    
}

const md_attributes_t* attributes(const Eval* eval) {
    (void)eval;
    return nullptr;
}

bool plan(Plan* out,  size_t cap,  md_bitset_t* out_pending,  const Eval* eval,  const Program* prog,  size_t num_frames,  md_allocator_i* alloc) {
    (void)out; (void)cap; (void)out_pending; (void)eval; (void)prog; (void)num_frames; (void)alloc;
    return false;
}

bool aggregate_request(Eval* eval,  uint64_t key,  str_t name,  const md_bitset_t* frames) {
    (void)eval; (void)key; (void)name; (void)frames;
    return false;
}

void aggregate_release(Eval* eval,  uint64_t key,  str_t name) {
    (void)eval; (void)key; (void)name;
    
}

bool aggregate_info(AggregateInfo* out,  const Eval* eval,  uint64_t key,  str_t name) {
    (void)out; (void)eval; (void)key; (void)name;
    return false;
}

bool aggregate_read(AggregateInfo* out,  float* values,  float* weights,  size_t num_values,  const Eval* eval,  uint64_t key,  str_t name) {
    (void)out; (void)values; (void)weights; (void)num_values; (void)eval; (void)key; (void)name;
    return false;
}

bool draw_ref(Drawing* out,  Ref ref,  int subidx,  const DrawCtx& ctx,  uint32_t flags) {
    (void)out; (void)ref; (void)subidx; (void)ctx; (void)flags;
    return false;
}

bool draw_text(Drawing* out,  str_t expr,  const DrawCtx& ctx,  uint32_t flags) {
    (void)out; (void)expr; (void)ctx; (void)flags;
    return false;
}

size_t value_text(char* buf,  size_t cap,  Ref ref,  int subidx,  const DrawCtx& ctx) {
    (void)buf; (void)cap; (void)ref; (void)subidx; (void)ctx;
    return 0;
}

// The language's lists, made once
static const std::vector<str_t>& keyword_list() {
    static const std::vector<str_t> list = [] {
        std::vector<str_t> v(md_vial_keywords(nullptr, 0));
        md_vial_keywords(v.data(), v.size());
        return v;
    }();
    return list;
}

size_t num_keywords(Language) { return keyword_list().size(); }

const str_t* keywords(Language) { return keyword_list().data(); }

size_t builtin_names(Language, str_t* out, size_t cap) { return md_vial_builtin_names(out, cap); }

static CompletionKind kind_of(md_vial_completion_kind_t k) {
    switch (k) {
    case MD_VIAL_COMPLETION_KEYWORD:   return CompletionKind::Keyword;
    case MD_VIAL_COMPLETION_FUNCTION:  return CompletionKind::Procedure;
    case MD_VIAL_COMPLETION_PROPERTY:  return CompletionKind::Procedure;
    case MD_VIAL_COMPLETION_CONSTANT:  return CompletionKind::Constant;
    case MD_VIAL_COMPLETION_VARIABLE:  return CompletionKind::Variable;
    case MD_VIAL_COMPLETION_PARAMETER: return CompletionKind::Parameter;
    case MD_VIAL_COMPLETION_VALUE:     return CompletionKind::Value;
    }
    return CompletionKind::Keyword;
}

static Completions completions_of(const md_vial_completions_t& c, md_allocator_i* alloc) {
    Completions r;
    r.range = range_of(c.span);
    r.prefix = c.prefix;
    r.count = c.count;
    r.items = c.count ? (Completion*)md_alloc(alloc, c.count * sizeof(Completion)) : nullptr;
    for (size_t i = 0; i < c.count; ++i) {
        r.items[i].label = c.items[i].label;
        r.items[i].text  = c.items[i].text;
        r.items[i].kind  = kind_of(c.items[i].kind);
        // A parameter of the call, or a value of the keyword: what fits where the cursor is
        r.items[i].fits  = c.items[i].kind == MD_VIAL_COMPLETION_PARAMETER || c.items[i].kind == MD_VIAL_COMPLETION_VALUE;
    }
    return r;
}

Completions complete(Language, str_t src, int cursor, const md_system_t* sys, md_allocator_i* alloc) {
    if (cursor < 0 || (size_t)cursor > src.len) return Completions{};
    md_vial_compile_desc_t d = {};
    d.sys = sys;
    md_vial_program_t* prog = md_vial_compile(src, &d, alloc);
    Completions r;
    md_vial_completions_t c;
    if (prog && md_vial_complete(&c, prog, (uint32_t)cursor, alloc)) r = completions_of(c, alloc);
    // What the completions hold is alloc's, not the program's
    if (prog) md_vial_program_free(prog);
    return r;
}

str_t attr_path(Language, const md_attributes_t* table, const md_attribute_t* attr) {
    // attr(path) reads <run>/<path>: the run is the prefix of attr's path with a frame axis, <run>/time
    if (!table || !attr) return str_t{};
    const str_t path = attr->path;
    char buf[512];
    for (size_t i = 1; i < path.len && i + 6 < sizeof(buf); ++i) {
        if (path.ptr[i] != '/') continue;
        memcpy(buf, path.ptr, i);
        memcpy(buf + i, "/time", 5);
        const md_attribute_t* t = md_attributes_find(table, str_t{buf, i + 5});
        if (t && t != attr) return md_vial_attr_path(str_t{path.ptr, i}, path);
    }
    return str_t{};
}

bool name_valid(Language, str_t name) { return md_vial_name_valid(name); }

FilterResult filter(Language lang,  md_bitset_t* dst,  str_t expr,  const FilterCtx& ctx,  md_allocator_i* alloc) {
    (void)lang; (void)dst; (void)expr; (void)ctx; (void)alloc;
    return not_yet();
}

FilterResult filter_sets(Language lang,  md_array(md_bitset_t)* dst,  str_t expr,  const FilterCtx& ctx,  md_allocator_i* alloc) {
    (void)lang; (void)dst; (void)expr; (void)ctx; (void)alloc;
    return not_yet();
}

FilterResult filter_bonds(Language lang,  md_array(md_atom_pair_t)* dst,  str_t expr,  const FilterCtx& ctx,  md_allocator_i* alloc) {
    (void)lang; (void)dst; (void)expr; (void)ctx; (void)alloc;
    return not_yet();
}

Completions filter_complete(Language lang,  str_t expr,  int cursor,  const FilterCtx& ctx,  md_allocator_i* alloc) {
    (void)lang; (void)expr; (void)cursor; (void)ctx; (void)alloc;
    return Completions{};
}

}  // namespace v2
}  // namespace script
