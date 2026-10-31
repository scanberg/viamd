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

// An aggregate asked for over chosen frames (a plot period): recomputed from the kept per-frame values whenever its
// frames change or a run finishes
struct AggSlot {
    uint64_t        key;
    char            name[64];
    md_bitset_t     frames;
    uint64_t        version;
    uint64_t        seen_runs;  // the runs it was computed after
    bool            computed;
    md_array(float) values;
};

struct EvalV2 {
    md_vial_eval_t*     eval;
    const Program*      prog;
    uint64_t            fingerprint;
    size_t              num_frames;
    volatile uint64_t   running;    // eval_frames is running
    volatile uint64_t   runs;       // runs finished
    md_array(AggSlot)   slots;
};

namespace v2 {



static uint64_t next_version() {
    static volatile uint64_t counter = 0;
    return __atomic_add_fetch(&counter, 1, __ATOMIC_RELAXED);
}

static inline ProgramV2* pv(const Program* p) { return p && p->lang == Language::V2 ? p->v2 : nullptr; }
static inline EvalV2* ev(const Eval* e) { return e && e->lang == Language::V2 ? e->v2 : nullptr; }
static inline md_vial_eval_t* ve(const Eval* e) { EvalV2* x = ev(e); return x ? x->eval : nullptr; }
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

Eval* eval_create(const Eval* prev, const Program* prog, const EvalDesc& desc, md_allocator_i* alloc) {
    if (!v2::valid(prog) || !desc.sys) return nullptr;
    md_vial_eval_desc_t d = {};
    d.sys = desc.sys;
    d.run = desc.run;
    d.num_threads = desc.num_threads;
    d.parallel_for = desc.parallel_for;
    d.parallel_user = desc.parallel_user;
    // What an earlier evaluation of v2 computed alike is taken over
    md_vial_eval_t* v = ve(prev) ? md_vial_eval_create_from(ve(prev), vp(prog), &d, alloc) : md_vial_eval_create(vp(prog), &d, alloc);
    if (!v) return nullptr;
    Eval* e = (Eval*)md_alloc(alloc, sizeof(Eval));
    memset(e, 0, sizeof(*e));
    e->lang = Language::V2;
    e->alloc = alloc;
    EvalV2* x = (EvalV2*)md_alloc(alloc, sizeof(EvalV2));
    memset(x, 0, sizeof(*x));
    x->eval = v;
    x->prog = prog;
    x->fingerprint = v2::fingerprint(prog);
    x->num_frames = desc.num_frames;
    e->v2 = x;
    // The periods' requests come along, to be computed for this evaluation
    if (EvalV2* px = ev(prev)) {
        for (size_t i = 0; i < md_array_size(px->slots); ++i) {
            AggSlot s = {};
            s.key = px->slots[i].key;
            memcpy(s.name, px->slots[i].name, sizeof(s.name));
            s.frames = md_bitset_copy(px->slots[i].frames, alloc);
            s.version = next_version();
            md_array_push(x->slots, s, alloc);
        }
    }
    return e;
}

void eval_free(Eval* e) {
    EvalV2* x = ev(e);
    if (!x) return;
    md_allocator_i* alloc = e->alloc;
    md_vial_eval_free(x->eval);
    for (size_t i = 0; i < md_array_size(x->slots); ++i) {
        md_bitset_free(&x->slots[i].frames, alloc);
        md_array_free(x->slots[i].values, alloc);
    }
    md_array_free(x->slots, alloc);
    md_free(alloc, x, sizeof(EvalV2));
    md_free(alloc, e, sizeof(Eval));
}

uint64_t eval_fingerprint(const Eval* e) { return ev(e) ? ev(e)->fingerprint : 0; }

size_t eval_num_outputs(const Eval* e) { return ev(e) ? v2::num_outputs(ev(e)->prog) : 0; }

uint64_t eval_output_fingerprint(const Eval* e, str_t name) { return ev(e) ? v2::output_fingerprint(ev(e)->prog, name) : 0; }

size_t eval_num_frames(const Eval* e) { return ev(e) ? ev(e)->num_frames : 0; }

md_bitset_t eval_output_frames(const Eval* e, str_t name, md_allocator_i* alloc) {
    EvalV2* x = ev(e);
    if (!x || !x->num_frames) return md_bitset_empty_set();
    md_allocator_i* temp = md_get_heap_allocator();
    uint8_t* held = (uint8_t*)md_alloc(temp, x->num_frames);
    const size_t n = md_vial_eval_output_frames(x->eval, name, held, x->num_frames);
    md_bitset_t r = md_bitset_empty_set();
    if (n == x->num_frames) {
        uint32_t* idx = (uint32_t*)md_alloc(temp, n * sizeof(uint32_t));
        size_t k = 0;
        for (size_t f = 0; f < n; ++f) if (held[f]) idx[k++] = (uint32_t)f;
        r = md_bitset_from_indices(idx, k, alloc);
        md_free(temp, idx, n * sizeof(uint32_t));
    }
    md_free(temp, held, x->num_frames);
    return r;
}

md_bitset_t eval_pending_frames(Eval* e, md_allocator_i* alloc) {
    // The frames some output does not hold yet; nothing when every output is held, as when all are taken over
    EvalV2* x = ev(e);
    if (!x || !x->num_frames) return md_bitset_empty_set();
    md_allocator_i* temp = md_get_heap_allocator();
    uint8_t* held = (uint8_t*)md_alloc(temp, x->num_frames);
    uint8_t* pending = (uint8_t*)md_alloc(temp, x->num_frames);
    memset(pending, 0, x->num_frames);
    const size_t no = v2::num_outputs(x->prog);
    const str_t* names = v2::output_names(x->prog);
    for (size_t i = 0; i < no; ++i) {
        if (md_vial_eval_output_frames(x->eval, names[i], held, x->num_frames) != x->num_frames) continue;
        for (size_t f = 0; f < x->num_frames; ++f) pending[f] |= !held[f];
    }
    uint32_t* idx = (uint32_t*)md_alloc(temp, x->num_frames * sizeof(uint32_t));
    size_t k = 0;
    for (size_t f = 0; f < x->num_frames; ++f) if (pending[f]) idx[k++] = (uint32_t)f;
    md_bitset_t r = md_bitset_from_indices(idx, k, alloc);
    md_free(temp, idx, x->num_frames * sizeof(uint32_t));
    md_free(temp, pending, x->num_frames);
    md_free(temp, held, x->num_frames);
    return r;
}

bool eval_frames(Eval* e, const Program*, const md_system_t*, str_t, const uint32_t*, size_t) {
    // Every frame in one run, spread over the evaluation's threads: a call while one runs does nothing
    EvalV2* x = ev(e);
    if (!x) return false;
    uint64_t idle = 0;
    if (!__atomic_compare_exchange_n(&x->running, &idle, 1, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) return true;
    const bool ok = md_vial_eval_run(x->eval);
    __atomic_add_fetch(&x->runs, 1, __ATOMIC_RELEASE);
    __atomic_store_n(&x->running, 0, __ATOMIC_RELEASE);
    return ok;
}

bool eval_whole(const Eval* e) { return ev(e) != nullptr; }

float eval_progress(const Eval* e) {
    EvalV2* x = ev(e);
    if (!x) return 1.0f;
    if (__atomic_load_n(&x->running, __ATOMIC_ACQUIRE)) return (float)md_vial_eval_progress(x->eval);
    return __atomic_load_n(&x->runs, __ATOMIC_ACQUIRE) ? 1.0f : 0.0f;
}

void eval_interrupt(Eval* e) { if (ve(e)) md_vial_eval_interrupt(ve(e)); }

void eval_reset_interrupt(Eval* e) { if (ve(e)) md_vial_eval_reset_interrupt(ve(e)); }

const md_attributes_t* attributes(const Eval* e) { return ve(e) ? md_vial_eval_attributes(ve(e)) : nullptr; }

bool plan(Plan* out, size_t cap, md_bitset_t* out_pending, const Eval* eval, const Program* prog, size_t num_frames, md_allocator_i* alloc) {
    ProgramV2* x = pv(prog);
    if (!v2::valid(prog) || cap < md_array_size(x->outputs)) return false;
    const md_vial_output_t* o;
    const size_t n = md_vial_outputs(x->prog, &o);
    md_allocator_i* temp = md_get_heap_allocator();
    md_vial_plan_t* tmp = (md_vial_plan_t*)md_alloc(temp, (n ? n : 1) * sizeof(md_vial_plan_t));
    const bool ok = md_vial_eval_plan(tmp, n, ve(eval), x->prog, num_frames);
    bool any_pending = false;
    if (ok) {
        // One per output plots see: md_vial_outputs less the selections
        for (size_t i = 0; i < md_array_size(x->outputs); ++i) {
            const md_vial_plan_t& t = tmp[x->output_idx[i]];
            out[i].status      = t.status == MD_VIAL_PLAN_HELD ? Status::Held : t.status == MD_VIAL_PLAN_CHANGED ? Status::Changed : Status::New;
            out[i].from        = t.from;
            out[i].num_frames  = t.num_frames;
            out[i].num_pending = t.num_pending;
            any_pending = any_pending || t.num_pending;
        }
    }
    md_free(temp, tmp, (n ? n : 1) * sizeof(md_vial_plan_t));
    // A run evaluates every frame for what it does not hold
    if (ok && out_pending) *out_pending = any_pending ? md_bitset_range(0, (uint32_t)num_frames) : md_bitset_empty_set();
    (void)alloc;
    return ok;
}

static AggSlot* find_slot(EvalV2* x, uint64_t key, str_t name) {
    for (size_t i = 0; i < md_array_size(x->slots); ++i) {
        AggSlot* s = &x->slots[i];
        if (s->key == key && str_eq(str_from_cstr(s->name), name)) return s;
    }
    return nullptr;
}

bool aggregate_request(Eval* e, uint64_t key, str_t name, const md_bitset_t* frames) {
    EvalV2* x = ev(e);
    if (!x || !frames || name.len >= sizeof(AggSlot::name) || !(v2::output_flags(x->prog, name) & OUTPUT_DISTRIBUTION)) return false;
    AggSlot* s = find_slot(x, key, name);
    if (!s) {
        AggSlot n = {};
        n.key = key;
        memcpy(n.name, name.ptr, name.len);
        n.frames = md_bitset_empty_set();
        md_array_push(x->slots, n, e->alloc);
        s = md_array_last(x->slots);
        s->version = next_version();
    }
    if (!md_bitset_equal(s->frames, *frames)) {
        md_bitset_free(&s->frames, e->alloc);
        s->frames = md_bitset_copy(*frames, e->alloc);
        s->computed = false;
        s->version = next_version();
    }
    return true;
}

void aggregate_release(Eval* e, uint64_t key, str_t name) {
    EvalV2* x = ev(e);
    if (!x) return;
    for (size_t i = 0; i < md_array_size(x->slots); ++i) {
        AggSlot* s = &x->slots[i];
        if (s->key != key || !str_eq(str_from_cstr(s->name), name)) continue;
        md_bitset_free(&s->frames, e->alloc);
        md_array_free(s->values, e->alloc);
        md_array_swap_back_and_pop(x->slots, i);
        return;
    }
}

// The slot computed for what the evaluation holds now: over its frames, once the run has finished
static bool slot_update(EvalV2* x, AggSlot* s, md_allocator_i* alloc) {
    const uint64_t runs = __atomic_load_n(&x->runs, __ATOMIC_ACQUIRE);
    if (__atomic_load_n(&x->running, __ATOMIC_ACQUIRE)) return s->computed;
    if (s->computed && s->seen_runs == runs) return true;
    // Kept per frame, so recomputed over the period's frames: nothing is evaluated again
    uint8_t* mask = (uint8_t*)md_alloc(md_get_heap_allocator(), x->num_frames ? x->num_frames : 1);
    memset(mask, 0, x->num_frames);
    md_bitset_iter_t it = md_bitset_iter(s->frames);
    while (md_bitset_iter_next(&it)) if (it.idx < x->num_frames) mask[it.idx] = 1;
    const size_t n = md_vial_eval_aggregate(x->eval, str_from_cstr(s->name), mask, nullptr, 0);
    bool ok = n > 0;
    if (ok) {
        md_array_resize(s->values, n, alloc);
        ok = md_vial_eval_aggregate(x->eval, str_from_cstr(s->name), mask, s->values, n) == n;
    }
    md_free(md_get_heap_allocator(), mask, x->num_frames ? x->num_frames : 1);
    if (ok) {
        s->computed = true;
        s->seen_runs = runs;
        s->version = next_version();
    }
    return ok;
}

static AggregateInfo slot_info(const EvalV2* x, const AggSlot* s) {
    AggregateInfo a;
    a.version = s->version;
    const uint32_t want = (uint32_t)md_bitset_count(s->frames);
    a.num_frames = s->computed ? want : 0;
    a.num_pending = s->computed ? 0 : want;
    a.num_values = s->computed ? md_array_size(s->values) : 0;
    (void)x;
    return a;
}

bool aggregate_info(AggregateInfo* out, const Eval* e, uint64_t key, str_t name) {
    EvalV2* x = ev(e);
    AggSlot* s = x ? find_slot(x, key, name) : nullptr;
    if (!s) return false;
    slot_update(x, s, e->alloc);
    *out = slot_info(x, s);
    return true;
}

bool aggregate_read(AggregateInfo* out, float* values, float* weights, size_t num_values, const Eval* e, uint64_t key, str_t name) {
    EvalV2* x = ev(e);
    AggSlot* s = x ? find_slot(x, key, name) : nullptr;
    if (!s || !slot_update(x, s, e->alloc)) return false;
    if (out) *out = slot_info(x, s);
    const size_t n = md_array_size(s->values);
    if (num_values != n) return false;
    if (values) memcpy(values, s->values, n * sizeof(float));
    // Every bin counts alike when bins are merged for drawing
    for (size_t i = 0; weights && i < n; ++i) weights[i] = 1.0f;
    return true;
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
