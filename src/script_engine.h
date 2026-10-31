#pragma once

// The script, behind one interface. viamd compiles, evaluates, plots, draws, completes and filters through these
// functions and their types alone; no other file names md_script_* or md_filter*. A program knows its language, and
// every function that takes one does what that language does: v1 is mdlib's md_script and md_filter
// (script_engine_v1.cpp), v2 is VIAL (script_engine_v2.cpp), and script_engine.cpp calls the one of the handle's or
// the call's language. A workspace's language picks which one compile and the filters use.
//
// The types here are viamd's own, not either language's: what the call sites read is the same whichever language
// stands behind them. Ranges are byte offsets into the source. Strings are views into what the program, the
// evaluation or the given allocator holds, and live as long as it does.

#include <core/md_str.h>
#include <core/md_bitset.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>
#include <core/md_unit.h>

#include <stdint.h>
#include <stddef.h>

struct md_allocator_i;
struct md_system_t;
struct md_system_state_t;
struct md_attribute_t;
struct md_attributes_t;
struct md_atom_pair_t;
struct md_writer_i;

namespace script {

enum class Language : uint32_t {
    V1 = 1,
    V2 = 2,
};

// Byte offsets into the source, [beg, end)
struct Range {
    int beg = 0;
    int end = 0;
};

// A part of a compiled program that can be drawn: an output, or a subexpression the editor marks. It names the part
// within the one compilation that made it and holds no pointer, so it can be kept, copied and passed between threads;
// against any other program it resolves to nothing. A zeroed one refers to nothing.
struct Ref {
    uint64_t program = 0;   // the compilation that made it, 0 for none
    uint32_t node = 0;
};

inline bool ref_empty(Ref r) { return r.program == 0; }
inline bool ref_equal(Ref a, Ref b) { return a.program == b.program && a.node == b.node; }

// ### COMPILE ###

struct Program;     // a compilation, always made, also of broken text; read-only once made

// A stored selection the script may name
struct NamedSet {
    str_t       name;
    md_bitset_t atoms;
};

struct CompileDesc {
    const md_system_t* sys = nullptr;
    const NamedSet*    sets = nullptr;
    size_t             num_sets = 0;
    str_t              run = {};        // v2: the run it is evaluated over, which attr() reads and frames are counted in
};

// What compiling said about a range of the source. context: atoms the message is about, or NULL.
struct Diagnostic {
    Range              range;
    str_t              text;
    const md_bitset_t* context = nullptr;
};

// A part of the source the editor marks for hover: its range, how deeply nested, a line of text, and what it draws
struct Mark {
    Range range;
    int   depth = 0;
    str_t text;
    Ref   ref;
};

// Empty src is not compiled: the program has the sets' names alone, no errors, and is not valid
Program* compile(Language lang, str_t src, const CompileDesc& desc, md_allocator_i* alloc);
void     program_free(Program* prog);

Language language(const Program* prog);
bool     valid(const Program* prog);            // compiled without errors; false for NULL
uint64_t fingerprint(const Program* prog);      // changes when what it computes does

size_t            num_errors(const Program* prog);
const Diagnostic* errors(const Program* prog);
size_t            num_warnings(const Program* prog);
const Diagnostic* warnings(const Program* prog);
size_t            num_marks(const Program* prog);
const Mark*       marks(const Program* prog);

// The errors and warnings, as a compiler prints them, with the source lines they are about; nothing if there are none.
// To see why a program did not compile, in a test or a tool: write_diagnostics(md_writer_stdout(), prog).
void write_diagnostics(md_writer_i* w, const Program* prog);

// Every name the program binds, for picking one that is not taken
size_t       num_identifiers(const Program* prog);
const str_t* identifiers(const Program* prog);

// ### OUTPUTS ###
// The values a program publishes, by name, in the order it binds them

enum OutputFlags : uint32_t {
    OUTPUT_NONE         = 0,
    OUTPUT_TEMPORAL     = 0x0001,   // a value per frame
    OUTPUT_DISTRIBUTION = 0x0002,   // an aggregate over frames: bins and weights
    OUTPUT_VOLUME       = 0x0004,   // an aggregate over frames: a grid
    OUTPUT_PERIODIC     = 0x0008,   // temporal values that wrap, such as angles; see script/<name>/period and /seam
    OUTPUT_SDF          = 0x0010,
};

size_t       num_outputs(const Program* prog);
const str_t* output_names(const Program* prog);
uint32_t     output_flags(const Program* prog, str_t name);         // OUTPUT_*; 0 for none of that name
Range        output_range(const Program* prog, str_t name);         // where it is named; empty for none
uint64_t     output_fingerprint(const Program* prog, str_t name);   // what it is computed from; 0 for none
Ref          output_ref(const Program* prog, str_t name);           // empty for none

// Whether ref was made by this compilation of prog, and so resolves in it
bool  ref_valid(const Program* prog, Ref ref);
str_t ref_name(const Program* prog, Ref ref);       // the name of what ref refers to, empty if none
int   ref_dim(const Program* prog, Ref ref);        // its outer extent, -1 when not known, 0 if it does not resolve

// ### EVALUATE ###
// An evaluation holds a program's outputs over the frames it has evaluated, published as attributes under script/
// (see attributes()). Made from an earlier one, it takes over every output computed alike, with the frames it holds.

struct Eval;

// v2 spreads frames over threads itself: through this, which calls fn(arg, i) for every i in [0, count) and returns when
// every call has returned. viamd passes its task pool.
typedef void (*ParallelFor)(void* user, size_t count, void (*fn)(void* arg, size_t i), void* arg);

struct EvalDesc {
    size_t             num_frames = 0;
    // v2: the system and run it evaluates over, and its threads
    const md_system_t* sys = nullptr;
    str_t              run = {};
    size_t             num_threads = 0;
    ParallelFor        parallel_for = nullptr;
    void*              parallel_user = nullptr;
};

// prev may be NULL; it is not changed and may be freed after. prog must outlive the evaluation.
Eval* eval_create(const Eval* prev, const Program* prog, const EvalDesc& desc, md_allocator_i* alloc);
void  eval_free(Eval* eval);

uint64_t eval_fingerprint(const Eval* eval);        // the fingerprint of the program it was made from
size_t   eval_num_outputs(const Eval* eval);
uint64_t eval_output_fingerprint(const Eval* eval, str_t name);
size_t   eval_num_frames(const Eval* eval);

// The frames an output holds, and the frames that remain to be evaluated; copies allocated with alloc
md_bitset_t eval_output_frames(const Eval* eval, str_t name, md_allocator_i* alloc);
md_bitset_t eval_pending_frames(Eval* eval, md_allocator_i* alloc);

// Evaluates frames of the run. Safe from several threads at once on one evaluation, also with overlapping frames.
bool eval_frames(Eval* eval, const Program* prog, const md_system_t* sys, str_t run, const uint32_t* frames, size_t num_frames);

// Whether one call of eval_frames evaluates every pending frame, whichever frames it is given, spreading them over the
// threads of EvalDesc itself (v2): call it once, from one task. Otherwise (v1) the pending frames are split among calls.
bool eval_whole(const Eval* eval);

// How far an evaluation is, 0 to 1, while eval_frames runs; 1 when it has nothing pending
float eval_progress(const Eval* eval);

// Stops an evaluation between frames; it stays stopped until reset
void eval_interrupt(Eval* eval);
void eval_reset_interrupt(Eval* eval);

// The outputs as attributes: script/<name>, with /range, /mean, /variance, /extent, /labels beside a temporal one,
// and /weight, /bin beside a distribution. Created with the evaluation, written in place as frames are evaluated.
const md_attributes_t* attributes(const Eval* eval);

// ### WHAT AN EDIT CHANGES ###
// How a program's outputs stand against an evaluation of an earlier one: what eval_create(eval, prog) would take
// over, and what would remain to be evaluated. Nothing is evaluated to find out.

enum class Status : uint32_t {
    New = 0,        // nothing to take over
    Changed,        // an output of its name is computed from something else
    Held,           // taken over from an output computed alike, with the frames it holds
};

struct Plan {
    Status   status = Status::New;
    str_t    from;              // Held: what it is taken over from, named otherwise when it was renamed
    uint32_t num_frames = 0;    // the frames it would hold
    uint32_t num_pending = 0;   // the frames it would have evaluated, for itself and its aggregates
};

// One plan per output of prog, in its order, and the frames that would remain, to out_pending when not NULL.
// eval may be NULL. False, and nothing written, if prog is not valid or cap is less than its number of outputs.
bool plan(Plan* out, size_t cap, md_bitset_t* out_pending, const Eval* eval, const Program* prog, size_t num_frames, md_allocator_i* alloc);

// ### AGGREGATES OVER CHOSEN FRAMES ###
// A distribution or a volume over the frames of a plot period, asked for under a key of the caller's. It moves toward
// the frames asked for as frames are evaluated, and always holds the aggregate over the frames it holds.

struct AggregateInfo {
    uint64_t version = 0;       // changes whenever what it holds does, never repeats, never 0
    uint32_t num_frames = 0;    // the frames it holds
    uint32_t num_pending = 0;   // the frames it has yet to add or take out
    size_t   num_values = 0;    // bins, or voxels
};

bool aggregate_request(Eval* eval, uint64_t key, str_t name, const md_bitset_t* frames);
void aggregate_release(Eval* eval, uint64_t key, str_t name);
bool aggregate_info(AggregateInfo* out, const Eval* eval, uint64_t key, str_t name);
bool aggregate_read(AggregateInfo* out, float* values, float* weights, size_t num_values, const Eval* eval, uint64_t key, str_t name);

// ### DRAWING ###
// What a part of a program looks like at a frame: viamd's own arrays, allocated with the drawing's allocator.

struct Vertex {
    vec3_t   pos;
    uint32_t color;
    vec3_t   normal;
    uint32_t picking_idx;
};

struct Text {
    vec3_t    pos;
    str_t     str;      // always set, drawable as it stands
    double    value;    // when unit is not none: the quantity, to format in units of the viewer's choosing
    md_unit_t unit;
};

struct Drawing {
    md_allocator_i*          alloc = nullptr;
    md_array(Vertex)         points = nullptr;
    md_array(Vertex)         lines = nullptr;
    md_array(Vertex)         triangles = nullptr;
    md_array(Text)           text = nullptr;
    struct {
        md_array(mat4_t)      matrices = nullptr;
        md_array(md_bitset_t) structures = nullptr;
        float                 extent = 0.0f;
    } sdf;
    md_array(md_bitset_t)    structures = nullptr;
    md_array(md_bitset_t)    groups = nullptr;     // the atoms of each element drawn, in order
    md_bitset_t              atom_mask = {};       // every atom drawn
};

enum DrawFlags : uint32_t {
    DRAW_DEFAULT  = 0,  // everything
    DRAW_GEOMETRY = 1,
    DRAW_ATOMS    = 2,
    DRAW_SDF      = 4,
    DRAW_TEXT     = 8,
};

struct DrawCtx {
    const Program*           prog = nullptr;
    const md_system_t*       sys = nullptr;
    const md_system_state_t* state = nullptr;
    // v2: an evaluation of prog, for values over the trajectory, and the run, whose frame state->frame is
    const Eval*              eval = nullptr;
    str_t                    run = {};
};

void drawing_init(Drawing* d, md_allocator_i* alloc);
void drawing_clear(Drawing* d);

// What ref refers to, resolved against ctx.prog; subidx picks one element of an array, -1 all. False if it does not resolve.
bool draw_ref(Drawing* out, Ref ref, int subidx, const DrawCtx& ctx, uint32_t flags);
// What an expression written in the program's language draws
bool draw_text(Drawing* out, str_t expr, const DrawCtx& ctx, uint32_t flags);
// The value of what ref refers to in a few lines, at ctx.state; returns its length, 0 if nothing
size_t value_text(char* buf, size_t cap, Ref ref, int subidx, const DrawCtx& ctx);

// ### EDITOR ###

enum class CompletionKind : uint32_t {
    Keyword = 1,
    Procedure,
    Constant,
    Variable,
    Parameter,
    Value,
};

struct Completion {
    str_t          label;   // what to show, and to match what has been typed against
    str_t          text;    // what replaces the range
    CompletionKind kind = CompletionKind::Keyword;
    bool           fits = false;    // has the type expected at the cursor
};

struct Completions {
    Range       range;      // the bytes a completion replaces
    str_t       prefix;     // what has been typed of it
    size_t      count = 0;
    Completion* items = nullptr;
};

size_t       num_keywords(Language lang);
const str_t* keywords(Language lang);
size_t       builtin_names(Language lang, str_t* out, size_t cap);    // how many there are; writes up to cap

// What can be written at the byte offset cursor of src; everything from alloc, a temporary arena
Completions complete(Language lang, str_t src, int cursor, const md_system_t* sys, md_allocator_i* alloc);

// The path attr() reads attr by, as completion offers it; empty if attr() cannot read it
str_t attr_path(Language lang, const md_attributes_t* table, const md_attribute_t* attr);

bool name_valid(Language lang, str_t name);

// ### FILTERS ###
// An expression in the workspace's language that picks atoms, evaluated within a context

struct FilterCtx {
    const md_system_t*       sys = nullptr;     // required
    const md_system_state_t* state = nullptr;   // NULL for the system's reference state
    const Program*           prog = nullptr;    // a program whose names it may use; may be NULL
    const md_bitset_t*       atoms = nullptr;   // the atoms it is evaluated within; NULL for the whole system
    str_t                    run = {};          // v2: the run attr() reads
};

struct FilterResult {
    bool  success = false;      // it picked atoms (possibly none)
    bool  is_dynamic = false;   // what it picks depends on the coordinates
    Range error_range;          // what the first error is about; empty when none, or about the whole
    char  error[256] = {};      // why it failed; empty on success
};

FilterResult filter(Language lang, md_bitset_t* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
FilterResult filter_sets(Language lang, md_array(md_bitset_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
FilterResult filter_bonds(Language lang, md_array(md_atom_pair_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
Completions  filter_complete(Language lang, str_t expr, int cursor, const FilterCtx& ctx, md_allocator_i* alloc);

}  // namespace script
