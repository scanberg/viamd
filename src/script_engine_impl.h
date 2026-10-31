#pragma once

// What the two implementations of script_engine.h share: the handles, which know their language, and each language's
// functions, script::v1:: and script::v2::, with the interface's signatures. script_engine.cpp calls the one of the
// handle's or the call's language.

#include "script_engine.h"

struct md_script_ir_t;
struct md_script_eval_t;

namespace script {

struct ProgramV2;   // script_engine_v2.cpp
struct EvalV2;

struct Program {
    Language                lang;
    md_allocator_i*         alloc;
    md_array(Diagnostic)    errs;
    md_array(Diagnostic)    warns;
    md_array(Mark)          marks;
    md_script_ir_t*         v1;
    ProgramV2*              v2;
};

struct Eval {
    Language                lang;
    md_allocator_i*         alloc;
    md_script_eval_t*       v1;
    EvalV2*                 v2;
};

namespace v1 {
Program* compile(Language lang, str_t src, const CompileDesc& desc, md_allocator_i* alloc);
void program_free(Program* prog);
Language language(const Program* prog);
bool valid(const Program* prog);
uint64_t fingerprint(const Program* prog);
size_t num_errors(const Program* prog);
const Diagnostic* errors(const Program* prog);
size_t num_warnings(const Program* prog);
const Diagnostic* warnings(const Program* prog);
size_t num_marks(const Program* prog);
const Mark* marks(const Program* prog);
void write_diagnostics(md_writer_i* w, const Program* prog);
size_t num_identifiers(const Program* prog);
const str_t* identifiers(const Program* prog);
size_t num_outputs(const Program* prog);
const str_t* output_names(const Program* prog);
uint32_t output_flags(const Program* prog, str_t name);
Range output_range(const Program* prog, str_t name);
uint64_t output_fingerprint(const Program* prog, str_t name);
Ref output_ref(const Program* prog, str_t name);
bool ref_valid(const Program* prog, Ref ref);
str_t ref_name(const Program* prog, Ref ref);
int ref_dim(const Program* prog, Ref ref);
Eval* eval_create(const Eval* prev, const Program* prog, const EvalDesc& desc, md_allocator_i* alloc);
void eval_free(Eval* eval);
uint64_t eval_fingerprint(const Eval* eval);
size_t eval_num_outputs(const Eval* eval);
uint64_t eval_output_fingerprint(const Eval* eval, str_t name);
size_t eval_num_frames(const Eval* eval);
md_bitset_t eval_output_frames(const Eval* eval, str_t name, md_allocator_i* alloc);
md_bitset_t eval_pending_frames(Eval* eval, md_allocator_i* alloc);
bool eval_frames(Eval* eval, const Program* prog, const md_system_t* sys, str_t run, const uint32_t* frames, size_t num_frames);
bool eval_whole(const Eval* eval);
float eval_progress(const Eval* eval);
EvalErrors eval_errors(const Eval* eval);
void eval_interrupt(Eval* eval);
void eval_reset_interrupt(Eval* eval);
const md_attributes_t* attributes(const Eval* eval);
bool plan(Plan* out, size_t cap, md_bitset_t* out_pending, const Eval* eval, const Program* prog, size_t num_frames, md_allocator_i* alloc);
bool aggregate_request(Eval* eval, uint64_t key, str_t name, const md_bitset_t* frames);
void aggregate_release(Eval* eval, uint64_t key, str_t name);
bool aggregate_info(AggregateInfo* out, const Eval* eval, uint64_t key, str_t name);
bool aggregate_read(AggregateInfo* out, float* values, float* weights, size_t num_values, const Eval* eval, uint64_t key, str_t name);
void drawing_init(Drawing* d, md_allocator_i* alloc);
void drawing_clear(Drawing* d);
bool draw_ref(Drawing* out, Ref ref, int subidx, const DrawCtx& ctx, uint32_t flags);
bool draw_text(Drawing* out, str_t expr, const DrawCtx& ctx, uint32_t flags);
size_t value_text(char* buf, size_t cap, Ref ref, int subidx, const DrawCtx& ctx);
size_t num_keywords(Language lang);
const str_t* keywords(Language lang);
size_t builtin_names(Language lang, str_t* out, size_t cap);
Completions complete(Language lang, str_t src, int cursor, const md_system_t* sys, md_allocator_i* alloc);
str_t attr_path(Language lang, const md_attributes_t* table, const md_attribute_t* attr);
bool name_valid(Language lang, str_t name);
FilterResult filter(Language lang, md_bitset_t* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
FilterResult filter_sets(Language lang, md_array(md_bitset_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
FilterResult filter_bonds(Language lang, md_array(md_atom_pair_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
Completions filter_complete(Language lang, str_t expr, int cursor, const FilterCtx& ctx, md_allocator_i* alloc);
}  // namespace v1

namespace v2 {
Program* compile(Language lang, str_t src, const CompileDesc& desc, md_allocator_i* alloc);
void program_free(Program* prog);
Language language(const Program* prog);
bool valid(const Program* prog);
uint64_t fingerprint(const Program* prog);
size_t num_errors(const Program* prog);
const Diagnostic* errors(const Program* prog);
size_t num_warnings(const Program* prog);
const Diagnostic* warnings(const Program* prog);
size_t num_marks(const Program* prog);
const Mark* marks(const Program* prog);
void write_diagnostics(md_writer_i* w, const Program* prog);
size_t num_identifiers(const Program* prog);
const str_t* identifiers(const Program* prog);
size_t num_outputs(const Program* prog);
const str_t* output_names(const Program* prog);
uint32_t output_flags(const Program* prog, str_t name);
Range output_range(const Program* prog, str_t name);
uint64_t output_fingerprint(const Program* prog, str_t name);
Ref output_ref(const Program* prog, str_t name);
bool ref_valid(const Program* prog, Ref ref);
str_t ref_name(const Program* prog, Ref ref);
int ref_dim(const Program* prog, Ref ref);
Eval* eval_create(const Eval* prev, const Program* prog, const EvalDesc& desc, md_allocator_i* alloc);
void eval_free(Eval* eval);
uint64_t eval_fingerprint(const Eval* eval);
size_t eval_num_outputs(const Eval* eval);
uint64_t eval_output_fingerprint(const Eval* eval, str_t name);
size_t eval_num_frames(const Eval* eval);
md_bitset_t eval_output_frames(const Eval* eval, str_t name, md_allocator_i* alloc);
md_bitset_t eval_pending_frames(Eval* eval, md_allocator_i* alloc);
bool eval_frames(Eval* eval, const Program* prog, const md_system_t* sys, str_t run, const uint32_t* frames, size_t num_frames);
bool eval_whole(const Eval* eval);
float eval_progress(const Eval* eval);
EvalErrors eval_errors(const Eval* eval);
void eval_interrupt(Eval* eval);
void eval_reset_interrupt(Eval* eval);
const md_attributes_t* attributes(const Eval* eval);
bool plan(Plan* out, size_t cap, md_bitset_t* out_pending, const Eval* eval, const Program* prog, size_t num_frames, md_allocator_i* alloc);
bool aggregate_request(Eval* eval, uint64_t key, str_t name, const md_bitset_t* frames);
void aggregate_release(Eval* eval, uint64_t key, str_t name);
bool aggregate_info(AggregateInfo* out, const Eval* eval, uint64_t key, str_t name);
bool aggregate_read(AggregateInfo* out, float* values, float* weights, size_t num_values, const Eval* eval, uint64_t key, str_t name);
bool draw_ref(Drawing* out, Ref ref, int subidx, const DrawCtx& ctx, uint32_t flags);
bool draw_text(Drawing* out, str_t expr, const DrawCtx& ctx, uint32_t flags);
size_t value_text(char* buf, size_t cap, Ref ref, int subidx, const DrawCtx& ctx);
size_t num_keywords(Language lang);
const str_t* keywords(Language lang);
size_t builtin_names(Language lang, str_t* out, size_t cap);
Completions complete(Language lang, str_t src, int cursor, const md_system_t* sys, md_allocator_i* alloc);
str_t attr_path(Language lang, const md_attributes_t* table, const md_attribute_t* attr);
bool name_valid(Language lang, str_t name);
FilterResult filter(Language lang, md_bitset_t* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
FilterResult filter_sets(Language lang, md_array(md_bitset_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
FilterResult filter_bonds(Language lang, md_array(md_atom_pair_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc);
Completions filter_complete(Language lang, str_t expr, int cursor, const FilterCtx& ctx, md_allocator_i* alloc);
}  // namespace v2

}  // namespace script
