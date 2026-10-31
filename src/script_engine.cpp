// The script interface: each call goes to the implementation of its language (script_engine_impl.h).

#include "script_engine_impl.h"

namespace script {

static inline Language lang_of(const Program* p) { return p ? p->lang : Language::V1; }
static inline Language lang_of(const Eval* e) { return e ? e->lang : Language::V1; }

Program* compile(Language lang, str_t src, const CompileDesc& desc, md_allocator_i* alloc) {
    return lang == Language::V2 ? v2::compile(lang, src, desc, alloc) : v1::compile(lang, src, desc, alloc);
}
void program_free(Program* prog) {
    if (lang_of(prog) == Language::V2) v2::program_free(prog);
    else v1::program_free(prog);
}
Language language(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::language(prog) : v1::language(prog);
}
bool valid(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::valid(prog) : v1::valid(prog);
}
uint64_t fingerprint(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::fingerprint(prog) : v1::fingerprint(prog);
}
size_t num_errors(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::num_errors(prog) : v1::num_errors(prog);
}
const Diagnostic* errors(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::errors(prog) : v1::errors(prog);
}
size_t num_warnings(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::num_warnings(prog) : v1::num_warnings(prog);
}
const Diagnostic* warnings(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::warnings(prog) : v1::warnings(prog);
}
size_t num_marks(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::num_marks(prog) : v1::num_marks(prog);
}
const Mark* marks(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::marks(prog) : v1::marks(prog);
}
void write_diagnostics(md_writer_i* w, const Program* prog) {
    if (lang_of(prog) == Language::V2) v2::write_diagnostics(w, prog);
    else v1::write_diagnostics(w, prog);
}
size_t num_identifiers(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::num_identifiers(prog) : v1::num_identifiers(prog);
}
const str_t* identifiers(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::identifiers(prog) : v1::identifiers(prog);
}
size_t num_outputs(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::num_outputs(prog) : v1::num_outputs(prog);
}
const str_t* output_names(const Program* prog) {
    return lang_of(prog) == Language::V2 ? v2::output_names(prog) : v1::output_names(prog);
}
uint32_t output_flags(const Program* prog, str_t name) {
    return lang_of(prog) == Language::V2 ? v2::output_flags(prog, name) : v1::output_flags(prog, name);
}
Range output_range(const Program* prog, str_t name) {
    return lang_of(prog) == Language::V2 ? v2::output_range(prog, name) : v1::output_range(prog, name);
}
uint64_t output_fingerprint(const Program* prog, str_t name) {
    return lang_of(prog) == Language::V2 ? v2::output_fingerprint(prog, name) : v1::output_fingerprint(prog, name);
}
Ref output_ref(const Program* prog, str_t name) {
    return lang_of(prog) == Language::V2 ? v2::output_ref(prog, name) : v1::output_ref(prog, name);
}
bool ref_valid(const Program* prog, Ref ref) {
    return lang_of(prog) == Language::V2 ? v2::ref_valid(prog, ref) : v1::ref_valid(prog, ref);
}
str_t ref_name(const Program* prog, Ref ref) {
    return lang_of(prog) == Language::V2 ? v2::ref_name(prog, ref) : v1::ref_name(prog, ref);
}
int ref_dim(const Program* prog, Ref ref) {
    return lang_of(prog) == Language::V2 ? v2::ref_dim(prog, ref) : v1::ref_dim(prog, ref);
}
Eval* eval_create(const Eval* prev, const Program* prog, const EvalDesc& desc, md_allocator_i* alloc) {
    return lang_of(prog) == Language::V2 ? v2::eval_create(prev, prog, desc, alloc) : v1::eval_create(prev, prog, desc, alloc);
}
void eval_free(Eval* eval) {
    if (lang_of(eval) == Language::V2) v2::eval_free(eval);
    else v1::eval_free(eval);
}
uint64_t eval_fingerprint(const Eval* eval) {
    return lang_of(eval) == Language::V2 ? v2::eval_fingerprint(eval) : v1::eval_fingerprint(eval);
}
size_t eval_num_outputs(const Eval* eval) {
    return lang_of(eval) == Language::V2 ? v2::eval_num_outputs(eval) : v1::eval_num_outputs(eval);
}
uint64_t eval_output_fingerprint(const Eval* eval, str_t name) {
    return lang_of(eval) == Language::V2 ? v2::eval_output_fingerprint(eval, name) : v1::eval_output_fingerprint(eval, name);
}
size_t eval_num_frames(const Eval* eval) {
    return lang_of(eval) == Language::V2 ? v2::eval_num_frames(eval) : v1::eval_num_frames(eval);
}
md_bitset_t eval_output_frames(const Eval* eval, str_t name, md_allocator_i* alloc) {
    return lang_of(eval) == Language::V2 ? v2::eval_output_frames(eval, name, alloc) : v1::eval_output_frames(eval, name, alloc);
}
md_bitset_t eval_pending_frames(Eval* eval, md_allocator_i* alloc) {
    return lang_of(eval) == Language::V2 ? v2::eval_pending_frames(eval, alloc) : v1::eval_pending_frames(eval, alloc);
}
bool eval_frames(Eval* eval, const Program* prog, const md_system_t* sys, str_t run, const uint32_t* frames, size_t num_frames) {
    return lang_of(prog) == Language::V2 ? v2::eval_frames(eval, prog, sys, run, frames, num_frames) : v1::eval_frames(eval, prog, sys, run, frames, num_frames);
}
bool eval_whole(const Eval* eval) {
    return lang_of(eval) == Language::V2 ? v2::eval_whole(eval) : v1::eval_whole(eval);
}
float eval_progress(const Eval* eval) {
    return lang_of(eval) == Language::V2 ? v2::eval_progress(eval) : v1::eval_progress(eval);
}
EvalErrors eval_errors(const Eval* eval) {
    return lang_of(eval) == Language::V2 ? v2::eval_errors(eval) : v1::eval_errors(eval);
}
void eval_interrupt(Eval* eval) {
    if (lang_of(eval) == Language::V2) v2::eval_interrupt(eval);
    else v1::eval_interrupt(eval);
}
void eval_reset_interrupt(Eval* eval) {
    if (lang_of(eval) == Language::V2) v2::eval_reset_interrupt(eval);
    else v1::eval_reset_interrupt(eval);
}
const md_attributes_t* attributes(const Eval* eval) {
    return lang_of(eval) == Language::V2 ? v2::attributes(eval) : v1::attributes(eval);
}
bool plan(Plan* out, size_t cap, md_bitset_t* out_pending, const Eval* eval, const Program* prog, size_t num_frames, md_allocator_i* alloc) {
    return lang_of(prog) == Language::V2 ? v2::plan(out, cap, out_pending, eval, prog, num_frames, alloc) : v1::plan(out, cap, out_pending, eval, prog, num_frames, alloc);
}
bool aggregate_request(Eval* eval, uint64_t key, str_t name, const md_bitset_t* frames) {
    return lang_of(eval) == Language::V2 ? v2::aggregate_request(eval, key, name, frames) : v1::aggregate_request(eval, key, name, frames);
}
void aggregate_release(Eval* eval, uint64_t key, str_t name) {
    if (lang_of(eval) == Language::V2) v2::aggregate_release(eval, key, name);
    else v1::aggregate_release(eval, key, name);
}
bool aggregate_info(AggregateInfo* out, const Eval* eval, uint64_t key, str_t name) {
    return lang_of(eval) == Language::V2 ? v2::aggregate_info(out, eval, key, name) : v1::aggregate_info(out, eval, key, name);
}
bool aggregate_read(AggregateInfo* out, float* values, float* weights, size_t num_values, const Eval* eval, uint64_t key, str_t name) {
    return lang_of(eval) == Language::V2 ? v2::aggregate_read(out, values, weights, num_values, eval, key, name) : v1::aggregate_read(out, values, weights, num_values, eval, key, name);
}
void drawing_init(Drawing* d, md_allocator_i* alloc) { v1::drawing_init(d, alloc); }
void drawing_clear(Drawing* d) { v1::drawing_clear(d); }
bool draw_ref(Drawing* out, Ref ref, int subidx, const DrawCtx& ctx, uint32_t flags) {
    return lang_of(ctx.prog) == Language::V2 ? v2::draw_ref(out, ref, subidx, ctx, flags) : v1::draw_ref(out, ref, subidx, ctx, flags);
}
bool draw_text(Drawing* out, str_t expr, const DrawCtx& ctx, uint32_t flags) {
    return lang_of(ctx.prog) == Language::V2 ? v2::draw_text(out, expr, ctx, flags) : v1::draw_text(out, expr, ctx, flags);
}
size_t value_text(char* buf, size_t cap, Ref ref, int subidx, const DrawCtx& ctx) {
    return lang_of(ctx.prog) == Language::V2 ? v2::value_text(buf, cap, ref, subidx, ctx) : v1::value_text(buf, cap, ref, subidx, ctx);
}
size_t num_keywords(Language lang) {
    return lang == Language::V2 ? v2::num_keywords(lang) : v1::num_keywords(lang);
}
const str_t* keywords(Language lang) {
    return lang == Language::V2 ? v2::keywords(lang) : v1::keywords(lang);
}
size_t builtin_names(Language lang, str_t* out, size_t cap) {
    return lang == Language::V2 ? v2::builtin_names(lang, out, cap) : v1::builtin_names(lang, out, cap);
}
Completions complete(Language lang, str_t src, int cursor, const md_system_t* sys, md_allocator_i* alloc) {
    return lang == Language::V2 ? v2::complete(lang, src, cursor, sys, alloc) : v1::complete(lang, src, cursor, sys, alloc);
}
str_t attr_path(Language lang, const md_attributes_t* table, const md_attribute_t* attr) {
    return lang == Language::V2 ? v2::attr_path(lang, table, attr) : v1::attr_path(lang, table, attr);
}
bool name_valid(Language lang, str_t name) {
    return lang == Language::V2 ? v2::name_valid(lang, name) : v1::name_valid(lang, name);
}
FilterResult filter(Language lang, md_bitset_t* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc) {
    return lang == Language::V2 ? v2::filter(lang, dst, expr, ctx, alloc) : v1::filter(lang, dst, expr, ctx, alloc);
}
FilterResult filter_sets(Language lang, md_array(md_bitset_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc) {
    return lang == Language::V2 ? v2::filter_sets(lang, dst, expr, ctx, alloc) : v1::filter_sets(lang, dst, expr, ctx, alloc);
}
FilterResult filter_bonds(Language lang, md_array(md_atom_pair_t)* dst, str_t expr, const FilterCtx& ctx, md_allocator_i* alloc) {
    return lang == Language::V2 ? v2::filter_bonds(lang, dst, expr, ctx, alloc) : v1::filter_bonds(lang, dst, expr, ctx, alloc);
}
Completions filter_complete(Language lang, str_t expr, int cursor, const FilterCtx& ctx, md_allocator_i* alloc) {
    return lang == Language::V2 ? v2::filter_complete(lang, expr, cursor, ctx, alloc) : v1::filter_complete(lang, expr, cursor, ctx, alloc);
}

}  // namespace script
