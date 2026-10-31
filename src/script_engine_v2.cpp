// The script interface (script_engine.h) for v2: VIAL, mdlib's md_vial.

#include "script_engine_impl.h"

#include <md_vial.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_array.h>

#include <string.h>
#include <stdio.h>

namespace script {
namespace v2 {

static FilterResult not_yet() {
    FilterResult r;
    snprintf(r.error, sizeof(r.error), "not yet in VIAL");
    return r;
}

Program* compile(Language lang,  str_t src,  const CompileDesc& desc,  md_allocator_i* alloc) {
    (void)lang; (void)src; (void)desc; (void)alloc;
    return nullptr;
}

void program_free(Program* prog) {
    (void)prog;
    
}

Language language(const Program* prog) {
    (void)prog;
    return Language::V2;
}

bool valid(const Program* prog) {
    (void)prog;
    return false;
}

uint64_t fingerprint(const Program* prog) {
    (void)prog;
    return 0;
}

size_t num_errors(const Program* prog) {
    (void)prog;
    return 0;
}

const Diagnostic* errors(const Program* prog) {
    (void)prog;
    return nullptr;
}

size_t num_warnings(const Program* prog) {
    (void)prog;
    return 0;
}

const Diagnostic* warnings(const Program* prog) {
    (void)prog;
    return nullptr;
}

size_t num_marks(const Program* prog) {
    (void)prog;
    return 0;
}

const Mark* marks(const Program* prog) {
    (void)prog;
    return nullptr;
}

void write_diagnostics(md_writer_i* w,  const Program* prog) {
    (void)w; (void)prog;
    
}

size_t num_identifiers(const Program* prog) {
    (void)prog;
    return 0;
}

const str_t* identifiers(const Program* prog) {
    (void)prog;
    return nullptr;
}

size_t num_outputs(const Program* prog) {
    (void)prog;
    return 0;
}

const str_t* output_names(const Program* prog) {
    (void)prog;
    return nullptr;
}

uint32_t output_flags(const Program* prog,  str_t name) {
    (void)prog; (void)name;
    return 0;
}

Range output_range(const Program* prog,  str_t name) {
    (void)prog; (void)name;
    return Range{};
}

uint64_t output_fingerprint(const Program* prog,  str_t name) {
    (void)prog; (void)name;
    return 0;
}

Ref output_ref(const Program* prog,  str_t name) {
    (void)prog; (void)name;
    return Ref{};
}

bool ref_valid(const Program* prog,  Ref ref) {
    (void)prog; (void)ref;
    return false;
}

str_t ref_name(const Program* prog,  Ref ref) {
    (void)prog; (void)ref;
    return str_t{};
}

int ref_dim(const Program* prog,  Ref ref) {
    (void)prog; (void)ref;
    return 0;
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

size_t num_keywords(Language lang) {
    (void)lang;
    return 0;
}

const str_t* keywords(Language lang) {
    (void)lang;
    return nullptr;
}

size_t builtin_names(Language lang,  str_t* out,  size_t cap) {
    (void)lang; (void)out; (void)cap;
    return 0;
}

Completions complete(Language lang,  str_t src,  int cursor,  const md_system_t* sys,  md_allocator_i* alloc) {
    (void)lang; (void)src; (void)cursor; (void)sys; (void)alloc;
    return Completions{};
}

str_t attr_path(Language lang,  const md_attributes_t* table,  const md_attribute_t* attr) {
    (void)lang; (void)table; (void)attr;
    return str_t{};
}

bool name_valid(Language lang,  str_t name) {
    (void)lang; (void)name;
    return false;
}

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
