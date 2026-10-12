// Throwaway: compiles and evaluates the script of a .via workspace on the aspirin data, with mdlib.
#include <md_gro.h>
#include <md_xtc.h>
#include <md_edr.h>
#include <md_script.h>
#include <md_util.h>
#include <md_system.h>
#include <md_attributes.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <fstream>
#include <sstream>

int main(int argc, char** argv) {
    if (argc < 3) return 1;
    std::ifstream in(argv[2]);
    std::stringstream ss; ss << in.rdbuf();
    std::string via = ss.str();
    size_t a = via.find("Text=\"\"\"");
    size_t b = via.find("\"\"\"", a + 8);
    std::string src = via.substr(a + 8, b - (a + 8));

    std::string dir = argv[1];
    std::string gro = dir + "/aspirin-phospholipase.gro", xtc = dir + "/aspirin-phospholipase.xtc", edr = dir + "/aspirin-phospholipase.edr";
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {};
    sys.alloc = arena;
    md_system_state_t st = {};
    st.alloc = arena;
    if (!md_gro_system_init_from_file(&sys, &st, str_from_cstr(gro.c_str()))) { printf("gro failed\n"); return 2; }
    md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL);
    const str_t run = STR_LIT("run/aspirin-phospholipase");
    if (!md_xtc_system_publish_run(&sys, str_from_cstr(xtc.c_str()), run, MD_RUN_FLAG_DISABLE_CACHE_WRITE)) { printf("xtc failed\n"); return 3; }
    if (!md_edr_system_supplement_from_file(&sys, str_from_cstr(edr.c_str()), run)) { printf("edr failed\n"); return 4; }

    char buf[512];
    const md_attribute_t* time = md_attributes_find(&sys.attributes, md_run_path(buf, sizeof(buf), run, STR_LIT("time")));
    const uint32_t F = time ? time->format.shape[0] : 0;
    printf("frames %u\n", F);

    md_script_ir_t* ir = md_script_ir_create(arena);
    bool ok = md_script_ir_compile_from_source(ir, str_t{src.c_str(), src.size()}, &sys, NULL);
    printf("compile %d valid %d errors %zu\n", ok, md_script_ir_valid(ir), md_script_ir_num_errors(ir));
    const md_log_token_t* errs = md_script_ir_errors(ir);
    for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) printf("  error: %.*s\n", (int)errs[i].text.len, errs[i].text.ptr);
    if (!md_script_ir_valid(ir)) return 5;

    md_script_eval_t* ev = md_script_eval_create(F, ir, arena);
    printf("eval %d\n", md_script_eval_frame_range(ev, ir, &sys, run, 0, F));
    size_t n = md_script_ir_property_count(ir);
    const str_t* names = md_script_ir_property_names(ir);
    const uint32_t frames[] = {0, 4, 30, 120, 200, 288, 330, 450, 600};
    for (size_t i = 0; i < n; ++i) {
        std::string path = "script/" + std::string(names[i].ptr, names[i].len);
        const md_attribute_t* at = md_attributes_find(md_script_eval_attributes(ev), str_from_cstr(path.c_str()));
        printf("%-14s", path.c_str() + 7);
        if (!at) { printf(" (no attribute)\n"); continue; }
        for (uint32_t f : frames) {
            double v[4] = {};
            md_attribute_slice_t row = md_attribute_slice_1(f);
            size_t k = md_attribute_extract_f64(v, 4, at, row, md_unit_none());
            printf(" %9.2f", k ? v[0] : -9999.0);
        }
        printf("\n");
    }
    return 0;
}
