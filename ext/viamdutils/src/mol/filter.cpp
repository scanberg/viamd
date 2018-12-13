#include "filter.h"
#include <core/common.h>
#include <core/hash.h>
#include <core/log.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>

namespace filter {
static DynamicArray<FilterCommand> filter_commands;

static bool is_modifier(CString str) {
    if (compare(str, "and", true)) return true;
    if (compare(str, "or", true)) return true;
    return false;
}

static bool is_keyword(CString str) {
    if (compare(str, "and", true)) return true;
    if (compare(str, "or", true)) return true;
    if (compare(str, "not", true)) return true;
    return false;
}

int32 count_parentheses(CString str) {
    int beg_parentheses_count = 0;
    int end_parentheses_count = 0;

    for (int64 i = 0; i < str.count; i++) {
        if (str[i] == '(') beg_parentheses_count++;
        if (str[i] == ')') end_parentheses_count++;
    }
    return beg_parentheses_count - end_parentheses_count;
}

CString extract_parenthesis(CString str) {
    const char* beg = str.beg();
    while (beg != str.end() && *beg != '(') beg++;
    if (beg == str.end()) return {};

    const char* end = beg + 1;
    int count = 1;
    while (end++ != str.end() && count > 0) {
        if (*end == '(') count++;
        if (*end == ')') count--;
    }
    if (end == str.end()) return {};
    return CString(beg, end);
}

DynamicArray<CString> extract_chunks(CString str) {
    DynamicArray<CString> chunks;

    const char* beg = str.beg();
    while (beg != str.end()) {
        if (*beg == '(') {
            CString par = extract_parenthesis(CString(beg, str.end()));
            chunks.push_back({par.beg(), par.end()});  // Exclude actual parentheses
            beg = par.end();
        } else if (*beg != ' ') {
            const char* end = beg;
            while (end != str.end() && *end != ' ') end++;
            chunks.push_back(CString(beg, end));
            beg = end;
        } else
            beg++;
    }

    DynamicArray<CString> big_chunks;

    CString* chunk = chunks.beg();
    while (chunk != chunks.end()) {
        if (chunk->front() == '(') {
            big_chunks.push_back(*chunk);
            chunk++;
        } else if (is_keyword(*chunk)) {
            big_chunks.push_back(*chunk);
            chunk++;
        } else {
            CString* beg_chunk = chunk;
            CString* end_chunk = chunk + 1;
            while (end_chunk != chunks.end() && !is_modifier(*end_chunk) && end_chunk->front() != '(') end_chunk++;
            big_chunks.push_back({beg_chunk->beg(), (end_chunk - 1)->end()});
            chunk = end_chunk;
        }
    }

    return big_chunks;
}

FilterCommand* find_filter_command(CString command) {
    for (auto& f : filter_commands) {
        if (compare(command, f.keyword)) return &f;
    }
    return nullptr;
}

void combine_mask_and(Array<bool> dst, Array<bool> src_a, Array<bool> src_b, bool state_not) {
    if (state_not) {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] & !src_b[i];
    } else {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] & src_b[i];
    }
}

void combine_mask_or(Array<bool> dst, Array<bool> src_a, Array<bool> src_b, bool state_not) {
    if (state_not) {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] | !src_b[i];
    } else {
        for (int i = 0; i < dst.count; i++) dst[i] = src_a[i] | src_b[i];
    }
}

bool internal_filter_mask(Array<bool> mask, const MoleculeDynamic& dyn, CString filter) {
    DynamicArray<CString> chunks = extract_chunks(filter);
    DynamicArray<bool> chunk_mask(mask.count);

    bool state_and = true;
    bool state_or = false;
    bool state_not = false;

    for (const auto& chunk : chunks) {
        if (compare(chunk, "and", true)) {
            state_and = true;
        } else if (compare(chunk, "or", true)) {
            state_or = true;
        } else if (compare(chunk, "not", true)) {
            state_not = true;
        } else {
            if (chunk.front() == '(') {
                ASSERT(chunk.back() == ')');
                if (!internal_filter_mask(chunk_mask, dyn, CString(chunk.beg() + 1, chunk.end() - 1))) return false;
            } else {
                auto tokens = ctokenize(chunk);
                auto cmd = find_filter_command(tokens[0]);
                if (!cmd) {
                    StringBuffer<32> buf = tokens[0];
                    LOG_ERROR("Could not match command: '%s'\n", buf.beg());
                    return false;
                }

                auto args = tokens.sub_array(1);

                while (args.count > 0 && compare(args[0], "not", true)) {
                    state_not = !state_not;
                    args = args.sub_array(1);
                }

                if (!cmd->func(chunk_mask, dyn, args)) {
                    StringBuffer<32> buf = tokens[0];
                    LOG_ERROR("Could not parse command: '%s' with arguments: ", buf.beg());
                    for (const auto& arg : args) {
                        buf = arg;
                        printf("'%s'", buf.beg());
                    }
                    return false;
                }
            }

            if (state_and)
                combine_mask_and(mask, mask, chunk_mask, state_not);
            else if (state_or)
                combine_mask_or(mask, mask, chunk_mask, state_not);

            state_and = false;
            state_or = false;
            state_not = false;
        }

        if (state_and && state_or) {
            LOG_ERROR("Cannot use both 'and' and 'or' to combine filter options\n");
            return false;
        }
    }

    return true;
}

bool compute_filter_mask(Array<bool> mask, const MoleculeDynamic& dyn, CString filter) {
    ASSERT(dyn.molecule.atom.count == mask.count);

    if (count_parentheses(filter) != 0) {
        LOG_ERROR("Unmatched parentheses\n");
        return false;
    }

    memset(mask.data, 1, mask.count);
    return internal_filter_mask(mask, dyn, filter);
}

void filter_colors(Array<uint32> colors, Array<bool> mask) {
    ASSERT(colors.count == mask.count);
    for (int i = 0; i < colors.count; i++) {
        if (mask[i])
            colors[i] |= 0xff000000;
        else
            colors[i] &= ~0xff000000;
    }
}

void initialize() {

    /*
            all
            water
            aminoacid
            backbone?
            protein

            name
            element
            atomicnumber
            atom
            residue
            resname
            resid
            chain
            chainid
    */

    auto filter_amino_acid = [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString>) {
        memset(mask.data, 0, mask.count);
        for (const auto& res : dyn.molecule.residues) {
            if (is_amino_acid(res)) {
                memset(mask.data + res.atom_idx.beg, 1, res.atom_idx.end - res.atom_idx.beg);
            }
        }
        return true;
    };

    filter_commands.push_back({"all", [](Array<bool> mask, const MoleculeDynamic&, Array<const CString>) {
                                   memset(mask.data, 1, mask.count);
                                   return true;
                               }});
    filter_commands.push_back({"water", [](Array<bool>, const MoleculeDynamic&, Array<const CString>) { return true; }});  // NOT DONE
    filter_commands.push_back({"aminoacid", filter_amino_acid});
    filter_commands.push_back({"backbone", [](Array<bool>, const MoleculeDynamic&, Array<const CString>) { return true; }});  // NOT DONE
    filter_commands.push_back({"protein", filter_amino_acid});
    filter_commands.push_back({"dna", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString>) {
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (const auto& res : dyn.molecule.residues) {
                                       if (is_dna(res)) {
                                           memset(mask.data + res.atom_idx.beg, 1, res.atom_idx.end - res.atom_idx.beg);
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"name", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   if (args.count == 0) return false;

                                   for (int i = 0; i < dyn.molecule.atom.count; i++) {
                                       mask[i] = false;
                                       for (const auto& arg : args) {
                                           if (compare(dyn.molecule.atom.labels[i], arg)) {
                                               mask[i] = true;
                                               break;
                                           }
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"label", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   return find_filter_command("label")->func(mask, dyn, args);
                               }});

    filter_commands.push_back({"element", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   Array<Element> elements = {(Element*)(TMP_MALLOC(args.count * sizeof(Element))), args.count};
                                   for (int i = 0; i < elements.count; i++) {
                                       elements[i] = element::get_from_string(args[i]);
                                       if (elements[i] == Element::Unknown) return false;
                                   }

                                   for (int i = 0; i < dyn.molecule.atom.count; i++) {
                                       mask[i] = false;
                                       for (const auto& ele : elements) {
                                           if (dyn.molecule.atom.elements[i] == ele) {
                                               mask[i] = true;
                                               break;
                                           }
                                       }
                                   }
                                   TMP_FREE(elements.data);
                                   return true;
                               }});

    filter_commands.push_back({"atomicnumber", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   for (int i = 0; i < dyn.molecule.atom.count; i++) {
                                       int atomnr = (int)dyn.molecule.atom.elements[i];
                                       mask[i] = false;
                                       for (auto range : ranges) {
                                           if (range.x <= atomnr && atomnr <= range.y) {
                                               mask[i] = true;
                                               break;
                                           }
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"atom", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (auto range : ranges) {
                                       range.x = math::clamp(range.x - 1, 0, (int32)dyn.molecule.atom.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)dyn.molecule.atom.count - 1);
                                       if (range.x == range.y)
                                           mask[range.x] = true;
                                       else
                                           memset(mask.data + range.x, 1, range.y - range.x + 1);
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"residue", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (auto range : ranges) {
                                       range.x = math::clamp(range.x, 0, (int32)dyn.molecule.residues.count - 1);
                                       range.y = math::clamp(range.y, 0, (int32)dyn.molecule.residues.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           const auto beg = dyn.molecule.residues[i].atom_idx.beg;
                                           const auto end = dyn.molecule.residues[i].atom_idx.end;
                                           memset(mask.data + beg, 1, end - beg);
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"resname", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   memset(mask.data, 0, mask.count);
                                   for (int i = 0; i < args.count; i++) {
                                       for (const auto& res : dyn.molecule.residues) {
                                           if (compare(args[i], res.name)) {
                                               const auto beg = res.atom_idx.beg;
                                               const auto end = res.atom_idx.end;
                                               memset(mask.data + beg, 1, end - beg);
                                           }
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"resid", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (auto range : ranges) {
                                       range.x = math::clamp(range.x - 1, 0, (int32)dyn.molecule.residues.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)dyn.molecule.residues.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           const auto beg = dyn.molecule.residues[i].atom_idx.beg;
                                           const auto end = dyn.molecule.residues[i].atom_idx.end;
                                           memset(mask.data + beg, 1, end - beg);
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"chain", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   DynamicArray<IntRange> ranges;
                                   if (!extract_ranges(&ranges, args)) return false;
                                   memset(mask.data, 0, mask.size_in_bytes());
                                   for (auto range : ranges) {
                                       range.x = math::clamp(range.x - 1, 0, (int32)dyn.molecule.residues.count - 1);
                                       range.y = math::clamp(range.y - 1, 0, (int32)dyn.molecule.residues.count - 1);
                                       for (int i = range.x; i <= range.y; i++) {
                                           Chain chain = get_chain(dyn.molecule, (ChainIdx)i);
                                           const auto beg = get_atom_beg_idx(dyn.molecule, chain);
                                           const auto end = get_atom_end_idx(dyn.molecule, chain);
                                           memset(mask.data + beg, 1, end - beg);
                                       }
                                   }
                                   return true;
                               }});

    filter_commands.push_back({"chainid", [](Array<bool> mask, const MoleculeDynamic& dyn, Array<const CString> args) {
                                   memset(mask.data, 0, mask.count);
                                   for (int i = 0; i < args.count; i++) {
                                       for (const auto& chain : dyn.molecule.chains) {
                                           if (compare(args[i], chain.id)) {
                                               const auto beg = get_atom_beg_idx(dyn.molecule, chain);
                                               const auto end = get_atom_end_idx(dyn.molecule, chain);
                                               memset(mask.data + beg, 1, end - beg);
                                           }
                                       }
                                   }
                                   return true;
                               }});
}

void shutdown() {}

}  // namespace filter
