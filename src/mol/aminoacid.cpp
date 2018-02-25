#include "mol/aminoacid.h"
#include <array>
#include <string.h>
#include <ctype.h>

namespace aminoacid {

static constexpr std::array<const char*, 26> names = {{"Unknown",
                                                       "Alanine",
                                                       "Arginine",
                                                       "Asparagine",
                                                       "Aspartic acid",
                                                       "Cysteine",
                                                       "Glutamic acid",
                                                       "Glutamine",
                                                       "Glycine",
                                                       "Histidine",
                                                       "Isoleucine",
                                                       "Luecine",
                                                       "Lysine",
                                                       "Methionine",
                                                       "Phenylalanine",
                                                       "Proline",
                                                       "Serine",
                                                       "Threonine",
                                                       "Tryptophan",
                                                       "Tyrosine",
                                                       "Valine",
                                                       "Selencysteine",
                                                       "Pyrrolysine",
                                                       "Asparagine or Aspartic acid",
                                                       "Glutamine or glutamic acid",
                                                       "Leucine or Isloleucine"}};

static constexpr std::array<const char*, 26> symbols = {{"XAA", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
                                                         "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASC", "GLX", "XLE"}};

constexpr const char* name(AminoAcid amino) { return names[static_cast<int>(amino)]; }
constexpr const char* symbol(AminoAcid amino) { return symbols[static_cast<int>(amino)]; }

AminoAcid get_from_string(const char* cstr) {

    if (!cstr[0] || !cstr[1] || !cstr[2]) return AminoAcid::Unknown;
    const char seq[3] = {toupper(cstr[0]), toupper(cstr[1]), toupper(cstr[2])};

    for (size_t i = 0; i < symbols.size(); i++) {
        if (symbols[i][0] == seq[0] &&
			symbols[i][1] == seq[1] &&
			symbols[i][2] == seq[2])
			return static_cast<AminoAcid>(i);
    }
    return AminoAcid::Unknown;
}
}  // namespace aminoacid