#include "mol/aminoacid.h"
#include <string.h>
#include <ctype.h>

namespace aminoacid {

static constexpr unsigned int NUM_AMINO_ACIDS = 26;

static constexpr const char* names[NUM_AMINO_ACIDS] = {"Unknown",
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
                                                       "Leucine or Isloleucine"};

static constexpr const char* symbols[NUM_AMINO_ACIDS] = {"XAA", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
                                                         "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASC", "GLX", "XLE"};

constexpr const char* name(AminoAcid amino) { return names[(int)amino]; }
constexpr const char* symbol(AminoAcid amino) { return symbols[(int)amino]; }

AminoAcid get_from_string(CString cstr) {

	// Skip leading numbers and crap
	while (cstr.count > 0 && !isalpha(*cstr.beg())) {
		cstr.data++;
		cstr.count--;
	};

	if (cstr.count != 3) return AminoAcid::Unknown;
    const char seq[3] = {(char)toupper(cstr[0]), (char)toupper(cstr[1]), (char)toupper(cstr[2])};

    for (unsigned int i = 0; i < NUM_AMINO_ACIDS; i++) {
        if (symbols[i][0] == seq[0] && symbols[i][1] == seq[1] && symbols[i][2] == seq[2]) return (AminoAcid)i;
    }
    return AminoAcid::Unknown;
}
}  // namespace aminoacid
