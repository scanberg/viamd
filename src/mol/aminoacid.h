#pragma once

#include <core/string_types.h>

enum class AminoAcid : unsigned char {
	Unknown = 0,
	Ala,
	Arg,
	Asn,
	Asp,
	Cys,
	Gln,
	Glu,
	Gly,
	His,
	Ile,
	Leu,
	Lys,
	Met,
	Phe,
	Pro,
	Ser,
	Thr,
	Trp,
	Tyr,
	Val,
	Sec,
	Pyl,
	Asx,
	Glx,
	Xle
};

namespace aminoacid {

// TODO: Add hydrophobicity and other animo-acid parameters
// Amino Acid functions
constexpr const char* name(AminoAcid amino);
constexpr const char* symbol(AminoAcid amino);
AminoAcid get_from_string(CString cstr);
}  // namespace animoacid