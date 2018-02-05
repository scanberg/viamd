/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 * Version 0.9
 *
 * Copyright (c) 2012-2016 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/

#include <mol/aminoacid.h>
#include <array>
#include <string.h>
#include <ctype.h>

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

static constexpr std::array<const char*, 26> symbols = {
    {"XAA", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
     "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASC", "GLX", "XLE"}};

constexpr const char* name(AminoAcid amino) { return names[static_cast<int>(amino)]; }
constexpr const char* symbol(AminoAcid amino) { return symbols[static_cast<int>(amino)]; }

AminoAcid getFromString(const char* cstr, int length) {
    if (length == -1)
        length = strlen(cstr);

    if (length >= 3) {
        char seq[3] = {(char)toupper(cstr[0]), (char)toupper(cstr[1]), (char)toupper(cstr[2])};
        for (size_t i = 0; i < symbols.size(); i++) {
            if (strncmp(seq,symbols[i], 3) == 0)
                return static_cast<AminoAcid>(i);
        }
    }
    return AminoAcid::Unknown;
}