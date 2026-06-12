# Molden Format: Edge Cases and Implementation Notes

## Overview

This document describes edge cases, format variations, and implementation considerations for Molden file support in VIAMD. The Molden format lacks a formal specification and is primarily documented through examples from various quantum chemistry programs.

## 1. Format Variations by Program

### 1.1 GAMESS
- Uses Cartesian basis functions by default
- Includes `[5D]` or `[5D10F]` tags for spherical harmonics
- Symmetry labels are typically single characters (A, B, E, T)
- May include `[GEOCONV]` section with optimization history

### 1.2 Gaussian
- Can output in Cartesian or Spherical format
- Uses descriptive symmetry labels (Ag, Bu, etc.)
- May include normalization factors different from standard
- Sometimes includes `[FR-COORD]` for fractional coordinates

### 1.3 ORCA
- Typically uses spherical harmonics for d and f orbitals
- Clear section delimiters
- Consistent formatting
- Often includes `[MO]` energies in eV (needs conversion)

### 1.4 Molpro
- May use non-standard section ordering
- Custom section names possible
- Different coefficient ordering for some basis sets

## 2. Section-Specific Edge Cases

### 2.1 [Atoms] Section

#### Standard Format
```
[Atoms] (AU|Angs)
element_symbol  atom_index  atomic_number  x  y  z
```

#### Edge Cases

**A. Missing Unit Specification**
```
[Atoms]
C   1   6   0.000000   0.000000   0.000000
```
- **Solution**: Default to Angstroms (most common)
- **Detection**: Check if `(AU|Angs)` is present on header line

**B. Omitted Atomic Numbers**
```
[Atoms] Angs
C   1      0.000000   0.000000   0.000000
H   2      1.089000   0.000000   0.000000
```
- **Solution**: Infer from element symbol using lookup table
- **Implementation**: `util::element_symbol_to_atomic_number()`

**C. Non-Sequential Atom Indices**
```
[Atoms] Angs
C   1   6   0.0   0.0   0.0
H   3   1   1.0   0.0   0.0  // Index 2 is missing
```
- **Solution**: Use index as-is, don't assume consecutive
- **Validation**: Check for duplicate indices

**D. Varying Whitespace**
```
[Atoms] Angs
C    1    6    0.000    0.000    0.000
  H  2    1  1.089    0.000    0.000
```
- **Solution**: Split on any whitespace, trim leading/trailing
- **Implementation**: Use flexible tokenization

### 2.2 [GTO] Section

#### Standard Format
```
[GTO]
atom_index  0
shell_type  num_primitives  [scale_factor]
exponent  coefficient
...
blank_line
```

#### Edge Cases

**A. SP Shells (Combined S and P)**
```
1  0
sp  3  1.0
13.01  0.03975  0.02388
1.962  0.38075  0.15680
0.4446 0.66750  0.48750
```
- **Challenge**: Two coefficient columns
- **Solution**: Store both coefficients in `PrimitiveGaussian::coefficient_sp`
- **Expansion**: May need to split into separate S and P shells for processing
- **Basis count**: Counts as 1 (S) + 3 (P) = 4 basis functions

**B. Missing Scale Factor**
```
2  0
s  3
13.01  0.03975
1.962  0.38075
```
- **Solution**: Default to 1.0
- **Common**: Most files omit when scale factor is 1.0

**C. The Mysterious "0" After Atom Index**
```
1  0  // What does the 0 mean?
```
- **Observation**: Always present, meaning unclear
- **Solution**: Parse and ignore, likely historical artifact
- **Speculation**: May have indicated basis set type in older versions

**D. Blank Line Delimiters**
```
1  0
s  3  1.0
...

2  0  // Double blank lines sometimes appear
s  3  1.0
```
- **Solution**: Skip all blank lines between atom basis sets
- **Parser state**: Reset to "expecting atom header" on blank line

### 2.3 [MO] Section

#### Standard Format
```
[MO]
Ene= -11.234
Spin= Alpha
Occup= 2.0
Sym= A1g
0.1234
-0.5678
...
```

#### Edge Cases

**A. Missing Symmetry Labels**
```
[MO]
Ene= -11.234
Spin= Alpha
Occup= 2.0
0.1234  // No Sym= line
```
- **Solution**: Leave `symmetry` field empty
- **Common**: Many programs omit symmetry for restricted calculations
- **Validation**: Don't require Sym= field

**B. Coefficient Count Mismatch**
```
[MO]
Ene= -11.234
Spin= Alpha
Occup= 2.0
0.1234
0.5678
// Missing coefficients!
[MO]  // Next orbital starts
```
- **Detection**: Count coefficients, compare to total basis functions
- **Solution**: Mark orbital as invalid, log warning
- **Recovery**: Continue parsing remaining orbitals

**C. Multiple Coefficients Per Line**
Some programs write multiple coefficients on one line:
```
0.1234 -0.5678 0.9012
0.3456 -0.7890
```
- **Solution**: Split each line and collect all coefficients
- **Format**: Variable number of coefficients per line

**D. Different Energy Units**
```
Ene= -11.234  // Hartrees (a.u.)
Ene= -305.6   // eV (needs conversion)
```
- **Detection**: Large negative values likely eV
- **Solution**: Add optional energy unit field
- **Conversion**: 1 Hartree = 27.2114 eV
- **Phase 1**: Document issue, assume Hartrees

**E. Virtual vs Occupied Orbitals**
```
Occup= 2.0   // Doubly occupied
Occup= 1.0   // Singly occupied (radicals)
Occup= 0.0   // Virtual (unoccupied)
Occup= 1.5   // Fractional (some methods)
```
- **Validation**: Check range [0.0, 2.0]
- **Common**: 0.0, 1.0, 2.0 for closed shell
- **Unusual**: Fractional occupations from ensemble methods

### 2.4 [5D], [5D7F], [5D10F] Tags

These tags indicate spherical harmonic basis functions:

```
[5D]      // 5 d functions instead of 6 Cartesian
[5D7F]    // 5 d and 7 f functions
[5D10F]   // 5 d and 10 f (Cartesian f)
```

- **Impact**: Changes basis function count
- **Implementation**: Set `BasisFormat` accordingly
- **Default**: Cartesian if not specified

## 3. Basis Function Ordering

### 3.1 Cartesian Ordering

| Shell | Functions | Ordering |
|-------|-----------|----------|
| s     | 1         | s |
| p     | 3         | px, py, pz |
| d     | 6         | d(x²), d(y²), d(z²), d(xy), d(xz), d(yz) |
| f     | 10        | f(x³), f(y³), f(z³), f(x²y), f(x²z), f(y²x), f(y²z), f(z²x), f(z²y), f(xyz) |

**Note**: Exact ordering may vary by program!

### 3.2 Spherical Harmonic Ordering

| Shell | Functions | Ordering (m quantum number) |
|-------|-----------|------------------------------|
| s     | 1         | m=0 |
| p     | 3         | m=-1,0,+1 |
| d     | 5         | m=-2,-1,0,+1,+2 |
| f     | 7         | m=-3,-2,-1,0,+1,+2,+3 |

## 4. Multiple Geometry Support

Some files contain multiple geometries:

```
[GEOMETRIES] XYZ
3
Comment for geometry 1
C 0.0 0.0 0.0
H 1.0 0.0 0.0
H 0.0 1.0 0.0
3
Comment for geometry 2
C 0.1 0.0 0.0
...
```

- **Phase 1**: Not implemented, use first geometry only
- **Future**: Support as trajectory-like data
- **Alternative**: Some use multiple `[Atoms]` sections

## 5. Character Encoding and Line Endings

- **Line endings**: Unix (LF), Windows (CRLF), or old Mac (CR)
- **Solution**: Handle all line ending types
- **Encoding**: Usually ASCII, sometimes UTF-8
- **Special characters**: Avoid in critical fields (coefficients, coordinates)

## 6. Numerical Precision Issues

### 6.1 Exponent Ranges
```
// Very large exponents (core orbitals)
1234567.89  0.00001

// Very small exponents (diffuse functions)  
0.000012  0.56789
```
- **Solution**: Use double precision
- **Validation**: Check for extremely small or negative exponents

### 6.2 Scientific Notation
```
1.234E+03  5.678E-04
1.234D+03  5.678D-04  // Fortran-style (D instead of E)
```
- **Solution**: Support both E and D notation
- **Parsing**: Replace 'D' with 'E' before parsing

### 6.3 Missing Decimal Points
Some Fortran programs output:
```
1234567    0.56789
0.000012   123456
```
- **Challenge**: Hard to determine if integer or decimal
- **Context-dependent**: Exponents are typically larger
- **Solution**: Type-aware parsing based on field position

## 7. Validation Checklist

### Critical Validations
- [ ] Atom indices are positive integers
- [ ] Atomic numbers match element symbols (if both present)
- [ ] Coordinate units are specified or assumed
- [ ] Shell types are valid (s, p, d, f, g, sp)
- [ ] Exponents are positive
- [ ] Coefficient count matches basis function count for each MO
- [ ] Occupation numbers are in range [0.0, 2.0]
- [ ] No duplicate atom indices

### Warning-Level Issues
- [ ] Missing symmetry labels (optional)
- [ ] Non-sequential atom indices
- [ ] Scale factors != 1.0
- [ ] Fractional occupation numbers
- [ ] Unusual basis set sizes

## 8. Test Data Requirements

For robust implementation, test with:

1. **Minimal example**: Small molecule (H₂O) with minimal basis (STO-3G)
2. **SP shells**: Molecule with SP shell basis sets
3. **Cartesian basis**: 6d, 10f functions
4. **Spherical basis**: 5d, 7f functions with [5D7F] tag
5. **Open shell**: Radical with alpha/beta orbitals
6. **Large system**: Protein or complex with 100+ atoms
7. **Various programs**: GAMESS, Gaussian, ORCA output
8. **Edge cases**: Missing fields, non-standard formatting

## 9. Implementation Priority

### Phase 1 (Current) - Data Structures ✓
- Define all data structures
- Document edge cases
- Create validation utilities

### Phase 2 - Basic Parser
- Parse [Atoms] section
- Parse [GTO] section (Cartesian only)
- Parse [MO] section
- Handle most common edge cases

### Phase 3 - Robust Parser
- Support SP shells
- Support spherical harmonics
- Handle all documented edge cases
- Comprehensive error reporting

### Phase 4 - Advanced Features
- Multiple geometries
- Frequency data
- Integration with VIAMD visualization

## 10. References

### Documentation Sources
- GAMESS manual Appendix on Molden format
- Examples from quantum chemistry output files
- Community knowledge (mailing lists, forums)
- Source code of programs that read Molden files

### Related Formats
- Gaussian fchk (formatted checkpoint)
- Gaussian cube (for volumetric data)
- WFN/WFX (wavefunction files)
- These may be better suited for certain use cases

## 11. Error Handling Strategy

### Parse Errors
- **Syntax errors**: Invalid section headers, malformed lines
- **Response**: Skip invalid lines, log warning, continue parsing

### Validation Errors  
- **Data inconsistency**: Coefficient count mismatch
- **Response**: Mark data as incomplete, allow partial loading

### Critical Errors
- **Missing required sections**: No [Atoms] section
- **Response**: Abort parsing, return error

### Error Reporting
- Provide line numbers for parse errors
- Descriptive error messages
- Validation summary (warnings + errors)

## 12. Performance Considerations

### Large Files
- Some files can be 100+ MB (large basis sets, many orbitals)
- **Solution**: Stream parsing, don't load entire file into memory
- **Optimization**: Reserve vector capacity when sizes are known

### Memory Usage
- MO coefficients can be large (N_basis × N_orbitals doubles)
- **Estimate**: 1000 basis functions × 1000 orbitals × 8 bytes = 8 MB
- **Solution**: Acceptable for modern systems, but consider lazy loading

### Parsing Speed
- Text parsing is typically not the bottleneck
- **Focus**: Correctness over speed in Phase 1
- **Optimization**: Profile before optimizing in later phases
