#pragma once

enum class Element : unsigned char {
    Unknown = 0,
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe,
    Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn,
    Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W,
    Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf,
    Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Uut, Fl, Mc, Lv, Ts, Og
};

// Element functions
constexpr const char*   name(Element symbol);
constexpr const char*   symbol(Element symbol);
constexpr unsigned int  color(Element symbol);
constexpr float         vdw_radius(Element symbol);
constexpr float         covalent_radius(Element symbol);

Element get_from_string(const char* cstr, int length);
