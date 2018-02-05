#include <mol/element.h>
#include <array>

static constexpr std::array<const char*, 119> names = {
    {"Unknown",     "Hydrogen",     "Helium",       "Lithium",     "Beryllium",   "Boron",         "Carbon",     "Nitrogen",   "Oxygen",
     "Fluorine",    "Neon",         "Sodium",       "Magnesium",   "Aluminium",   "Silicon",       "Phosphorus", "Sulfur",     "Chlorine",
     "Argon",       "Potassium",    "Calcium",      "Scandium",    "Titanium",    "Vanadium",      "Chromium",   "Manganese",  "Iron",
     "Cobalt",      "Nickel",       "Copper",       "Zinc",        "Gallium",     "Germanium",     "Arsenic",    "Selenium",   "Bromine",
     "Krypton",     "Rubidium",     "Strontium",    "Yttrium",     "Zirconium",   "Niobium",       "Molybdenum", "Technetium", "Ruthenium",
     "Rhodium",     "Palladium",    "Silver",       "Cadmium",     "Indium",      "Tin",           "Antimony",   "Tellurium",  "Iodine",
     "Xenon",       "Caesium",      "Barium",       "Lanthanum",   "Cerium",      "Praseodymium",  "Neodymium",  "Promethium", "Samarium",
     "Europium",    "Gadolinium",   "Terbium",      "Dysprosium",  "Holmium",     "Erbium",        "Thulium",    "Ytterbium",  "Lutetium",
     "Hafnium",     "Tantalum",     "Tungsten",     "Rhenium",     "Osmium",      "Iridium",       "Platinum",   "Gold",       "Mercury",
     "Thallium",    "Lead",         "Bismuth",      "Polonium",    "Astatine",    "Radon",         "Francium",   "Radium",     "Actinium",
     "Thorium",     "Protactinium", "Uranium",      "Neptunium",   "Plutonium",   "Americium",     "Curium",     "Berkelium",  "Californium",
     "Einsteinium", "Fermium",      "Mendelevium",  "Nobelium",    "Lawrencium",  "Rutherfordium", "Dubnium",    "Seaborgium", "Bohrium",
     "Hassium",     "Meitnerium",   "Darmstadtium", "Roentgenium", "Copernicium", "Nihonium",      "Flerovium",  "Moscovium",  "Livermorium",
     "Tennessine",  "Oganesson"}};

static constexpr std::array<const char*, 119> symbols = {
    {"Xx", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
     "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
     "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
     "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
     "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"}};

// http://dx.doi.org/10.1039/b801115j
static constexpr std::array<float, 119> covalent_radii = {
    {2.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53,
     1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45,
     1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,
     1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80,
     1.69, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60}};

// https://dx.doi.org/10.1021/jp8111556
static constexpr std::array<float, 119> vdw_radii = {
    {0.00, 1.10, 1.40, 1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, 2.30, 2.15, 2.05,
     2.05, 2.05, 2.05, 2.00, 2.00, 2.00, 2.10, 1.87, 2.11, 1.85, 1.90, 1.83, 2.02, 3.03, 2.49, 2.40, 2.30, 2.15, 2.10, 2.05, 2.05, 2.00, 2.05, 2.10,
     2.20, 2.20, 1.93, 2.17, 2.06, 1.98, 2.16, 3.43, 2.68, 2.50, 2.48, 2.47, 2.45, 2.43, 2.42, 2.40, 2.38, 2.37, 2.35, 2.33, 2.32, 2.30, 2.28, 2.27,
     2.25, 2.20, 2.10, 2.05, 2.00, 2.00, 2.05, 2.10, 2.05, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, 2.00, 2.40, 2.00, 2.30, 2.00, 2.00, 2.00,
     2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00}};

// http://jmol.sourceforge.net/jscolors/
static constexpr std::array<unsigned int, 119> colors = {
    {0x000000FF, 0xFFFFFFFF, 0xD9FFFFFF, 0xB22222FF, 0xC2FF00FF, 0xFFB5B5FF, 0xB0B0B0FF, 0x8F8FFFFF, 0xF00000FF, 0x90E050FF, 0xB3E3F5FF, 0xAB5CF2FF,
     0x8AFF00FF, 0x808090FF, 0xF0C8A0FF, 0xFFA500FF, 0xFFC832FF, 0x1FF01FFF, 0x80D1E3FF, 0x8F40D4FF, 0x808090FF, 0xE6E6E6FF, 0x808090FF, 0xA6A6ABFF,
     0x808090FF, 0x808090FF, 0xFFA500FF, 0xF090A0FF, 0xA52A2AFF, 0xA52A2AFF, 0xA52A2AFF, 0xC28F8FFF, 0x668F8FFF, 0xBD80E3FF, 0xFFA100FF, 0xA52A2AFF,
     0x5CB8D1FF, 0x702EB0FF, 0x00FF00FF, 0x94FFFFFF, 0x94E0E0FF, 0x73C2C9FF, 0x54B5B5FF, 0x3B9E9EFF, 0x248F8FFF, 0x0A7D8CFF, 0x006985FF, 0x808090FF,
     0xFFD98FFF, 0xA67573FF, 0x668080FF, 0x9E63B5FF, 0xD47A00FF, 0x940094FF, 0x429EB0FF, 0x57178FFF, 0xFFA500FF, 0x70D4FFFF, 0xFFFFC7FF, 0xD9FFC7FF,
     0xC7FFC7FF, 0xA3FFC7FF, 0x8FFFC7FF, 0x61FFC7FF, 0x45FFC7FF, 0x30FFC7FF, 0x1FFFC7FF, 0x00FF9CFF, 0x00E675FF, 0x00D452FF, 0x00BF38FF, 0x00AB24FF,
     0x4DC2FFFF, 0x4DA6FFFF, 0x2194D6FF, 0x267DABFF, 0x266696FF, 0x175487FF, 0xD0D0E0FF, 0xFFD123FF, 0xB8B8D0FF, 0xA6544DFF, 0x575961FF, 0x9E4FB5FF,
     0xAB5C00FF, 0x754F45FF, 0x428296FF, 0x420066FF, 0x007D00FF, 0x70ABFAFF, 0x00BAFFFF, 0x00A1FFFF, 0x008FFFFF, 0x0080FFFF, 0x006BFFFF, 0x545CF2FF,
     0x785CE3FF, 0x8A4FE3FF, 0xA136D4FF, 0xB31FD4FF, 0xB31FBAFF, 0xB30DA6FF, 0xBD0D87FF, 0xC70066FF, 0xCC0059FF, 0xD1004FFF, 0xD90045FF, 0xE00038FF,
     0xE6002EFF, 0xEB0026FF, 0xF00022FF, 0xF60020FF, 0xF8001EFF, 0xFA001CFF, 0xFC001AFF, 0xFD0018FF, 0xFE0016FF, 0xFF0014FF, 0xFF0012FF}};

constexpr const char* name(Element symbol) { return names[static_cast<int>(symbol)]; }
constexpr const char* symbol(Element symbol) { return symbols[static_cast<int>(symbol)]; }
constexpr unsigned int color(Element symbol) { return colors[static_cast<int>(symbol)]; }
constexpr float vdw_radius(Element symbol) { return vdw_radii[static_cast<int>(symbol)]; }
constexpr float covalent_radius(Element symbol) { return covalent_radii[static_cast<int>(symbol)]; }

Element get_from_string(const char* cstr, int length = -1) {
    if (length == -1)
        length = strlen(cstr);

    if (length == 0) return Element::Unknown;
    //length = length < 3 ? length : 3;

    const char* beg = cstr;
    const char* end = cstr + length;
    while (!isupper(*beg) && beg != end) beg++;
    if (beg == end) return Element::Unknown;
    const char* tmp = beg + 1;
    while (tmp != end && islower(*tmp)) tmp++;
    end = tmp;
    if (end - beg > 3) return Element::Unknown;

    for (size_t i = 0; i < symbols.size(); i++) {
        if (strncmp(beg, symbols[i], strlen(symbols[i])) == 0) {
            return static_cast<Element>(i);
        }
    }
    return Element::Unknown;
}
