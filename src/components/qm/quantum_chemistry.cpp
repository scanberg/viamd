#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd_event.h>
#include <viamd.h>
#include <representation.h>
#include <task_system.h>
#include <color_utils.h>

#include <gfx/volumerender_utils.h>
#include <gfx/isomesh_utils.h>
#include <iso_meshes.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <md_gto.h>
#include <md_gto_iso.h>
#include <md_util.h>
#include <md_topo.h>
#include <md_vector_graphics.h>
#include <md_flow.h>
#include <md_qm_ct.h>
#include <md_fragment.h>
#include <md_csv.h>
#include <md_xvg.h>
#include <core/md_vec_math.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

#include <serialization_utils.h>

#include <imgui_widgets.h>
#include <imgui_internal.h>

#include "flow_draw.h"
#include "fragment_tree.h"

#include <implot_widgets.h>
#include <implot_internal.h>

#include <algorithm>
#include <atomic>
#include <app/IconsFontAwesome6.h>

#define BLK_DIM 8
#define ROTATORY_STRENGTH_IN_CGS 471.443648175
#define EXTINCTION_COEFFICIENT_FROM_BETA 19.603697575813566
#define ROTATORY_STRENGTH_TO_DELTA_EPSILON 0.01386075702557538652

#define HARTREE_TO_KJ_PER_MOL 2625.4996394799
#define HARTREE_TO_EV 27.2114079527
#define EV_TO_HARTREE (1.0 / HARTREE_TO_EV)

#define AU_TO_GM 40479.02797814119
#define OSCILLATOR_STRENGTH_TO_EPSILON 1054.9516366171872

#define DEFAULT_SAMPLES_PER_ANGSTROM 8
#define DEFAULT_GTO_CUTOFF_VALUE 1.0e-6
#define MAX_NTO_GROUPS 16
#define MAX_NTO_LAMBDAS 3
// Upper bound on the NTO pairs the flow diagram keeps LABEL storage for. The lambda axis of
// vlx/rsp/nto/lambda is padded to the widest row in the file, so this caps what is drawn rather
// than what exists - a file with more pairs than this simply stops at this many columns.
#define MAX_NTO_LABELS 32
#define NTO_LAMBDA_CUTOFF_VALUE 0.1

#define FORCE_CPU_PATH 0

#define USE_AABB_GRID 0

// Resolution for broadened plots
#define NUM_SAMPLES 1024

#define IM_GREEN ImVec4(0, 1, 0, 1)
#define IM_RED ImVec4(1, 0, 0, 1)
#define IM_YELLOW ImVec4(1, 1, 0.5, 0.3)
#define IM_BLUE ImVec4(0.5, 0.5, 1, 0.3)

#define COLOR_PEAK_SELECTED ImVec4(0.5, 0.5, 1.0, 1.0)
#define COLOR_PEAK_HOVER    ImVec4(1.0, 1.0, 0.5, 1.0)
// COLOR_PEAK_HOVER as a table row background. Kept translucent: the row text and the selectable's
// own fill both draw on top of it, and an opaque tint makes either unreadable.
#define COLOR_ROW_HIGHLIGHT ImVec4(1.0, 1.0, 0.5, 0.25)

#define COLOR_TABLE_SELECTED IM_BLUE
#define COLOR_TABLE_HOVERED  IM_YELLOW

#define U32_MAGENTA IM_COL32(255, 0, 255, 255)

// The attachment/detachment pair. The teal is VeloxChem's own accent colour, kept because the
// transition analysis has been read against it in print; the magenta is its complement, so the two
// densities separate on either side of neutral rather than by lightness alone.
#define U32_ATTACHMENT  IM_COL32(0, 162, 135, 191)
#define VEC4_ATTACHMENT {0, 162.0f/255.0f, 135.0f/255.0f, 0.75f}
#define VEC4_DETACHMENT {162.0f/255.0f, 35.0f/255.0f, 135.0f/255.0f, 0.75f}

constexpr uint64_t interaction_surface_nto = HASH_STR_LIT64("interaction surface nto");
constexpr uint64_t interaction_surface_orb = HASH_STR_LIT64("interaction surface orb");
constexpr uint64_t interaction_surface_flow = HASH_STR_LIT64("interaction surface flow");

constexpr PickingDomainID PickingDomain_CriticalPoints = HASH_STR_LIT64("Picking Domain Critical Points");

// How a critical point of the electron density is shown, by type (md_topo_critical_point_type_t). Named as in
// Bader's QTAIM (and AIMAll): attractor, bond, ring and cage critical point, with the rank and signature of the
// Hessian there, (3,-3) to (3,+3); beside it the topological name md_topo uses (maximum, split saddle, join
// saddle, minimum). An attractor is not called nuclear: non-nuclear attractors exist. The colours are Okabe and
// Ito's, told apart with each common colour vision deficiency, and none of them is the selection's blue or the
// highlight's yellow: those mark a point with a ring around it, so a marked point keeps the colour of its type.
// Size tells the type too: the attractors, at the nuclei, are the largest.
struct CriticalPointStyle {
    const char* name;       // QTAIM, as a title
    const char* label;      // QTAIM, short, for the table and the legend
    const char* plural;
    const char* signature;  // (rank, signature)
    const char* topology;   // md_topo's name
    const char* topology_plural;
    vec4_t      color;      // display referred: the overlay goes over the image after the display transform
    float       radius;     // Å, at a point size of 1
};

static const CriticalPointStyle critical_point_styles[5] = {
    { "Undefined",            "Undefined", "Undefined",  "",       "undefined",    "undefined",     {0.60f, 0.60f, 0.60f, 1.0f}, 0.09f },
    { "Attractor",            "Attractor", "Attractors", "(3,-3)", "maximum",      "maxima",        {0.80f, 0.47f, 0.65f, 1.0f}, 0.12f },   // MD_TOPO_MAXIMUM: reddish purple
    { "Bond critical point",  "BCP",       "BCPs",       "(3,-1)", "split saddle", "split saddles", {0.90f, 0.62f, 0.00f, 1.0f}, 0.09f },   // MD_TOPO_SPLIT_SADDLE: orange
    { "Cage critical point",  "CCP",       "CCPs",       "(3,+3)", "minimum",      "minima",        {0.00f, 0.45f, 0.70f, 1.0f}, 0.09f },   // MD_TOPO_MINIMUM: blue
    { "Ring critical point",  "RCP",       "RCPs",       "(3,+1)", "join saddle",  "join saddles",  {0.00f, 0.62f, 0.45f, 1.0f}, 0.09f },   // MD_TOPO_JOIN_SADDLE: bluish green
};

static inline const CriticalPointStyle& critical_point_style(int type) {
    return critical_point_styles[(type >= 0 && type < (int)MD_COUNTOF(critical_point_styles)) ? type : 0];
}

// The order the types are listed and sorted in: attractors, bond, ring and cage critical points
static constexpr int critical_point_type_order[4] = { MD_TOPO_MAXIMUM, MD_TOPO_SPLIT_SADDLE, MD_TOPO_JOIN_SADDLE, MD_TOPO_MINIMUM };

static inline int critical_point_type_rank(int type) {
    for (int k = 0; k < 4; ++k) if (critical_point_type_order[k] == type) return k;
    return 4;
}

// The type's marker as a swatch of the text's height, for a legend or a table cell
static void critical_point_swatch(int type) {
    const vec4_t c = critical_point_style(type).color;
    const float h = ImGui::GetTextLineHeight();
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const ImVec2 center = ImVec2(p.x + 0.5f * h, p.y + 0.5f * h);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->AddCircleFilled(center, 0.40f * h, IM_COL32(26, 26, 26, 230));
    dl->AddCircleFilled(center, 0.31f * h, ImGui::ColorConvertFloat4ToU32(ImVec4(c.x, c.y, c.z, 1.0f)));
    ImGui::Dummy(ImVec2(h, h));
}

// --- Export of the critical points (QuantumChemistry::critical_points_export)
//
// Every number is written at the precision it was computed in, so that it reads back to the same bits: the
// analytic search's points come in double (md_topo_extremum_graph_t::vertices_f64: positions, density and the
// Hessian there), written with 17 significant digits; the volume pipeline's in float, with 9. Lengths in bohr,
// the units of the computation, the density in au (e bohr^-3), its second derivatives in e bohr^-5.
struct CriticalPointsExport {
    const md_topo_extremum_graph_t* graph = nullptr;
    bool   analytic  = false;   // else the volume pipeline
    double rho_min   = 0.0;     // analytic: the threshold of the search
    bool   complete  = false;   // analytic: every critical point above the threshold is there
    bool   cancelled = false;
    md_str_t  field     = {};      // the density's attribute path
};

// The Morse index (negative eigenvalues of the Hessian): TTK's CellDimension of a critical point
static int critical_point_morse_index(int type) {
    switch (type) {
    case MD_TOPO_MAXIMUM:      return 3;
    case MD_TOPO_SPLIT_SADDLE: return 2;
    case MD_TOPO_JOIN_SADDLE:  return 1;
    case MD_TOPO_MINIMUM:      return 0;
    default:                   return -1;
    }
}

// A point's position and density in double, from vertices_f64 where the graph has it
static void critical_point_xyzv(double out[4], const md_topo_extremum_graph_t& g, size_t i) {
    if (g.vertices_f64) {
        out[0] = g.vertices_f64[i].x; out[1] = g.vertices_f64[i].y; out[2] = g.vertices_f64[i].z; out[3] = g.vertices_f64[i].value;
    } else {
        out[0] = g.vertices[i].x; out[1] = g.vertices[i].y; out[2] = g.vertices[i].z; out[3] = g.vertices[i].value;
    }
}

// A number as JSON and VTK read it back exactly: 17 significant digits of a double, 9 of a float. Not finite
// (which neither format has a literal for in common): JSON null, VTK 0.
static void write_number(md_file_t file, double x, bool f64, bool json) {
    if (!isfinite(x)) {
        md_file_printf(file, json ? "null" : "0");
    } else {
        md_file_printf(file, f64 ? "%.17g" : "%.9g", x);
    }
}

static void critical_points_write_json(md_file_t f, const CriticalPointsExport& ex) {
    const md_topo_extremum_graph_t& g = *ex.graph;
    const bool f64 = g.vertices_f64 != nullptr;
    md_file_printf(f, "{\n");
    md_file_printf(f, "  \"format\": \"viamd-critical-points\",\n");
    md_file_printf(f, "  \"format_version\": 1,\n");
    md_file_printf(f, "  \"description\": \"Critical points of the electron density and the paths between them, as VIAMD found them (QTAIM terms; md_topo's in topology)\",\n");
    md_file_printf(f, "  \"units\": {\"length\": \"bohr\", \"density\": \"e bohr^-3\", \"hessian\": \"e bohr^-5\"},\n");
    md_file_printf(f, "  \"precision\": \"%s\",\n", f64 ? "float64" : "float32");
    md_file_printf(f, "  \"field\": \"" MD_STR_FMT "\",\n", MD_STR_ARG(ex.field));
    md_file_printf(f, "  \"method\": \"%s\",\n", ex.analytic ? "analytic" : "volume");
    if (ex.analytic) {
        md_file_printf(f, "  \"rho_min\": ");
        write_number(f, ex.rho_min, true, true);
        md_file_printf(f, ",\n  \"complete\": %s,\n", ex.complete ? "true" : "false");
        md_file_printf(f, "  \"cancelled\": %s,\n", ex.cancelled ? "true" : "false");
    }
    md_file_printf(f, "  \"critical_points\": [");
    for (size_t i = 0; i < g.num_vertices; ++i) {
        const int type = g.types[i];
        const CriticalPointStyle& s = critical_point_style(type);
        const int index = critical_point_morse_index(type);
        double p[4];
        critical_point_xyzv(p, g, i);
        md_file_printf(f, "%s\n    {\"index\": %zu, \"type\": \"%s\", \"name\": \"%s\", \"rank\": 3, \"signature\": %d, \"topology\": \"%s\", \"morse_index\": %d, \"position\": [",
                       i ? "," : "", i, s.label, s.name, index >= 0 ? 3 - 2 * index : 0, s.topology, index);
        for (int k = 0; k < 3; ++k) {
            if (k) md_file_printf(f, ", ");
            write_number(f, p[k], f64, true);
        }
        md_file_printf(f, "], \"rho\": ");
        write_number(f, p[3], f64, true);
        if (f64) {
            const md_topo_vert_f64_t& d = g.vertices_f64[i];
            const double* ev = d.hessian_eigenvalues;
            md_file_printf(f, ", \"laplacian\": ");
            write_number(f, ev[0] + ev[1] + ev[2], true, true);
            md_file_printf(f, ", \"hessian_eigenvalues\": [");
            for (int k = 0; k < 3; ++k) {
                if (k) md_file_printf(f, ", ");
                write_number(f, ev[k], true, true);
            }
            md_file_printf(f, "], \"hessian_eigenvectors\": [");
            for (int k = 0; k < 3; ++k) {
                md_file_printf(f, "%s[", k ? ", " : "");
                for (int j = 0; j < 3; ++j) {
                    if (j) md_file_printf(f, ", ");
                    write_number(f, d.hessian_eigenvectors[k][j], true, true);
                }
                md_file_printf(f, "]");
            }
            md_file_printf(f, "]");
            if (type == MD_TOPO_SPLIT_SADDLE) {
                // lambda1 / lambda2 - 1 of the two negative eigenvalues, |lambda1| >= |lambda2|
                md_file_printf(f, ", \"ellipticity\": ");
                write_number(f, ev[0] / ev[1] - 1.0, true, true);
            }
        }
        md_file_printf(f, "}");
    }
    md_file_printf(f, "\n  ],\n");
    md_file_printf(f, "  \"paths\": [");
    for (size_t i = 0; i < g.num_edges; ++i) {
        const md_topo_edge_t e = g.edges[i];
        const bool bond = e.from < g.num_vertices && g.types[e.from] == MD_TOPO_SPLIT_SADDLE;
        md_file_printf(f, "%s\n    {\"index\": %zu, \"type\": \"%s\", \"saddle\": %u, \"extremum\": %u}",
                       i ? "," : "", i, bond ? "bond path" : "RCP attractor path", e.from, e.to);
    }
    md_file_printf(f, "\n  ],\n");
    md_file_printf(f, "  \"notes\": [\n");
    md_file_printf(f, "    \"index: the point's index in this file, which paths refer to. position: bohr.\",\n");
    md_file_printf(f, "    \"rank, signature: QTAIM's (rank, signature) of the Hessian; morse_index: its number of negative eigenvalues.\",\n");
    md_file_printf(f, "    \"hessian_eigenvalues: ascending; hessian_eigenvectors[k]: the unit eigenvector of eigenvalue k. ellipticity: bond critical points only.\",\n");
    md_file_printf(f, "    \"paths: which points a gradient path joins (a bond path from a BCP to an attractor, an RCP attractor path between an RCP and a CCP), not the path's geometry.\"\n");
    md_file_printf(f, "  ]\n");
    md_file_printf(f, "}\n");
}

// VTK XML PolyData: a vertex cell per critical point, a line cell per path, from the saddle to the extremum.
// The arrays are named as the Topology ToolKit names those of a Morse-Smale complex where they mean the same
// (CellDimension, CellId; SeparatrixType, SeparatrixId, SourceId, DestinationId), and after QTAIM where TTK
// has none. A cell array is per cell, vertices first: -1 for a vertex where it is a path's.
static void critical_points_write_vtp(md_file_t f, const CriticalPointsExport& ex) {
    const md_topo_extremum_graph_t& g = *ex.graph;
    const bool f64 = g.vertices_f64 != nullptr;
    const size_t n = g.num_vertices, m = g.num_edges;
    const char* real = f64 ? "Float64" : "Float32";

    md_file_printf(f, "<?xml version=\"1.0\"?>\n");
    md_file_printf(f, "<!-- Critical points of the electron density (" MD_STR_FMT ") and the paths between them, as VIAMD found them\n", MD_STR_ARG(ex.field));
    md_file_printf(f, "     (%s). Lengths in bohr, the density in au (e bohr^-3), the Hessian in e bohr^-5; %s.\n",
                   ex.analytic ? "the analytic search" : "the volume pipeline", f64 ? "every number in double, to 17 digits" : "in float, to 9 digits");
    if (ex.analytic) {
        md_file_printf(f, "     The search: rho_min %.17g, %s.\n", ex.rho_min, ex.cancelled ? "cancelled" : (ex.complete ? "complete" : "incomplete"));
    }
    md_file_printf(f, "     CellDimension: the Morse index (3 attractor, 2 BCP, 1 RCP, 0 CCP); Signature: QTAIM's (3, Signature).\n");
    md_file_printf(f, "     SeparatrixType: 2 a bond path, 0 an RCP attractor path; the lines join their points, they are not the paths' geometry. -->\n");
    md_file_printf(f, "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    md_file_printf(f, "  <PolyData>\n");
    md_file_printf(f, "    <Piece NumberOfPoints=\"%zu\" NumberOfVerts=\"%zu\" NumberOfLines=\"%zu\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n", n, n, m);

    // --- per point
    md_file_printf(f, "      <PointData Scalars=\"Density\">\n");
    md_file_printf(f, "        <DataArray type=\"%s\" Name=\"Density\" format=\"ascii\">\n", real);
    for (size_t i = 0; i < n; ++i) {
        double p[4];
        critical_point_xyzv(p, g, i);
        md_file_printf(f, "          ");
        write_number(f, p[3], f64, false);
        md_file_printf(f, "\n");
    }
    md_file_printf(f, "        </DataArray>\n");
    md_file_printf(f, "        <DataArray type=\"Int32\" Name=\"CellDimension\" format=\"ascii\">\n");
    for (size_t i = 0; i < n; ++i) md_file_printf(f, "          %d\n", critical_point_morse_index(g.types[i]));
    md_file_printf(f, "        </DataArray>\n");
    md_file_printf(f, "        <DataArray type=\"Int32\" Name=\"Signature\" format=\"ascii\">\n");
    for (size_t i = 0; i < n; ++i) {
        const int index = critical_point_morse_index(g.types[i]);
        md_file_printf(f, "          %d\n", index >= 0 ? 3 - 2 * index : 0);
    }
    md_file_printf(f, "        </DataArray>\n");
    md_file_printf(f, "        <DataArray type=\"Int64\" Name=\"CellId\" format=\"ascii\">\n");
    for (size_t i = 0; i < n; ++i) md_file_printf(f, "          %zu\n", i);
    md_file_printf(f, "        </DataArray>\n");
    if (f64) {
        md_file_printf(f, "        <DataArray type=\"Float64\" Name=\"Laplacian\" format=\"ascii\">\n");
        for (size_t i = 0; i < n; ++i) {
            const double* ev = g.vertices_f64[i].hessian_eigenvalues;
            md_file_printf(f, "          ");
            write_number(f, ev[0] + ev[1] + ev[2], true, false);
            md_file_printf(f, "\n");
        }
        md_file_printf(f, "        </DataArray>\n");
        md_file_printf(f, "        <DataArray type=\"Float64\" Name=\"HessianEigenvalues\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (size_t i = 0; i < n; ++i) {
            md_file_printf(f, "          ");
            for (int k = 0; k < 3; ++k) {
                if (k) md_file_printf(f, " ");
                write_number(f, g.vertices_f64[i].hessian_eigenvalues[k], true, false);
            }
            md_file_printf(f, "\n");
        }
        md_file_printf(f, "        </DataArray>\n");
        for (int k = 0; k < 3; ++k) {
            md_file_printf(f, "        <DataArray type=\"Float64\" Name=\"HessianEigenvector%d\" NumberOfComponents=\"3\" format=\"ascii\">\n", k);
            for (size_t i = 0; i < n; ++i) {
                md_file_printf(f, "          ");
                for (int j = 0; j < 3; ++j) {
                    if (j) md_file_printf(f, " ");
                    write_number(f, g.vertices_f64[i].hessian_eigenvectors[k][j], true, false);
                }
                md_file_printf(f, "\n");
            }
            md_file_printf(f, "        </DataArray>\n");
        }
    }
    md_file_printf(f, "      </PointData>\n");

    // --- per cell: the vertices, then the lines
    struct { const char* name; } cell_arrays[] = { {"SeparatrixType"}, {"SeparatrixId"}, {"SourceId"}, {"DestinationId"} };
    md_file_printf(f, "      <CellData>\n");
    for (int a = 0; a < 4; ++a) {
        md_file_printf(f, "        <DataArray type=\"Int64\" Name=\"%s\" format=\"ascii\">\n", cell_arrays[a].name);
        for (size_t i = 0; i < n; ++i) md_file_printf(f, "          -1\n");
        for (size_t i = 0; i < m; ++i) {
            const md_topo_edge_t e = g.edges[i];
            const bool bond = e.from < n && g.types[e.from] == MD_TOPO_SPLIT_SADDLE;
            const int64_t v = a == 0 ? (bond ? 2 : 0) : a == 1 ? (int64_t)i : a == 2 ? (int64_t)e.from : (int64_t)e.to;
            md_file_printf(f, "          %lld\n", (long long)v);
        }
        md_file_printf(f, "        </DataArray>\n");
    }
    md_file_printf(f, "      </CellData>\n");

    md_file_printf(f, "      <Points>\n");
    md_file_printf(f, "        <DataArray type=\"%s\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n", real);
    for (size_t i = 0; i < n; ++i) {
        double p[4];
        critical_point_xyzv(p, g, i);
        md_file_printf(f, "          ");
        for (int k = 0; k < 3; ++k) {
            if (k) md_file_printf(f, " ");
            write_number(f, p[k], f64, false);
        }
        md_file_printf(f, "\n");
    }
    md_file_printf(f, "        </DataArray>\n");
    md_file_printf(f, "      </Points>\n");

    md_file_printf(f, "      <Verts>\n");
    md_file_printf(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for (size_t i = 0; i < n; ++i) md_file_printf(f, "          %zu\n", i);
    md_file_printf(f, "        </DataArray>\n");
    md_file_printf(f, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for (size_t i = 0; i < n; ++i) md_file_printf(f, "          %zu\n", i + 1);
    md_file_printf(f, "        </DataArray>\n");
    md_file_printf(f, "      </Verts>\n");
    md_file_printf(f, "      <Lines>\n");
    md_file_printf(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for (size_t i = 0; i < m; ++i) md_file_printf(f, "          %u %u\n", g.edges[i].from, g.edges[i].to);
    md_file_printf(f, "        </DataArray>\n");
    md_file_printf(f, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for (size_t i = 0; i < m; ++i) md_file_printf(f, "          %zu\n", 2 * (i + 1));
    md_file_printf(f, "        </DataArray>\n");
    md_file_printf(f, "      </Lines>\n");
    md_file_printf(f, "    </Piece>\n");
    md_file_printf(f, "  </PolyData>\n");
    md_file_printf(f, "</VTKFile>\n");
}

// This is the internal storage order of volumes and textures related to NTOs
enum NTO {
    NTO_Attachment,
    NTO_Detachment,
};

// Which spin channel an orbital belongs to. The component's own two valued vocabulary, not a
// reader's: alpha is 0 and beta is 1, so the value indexes the two element label arrays and the two
// homo_idx/lumo_idx slots directly, and es_path turns it into the attribute to read.
//
// ElectronicStructureSpin answers a different question - WHICH DENSITY a representation draws, with
// Total and Difference beside the two channels - and its Alpha and Beta are not 0 and 1. Reusing it
// here would leave every "which of the two orbital sets" site doing arithmetic to get an index back.
enum SpinChannel {
    SPIN_ALPHA = 0,
    SPIN_BETA  = 1,
};

static const char* electronic_structure_value_mode_str(ElectronicStructureSource source, ElectronicStructureSpin spin, bool use_magnitude) {
    switch (source) {
    case ElectronicStructureSource::MolecularOrbital:
    case ElectronicStructureSource::NaturalTransitionOrbital:
        return use_magnitude ? "Magnitude" : "Signed";
    case ElectronicStructureSource::TransitionDensity:
        return "Density";
    case ElectronicStructureSource::ElectronDensity:
        return spin == ElectronicStructureSpin::Difference ? (use_magnitude ? "Magnitude" : "Signed") : "Density";
    default:
        return "";
    }
}

// What the Charge Transfer diagram shows.
//
// The two routes the prototype had besides these are gone. 'Groups only' split the group to group
// flow proportionally from its two marginals, which is a guess. 'Atoms via NTOs' routed atoms
// through the NTO pairs with one product per pair; that drops the interference between pairs, and
// when two local excitations mix into delocalized NTOs it shows charge transfer that is not there
// (a CIS state of test_data/vlx/mol.h5 with three comparable pairs: 19% amine <-> phenyl transfer
// shown, 0.5% real). The fragment route computes the transfer exactly instead (md_qm_ct.h).
enum flow_mode_t {
    FLOW_MODE_FRAGMENTS = 0,    // hole fragment -> particle fragment: the exact transfer matrix
    FLOW_MODE_ORBITALS,         // occupied MO -> NTO pair -> virtual MO
    FLOW_MODE_COUNT,
};

static const char* flow_mode_str[] = { "Fragments", "Orbitals" };

enum x_unit_t {
    X_UNIT_EV,
    X_UNIT_NM,
    X_UNIT_CM_INVERSE,
    X_UNIT_HARTREE,
    X_UNIT_COUNT,
};

static const char* x_unit_full_str[] = {"Energy (eV)", "Wavelength (nm)", (const char*)u8"Wavenumber (cm⁻¹)", "Energy (au)"};
static const char* x_unit_label_str[] = {"Energy", "Wavelength", "Wavenumber", "Energy"};
static const char* x_unit_short_str[] = {"eV", "nm", (const char*)u8"cm⁻¹", "au"};

enum broadening_mode_t {
    BROADENING_MODE_GAUSSIAN,
    BROADENING_MODE_LORENTZIAN,
    BROADENING_MODE_COUNT,
};

static const char* broadening_mode_str[] = { "Gaussian", "Lorentzian" };

// Attribute charge contributions of an AO-basis density matrix D to atoms or groups
// using a Mulliken-style partitioning: q_g = sum_{mu in g} sum_nu D[mu,nu] * S[mu,nu].
// Both D and S are symmetric (num_ao x num_ao) row-major matrices.
// The total tr(D*S) is exactly preserved: sum_g out_charge[g] == tr(D*S).
static void attribute_charge_density(double out_charge[], const int* ao_to_idx, const double* D, const double* S, size_t num_ao) {
	MD_ASSERT(D);
	MD_ASSERT(S);
	MD_ASSERT(ao_to_idx);

	for (size_t mu = 0; mu < num_ao; ++mu) {
		const int idx = ao_to_idx[mu];
		const double* D_row = &D[mu * num_ao];
		const double* S_row = &S[mu * num_ao];
		double s = 0.0;
		for (size_t nu = 0; nu < num_ao; ++nu) {
			s += D_row[nu] * S_row[nu];
		}
		out_charge[idx] += s;
	}
}

// Mulliken partitioning of a SINGLE orbital vector, which is the same quantity as above for the
// rank-one density D = C C^T. Written out rather than routed through attribute_charge_density so
// it does not have to materialize an N x N matrix per orbital - with 32 NTOs and a few hundred
// AOs, forming those densities is most of the work and none of the answer.
//
//   q_g = sum_{mu in g} C[mu] * (S C)[mu],  and  sum_g q_g = C^T S C = 1 for a normalized orbital.
static void attribute_orbital_density(double out_charge[], const int* ao_to_idx, const double* C, const double* S, size_t num_ao) {
	MD_ASSERT(C);
	MD_ASSERT(S);
	MD_ASSERT(ao_to_idx);

	for (size_t mu = 0; mu < num_ao; ++mu) {
		const double* S_row = &S[mu * num_ao];
		double sc = 0.0;
		for (size_t nu = 0; nu < num_ao; ++nu) {
			sc += S_row[nu] * C[nu];
		}
		out_charge[ao_to_idx[mu]] += C[mu] * sc;
	}
}

// Column identifiers for the Orbital Grid's "Molecular Orbitals" table. File scope because the
// active sort column is remembered in Orb across frames, and because the order it produces is
// what the orbital grid takes its span from.
enum OrbitalColumn {
    OrbitalColumn_Idx,
    OrbitalColumn_OccAlpha,
    OrbitalColumn_OccBeta,
    OrbitalColumn_EneAlpha,
    OrbitalColumn_EneBeta,
};

// The panels' isosurfaces as meshes (IsosurfaceRendering_Mesh): a field per panel, all extracted together
// (iso_meshes.h), and the running mean of each panel's image
struct PanelMeshes {
    static constexpr int MAX = (int)ISO_MESHES_MAX_SLOTS;
    md_gto_iso_field_t* field[MAX] = {};
    uint64_t  field_key[MAX] = {};
    int       slot[MAX] = {};       // the panel's slot in 'meshes', -1 for none
    IsoMeshes meshes = {};
    isomesh::Target target[MAX] = {};
    bool      converged = true;     // every panel's image is final: the panels need not be redrawn
};

static void panel_meshes_free(PanelMeshes* pm, ApplicationState* state) {
    iso_meshes_free(&pm->meshes, state);
    for (int i = 0; i < PanelMeshes::MAX; ++i) {
        iso_field_release(pm->field[i]);
        isomesh::target_free(&pm->target[i]);
    }
    *pm = {};
}

// Where the atoms are: a field is of one geometry
static uint64_t panel_positions_hash(const ApplicationState& state) {
    const md_system_state_t& st = state.mold.state;
    return (st.xyz && st.num_atoms) ? md_hash64(st.xyz, sizeof(st.xyz[0]) * st.num_atoms, 0x9a7e1) : 0x9a7e1;
}

// Panel i's field, made by 'make' when its key changed (0: none)
template <typename F>
static void panel_meshes_field(PanelMeshes* pm, ApplicationState* state, int i, uint64_t key, F make) {
    if (pm->field_key[i] == key) return;
    // The extraction state holds the fields: let it go before one of them does
    iso_meshes_free(&pm->meshes, state);
    iso_field_release(pm->field[i]);
    pm->field[i] = key ? make() : nullptr;
    pm->field_key[i] = key;
}

// The meshes of the first 'count' panels at one isovalue
static void panel_meshes_update(PanelMeshes* pm, ApplicationState* state, int count, double iso, uint32_t sheets) {
    IsoMeshSlot slots[PanelMeshes::MAX];
    uint32_t n = 0;
    uint64_t key = 0x5107;
    for (int i = 0; i < PanelMeshes::MAX; ++i) {
        pm->slot[i] = -1;
        if (i < count && pm->field[i] && iso > 0.0) {
            pm->slot[i] = (int)n;
            slots[n++] = {pm->field[i], iso, sheets};
            key = md_hash64(&pm->field_key[i], sizeof(pm->field_key[i]), key);
        }
    }
    iso_meshes_update(&pm->meshes, state, slots, n, iso_meshes_cell_size(VolumeResolution::Mid), key);
}

// Panel i's surfaces: the positive sheet in colors[0], the negative one in colors[1]
static size_t panel_meshes_surfaces(isomesh::Surface out[2], const PanelMeshes& pm, int i, const vec4_t colors[2]) {
    const int k = pm.slot[i];
    if (k < 0 || !pm.meshes.ok) return 0;
    size_t n = 0;
    for (int s = 0; s < 2; ++s) {
        if (pm.meshes.index_count[k][s] == 0) continue;
        isomesh::Surface& surf = out[n++];
        surf = {};
        surf.mesh = &pm.meshes.mesh;
        surf.index_offset = pm.meshes.index_offset[k][s];
        surf.index_count  = pm.meshes.index_count[k][s];
        surf.model = iso_meshes_model_matrix();
        surf.color = colors[s];
    }
    return n;
}

// Fills row_to_mo with the mo indices in the order the table lists them, and mo_to_row with its
// inverse. The grid shows a span of table ROWS rather than a span of mo indices, so both views
// agree by construction whatever the sort is - see draw_orb_window.
//
// The default, mo index descending, is the conventional orbital diagram order with the highest
// energy on top, and is exactly what this table showed before it became sortable.
//
// sum_occupancy: in the restricted layouts a single "Occupancy" column shows alpha + beta, and
// the sort key has to be the value on screen, not the half of it that alpha contributes.
static void orbital_table_build_order(int* row_to_mo, int* mo_to_row, int num_mos, int sort_column, bool ascending,
                                      const double* occ_alpha, const double* occ_beta,
                                      const double* ene_alpha, const double* ene_beta, bool sum_occupancy)
{
    for (int i = 0; i < num_mos; ++i) {
        row_to_mo[i] = i;
    }

    // A column whose data never arrived falls back to the index, which keeps the order total and
    // the comparator a strict weak ordering rather than sorting on a null dereference.
    auto key_of = [&](int n) -> double {
        switch (sort_column) {
        case OrbitalColumn_OccAlpha: return occ_alpha ? ((sum_occupancy && occ_beta) ? occ_alpha[n] + occ_beta[n] : occ_alpha[n]) : (double)n;
        case OrbitalColumn_OccBeta:  return occ_beta  ? occ_beta[n]  : (double)n;
        case OrbitalColumn_EneAlpha: return ene_alpha ? ene_alpha[n] : (double)n;
        case OrbitalColumn_EneBeta:  return ene_beta  ? ene_beta[n]  : (double)n;
        case OrbitalColumn_Idx:
        default:                     return (double)n;
        }
    };

    std::sort(row_to_mo, row_to_mo + num_mos, [&](int a, int b) {
        const double ka = key_of(a);
        const double kb = key_of(b);
        if (ka != kb) {
            return ascending ? (ka < kb) : (ka > kb);
        }
        // Ties are the common case, not the exception - occupancies are mostly 2, 1 or 0, and
        // symmetry makes energies degenerate - so the mo index breaks them, in the same direction
        // as the sort so that a degenerate run reads the way the index column would.
        return ascending ? (a < b) : (a > b);
    });

    for (int r = 0; r < num_mos; ++r) {
        mo_to_row[row_to_mo[r]] = r;
    }
}

// What the charge transfer flow graph procedures know about a calculation, and the whole of it.
// Gathered by the caller - the only place that has to know whether this came from an attribute
// table, a reader, or somewhere else entirely. Past that point there is no system and no reader
// here: arrays, counts and indices, and that is the entire contract.
//
// Row major throughout. S is {A,A} and C_mo is {M,A}, so orbital m starts at C_mo + m * num_ao. The
// NTO blocks are {K,A} for ONE excited state, with nto_lambda[k] the weight of pair k. That axis is
// PADDED - a state has as many pairs as it has - and the pad carries weight zero, which the lambda
// cutoff every consumer already applies discards in exactly the same place. ao_to_atom is one atom
// index per AO.
//
// Any pointer may be null and any count zero: a file without response data has no NTOs, and a
// builder checks what it uses rather than assuming a complete set.
struct FlowData {
    size_t        num_ao      = 0;
    size_t        num_mo      = 0;
    int           homo_idx    = -1;      // -1 when the channel holds no electrons

    const double* S           = nullptr; // {A,A} AO overlap
    const double* C_mo        = nullptr; // {M,A} MO coefficients
    const int*    ao_to_atom  = nullptr; // [A]

    size_t        num_nto     = 0;       // K, padded
    const double* nto_lambda  = nullptr; // [K]
    const double* C_nto_hole  = nullptr; // {K,A}
    const double* C_nto_part  = nullptr; // {K,A}

    // The orbitals in the Lowdin basis, {M,A} like C_mo; null until the job computing them is done
    const double* C_lowdin    = nullptr;
};

// The density whose critical points are reported. QTAIM is defined on the total density; a
// calculation that publishes a single spin channel publishes only the alpha one, and that one is then
// the total.
static md_str_t critical_points_density_path(const md_system_t& sys) {
    return md_attributes_find(&sys.attributes, es_path::total_density) ? es_path::total_density : es_path::alpha_density;
}

// One run of the analytic critical point search (md_topo_compute_extremum_graph_gto) on the task pool.
// Everything the run reads - basis, geometry, density matrix - is copied in here before it starts, so
// it never touches the attribute table, the component's arena or the system while the main thread
// carries on, and a system freed or reloaded underneath it takes nothing away from it. The main thread
// reads graph and info only once 'done' is set, and frees the job only after waiting for the task.
struct GtoCriticalPointJob {
    md_allocator_i*    arena    = nullptr;   // owns the input copies below
    md_gto_basis_t     basis    = {};
    float*             atom_xyz = nullptr;   // [3 * num_atoms] Bohr
    double*            density  = nullptr;   // [num_ao * num_ao]
    md_topo_gto_desc_t desc     = {};

    md_topo_extremum_graph_t graph   = {};   // heap, filled by the run
    md_topo_gto_info_t       info    = {};
    bool                     ok      = false;
    double                   seconds = 0.0;

    uint64_t          hash    = 0;           // the inputs this run was started for
    md_tick_t         start   = 0;
    bool              dropped = false;       // main thread only: superseded, never to be shown
    volatile int32_t  cancel  = 0;
    volatile float    progress = 0.0f;       // written by the run: 0 to 1, an estimate (md_topo_gto_desc_t)
    std::atomic<bool> done    = false;
    task_system::Job  task;
#if MD_ENABLE_GPU
    md_gpu_stream_t   stream  = nullptr;   // own stream: the sweep runs on the GPU, used by the task only
#endif
};

// Starts a search of the system as it is now. Returns null, with the reason logged, when the system
// does not carry what the search needs.
static GtoCriticalPointJob* gto_critical_point_job_start(ApplicationState& app, double rho_min, bool use_gpu, uint64_t hash) {
    const md_system_t& sys = app.mold.sys;
    const md_system_state_t& state = app.mold.state;
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        MD_LOG_ERROR("Critical points: the system carries no GTO basis");
        return nullptr;
    }
    const size_t num_ao = md_gto_basis_num_ao(&basis);

    size_t num_atoms = 0;
    const vec3_t* atom_pos = basis_atom_positions_extract(&num_atoms, temp, sys, state);
    if (!atom_pos || num_atoms != md_gto_basis_num_atoms(&basis)) {
        MD_LOG_ERROR("Critical points: no positions for the atoms the basis is centred on");
        return nullptr;
    }

    const md_str_t path = critical_points_density_path(sys);
    size_t dim = 0;
    const double* density = density_matrix_extract(&dim, temp, sys, path, nullptr);
    if (!density || dim != num_ao) {
        MD_LOG_ERROR("Critical points: no density matrix matching the basis (%zu AOs) at '" MD_STR_FMT "'", num_ao, MD_STR_ARG(path));
        return nullptr;
    }

    GtoCriticalPointJob* job = new GtoCriticalPointJob();
    job->arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_allocator_i* alloc = job->arena;

    job->basis.num_shells     = basis.num_shells;
    job->basis.num_primitives = basis.num_primitives;
    job->basis.shells = (md_gto_shell_t*)md_alloc(alloc, sizeof(md_gto_shell_t) * basis.num_shells);
    job->basis.alpha  = (float*)md_alloc(alloc, sizeof(float) * basis.num_primitives);
    job->basis.coeff  = (float*)md_alloc(alloc, sizeof(float) * basis.num_primitives);
    md_memcpy(job->basis.shells, basis.shells, sizeof(md_gto_shell_t) * basis.num_shells);
    md_memcpy(job->basis.alpha,  basis.alpha,  sizeof(float) * basis.num_primitives);
    md_memcpy(job->basis.coeff,  basis.coeff,  sizeof(float) * basis.num_primitives);

    job->atom_xyz = (float*)md_alloc(alloc, sizeof(float) * 3 * num_atoms);
    for (size_t i = 0; i < num_atoms; ++i) {
        job->atom_xyz[3 * i + 0] = atom_pos[i].x;
        job->atom_xyz[3 * i + 1] = atom_pos[i].y;
        job->atom_xyz[3 * i + 2] = atom_pos[i].z;
    }

    job->density = (double*)md_alloc(alloc, sizeof(double) * num_ao * num_ao);
    md_memcpy(job->density, density, sizeof(double) * num_ao * num_ao);

    job->desc.basis              = &job->basis;
    job->desc.atom_xyz           = job->atom_xyz;
    job->desc.atom_xyz_stride    = sizeof(float) * 3;
    job->desc.density_matrix     = job->density;
    job->desc.rho_min            = rho_min;
    job->desc.trace_separatrices = true;
    job->desc.cancel             = &job->cancel;
    job->desc.executor           = task_system::executor();   // its loops on the pool, stopped with the job
    job->desc.progress           = &job->progress;

    job->graph.alloc = md_get_heap_allocator();
    job->hash  = hash;
    job->start = md_tick_now();
#if MD_ENABLE_GPU
    // The octree sweep goes to the GPU when there is one; the CPU keeps the double precision parts.
    if (use_gpu && app.gpu_device) job->stream = md_gpu_stream_create(app.gpu_device, MD_GPU_STREAM_COMPUTE, "Critical points");
#endif

    // Hidden as a task: the task overlay's bar for it would show no fraction, and its close button
    // interrupts a task through a flag this run never sees. The overlay shows the run's own progress
    // instead (EventType_ViamdTaskProgressRequest), as does the Critical Points window, each with a cancel button.
    job->task.start({}, [job]() {
#if MD_ENABLE_GPU
        if (job->stream) {
            job->ok = md_topo_compute_extremum_graph_gto_gpu(&job->graph, &job->info, &job->desc, job->stream);
        } else
#endif
        job->ok = md_topo_compute_extremum_graph_gto(&job->graph, &job->info, &job->desc);
        job->seconds = md_tick_to_seconds(md_tick_now() - job->start);
        job->done.store(true, std::memory_order_release);
    });
    return job;
}

static void gto_critical_point_job_destroy(GtoCriticalPointJob* job) {
    if (!job) return;
    job->cancel = 1;
    job->task.wait();
#if MD_ENABLE_GPU
    if (job->stream) md_gpu_stream_destroy(job->stream);
#endif
    md_topo_extremum_graph_free(&job->graph);
    if (job->arena) md_arena_allocator_destroy(job->arena);
    delete job;
}

// The orbitals in the Lowdin basis (md_qm_lowdin_orbitals), which the charge transfer numbers are
// computed in. One O(M^3) computation per system, so a job rather than a stall: for a system near
// the top of what VIAMD is meant for (~10k AOs) it is a matter of minutes on one core, spread over
// the task system's threads (products and eigendecomposition alike). Every excited state then
// reuses the result, and the overlap and the AO to atom map it
// keeps with it, so a state's own job (CtStateJob) needs nothing from the attribute table.
struct LowdinJob {
    double*   C        = nullptr;   // {M,A} input, owned
    double*   C_lowdin = nullptr;   // {M,A} output, owned
    double*   S        = nullptr;   // {A,A} owned
    uint32_t* ao_atom  = nullptr;   // [A] QM atom of each AO, num_atoms for none
    size_t    num_mo    = 0;
    size_t    num_ao    = 0;
    size_t    num_atoms = 0;        // QM atoms
    bool      ok        = false;
    std::atomic<bool> done = false;
    task_system::Job task;
};

static void lowdin_job_destroy(LowdinJob* job) {
    if (!job) return;
    job->task.cancel_and_wait();
    md_allocator_i* heap = md_get_heap_allocator();
    if (job->C)        md_free(heap, job->C,        sizeof(double) * job->num_mo * job->num_ao);
    if (job->C_lowdin) md_free(heap, job->C_lowdin, sizeof(double) * job->num_mo * job->num_ao);
    if (job->S)        md_free(heap, job->S,        sizeof(double) * job->num_ao * job->num_ao);
    if (job->ao_atom)  md_free(heap, job->ao_atom,  sizeof(uint32_t) * job->num_ao);
    delete job;
}

// Null when the system carries no orbitals, or no overlap or basis matching them
static LowdinJob* lowdin_job_start(const md_system_t& sys, const md_gto_basis_t& basis, size_t num_qm_atoms) {
    size_t num_mo = 0, num_ao = 0;
    if (!es_orbital_extent(sys, &num_mo, &num_ao) || num_mo == 0 || num_ao == 0) return nullptr;
    const md_attribute_t* coeff = md_attributes_find(&sys.attributes, es_path::alpha_coefficient);
    const md_attribute_t* ovl   = md_attributes_find(&sys.attributes, es_path::overlap);
    if (!coeff || !ovl || ovl->format.rank != 2 || ovl->format.shape[0] != num_ao || ovl->format.shape[1] != num_ao) {
        MD_LOG_ERROR("Charge transfer: the system has no overlap matrix matching its orbitals");
        return nullptr;
    }
    if (basis.num_shells == 0 || md_gto_basis_num_ao(&basis) != num_ao) {
        MD_LOG_ERROR("Charge transfer: the basis does not match the orbitals");
        return nullptr;
    }

    md_allocator_i* heap = md_get_heap_allocator();
    LowdinJob* job = new LowdinJob();
    job->num_mo    = num_mo;
    job->num_ao    = num_ao;
    job->num_atoms = num_qm_atoms;
    job->C        = (double*)  md_alloc(heap, sizeof(double) * num_mo * num_ao);
    job->C_lowdin = (double*)  md_alloc(heap, sizeof(double) * num_mo * num_ao);
    job->S        = (double*)  md_alloc(heap, sizeof(double) * num_ao * num_ao);
    job->ao_atom  = (uint32_t*)md_alloc(heap, sizeof(uint32_t) * num_ao);
    const bool read =
        job->C && job->C_lowdin && job->S && job->ao_atom &&
        md_attribute_extract_f64(job->C, num_mo * num_ao, coeff, md_attribute_slice_all(), md_unit_none()) == num_mo * num_ao &&
        md_attribute_extract_f64(job->S, num_ao * num_ao, ovl,   md_attribute_slice_all(), md_unit_none()) == num_ao * num_ao &&
        md_gto_basis_ao_to_atom(job->ao_atom, &basis) == num_ao;
    if (!read) {
        MD_LOG_ERROR("Charge transfer: could not read the orbitals and the overlap");
        job->done.store(true, std::memory_order_release);
        return job;
    }
    for (size_t i = 0; i < num_ao; ++i) {
        if (job->ao_atom[i] >= num_qm_atoms) job->ao_atom[i] = (uint32_t)num_qm_atoms;
    }
    job->task.start(MD_STR_LIT("Lowdin basis"), [job]() {
        job->ok = md_qm_lowdin_orbitals(job->C_lowdin, job->C, job->num_mo, job->num_ao, task_system::executor());
        job->done.store(true, std::memory_order_release);
    });
    return job;
}

static bool lowdin_job_ready(const LowdinJob* job) {
    return job && job->done.load(std::memory_order_acquire) && job->ok;
}

// One excited state's charge transfer numbers per QM ATOM, from its NTO pairs (md_qm_ct.h): the pairs in
// the Lowdin basis, their K x K blocks per atom (which fragments sum, and the per ribbon orbital pairs are
// computed from), and the atom to atom transfer matrix. Everything coarser is a sum over these, done on
// the main thread when the diagram is built.
struct CtStateJob {
    int32_t state = -1;
    const LowdinJob* lowdin = nullptr;  // outlives this job: destroyed after it
    size_t  K = 0;
    size_t  num_atoms = 0;
    double* lambda = nullptr;           // [K]
    double* hole   = nullptr;           // {K,A}, AO basis in; Lowdin basis after the run
    double* part   = nullptr;
    double* hole_blocks = nullptr;      // {N,K,K}
    double* part_blocks = nullptr;
    double* omega  = nullptr;           // {N,N}, normalized to sum 1, the tail below 1e-6 of it dropped
    double* weight[2] = {};             // [N] per atom: hole (row sums), particle (column sums)
    float   nto_dropped  = 0.0f;        // share of the excitation in the pairs left out
    float   link_dropped = 0.0f;        // share in the atom pairs below the floor
    bool    ok = false;
    std::atomic<bool> done = false;
    task_system::Job task;
};

static void ct_state_job_destroy(CtStateJob* job) {
    if (!job) return;
    job->task.cancel_and_wait();
    md_allocator_i* heap = md_get_heap_allocator();
    const size_t A = job->lowdin ? job->lowdin->num_ao : 0;
    const size_t N = job->num_atoms, K = job->K;
    if (job->lambda)      md_free(heap, job->lambda,      sizeof(double) * K);
    if (job->hole)        md_free(heap, job->hole,        sizeof(double) * K * A);
    if (job->part)        md_free(heap, job->part,        sizeof(double) * K * A);
    if (job->hole_blocks) md_free(heap, job->hole_blocks, sizeof(double) * N * K * K);
    if (job->part_blocks) md_free(heap, job->part_blocks, sizeof(double) * N * K * K);
    if (job->omega)       md_free(heap, job->omega,       sizeof(double) * N * N);
    for (int c = 0; c < 2; ++c) if (job->weight[c]) md_free(heap, job->weight[c], sizeof(double) * N);
    delete job;
}

static bool ct_state_job_ready(const CtStateJob* job) {
    return job && job->done.load(std::memory_order_acquire) && job->ok;
}

// A state's NTO pairs, read on the main thread from the attribute table into heap arrays the caller owns:
// every pair above round off of the whole (with the cross terms in, a small pair can still move a large
// share between two fragments). nto_dropped is the share of the excitation left out.
struct NtoPairs {
    size_t  K = 0;
    double* lambda = nullptr;   // [K]
    double* hole   = nullptr;   // {K,A}
    double* part   = nullptr;   // {K,A}
    float   nto_dropped = 0.0f;
};

static bool nto_pairs_read(NtoPairs* out, const md_system_t& sys, size_t A, int32_t state) {
    *out = {};
    if (state < 0 || A == 0) return false;
    const md_attribute_t* lambda_attr = md_attributes_find(&sys.attributes, es_path::nto_lambda);
    const md_attribute_t* hole_attr   = md_attributes_find(&sys.attributes, es_path::nto_hole);
    const md_attribute_t* part_attr   = md_attributes_find(&sys.attributes, es_path::nto_particle);
    if (!lambda_attr || !hole_attr || !part_attr) return false;

    const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)state);
    const size_t L = md_attribute_slice_count(lambda_attr, slice);
    if (L == 0 || md_attribute_slice_count(hole_attr, slice) < L * A || md_attribute_slice_count(part_attr, slice) < L * A) return false;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    double* lam  = md_temp_alloc_array(temp, double, L);
    double* hole = md_temp_alloc_array(temp, double, L * A);
    double* part = md_temp_alloc_array(temp, double, L * A);
    if (md_attribute_extract_f64(lam,  L,     lambda_attr, slice, md_unit_none()) < L ||
        md_attribute_extract_f64(hole, L * A, hole_attr,   slice, md_unit_none()) < L * A ||
        md_attribute_extract_f64(part, L * A, part_attr,   slice, md_unit_none()) < L * A) {
        return false;
    }

    double lambda_all = 0.0;
    for (size_t k = 0; k < L; ++k) lambda_all += MD_MAX(0.0, lam[k]);
    if (lambda_all <= 0.0) return false;
    size_t K = 0;
    for (size_t k = 0; k < L; ++k) K += (lam[k] >= 1.0e-6 * lambda_all) ? 1 : 0;
    if (K == 0) return false;

    md_allocator_i* heap = md_get_heap_allocator();
    out->K = K;
    out->lambda = (double*)md_alloc(heap, sizeof(double) * K);
    out->hole   = (double*)md_alloc(heap, sizeof(double) * K * A);
    out->part   = (double*)md_alloc(heap, sizeof(double) * K * A);
    double lambda_kept = 0.0;
    for (size_t k = 0, j = 0; k < L; ++k) {
        if (!(lam[k] >= 1.0e-6 * lambda_all)) continue;
        out->lambda[j] = lam[k];
        lambda_kept += lam[k];
        md_memcpy(out->hole + j * A, hole + k * A, sizeof(double) * A);
        md_memcpy(out->part + j * A, part + k * A, sizeof(double) * A);
        j++;
    }
    out->nto_dropped = (float)(1.0 - lambda_kept / lambda_all);
    return true;
}

// The hole and particle densities of one excited state for the charge transfer window's side views, from its NTO
// pairs: D_hole = sum_k lambda_k h_k h_k^T, D_particle = sum_k lambda_k p_k p_k^T, in the AO basis. The same pairs the
// diagram partitions, so the views show what its columns count (for TDA these are the detachment and attachment
// densities; for RPA they differ slightly, which is why the views say "hole" and "particle"). The pairs are read on the
// main thread once they are resident (es_nto_ready), the matrices made on a worker: switching state never waits.
struct CtDensityJob {
    int32_t state = -1;
    size_t  A = 0;
    NtoPairs pairs = {};
    double* D[2] = {};          // {A,A} hole, particle
    bool    ok = false;
    std::atomic<bool> done = false;
    task_system::Job task;
};

static void ct_density_job_destroy(CtDensityJob* job) {
    if (!job) return;
    job->task.cancel_and_wait();
    md_allocator_i* heap = md_get_heap_allocator();
    const size_t K = job->pairs.K, A = job->A;
    if (job->pairs.lambda) md_free(heap, job->pairs.lambda, sizeof(double) * K);
    if (job->pairs.hole)   md_free(heap, job->pairs.hole,   sizeof(double) * K * A);
    if (job->pairs.part)   md_free(heap, job->pairs.part,   sizeof(double) * K * A);
    for (int i = 0; i < 2; ++i) if (job->D[i]) md_free(heap, job->D[i], sizeof(double) * A * A);
    delete job;
}

// The number of AOs the NTOs are written in: the last axis of the {S,L,A} coefficients
static size_t nto_num_ao(const md_system_t& sys) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, es_path::nto_hole);
    if (!attr || attr->format.rank < 1) return 0;
    return attr->format.shape[attr->format.rank - 1];
}

// Null when the state's pairs cannot be read (yet: ask es_nto_ready first)
static CtDensityJob* ct_density_job_start(const md_system_t& sys, int32_t state) {
    const size_t A = nto_num_ao(sys);
    if (state < 0 || A == 0) return nullptr;
    CtDensityJob* job = new CtDensityJob();
    job->state = state;
    job->A = A;
    if (!nto_pairs_read(&job->pairs, sys, A, state)) {
        delete job;
        return nullptr;
    }
    job->task.start({}, [job]() {
        const size_t A = job->A, K = job->pairs.K;
        md_allocator_i* heap = md_get_heap_allocator();
        const double* src[2] = { job->pairs.hole, job->pairs.part };
        for (int c = 0; c < 2; ++c) {
            double* D = (double*)md_alloc(heap, sizeof(double) * A * A);
            if (!D) { job->done.store(true, std::memory_order_release); return; }
            md_memset(D, 0, sizeof(double) * A * A);
            job->D[c] = D;
            for (size_t k = 0; k < K; ++k) {
                if (task_system::cancelled()) { job->done.store(true, std::memory_order_release); return; }
                const double lam = MD_MAX(0.0, job->pairs.lambda[k]);
                const double* v = src[c] + k * A;
                for (size_t i = 0; i < A; ++i) {
                    const double vi = lam * v[i];
                    if (vi == 0.0) continue;
                    double* row = D + i * A;
                    for (size_t j = 0; j < A; ++j) row[j] += vi * v[j];
                }
            }
        }
        job->ok = true;
        job->done.store(true, std::memory_order_release);
    });
    return job;
}

// One ribbon's orbital pair for the side views: a block NTO (md_qm_ct_block_nto), which the main thread forms
// in the Lowdin basis from the state's charge transfer data (K x K work and a pass over K x A), and this takes
// back to AO coefficients, C C~^T x, for the evaluation on the grid. O(M A) per vector: a worker's.
struct CtOrbitalJob {
    uint64_t key = 0;
    const LowdinJob* lowdin = nullptr;  // outlives this job
    double*  x[2] = {};                 // [A] hole and particle in the Lowdin basis, owned
    double*  c[2] = {};                 // [A] the same in AO coefficients, owned
    char     info[2][64] = {};
    char     sub[96] = {};
    bool     ok = false;
    std::atomic<bool> done = false;
    task_system::Job task;
};

static void ct_orbital_job_destroy(CtOrbitalJob* job) {
    if (!job) return;
    job->task.cancel_and_wait();
    md_allocator_i* heap = md_get_heap_allocator();
    const size_t A = job->lowdin ? job->lowdin->num_ao : 0;
    for (int i = 0; i < 2; ++i) {
        if (job->x[i]) md_free(heap, job->x[i], sizeof(double) * A);
        if (job->c[i]) md_free(heap, job->c[i], sizeof(double) * A);
    }
    delete job;
}

// x[0], x[1] are taken over: allocated by the caller from the heap, A doubles each
static CtOrbitalJob* ct_orbital_job_start(const LowdinJob* lowdin, uint64_t key, double* x_hole, double* x_part) {
    md_allocator_i* heap = md_get_heap_allocator();
    CtOrbitalJob* job = new CtOrbitalJob();
    job->key    = key;
    job->lowdin = lowdin;
    job->x[0]   = x_hole;
    job->x[1]   = x_part;
    const size_t A = lowdin->num_ao;
    job->c[0] = (double*)md_alloc(heap, sizeof(double) * A);
    job->c[1] = (double*)md_alloc(heap, sizeof(double) * A);
    job->task.start(MD_STR_LIT("Ribbon orbitals"), [job]() {
        const LowdinJob* lw = job->lowdin;
        for (int i = 0; i < 2; ++i) {
            if (task_system::cancelled()) { job->done.store(true, std::memory_order_release); return; }
            md_qm_lowdin_inverse_transform(job->c[i], job->x[i], 1, lw->C, lw->C_lowdin, lw->num_mo, lw->num_ao);
        }
        job->ok = true;
        job->done.store(true, std::memory_order_release);
    });
    return job;
}

// The NTO pairs are read here, on the main thread, from the attribute table; the rest runs on a worker.
static CtStateJob* ct_state_job_start(const md_system_t& sys, const LowdinJob* lowdin, int32_t state) {
    if (!lowdin_job_ready(lowdin) || state < 0) return nullptr;
    NtoPairs pairs = {};
    if (!nto_pairs_read(&pairs, sys, lowdin->num_ao, state)) return nullptr;

    CtStateJob* job = new CtStateJob();
    job->state       = state;
    job->lowdin      = lowdin;
    job->K           = pairs.K;
    job->num_atoms   = lowdin->num_atoms;
    job->lambda      = pairs.lambda;
    job->hole        = pairs.hole;
    job->part        = pairs.part;
    job->nto_dropped = pairs.nto_dropped;

    job->task.start({}, [job]() {
        const LowdinJob* lw = job->lowdin;
        const size_t A = lw->num_ao, M = lw->num_mo, N = job->num_atoms, K = job->K;
        md_allocator_i* heap = md_get_heap_allocator();

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        double* ht = md_temp_alloc_array(temp, double, K * A);
        double* pt = md_temp_alloc_array(temp, double, K * A);
        md_qm_lowdin_transform(ht, job->hole, K, lw->S, lw->C, lw->C_lowdin, M, A);
        if (task_system::cancelled()) { job->done.store(true, std::memory_order_release); return; }
        md_qm_lowdin_transform(pt, job->part, K, lw->S, lw->C, lw->C_lowdin, M, A);
        md_memcpy(job->hole, ht, sizeof(double) * K * A);
        md_memcpy(job->part, pt, sizeof(double) * K * A);

        job->hole_blocks = (double*)md_alloc(heap, sizeof(double) * N * K * K);
        job->part_blocks = (double*)md_alloc(heap, sizeof(double) * N * K * K);
        job->omega       = (double*)md_alloc(heap, sizeof(double) * N * N);
        job->weight[0]   = (double*)md_alloc(heap, sizeof(double) * N);
        job->weight[1]   = (double*)md_alloc(heap, sizeof(double) * N);
        if (!job->hole_blocks || !job->part_blocks || !job->omega || !job->weight[0] || !job->weight[1]) {
            MD_LOG_ERROR("Charge transfer: out of memory (%zu atoms, %zu NTO pairs)", N, K);
            job->done.store(true, std::memory_order_release);
            return;
        }
        md_qm_ct_group_blocks(job->hole_blocks, job->hole, K, A, lw->ao_atom, N);
        md_qm_ct_group_blocks(job->part_blocks, job->part, K, A, lw->ao_atom, N);
        md_qm_ct_omega(job->omega, job->lambda, K, job->hole_blocks, N, job->part_blocks, N);

        double total = 0.0;
        for (size_t i = 0; i < N * N; ++i) total += job->omega[i];
        if (total <= 0.0) { job->done.store(true, std::memory_order_release); return; }

        // The atom to atom tail is not drawn below a floor. Dropping and renormalizing what is left keeps
        // the graph exactly conserved, and what was dropped is reported.
        const double floor = 1.0e-6 * total;
        double kept = 0.0;
        for (size_t i = 0; i < N * N; ++i) {
            if (job->omega[i] < floor) job->omega[i] = 0.0;
            kept += job->omega[i];
        }
        job->link_dropped = (float)(1.0 - kept / total);
        md_memset(job->weight[0], 0, sizeof(double) * N);
        md_memset(job->weight[1], 0, sizeof(double) * N);
        for (size_t a = 0; a < N; ++a) {
            for (size_t b = 0; b < N; ++b) {
                const double w = job->omega[a * N + b] / kept;
                job->omega[a * N + b] = w;
                job->weight[0][a] += w;
                job->weight[1][b] += w;
            }
        }
        job->ok = true;
        job->done.store(true, std::memory_order_release);
    });
    return job;
}

// One state's entry in the state strip: its charge transfer between the automatic UNITS (a fixed level,
// so states compare), the distance between its hole and particle centroids, and how many NTO pairs it
// takes (the participation ratio of the lambdas). A light job: the transfer matrix over units, not atoms.
struct StripJob {
    int32_t state = -1;
    const LowdinJob* lowdin = nullptr;  // outlives this job
    NtoPairs pairs = {};
    const uint32_t* ao_unit = nullptr;  // [A], owned by the caller and kept until the job is destroyed
    size_t  num_units = 0;
    const double* atom_xyz = nullptr;   // [N*3] Angstrom, likewise
    float   ct = 0.0f;
    float   distance = 0.0f;
    float   pr = 0.0f;
    bool    ok = false;
    std::atomic<bool> done = false;
    task_system::Job task;
};

static void strip_job_destroy(StripJob* job) {
    if (!job) return;
    job->task.cancel_and_wait();
    md_allocator_i* heap = md_get_heap_allocator();
    const size_t A = job->lowdin ? job->lowdin->num_ao : 0, K = job->pairs.K;
    if (job->pairs.lambda) md_free(heap, job->pairs.lambda, sizeof(double) * K);
    if (job->pairs.hole)   md_free(heap, job->pairs.hole,   sizeof(double) * K * A);
    if (job->pairs.part)   md_free(heap, job->pairs.part,   sizeof(double) * K * A);
    delete job;
}

static StripJob* strip_job_start(const md_system_t& sys, const LowdinJob* lowdin, int32_t state,
                                 const uint32_t* ao_unit, size_t num_units, const double* atom_xyz) {
    if (!lowdin_job_ready(lowdin) || !ao_unit || num_units == 0) return nullptr;
    StripJob* job = new StripJob();
    job->state     = state;
    job->lowdin    = lowdin;
    job->ao_unit   = ao_unit;
    job->num_units = num_units;
    job->atom_xyz  = atom_xyz;
    if (!nto_pairs_read(&job->pairs, sys, lowdin->num_ao, state)) {
        job->done.store(true, std::memory_order_release);
        return job;
    }
    job->task.start({}, [job]() {
        const LowdinJob* lw = job->lowdin;
        const size_t A = lw->num_ao, M = lw->num_mo, N = lw->num_atoms, U = job->num_units, K = job->pairs.K;
        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        double* ht = md_temp_alloc_array(temp, double, K * A);
        double* pt = md_temp_alloc_array(temp, double, K * A);
        md_qm_lowdin_transform(ht, job->pairs.hole, K, lw->S, lw->C, lw->C_lowdin, M, A);
        if (task_system::cancelled()) { job->done.store(true, std::memory_order_release); return; }
        md_qm_lowdin_transform(pt, job->pairs.part, K, lw->S, lw->C, lw->C_lowdin, M, A);

        double* hb = md_temp_alloc_array(temp, double, U * K * K);
        double* pb = md_temp_alloc_array(temp, double, U * K * K);
        double* om = md_temp_alloc_array(temp, double, U * U);
        md_qm_ct_group_blocks(hb, ht, K, A, job->ao_unit, U);
        md_qm_ct_group_blocks(pb, pt, K, A, job->ao_unit, U);
        md_qm_ct_omega(om, job->pairs.lambda, K, hb, U, pb, U);
        double total = 0.0, local = 0.0;
        for (size_t a = 0; a < U; ++a) {
            for (size_t b = 0; b < U; ++b) total += om[a * U + b];
            local += om[a * U + a];
        }

        // Centroids from the per atom hole and particle, which decompose exactly per pair
        double rh[3] = {0, 0, 0}, rp[3] = {0, 0, 0}, sh = 0.0, sp = 0.0;
        if (job->atom_xyz) {
            for (size_t mu = 0; mu < A; ++mu) {
                const uint32_t atom = lw->ao_atom[mu];
                if (atom >= N) continue;
                double qh = 0.0, qp = 0.0;
                for (size_t k = 0; k < K; ++k) {
                    qh += job->pairs.lambda[k] * ht[k * A + mu] * ht[k * A + mu];
                    qp += job->pairs.lambda[k] * pt[k * A + mu] * pt[k * A + mu];
                }
                for (int d = 0; d < 3; ++d) {
                    rh[d] += qh * job->atom_xyz[atom * 3 + d];
                    rp[d] += qp * job->atom_xyz[atom * 3 + d];
                }
                sh += qh;
                sp += qp;
            }
        }
        double l1 = 0.0, l2 = 0.0;
        for (size_t k = 0; k < K; ++k) { l1 += job->pairs.lambda[k]; l2 += job->pairs.lambda[k] * job->pairs.lambda[k]; }

        job->ct = total > 0.0 ? (float)(1.0 - local / total) : 0.0f;
        if (sh > 0.0 && sp > 0.0) {
            double d2 = 0.0;
            for (int d = 0; d < 3; ++d) { const double x = rh[d] / sh - rp[d] / sp; d2 += x * x; }
            job->distance = (float)sqrt(d2);
        }
        job->pr = l2 > 0.0 ? (float)(l1 * l1 / l2) : 0.0f;
        job->ok = total > 0.0;
        job->done.store(true, std::memory_order_release);
    });
    return job;
}

// Everything below reads the SYSTEM's attribute table. There is no reader object here any more and
// no second parse of the file: whatever loaded the system published what it carried, and these
// windows are consumers of that table like any other - which is what lets a system published by
// some other quantum chemistry reader drive the same windows without this file knowing it exists.
//
// The paths are named in one place: es_path in viamd.h for the format neutral ones, and the vlx/
// literals below for what is still VeloxChem's own vocabulary waiting for a neutral name.
struct QuantumChemistry : viamd::EventHandler {
    QuantumChemistry() { viamd::event_system_register_handler(*this); }

    bool use_gpu_path = false;
#if MD_ENABLE_GPU
    // No GPU handles are kept here any more. The device, its stream, its pools and the evaluation
    // scratch belong to the application; the uploaded basis and the atom buffer belong to the
    // dataset. Everything this component does with them goes through the application's evaluators,
    // which read them from the ApplicationState it is already handed - so a stale copy taken at
    // load time is not something that can exist.
    //
    // The one exception is the critical point extraction below, which records its own kernel over
    // the scratch volume. It reads state.gpu_* at the point of use, for the same reason.
#endif

    // GL representations
    md_gl_rep_t gl_rep = {};

    int homo_idx[2] = {};
    int lumo_idx[2] = {};

    AABB aabb = {};
    OABB oabb = {};

    dvec3_t center_of_charge = {};

    // Up to date packed atom xyz in BOHR ready to be supplied for evaluation
    // The w component is padding for now. 16-bytes are used for better alignment (gpu buffer 1:1)
    vec4_t* atom_xyzw = nullptr;

    // This is the default view which is used as a reset view target
    ViewTransform default_view = {};

    // How many atoms the calculation covered, and nothing else about them.
    //
    // The QM ATOM DOMAIN itself (see QmAtoms in viamd.h) is NOT cached here. es_qm_atoms hands back
    // pointers into the attribute table's own storage, and md_attributes_replace frees the old
    // buffer before it writes the new one - so a supplemental load, which republishes qm/atom/*
    // onto an already loaded system, invalidates them. Keeping them alive would mean re-gathering
    // on a signal from a DIFFERENT attribute's version, which is a coupling that breaks silently.
    // Every site that reads the domain has the system in hand and hoists the pointers out of its
    // own loop anyway, so it gathers them itself: three binary searches, once per frame.
    //
    // The COUNT survives because two handlers have no system to ask - the menu event carries no
    // payload at all, and serialize_workspace takes only a serialization state. It is a plain
    // number, so nothing it names can be freed underneath it.
    size_t num_qm_atoms = 0;

    // The {M,A} shape of the coefficient attribute, cached beside it: the per frame window gate asks
    // for it and the answer only changes with the system.
    size_t num_mos = 0;

    // Excited states carrying NTOs. Cached for the same reason, and because the menu handler has no
    // ApplicationState to reach the table through.
    size_t num_excited = 0;

    // Version of the coefficient attribute the setup above last ran for. See the load event.
    uint64_t qm_version = 0;
    md_bitset_t* qm_mask = nullptr;

    struct Summary {
        bool show_window = false;
    } summary;

    struct CriticalPoints {
        bool show_window = false;       // its window (draw_critical_points_window)
        bool enabled = false;

        // How they are found. Analytic searches the GTO density itself and proves what it reports
        // (md_topo_compute_extremum_graph_gto); Volume samples the density onto a grid first and can
        // place a critical point no better than to the voxel it falls in.
        enum Method : int { Method_Analytic = 0, Method_Volume = 1 };
        int method = Method_Analytic;

        uint64_t raw_hash = 0;
        uint64_t simp_hash = 0;

        // Analytic: the threshold, the run in flight, the inputs wanted now and the inputs the shown
        // graph was computed for, and what that run reported about itself.
        float                rho_min     = 1.0e-4f;
        bool                 rho_min_editing = false; // its slider is held: no search is started for every value
        bool                 gto_use_gpu = true;      // only consulted when the app has a GPU device
        GtoCriticalPointJob* gto_job     = nullptr;
        uint64_t             gto_want    = 0;
        uint64_t             gto_have    = 0;
        md_topo_gto_info_t   gto_info    = {};
        bool                 gto_failed  = false;
        double               gto_seconds = 0.0;
        double               gto_rho_min = 0.0;       // the threshold the shown graph was searched with

        // Volume: the grid, and what is left out of the graph of its samples
        VolumeResolution vol_res            = VolumeResolution::Mid;
        bool             vol_per_frame      = false;  // sampled and extracted again every frame
        bool             vol_simplify       = true;
        float            vol_simp_threshold = 1.0E-4f;

#if MD_ENABLE_GPU
        md_topo_gpu_context_t* topo_ctx      = nullptr;
        uint32_t               topo_dims[3]  = {};
#endif
        md_topo_extremum_graph_t raw_graph  = { 0 };
        md_topo_extremum_graph_t simp_graph = { 0 };

        // How they are drawn (draw_critical_points)
        float point_scale = 1.0f;       // the markers' size, times critical_point_styles' radii
        bool  show_paths  = true;       // the bond paths and RCP attractor paths (the graph's edges)

        PickingRange picking_range = {};

        // Over the critical points, in the persistent allocator (initialized with the component)
        md_bitmask_t selection_mask = {};
        md_bitmask_t highlight_mask = {};
    } critical_points;

    struct Opt {
        int hovered = -1;
        int selected = -1;

        bool coord_modified = false;

        uint64_t energy_hash = 0;
    } opt;

    struct Scf {
        // Convergence history the SCF plot last fitted its axes to. See draw_summary_window.
        uint64_t plot_hash = 0;
    } scf;

    struct Orb {
        bool show_window = false;
        Volume   vol[16] = {};
        int      vol_mo_idx[16] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        SpinChannel vol_mo_type[16] = {};
        uint32_t iso_tex[16] = {};      // RGBA16F: the panel's isosurfaces, premultiplied HDR
        uint32_t out_tex[16] = {};      // RGBA8: the finished panel
        uint32_t out_fbo = 0;
        uint64_t render_hash = 0;       // of everything the panels were last rendered from: they are redrawn when it changes
        int num_x  = 3;
        int num_y  = 3;
        int mo_idx = -1;
        int scroll_to_idx = -1;

        // Active sort of the orbital table. Mo index descending is what the table showed before
        // it was sortable, so it stays the default.
        int  sort_column    = OrbitalColumn_Idx;
        bool sort_ascending = false;

        PanelMeshes mesh = {};

        IsoDesc iso = {
            .count = 2,
            .values = {0.05f, -0.05},
            .colors = {{0.f, 75.f/255.f, 135.f/255.f, 0.75f}, {1.0f, 205.f/255.f, 0.f, 0.75f}},
        };

        GBuffer gbuf = {};
        PickingSurface picking_surface = {};
        Camera camera = {};
        ViewTransform target;

        float distance_scale = 2.0f;
        bool show_coordinate_system_widget = true;
    } orb;

    struct Nto {
        bool show_window = false;
        Volume   vol[2] = {};
        uint32_t iso_tex[2] = {};       // RGBA16F: the panel's isosurfaces, premultiplied HDR
        uint32_t out_tex[2] = {};       // RGBA8: the finished panel
        uint32_t out_fbo = 0;
        uint64_t render_hash = 0;       // of everything the panels were last rendered from: they are redrawn when it changes
        PanelMeshes mesh = {};

        int sel_nto_idx = -1;

        md_grid_t grid = {0};

        md_array(uint32_t) atom_group_idx = nullptr;

        md_gl_rep_t gl_rep = {0};

        size_t transition_matrix_dim = 0;
        float* transition_matrix = nullptr;
        float* transition_density_hole = nullptr;
        float* transition_density_part = nullptr;

        struct {
            size_t count = 0;

            char   label[MAX_NTO_GROUPS][64] = {};
            vec4_t color[MAX_NTO_GROUPS] = {};
            int8_t hovered_index = -1;
        } group;

        double iso_val = 0.05;
        vec4_t col_pos = { 0.0f, 0.294f, 0.529f, 0.75f };
        vec4_t col_neg = { 1.0f, 0.804f, 0.0f,   0.75f };
        vec4_t col_den = { 1.0f, 1.0f,   1.0f,   0.75f };
        vec4_t col_att = VEC4_ATTACHMENT;
        vec4_t col_det = VEC4_DETACHMENT;

        struct {
            bool enabled = false;
            vec4_t colors[3] = {
                {0.0f, 1.0f, 1.0f, 1.0f},
                {1.0f, 0.0f, 1.0f, 1.0f}
            };
            float vector_scale = 1.0f;
            bool show_angle = false;
        } dipole;

        GBuffer gbuf = {};
        PickingSurface picking_surface = {};
        Camera camera = {};
		ViewTransform target = {};

        float distance_scale = 2.0f;
        bool link_attachment_detachment_density = false;
        bool show_coordinate_system_widget = true;
    } nto;

    // One excited state's attachment and detachment densities in the AO basis.
    //
    // Both are VIRTUAL attributes: every read rebuilds the matrix from the response solution
    // vector and the MO coefficients, which is by far the most expensive thing this component
    // does. Selecting a state asks for exactly the same two matrices twice in the same frame -
    // once for the Mulliken group attribution in update_nto_derived_data, once for the two
    // volumes the Transition Analysis window evaluates - and the export window asks for a third
    // copy. Rebuilt once here and held until the state or the dataset changes.
    //
    // Allocated from 'arena'; reset_data() rewinds that arena, so the pointers are dropped rather
    // than freed there. The difference density is not cached: it is one subtraction away from the
    // two that are, and only the export window asks for it.
    struct TransitionDensityCache {
        int     state_idx = -1;     // which excited state, -1 when nothing is cached
        size_t  dim       = 0;      // atomic orbitals per side
        double* attach    = nullptr;
        double* detach    = nullptr;
    } td_cache;

    // The charge transfer window's side views: the molecule with the hole density beside the diagram's hole column and
    // the particle density beside its particle column, one camera for both. Optional (they take room).
    //
    // Nothing here waits. A state's two densities are made from its NTO pairs on a worker (CtDensityJob) and evaluated
    // into volumes on the GPU, at half the resolution of the main grid (the views are small). The volumes are kept in a
    // cache of a set budget, least recently used out first, and while nothing else is asked the cache is filled with
    // the other states, the selected one's neighbours first: switching state shows what is cached at once, and what is
    // not yet keeps the last picture until it is.
    struct CtViews {
        bool   show = false;
        float  width_frac = 0.26f;      // of the window's canvas, each side
        double iso_val = 0.05;          // drawn at iso_val^2, as the transition analysis draws its densities
        int    budget_mb = 256;
        md_grid_t grid = {};
        bool   initialized = false;

        // A state's two densities, or one ribbon's orbital pair (block NTO) or the state's own NTO pair: the
        // volumes the views show, under one budget
        struct Entry {
            uint64_t key = 0;           // ctv_density_key, or an orbital key (top bit set)
            int32_t  state = -1;
            Volume   vol[2] = {};       // hole, particle
            uint64_t v0[2] = {};        // their data versions when made: filled once both have moved on
            uint64_t used = 0;
            bool     failed = false;    // nothing to show (no NTO pairs, or the evaluation failed)
            bool     orbital = false;   // signed: drawn at +-iso rather than at iso^2
            char     info[2][64] = {};  // per view, beside its caption: what is shown
            char     sub[96] = {};      // under both
        };
        md_array(Entry) cache = nullptr;    // heap
        uint64_t tick = 0;
        CtDensityJob* job = nullptr;
        int32_t shown = -1;             // the state whose volumes the views show

        // Ribbon orbitals. Hovering a ribbon asks for its orbital pair once the cursor has rested on it for
        // hover_delay (sweeping across the diagram asks for nothing); clicking one pins it, so that the views
        // keep it while the cursor is off the diagram (turning them, for one). The orbital route's NTO pairs
        // likewise, from the middle column.
        CtOrbitalJob* orb_job = nullptr;
        uint64_t orb_shown = 0;         // the orbital entry last shown
        uint64_t hover_cand = 0;        // the orbital under the cursor, and since when
        double   hover_t0 = 0.0;
        float    hover_delay = 0.15f;
        struct Pin {
            bool     on = false;
            int32_t  state = -1;
            uint64_t src_key = 0, dst_key = 0;  // the ribbon, by its ends (fragment route)
            uint32_t pair = 0;
            int32_t  nto = -1;                  // the NTO pair (orbital route)
        } pin;

        GBuffer  gbuf = {};
        uint32_t iso_tex[2] = {};
        uint32_t out_tex[2] = {};
        uint32_t out_fbo = 0;
        uint64_t render_hash = 0;
        PickingSurface picking_surface = {};
        Camera camera = {};
        ViewTransform target = {};
        md_gl_rep_t gl_rep = {};
    } ctv;

    void ctv_free() {
        ct_density_job_destroy(ctv.job);
        ctv.job = nullptr;
        ct_orbital_job_destroy(ctv.orb_job);
        ctv.orb_job = nullptr;
        for (size_t i = 0; i < md_array_size(ctv.cache); ++i) {
            for (int c = 0; c < 2; ++c) if (ctv.cache[i].vol[c].tex_id) gl::free_texture(&ctv.cache[i].vol[c].tex_id);
        }
        md_array_free(ctv.cache, md_get_heap_allocator());
        for (int i = 0; i < 2; ++i) {
            if (ctv.iso_tex[i]) gl::free_texture(&ctv.iso_tex[i]);
            if (ctv.out_tex[i]) gl::free_texture(&ctv.out_tex[i]);
        }
        if (ctv.out_fbo) glDeleteFramebuffers(1, &ctv.out_fbo);
        if (ctv.gbuf.fbo) gbuffer_free(&ctv.gbuf);
        md_gl_rep_destroy(ctv.gl_rep);
        const CtViews prefs = { .show = ctv.show, .width_frac = ctv.width_frac, .iso_val = ctv.iso_val, .budget_mb = ctv.budget_mb };
        ctv = CtViews{};
        ctv.show       = prefs.show;        // preferences, not data
        ctv.width_frac = prefs.width_frac;
        ctv.iso_val    = prefs.iso_val;
        ctv.budget_mb  = prefs.budget_mb;
    }

    // Charge Transfer Analysis: the layered flow diagram (docs/transition_flow_design.md).
    //
    // Separate window, separate state, but the SAME group definitions as the NTO window - asking
    // someone to define their groups twice in two windows would be worse than any coupling this
    // introduces.
    struct Flow {
        bool show_window = false;

        flow_mode_t mode = FLOW_MODE_FRAGMENTS;

        md_flow_graph_t  graph  = {};
        md_flow_cut_t    cut    = {};
        md_flow_layout_t layout = {};
        md_flow_layout_params_t params = {};

        // Node labels are borrowed by the graph rather than owned, so they have to outlive it.
        // Group labels already live in nto.group.label; these two are ours.
        char nto_label[MAX_NTO_LABELS][40] = {};
        uint32_t nto_index[MAX_NTO_LABELS] = {};   // the state's NTO pair each middle node stands for

        // Second line: the dominant MO character of each NTO, e.g. "HOMO-1 -> LUMO". Empty when
        // no single pair dominates enough to name honestly.
        char nto_character[MAX_NTO_LABELS][40] = {};

        // Fixed-stride blocks, allocated once per file. The graph BORROWS these strings, so they
        // have to outlive it; rebuilding them per graph rebuild would leak into the arena on every
        // state change. Atom labels never change at all; MO labels are rewritten in place.
        char*  atom_labels = nullptr;
        size_t num_atom_labels = 0;
        char*  mo_labels = nullptr;
        size_t num_mo_labels = 0;

        // Most atoms contribute almost nothing to a given transition, so the diagram opens with a
        // threshold rather than with a hundred hairline ribbons. The slider is in the toolbar.
        float initial_threshold = 0.01f;

        // How hard to flatten the ORBITAL columns' heights. 0 = strictly proportional, 1 = nearly
        // uniform. Applied only to orbital columns: an atom column's proportions ARE the charge
        // transfer story, while an NTO column is a channel picker, and a channel drawn two pixels
        // tall cannot be read, labelled or clicked.
        float size_normalize = 0.75f;

        // NTO pairs below this share of the transition are dropped and the rest renormalized.
        // Distinct from the cut threshold, which hides nodes that are still in the graph.
        float lambda_cutoff = 1.0e-3f;

        // What the last Fragments build left out, so the window can say it. The pairs below a round
        // off floor of the excitation (interference makes every other pair count), and the atom to
        // atom entries too small to draw - between distant atoms there are many of them.
        float nto_dropped  = 0.0f;
        float link_dropped = 0.0f;

        // Rows per column up to which the diagram opens fragments on its own (md_flow_cut_t.max_rows)
        int max_rows = 12;

        // The system's fragments (sys->fragment, md_fragment.h) as a tree to show (fragment_tree.h), and the
        // two atom index spaces: QM atoms (what the basis, the NTOs and the groups are over) and system atoms
        // (what fragments and selections are over). Made again when the topology changes (the fragments come
        // with the system's inference, which may finish after the window opens); the arena reset frees them.
        FragmentTree frags = {};
        uint64_t  frags_version = 0;       // the topology version the tree was made from
        uint32_t* qm_to_sys = nullptr;     // [num_qm_atoms]
        int32_t*  sys_to_qm = nullptr;     // [num_sys_atoms], -1 for an atom outside the calculation
        size_t    num_qm_atoms  = 0;
        size_t    num_sys_atoms = 0;
        bool      frags_built   = false;
        int       ct_failed_state = -2;    // the state whose NTO pairs could not be read, not to be retried

        // The state strip: one entry per excited state, filled one state at a time (StripJob)
        struct StripEntry {
            float ct = 0.0f, distance = 0.0f, pr = 0.0f;
            bool  done = false, ok = false;
        };
        md_array(StripEntry) strip = nullptr;
        uint32_t* ao_unit = nullptr;       // [A] the unit (leaf fragment, compacted) of each AO
        size_t    num_units = 0;
        double*   atom_xyz = nullptr;      // [N*3] the QM atoms, Angstrom
        bool      show_strip = true;
        bool      show_legend = true;
        float     bar_width = 0.03f;    // the layout's node thickness, graph units: thin bars, labels outside

        // The orbital pairs (block NTOs) of every visible ribbon, as shares of it: see FlowStripes
        md_array(uint32_t) stripe_offset = nullptr;
        md_array(float)    stripe_frac   = nullptr;
        uint64_t           stripe_hash   = 0;

        FlowDrawStyle style = {};

        // View transform, in graph units. Kept here rather than in the layout so that panning and
        // zooming never trigger a relayout.
        vec2_t pan  = {0, 0};
        float  zoom = 1.0f;

        // Last frame's hit, so the dimming does not flicker while the cursor crosses a seam.
        FlowHit hover = {};

        // Deferred because fitting needs a laid-out graph and a canvas size, neither of which the
        // rebuild has. Consumed on the next draw.
        bool fit_requested = true;

        // Rubber-band region select. The button is captured at press time rather than read live,
        // so a drag started with the left button stays an ADD even if the right one gets pressed
        // part way through.
        bool   region_active = false;
        bool   region_add    = true;
        ImVec2 region_start  = {0, 0};

        // ATOM selection is NOT held here. It lives in state.selection.selection_mask, and this
        // window reads it every frame - so a selection made in the 3D view, by a script, or in any
        // other window shows up here, and one made here shows up everywhere else. A private copy
        // would be a second source of truth that silently drifts from the first.
        //
        // ORBITAL selection - NTO pairs and canonical MOs - is a different quantity: it picks a
        // transition channel rather than atoms, so it has nowhere else to live and is held here,
        // by key, like expansion.
        md_array(uint64_t) selected_orbital_keys = nullptr;

        uint64_t data_hash = 0;
        bool built = false;
    } flow;

    // Shared by every excited state of the system; made on first use by the Charge Transfer window
    LowdinJob* lowdin = nullptr;
    // The selected state's numbers per atom; replaced when another state is selected
    CtStateJob* ct_job = nullptr;
    // The state strip's current state, one at a time
    StripJob* strip_job = nullptr;

    // Workspace group assignments, waiting for a molecule.
    //
    // load_workspace parses the ENTIRE workspace text before it loads any data, so when the
    // [QuantumChemistry] section is read there is no system and no atom count yet. Groups are
    // per-atom, so they have to be buffered and stamped on at the end of init_qm_from_system.
    //
    // The buffer lives on the heap, NOT on 'arena': reset_data() resets that arena, and a
    // workspace load goes through reset_data before it gets anywhere near init_qm_from_system.
    struct PendingGroups {
        bool   pending = false;
        size_t count = 0;
        char   label[MAX_NTO_GROUPS][64] = {};
        vec4_t color[MAX_NTO_GROUPS] = {};
        md_bitset_t atoms[MAX_NTO_GROUPS] = {};
        bool   atoms_initialized = false;
    } pending_groups;

    struct Rsp {
        bool show_window = false;
        bool show_export_window = false;
        int hovered = -1;
        int selected = -1;

        // Broadening parameter
        double broadening_fwhm = 0.123; // FWHM in eV, default value is 0.123 eV
        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        x_unit_t x_unit = X_UNIT_EV;

        uint64_t hash = 0; // Hash for tracking the unit of the x-axis, if it changes we need to recalculate the spectra (and shown peaks)
    } rsp;

    // Settings + cached results for the RIXS map.
    //
    // The struct owns its own arena ('alloc'), created at ViamdInitialize and released at
    // ViamdShutdown. The cached arrays below are the *only* thing allocated from it, and
    // compute_rixs_map() rebuilds all of them from scratch, so the arena is simply reset at the top
    // of that function: memory never exceeds the footprint of a single map, and no array can be
    // leaked by forgetting to free it individually.
    //
    // A vm arena is used rather than md_arena_allocator because reset() is a pointer rewind that
    // keeps the pages committed. The map is rebuilt on every settings change (including while a
    // slider is being dragged), and md_arena_allocator_reset() would release every page back to the
    // backing allocator only to re-acquire it on the next frame.
    //
    // NOTE: 'alloc' must survive the value-initialization of this struct in reset_data().
    struct Rixs {
        // --- Settings (mirrors the keyword arguments of spectrumplot.plot_rixs_map) ---

        // Include the elastic (Rayleigh) line as an additional stick per photon energy.
        bool plot_elastic_line = false;

        // Normalize the map by its global maximum ('global_max' in the python version).
        bool normalize = true;

        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        double broadening_fwhm_ev = 0.24;   // FWHM in eV applied along the energy loss axis
        double x_step_ev = 0.01;            // Grid spacing along the energy loss axis, in eV
        ImPlotColormap colormap = ImPlotColormap_Viridis;

        // Arena backing every pointer below. Owned by the module, see the comment above the struct.
        md_allocator_i* alloc = nullptr;

        // --- Cached results, rebuilt whenever 'hash' changes ---
        uint64_t hash  = 0;
        bool     valid = false;

        // The map, row-major [num_rows][num_cols].
        // NOTE: row 0 holds the HIGHEST photon energy. ImPlot::PlotHeatmap renders the first row at
        // bounds_max.y (top), so reversing the row order here reproduces matplotlib's origin='lower'.
        double* map = nullptr;

        double* grid      = nullptr;  // [num_cols] energy loss in eV
        double* photon_ev = nullptr;  // [num_rows] incoming photon energies in eV, ascending
        double* xas_x     = nullptr;  // [num_xas]  broadened XAS energies in eV
        double* xas_y     = nullptr;  // [num_xas]  broadened XAS intensities
        double* core_ev   = nullptr;  // [num_core] core excitation energies in eV, ascending
        double* core_f    = nullptr;  // [num_core] core oscillator strengths, sorted with core_ev

        size_t num_rows = 0;  // == number of incoming photon energies (P)
        size_t num_cols = 0;  // == number of grid samples along the energy loss axis (N)
        size_t num_xas  = 0;
        size_t num_core = 0;

        double map_min = 0.0, map_max = 0.0;
        double x_min_ev = 0.0, x_max_ev = 0.0;  // extent of the energy loss axis
        double y_min_ev = 0.0, y_max_ev = 0.0;  // extent of the photon energy axis (pixel edges)
        double xas_max  = 0.0;

        // Linked photon-energy axis, shared between the map and the XAS panel
        double y_link_min = 0.0, y_link_max = 1.0;

        int hovered_core  = -1;
        int selected_core = -1;

        // --- RIXS spectrum (a single incoming photon energy, i.e. one row of the map) ---
        // The spectrum is cheap enough to rebuild every frame, so nothing is cached here beyond the
        // selection state. The broadening settings above are shared with the map.
        int selected_photon = 0;   // index into the ascending photon energy list
        int hovered_final   = -1;
        int selected_final  = -1;
        uint64_t spectrum_hash = 0;  // only used to refit the axes when the plotted data changes

        // Set once the non-uniform photon spacing has been reported, to avoid spamming the log while
        // a settings slider is being dragged.
        bool warned_non_uniform = false;
    } rixs;

    // Settings + interaction state for the XPS plot.
    // No cached arrays: the broadened curve is generated on the fly by an ImPlotGetter and the peak
    // arrays are a few tens of doubles built into the frame temp arena, so there is nothing worth
    // keeping between frames.
    struct Xps {
        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        double broadening_fwhm_ev = 0.5;    // Matches the default of plot_xps_spectrum()

        bool show_sticks = true;
        bool invert_x    = false;   // XPS binding energy axes are often drawn decreasing to the right

        // false -> every core hole contributes a stick of height 1 (what VeloxChem's plot does)
        // true  -> stick height is the atom's share of the core MO
        // These differ only for delocalized holes, see the comment in draw_xps_plot().
        bool weight_by_contribution = false;

        // Which element is on screen. One at a time: a peak is now an ATOM rather than an element,
        // and core binding energies of different elements are hundreds of eV apart, so putting two
        // elements on one axis leaves each of them an unreadable spike. Zero means nothing has been
        // picked yet - the first group with entries is adopted on the first draw.
        md_element_t element = 0;

        // Index into the FLAT entry array, not into a group. -1 when the cursor is not on a peak.
        // There is deliberately no 'selected' counterpart: a peak draws as selected when its atom is
        // in the application selection mask, so the plot is a view of that mask and cannot drift
        // out of sync with the viewport or with a selection made anywhere else.
        int hovered = -1;

        // MO of the hovered peak, handed to the orbital grid window so it can highlight the matching
        // table row, or -1 when no peak is hovered. Written by draw_xps_plot, consumed and cleared by
        // draw_orb_window - see the comment at the top of that function for the ordering.
        int32_t highlight_mo_idx = -1;

        uint64_t hash = 0;  // Refit the axes when the settings change
    } xps;

    struct Vib {
        int hovered = -1;
        int selected = -1;

        // In case of resonance raman, this is the index of the selected electronic transition
        int external_frequency_index = 0;

        double broadening_fwhm = 20.0; // FWHM in cm^-1, default value is 20 cm^-1
        broadening_mode_t broadening_mode = BROADENING_MODE_LORENTZIAN;
        bool coord_modified = false;
        float displacement_amp_scl  = 1.0f;
        float displacement_freq_scl = 1.0f;
        bool invert_x = true;
        bool invert_y = false;
        bool displace_aos = false;

        // Time accumulation for vibrational mode visualization (pertubation of atoms)
        double t = 0;

        // Scaling applied to input frequencies (to account for inaccuracies of the chosen basis set)
        double freq_scaling_factor = 1.0;

        // Hash to track if the underlying peaks have changed
        uint64_t hash = 0;
    } vib;

    struct Export {
        char export_path[1024] = {};
        ElectronicStructureSource source = ElectronicStructureSource::MolecularOrbital;
        bool use_magnitude = false;
        ElectronicStructureSpin spin = ElectronicStructureSpin::Total;
        ElectronicStructureNtoComponent nto_component = ElectronicStructureNtoComponent::Particle;
        ElectronicStructureTransitionDensityComponent transition_density_component = ElectronicStructureTransitionDensityComponent::Attachment;
        VolumeResolution resolution = VolumeResolution::Mid;

        struct {
            SpinChannel type = SPIN_ALPHA;
            int idx = 0;
        } mo;

        struct {
            ElectronicStructureNtoComponent type = ElectronicStructureNtoComponent::Particle;
            int lambda_idx = 0;
            int idx = 1;
        } nto;

        bool use_obb = true;
        bool show_window = false;
    md_system_state_t export_state = {0};
    } export_state;

    // Arena for persistent allocations for this module, reset by reset_data() when the system goes
    // away. Nothing the attribute table owns is allocated from it - those pointers belong to the
    // system - so what lives here is this component's own derived state.
    md_allocator_i* arena = 0;
    ApplicationState* app = nullptr;    // for what has to be released with the GPU stream (the panels' meshes)

    // The GTO basis, built from the system's basis/ attributes. One shell per contracted radial
    // shell,
    // pure radial normalization coefficients (sph→cart factors baked in by gto.c).
    // Allocated from arena; zeroed when reset_data() is called.
    md_gto_basis_t basis = {};

    // The {M,A} shape of the coefficient attribute. Cached at load beside the QM atom domain,
    // because the gate below asks for it every frame and the answer only changes with the system.
    size_t num_molecular_orbitals() const {
        return num_mos;
    }

    size_t num_natural_transition_orbitals() const {
        return num_excited;
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t event_idx = 0; event_idx < num_events; ++event_idx) {
            const viamd::Event& e = events[event_idx];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;

                arena = md_arena_allocator_create(state.allocator.persistent, MD_MEGABYTES(4));
                app = &state;
                md_bitmask_init(&critical_points.selection_mask, 0, MD_BITSET_MAX_INDEX, state.allocator.persistent);
                md_bitmask_init(&critical_points.highlight_mask, 0, MD_BITSET_MAX_INDEX, state.allocator.persistent);

                // Reservation only, pages are committed on demand. The high-water mark is a single
                // RIXS map (num_photon_energies * num_grid_samples doubles, the latter capped at
                // 4096) plus the XAS panel, so this is generous by a wide margin.
                rixs.alloc = md_vm_arena_create(MD_MEGABYTES(512));

                int gl_major, gl_minor;
                glGetIntegerv(GL_MAJOR_VERSION, &gl_major);
                glGetIntegerv(GL_MINOR_VERSION, &gl_minor);
                use_gpu_path = (!FORCE_CPU_PATH && gl_major >= 4 && gl_minor >= 3);
#if MD_ENABLE_GPU
                md_gto_gpu_initialize(state.gpu_device);
                md_topo_gpu_initialize(state.gpu_device);
#endif
                break;
            }
            case viamd::EventType_ViamdShutdown:
                if (app) {
                    panel_meshes_free(&orb.mesh, app);
                    panel_meshes_free(&nto.mesh, app);
                }
#if MD_ENABLE_GPU
                gl::pbo_upload_shutdown();
                md_gto_gpu_shutdown();
                md_topo_gpu_shutdown();
                // Nothing is freed here: the device scratch and the dataset's uploaded basis are
                // both the application's, and it destroys them.
#endif
                gto_critical_point_job_destroy(critical_points.gto_job);
                critical_points.gto_job = nullptr;
                md_arena_allocator_destroy(arena);
                md_bitmask_free(&critical_points.selection_mask);
                md_bitmask_free(&critical_points.highlight_mask);
                if (rixs.alloc) {
                    md_vm_arena_destroy(rixs.alloc);
                    rixs.alloc = nullptr;
                }
                break;
            case viamd::EventType_ViamdFrameTick: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;

                // Here rather than in the window that started it, so a result lands whether or not
                // that window is open.
                critical_points_poll_gto();

                // What the SYSTEM holds, not what this component parsed. Orbitals are what every
                // window here is ultimately about, and the transition windows want their own data on
                // top of that - so a system published by any QM reader opens the ones it can fill
                // and none of the ones it cannot.
                if (num_qm_atoms > 0 && num_molecular_orbitals() > 0) {
                    // Before any window draws, so both transition windows see the same numbers and
                    // neither has to be open for the other to have data.
                    update_nto_derived_data(state.mold.sys);

                    update_critical_points(state);

                    draw_orb_window(state);
                    draw_summary_window(state);
                    draw_critical_points_window(state);
                    draw_rsp_window(state);
                    if (num_excited > 0) {
                        draw_nto_window(state);
                        draw_flow_window(state);
                    }
                    draw_export_window(state);
                }
                break;
            }
            case viamd::EventType_ViamdWindowDrawMenu:
                if (num_qm_atoms > 0 && num_molecular_orbitals() > 0) {
                    if (ImGui::BeginMenu("Quantum Chemistry")) {
                        ImGui::Checkbox("Summary", &summary.show_window);
                        ImGui::Checkbox("Critical Points", &critical_points.show_window);
                        ImGui::Checkbox("Response", &rsp.show_window);
                        ImGui::Checkbox("Orbital Grid", &orb.show_window);
                        if (num_excited > 0) {
                            ImGui::Checkbox("Transition Analysis", &nto.show_window);
                            ImGui::Checkbox("Charge Transfer Analysis", &flow.show_window);
                        }
                        ImGui::Checkbox("Export", &export_state.show_window);
                        ImGui::EndMenu();
                    }
                }
                break;
            case viamd::EventType_ViamdRenderTransparent: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
				const ApplicationState& state = *(ApplicationState*)e.payload;

                md_temp_scope_t temp = md_temp_begin_in(state.allocator.frame);
                MD_DEFER { md_temp_end(temp); };

                if (critical_points.enabled && critical_points.simp_graph.num_vertices > 0) {
                    draw_critical_points(state, temp);
                }

                //ApplicationState& state = *(ApplicationState*)e.payload;
                //draw_orb_volume(state);
                break;
            }
            case viamd::EventType_ViamdRenderKey: {
                // What the critical points are drawn from (RenderTransparent): they are selected and highlighted
                // by hovering and clicking them
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_RenderKey);
                RenderKey* rk = (RenderKey*)e.payload;
                if (critical_points.enabled && critical_points.simp_graph.num_vertices > 0) {
                    const md_topo_extremum_graph_t& g = critical_points.simp_graph;
                    const uint64_t counts[2] = {(uint64_t)g.num_vertices, (uint64_t)g.num_edges};
                    uint64_t h = md_hash64(counts, sizeof(counts), rk->key);
                    h = md_hash64(&g.vertices, sizeof(g.vertices), h);
                    h = md_bitset_hash64(md_bitmask_view(&critical_points.selection_mask), h);
                    h = md_bitset_hash64(md_bitmask_view(&critical_points.highlight_mask), h);
                    h = md_hash64(&rk->state->selection.color, sizeof(rk->state->selection.color), h);
                    h = md_hash64(&critical_points.point_scale, sizeof(critical_points.point_scale), h);
                    h = md_hash64(&critical_points.show_paths, sizeof(critical_points.show_paths), h);
                    h = md_hash64(&rk->state->visuals.background.color, sizeof(rk->state->visuals.background.color), h);
                    rk->key = h;
                }
                break;
            }
            case viamd::EventType_ViamdInteractionSurface: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_InteractionSurfaceEvent);
                MD_ASSERT(e.payload);
                InteractionSurfaceEvent* event = (InteractionSurfaceEvent*)e.payload;
                if (event->surface_id == interaction_surface_main) {
                    if (critical_points.enabled && arena) {
                        switch (event->kind) {
                        case InteractionSurfaceEventKind::Hover: [[fallthrough]];
                        case InteractionSurfaceEventKind::Click: {
                            md_bitmask_clear(&critical_points.highlight_mask);
                            if (event->hit.domain == PickingDomain_CriticalPoints) {
                                md_bitmask_or(&critical_points.highlight_mask, md_bitset_range((uint32_t)(event->hit.local_idx), (uint32_t)(event->hit.local_idx) + 1));
                            }

                            if (event->kind == InteractionSurfaceEventKind::Click && event->selection_mode == InteractionSelectionMode::None && event->hit.domain == 0) {
                                // A plain click on empty space clears the selection
                                md_bitmask_clear(&critical_points.selection_mask);
                            } else if (event->selection_mode == InteractionSelectionMode::Append) {
                                md_bitmask_or(&critical_points.selection_mask, md_bitmask_view(&critical_points.highlight_mask));
                            } else if (event->selection_mode == InteractionSelectionMode::Remove) {
                                if (event->hit.domain == PickingDomain_CriticalPoints) {
                                    md_bitmask_andnot(&critical_points.selection_mask, md_bitmask_view(&critical_points.highlight_mask));
                                } else if (event->hit.domain == 0) {
                                    md_bitmask_clear(&critical_points.selection_mask);
                                }
                            }
                            break;
                        }
                        case InteractionSurfaceEventKind::RegionSelect: {
                            // The selection is read (Remove) while the highlight is written: two masks
                            md_bitmask_t* region = &critical_points.highlight_mask;
                            md_bitmask_clear(region);
                            if (event->selection_mode == InteractionSelectionMode::Remove) {
                                // When removing, only consider currently selected items as candidates for region selection
                                md_bitset_iter_t it = md_bitmask_iter(&critical_points.selection_mask);
                                while (md_bitset_iter_next(&it)) {
                                    size_t idx = it.idx;
                                    vec3_t xyz = { critical_points.simp_graph.vertices[idx].x,
                                                   critical_points.simp_graph.vertices[idx].y,
                                                   critical_points.simp_graph.vertices[idx].z };
                                    xyz *= BOHR_TO_ANGSTROM;
                                    vec2_t coord = world_to_surface_project(xyz, event->world_to_clip, event->surface_size);
                                    if (coord.x >= event->region_min.x && coord.x <= event->region_max.x &&
                                        coord.y >= event->region_min.y && coord.y <= event->region_max.y)
                                    {
                                        md_bitmask_set_bit(region, (uint32_t)idx);
                                    }
                                }
                            } else {
                                for (size_t i = 0; i < critical_points.simp_graph.num_vertices; ++i) {
                                    vec3_t xyz = { critical_points.simp_graph.vertices[i].x,
                                                   critical_points.simp_graph.vertices[i].y,
                                                   critical_points.simp_graph.vertices[i].z };
                                    xyz *= BOHR_TO_ANGSTROM;
                                    vec2_t coord = world_to_surface_project(xyz, event->world_to_clip, event->surface_size);
                                    if (coord.x >= event->region_min.x && coord.x <= event->region_max.x &&
                                        coord.y >= event->region_min.y && coord.y <= event->region_max.y)
                                    {
                                        md_bitmask_set_bit(region, (uint32_t)i);
                                    }
                                }
                            }

                            if (event->region_phase == InteractionSurfaceEventPhase::Commit) {
                                // Merge highlight into selection
                                if (event->selection_mode == InteractionSelectionMode::Append) {
                                    md_bitmask_or(&critical_points.selection_mask, md_bitmask_view(&critical_points.highlight_mask));
                                }
                                else if (event->selection_mode == InteractionSelectionMode::Remove) {
                                    md_bitmask_andnot(&critical_points.selection_mask, md_bitmask_view(&critical_points.highlight_mask));
                                }
                                md_bitmask_clear(&critical_points.highlight_mask);
                            }
                            break;
                        }
                        default:
                            break;
                        }
                    }
                }
                break;
            }
            case viamd::EventType_ViamdPickingTooltipTextRequest: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_PickingTooltipTextRequest);
                MD_ASSERT(e.payload);
                PickingTooltipTextRequest* req = (PickingTooltipTextRequest*)e.payload;
                if (req->hit.domain == PickingDomain_CriticalPoints) {
                    uint32_t cp_idx = req->hit.local_idx;
                    if (cp_idx < critical_points.simp_graph.num_vertices) {
                        const md_topo_critical_point_type_t type = md_topo_vertex_type(&critical_points.simp_graph, cp_idx);
                        const md_topo_vert_t& v = critical_points.simp_graph.vertices[cp_idx];
                        // The title in lower case, as the tooltips have them ("bond critical point")
                        const CriticalPointStyle& style = critical_point_style(type);
                        char type_str[32];
                        const int type_len = snprintf(type_str, sizeof(type_str), "%s", style.name);
                        md_convert_to_lower(type_str, (size_t)MD_CLAMP(type_len, 0, (int)sizeof(type_str) - 1));
                        tooltip_title(req, "%s", type_str);
                        tooltip_row(req, MD_STR_LIT("type"), "%s %s" TOOLTIP_SEP "%s", style.label, style.signature, style.topology);
                        tooltip_row(req, MD_STR_LIT("density"), "%.6g au", v.value);
                        // The graph is in bohr
                        tooltip_row(req, MD_STR_LIT("position"), "%.3f, %.3f, %.3f " TOOLTIP_ANGSTROM,
                                    v.x * BOHR_TO_ANGSTROM, v.y * BOHR_TO_ANGSTROM, v.z * BOHR_TO_ANGSTROM);
                    }
                }
                break;
            }
            case viamd::EventType_ViamdLoadData: {
                 MD_ASSERT(e.payload_type == viamd::EventPayloadType_LoadData);
                 LoadDataPayload& payload = *(LoadDataPayload*)e.payload;
                 ApplicationState& state = *payload.app_state;

                 // Nothing is parsed here. The loader published this system's table before the
                 // event was broadcast, and there is no branch on WHICH loader did: a file type
                 // check would gate these windows on VeloxChem having been the reader, which is
                 // exactly the coupling this component no longer has.
                 //
                 // Keyed on the coefficient attribute's VERSION, not on merely being present: this
                 // event also fires for a trajectory or a supplemental file loaded on top of a
                 // system that already has orbitals, and re-running the setup then would throw away
                 // the user's cameras and selections. A version is bumped only when the attribute
                 // is actually written, so a reload re-initialises and nothing else does.
                 const md_attribute_t* coeff = md_attributes_find(&state.mold.sys.attributes, es_path::alpha_coefficient);
                 const uint64_t version = coeff ? md_attributes_version(&state.mold.sys.attributes, coeff->id) : 0;
                 if (version != 0 && version != qm_version) {
                     qm_version = version;
                     init_qm_from_system(state);
                 }
                 break;
            }
            case viamd::EventType_ViamdSystemFree: {
                reset_data();
                break;
            }
            case viamd::EventType_ViamdSerialize: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_SerializationState);
                serialize_workspace(*(viamd::serialization_state_t*)e.payload);
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_DeserializationState);
                viamd::deserialization_state_t& deser = *(viamd::deserialization_state_t*)e.payload;
                if (md_str_eq(viamd::section_header(deser), MD_STR_LIT("QuantumChemistry"))) {
                    deserialize_workspace(deser);
                }
                break;
            }
            case viamd::EventType_ViamdSystemStateChanged: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState& state = *(ApplicationState*)e.payload;

                // Nothing below reads the reader any more - the QM atom domain and the dipole both
                // come out of the table - so the gate is whether there are QM atoms to re-gather.
                // atom_xyzw was sized by init_qm_from_system against the domain as it stood then.
                // A supplemental load can republish qm/atom/* at a different length without
                // touching the coefficient attribute that gates re-init, so the two lengths
                // agreeing is something to CHECK here rather than an invariant to assert. A
                // disagreement skips the update; the load event that caused it re-initialises.
                QmAtoms qm = {};
                if (num_qm_atoms > 0 && es_qm_atoms(&qm, state.mold.sys) &&
                    md_array_size(atom_xyzw) == qm.count) {
                    // Update atom_xyz
					// Calculate nuclei charge weighted center of charge for later use in orbital centering
					dvec3_t nucl_dipole = {0};

					// One entry per QM atom, gathered from wherever that atom sits in the system.
					// No branch on whether a map exists: without one the QM atoms ARE the system's,
					// so the count matches and qm_to_system_atom is the identity.
					const size_t count = qm.count;
                    for (size_t i = 0; i < count; ++i) {
                        const size_t idx = qm_to_system_atom(qm, i);
                        atom_xyzw[i].x = (float)(state.mold.state.xyz[idx].x * ANGSTROM_TO_BOHR);
                        atom_xyzw[i].y = (float)(state.mold.state.xyz[idx].y * ANGSTROM_TO_BOHR);
                        atom_xyzw[i].z = (float)(state.mold.state.xyz[idx].z * ANGSTROM_TO_BOHR);

						int z = md_atom_atomic_number(&state.mold.sys.atom, idx);
						nucl_dipole.x += atom_xyzw[i].x * z;
						nucl_dipole.y += atom_xyzw[i].y * z;
						nucl_dipole.z += atom_xyzw[i].z * z;
                    }

                    // The occupations summed, rather than a count the reader kept beside them.
                    const double num_electrons = es_electron_count(state.mold.sys, es_path::alpha_occupation) +
                                                 es_electron_count(state.mold.sys, es_path::beta_occupation);

                    // The centroid of the ELECTRONIC charge: mu = sum_A Z_A R_A - <r_e> N_e, so
                    // <r_e> = (sum_A Z_A R_A - mu) / N_e. This used to divide the nuclear term by
                    // N_e and only then subtract the dipole, which subtracts a dipole from a
                    // position - dimensionally impossible, and disagreeing with the load path,
                    // which had it right. The two now compute one formula.
                    const double inv_ne = num_electrons > 0.0 ? 1.0 / num_electrons : 1.0;
                    center_of_charge = (nucl_dipole - dipole_ground_state(state.mold.sys)) * inv_ne;
                    
                    // Recalculate OABB and AABB
                    oabb.orientation = mat3_PCA(atom_xyzw, md_array_size(atom_xyzw));
                    calculate_bounds(oabb.min_ext.elem, oabb.max_ext.elem, atom_xyzw, md_array_size(atom_xyzw), oabb.orientation);
                    calculate_bounds(aabb.min_ext.elem, aabb.max_ext.elem, atom_xyzw, md_array_size(atom_xyzw));
                }

                for (size_t i = 0; i < md_array_size(state.representation.reps); ++i) {
                    if (state.representation.reps[i].type == RepresentationType::ElectronicStructure) {
                        flag_representation_as_dirty(&state.representation.reps[i]);
                    }
                }
                
                break;
            }
            case viamd::EventType_ViamdPickingRangeReserve: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_PickingSpace);
                PickingSpace* space = (PickingSpace*)e.payload;
                if (critical_points.enabled) {
                    size_t num_cp = md_topo_num_critical_points(&critical_points.simp_graph);
                    picking_range_reserve(&critical_points.picking_range, space, PickingDomain_CriticalPoints, num_cp);
                }
                break;
            }
            case viamd::EventType_ViamdTaskProgressRequest: {
                // The analytic critical point search, for the task overlay: once it has run for a moment, with the
                // flag its cancel button raises (collected as usual: what was found so far is shown, marked)
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_TaskProgressRequest);
                TaskProgressRequest* req = (TaskProgressRequest*)e.payload;
                GtoCriticalPointJob* job = critical_points.gto_job;
                if (job && !job->cancel && !job->dropped && md_tick_to_seconds(md_tick_now() - job->start) >= 0.25) {
                    task_progress_add(req, "Critical points", job->progress, &job->cancel);
                }
                break;
            }
            case viamd::EventType_ViamdViewFit: {
                MD_ASSERT(e.payload_type == viamd::EventPayloadType_ViewFitRequest);
                ViewFitRequest* req = (ViewFitRequest*)e.payload;
                if (req && req->surface_id == interaction_surface_main) {
                    const md_bitmask_t* bf = nullptr;
                    switch (req->round) {
                    case ViewFitRound_Highlight:
                        bf = &critical_points.highlight_mask; break;
                    case ViewFitRound_Selection:
                        bf = &critical_points.selection_mask; break;
					case ViewFitRound_Visible: [[fallthrough]];
                    default:
                        break;
                    }

                    // Add highlighted atoms to the mask used for view fitting, so that they are included in the computed optimal view
                    if (critical_points.enabled) {
                        size_t popcount = bf ? md_bitmask_count(bf) : 0;
                        if (popcount > 0) {
                            vec4_t* dst_xyzw = md_array_extend(req->xyzw, popcount, req->alloc);
                            if (dst_xyzw) {
                                md_bitset_iter_t it = md_bitmask_iter(bf);
                                while (md_bitset_iter_next(&it)) {
                                    size_t idx = it.idx;
                                    dst_xyzw->x = (float)(critical_points.simp_graph.vertices[idx].x * BOHR_TO_ANGSTROM);
                                    dst_xyzw->y = (float)(critical_points.simp_graph.vertices[idx].y * BOHR_TO_ANGSTROM);
                                    dst_xyzw->z = (float)(critical_points.simp_graph.vertices[idx].z * BOHR_TO_ANGSTROM);
                                    dst_xyzw->w = 1.0f;
                                    dst_xyzw++;
                                }
                            }
                        }
                    }
                }
                break;
            }
            default:
                break;
            }
        }
    }

    // The critical points over the image, in the overlay: neither hidden by the atoms around them nor hiding them
    // except where a marker is drawn, and picked wherever it is seen. Back to front: the paths as thin strips,
    // then a disc per point facing the camera, in its type's colour (critical_point_styles) inside an outline
    // that stands out from the background. The highlighted point is drawn larger, inside a ring of the
    // highlight's colour, a selected one inside a ring of the selection's; while any are selected, the others
    // keep only part of their colour, as atoms do. Sizes are in Å, so they scale with the view as the atoms do
    // (and an orthographic view does not make them huge, as points of a size by depth did).
    void draw_critical_points(const ApplicationState& state, md_temp_scope_t temp) {
        const md_topo_extremum_graph_t& g = critical_points.simp_graph;
        const size_t num_cp = g.num_vertices;

        immediate::Scope scope(state.gfx.overlay, "qm_critical_points");
        immediate::set_picking_base_idx(scope, critical_points.picking_range.beg);

        const mat4_t view = camera_world_to_view_matrix(state.view.camera);
        const vec3_t right = {view[0][0], view[1][0], view[2][0]};
        const vec3_t up    = {view[0][1], view[1][1], view[2][1]};
        const vec3_t back  = {view[0][2], view[1][2], view[2][2]};    // towards the camera

        const vec3_t bg = state.visuals.background.color;
        const bool dark_bg = 0.2126f * bg.x + 0.7152f * bg.y + 0.0722f * bg.z < 0.5f;
        const vec4_t outline_color = dark_bg ? vec4_t{0.92f, 0.92f, 0.92f, 0.90f} : vec4_t{0.10f, 0.10f, 0.10f, 0.90f};
        const vec4_t bond_path_color = dark_bg ? vec4_t{0.80f, 0.80f, 0.80f, 0.85f} : vec4_t{0.25f, 0.25f, 0.25f, 0.85f};
        vec4_t rcp_path_color = critical_point_style(MD_TOPO_JOIN_SADDLE).color;
        rcp_path_color.w = 0.75f;
        const float scale = critical_points.point_scale;

        auto pos = [&g](size_t i) {
            return vec3_t{g.vertices[i].x, g.vertices[i].y, g.vertices[i].z} * (float)BOHR_TO_ANGSTROM;
        };

        // The paths: the bond paths, from a bond critical point to the attractors it joins, and the RCP attractor
        // paths, from a ring critical point to the cage critical points beside it (md_topo traces each from the
        // saddle, so a graph's edge runs from it)
        if (critical_points.show_paths && g.num_edges > 0) {
            immediate::Vertex* v = md_temp_alloc_array(temp, immediate::Vertex, 6 * g.num_edges);
            size_t n = 0;
            for (size_t i = 0; i < g.num_edges; ++i) {
                const uint32_t i0 = g.edges[i].from, i1 = g.edges[i].to;
                if (i0 >= num_cp || i1 >= num_cp) continue;
                const bool bond = md_topo_vertex_type(&g, i0) == MD_TOPO_SPLIT_SADDLE;
                const vec3_t p0 = pos(i0), p1 = pos(i1);
                const vec3_t side = vec3_cross(p1 - p0, back);
                const float len = vec3_length(side);
                if (len < 1.0e-6f) continue;    // along the line of sight: nothing to see
                const vec3_t d = side * ((bond ? 0.018f : 0.012f) * scale / len);
                const uint32_t col = convert_color(bond ? bond_path_color : rcp_path_color);
                const vec3_t q[4] = {p0 - d, p0 + d, p1 + d, p1 - d};
                const int tri[6] = {0, 1, 2, 0, 2, 3};
                for (int k = 0; k < 6; ++k) {
                    v[n].coord = q[tri[k]];
                    v[n].color = col;
                    v[n].normal = back;
                    v[n].picking_idx = 0xFFFFFFFF;   // a path is not a point to pick
                    n++;
                }
            }
            immediate::triangles(scope, v, n);
        }

        // The points, farthest first
        uint32_t* order = md_temp_alloc_array(temp, uint32_t, num_cp);
        float* depth = md_temp_alloc_array(temp, float, num_cp);
        for (size_t i = 0; i < num_cp; ++i) {
            order[i] = (uint32_t)i;
            depth[i] = vec3_dot(pos(i), back);
        }
        std::sort(order, order + num_cp, [depth](uint32_t a, uint32_t b) { return depth[a] < depth[b]; });

        constexpr int SEG = 16;
        vec2_t circle[SEG + 1];
        for (int k = 0; k <= SEG; ++k) {
            const float a = (float)(2.0 * MD_PI * k / SEG);
            circle[k] = {cosf(a), sinf(a)};
        }

        const vec4_t selection_color = state.selection.color.selection.outline;
        const vec4_t highlight_color = state.selection.color.highlight.outline;
        const bool any_selected = md_bitmask_count(&critical_points.selection_mask) > 0;
        const float saturation = any_selected ? state.selection.color.saturation : 1.0f;

        // At most four discs a point: the highlight's ring, the selection's, the outline, the fill
        immediate::Vertex* v = md_temp_alloc_array(temp, immediate::Vertex, num_cp * 4 * SEG * 3);
        size_t n = 0;
        auto disc = [&](vec3_t c, float r, vec4_t color, uint32_t picking_idx) {
            const uint32_t col = convert_color(color);
            const vec3_t u = right * r, w = up * r;
            for (int k = 0; k < SEG; ++k) {
                const vec3_t e[3] = {c, c + u * circle[k].x + w * circle[k].y, c + u * circle[k + 1].x + w * circle[k + 1].y};
                for (int j = 0; j < 3; ++j) {
                    v[n].coord = e[j];
                    v[n].color = col;
                    v[n].normal = back;
                    v[n].picking_idx = picking_idx;
                    n++;
                }
            }
        };
        for (size_t j = 0; j < num_cp; ++j) {
            const uint32_t i = order[j];
            const CriticalPointStyle& style = critical_point_style(md_topo_vertex_type(&g, i));
            const bool selected    = md_bitmask_test_bit(&critical_points.selection_mask, i);
            const bool highlighted = md_bitmask_test_bit(&critical_points.highlight_mask, i);
            const float r = style.radius * scale * (highlighted ? 1.3f : 1.0f);
            const vec3_t c = pos(i);
            if (highlighted) disc(c, r * (selected ? 1.75f : 1.4f), highlight_color, i);
            if (selected)    disc(c, r * 1.4f, selection_color, i);
            disc(c, r, outline_color, i);
            const vec4_t fill = (any_selected && !selected) ? tint_color(style.color, style.color, 0.0f, saturation) : style.color;
            disc(c, r * 0.78f, fill, i);
        }
        immediate::triangles(scope, v, n);
    }

    void reset_data() {
        // First: the run holds nothing of the system's, but its result would land in a graph that
        // belongs to the next one.
        gto_critical_point_job_destroy(critical_points.gto_job);
        critical_points.gto_job = nullptr;
#if MD_ENABLE_GPU
        // The topo context is the only GPU object this component owns. The uploaded basis goes with
        // the system it was built from, released by free_system_data.
        md_topo_gpu_context_destroy(critical_points.topo_ctx);
#endif
        //md_gl_mol_destroy(gl_mol);
        md_gl_rep_destroy(gl_rep);
        gl_rep = {};
        // Zeroed with nto below, which would leave its slot taken for good: one for each system loaded
        md_gl_rep_destroy(nto.gl_rep);
        ctv_free();
        num_qm_atoms = 0;
        num_mos = 0;
        num_excited = 0;
        qm_version = 0;
        // The charge transfer jobs read the arena (the strip's AO to unit map) and the Lowdin job: they go first
        ct_state_job_destroy(ct_job);
        ct_job = nullptr;
        strip_job_destroy(strip_job);
        strip_job = nullptr;
        lowdin_job_destroy(lowdin);
        lowdin = nullptr;
        md_arena_allocator_reset(arena);
        md_topo_extremum_graph_free(&critical_points.raw_graph);
        md_topo_extremum_graph_free(&critical_points.simp_graph);
        // The masks are made once, at initialisation, in the persistent allocator: carried across the
        // value-initialisation and emptied. Zeroed, a mask has an empty domain and drops every index
        // set in it, which left hovering, clicking and region-selecting critical points doing nothing
        // once a system had been freed (any load after the first).
        // The window and how the points are drawn are the user's, not the system's: carried across too.
        const md_bitmask_t cp_selection = critical_points.selection_mask;
        const md_bitmask_t cp_highlight = critical_points.highlight_mask;
        const bool  cp_show_window = critical_points.show_window;
        const float cp_point_scale = critical_points.point_scale;
        const bool  cp_show_paths  = critical_points.show_paths;
        critical_points = {};
        critical_points.selection_mask = cp_selection;
        critical_points.highlight_mask = cp_highlight;
        critical_points.show_window    = cp_show_window;
        critical_points.point_scale    = cp_point_scale;
        critical_points.show_paths     = cp_show_paths;
        md_bitmask_clear(&critical_points.selection_mask);
        md_bitmask_clear(&critical_points.highlight_mask);
        critical_points.raw_graph.alloc = arena;
        critical_points.simp_graph.alloc = arena;
        basis = {};  // arena reset above invalidates allocations; zero the struct
        // Same as basis: the arena reset above released these, so the struct is zeroed rather
        // than freed through td_cache_free().
        td_cache = TransitionDensityCache{};
        atom_xyzw = nullptr;
        if (app) {
            panel_meshes_free(&orb.mesh, app);
            panel_meshes_free(&nto.mesh, app);
        }
        orb = QuantumChemistry::Orb{};
        nto = QuantumChemistry::Nto{};
        // The arena reset above already released everything these held.
        flow = QuantumChemistry::Flow{};
        rsp = QuantumChemistry::Rsp{};
        // The rixs arena outlives the VLX object, so carry it across the value-initialization.
        // Rewinding it first drops the cached arrays that the zeroed struct would otherwise orphan.
        md_allocator_i* rixs_alloc = rixs.alloc;
        reset_rixs_cache(rixs);
        rixs = QuantumChemistry::Rixs{};
        rixs.alloc = rixs_alloc;
        vib = QuantumChemistry::Vib{};
        opt = QuantumChemistry::Opt{};
    }

    void update_nto_group_colors(const md_system_t& sys) {
        QmAtoms qm = {};
        if (!es_qm_atoms(&qm, sys)) return;

        md_temp_scope_t temp = md_temp_begin();
		MD_DEFER { md_temp_end(temp); };
        size_t num_atoms = qm.count;
        uint32_t* colors = (uint32_t*)md_temp_alloc_array(temp, uint32_t, num_atoms);
        // Group 0 is "not in any group", so it keeps the element colour the molecule would have
        // anyway. Stamping the unassigned grey over everything would make an ungrouped system look
        // like a system whose groups all happen to be grey.
        const md_element_t* atom_z = qm.atomic_number;
        for (size_t i = 0; i < num_atoms; ++i) {
            const size_t group_idx = nto.atom_group_idx[i];
            uint32_t color = U32_MAGENTA;
            if (group_idx == 0) {
                color = atom_z ? md_util_element_cpk_color(atom_z[i]) : U32_MAGENTA;
            } else if (group_idx < nto.group.count) {
                color = convert_color(nto.group.color[group_idx]);
            }
            colors[i] = color;
        }

        md_gl_rep_set_atom_colors(nto.gl_rep, 0, (uint32_t)num_atoms, colors, 0);
        nto.render_hash = 0;
    }

    // Everything this component needs about a quantum chemistry system, out of the system's own
    // attribute table. No reader appears in it, which is the whole point: it runs for ANY loaded
    // system that carries orbitals, whichever producer filled the table.
    void init_qm_from_system(ApplicationState& state) {
        // The system's attribute table is already populated by whichever loader read the
        // file - md_vlx_system_init_from_file and md_vlx_system_supplement_from_file for
        // an h5, but the loader is not named anywhere below. This component is a consumer
        // of that table and never its producer, which is what lets a system published by
        // some other QM reader drive the same windows.

        // The basis comes out of that table too: an evaluator needs the system, not a
        // reader.
        md_gto_basis_free(&basis, arena);
        if (!md_gto_basis_extract_attributes(&basis, &state.mold.sys.attributes, arena)) {
            MD_LOG_ERROR("Failed to build a GTO basis from the system's attributes");
            return;
        }
        // The QM atom domain, out of the table the load path filled. Standalone or
        // supplemental was decided there, by which entry point the loader called, and
        // the atomic numbers were checked before the flag was ever set
        // (md_vlx_system_is_file_supplemental) - so there is nothing left to validate
        // or to publish here.
        QmAtoms qm = {};
        if (!es_qm_atoms(&qm, state.mold.sys)) {
            MD_LOG_ERROR("The system carries no QM atom domain");
            return;
        }
        // The one number kept: see num_qm_atoms above for why the pointers are not.
        num_qm_atoms = qm.count;
        const size_t num_atoms = qm.count;

        // A new dataset: whatever state index was cached names a different state now.
        td_cache_free();

        num_mos = 0;
        es_orbital_extent(state.mold.sys, &num_mos, nullptr);
        num_excited = es_excited_state_count(state.mold.sys);

        // From the published occupations rather than from the vlx object: the frontier
        // is one scan of a column already in the table, and reading it there is what a
        // system published by some other QM reader will need.
        //
        // It also answers two cases md_vlx's own indices got wrong. A channel holding
        // no electrons at all - beta in a fully polarised open shell - reported HOMO 0
        // and labelled the first orbital; it now reports -1 and labels nothing. A fully
        // occupied set reported 0 for both; LUMO is now num_orbitals, out of range,
        // because there is no LUMO to point at.
        OrbitalFrontier frontier = {};
        es_orbital_frontier(&frontier, state.mold.sys, es_path::alpha_occupation);
        homo_idx[0] = frontier.homo_idx;
        lumo_idx[0] = frontier.lumo_idx;

        es_orbital_frontier(&frontier, state.mold.sys, es_path::beta_occupation);
        homo_idx[1] = frontier.homo_idx;
        lumo_idx[1] = frontier.lumo_idx;

        size_t num_colors = state.mold.sys.atom.count;
        uint32_t* colors = (uint32_t*)md_vm_arena_push(state.allocator.frame, num_colors * sizeof(uint32_t));
        color_atoms_type(colors, num_colors, state.mold.sys);

        if (qm.system_index) {
            md_bitset_t mask = md_bitset_from_indices((uint32_t*)qm.system_index, qm.count, state.allocator.frame);
            filter_colors(colors, num_colors, mask);
            reset_view(&default_view, state.mold.state, &mask);
        } else {
            reset_view(&default_view, state.mold.state);
        }

        gl_rep = md_gl_rep_create(state.mold.gl_mol);
        md_gl_rep_set_atom_colors(gl_rep, 0, (uint32_t)num_colors, colors, 0);
        orb.render_hash = 0;

					const md_element_t* atom_z = qm.atomic_number;
					dvec3_t nucl_dipole = { 0, 0, 0 };

        md_array_resize(atom_xyzw, num_atoms, arena);
        const dvec3_t* atom_vlx = qm.coordinate;
        for (size_t i = 0; i < num_atoms; ++i) {
            dvec3_t xyz = atom_vlx[i] * ANGSTROM_TO_BOHR;
            atom_xyzw[i] = vec4_set((float)xyz.x, (float)xyz.y, (float)xyz.z, 1.0f);
						nucl_dipole += xyz * atom_z[i];
        }

					const dvec3_t ground_state_dipole = dipole_ground_state(state.mold.sys);

        const double num_electrons = es_electron_count(state.mold.sys, es_path::alpha_occupation) +
                                     es_electron_count(state.mold.sys, es_path::beta_occupation);
        const double inv_ne = num_electrons > 0.0 ? 1.0 / num_electrons : 1.0;
        center_of_charge = (nucl_dipole - ground_state_dipole) * inv_ne;

        oabb.orientation = mat3_PCA(atom_xyzw, md_array_size(atom_xyzw));
        calculate_bounds(oabb.min_ext.elem, oabb.max_ext.elem, atom_xyzw, md_array_size(atom_xyzw), oabb.orientation);
        calculate_bounds(aabb.min_ext.elem, aabb.max_ext.elem, atom_xyzw, md_array_size(atom_xyzw));

        // NTO

        // NTOs specifically, not excited states in general: everything set up below is
        // the transition analysis, which has nothing to show without them.
        if (num_excited > 0) {
            //nto.show_window = true;

            picking_surface_init(&nto.picking_surface, interaction_surface_nto);
						vec3_t center = (oabb.min_ext + oabb.max_ext) * 0.5f * BOHR_TO_ANGSTROM;
						vec3_t half_ext = (oabb.max_ext - oabb.min_ext) * 0.5f * BOHR_TO_ANGSTROM;
            nto.target = compute_optimal_view(center, half_ext, oabb.orientation, nto.distance_scale);
            nto.camera = nto.target;

						nto.atom_group_idx = md_array_create(uint32_t, num_atoms, arena);
            md_memset(nto.atom_group_idx, 0, sizeof(uint32_t) * num_atoms);

            snprintf(nto.group.label[0], sizeof(nto.group.label[0]), "Unassigned");
            nto.group.color[0] = vec4_t{ 0.25f, 0.25f, 0.25f, 1.0f };

            for (int i = 1; i < (int)MD_COUNTOF(nto.group.color); ++i) {
                ImVec4 color = ImPlot::GetColormapColor(i - 1, ImPlotColormap_Deep);
                nto.group.color[i] = vec_cast(color);
                snprintf(nto.group.label[i], sizeof(nto.group.label[i]), "Group %i", i);
            }

            // No groups until the user makes one. Slot 0 is the "everything else"
            // bucket every atom starts in; it is not a group the user chose and is
            // never drawn as one. Splitting the molecule in half, or stamping a
            // hardcoded mapping onto a file with a particular name, put invented
            // structure in front of the user and then invited them to read charge
            // transfer off it.
            nto.group.count = 1;
            nto.gl_rep = md_gl_rep_create(state.mold.gl_mol);
            update_nto_group_colors(state.mold.sys);

            // Callculate ballpark scaling factor for dipole vectors
            vec3_t extent = aabb.max_ext - aabb.min_ext;
            float max_ext = MD_MAX(extent.x, MD_MAX(extent.y, extent.z));
            float max_len = 0;
            {
                md_temp_scope_t temp = md_temp_begin();
                MD_DEFER { md_temp_end(temp); };
                max_len = (float)MD_MAX(dipole_max_length(temp, state.mold.sys, MD_STR_LIT("dipole/electric_transition/vector")),
                                     dipole_max_length(temp, state.mold.sys, MD_STR_LIT("dipole/magnetic_transition/vector")));
            }
            if (max_len > 0.0f) {
                nto.dipole.vector_scale = MD_CLAMP((max_ext * 0.75f) / max_len, 0.1f, 10.0f);
            }

#if USE_AABB_GRID
            init_grid(&nto.grid, mat3_ident(), aabb.min_ext, aabb.max_ext, DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM);
#else
            init_grid(&nto.grid, oabb.orientation, oabb.min_ext, oabb.max_ext, DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM);
#endif
            nto.target = default_view;
            nto.camera = nto.target;

            // The charge transfer window's side views: the same view to start from, the molecule in its element
            // colours, a grid at half the density of the main one. What a previous system left goes first.
            ctv_free();
            picking_surface_init(&ctv.picking_surface, interaction_surface_flow);
            ctv.target = default_view;
            ctv.camera = ctv.target;
            ctv.gl_rep = md_gl_rep_create(state.mold.gl_mol);
            {
                QmAtoms qm = {};
                if (es_qm_atoms(&qm, state.mold.sys) && qm.count > 0) {
                    md_temp_scope_t temp = md_temp_begin();
                    MD_DEFER { md_temp_end(temp); };
                    uint32_t* colors = md_temp_alloc_array(temp, uint32_t, qm.count);
                    for (size_t i = 0; i < qm.count; ++i) colors[i] = qm.atomic_number ? md_util_element_cpk_color(qm.atomic_number[i]) : U32_MAGENTA;
                    md_gl_rep_set_atom_colors(ctv.gl_rep, 0, (uint32_t)qm.count, colors, 0);
                }
            }
#if USE_AABB_GRID
            init_grid(&ctv.grid, mat3_ident(), aabb.min_ext, aabb.max_ext, 0.5 * DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM);
#else
            init_grid(&ctv.grid, oabb.orientation, oabb.min_ext, oabb.max_ext, 0.5 * DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM);
#endif
            ctv.initialized = true;
        }

        // RSP. Oscillator strengths are what a linear response calculation produces and
        // what the absorption plot draws, so their presence is the question - a reader
        // with no response type of its own answers it just by publishing them.
        if (md_attributes_find(&state.mold.sys.attributes, MD_STR_LIT("vlx/rsp/oscillator_strength"))) {
            rsp.hovered = -1;
            rsp.selected = 0;
        }

        // OPT. Opens on the last step, which is the converged geometry.
        size_t num_opt_steps = 0;
        if (attribute_series_f64(&num_opt_steps, state.mold.sys, MD_STR_LIT("vlx/opt/energy")) && num_opt_steps > 0) {
            opt.selected = (int)(num_opt_steps - 1);
        }

        // ORB
        //orb.show_window = true;
        picking_surface_init(&orb.picking_surface, interaction_surface_orb);
        orb.target = default_view;
        orb.camera = orb.target;
        // homo_idx is -1 when the file carries no occupations, which is not an orbital to open on.
        const int initial_mo = (homo_idx[0] >= 0 && (size_t)homo_idx[0] < num_mos) ? homo_idx[0] : 0;
        orb.mo_idx = initial_mo;
        orb.scroll_to_idx = initial_mo;

        // Export
        export_state.mo.idx = initial_mo;

        // CriPoAl
        md_bitmask_clear(&critical_points.selection_mask);
        md_bitmask_clear(&critical_points.highlight_mask);

        // Last, so it overrides the placeholder grouping that init sets up.
        apply_pending_groups(state.mold.sys);

        MD_LOG_INFO("Successfully initialized quantum chemistry data");
    }

    // Every one of these forwards to the application's evaluator. What used to live here - packing
    // coefficients for the device, launching the kernel, queueing the readback, and the GL fallback
    // beside it - was never specific to this component; it is the same work for any system that
    // publishes a basis, and it now lives once in viamd.cpp where the GPU resources do.
    //
    // What is left is the translation from this component's vocabulary to an attribute path and a
    // slice, which is the only part that was ever ours.

    bool evaluate_nto(ApplicationState& state, uint32_t vol_tex, const md_grid_t& grid, size_t nto_idx, size_t lambda_idx, ElectronicStructureNtoComponent type, md_gto_op_t op) {
        const md_str_t path = (type == ElectronicStructureNtoComponent::Hole) ? es_path::nto_hole : es_path::nto_particle;
        const md_attribute_slice_t slice = md_attribute_slice_2((uint32_t)nto_idx, (uint32_t)lambda_idx);
        return orbital_evaluate(&state, vol_tex, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, op, DEFAULT_GTO_CUTOFF_VALUE);
    }

    bool evaluate_mo(ApplicationState& state, uint32_t vol_tex, const md_grid_t& grid, SpinChannel mo_type, size_t mo_idx, md_gto_op_t op) {
        const md_str_t path = (mo_type == SPIN_BETA) ? es_path::beta_coefficient : es_path::alpha_coefficient;
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)mo_idx);
        return orbital_evaluate(&state, vol_tex, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, op, DEFAULT_GTO_CUTOFF_VALUE);
    }

    // Releases what the cache holds. Only valid while the arena it was allocated from is alive -
    // reset_data() rewinds that arena and zeroes the struct instead of calling this.
    void td_cache_free() {
        if (td_cache.attach) {
            md_free(arena, td_cache.attach, sizeof(double) * td_cache.dim * td_cache.dim);
        }
        if (td_cache.detach) {
            md_free(arena, td_cache.detach, sizeof(double) * td_cache.dim * td_cache.dim);
        }
        td_cache = TransitionDensityCache{};
    }

    // Both matrices for one excited state, rebuilding them only when the state asked for is not
    // the one already held. False when the system carries no transition densities, or when the two
    // do not agree on their dimension - the callers then draw and attribute nothing, which is what
    // they did when the extract itself failed.
    bool td_cache_ensure(const md_system_t& sys, int state_idx) {
        if (state_idx < 0) {
            return false;
        }
        if (td_cache.state_idx == state_idx && td_cache.attach && td_cache.detach) {
            return true;
        }

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)state_idx);

        size_t a_dim = 0;
        size_t d_dim = 0;
        const double* attach = density_matrix_extract(&a_dim, temp, sys, es_path::attachment_density, &slice);
        const double* detach = density_matrix_extract(&d_dim, temp, sys, es_path::detachment_density, &slice);
        if (!attach || !detach || a_dim == 0 || a_dim != d_dim) {
            // Whatever is held is for a state that is no longer the one being asked about, so it
            // is dropped rather than left to answer the next call.
            td_cache_free();
            return false;
        }

        if (td_cache.dim != a_dim || !td_cache.attach || !td_cache.detach) {
            td_cache_free();
            const size_t bytes = sizeof(double) * a_dim * a_dim;
            td_cache.attach = (double*)md_alloc(arena, bytes);
            td_cache.detach = (double*)md_alloc(arena, bytes);
            if (!td_cache.attach || !td_cache.detach) {
                td_cache_free();
                return false;
            }
            td_cache.dim = a_dim;
        }

        const size_t bytes = sizeof(double) * a_dim * a_dim;
        md_memcpy(td_cache.attach, attach, bytes);
        md_memcpy(td_cache.detach, detach, bytes);
        td_cache.state_idx = state_idx;
        return true;
    }

    bool evaluate_transition_density(ApplicationState& state, uint32_t vol_tex, const md_grid_t& grid, size_t state_idx, ElectronicStructureTransitionDensityComponent component) {
        // Out of the cache rather than out of the attribute: see TransitionDensityCache. The
        // difference is the one component the table would have reconstructed a third matrix for,
        // and it is the subtraction the provider performs anyway.
        if (!td_cache_ensure(state.mold.sys, (int)state_idx)) {
            return false;
        }

        const size_t dim = td_cache.dim;
        const double* density_matrix = td_cache.attach;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        switch (component) {
        case ElectronicStructureTransitionDensityComponent::Detachment:
            density_matrix = td_cache.detach;
            break;
        case ElectronicStructureTransitionDensityComponent::Difference: {
            double* diff = md_temp_alloc_array(temp, double, dim * dim);
            for (size_t i = 0; i < dim * dim; ++i) {
                diff[i] = td_cache.attach[i] - td_cache.detach[i];
            }
            density_matrix = diff;
            break;
        }
        case ElectronicStructureTransitionDensityComponent::Attachment:
        default:
            break;
        }

        return density_matrix_evaluate(&state, vol_tex, grid, density_matrix, dim, MD_GTO_OP_SET);
    }

    bool evaluate_electron_density(ApplicationState& state, uint32_t vol_tex, const md_grid_t& grid, ElectronicStructureSpin spin, md_gto_op_t op) {
        return density_evaluate(&state, vol_tex, grid, es_electron_density_path(spin), nullptr, op);
    }

    static inline ImVec4 make_highlight_color(const ImVec4& color, float factor = 0.2f) {
        // Ensure the factor is not too high, to avoid over-brightening
        factor = MD_CLAMP(factor, 0.0f, 1.0f);

        // Increase the brightness of each component by the factor
        float r = color.x + factor;
        float g = color.y + factor;
        float b = color.z + factor;

        // Clamp the values to the range [0, 1]
        r = ImClamp(r, 0.0f, 1.0f);
        g = ImClamp(g, 0.0f, 1.0f);
        b = ImClamp(b, 0.0f, 1.0f);

        // Return the new highlighted color with the same alpha
        return ImVec4(r, g, b, color.w);
    }

    static inline vec2_t vec2_from_imvec2(ImVec2 v) {
        return { v.x, v.y };
    }

    static inline md_vg_color_t vg_color_from_imvec4(const ImVec4& color) {
        return md_vg_color_from_vec4({ color.x, color.y, color.z, color.w });
    }

    static inline ImVec2 aligned_text_pos(const char* text, ImVec2 pos, ImVec2 alignment = { 0,0 }) {
        return pos - ImGui::CalcTextSize(text) * alignment;
    }

    static inline void compute_vertical_sankey_bezier(ImVec2* p1, ImVec2* p2, ImVec2* p3, ImVec2* p4, ImVec2 source_pos, ImVec2 dest_pos, float thickness) {
        MD_ASSERT(p1);
        MD_ASSERT(p2);
        MD_ASSERT(p3);
        MD_ASSERT(p4);

        *p1 = ImVec2(source_pos.x + thickness * 0.5f, source_pos.y);
        *p4 = ImVec2(dest_pos.x + thickness * 0.5f, dest_pos.y);

        const float dist = fabsf(p1->y - p4->y);
        const float curve_offset = dist * 0.25f;

        *p2 = ImVec2(source_pos.x + thickness * 0.5f, source_pos.y - curve_offset);
        *p3 = ImVec2(dest_pos.x + thickness * 0.5f, dest_pos.y + curve_offset);
    }

    // Draws a vertical flow based on cubic bezier and returns the abs delta from the mouse cursor to this line
    static inline ImVec2 draw_vertical_sankey_flow(ImDrawList* draw_list, ImVec2 source_pos, ImVec2 dest_pos, float thickness, ImU32 flow_color) {
        ImVec2 p1 = {};
        ImVec2 p2 = {};
        ImVec2 p3 = {};
        ImVec2 p4 = {};
        compute_vertical_sankey_bezier(&p1, &p2, &p3, &p4, source_pos, dest_pos, thickness);

        const int num_segments = 50;

        // Define the color and thickness for the flow
        //ImU32 flow_color = IM_COL32(100, 149, 237, 255);  // Cornflower Blue
        // ImBezierCubicCalc Use this for calculating mouse distance to curve
        // Draw the Bezier curve representing the flow
        draw_list->AddBezierCubic(p1, p2, p3, p4, flow_color, thickness, num_segments);

        ImVec2 mouse_pos = ImGui::GetMousePos();
        float clamped_y = ImClamp(mouse_pos.y, ImMin(p1.y, p4.y), ImMax(p1.y, p4.y));

        ImVec2 point_on_bezier = ImBezierCubicClosestPointCasteljau(p1, p2, p3, p4, mouse_pos, 1.0e-3f);
        ImVec2 delta = mouse_pos - point_on_bezier;
        return {ImAbs(delta.x), ImAbs(clamped_y - mouse_pos.y)};
    }

    static inline void draw_aligned_text(ImDrawList* draw_list, const char* text, ImVec2 pos, ImVec2 alignment = { 0,0 }) {
        draw_list->AddText(aligned_text_pos(text, pos, alignment), ImGui::ColorConvertFloat4ToU32({ 0,0,0,1 }), text);
    }


    // --- Ribbon label placement -----------------------------------------------------------
    //
    // A ribbon may carry a percentage label. Labels are placed at discrete slots ALONG the
    // ribbon's own centre-line bezier, in order of decreasing flow magnitude, and a label that
    // finds no free slot is either left overlapping or dropped - it is never pushed somewhere
    // else.
    //
    // The cost is bounded by construction: num_labels * num_slots rectangle tests, one pass, no
    // iteration towards a fixed point. The previous scheme - repeatedly pushing every pair of
    // overlapping labels apart along their curve directions until nothing overlapped - had no
    // such bound and could not terminate at all whenever two labels shared (or mirrored) a
    // direction: they then translated in lockstep and stayed overlapped forever. Straight-up
    // parallel ribbons (every self-flow is vertical) make that the common case rather than the
    // pathological one. It also had nothing tying a label to its ribbon or to the plot area, so
    // even when it did converge it could converge to labels off screen.

    struct FlowLabel {
        ImVec2 bezier[4];   // ribbon centre-line
        float  value;       // flow fraction; doubles as priority
        ImVec2 pos;         // out: label centre
        bool   visible;     // out
    };

    // drop_unplaceable: a label with no free slot is omitted rather than drawn on top of another.
    static inline void layout_flow_labels(FlowLabel* labels, size_t num_labels, ImVec2 text_size, ImVec2 padding, bool drop_unplaceable, md_temp_scope_t temp) {
        if (num_labels == 0) return;

        // Slot positions along the curve, ideal (the middle) first, alternating outwards. Kept
        // well inside [0,1] so a displaced label still sits on its own ribbon, between the bars.
        static const float slot_t[] = { 0.5f, 0.40f, 0.60f, 0.31f, 0.69f, 0.23f, 0.77f, 0.16f, 0.84f };

        const ImVec2 half = text_size * 0.5f + padding;

        int* order = md_temp_alloc_array(temp, int, num_labels);
        for (size_t i = 0; i < num_labels; ++i) {
            order[i] = (int)i;
        }
        // Biggest flow gets first pick of the slots, so what gets dropped is what matters least.
        std::sort(order, order + num_labels, [labels](int a, int b) {
            return labels[a].value > labels[b].value;
        });

        ImRect* placed = md_temp_alloc_array(temp, ImRect, num_labels);
        size_t num_placed = 0;

        for (size_t i = 0; i < num_labels; ++i) {
            FlowLabel& label = labels[order[i]];
            label.pos = ImBezierCubicCalc(label.bezier[0], label.bezier[1], label.bezier[2], label.bezier[3], slot_t[0]);
            label.visible = false;

            for (size_t s = 0; s < MD_COUNTOF(slot_t); ++s) {
                const ImVec2 pos = ImBezierCubicCalc(label.bezier[0], label.bezier[1], label.bezier[2], label.bezier[3], slot_t[s]);
                const ImRect rect = { pos - half, pos + half };

                bool slot_free = true;
                for (size_t p = 0; p < num_placed; ++p) {
                    if (rect.Overlaps(placed[p])) {
                        slot_free = false;
                        break;
                    }
                }

                if (slot_free) {
                    label.pos = pos;
                    label.visible = true;
                    placed[num_placed++] = rect;
                    break;
                }
            }

            if (!label.visible && !drop_unplaceable) {
                // Every slot is taken, but the caller asked for all the labels. Leave this one at
                // its ideal position and let it overlap.
                label.visible = true;
                placed[num_placed++] = ImRect{ label.pos - half, label.pos + half };
            }
        }
    }

    static inline void vg_sankey_diagram(md_vg_scene_t* scene, ImRect area, Nto* nto, bool hide_text_overlaps) {
        if (!scene) return;
        if (!nto->transition_density_hole) return;
        if (!nto->transition_density_part) return;

        md_temp_scope_t temp = md_temp_begin_avoid(scene->alloc);
        MD_DEFER { md_temp_end(temp); };

        md_vg_add_rect_filled(scene, vec2_from_imvec2(area.Min), vec2_from_imvec2(area.Max), MD_VG_COLOR_RGB(255, 255, 255), 0.0f);

        ImRect plot_area = area;
        const float plot_percent = 0.8f;
        plot_area.Expand({ area.GetWidth() * -(1 - plot_percent), area.GetHeight() * -(1 - plot_percent) });

        const float bar_height = plot_area.GetHeight() * 0.05f;
        int num_bars = 0;
        for (size_t i = 0; i < nto->group.count; ++i) {
            if (nto->transition_density_hole[i] != 0.0f) {
                num_bars += 1;
            }
        }
        const int num_gaps = MD_MAX(num_bars - 1, 0);
        const float gap_size = plot_area.GetWidth() * 0.05f;
        const float bars_avail_width = plot_area.GetWidth() - gap_size * num_gaps;

        float hole_sum = 0.0f;
        float part_sum = 0.0f;
        for (size_t i = 0; i < nto->group.count; ++i) {
            hole_sum += nto->transition_density_hole[i];
            part_sum += nto->transition_density_part[i];
        }
        hole_sum = MD_MAX(hole_sum, 1.0e-6f);
        part_sum = MD_MAX(part_sum, 1.0e-6f);

        float* hole_percentages = md_temp_alloc_array(temp, float, nto->group.count);
        float* part_percentages = md_temp_alloc_array(temp, float, nto->group.count);
        for (size_t i = 0; i < nto->group.count; ++i) {
            hole_percentages[i] = nto->transition_density_hole[i] / hole_sum;
            part_percentages[i] = nto->transition_density_part[i] / part_sum;
        }

        float* start_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float cur_bottom_pos = plot_area.Min.x;
        for (size_t i = 0; i < nto->group.count; ++i) {
            start_positions[i] = cur_bottom_pos;
            cur_bottom_pos += bars_avail_width * hole_percentages[i];
            if (hole_percentages[i] != 0.0f) {
                cur_bottom_pos += gap_size;
            }
        }

        float* end_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float* sub_end_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float cur_pos = plot_area.Min.x;
        for (size_t end_i = 0; end_i < nto->group.count; ++end_i) {
            end_positions[end_i] = cur_pos;
            sub_end_positions[end_i] = cur_pos;
            cur_pos += bars_avail_width * part_percentages[end_i];
            if (hole_percentages[end_i] != 0.0f) {
                cur_pos += gap_size;
            }
        }

        FlowLabel* flow_labels = md_temp_alloc_array(temp, FlowLabel, nto->group.count * nto->group.count);
        size_t num_flow_labels = 0;

        const float font_size = ImGui::GetFontSize();
        const md_vg_color_t text_color = MD_VG_COLOR_RGB(0, 0, 0);
        const md_vg_color_t outline_color = MD_VG_COLOR_RGBA(0, 0, 0, 128);

        for (size_t beg_i = 0; beg_i < nto->group.count; ++beg_i) {
            if (hole_percentages[beg_i] == 0.0f) continue;

            ImVec2 start_pos = { start_positions[beg_i], plot_area.Max.y - bar_height + 0.1f * bar_height };
            ImVec4 start_col = vec_cast(nto->group.color[beg_i]);
            start_col.w = 0.5f;

            for (size_t end_i = 0; end_i < nto->group.count; ++end_i) {
                const float percentage = nto->transition_matrix[end_i * nto->group.count + beg_i];
                if (percentage <= 0.0f) continue;

                const float width = bars_avail_width * percentage;
                ImVec2 end_pos = { sub_end_positions[end_i], plot_area.Min.y + bar_height - 0.1f * bar_height };
                ImVec4 end_col = vec_cast(nto->group.color[end_i]);

                ImVec2 p1 = {};
                ImVec2 p2 = {};
                ImVec2 p3 = {};
                ImVec2 p4 = {};
                compute_vertical_sankey_bezier(&p1, &p2, &p3, &p4, start_pos, end_pos, width);

                md_vg_style_t flow_style = md_vg_style_stroke(vg_color_from_imvec4(start_col), width);
                flow_style.line_cap = MD_VG_LINE_CAP_BUTT;
                flow_style.line_join = MD_VG_LINE_JOIN_ROUND;
                if (end_i != beg_i) {
                    flow_style.stroke = md_vg_paint_linear_gradient(vec2_from_imvec2({ start_pos.x, start_pos.y }), vec2_from_imvec2({ start_pos.x, end_pos.y }), vg_color_from_imvec4(start_col), vg_color_from_imvec4(end_col));
                }
                md_vg_add_bezier_cubic_styled(scene, vec2_from_imvec2(p1), vec2_from_imvec2(p2), vec2_from_imvec2(p3), vec2_from_imvec2(p4), &flow_style);

                char label[32];
                snprintf(label, sizeof(label), "%3.2f%%", percentage * 100.0f);
                const ImVec2 label_size = ImGui::CalcTextSize(label);
                // A ribbon narrower than its own label cannot carry one at all.
                if (width > label_size.x) {
                    FlowLabel& fl = flow_labels[num_flow_labels++];
                    fl.bezier[0] = p1;
                    fl.bezier[1] = p2;
                    fl.bezier[2] = p3;
                    fl.bezier[3] = p4;
                    fl.value     = percentage;
                }

                start_pos.x += width;
                sub_end_positions[end_i] += width;
            }
        }

        layout_flow_labels(flow_labels, num_flow_labels, ImGui::CalcTextSize("99.99%"), { 2.0f, 1.0f }, hide_text_overlaps, temp);

        for (size_t i = 0; i < num_flow_labels; ++i) {
            if (!flow_labels[i].visible) continue;

            char label[16];
            snprintf(label, sizeof(label), "%3.2f%%", flow_labels[i].value * 100.0f);
            const ImVec2 pos = aligned_text_pos(label, flow_labels[i].pos, { 0.5f, 0.5f });
            md_vg_add_text(scene, vec2_from_imvec2(pos), font_size, text_color, md_str_from_cstr(label));
        }

        bool* show_start_text = md_temp_alloc_array(temp, bool, nto->group.count);
        bool* show_end_text = md_temp_alloc_array(temp, bool, nto->group.count);
        int* size_order = md_temp_alloc_array(temp, int, nto->group.count);
        float* text_sizes = md_temp_alloc_array(temp, float, nto->group.count);
        ImRect* rectangles = md_temp_alloc_array(temp, ImRect, nto->group.count);

        for (size_t i = 0; i < nto->group.count; ++i) {
            show_start_text[i] = true;
            show_end_text[i] = true;
        }

        if (hide_text_overlaps) {
            for (int i = 0; i < (int)nto->group.count; ++i) {
                size_order[i] = i;
                text_sizes[i] = MD_MAX(ImGui::CalcTextSize(nto->group.label[i]).x, ImGui::CalcTextSize("99.99%").x);
                const float midpoint = start_positions[i] + bars_avail_width * hole_percentages[i] * 0.5f;
                rectangles[i] = { midpoint - text_sizes[i] * 0.5f, 0.0f, midpoint + text_sizes[i] * 0.5f, 1.0f };
            }

            std::sort(size_order, size_order + nto->group.count, [values = hole_percentages](int a, int b) {
                return values[a] > values[b];
            });

            for (int i = 1; i < (int)nto->group.count; ++i) {
                for (int j = i - 1; j >= 0; --j) {
                    if (hole_percentages[size_order[i]] == 0.0f) {
                        show_start_text[size_order[i]] = false;
                        break;
                    }
                    if (show_start_text[size_order[j]] && rectangles[size_order[i]].Overlaps(rectangles[size_order[j]])) {
                        show_start_text[size_order[i]] = false;
                        break;
                    }
                }
            }

            for (int i = 0; i < (int)nto->group.count; ++i) {
                size_order[i] = i;
                const float midpoint = end_positions[i] + bars_avail_width * part_percentages[i] * 0.5f;
                rectangles[i] = { midpoint - text_sizes[i] * 0.5f, 0.0f, midpoint + text_sizes[i] * 0.5f, 1.0f };
            }

            std::sort(size_order, size_order + nto->group.count, [values = part_percentages](int a, int b) {
                return values[a] > values[b];
            });

            for (int i = 1; i < (int)nto->group.count; ++i) {
                for (int j = i - 1; j >= 0; --j) {
                    if (hole_percentages[size_order[i]] == 0.0f) {
                        show_end_text[size_order[i]] = false;
                        break;
                    }
                    if (show_end_text[size_order[j]] && rectangles[size_order[i]].Overlaps(rectangles[size_order[j]])) {
                        show_end_text[size_order[i]] = false;
                        break;
                    }
                }
            }
        }

        for (size_t i = 0; i < nto->group.count; ++i) {
            if (hole_percentages[i] == 0.0f) continue;

            const ImVec4 bar_color = vec_cast(nto->group.color[i]);

            char start_label[16];
            char end_label[16];
            snprintf(start_label, sizeof(start_label), "%3.2f%%", hole_percentages[i] * 100.0f);
            snprintf(end_label, sizeof(end_label), "%3.2f%%", part_percentages[i] * 100.0f);

            ImVec2 start_p0 = { start_positions[i], plot_area.Max.y - bar_height };
            ImVec2 start_p1 = { start_positions[i] + bars_avail_width * hole_percentages[i], plot_area.Max.y };
            ImVec2 start_midpoint = { (start_p0.x + start_p1.x) * 0.5f, start_p1.y };

            ImVec2 end_p0 = { end_positions[i], plot_area.Min.y };
            ImVec2 end_p1 = { end_positions[i] + bars_avail_width * part_percentages[i], plot_area.Min.y + bar_height };
            ImVec2 end_midpoint = { (end_p0.x + end_p1.x) * 0.5f, end_p0.y };

            md_vg_add_rect_filled(scene, vec2_from_imvec2(start_p0), vec2_from_imvec2(start_p1), vg_color_from_imvec4(bar_color), 0.0f);
            md_vg_add_rect(scene, vec2_from_imvec2(start_p0), vec2_from_imvec2(start_p1), outline_color, 0.0f, 1.0f);
            if (show_start_text[i]) {
                md_vg_add_text(scene, vec2_from_imvec2(aligned_text_pos(nto->group.label[i], start_midpoint, { 0.5f, -0.2f })), font_size, text_color, md_str_from_cstr(nto->group.label[i]));
                md_vg_add_text(scene, vec2_from_imvec2(aligned_text_pos(start_label, start_midpoint, { 0.5f, -1.2f })), font_size, text_color, md_str_from_cstr(start_label));
            }

            md_vg_add_rect_filled(scene, vec2_from_imvec2(end_p0), vec2_from_imvec2(end_p1), vg_color_from_imvec4(bar_color), 0.0f);
            md_vg_add_rect(scene, vec2_from_imvec2(end_p0), vec2_from_imvec2(end_p1), outline_color, 0.0f, 1.0f);
            if (show_end_text[i]) {
                md_vg_add_text(scene, vec2_from_imvec2(aligned_text_pos(nto->group.label[i], end_midpoint, { 0.5f, 1.2f })), font_size, text_color, md_str_from_cstr(nto->group.label[i]));
                md_vg_add_text(scene, vec2_from_imvec2(aligned_text_pos(end_label, end_midpoint, { 0.5f, 2.2f })), font_size, text_color, md_str_from_cstr(end_label));
            }
        }
    }

    static inline bool export_transition_diagram(ImVec2 size, Nto* nto, bool hide_text_overlaps) {
        if (!nto || !nto->transition_density_hole || !nto->transition_density_part || !nto->transition_matrix) {
            MD_LOG_ERROR("Transition diagram export failed: transition diagram data is not initialized");
            return false;
        }

        char path_buf[1024] = {};
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, MD_STR_LIT("svg"))) {
            return false;
        }

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        md_vg_scene_t scene = {};
        md_vg_scene_init(&scene, { size.x, size.y }, md_temp_allocator(temp));
        vg_sankey_diagram(&scene, { 0.0f, 0.0f, size.x, size.y }, nto, hide_text_overlaps);

        const bool success = md_vg_scene_write_svg_file(&scene, md_str_from_cstr(path_buf));
        if (success) {
            MD_LOG_INFO("Exported transition diagram to '%s'", path_buf);
        } else {
            MD_LOG_ERROR("Failed to export transition diagram to '%s'", path_buf);
        }
        return success;
    }

    static inline void im_sankey_diagram(ApplicationState* state, ImRect area, Nto* nto, bool hide_text_overlaps) {
        if (!nto->transition_density_hole) return;
        if (!nto->transition_density_part) return;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        /*
        * A sankey diagram needs to implement bezier curves, that are connecting two points
        */
        ImDrawList* draw_list = ImGui::GetWindowDrawList();

        //Draw background
        //draw_list->AddRectFilled(area.Min, area.Max, ImGui::ColorConvertFloat4ToU32({ 1,1,1,1 }));
        //draw_list->AddRect(area.Min, area.Max, ImGui::ColorConvertFloat4ToU32({ 0,0,0,1 }));

        ImRect plot_area = area;
        const float plot_percent = 0.8f;
        plot_area.Expand({ area.GetWidth() * -(1 - plot_percent), area.GetHeight() * -(1 - plot_percent) });
        //draw_list->AddRect(plot_area.Min, plot_area.Max, ImGui::ColorConvertFloat4ToU32({ 1,0,0,1 })); //Use this to draw debug of plot area

        ////The given data
        //const char* names[] = { "THIO", "QUIN" };
        //float initial_percentages[2] = { 0.3, 0.7 };
        //float transitions[2][2] = {
        //    /*
        //    *   A B   col -->
        //    * A
        //    * B
        //    * |
        //    * |
        //    * v
        //    * row
        //    */
        //    {0.22, 0.78}, //Read as "Starting from A, 22% of A_Start goes to A_End, and 78% goes to B_End
        //    {0.0, 1.0}
        //};

        //Bar definitions
        const float bar_height = plot_area.GetHeight() * 0.05f;
        int num_bars = 0;
        for (size_t i = 0; i < nto->group.count; i++) {
            if (nto->transition_density_hole[i] != 0.0f) {
                num_bars++;
            }
        }
        int num_gaps = num_bars - 1;
        float gap_size = plot_area.GetWidth() * 0.05f;
        float bars_avail_width = plot_area.GetWidth() - gap_size * num_gaps;

        //Calculate bar percentages
        float hole_sum = 0;
        float part_sum = 0;
        for (size_t i = 0; i < nto->group.count; i++) {
            hole_sum += nto->transition_density_hole[i];
            part_sum += nto->transition_density_part[i];
        }

        // Prevent division by zero
        hole_sum = MD_MAX(hole_sum, 1.0e-6f);
        part_sum = MD_MAX(part_sum, 1.0e-6f);

        float* hole_percentages = md_temp_alloc_array(temp, float, nto->group.count);
        float* part_percentages = md_temp_alloc_array(temp, float, nto->group.count);
        for (size_t i = 0; i < nto->group.count; i++) {
            hole_percentages[i] = nto->transition_density_hole[i] / hole_sum;
            part_percentages[i] = nto->transition_density_part[i] / part_sum;
        }

        //Calculate start positions
        float* start_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float cur_bottom_pos = plot_area.Min.x;
        for (size_t i = 0; i < nto->group.count; i++) {
            start_positions[i] = cur_bottom_pos;
            cur_bottom_pos += bars_avail_width * hole_percentages[i];
            if (hole_percentages[i] != 0.0) {
                cur_bottom_pos += gap_size;
            }
        }

        //Calculate end positions
        float* end_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float* sub_end_positions = md_temp_alloc_array(temp, float, nto->group.count);
        float cur_pos = plot_area.Min.x;
        for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
            end_positions[end_i] = cur_pos;
            sub_end_positions[end_i] = cur_pos;
            cur_pos += bars_avail_width * part_percentages[end_i];
            if (hole_percentages[end_i] != 0.0) {
                cur_pos += gap_size;
            }
        }

        float  hovered_flow_width = 0.0f;
        ImVec2 hovered_flow_beg = {};
        ImVec2 hovered_flow_end = {};
        ImVec2 mouse_pos = ImGui::GetMousePos();

        char mouse_label[32] = { 0 };

        const float min_test_width = 1.0f;

        //Draw curves
        FlowLabel* flow_labels = md_temp_alloc_array(temp, FlowLabel, nto->group.count * nto->group.count);
        size_t num_flow_labels = 0;

        for (size_t beg_i = 0; beg_i < nto->group.count; beg_i++) {
            if (hole_percentages[beg_i] != 0.0) {

                ImVec2 start_pos = { start_positions[beg_i], plot_area.Max.y - bar_height + 0.1f * bar_height };
                ImVec4 start_col = vec_cast(nto->group.color[beg_i]);

                start_col.w = 0.5f;
                for (size_t end_i = 0; end_i < nto->group.count; end_i++) {
                    float percentage = nto->transition_matrix[end_i * nto->group.count + beg_i];
                    ImVec4 end_col = vec_cast(nto->group.color[end_i]);

                    if (percentage > 0.0f) {
                        float width = bars_avail_width * percentage;
                        ImVec2 end_pos = { sub_end_positions[end_i], plot_area.Min.y + bar_height - 0.1f * bar_height };

                        int vert_beg = draw_list->VtxBuffer.Size;
                        ImVec2 mouse_delta = draw_vertical_sankey_flow(draw_list, start_pos, end_pos, width, ImGui::ColorConvertFloat4ToU32(start_col));
                        int vert_end = draw_list->VtxBuffer.Size;

                        // The width test here is clamped to be atleast some pixel in width, otherwise it would be impossible to hit thin flows
                        bool flow_hovered = mouse_delta.x < MD_MAX(width * 0.5f, min_test_width) && mouse_delta.y < 1.0e-4;

                        // Apply a gradient based on y-value from start color to end color if they belong to different groups
                        if (end_i != beg_i) {
                            ImVec2 grad_p0 = { start_pos.x, start_pos.y };
                            ImVec2 grad_p1 = { start_pos.x, end_pos.y };
                            ImGui::ShadeVertsLinearColorGradientKeepAlpha(draw_list, vert_beg, vert_end, grad_p0, grad_p1, ImGui::ColorConvertFloat4ToU32(start_col), ImGui::ColorConvertFloat4ToU32(end_col));
                        }

                        char label[32];
                        snprintf(label, sizeof(label), "%3.2f%%", percentage * 100);
                        const ImVec2 label_size = ImGui::CalcTextSize(label);

                        // @NOTE: The comparison for width here is to ensure that the narrowest flow is always the one hovered (when multiple flows are overlapping)
                        if (flow_hovered && (hovered_flow_width == 0.0f || width < hovered_flow_width)) {
                            hovered_flow_width = width;
                            hovered_flow_beg = start_pos;
                            hovered_flow_end = end_pos;
                            md_memcpy(mouse_label, label, sizeof(label));
                        }

                        // A ribbon narrower than its own label cannot carry one at all.
                        if (width > label_size.x) {
                            FlowLabel& fl = flow_labels[num_flow_labels++];
                            compute_vertical_sankey_bezier(&fl.bezier[0], &fl.bezier[1], &fl.bezier[2], &fl.bezier[3], start_pos, end_pos, width);
                            fl.value = percentage;
                        }

                        start_pos.x += width;
                        sub_end_positions[end_i] += width;
                    }
                }
            }
        }

        // Place the ribbon labels. Bounded, single pass - see layout_flow_labels.
        layout_flow_labels(flow_labels, num_flow_labels, ImGui::CalcTextSize("99.99%"), { 2.0f, 1.0f }, hide_text_overlaps, temp);

        for (size_t i = 0; i < num_flow_labels; ++i) {
            if (!flow_labels[i].visible) continue;

            char label[16];
            snprintf(label, sizeof(label), "%3.2f%%", flow_labels[i].value * 100.0f);
            draw_aligned_text(draw_list, label, flow_labels[i].pos, { 0.5f, 0.5f });
        }

        // Some flow is hovered, So we render a new flow on top of the existing ones to ensure that this is rendered on top
        if (hovered_flow_width > 0.0f) {
            draw_vertical_sankey_flow(draw_list, hovered_flow_beg, hovered_flow_end, hovered_flow_width, ImColor(1.0f, 1.0f, 1.0f, 0.2f));
        }

        //Calculate text visibility
        bool* show_start_text  = md_temp_alloc_array(temp, bool,   nto->group.count);
        bool* show_end_text    = md_temp_alloc_array(temp, bool,   nto->group.count);
        int*  size_order       = md_temp_alloc_array(temp, int,    nto->group.count);
        float* text_sizes      = md_temp_alloc_array(temp, float,  nto->group.count);
        ImRect* rectangles     = md_temp_alloc_array(temp, ImRect, nto->group.count);

        for (size_t i = 0; i < nto->group.count; i++) {
            show_start_text[i] = true;
            show_end_text[i] = true;
        }
        
        if (hide_text_overlaps) {
            //Start texts
            for (int i = 0; i < (int)nto->group.count; i++) {
                size_order[i] = i;
                text_sizes[i] = MD_MAX(ImGui::CalcTextSize(nto->group.label[i]).x, ImGui::CalcTextSize("99.99%").x);

                float midpoint = start_positions[i] + bars_avail_width * hole_percentages[i] * 0.5f;
                rectangles[i] = { midpoint - text_sizes[i] * 0.5f, 0, midpoint + text_sizes[i] * 0.5f, 1};
            }
        
            std::sort(size_order, size_order + nto->group.count, [values = hole_percentages](int a, int b) {
                return values[a] > values[b];  // Sort in descending order
            });

            for (int i = 1; i < (int)nto->group.count; i++) { //First one is always drawn
                for (int j = i - 1; j >= 0 ; j--) {//The bigger ones
                    if (hole_percentages[size_order[i]] == 0.0f) {
                        show_start_text[size_order[i]] = false;
                        break;
                    }
                    else if (show_start_text[size_order[j]] && rectangles[size_order[i]].Overlaps(rectangles[size_order[j]])) {
                        show_start_text[size_order[i]] = false;
                        break;
                    }
                }
            }

            //End texts
            for (int i = 0; i < (int)nto->group.count; i++) {
                size_order[i] = i;

                float midpoint = end_positions[i] + bars_avail_width * part_percentages[i] * 0.5f;
                rectangles[i] = { midpoint - text_sizes[i] * 0.5f, 0, midpoint + text_sizes[i] * 0.5f, 1 };
            }

            std::sort(size_order, size_order + nto->group.count, [values = part_percentages](int a, int b) {
                return values[a] > values[b];  // Sort in descending order
            });

            for (int i = 1; i < (int)nto->group.count; i++) { //First one is always drawn
                for (int j = i - 1; j >= 0; j--) {//The bigger ones
                    if (hole_percentages[size_order[i]] == 0.0f) {
                        show_end_text[size_order[i]] = false;
                        break;
                    }
                    else if (show_end_text[size_order[j]] && rectangles[size_order[i]].Overlaps(rectangles[size_order[j]])) {
                        show_end_text[size_order[i]] = false;
                        break;
                    }
                }
            }
        }

        //Draw bars
        for (size_t i = 0; i < nto->group.count; i++) {
            if (hole_percentages[i] != 0.0f) {
                ImVec4 bar_color = vec_cast(nto->group.color[i]);

                char start_label1[16];
                char end_label1[16]; 
                snprintf(start_label1, sizeof(start_label1), "%3.2f%%", hole_percentages[i] * 100);
                snprintf(end_label1, sizeof(end_label1), "%3.2f%%", part_percentages[i] * 100);

                //Calculate start
                ImVec2 start_p0 = { start_positions[i], plot_area.Max.y - bar_height };
                ImVec2 start_p1 = { start_positions[i] + bars_avail_width * hole_percentages[i], plot_area.Max.y };
                ImVec2 start_midpoint = { (start_p0.x + start_p1.x) * 0.5f, start_p1.y };
                ImRect start_bar = ImRect{ start_p0, start_p1 };

                //Calculate end
                ImVec2 end_p0 = { end_positions[i], plot_area.Min.y };
                ImVec2 end_p1 = { end_positions[i] + bars_avail_width * part_percentages[i], plot_area.Min.y + bar_height };
                ImRect end_bar = ImRect{ end_p0, end_p1 };
                ImVec2 end_midpoint = { (end_p0.x + end_p1.x) * 0.5f, end_p0.y };

                bool index_hovered = false;
                if (start_bar.Contains(mouse_pos)) {
                    md_memcpy(mouse_label, start_label1, sizeof(start_label1));
                    index_hovered = true;
                }
                else if (end_bar.Contains(mouse_pos)) {
                    md_memcpy(mouse_label, end_label1, sizeof(end_label1));
                    index_hovered = true;
                }

                if (index_hovered) {
                    nto->group.hovered_index = (int8_t)i;
                    bar_color = make_highlight_color(bar_color);
                    for (size_t atom_i = 0; atom_i < md_array_size(nto->atom_group_idx); atom_i++) {
                        if (nto->atom_group_idx[atom_i] == (uint32_t)i) {
                            md_bitmask_set_bit(&state->selection.highlight_mask, (uint32_t)atom_i);
                        }
                    }
                }

                //Draw start
                draw_list->AddRectFilled(start_p0, start_p1, ImGui::ColorConvertFloat4ToU32(bar_color));
                draw_list->AddRect(start_p0, start_p1, ImGui::ColorConvertFloat4ToU32({ 0,0,0,0.5 }));
                if (show_start_text[i]) {
                    draw_aligned_text(draw_list, nto->group.label[i], start_midpoint, { 0.5, -0.2 });
                    draw_aligned_text(draw_list, start_label1, start_midpoint, { 0.5, -1.2 });
                }

                //Draw end
                draw_list->AddRectFilled(end_p0, end_p1, ImGui::ColorConvertFloat4ToU32(bar_color));
                draw_list->AddRect(end_p0, end_p1, ImGui::ColorConvertFloat4ToU32({ 0,0,0,0.5 }));
                if (show_end_text[i]) {
                    draw_aligned_text(draw_list, nto->group.label[i], end_midpoint, { 0.5, 1.2 });
                    draw_aligned_text(draw_list, end_label1, end_midpoint, { 0.5, 2.2 });
                }
            }
        }

        size_t mouse_label_len = strnlen(mouse_label, sizeof(mouse_label));
        if (mouse_label_len > 0) {
            ImVec2 offset = { 15, 15 };
            ImVec2 pos = ImGui::GetMousePos() + offset;
            draw_list->AddText(pos, IM_COL32_BLACK, mouse_label, mouse_label + mouse_label_len);
        }
    }

    static inline double convert_value_from_au(double value, x_unit_t unit) {
        switch (unit) {
        case X_UNIT_EV:
            return value * HARTREE_TO_EV;
        case X_UNIT_NM:
            return 45.56331418628337 / value;
        case X_UNIT_CM_INVERSE:
            return value * 219479.86946633097;
        case X_UNIT_HARTREE:
            return value;
        default:
            MD_ASSERT(false); // Should not happen
            return value;
        }
    }

    static inline double convert_value_to_au(double value, x_unit_t unit) {
        switch (unit) {
        case X_UNIT_EV:
            return value * EV_TO_HARTREE;
        case X_UNIT_NM:
            return 45.56331418628337 / value;
        case X_UNIT_CM_INVERSE:
            return value / 219479.86946633097;
        case X_UNIT_HARTREE:
            return value;
        default:
            MD_ASSERT(false); // Should not happen
            return value;
        }
    }

    static inline void convert_values_from_au(double* out_values, const double* in_values, size_t num_values, x_unit_t unit) {
        if (!out_values || !in_values || num_values == 0) return;

        switch (unit) {
        case X_UNIT_EV:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i] * HARTREE_TO_EV;
            }
            break;
        case X_UNIT_NM:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = 45.56331418628337 / in_values[i];
            }
            break;
        case X_UNIT_CM_INVERSE:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i] * 219479.86946633097;
            }
            break;
        case X_UNIT_HARTREE:
            for (size_t i = 0; i < num_values; ++i) {
                out_values[i] = in_values[i];
            }
            break;
        default:
            MD_ASSERT(false); // Should not happen
            break;
        }
    }

    static inline double sigma_from_fwhm(double fwhm) {
        return fwhm * 0.21233045007200476068;
    }

    // Evaluates lorentzian broadening around a point x, given a set of peaks and their strengths.
    // This is a NON normalized version.
    // To normalize this, just multiply with the normalization factor
	static inline double lorentzian_base(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
		const double gamma = fwhm * 0.5;
		const double gamma2 = gamma * gamma;

		double y = 0;
		for (size_t k = 0; k < num_peaks; k++) {
			double dx = x - peaks_x[k];
			y += peaks_y[k] * gamma / fma(dx, dx, gamma2);
		}
		return y;
	}

    static inline double lorentzian_abs(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double gamma = fwhm * 0.5;
		const double gamma2 = gamma * gamma;
        double y = 0.0;
		for (size_t k = 0; k < num_peaks; k++) {
			double dx = x - peaks_x[k];
			y += peaks_y[k] / peaks_x[k] * gamma / fma(dx, dx, gamma2);
		}
        return y / MD_PI;
    }

    static inline double lorentzian_ecd(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return lorentzian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    static inline double lorentzian_vib(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return (1.0 / MD_PI) * lorentzian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    static inline double lorentzian_tpa(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return x * x * lorentzian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    static inline double lorentzian_xps(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return (1.0 / MD_PI) * lorentzian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    // Evaluates gaussian broadening around a point x, given a set of peaks and their strengths.
    // This is a NON normalized version.
    // To normalize this, just multiply with the normalization factor
	static inline double gaussian_base(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
		const double sigma = sigma_from_fwhm(fwhm);
		const double factor = -0.5 / (sigma * sigma);

        double y = 0.0;
        for (size_t k = 0; k < num_peaks; k++) {
            double dx = x - peaks_x[k];
            y += peaks_y[k] * exp(factor * dx * dx);
        }
		return y;
	}

    static inline double gaussian_abs(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double sigma = sigma_from_fwhm(fwhm);
		const double factor = -0.5 / (sigma * sigma);
        double y = 0.0;
        for (size_t k = 0; k < num_peaks; k++) {
            double dx = x - peaks_x[k];
            y += peaks_y[k] / peaks_x[k] * exp(factor * dx * dx);
        }
		return y / (sigma * sqrt(2.0 * MD_PI));
    }

    static inline double gaussian_ecd(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double sigma = sigma_from_fwhm(fwhm);
		const double y = gaussian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
		return y * MD_PI / (sigma * sqrt(2.0 * MD_PI));
    }

    static inline double gaussian_vib(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double sigma = sigma_from_fwhm(fwhm);
		const double y = gaussian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
		return y * sqrt(2.0) / (sigma * sqrt(MD_PI));
    }

    static inline double gaussian_xps(double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        const double sigma = sigma_from_fwhm(fwhm);
		const double y = gaussian_base(x, peaks_x, peaks_y, num_peaks, fwhm);
		return y / (sigma * sqrt(2.0 * MD_PI));
    }


    //Constructs plot limits from peaks
    static inline ImPlotRect get_plot_limits(const double* x_samples, const double* y_peaks, size_t num_peaks, size_t num_samples, double ext_fac = 0.1) {
        ImPlotRect lim = { MD_MIN(x_samples[0],x_samples[num_samples - 1]), MD_MAX(x_samples[0],x_samples[num_samples - 1]),0,0};
        for (size_t i = 0; i < num_peaks; i++) {
            //Use Contains to check if values are within the limits, or if they should extend the limits
            if (lim.Y.Max < y_peaks[i]) {
                lim.Y.Max = y_peaks[i];
            }
            else if (lim.Y.Min > y_peaks[i]) {
                lim.Y.Min = y_peaks[i];
            }
        }

        double height = lim.Y.Max - lim.Y.Min;
        double width = lim.X.Max - lim.X.Min;

        lim.Y.Max += height * ext_fac;
        lim.Y.Min -= height * ext_fac;
        lim.X.Max += width * ext_fac;
        lim.X.Min -= width * ext_fac;

        //The y limits needs to be symmetric so that the spectra is not clipped. For example, if the y peak max is a positive value but the y spectra value is negative, we get issues otherwise. This ensures space for all the data.
        //This is needed because Y2 scaling is based on the maximum peak value
        double abs_max_y = MD_MAX(fabs(lim.Y.Min), fabs(lim.Y.Max));
        lim.Y.Min = -abs_max_y;
        lim.Y.Max = abs_max_y;
        return lim;
    }

    //converts x and y peaks into pixel points in context of the current plot. Use between BeginPlot() and EndPlot()
    static inline void peaks_to_pixels(ImVec2* pixel_peaks, const double* x_peaks, const double* y_peaks, size_t num_peaks) {
        for (size_t i = 0; i < num_peaks; i++) {
            pixel_peaks[i] = ImPlot::PlotToPixels(ImPlotPoint{ x_peaks[i], y_peaks[i] });
        }
    }

    static inline int get_hovered_plot_peak(const double* xs, const double* ys, size_t count, double distance_threshold = 10.0) {
        ImPlotPoint mouse_pos = ImPlot::GetPlotMousePos();
        int closest_idx = -1;
        double closest_dist = DBL_MAX;
        for (size_t i = 0; i < count; ++i) {
            double dx = mouse_pos.x - xs[i];
            double dy = MD_CLAMP(mouse_pos.y, 0.0, ys[i]) - mouse_pos.y;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist < distance_threshold && dist < closest_dist) {
                closest_idx = (int)i;
                closest_dist = dist;
            }
        }
        return closest_idx;
    }

    static inline int get_hovered_plot_peak(const double* values, size_t count, double x_base = 0.0, double x_step = 1.0, double distance_threshold = 10.0) {
        ImPlotPoint mouse_pos = ImGui::GetMousePos();
        int closest_idx = -1;
        double closest_dist = DBL_MAX;
        ImVec2 zero = ImPlot::PlotToPixels(0, 0);
        for (size_t i = 0; i < count; ++i) {
            ImVec2 point = ImPlot::PlotToPixels(x_base + i * x_step, values[i]);
            double dx = mouse_pos.x - point.x;
            double dy = MD_CLAMP(mouse_pos.y, MD_MIN(zero.y, point.y), MD_MAX(zero.y, point.y)) - mouse_pos.y;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist < distance_threshold && dist < closest_dist) {
                closest_idx = (int)i;
                closest_dist = dist;
            }
        }
        return closest_idx;
    }

    // Returns peak index closest to mouse pixel position, assumes that x-values are sorted.
    static inline int get_hovered_peak(const ImVec2 mouse_pos, const ImVec2* pixel_peaks, const ImVec2* pixel_points, size_t num_peaks, bool y_flipped = false, double proxy_distance = 10.0) {
        int closest_idx = 0;
        double x = mouse_pos.x;
        double y = mouse_pos.y;
        double y_max = 0;
        double y_min = 0;
        double distance_x = 0;
        double distance_y = 0;
        double closest_distance = 0;
        double pixel_y0 = ImPlot::PlotToPixels(0, 0).y;

        //Keep in mind that pixel y is 0 at the top, so you flip the comparison compared to plot y. The code below still seems to work as intended though.

        for (size_t i = 0; i < num_peaks; i++) {
            y_max = MD_MAX(pixel_peaks[i].y, pixel_y0);
            y_min = MD_MIN(pixel_peaks[i].y, pixel_y0);

            //Check if the y location is within the range of y_min,ymax
            if (y > y_max) {
                distance_y = fabs(y - y_max);
                if (y_flipped) {
                    distance_y = fabs(y - pixel_points[i].y) < distance_y ? fabs(y - pixel_points[i].y) : distance_y;
                }
            }
            else if (y < y_min) {
                distance_y = fabs(y - y_min);
                if (!y_flipped) {
                    distance_y = fabs(y - pixel_points[i].y) < distance_y ? fabs(y - pixel_points[i].y) : distance_y;
                }
            }
            else {
                distance_y = 0;
            }

            distance_x = fabs(pixel_peaks[i].x - x);

            // We need a special case for i == 0 to set a reference for comparison, so that closest_distance does not start on 0;
            if (i == 0 && distance_y == 0 ){
                closest_distance = distance_x;
                closest_idx = 0;

            }
            else if (i == 0) {
                closest_distance = sqrt(pow(distance_x, 2) + pow(distance_y, 2)); // Is there a better function for doing this? Is this expensive?
                closest_idx = 0;
            }
            else if (distance_y == 0 && distance_x < closest_distance) {
                closest_distance = fabs(pixel_peaks[i].x - x);
                closest_idx = (int)i;
            }
            else if (sqrt(pow(distance_x, 2) + pow(distance_y, 2)) < closest_distance) {
                closest_distance = sqrt(pow(distance_x, 2) + pow(distance_y, 2));
                closest_idx = (int)i;
                //ImPlot::Annotation()
            }
            else if (distance_x > closest_distance){
                // We are now so far away that a closer bar will not occur, no matter the y value.
                break;
            }
        }
        return closest_distance < proxy_distance ? closest_idx : -1;
    }

    // Draws a bar within a plot at x = x, and from y = 0 to y = y, with specified width and color.
    // It also draws a point at the top with a radius
    static inline void draw_bar(double x, double y, double bar_width_in_pixels, double point_radius_in_pixels, ImVec4 color) {
        ImVec2 rect_min = ImPlot::PlotToPixels(x, 0, IMPLOT_AUTO, IMPLOT_AUTO);
        ImVec2 rect_max = ImPlot::PlotToPixels(x, y, IMPLOT_AUTO, IMPLOT_AUTO);

        rect_min = ImVec2{rect_min.x - (float)bar_width_in_pixels * 0.5f, rect_min.y};
        rect_max = ImVec2{rect_max.x + (float)bar_width_in_pixels * 0.5f, rect_max.y};
        ImVec2 pos = ImPlot::PlotToPixels(x, y, IMPLOT_AUTO, IMPLOT_AUTO);
        float  rad = (float)point_radius_in_pixels;

        ImPlot::PushPlotClipRect();
        ImDrawList& DrawList = *ImPlot::GetPlotDrawList();

        ImU32 col32 = ImGui::ColorConvertFloat4ToU32(color);
        DrawList.AddRectFilled(rect_min, rect_max, col32);
        DrawList.AddCircleFilled(pos, rad, col32);

        ImPlot::PopPlotClipRect();
    }

    // Compute a simple display multiplier to map oscillator-strength values (f) to eps for axis-syncing.
    // Strategy: Match peak magnitudes: multiplier = max(abs(eps_samples))/max(abs(f_peaks)).
    double compute_display_multiplier_from_peaks(const double* eps_samples, size_t eps_count, const double* peak_values, size_t peak_count) {
        double eps_max = 0.0;
        for (size_t i = 0; i < eps_count; ++i) {
            eps_max = MD_MAX(eps_max, fabs(eps_samples[i]));
        }
        double peak_max = 0.0;
        for (size_t i = 0; i < peak_count; ++i) {
            peak_max = MD_MAX(peak_max, fabs(peak_values[i]));
        }
        if (peak_max <= 0.0 || !isfinite(peak_max) || !isfinite(eps_max)) return 1.0;

        return eps_max / peak_max;
    }

    enum PlotPeaksFlags_ {
        PlotPeaksFlags_None = 0,
        PlotPeaksFlags_Bars = 1 << 0,
        PlotPeaksFlags_Points = 1 << 1,
        PlotPeaksFlags_All = PlotPeaksFlags_Bars | PlotPeaksFlags_Points
    };

    typedef uint32_t PlotPeaksFlags;

    // This is a helper function to plot peaks and handle hovering
    // Only valid to call within a plot context
    void plot_peaks(const char* title, const double x_values[], const double y_values[], size_t num_peaks, int& selected, int& hovered, PlotPeaksFlags flags = PlotPeaksFlags_All) {
        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        const double bar_width_in_pixels = 2.0;
        const double point_radius_in_pixels = 3.0;

        // If x values are not given we implicitly assume that they correspond to 1..num_peaks
        if (!x_values) {
            double* x_vals = md_temp_alloc_array(temp, double, num_peaks);
            for (size_t i = 0; i < num_peaks; i++) {
                x_vals[i] = (double)(i + 1);
            }
            x_values = x_vals;
        }

        // An item hidden through the legend draws nothing, so it must not claim a hover or take a
        // click either. ImPlot only resolves that inside BeginItem, which happens further down -
        // after both would already have run - so the item is looked up directly instead. It is
        // absent on the very first frame it is drawn, which BeginItem will treat as shown.
        const ImPlotItem* item = ImPlot::GetItem(title);
        const bool item_visible = (item == nullptr) || item->Show;

        // Dropped whenever the cursor is inside the plot, hidden or not, so that a hover taken
        // before the item was toggled off cannot keep driving a tooltip in the caller.
        if (ImPlot::IsPlotHovered()) {
            hovered = -1;
        }

        if (item_visible && ImPlot::IsPlotHovered()) {
            // Calculate the hovered index as the closest peak or 'bar' to the mouse position
            // This comparison needs to be done in pixel space

            const ImVec2 mouse_pos = ImGui::GetMousePos();
            const double pixel_tolerance = 7.0;
            const double min_line_d2  = pow(bar_width_in_pixels * 0.5 + pixel_tolerance, 2);  // 7 pixels tolerance
            const double min_point_d2 = pow(point_radius_in_pixels + pixel_tolerance, 2);     // 7 pixels tolerance

            double min_dist = DBL_MAX;
            for (size_t i = 0; i < num_peaks; i++) {
                ImVec2 min_raw = ImPlot::PlotToPixels(x_values[i], 0.0, IMPLOT_AUTO, IMPLOT_AUTO);
                ImVec2 max_raw = ImPlot::PlotToPixels(x_values[i], y_values[i], IMPLOT_AUTO, IMPLOT_AUTO);

                // The axis may be flipped, so we need to determine min and max correctly
                ImVec2 min = ImMin(min_raw, max_raw);
                ImVec2 max = ImMax(min_raw, max_raw);

                if (flags & PlotPeaksFlags_Bars) {
                    // Calculate distance to the line which forms the bar (1D clamp trick)
                    ImVec2 line_dx = mouse_pos - ImClamp(mouse_pos, min, max);
                    double line_d2 = line_dx.x * line_dx.x + line_dx.y * line_dx.y;

                    if (line_d2 < min_line_d2 && line_d2 < min_dist) {
                        hovered = (int)i;
                        min_dist = line_d2;
                    }
                }

                if (flags & PlotPeaksFlags_Points) {
                    // Calculate distance to the point at the top of the bar
                    ImVec2 point_dx = mouse_pos - max_raw;
                    double point_d2 = point_dx.x * point_dx.x + point_dx.y * point_dx.y;

                    if (point_d2 < min_point_d2 && point_d2 < min_dist) {
                        hovered = (int)i;
                        min_dist = point_d2;
                    }
                }
            }
        }

        // Update selected peak on click
        if (item_visible && ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) &&
            ImPlot::IsPlotHovered()) {
            selected = hovered == selected ? -1 : hovered;
        }

        if (ImPlot::BeginItem(title)) {
            ImPlot::PushPlotClipRect();
            ImDrawList& DrawList = *ImPlot::GetPlotDrawList();
            const float rad = (float)point_radius_in_pixels;

            for (int i = 0; i < (int)num_peaks; i++) {
                const ImPlotNextItemData& s = ImPlot::GetItemData();
                ImVec4 color = (i == hovered) ? COLOR_PEAK_HOVER : (i == selected) ? COLOR_PEAK_SELECTED : s.Colors[ImPlotCol_Fill];

                const double x = x_values[i];
                const double y = y_values[i];

                ImVec2 rect_min = ImPlot::PlotToPixels(x, 0, IMPLOT_AUTO, IMPLOT_AUTO);
                ImVec2 rect_max = ImPlot::PlotToPixels(x, y, IMPLOT_AUTO, IMPLOT_AUTO);

                rect_min = ImVec2{rect_min.x - (float)bar_width_in_pixels * 0.5f, rect_min.y};
                rect_max = ImVec2{rect_max.x + (float)bar_width_in_pixels * 0.5f, rect_max.y};
                ImVec2 pos = ImPlot::PlotToPixels(x, y, IMPLOT_AUTO, IMPLOT_AUTO);

                ImU32 col32 = ImGui::ColorConvertFloat4ToU32(color);
                if (flags & PlotPeaksFlags_Bars) {
                    DrawList.AddRectFilled(rect_min, rect_max, col32);
                }

                if (flags & PlotPeaksFlags_Points) {
                    DrawList.AddCircleFilled(pos, rad, col32);
                }

                if (ImPlot::FitThisFrame()) {
                    ImPlot::FitPoint(ImPlotPoint(x, 0.0));
                    ImPlot::FitPoint(ImPlotPoint(x, y));
                }
            }
            ImPlot::PopPlotClipRect();

            ImPlot::EndItem();
        }
    }

    // =================================================================================================
    // RIXS map & spectrum
    //
    // ImGui/ImPlot port of VeloxChem's spectrumplot.plot_rixs_map(), plus a 1D spectrum widget that
    // shows a single row of that map (one incoming photon energy).
    // The implementation is split into a pure data pass (compute_rixs_map) and a drawing pass
    // (draw_rixs_map) so that the map can be rebuilt only when the settings or the input actually
    // change. Both take the input data through RixsMapInput and know nothing about where it came
    // from; rixs_map_input_gather() reads it out of the system's attribute table.
    // =================================================================================================

    // Raw inputs, mirroring the 'rixs_results' dict consumed by plot_rixs_map().
    // Dimensions:
    //   P = num_photon_energies (incoming photons)
    //   F = num_final_states    (valence / final states)
    //   C = num_core_states     (core-excited states, only used for the XAS panel)
    // The 2D arrays are row-major [F][P], i.e. element (f, p) is at [f * P + p].
    struct RixsMapInput {
        const double* photon_energies_au     = nullptr;  // [P], a.u.
        const double* elastic_cross_sections = nullptr;  // [P], a.u.    (optional)
        const double* cross_sections         = nullptr;  // [F][P], a.u.
        const double* energy_losses_au       = nullptr;  // [F][P], a.u.
        const double* core_eigenvalues_au    = nullptr;  // [C], a.u.    (optional)
        const double* core_osc_strengths     = nullptr;  // [C]          (optional)
        double gamma_fwhm_ev = 0.0;                      // core-hole lifetime broadening, eV
        size_t num_photon_energies = 0;
        size_t num_final_states    = 0;
        size_t num_core_states     = 0;
    };

    // Collects the RIXS input out of the system's attribute table, under vlx/rsp/rixs/.
    // Returns false if the system carries no RIXS data or carries it incompletely.
    //
    // The extents come from the attributes' OWN SHAPES rather than from three published counts: a
    // {F,P} matrix cannot disagree with its own dimensions, and the one thing a separate count could
    // add is a way for it to. The photon axis is checked across the three quantities that share it
    // for the same reason - siblings in one group share an index space, and a file where they do not
    // is not a map anyone can draw.
    static bool rixs_map_input_gather(RixsMapInput& out, const md_system_t& sys) {
        if (response_type(sys) != RSP_RIXS) return false;

        size_t num_photon   = 0;
        size_t num_elastic  = 0;
        size_t num_core     = 0;
        size_t num_core_osc = 0;

        out.photon_energies_au     = attribute_series_f64(&num_photon,   sys, MD_STR_LIT("vlx/rsp/rixs/photon_energy"));
        out.elastic_cross_sections = attribute_series_f64(&num_elastic,  sys, MD_STR_LIT("vlx/rsp/rixs/elastic_cross_section"));
        out.core_eigenvalues_au    = attribute_series_f64(&num_core,     sys, MD_STR_LIT("vlx/rsp/rixs/core_energy"));
        out.core_osc_strengths     = attribute_series_f64(&num_core_osc, sys, MD_STR_LIT("vlx/rsp/rixs/core_oscillator_strength"));
        out.gamma_fwhm_ev          = attribute_scalar_f64(sys, MD_STR_LIT("vlx/rsp/rixs/gamma_fwhm"));

        size_t cs_rows = 0, cs_cols = 0;
        size_t el_rows = 0, el_cols = 0;
        out.cross_sections   = attribute_matrix_f64(&cs_rows, &cs_cols, sys, MD_STR_LIT("vlx/rsp/rixs/cross_section"));
        out.energy_losses_au = attribute_matrix_f64(&el_rows, &el_cols, sys, MD_STR_LIT("vlx/rsp/rixs/energy_loss"));

        // The map itself needs the photon energies, the cross sections and the energy losses.
        if (!out.photon_energies_au || !out.cross_sections || !out.energy_losses_au) return false;
        if (num_photon == 0 || cs_rows == 0) return false;
        if (cs_cols != num_photon || el_rows != cs_rows || el_cols != num_photon) return false;

        out.num_photon_energies = num_photon;
        out.num_final_states    = cs_rows;

        // The XAS panel is optional and is dropped whole when the core data is missing or does not
        // agree with itself - half of it would draw oscillator strengths against the wrong energies.
        if (!out.core_eigenvalues_au || !out.core_osc_strengths || num_core != num_core_osc) {
            out.core_eigenvalues_au = nullptr;
            out.core_osc_strengths  = nullptr;
            num_core = 0;
        }
        out.num_core_states = num_core;

        // Likewise the elastic line, which is drawn over the photon axis.
        if (num_elastic != num_photon) {
            out.elastic_cross_sections = nullptr;
        }

        return true;
    }

    // Drops every cached array by rewinding the arena. The pointers are nulled rather than freed
    // individually: the arena owns nothing else, so a reset is sufficient and exhaustive.
    static void reset_rixs_cache(Rixs& r) {
        if (r.alloc) {
            md_vm_arena_reset(r.alloc);
        }
        r.map = r.grid = r.photon_ev = nullptr;
        r.xas_x = r.xas_y = r.core_ev = r.core_f = nullptr;
        r.num_rows = r.num_cols = r.num_xas = r.num_core = 0;
        r.valid = false;
    }

    // Evaluates the configured broadening kernel at x for a set of sticks.
    // Both kernels use exactly the same normalization as lorentzian_xps / gaussian_xps in
    // VeloxChem's spectrumplot.py, so the resulting intensities are directly comparable.
    static inline double rixs_broaden(broadening_mode_t mode, double x, const double* peaks_x, const double* peaks_y, size_t num_peaks, double fwhm) {
        return (mode == BROADENING_MODE_LORENTZIAN)
            ? lorentzian_xps(x, peaks_x, peaks_y, num_peaks, fwhm)
            : gaussian_xps  (x, peaks_x, peaks_y, num_peaks, fwhm);
    }

    // Rebuilds the cached RIXS map, XAS curve and axis extents from 'in' using the settings in 'r'.
    // This is the expensive part; it is O(P * N * F) and is therefore guarded by a hash in draw_rixs_map().
    void compute_rixs_map(Rixs& r, const RixsMapInput& in) {
        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        // Every array below is rebuilt from scratch, so nothing from the previous map needs to
        // survive. Rewinding here keeps the footprint at exactly one map regardless of how many
        // times the settings are changed. Also nulls the pointers, which matters for the early
        // returns further down: they leave the cache invalid but must not leave it dangling.
        reset_rixs_cache(r);

        md_allocator_i* alloc = r.alloc;
        MD_ASSERT(alloc);

        const size_t P = in.num_photon_energies;
        const size_t F = in.num_final_states;

        const double* x_src = in.energy_losses_au;
        if (P == 0 || F == 0 || !in.photon_energies_au || !in.cross_sections || !x_src) return;

        // Include the elastic line as one extra stick per photon energy, if requested and available.
        const bool use_elastic = r.plot_elastic_line && in.elastic_cross_sections != nullptr;
        const size_t num_sticks = F + (use_elastic ? 1 : 0);

        // ---- Sort the incoming photon energies ascending (np.argsort) -------------------------------
        int* sort_idx = md_temp_alloc_array(temp, int, P);
        for (size_t i = 0; i < P; ++i) sort_idx[i] = (int)i;
        std::sort(sort_idx, sort_idx + P, [&](int a, int b) {
            return in.photon_energies_au[a] < in.photon_energies_au[b];
        });

        md_array_resize(r.photon_ev, P, alloc);
        for (size_t i = 0; i < P; ++i) {
            r.photon_ev[i] = in.photon_energies_au[sort_idx[i]] * HARTREE_TO_EV;
        }

        // ---- Gather the sticks per photon energy, column-contiguous ---------------------------------
        // stick_x[i * num_sticks + s], stick_y[i * num_sticks + s]
        double* stick_x = md_temp_alloc_array(temp, double, P * num_sticks);
        double* stick_y = md_temp_alloc_array(temp, double, P * num_sticks);

        double x_lo =  DBL_MAX;
        double x_hi = -DBL_MAX;

        for (size_t i = 0; i < P; ++i) {
            const size_t p = (size_t)sort_idx[i];
            double* sx = stick_x + i * num_sticks;
            double* sy = stick_y + i * num_sticks;

            for (size_t f = 0; f < F; ++f) {
                const double x = x_src[f * P + p] * HARTREE_TO_EV;
                sx[f] = x;
                sy[f] = in.cross_sections[f * P + p];
                x_lo = MD_MIN(x_lo, x);
                x_hi = MD_MAX(x_hi, x);
            }

            if (use_elastic) {
                // Elastic scattering: zero energy loss.
                const double x = 0.0;
                sx[F] = x;
                sy[F] = in.elastic_cross_sections[p];
                x_lo = MD_MIN(x_lo, x);
                x_hi = MD_MAX(x_hi, x);
            }
        }

        if (!(x_lo <= x_hi)) return;  // catches NaN / empty input

        // ---- Build the grid along the energy loss axis -------------------------------------------
        // The python version uses np.arange(val_min, val_max, xstep). Here the sample count is derived
        // from the same step but clamped, so that a small step cannot blow up the allocation.
        const double pad = 0.5;
        const double val_min = x_lo - pad;
        const double val_max = x_hi + pad;

        const double step = (r.x_step_ev > 0.0) ? r.x_step_ev : 0.01;
        size_t N = (size_t)((val_max - val_min) / step + 0.5);
        N = MD_CLAMP(N, (size_t)2, (size_t)4096);

        md_array_resize(r.grid, N, alloc);
        for (size_t j = 0; j < N; ++j) {
            r.grid[j] = val_min + (val_max - val_min) * (double)j / (double)(N - 1);
        }

        // ---- Broaden each photon-energy column into a row of the map ---------------------------------
        md_array_resize(r.map, P * N, alloc);

        double map_lo =  DBL_MAX;
        double map_hi = -DBL_MAX;

        for (size_t i = 0; i < P; ++i) {
            const double* sx = stick_x + i * num_sticks;
            const double* sy = stick_y + i * num_sticks;

            // Reverse the row order so that the lowest photon energy ends up at the bottom of the
            // heatmap, matching matplotlib's origin='lower'.
            double* row = r.map + (P - 1 - i) * N;

            for (size_t j = 0; j < N; ++j) {
                const double y = rixs_broaden(r.broadening_mode, r.grid[j], sx, sy, num_sticks, r.broadening_fwhm_ev);
                row[j] = y;
                map_lo = MD_MIN(map_lo, y);
                map_hi = MD_MAX(map_hi, y);
            }
        }

        // ---- Normalize to the global maximum ---------------------------------------------------------
        if (r.normalize && map_hi > 0.0) {
            const double scl = 1.0 / map_hi;
            for (size_t k = 0; k < P * N; ++k) {
                r.map[k] *= scl;
            }
            map_lo *= scl;
            map_hi  = 1.0;
        }

        r.num_rows = P;
        r.num_cols = N;
        r.map_min  = map_lo;
        r.map_max  = map_hi;
        r.x_min_ev = val_min;
        r.x_max_ev = val_max;

        // ---- Photon energy axis extent (pixel edges, half a step beyond the outermost centers) -------
        double photon_step = 1.0;
        if (P > 1) {
            photon_step = (r.photon_ev[P - 1] - r.photon_ev[0]) / (double)(P - 1);

            // The heatmap assumes a uniform spacing; warn instead of asserting like the python version.
            for (size_t i = 1; i < P; ++i) {
                const double d = r.photon_ev[i] - r.photon_ev[i - 1];
                if (fabs(d - photon_step) > 1.0e-5 * MD_MAX(1.0, fabs(photon_step))) {
                    if (!r.warned_non_uniform) {
                        r.warned_non_uniform = true;
                        MD_LOG_INFO("RIXS map: incoming photon energies are not uniformly spaced, the map will be distorted.");
                    }
                    break;
                }
            }
        }
        r.y_min_ev = r.photon_ev[0]     - 0.5 * photon_step;
        r.y_max_ev = r.photon_ev[P - 1] + 0.5 * photon_step;

        // ---- XAS side panel --------------------------------------------------------------------------
        const size_t C = in.num_core_states;
        if (C > 0 && in.core_eigenvalues_au && in.core_osc_strengths) {
            int* xas_sort = md_temp_alloc_array(temp, int, C);
            for (size_t i = 0; i < C; ++i) xas_sort[i] = (int)i;
            std::sort(xas_sort, xas_sort + C, [&](int a, int b) {
                return in.core_eigenvalues_au[a] < in.core_eigenvalues_au[b];
            });

            md_array_resize(r.core_ev, C, alloc);
            md_array_resize(r.core_f,  C, alloc);
            for (size_t i = 0; i < C; ++i) {
                r.core_ev[i] = in.core_eigenvalues_au[xas_sort[i]] * HARTREE_TO_EV;
                r.core_f[i]  = in.core_osc_strengths[xas_sort[i]];
            }
            r.num_core = C;

            const double xas_min = r.core_ev[0]     - 1.0;
            const double xas_max_e = r.core_ev[C - 1] + 1.0;

            size_t M = (size_t)((xas_max_e - xas_min) / 0.01 + 0.5);
            M = MD_CLAMP(M, (size_t)2, (size_t)4096);

            md_array_resize(r.xas_x, M, alloc);
            md_array_resize(r.xas_y, M, alloc);

            // The XAS is broadened with the core-hole lifetime (gamma), not with the map broadening.
            const double xas_fwhm = (in.gamma_fwhm_ev > 0.0) ? in.gamma_fwhm_ev : r.broadening_fwhm_ev;

            double y_hi = 0.0;
            for (size_t j = 0; j < M; ++j) {
                const double x = xas_min + (xas_max_e - xas_min) * (double)j / (double)(M - 1);
                const double y = rixs_broaden(r.broadening_mode, x, r.core_ev, r.core_f, C, xas_fwhm);
                r.xas_x[j] = x;
                r.xas_y[j] = y;
                y_hi = MD_MAX(y_hi, y);
            }

            // The sticks are drawn on the same axis, so make sure they are not clipped either.
            for (size_t i = 0; i < C; ++i) {
                y_hi = MD_MAX(y_hi, r.core_f[i]);
            }
            r.num_xas = M;
            r.xas_max = y_hi;
        }

        r.valid = true;
    }

    // Draws the XAS side panel content: a filled broadened curve plotted as x = f(y), plus the
    // core-excitation sticks as horizontal bars. Only valid inside a plot context.
    void plot_xas_sidebar(const Rixs& r, int& hovered, ImVec4 curve_col, ImVec4 stick_col) {
        ImDrawList& draw_list = *ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        MD_DEFER { ImPlot::PopPlotClipRect(); };

        // Fill between x = 0 and the curve. ImPlot::PlotShaded only shades along y, so the horizontal
        // equivalent of matplotlib's fill_betweenx is drawn manually, one quad per segment.
        if (r.num_xas > 1) {
            const ImU32 fill_col = ImGui::ColorConvertFloat4ToU32(ImVec4(curve_col.x, curve_col.y, curve_col.z, curve_col.w * 0.25f));
            for (size_t j = 0; j + 1 < r.num_xas; ++j) {
                const ImVec2 a = ImPlot::PlotToPixels(0.0,          r.xas_x[j],     IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 b = ImPlot::PlotToPixels(r.xas_y[j],   r.xas_x[j],     IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 c = ImPlot::PlotToPixels(r.xas_y[j+1], r.xas_x[j + 1], IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 d = ImPlot::PlotToPixels(0.0,          r.xas_x[j + 1], IMPLOT_AUTO, IMPLOT_AUTO);
                draw_list.AddQuadFilled(a, b, c, d, fill_col);
            }
        }

        // Sticks (matplotlib barh), drawn as horizontal bars of a fixed pixel thickness.
        const float bar_thickness_in_pixels = 2.0f;
        hovered = -1;

        if (ImPlot::IsPlotHovered() && r.num_core > 0) {
            const ImVec2 mouse_pos = ImGui::GetMousePos();
            const float  tolerance = 7.0f;
            double min_dist = DBL_MAX;

            for (size_t i = 0; i < r.num_core; ++i) {
                const ImVec2 p0 = ImPlot::PlotToPixels(0.0,          r.core_ev[i], IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 p1 = ImPlot::PlotToPixels(r.core_f[i],  r.core_ev[i], IMPLOT_AUTO, IMPLOT_AUTO);
                const ImVec2 lo = ImMin(p0, p1);
                const ImVec2 hi = ImMax(p0, p1);

                const ImVec2 dx = mouse_pos - ImClamp(mouse_pos, lo, hi);
                const double d2 = dx.x * dx.x + dx.y * dx.y;
                if (d2 < (tolerance * tolerance) && d2 < min_dist) {
                    hovered = (int)i;
                    min_dist = d2;
                }
            }
        }

        for (size_t i = 0; i < r.num_core; ++i) {
            const ImVec4 col = ((int)i == hovered)            ? COLOR_PEAK_HOVER
                             : ((int)i == r.selected_core)   ? COLOR_PEAK_SELECTED
                                                             : stick_col;
            ImVec2 p0 = ImPlot::PlotToPixels(0.0,         r.core_ev[i], IMPLOT_AUTO, IMPLOT_AUTO);
            ImVec2 p1 = ImPlot::PlotToPixels(r.core_f[i], r.core_ev[i], IMPLOT_AUTO, IMPLOT_AUTO);
            p0.y -= bar_thickness_in_pixels * 0.5f;
            p1.y += bar_thickness_in_pixels * 0.5f;
            draw_list.AddRectFilled(ImMin(p0, p1), ImMax(p0, p1), ImGui::ColorConvertFloat4ToU32(col));
        }
    }

    // The settings for the whole RIXS section. Drawn once, above both plots: the spectrum is a single
    // row of the map, so they share the broadening and the elastic line and must not drift apart.
    void draw_rixs_settings(const RixsMapInput& in, Rixs& r) {
        const float avail_width = ImGui::GetContentRegionAvail().x;
        ImGui::PushItemWidth(MD_MIN(avail_width, 200.0f));
        MD_DEFER { ImGui::PopItemWidth(); };

        static const double broadening_min = 0.01;
        static const double broadening_max = 2.0;
        ImGui::SliderScalar((const char*)u8"Broadening FWHM (eV)", ImGuiDataType_Double, &r.broadening_fwhm_ev, &broadening_min, &broadening_max);
        ImGui::Combo("Broadening mode", (int*)(&r.broadening_mode), broadening_mode_str, BROADENING_MODE_COUNT);

        if (in.elastic_cross_sections) {
            ImGui::Checkbox("Elastic line", &r.plot_elastic_line);
            ImGui::SameLine();
        }
        // Only affects the map, the spectrum is drawn on its own auto fitted axes.
        ImGui::Checkbox("Normalize map", &r.normalize);
    }

    // The RIXS section: one set of settings, then the spectrum and the map, each in its own tree node
    // so either can be collapsed away.
    void draw_rixs_section(const RixsMapInput& in, Rixs& r) {
        // The enclosing tree node uses NoTreePushOnOpen, so without this everything drawn here would
        // share the ID scope of the parent window - and collide with the identically labelled
        // broadening widgets of the XPS section, which sits under the same kind of node.
        ImGui::PushID("rixs");
        MD_DEFER { ImGui::PopID(); };

        const ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_NoTreePushOnOpen;

        draw_rixs_settings(in, r);

        ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
        if (ImGui::TreeNodeEx("Spectrum", flags)) {
            draw_rixs_spectrum(in, r);
        }

        ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
        if (ImGui::TreeNodeEx("Map", flags)) {
            draw_rixs_map(in, r);
        }
    }

    // Standalone RIXS spectrum widget: the energy loss spectrum for one incoming photon energy, which
    // is exactly one row of the RIXS map. The photon energy is picked from a combo listing every
    // incoming energy of the calculation, the broadening is shared with the map.
    // Call from inside an ImGui window.
    void draw_rixs_spectrum(const RixsMapInput& in, Rixs& r, ImVec2 size = ImVec2(-1.0f, 350.0f)) {
        const size_t P = in.num_photon_energies;
        const size_t F = in.num_final_states;
        if (P == 0 || F == 0 || !in.photon_energies_au || !in.cross_sections || !in.energy_losses_au) return;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        constexpr int num_broadened_samples = 2048;

        // Photon energies ascending, so the combo lists them in order and the indices agree with the
        // row order of the map.
        int* sort_idx = md_temp_alloc_array(temp, int, P);
        for (size_t i = 0; i < P; ++i) sort_idx[i] = (int)i;
        std::sort(sort_idx, sort_idx + P, [&](int a, int b) {
            return in.photon_energies_au[a] < in.photon_energies_au[b];
        });

        double* photon_ev = md_temp_alloc_array(temp, double, P);
        for (size_t i = 0; i < P; ++i) {
            photon_ev[i] = in.photon_energies_au[sort_idx[i]] * HARTREE_TO_EV;
        }

        // The selection is kept across data sets, so it has to be validated against the current one.
        r.selected_photon = MD_CLAMP(r.selected_photon, 0, (int)P - 1);

        // ---- Which row of the map to show ------------------------------------------------------------
        // Not a setting of the section: it selects the data, the broadening settings above apply to
        // both plots and are drawn once by draw_rixs_settings().
        {
            const float avail_width = ImGui::GetContentRegionAvail().x;
            ImGui::PushItemWidth(MD_MIN(avail_width, 200.0f));
            MD_DEFER { ImGui::PopItemWidth(); };

            char preview[32];
            snprintf(preview, sizeof(preview), "%.4f eV", photon_ev[r.selected_photon]);
            if (ImGui::BeginCombo("Photon energy", preview)) {
                for (int i = 0; i < (int)P; ++i) {
                    char item[32];
                    snprintf(item, sizeof(item), "%.4f eV", photon_ev[i]);
                    if (ImGui::Selectable(item, i == r.selected_photon)) {
                        r.selected_photon = i;
                        r.hovered_final   = -1;
                        r.selected_final  = -1;
                    }
                    if (i == r.selected_photon) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }
        }

        // ---- Sticks for the selected photon energy ---------------------------------------------------
        // Column p of the [F][P] arrays, gathered into the tightly packed arrays the broadening kernels
        // and plot_peaks expect. A handful of doubles, from the frame temp arena.
        const size_t p = (size_t)sort_idx[r.selected_photon];
        const bool   use_elastic = r.plot_elastic_line && in.elastic_cross_sections != nullptr;
        const size_t num_peaks   = F + (use_elastic ? 1 : 0);

        double* peaks_x = md_temp_alloc_array(temp, double, num_peaks);
        double* peaks_y = md_temp_alloc_array(temp, double, num_peaks);

        double x_lo =  DBL_MAX;
        double x_hi = -DBL_MAX;
        for (size_t f = 0; f < F; ++f) {
            peaks_x[f] = in.energy_losses_au[f * P + p] * HARTREE_TO_EV;
            peaks_y[f] = in.cross_sections[f * P + p];
            x_lo = MD_MIN(x_lo, peaks_x[f]);
            x_hi = MD_MAX(x_hi, peaks_x[f]);
        }
        if (use_elastic) {
            // Elastic scattering: zero energy loss.
            peaks_x[F] = 0.0;
            peaks_y[F] = in.elastic_cross_sections[p];
            x_lo = MD_MIN(x_lo, peaks_x[F]);
            x_hi = MD_MAX(x_hi, peaks_x[F]);
        }

        if (!(x_lo <= x_hi)) return;  // catches NaN / empty input

        // Same padding as the map, so both cover the same range for a given photon energy.
        const double pad = 0.5;
        const double x_min = x_lo - pad;
        const double x_max = x_hi + pad;

        // ---- Refit whenever the plotted data changes -------------------------------------------------
        uint64_t hash = md_hash64(&r.selected_photon,   sizeof(r.selected_photon),   0);
        hash = md_hash64(&r.plot_elastic_line,  sizeof(r.plot_elastic_line),  hash);
        hash = md_hash64(&r.broadening_mode,    sizeof(r.broadening_mode),    hash);
        hash = md_hash64(&r.broadening_fwhm_ev, sizeof(r.broadening_fwhm_ev), hash);
        hash = md_hash64(&in,                   sizeof(in),                   hash);

        if (hash != r.spectrum_hash) {
            r.spectrum_hash = hash;
            ImPlot::SetNextAxesToFit();
        }

        if (ImPlot::BeginPlot("RIXS spectrum", size)) {
            // The energy-loss convention draws the axis reversed (high loss on the left), as in the map.
            ImPlot::SetupAxis(ImAxis_X1, "Energy loss [eV]", ImPlotAxisFlags_Invert);
            ImPlot::SetupAxis(ImAxis_Y1, "Intensity [a.u.]");
            ImPlot::SetupAxis(ImAxis_Y2, "Cross section [a.u.]", ImPlotAxisFlags_AuxDefault);
            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
            ImPlot::SetupFinish();

            BroadenedCurve curve = {
                peaks_x, peaks_y, num_peaks, x_min, x_max,
                r.broadening_fwhm_ev, num_broadened_samples, r.broadening_mode,
            };
            ImPlot::SetAxis(ImAxis_Y1);
            ImPlot::PlotLineG("Broadened Spectrum", broadened_curve_getter, &curve, num_broadened_samples);

            ImPlot::SetAxis(ImAxis_Y2);
            plot_peaks("Cross Section", peaks_x, peaks_y, num_peaks, r.selected_final, r.hovered_final);

            if (ImPlot::IsPlotHovered() && r.hovered_final != -1) {
                const int i = r.hovered_final;
                if (ImGui::BeginTooltip()) {
                    if (use_elastic && (size_t)i == F) {
                        ImGui::TextUnformatted("Elastic line");
                    } else {
                        ImGui::Text("Final state %i", i + 1);
                    }
                    ImGui::Text("Energy loss: %.4f eV", peaks_x[i]);
                    ImGui::Text("Cross section: %.4g", peaks_y[i]);
                    ImGui::EndTooltip();
                }
            }

            ImPlot::EndPlot();
        }
    }

    // RIXS map widget: 2D map, XAS side panel and colorbar. The settings it reads are drawn by
    // draw_rixs_settings().
    // Call from inside an ImGui window. 'size' is the total size of the plot row; pass -1 for the
    // width to use the available content region.
    void draw_rixs_map(const RixsMapInput& in, Rixs& r, ImVec2 size = ImVec2(-1.0f, 350.0f)) {
        // ---- Rebuild the cache only when something actually changed -----------------------------------
        uint64_t hash = md_hash64(&r.plot_elastic_line, sizeof(r.plot_elastic_line), 0);
        hash = md_hash64(&r.normalize,          sizeof(r.normalize),          hash);
        hash = md_hash64(&r.broadening_mode,    sizeof(r.broadening_mode),    hash);
        hash = md_hash64(&r.broadening_fwhm_ev, sizeof(r.broadening_fwhm_ev), hash);
        hash = md_hash64(&r.x_step_ev,          sizeof(r.x_step_ev),          hash);
        hash = md_hash64(&in,                   sizeof(in),                   hash);

        const bool refit = (hash != r.hash);
        if (refit) {
            r.hash = hash;
            compute_rixs_map(r, in);
        }

        if (!r.valid) {
            ImGui::TextUnformatted("No RIXS map could be computed from the loaded data.");
            return;
        }

        // ---- Layout: map | XAS | colorbar --------------------------------------------------------------
        const float spacing       = ImGui::GetStyle().ItemSpacing.x;
        const float colorbar_w    = 90.0f;
        const float height        = (size.y > 0.0f) ? size.y : 350.0f;
        const float total_w       = ((size.x > 0.0f) ? size.x : ImGui::GetContentRegionAvail().x) - colorbar_w - 2.0f * spacing;
        const bool  has_xas       = r.num_xas > 1;
        // Same 4 : 1.25 width ratio as the matplotlib gridspec.
        const float map_w         = has_xas ? total_w * (4.0f / 5.25f) : total_w;
        const float xas_w         = total_w - map_w;

        const char* x_label = "Energy loss [eV]";

        // The energy-loss convention draws the axis reversed (high loss on the left).
        const ImPlotAxisFlags x_flags = ImPlotAxisFlags_Invert;

        // Only the heatmap and the colorbar are drawn with the map colormap. The XAS panel is a
        // regular line plot and picks up the default item colors, like every other spectrum.
        ImPlot::PushColormap(r.colormap);

        // PlotHeatmap only rasterizes the cells it was given, so anything outside the data would fall
        // through to the window background. Painting the canvas in the color a zero sample maps to
        // makes the map read as an infinite field of zeros instead of a floating tile.
        const double map_range = MD_MAX(r.map_max - r.map_min, DBL_EPSILON);
        const float  zero_t    = ImClamp((float)((0.0 - r.map_min) / map_range), 0.0f, 1.0f);
        ImPlot::PushStyleColor(ImPlotCol_PlotBg, ImPlot::SampleColormap(zero_t, r.colormap));

        ImPlot::SetNextAxisLinks(ImAxis_Y1, &r.y_link_min, &r.y_link_max);
        if (ImPlot::BeginPlot("RIXS map", ImVec2(map_w, height), ImPlotFlags_NoLegend)) {
            ImPlot::SetupAxis(ImAxis_X1, x_label, x_flags);
            ImPlot::SetupAxis(ImAxis_Y1, "Photon energy [eV]");
            ImPlot::SetupAxesLimits(r.x_min_ev, r.x_max_ev, r.y_min_ev, r.y_max_ev, refit ? ImPlotCond_Always : ImPlotCond_Once);

            // Keep panning and zooming inside the data. Combined with the background above, there is
            // no way to end up looking at a region the map says nothing about.
            ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, r.x_min_ev, r.x_max_ev);
            ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, r.y_min_ev, r.y_max_ev);
            ImPlot::SetupFinish();

            ImPlot::PlotHeatmap("##rixs_map", r.map, (int)r.num_rows, (int)r.num_cols,
                                r.map_min, r.map_max, nullptr,
                                ImPlotPoint(r.x_min_ev, r.y_min_ev),
                                ImPlotPoint(r.x_max_ev, r.y_max_ev));

            // The row currently shown in the spectrum. White reads against the whole colormap, which
            // a colormap color would not.
            if (r.selected_photon >= 0 && (size_t)r.selected_photon < r.num_rows) {
                const double y = r.photon_ev[r.selected_photon];
                ImPlot::SetNextLineStyle(ImVec4(1.0f, 1.0f, 1.0f, 0.85f), 1.5f);
                ImPlot::PlotInfLines("##selected_photon", &y, 1, ImPlotInfLinesFlags_Horizontal);
            }

            // Clicking anywhere in the map moves the spectrum to the closest incoming photon energy.
            // Released rather than clicked, and dragging excluded, so that panning the plot does not
            // change the selection.
            if (ImPlot::IsPlotHovered() && ImGui::IsMouseReleased(ImGuiMouseButton_Left) &&
                !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) && r.num_rows > 0) {
                const double y = ImPlot::GetPlotMousePos().y;

                int    closest   = 0;
                double closest_d = DBL_MAX;
                for (size_t i = 0; i < r.num_rows; ++i) {
                    const double d = fabs(r.photon_ev[i] - y);
                    if (d < closest_d) {
                        closest_d = d;
                        closest   = (int)i;
                    }
                }

                if (closest != r.selected_photon) {
                    r.selected_photon = closest;
                    r.hovered_final   = -1;
                    r.selected_final  = -1;
                }
            }

            if (ImPlot::IsPlotHovered()) {
                const ImPlotPoint p = ImPlot::GetPlotMousePos();

                // Map the cursor back onto a cell to read out the intensity.
                const double u = (p.x - r.x_min_ev) / (r.x_max_ev - r.x_min_ev);
                const double v = (p.y - r.y_min_ev) / (r.y_max_ev - r.y_min_ev);
                if (u >= 0.0 && u < 1.0 && v >= 0.0 && v < 1.0) {
                    const size_t col = (size_t)(u * (double)r.num_cols);
                    const size_t row = r.num_rows - 1 - (size_t)(v * (double)r.num_rows);
                    ImGui::BeginTooltip();
                    ImGui::Text("Energy loss: %.3f eV", p.x);
                    ImGui::Text("Photon energy: %.3f eV", p.y);
                    ImGui::Text("Intensity: %.4g", r.map[row * r.num_cols + col]);
                    ImGui::EndTooltip();
                }
            }
            ImPlot::EndPlot();
        }
        ImPlot::PopStyleColor();
        ImPlot::PopColormap();

        if (has_xas) {
            ImGui::SameLine();
            ImPlot::SetNextAxisLinks(ImAxis_Y1, &r.y_link_min, &r.y_link_max);
            const ImPlotFlags xas_flags = ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText;
            if (ImPlot::BeginPlot("XAS", ImVec2(xas_w, height), xas_flags)) {
                ImPlot::SetupAxis(ImAxis_X1, "Intensity",  ImPlotAxisFlags_NoGridLines);
                // The photon energy axis is shared with the map, so it is not labelled again here.
                ImPlot::SetupAxis(ImAxis_Y1, nullptr, ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoLabel);
                ImPlot::SetupAxesLimits(0.0, MD_MAX(r.xas_max, 1.0e-12) * 1.1, r.y_min_ev, r.y_max_ev, refit ? ImPlotCond_Always : ImPlotCond_Once);
                // The photon energy axis is linked to the map, so it carries the same constraint.
                ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, r.y_min_ev, r.y_max_ev);
                ImPlot::SetupFinish();

                // The first two colors of the current colormap, as two ordinary plot items would get.
                // NOTE: GetColormapColor, not NextColormapColor: the latter advances a counter that
                // lives in the plot and is never reset per frame, so calling it for items that are not
                // registered with ImPlot walks through the colormap one step per frame (flickering).
                const ImVec4 curve_col = ImPlot::GetColormapColor(0);
                const ImVec4 stick_col = ImPlot::GetColormapColor(1);

                plot_xas_sidebar(r, r.hovered_core, curve_col, stick_col);

                ImPlot::SetNextLineStyle(curve_col);
                ImPlot::PlotLine("##xas", r.xas_y, r.xas_x, (int)r.num_xas);

                if (r.hovered_core != -1) {
                    ImGui::BeginTooltip();
                    ImGui::Text("Core state %i", r.hovered_core + 1);
                    ImGui::Text("Energy: %.3f eV", r.core_ev[r.hovered_core]);
                    ImGui::Text("Osc. strength: %.4g", r.core_f[r.hovered_core]);
                    ImGui::EndTooltip();
                }

                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left) && ImPlot::IsPlotHovered()) {
                    r.selected_core = (r.hovered_core == r.selected_core) ? -1 : r.hovered_core;
                }

                ImPlot::EndPlot();
            }
        }

        ImGui::SameLine();
        ImPlot::PushColormap(r.colormap);
        ImPlot::ColormapScale("Intensity [arb. u.]", r.map_min, r.map_max, ImVec2(colorbar_w, height),
                              r.normalize ? "%.2f" : "%.1e");
        ImPlot::PopColormap();
    }

    // =================================================================================================
    // XPS
    //
    // ImGui/ImPlot port of VeloxChem's spectrumplot.plot_xps_spectrum(). Like the python version this
    // draws one element at a time, picked from a combo box listing only the elements that actually
    // have core hole entries.
    //
    // Within that element every peak is one ATOM, so the plot doubles as a way to find the atom a
    // binding energy belongs to: hovering a peak highlights its atom in the viewport and clicking
    // adds it to the selection. That is also why there is no CPK coloring any more - the colour used
    // to encode the element, which the combo box now carries, and per atom CPK would encode nothing
    // at all since every peak on screen shares one element.
    //
    // Two ImPlot items, matching the other spectra in the Response window: "Broadened Spectrum" on Y1
    // and "Peaks" on Y2, each with its own legend entry and show/hide state.
    // =================================================================================================

    // Per-element input to the broadening getter. Lives on the stack for the duration of one
    // PlotLineG call, which consumes it immediately.
    // A broadened curve sampled on the fly by ImPlot, shared by the XPS and RIXS spectra.
    struct BroadenedCurve {
        const double*     peaks_x;
        const double*     peaks_y;
        size_t            num_peaks;
        double            x_min;
        double            x_max;
        double            fwhm;
        int               num_samples;
        broadening_mode_t mode;
    };

    static ImPlotPoint broadened_curve_getter(int idx, void* user_data) {
        const BroadenedCurve* c = (const BroadenedCurve*)user_data;
        const double x = c->x_min + (c->x_max - c->x_min) * idx / (double)(c->num_samples - 1);
        const double y = (c->mode == BROADENING_MODE_LORENTZIAN)
            ? lorentzian_xps(x, c->peaks_x, c->peaks_y, c->num_peaks, c->fwhm)
            : gaussian_xps  (x, c->peaks_x, c->peaks_y, c->num_peaks, c->fwhm);
        return ImPlotPoint{ x, y };
    }

    // ---- XPS ------------------------------------------------------------------------------------
    //
    // The core-hole table as it sits in the system's attribute table: six sibling columns over one
    // {C} index space. A record is not one attribute - a value has ONE type - so what the file
    // carries as an array of structs is read back as six pointers, which is also why the plot below
    // hands ImPlot contiguous arrays instead of a base pointer and a struct stride.
    //
    // The pointers are into the table's own storage and are invalidated by the next attribute create
    // or remove; nothing here holds one past the frame it gathered it in.
    struct XpsTable {
        size_t         count             = 0;
        const double*  ionization_energy = nullptr;  // [C], eV
        const double*  contribution      = nullptr;  // [C], the atom's share of the core MO
        const int32_t* atom_index        = nullptr;  // [C], into the QM ATOM DOMAIN; -1 when unassigned
        const int32_t* mo_index          = nullptr;  // [C], into the molecular orbitals
        const uint8_t* element           = nullptr;  // [C], atomic number
        const uint8_t* is_delocalized    = nullptr;  // [C], one byte per entry
    };

    // False when the system carries no XPS at all. The ionization energies and the elements are what
    // everything below indexes and groups by, so a table missing either is not one that can be
    // drawn; the other four are optional and are checked where they are read.
    static bool xps_table_gather(XpsTable* out, const md_system_t& sys) {
        size_t num_energy = 0;
        size_t num_element = 0;
        const double*  energy  = (const double*) attribute_column(&num_energy,  sys, MD_STR_LIT("vlx/xps/ionization_energy"), MD_ATTRIBUTE_TYPE_F64);
        const uint8_t* element = (const uint8_t*)attribute_column(&num_element, sys, MD_STR_LIT("vlx/xps/element"),           MD_ATTRIBUTE_TYPE_U8);
        if (!energy || !element || num_energy == 0 || num_element != num_energy) {
            return false;
        }

        // A sibling whose length disagrees with the index space is dropped rather than trusted: it
        // would be read at an index that means something else in it.
        auto column_of = [&](md_str_t path, md_attribute_type_t type) -> const void* {
            size_t count = 0;
            const void* data = attribute_column(&count, sys, path, type);
            return (data && count == num_energy) ? data : nullptr;
        };

        out->count             = num_energy;
        out->ionization_energy = energy;
        out->element           = element;
        out->contribution      = (const double*) column_of(MD_STR_LIT("vlx/xps/contribution"),   MD_ATTRIBUTE_TYPE_F64);
        out->atom_index        = (const int32_t*)column_of(MD_STR_LIT("vlx/xps/atom_index"),     MD_ATTRIBUTE_TYPE_I32);
        out->mo_index          = (const int32_t*)column_of(MD_STR_LIT("vlx/xps/mo_index"),       MD_ATTRIBUTE_TYPE_I32);
        out->is_delocalized    = (const uint8_t*)column_of(MD_STR_LIT("vlx/xps/is_delocalized"), MD_ATTRIBUTE_TYPE_U8);
        return true;
    }

    // One element's states: an OFFSET and a LENGTH into the columns above, not a copy of them.
    struct XpsGroup {
        uint8_t element = 0;
        size_t  offset  = 0;
        size_t  count   = 0;
    };

    // The per element runs, DERIVED rather than published. The entries are laid out as contiguous
    // runs of equal element - the reader sorts them by (element, ionization energy) and md_vlx.h
    // says so - so one scan of the element column reproduces the grouping exactly, and publishing
    // the runs as well would be a second representation of one fact with nothing keeping the two in
    // agreement.
    //
    // Returns the total number of runs and writes at most cap, so pass cap 0 to count.
    static size_t xps_groups_gather(XpsGroup out[], size_t cap, const XpsTable& tbl) {
        size_t num = 0;
        size_t i = 0;
        while (i < tbl.count) {
            const uint8_t element = tbl.element[i];
            size_t j = i;
            while (j < tbl.count && tbl.element[j] == element) ++j;
            if (out && num < cap) {
                out[num] = XpsGroup{ element, i, j - i };
            }
            num += 1;
            i = j;
        }
        return num;
    }

    // Draws the stick spectrum of one element as an ImPlot item under 'label'.
    // 'hovered' is an index into the flat entry array and 'base' is the offset of this group within
    // it, so a hover survives being compared against entries from any group.
    // 'stick_selected' is per entry and comes from the application selection mask - the plot does not
    // own a selection of its own.
    // Returns false if the item is hidden through the legend, in which case nothing was drawn.
    static bool plot_xps_sticks(const char* label, const double* stick_x, const double* stick_y,
                                const bool* stick_selected, const uint8_t* stick_delocalized,
                                size_t count, size_t base, int& hovered) {
        const float bar_width_in_pixels    = 2.0f;
        const float point_radius_in_pixels = 3.0f;

        if (!ImPlot::BeginItem(label)) {
            // Hidden via the legend. Clear any hover that belonged to this element so a stale
            // tooltip does not survive the toggle.
            if (hovered >= (int)base && hovered < (int)(base + count)) hovered = -1;
            return false;
        }
        MD_DEFER { ImPlot::EndItem(); };

        ImPlot::PushPlotClipRect();
        MD_DEFER { ImPlot::PopPlotClipRect(); };

        ImDrawList& draw_list = *ImPlot::GetPlotDrawList();
        const ImPlotNextItemData& s = ImPlot::GetItemData();
        const ImVec4 item_color = s.Colors[ImPlotCol_Fill];

        const bool  plot_hovered = ImPlot::IsPlotHovered();
        const ImVec2 mouse_pos   = ImGui::GetMousePos();
        const double min_line_d2  = pow(bar_width_in_pixels * 0.5 + 7.0, 2);
        const double min_point_d2 = pow(point_radius_in_pixels + 7.0, 2);
        double min_dist = DBL_MAX;

        for (size_t i = 0; i < count; ++i) {
            const double x = stick_x[i];
            const double y = stick_y[i];

            const ImVec2 p0_raw = ImPlot::PlotToPixels(x, 0.0, IMPLOT_AUTO, IMPLOT_AUTO);
            const ImVec2 p1_raw = ImPlot::PlotToPixels(x, y,   IMPLOT_AUTO, IMPLOT_AUTO);
            const ImVec2 p_min  = ImMin(p0_raw, p1_raw);
            const ImVec2 p_max  = ImMax(p0_raw, p1_raw);

            if (plot_hovered) {
                const ImVec2 line_dx = mouse_pos - ImClamp(mouse_pos, p_min, p_max);
                const double line_d2 = line_dx.x * line_dx.x + line_dx.y * line_dx.y;
                if (line_d2 < min_line_d2 && line_d2 < min_dist) {
                    hovered  = (int)(base + i);
                    min_dist = line_d2;
                }
                const ImVec2 point_dx = mouse_pos - p1_raw;
                const double point_d2 = point_dx.x * point_dx.x + point_dx.y * point_dx.y;
                if (point_d2 < min_point_d2 && point_d2 < min_dist) {
                    hovered  = (int)(base + i);
                    min_dist = point_d2;
                }
            }

            const int    flat_idx = (int)(base + i);
            const ImVec4 color    = (flat_idx == hovered)   ? COLOR_PEAK_HOVER
                                  : (stick_selected[i])     ? COLOR_PEAK_SELECTED
                                                            : item_color;
            const ImU32  col32    = ImGui::ColorConvertFloat4ToU32(color);

            draw_list.AddRectFilled(ImVec2{ p_min.x - bar_width_in_pixels * 0.5f, p_min.y },
                                    ImVec2{ p_max.x + bar_width_in_pixels * 0.5f, p_max.y }, col32);
            draw_list.AddCircleFilled(p1_raw, point_radius_in_pixels, col32);

            // A delocalized core hole is not attributable to a single atom, so mark it.
            if (stick_delocalized && stick_delocalized[i]) {
                draw_list.AddCircle(p1_raw, point_radius_in_pixels + 2.5f, col32, 0, 1.5f);
            }

            if (ImPlot::FitThisFrame()) {
                ImPlot::FitPoint(ImPlotPoint(x, 0.0));
                ImPlot::FitPoint(ImPlotPoint(x, y));
            }
        }

        return true;
    }

    // Maps an atom index in the QM ATOM DOMAIN onto an index into the loaded system. The two differ
    // when the calculation covers only a subset of the system. Returns -1 when the entry names no
    // atom at all - a delta-SCF entry can arrive without one - or when the mapping lands outside the
    // system, in which case there is nothing to highlight.
    // 'qm' is passed rather than gathered: this is called once per peak inside the plot's own
    // loops, and the domain does not change between two peaks of the same frame.
    static int32_t xps_atom_index(const QmAtoms& qm, int32_t qm_atom_index, const ApplicationState& state) {
        if (qm_atom_index < 0) return -1;
        if ((size_t)qm_atom_index >= qm.count) return -1;
        const size_t idx = qm_to_system_atom(qm, (size_t)qm_atom_index);
        return (idx < state.mold.sys.atom.count) ? (int32_t)idx : -1;
    }

    void draw_xps_plot(ApplicationState& state, ImVec2 size = ImVec2(-1.0f, 350.0f)) {
        XpsTable tbl = {};
        if (!xps_table_gather(&tbl, state.mold.sys)) return;

        // The QM atom domain, for mapping an entry's atom onto a system atom. Gathered once here
        // and handed down; the peaks below index it, they do not each look it up.
        QmAtoms qm = {};
        if (!es_qm_atoms(&qm, state.mold.sys)) return;

        // One run per element that carries a core hole, which is at most one per element there is.
        // A fixed array rather than an allocation: the gather is a scan of a column already in the
        // table and the result is used within this frame only.
        XpsGroup groups[128];
        const size_t num_groups  = MD_MIN(xps_groups_gather(groups, MD_COUNTOF(groups), tbl), MD_COUNTOF(groups));
        const size_t num_entries = tbl.count;
        if (num_groups == 0) return;

        // The enclosing tree node uses NoTreePushOnOpen, so this section would otherwise draw into the
        // parent window's ID scope. Two collisions follow from that: ImPlot derives the plot ID from
        // Window->GetID(title), so BeginPlot("XPS") lands on the exact ID of the TreeNodeEx("XPS")
        // header - which makes ImPlot report the plot as held while the header is dragged, panning the
        // axes - and the broadening widgets below collide with the identically labelled ones in the
        // RIXS section when a file carries both.
        ImGui::PushID("xps");
        MD_DEFER { ImGui::PopID(); };

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        constexpr int num_broadened_samples = 2048;

        // Resolves the selected element, adopting the first non empty group when the selection is
        // unset or no longer present. Called again after the combo so a change lands on the same
        // frame it is made rather than one frame later.
        auto resolve_group = [&]() -> const XpsGroup* {
            for (size_t i = 0; i < num_groups; ++i) {
                if (groups[i].element == xps.element && groups[i].count > 0) return &groups[i];
            }
            for (size_t i = 0; i < num_groups; ++i) {
                if (groups[i].count > 0) {
                    xps.element = groups[i].element;
                    return &groups[i];
                }
            }
            return nullptr;
        };

        const XpsGroup* grp = resolve_group();
        if (!grp) return;

        // ---- Settings ---------------------------------------------------------------------------
        const float avail_width = ImGui::GetContentRegionAvail().x;
        ImGui::PushItemWidth(MD_MIN(avail_width, 200.0f));
        MD_DEFER { ImGui::PopItemWidth(); };

        // The group list is exactly the set of elements that were computed, so listing it offers no
        // element that would draw an empty plot.
        {
            const md_str_t cur_sym = md_util_element_symbol(grp->element);
            char preview[8];
            snprintf(preview, sizeof(preview), "%.*s", (int)cur_sym.len, cur_sym.ptr);
            if (ImGui::BeginCombo("Element", preview)) {
                for (size_t g = 0; g < num_groups; ++g) {
                    const XpsGroup* it = &groups[g];
                    if (it->count == 0) continue;
                    const md_str_t sym = md_util_element_symbol(it->element);
                    char item[32];
                    snprintf(item, sizeof(item), "%.*s (%i)", (int)sym.len, sym.ptr, (int)it->count);
                    if (ImGui::Selectable(item, it->element == xps.element)) {
                        xps.element = it->element;
                        xps.hovered = -1;
                    }
                }
                ImGui::EndCombo();
            }
            // NOTE: no tooltip on the combo. IsItemHovered() after EndCombo() reports on whatever the
            // popup left in LastItemData, not on the combo frame. The bracketed count is the number
            // of core hole states computed for that element, one per atom.
            grp = resolve_group();
            if (!grp) return;
        }

        static const double fwhm_min = 0.05;
        static const double fwhm_max = 5.0;
        ImGui::SliderScalar((const char*)u8"Broadening FWHM (eV)", ImGuiDataType_Double, &xps.broadening_fwhm_ev, &fwhm_min, &fwhm_max);
        ImGui::Combo("Broadening mode", (int*)(&xps.broadening_mode), broadening_mode_str, BROADENING_MODE_COUNT);

        ImGui::Checkbox("Sticks", &xps.show_sticks);
        ImGui::SameLine();
        ImGui::Checkbox("Invert energy axis", &xps.invert_x);
        ImGui::SameLine();
        ImGui::Checkbox("Weight by contribution", &xps.weight_by_contribution);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Off: every core hole contributes a unit peak (what VeloxChem plots).\n"
                              "On:  each peak is scaled by the atom's share of the core orbital, so a\n"
                              "     hole delocalized over N atoms sums to one instead of N.");
        }

        // ---- Per peak arrays ---------------------------------------------------------------------
        // The broadening kernels take tightly packed arrays, so the strided entry fields are gathered
        // here. A handful of doubles for one element, from the frame temp arena.
        // The energies are already contiguous over the run, so only the two DERIVED arrays are
        // built: the stick heights and the per peak selection. A handful of doubles for one element,
        // from the frame temp arena.
        const double* peaks_x        = tbl.ionization_energy + grp->offset;
        double*       peaks_y        = md_temp_alloc_array(temp, double, grp->count);
        bool*         peaks_selected = md_temp_alloc_array(temp, bool,   grp->count);
        for (size_t i = 0; i < grp->count; ++i) {
            const size_t e = grp->offset + i;
            // Weighting by contribution needs the column; without it every hole is a unit peak,
            // which is what VeloxChem plots anyway.
            peaks_y[i] = (xps.weight_by_contribution && tbl.contribution) ? tbl.contribution[e] : 1.0;
            const int32_t atom = xps_atom_index(qm, tbl.atom_index ? tbl.atom_index[e] : -1, state);
            peaks_selected[i]  = (atom >= 0) && md_bitmask_test_bit(&state.selection.selection_mask, (uint64_t)atom);
        }

        // ---- Energy range ------------------------------------------------------------------------
        // Over the selected element only. Core binding energies of different elements sit hundreds of
        // eV apart, so a range spanning all of them would compress the peaks actually on screen into
        // a single unreadable spike.
        double x_min = DBL_MAX;
        double x_max = -DBL_MAX;
        for (size_t i = 0; i < grp->count; ++i) {
            x_min = MD_MIN(x_min, peaks_x[i]);
            x_max = MD_MAX(x_max, peaks_x[i]);
        }
        const double pad = 5.0;
        x_min = MD_MAX(0.0, x_min - pad);
        x_max = x_max + pad;

        uint64_t hash = md_hash64(&xps.broadening_fwhm_ev, sizeof(xps.broadening_fwhm_ev), xps.broadening_mode);
        hash = md_hash64(&xps.weight_by_contribution, sizeof(xps.weight_by_contribution), hash);
        hash = md_hash64(&xps.element, sizeof(xps.element), hash);
        hash = md_hash64(peaks_x, grp->count * sizeof(double), hash);
        hash = md_hash64(peaks_y, grp->count * sizeof(double), hash);
        const bool refit = (hash != xps.hash);
        xps.hash = hash;

        if (refit) {
            ImPlot::SetNextAxesToFit();
        }

        const ImPlotAxisFlags x_flags = xps.invert_x ? ImPlotAxisFlags_Invert : ImPlotAxisFlags_None;

        // NoMenus: right click removes an atom from the selection here, and ImPlot would otherwise
        // open its context menu on the very same event. The ramachandran plot gives up its context
        // menu for the same reason.
        if (ImPlot::BeginPlot("XPS", size, ImPlotFlags_NoMenus)) {
            ImPlot::SetupAxis(ImAxis_X1, "Binding energy (eV)", x_flags);
            ImPlot::SetupAxis(ImAxis_Y1, "Intensity (a.u.)");
            if (xps.show_sticks) {
                ImPlot::SetupAxis(ImAxis_Y2, "Peak intensity", ImPlotAxisFlags_AuxDefault);
            }
            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
            ImPlot::SetupFinish();

            // Recomputed from scratch every frame by plot_xps_sticks. Cleared unconditionally so that
            // hiding the sticks - through the checkbox or the legend - cannot leave a stale hover
            // driving the tooltip and the highlight.
            xps.hovered = -1;

            BroadenedCurve curve = {
                peaks_x, peaks_y, grp->count, x_min, x_max,
                xps.broadening_fwhm_ev, num_broadened_samples, xps.broadening_mode,
            };
            ImPlot::SetAxis(ImAxis_Y1);
            ImPlot::PlotLineG("Broadened Spectrum", broadened_curve_getter, &curve, num_broadened_samples);

            if (xps.show_sticks) {
                ImPlot::SetAxis(ImAxis_Y2);
                plot_xps_sticks("Peaks", peaks_x, peaks_y, peaks_selected,
                                tbl.is_delocalized ? tbl.is_delocalized + grp->offset : nullptr,
                                grp->count, grp->offset, xps.hovered);
            }

            // ---- Interaction ---------------------------------------------------------------------
            const bool plot_hovered = ImPlot::IsPlotHovered();
            // An index into the FLAT columns, which is why plot_xps_sticks was handed the group's
            // offset as its base: a hover survives being compared against entries from any group.
            const int  hov = (xps.hovered >= 0 && xps.hovered < (int)num_entries) ? xps.hovered : -1;
            const int32_t hov_atom = (hov >= 0) ? xps_atom_index(qm, tbl.atom_index ? tbl.atom_index[hov] : -1, state) : -1;

            // While the cursor is inside the plot this owns the highlight mask, which is the contract
            // every other hover provider in the application follows. Cleared unconditionally, so
            // moving off a peak but staying inside the plot drops the highlight rather than leaving
            // the last atom lit up.
            if (plot_hovered) {
                md_bitmask_clear(&state.selection.highlight_mask);
                if (hov_atom >= 0) {
                    md_bitmask_or(&state.selection.highlight_mask, md_bitset_range((uint32_t)((uint64_t)hov_atom), (uint32_t)((uint64_t)hov_atom) + 1));
                }
            }

            // Same contract towards the orbital grid window: only written while the cursor is inside
            // the plot, so moving off a peak drops the row highlight exactly as it drops the atom
            // highlight. A delocalized hole names one shared core MO, so every entry of that hole
            // lights the same row - which is the truth of it.
            // The reader range checks atom_index on load but not mo_index, so it is checked here.
            if (plot_hovered && hov >= 0 && tbl.mo_index && tbl.mo_index[hov] >= 0 &&
                (size_t)tbl.mo_index[hov] < num_molecular_orbitals()) {
                xps.highlight_mo_idx = tbl.mo_index[hov];
            }

            if (plot_hovered && hov >= 0) {
                const int32_t hov_atom_index = tbl.atom_index ? tbl.atom_index[hov] : -1;
                const md_str_t sym = md_util_element_symbol(tbl.element[hov]);
                if (ImGui::BeginTooltip()) {
                    // Atom index is shown 1-based to match the labels VeloxChem prints.
                    if (hov_atom_index >= 0) {
                        ImGui::Text("%.*s%i", (int)sym.len, sym.ptr, hov_atom_index + 1);
                    } else {
                        ImGui::Text("%.*s (unassigned atom)", (int)sym.len, sym.ptr);
                    }
                    ImGui::Text("Binding energy: %.3f eV", tbl.ionization_energy[hov]);
                    if (tbl.contribution) ImGui::Text("Contribution: %.4g", tbl.contribution[hov]);
                    if (tbl.mo_index)     ImGui::Text("MO index: %i", tbl.mo_index[hov]);
                    if (tbl.is_delocalized && tbl.is_delocalized[hov]) {
                        // The hole is shared with the other atoms carrying this same mo_index; only
                        // the atom this entry names is highlighted.
                        ImGui::TextUnformatted("Delocalized core hole");
                    }
                    ImGui::EndTooltip();
                }
            }

            // Left adds the atom to the selection, right removes it - the same convention the
            // ramachandran plot and the dataset window use. Drags are excluded so that panning or
            // box selecting does not also grab whatever sat under the cursor when the drag began.
            if (plot_hovered) {
                if (hov_atom >= 0) {
                    if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Left)) {
                        if (!ImGui::IsKeyDown(ImGuiMod_Shift)) {
                            md_bitmask_clear(&state.selection.selection_mask);
                        }
                        md_bitmask_or(&state.selection.selection_mask, md_bitset_range((uint32_t)((uint64_t)hov_atom), (uint32_t)((uint64_t)hov_atom) + 1));
                    }
                    else if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Right)) {
                        md_bitmask_andnot(&state.selection.selection_mask, md_bitset_range((uint32_t)((uint64_t)hov_atom), (uint32_t)((uint64_t)hov_atom) + 1));
                    }
                }
                else if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && !ImGui::IsMouseDragPastThreshold(ImGuiMouseButton_Right)) {
                    md_bitmask_clear(&state.selection.selection_mask);
                }
            }

            ImPlot::EndPlot();
        }
    }

    void set_atom_coordinates(ApplicationState& state, const dvec3_t* atom_coords) {
        QmAtoms qm = {};
        if (!es_qm_atoms(&qm, state.mold.sys)) return;

        const dvec3_t* coords = atom_coords ? atom_coords : qm.coordinate;
        if (!coords) return;

        // Through the map. These are QM atom coordinates and the destination is the SYSTEM's state,
        // so writing them at i would put QM atom i's position on system atom i - which is only the
        // same atom when the calculation covers the whole system. Stepping an optimisation on a
        // subset used to scramble the first N atoms of the molecule.
        for (size_t i = 0; i < qm.count; i++) {
            const size_t idx = qm_to_system_atom(qm, i);
            if (idx >= state.mold.state.num_atoms) continue;
            state.mold.state.xyz[idx].x = (float)coords[i].x;
            state.mold.state.xyz[idx].y = (float)coords[i].y;
            state.mold.state.xyz[idx].z = (float)coords[i].z;
        }
		viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemStateChanged, viamd::EventPayloadType_ApplicationState, &state);
        state.mold.dirty_gpu_buffers |= MolBit_ClearVelocity | MolBit_DirtyPosition;
    }

    // The ground state dipole vector, in the e*bohr it is stored in, or zero when absent. Read as
    // stored (md_unit_none) rather than converted: the centre of charge arithmetic below is in the
    // same units as the nuclear term it is subtracted from, and a conversion here would be a
    // silent change of what that subtraction means.
    static dvec3_t dipole_ground_state(const md_system_t& sys) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, MD_STR_LIT("dipole/ground_state/vector"));
        if (!attr) return {0, 0, 0};
        double v[3] = {0, 0, 0};
        if (md_attribute_extract_f64(v, MD_COUNTOF(v), attr, md_attribute_slice_all(), md_unit_none()) != MD_COUNTOF(v)) {
            return {0, 0, 0};
        }
        return {v[0], v[1], v[2]};
    }

    // The longest vector in a {N,3} dipole attribute, or 0 when it is absent. Only the magnitude is
    // wanted - it sets an arrow scale - so the vectors themselves never have to come back.
    static double dipole_max_length(md_temp_scope_t temp, const md_system_t& sys, md_str_t path) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        if (!attr || attr->format.components != 3) {
            return 0.0;
        }
        const size_t count = md_attribute_element_count(&attr->format);
        double* v = md_temp_alloc_array(temp, double, count);
        if (!v || md_attribute_extract_f64(v, count, attr, md_attribute_slice_all(), md_unit_none()) != count) {
            return 0.0;
        }
        double max_len = 0.0;
        for (size_t i = 0; i + 2 < count; i += 3) {
            max_len = MD_MAX(max_len, dvec3_length({v[i], v[i + 1], v[i + 2]}));
        }
        return max_len;
    }

    // One string out of the table, or empty when the path is absent. Pool strings are NUL
    // terminated, so str_ptr goes straight into a "%s".
    static md_str_t attribute_str(const md_system_t& sys, md_str_t path) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        if (!attr || attr->format.type != MD_ATTRIBUTE_TYPE_STR) return {};
        return md_attribute_str(&sys.attributes, attr, 0);
    }

    // Reads one rank 1 series of doubles out of the system's attribute table.
    //
    // Hands back the stored buffer rather than extracting to float, deliberately: the vlx/ series
    // are resident F64 and this panel prints total energies to ten decimals, which a float does not
    // carry. That is the opt-in fast path - it is only valid because the attribute is resident, and
    // the null return covers every case where it is not.
    static const double* attribute_series_f64(size_t* out_count, const md_system_t& sys, md_str_t path) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        const double* values = (const double*)md_attribute_view(attr, MD_ATTRIBUTE_TYPE_F64, 1, 1);
        if (values && out_count) *out_count = attr->format.shape[0];
        return values;
    }

    // True when the table carries this path at all. "Does this file have X" is a question about the
    // table and nothing else: a block a file does not contain publishes nothing, so there is no
    // separate "was it parsed" to ask alongside it.
    static bool attribute_exists(const md_system_t& sys, md_str_t path) {
        return md_attributes_find(&sys.attributes, path) != nullptr;
    }

    // One single valued attribute as a double, or 'fallback' when the path is absent. Anything that
    // is not single valued returns the fallback too: md_attribute_extract_f64 refuses a cap smaller
    // than the attribute, which is exactly the guard wanted here.
    static double attribute_scalar_f64(const md_system_t& sys, md_str_t path, double fallback = 0.0) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        if (!attr) return fallback;
        double value = 0.0;
        if (md_attribute_extract_f64(&value, 1, attr, md_attribute_slice_all(), md_unit_none()) != 1) return fallback;
        return value;
    }

    // The length of a rank 1 series without reading it. Zero when the path is absent, which is how
    // "how many optimisation steps / response frequencies / normal modes does this file carry" is
    // asked now that there is no reader to ask it of.
    static size_t attribute_series_count(const md_system_t& sys, md_str_t path) {
        size_t count = 0;
        return attribute_series_f64(&count, sys, path) ? count : 0;
    }

    // A series that is only handed out when it covers the index space the caller is going to read it
    // over. Presence alone is not enough: an optional column is indexed with a count that came from
    // somewhere else (the coefficient matrix, the frequency axis), and a file that publishes a
    // shorter column - or a reader that fills it differently - would otherwise be read past its end.
    // NULL when absent, not a plain F64 series, or shorter than 'required'.
    static const double* attribute_series_f64_covering(const md_system_t& sys, md_str_t path, size_t required) {
        size_t count = 0;
        const double* data = attribute_series_f64(&count, sys, path);
        return (data && count >= required) ? data : nullptr;
    }

    // A rank 2 {R,C} block of doubles, straight out of the table's storage. NULL when the path is
    // absent, computed, or not that shape. Same opt-in as attribute_series_f64 above and for the
    // same reason: these are resident F64 and the callers want them at that precision.
    static const double* attribute_matrix_f64(size_t* out_rows, size_t* out_cols, const md_system_t& sys, md_str_t path) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        const double* values = (const double*)md_attribute_view(attr, MD_ATTRIBUTE_TYPE_F64, 1, 2);
        if (!values) return nullptr;

        if (out_rows) *out_rows = attr->format.shape[0];
        if (out_cols) *out_cols = attr->format.shape[1];
        return values;
    }

    // One row of such a matrix. NULL when the row is out of range, which is what a caller holding a
    // selection index left over from a previously loaded file needs - rather than a read past the
    // end that looks like data.
    static const double* attribute_matrix_row_f64(size_t* out_cols, const md_system_t& sys, md_str_t path, size_t row) {
        size_t num_rows = 0;
        size_t num_cols = 0;
        const double* data = attribute_matrix_f64(&num_rows, &num_cols, sys, path);
        if (!data || row >= num_rows) return nullptr;
        if (out_cols) *out_cols = num_cols;
        return data + row * num_cols;
    }

    // A rank 2 {R,C} of THREE component values - one vector per (r,c). This is the shape a per atom
    // vector series over some other axis is published at: qm/atom/normal_mode is {M,N} displacements,
    // vlx/opt/coordinate is {P,N} geometries. Components and rank are separate things, which is why
    // this is not the matrix above with C*3 columns.
    static const dvec3_t* attribute_vec3_rows(size_t* out_rows, size_t* out_cols, const md_system_t& sys, md_str_t path) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        if (!attr || !attr->data) return nullptr;
        if (attr->format.type != MD_ATTRIBUTE_TYPE_F64) return nullptr;
        if (attr->format.rank != 2 || attr->format.components != 3) return nullptr;

        if (out_rows) *out_rows = attr->format.shape[0];
        if (out_cols) *out_cols = attr->format.shape[1];
        return (const dvec3_t*)attr->data;
    }

    // One row of it - one normal mode's displacements, one optimisation step's geometry. NULL when
    // the row is out of range, or when a row holds fewer than min_cols entries: every caller walks
    // the row over the QM atoms, so a row that does not cover them is not one to read.
    static const dvec3_t* attribute_vec3_row(const md_system_t& sys, md_str_t path, size_t row, size_t min_cols) {
        size_t num_rows = 0;
        size_t num_cols = 0;
        const dvec3_t* data = attribute_vec3_rows(&num_rows, &num_cols, sys, path);
        if (!data || row >= num_rows || num_cols < min_cols) return nullptr;
        return data + row * num_cols;
    }

    // A rank 1 series of 3 component values: dipole/<group>/vector is one per excited state.
    static const dvec3_t* attribute_vec3_series(size_t* out_count, const md_system_t& sys, md_str_t path) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        if (!attr || !attr->data) return nullptr;
        if (attr->format.type != MD_ATTRIBUTE_TYPE_F64) return nullptr;
        if (attr->format.rank != 1 || attr->format.components != 3) return nullptr;

        if (out_count) *out_count = attr->format.shape[0];
        return (const dvec3_t*)attr->data;
    }

    // A rank 1 column of a given STORED type, for the sibling columns a record was published as.
    // Typed rather than converted: an atom index is an index, and a widened double of one is a
    // number that has lost the only thing that made it an index.
    static const void* attribute_column(size_t* out_count, const md_system_t& sys, md_str_t path, md_attribute_type_t type) {
        const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
        if (!attr || !attr->data) return nullptr;
        if (attr->format.type != type) return nullptr;
        if (attr->format.rank != 1 || attr->format.components != 1) return nullptr;

        if (out_count) *out_count = attr->format.shape[0];
        return attr->data;
    }

    // ---- What KIND of calculation this is -------------------------------------------------------
    //
    // From vlx/rsp/type and vlx/opt/type, which the reader publishes as the lowercased enumerator
    // name. Text in the table and an enum here: nothing stored depends on the ordering of a header
    // this component no longer includes, and a second QM reader saying "rixs" about its own run
    // lands in the same branch without either side having heard of the other.
    //
    // An absent path is NONE, and every branch below treats that as "nothing of this kind to draw" -
    // which is also what a file with no response section at all gets, correctly.

    enum RspType {
        RSP_NONE = 0,
        RSP_LINEAR,
        RSP_CPP,
        RSP_C6,
        RSP_TPA,
        RSP_TPA_TRANSITION,
        RSP_RIXS,
    };

    static RspType response_type(const md_system_t& sys) {
        const md_str_t s = attribute_str(sys, MD_STR_LIT("vlx/rsp/type"));
        if (md_str_eq(s, MD_STR_LIT("linear")))         return RSP_LINEAR;
        if (md_str_eq(s, MD_STR_LIT("cpp")))            return RSP_CPP;
        if (md_str_eq(s, MD_STR_LIT("c6")))             return RSP_C6;
        if (md_str_eq(s, MD_STR_LIT("tpa")))            return RSP_TPA;
        if (md_str_eq(s, MD_STR_LIT("tpa_transition"))) return RSP_TPA_TRANSITION;
        if (md_str_eq(s, MD_STR_LIT("rixs")))           return RSP_RIXS;
        return RSP_NONE;
    }

    enum OptType {
        OPT_NONE = 0,
        OPT_GEOMETRY,
        OPT_CONSTRAINED,
        OPT_TS,
        OPT_IRC,
    };

    static OptType optimization_type(const md_system_t& sys) {
        const md_str_t s = attribute_str(sys, MD_STR_LIT("vlx/opt/type"));
        if (md_str_eq(s, MD_STR_LIT("geometry")))         return OPT_GEOMETRY;
        if (md_str_eq(s, MD_STR_LIT("constrained")))      return OPT_CONSTRAINED;
        if (md_str_eq(s, MD_STR_LIT("transition_state"))) return OPT_TS;
        if (md_str_eq(s, MD_STR_LIT("irc")))              return OPT_IRC;
        return OPT_NONE;
    }

    // Forgets the shown graph and whatever produced it, so the next frame starts from nothing. A run
    // in flight is told to stop and is dropped when the poll collects it.
    void critical_points_reset() {
        CriticalPoints& cp = critical_points;
        if (cp.gto_job) {
            cp.gto_job->cancel  = 1;
            cp.gto_job->dropped = true;
        }
        cp.gto_want  = 0;
        cp.gto_have  = 0;
        cp.raw_hash  = 0;
        cp.simp_hash = 0;
        md_topo_extremum_graph_free(&cp.raw_graph);
        md_topo_extremum_graph_free(&cp.simp_graph);
        // Indices into the graph just freed
        md_bitmask_clear(&cp.selection_mask);
        md_bitmask_clear(&cp.highlight_mask);
    }

    // Collects a finished analytic run: shown if it is still what is wanted, dropped otherwise.
    void critical_points_poll_gto() {
        CriticalPoints& cp = critical_points;
        GtoCriticalPointJob* job = cp.gto_job;
        if (!job || !job->done.load(std::memory_order_acquire)) return;
        cp.gto_job = nullptr;

        if (!job->dropped && job->hash == cp.gto_want) {
            md_topo_extremum_graph_free(&cp.raw_graph);
            md_topo_extremum_graph_free(&cp.simp_graph);
            if (job->graph.num_vertices > 0) {
                md_topo_extremum_graph_copy(&cp.raw_graph,  &job->graph);
                // Nothing to simplify: every vertex is a proven critical point, and two saddles
                // between the same pair of maxima are two saddles, not a sampling artefact.
                md_topo_extremum_graph_copy(&cp.simp_graph, &job->graph);
            }
            cp.gto_have    = job->hash;
            cp.gto_info    = job->info;
            cp.gto_failed  = !job->ok && !job->info.cancelled;
            cp.gto_seconds = job->seconds;
            cp.gto_rho_min = job->desc.rho_min;
        }
        gto_critical_point_job_destroy(job);
    }

    // The analytic search, once a frame while it is enabled, whether or not its window is open. A new run starts
    // whenever what the result depends on changes - the geometry the basis sits on, the threshold, the orbitals -
    // and a run made stale by such a change is stopped and dropped, never shown.
    void update_critical_points_gto(ApplicationState& state) {
        CriticalPoints& cp = critical_points;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        size_t num_atoms = 0;
        const vec3_t* atom_pos = basis_atom_positions_extract(&num_atoms, temp, state.mold.sys, state.mold.state);
        uint64_t hash = md_hash64(&cp.rho_min, sizeof(cp.rho_min), qm_version);
        hash = md_hash64(&cp.gto_use_gpu, sizeof(cp.gto_use_gpu), hash);
        if (atom_pos) hash = md_hash64(atom_pos, sizeof(vec3_t) * num_atoms, hash);

        // Every threshold is a new search, so none is started while its slider is held.
        if (!cp.rho_min_editing && hash != cp.gto_want) {
            // A point drawn where it was for other inputs would read as a result.
            critical_points_reset();
            cp.gto_want = hash;
        }
        if (!cp.gto_job && cp.gto_have != cp.gto_want) {
            cp.gto_job = gto_critical_point_job_start(state, cp.rho_min, cp.gto_use_gpu, cp.gto_want);
            if (!cp.gto_job) {
                // Logged by the start. Remembered, so it is not retried every frame.
                cp.gto_have    = cp.gto_want;
                cp.gto_info    = {};
                cp.gto_failed  = true;
                cp.gto_seconds = 0.0;
            }
        }
    }

    // The discrete pipeline, once a frame while it is enabled: the density sampled onto a grid, the critical
    // points of the samples (md_topo), then md_topo_simplify.
    void update_critical_points_volume(ApplicationState& state) {
        CriticalPoints& cp = critical_points;
        if (cp.vol_per_frame) {
            cp.raw_hash = 0;
        }

        // Update volume if required
        const double samples_per_unit_length = (float)(volume_resolution_samples_per_angstrom[(int)cp.vol_res] * BOHR_TO_ANGSTROM);
        md_grid_t grid = { 0 };
#if USE_OABB_GRID
        init_grid(&grid, oabb.orientation, oabb.center, oabb.extents, samples_per_unit_length);
#else
        init_grid(&grid, mat3_ident(), aabb.min_ext, aabb.max_ext, samples_per_unit_length);
#endif

        // Compute a hash of the current state
        size_t frame_idx = (size_t)(state.animation.frame + 0.5);
        uint64_t raw_hash = md_hash64(&grid, sizeof(grid), frame_idx);

        if (cp.raw_hash != raw_hash) {
            cp.raw_hash = raw_hash;
            md_topo_extremum_graph_free(&cp.raw_graph);

            // The density comes out of the attribute table, evaluated by the same code
            // path a representation uses. The only thing special here is the tail: the
            // topology pass reads what the density pass wrote, on the same stream, so
            // the volume is never read back to the host at all.
#if MD_ENABLE_GPU
            // Ensure the topo context exists with the correct volume dimensions
            uint32_t new_dims[3] = { (uint32_t)grid.dim[0], (uint32_t)grid.dim[1], (uint32_t)grid.dim[2] };
            if (!cp.topo_ctx ||
                cp.topo_dims[0] != new_dims[0] ||
                cp.topo_dims[1] != new_dims[1] ||
                cp.topo_dims[2] != new_dims[2])
            {
                md_topo_gpu_context_destroy(cp.topo_ctx);
                cp.topo_ctx = md_topo_gpu_context_create(state.gpu_device, new_dims[0], new_dims[1], new_dims[2]);
                cp.topo_dims[0] = new_dims[0];
                cp.topo_dims[1] = new_dims[1];
                cp.topo_dims[2] = new_dims[2];
            }

            if (!cp.topo_ctx) {
                MD_LOG_ERROR("Failed to create topo GPU context for electron density");
            } else if (density_evaluate_to_gpu_volume(&state, grid, critical_points_density_path(state.mold.sys), nullptr, MD_GTO_OP_SET)) {
                // Program order on the stream is the only dependency needed.
                md_topo_gpu_record(state.gpu_stream, cp.topo_ctx, state.gpu_volume, &grid, 0.0f);
                md_gpu_stream_sync(state.gpu_stream);
                md_topo_gpu_context_extract(&cp.raw_graph, cp.topo_ctx);
            } else {
                MD_LOG_ERROR("Failed to evaluate the electron density for critical points analysis");
            }
#else
            md_temp_scope_t temp = md_temp_begin();
            MD_DEFER { md_temp_end(temp); };
            uint32_t tex_id = 0;
            if (gl::init_texture_3D(&tex_id, grid.dim[0], grid.dim[1], grid.dim[2], GL_R32F)) {
                size_t num_atom_pos = 0;
                const vec3_t* atom_pos = basis_atom_positions_extract(&num_atom_pos, temp, state.mold.sys, state.mold.state);
                if (density_evaluate_gl(tex_id, grid, state.mold.sys, atom_pos, num_atom_pos, critical_points_density_path(state.mold.sys), nullptr, MD_GTO_OP_SET)) {
                    if (!md_topo_compute_extremum_graph_GPU(&cp.raw_graph, tex_id, &grid, 0.0f)) {
                        MD_LOG_ERROR("Failed to compute extremum graph for electron density");
                    }
                }
                glDeleteTextures(1, &tex_id);
            }
#endif
        }

        uint64_t simp_hash = raw_hash;
        if (cp.vol_simplify) {
            simp_hash = md_hash64(&cp.vol_simp_threshold, sizeof(cp.vol_simp_threshold), raw_hash);
        }

        if (cp.simp_hash != simp_hash) {
            cp.simp_hash = simp_hash;

            md_topo_extremum_graph_free(&cp.simp_graph);
            if (cp.vol_simplify && cp.raw_graph.num_vertices > 0) {
                md_topo_simplify(&cp.simp_graph, &cp.raw_graph, cp.vol_simp_threshold, true);
            } else {
                md_topo_extremum_graph_copy(&cp.simp_graph, &cp.raw_graph);
            }
        }
    }

    // Once a frame: the search of the selected method while enabled, whether or not the window is open, so what
    // is drawn follows the geometry and the orbitals.
    void update_critical_points(ApplicationState& state) {
        CriticalPoints& cp = critical_points;
        if (!cp.enabled) {
            // Switched off mid search: nothing is going to show the result.
            if (cp.gto_job) critical_points_reset();
            return;
        }
        if (cp.method == CriticalPoints::Method_Analytic) {
            update_critical_points_gto(state);
        } else {
            update_critical_points_volume(state);
        }
    }

    // Where the analytic search is: its progress and a cancel button while it runs, else its outcome and what it took
    void draw_critical_points_gto_status() {
        CriticalPoints& cp = critical_points;
        if (cp.gto_job) {
            const GtoCriticalPointJob* job = cp.gto_job;
            char buf[64];
            snprintf(buf, sizeof(buf), "Searching %.0f%% (%.0f s)", job->progress * 100.0f, md_tick_to_seconds(md_tick_now() - job->start));
            const float cancel_width = ImGui::CalcTextSize("Cancel").x + ImGui::GetStyle().FramePadding.x * 2.0f + ImGui::GetStyle().ItemSpacing.x;
            ImGui::ProgressBar(job->progress, ImVec2(ImGui::GetContentRegionAvail().x - cancel_width, 0), buf);
            ImGui::SetItemTooltip("How far the search is: an estimate, as the work of the search unfolds as it goes");
            ImGui::SameLine();
            if (ImGui::Button("Cancel")) {
                cp.gto_job->cancel = 1;   // collected as usual: what was found so far is shown, marked
            }
            return;
        }

        const ImVec4 amber = ImVec4(1.0f, 0.75f, 0.3f, 1.0f);
        const md_topo_gto_info_t& info = cp.gto_info;
        if (cp.gto_failed) {
            ImGui::TextColored(IM_RED, "The search could not run, see the log");
            return;
        }
        if (cp.gto_have == 0) {
            return;     // nothing searched yet
        }
        if (info.cancelled) {
            ImGui::TextColored(amber, "Cancelled: only what was found so far is shown");
            ImGui::SameLine();
            if (ImGui::SmallButton("Restart")) {
                cp.gto_have = 0;
            }
        } else if (info.complete) {
            ImGui::TextColored(IM_GREEN, "Complete");
            ImGui::SetItemTooltip("Every critical point above the threshold is shown. Each was proven to be the only\n"
                                  "critical point in the cube it was found in, and every other cube to hold none.");
        } else {
            ImGui::TextColored(amber, "Incomplete: %u region%s left undecided", info.num_clusters, info.num_clusters == 1 ? "" : "s");
            ImGui::SetItemTooltip("Everything outside these regions is decided. Inside one, the gradient is too close to\n"
                                  "zero to separate critical points at the smallest cube size - typically a nearly\n"
                                  "degenerate pair, such as a bond and a ring critical point about to annihilate.");
        }
        ImGui::SameLine();
        if (info.used_gpu) {
            ImGui::TextDisabled("%.2f s on the GPU, %llu cubes", cp.gto_seconds, (unsigned long long)info.num_box_evals);
        } else {
            ImGui::TextDisabled("%.2f s on %u CPU threads, %llu cubes", cp.gto_seconds, info.num_threads, (unsigned long long)info.num_box_evals);
        }
        if (ImGui::BeginItemTooltip()) {
            // where the time went, to see which part to work on
            ImGui::Text("Setup: %.1f ms", info.ms_setup);
            ImGui::Text("Octree sweep: %.1f ms", info.ms_sweep);
            if (info.used_gpu) {
                const double host = info.ms_sweep - info.ms_sweep_gpu_wait - info.ms_sweep_polish - info.ms_sweep_cpu;
                ImGui::BulletText("GPU: %.1f ms, %llu cubes in %u dispatches", info.ms_sweep_gpu_wait, (unsigned long long)info.num_gpu_box_evals, info.num_gpu_dispatches);
                ImGui::BulletText("Root polish (CPU, double): %.1f ms", info.ms_sweep_polish);
                ImGui::BulletText("Escalated subtrees (CPU, double): %.1f ms, %llu cubes from %u", info.ms_sweep_cpu,
                                  (unsigned long long)(info.num_box_evals - info.num_gpu_box_evals), info.num_escalated_boxes);
                ImGui::BulletText("Other host work: %.1f ms", host > 0.0 ? host : 0.0);
            }
            ImGui::Text("Separatrices: %.1f ms", info.ms_separatrices);
            ImGui::Text("Clusters: %.1f ms", info.ms_clusters);
            if (info.used_gpu && info.num_escalated_boxes > 0) {
                ImGui::TextDisabled("Escalated: cubes too close to a critical point for fp32, finished in double on the CPU.");
            }
            ImGui::EndTooltip();
        }
        if (!info.cancelled && !info.complete && info.num_clusters > 0 && ImGui::TreeNode("Undecided regions")) {
            for (uint32_t i = 0; i < info.num_clusters; ++i) {
                const md_topo_gto_cluster_t& c = info.clusters[i];
                const double x = 0.5 * (c.lo[0] + c.hi[0]) * BOHR_TO_ANGSTROM;
                const double y = 0.5 * (c.lo[1] + c.hi[1]) * BOHR_TO_ANGSTROM;
                const double z = 0.5 * (c.lo[2] + c.hi[2]) * BOHR_TO_ANGSTROM;
                if (c.degree_residual < 0.1f) {
                    ImGui::BulletText("(%.3f, %.3f, %.3f): index sum %+d", x, y, z, c.degree);
                } else {
                    ImGui::BulletText("(%.3f, %.3f, %.3f): index sum undetermined", x, y, z);
                }
                ImGui::SetItemTooltip("Centre of the region in \xC3\x85, %u cubes.\n"
                                      "The index sum is the degree of the density gradient on the region's boundary\n"
                                      "(computed numerically): the sum over the critical points inside of -1 for an\n"
                                      "attractor (maximum) or RCP (join saddle) and +1 for a BCP (split saddle) or CCP\n"
                                      "(minimum). Non-zero proves at least one critical point is in there; zero means\n"
                                      "none, or pairs that cancel.", c.num_boxes);
            }
            if (info.clusters_truncated) {
                ImGui::BulletText("...");
            }
            ImGui::TreePop();
        }
    }

    // The critical points of the density in a window of their own: what is searched and how, how the points are
    // drawn and the export on two lines, where the search is, the legend, and the rest for the table of the points
    // (hovered and selected with the 3D view, both ways). The search runs whether or not this window is open
    // (update_critical_points).
    void draw_critical_points_window(ApplicationState& state) {
        CriticalPoints& cp = critical_points;
        cp.rho_min_editing = false;
        if (!cp.show_window) {
            return;
        }
        ImGui::SetNextWindowSize({ 640, 600 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Critical Points", &cp.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            draw_critical_points_window_contents(state);
        }
        ImGui::End();
    }

    void draw_critical_points_window_contents(ApplicationState& state) {
        CriticalPoints& cp = critical_points;
        const float em = ImGui::GetFontSize();
        const md_topo_extremum_graph_t& graph = cp.simp_graph;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        // --- what is searched, and how
        ImGui::Checkbox("Enabled", &cp.enabled);
        ImGui::SetItemTooltip("Search the electron density for its critical points, and show them");
        ImGui::SameLine();
        static const char* method_str[] = { "Analytic", "Volume" };
        ImGui::SetNextItemWidth(6.5f * em);
        if (ImGui::Combo("##Method", &cp.method, method_str, (int)MD_COUNTOF(method_str))) {
            critical_points_reset();   // what is shown came from the other method
        }
        ImGui::SetItemTooltip("Analytic: searches the GTO density itself. Every critical point above the threshold is\n"
                              "found and placed to numerical precision, with a proof that none was missed.\n"
                              "Volume: samples the density onto a grid and finds critical points of the samples,\n"
                              "placed no better than to the voxel.");
        if (cp.method == CriticalPoints::Method_Analytic) {
            ImGui::SameLine();
            ImGui::SetNextItemWidth(8.0f * em);
            ImGui::SliderFloat("##Density threshold", &cp.rho_min, 1.0e-5f, 1.0e-1f, "\xCF\x81 > %.1e", ImGuiSliderFlags_Logarithmic);
            cp.rho_min_editing = ImGui::IsItemActive();
            ImGui::SetItemTooltip("The density threshold: critical points where the total electron density is below this (au)\n"
                                  "are neither searched for nor reported. A lower threshold searches a larger region and takes longer.");
#if MD_ENABLE_GPU
            if (state.gpu_device) {
                ImGui::SameLine();
                ImGui::Checkbox("GPU", &cp.gto_use_gpu);
                ImGui::SetItemTooltip("On: the octree sweep runs on the GPU in fp32 with rounding margins, and the CPU polishes\n"
                                      "the critical points in double. Off: everything runs on the CPU. Same result either way.");
            }
#endif
        } else {
            ImGui::SameLine();
            ImGui::SetNextItemWidth(5.5f * em);
            ImGui::Combo("##Resolution", (int*)&cp.vol_res, volume_resolution_str, (int)VolumeResolution::Count);
            ImGui::SetItemTooltip("The resolution of the grid the density is sampled on");
            ImGui::SameLine();
            ImGui::Checkbox("Per frame", &cp.vol_per_frame);
            ImGui::SetItemTooltip("Sample the density and extract the graph again every frame");
            ImGui::SameLine();
            ImGui::Checkbox("Simplify", &cp.vol_simplify);
            ImGui::SetItemTooltip("Leave out the critical points below a density threshold, and duplicate saddles");
            if (cp.vol_simplify) {
                ImGui::SameLine();
                ImGui::SetNextItemWidth(7.0f * em);
                ImGui::SliderFloat("##Simplification threshold", &cp.vol_simp_threshold, 0.0f, 1.0E-2f, "\xCF\x81 > %.1e", ImGuiSliderFlags_Logarithmic);
                ImGui::SetItemTooltip("The density below which critical points are left out (au)");
            }
        }

        // --- how they are drawn, and the export at the end of the line
        ImGui::SetNextItemWidth(6.0f * em);
        ImGui::SliderFloat("Size", &cp.point_scale, 0.25f, 4.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
        ImGui::SetItemTooltip("The size of the markers");
        ImGui::SameLine();
        ImGui::Checkbox("Paths", &cp.show_paths);
        ImGui::SetItemTooltip("The bond paths, from each bond critical point to the two attractors it joins, and the\n"
                              "RCP attractor paths, from each ring critical point to the cage critical points beside it");
        {
            const char* lbl = "Export...";
            const float w = ImGui::CalcTextSize(lbl).x + 2.0f * ImGui::GetStyle().FramePadding.x;
            ImGui::SameLine();
            const float avail = ImGui::GetContentRegionAvail().x;
            if (avail > w) ImGui::SetCursorPosX(ImGui::GetCursorPosX() + avail - w);
            ImGui::BeginDisabled(!cp.enabled || graph.num_vertices == 0 || cp.gto_job != nullptr);
            if (ImGui::Button(lbl)) {
                critical_points_export(state);
            }
            ImGui::EndDisabled();
            ImGui::SetItemTooltip("The critical points and the paths between them, every number at the precision it was\n"
                                  "computed in: positions (bohr), density and, from the analytic search, the Hessian there.\n"
                                  ".json for scripts, .vtp (VTK PolyData, as the Topology ToolKit writes a Morse-Smale\n"
                                  "complex) for ParaView");
        }

        if (!cp.enabled) {
            ImGui::TextDisabled("Not searched: tick Enabled to search the density for its critical points");
            return;
        }

        // --- where the search is
        if (cp.method == CriticalPoints::Method_Analytic) {
            draw_critical_points_gto_status();
        }

        // --- the counts, as the legend of the markers: two by two
        uint32_t vertex_count[MD_TOPO_NUM_TYPES] = { 0 };
        md_topo_count_vertex_types(vertex_count, &graph);
        if (ImGui::BeginTable("##Critical point legend", 2, ImGuiTableFlags_SizingStretchSame)) {
            for (int k = 0; k < (int)MD_COUNTOF(critical_point_type_order); ++k) {
                const int type = critical_point_type_order[k];
                const CriticalPointStyle& style = critical_point_style(type);
                ImGui::TableNextColumn();
                critical_point_swatch(type);
                ImGui::SameLine();
                ImGui::Text("%s %s  %u", style.plural, style.signature, vertex_count[type]);
                ImGui::SameLine();
                ImGui::TextDisabled("%s", style.topology_plural);
                ImGui::SetItemTooltip("%s %s" TOOLTIP_SEP "%s", style.name, style.signature, style.topology);
            }
            ImGui::EndTable();
        }

        if (graph.num_vertices == 0) {
            return;
        }

        // --- the table of the points, in the rest of the window
        enum CriticalPointColumn {
            CriticalPointColumn_Idx,
            CriticalPointColumn_Type,
            CriticalPointColumn_Value,
            CriticalPointColumn_Laplacian,
            CriticalPointColumn_X,
            CriticalPointColumn_Y,
            CriticalPointColumn_Z,
            CriticalPointColumn_Count
        };
        // The Laplacian where the producer has the Hessian (the analytic search)
        const md_topo_vert_f64_t* f64 = graph.vertices_f64;

        uint32_t* row_indices = md_temp_alloc_array(temp, uint32_t, graph.num_vertices);
        for (uint32_t i = 0; i < graph.num_vertices; ++i) {
            row_indices[i] = i;
        }

        static const ImGuiTableFlags table_flags = ImGuiTableFlags_RowBg |
                                                   ImGuiTableFlags_Borders |
                                                   ImGuiTableFlags_ScrollY |
                                                   ImGuiTableFlags_Sortable |
                                                   ImGuiTableFlags_SortMulti |
                                                   ImGuiTableFlags_Resizable |
                                                   ImGuiTableFlags_SizingStretchProp;

        if (ImGui::IsWindowHovered()) {
            md_bitmask_clear(&critical_points.highlight_mask);
        }

        if (ImGui::BeginTable("Critical Points Table", CriticalPointColumn_Count, table_flags, ImVec2(0.0f, ImGui::GetContentRegionAvail().y))) {
            ImGui::TableSetupColumn("Idx",          ImGuiTableColumnFlags_DefaultSort, 0.55f, CriticalPointColumn_Idx);
            ImGui::TableSetupColumn("Type",         ImGuiTableColumnFlags_None,        1.15f, CriticalPointColumn_Type);
            ImGui::TableSetupColumn("\xCF\x81 (au)", ImGuiTableColumnFlags_None,       1.0f,  CriticalPointColumn_Value);
            ImGui::TableSetupColumn("Laplacian",    f64 ? ImGuiTableColumnFlags_None : ImGuiTableColumnFlags_Disabled, 1.25f, CriticalPointColumn_Laplacian);
            ImGui::TableSetupColumn("X (\xC3\x85)", ImGuiTableColumnFlags_None,        0.9f,  CriticalPointColumn_X);
            ImGui::TableSetupColumn("Y (\xC3\x85)", ImGuiTableColumnFlags_None,        0.9f,  CriticalPointColumn_Y);
            ImGui::TableSetupColumn("Z (\xC3\x85)", ImGuiTableColumnFlags_None,        0.9f,  CriticalPointColumn_Z);
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableHeadersRow();

            if (ImGui::IsWindowHovered()) {
                md_bitmask_clear(&critical_points.highlight_mask);
            }

            if (ImGuiTableSortSpecs* sort_specs = ImGui::TableGetSortSpecs()) {
                if (sort_specs->SpecsCount > 0) {
                    auto laplacian = [f64](uint32_t i) {
                        return f64 ? f64[i].hessian_eigenvalues[0] + f64[i].hessian_eigenvalues[1] + f64[i].hessian_eigenvalues[2] : 0.0;
                    };
                    auto compare_rows = [&graph, sort_specs, &laplacian](uint32_t lhs_idx, uint32_t rhs_idx) {
                        const md_topo_vert_t& lhs = graph.vertices[lhs_idx];
                        const md_topo_vert_t& rhs = graph.vertices[rhs_idx];

                        for (int spec_idx = 0; spec_idx < sort_specs->SpecsCount; ++spec_idx) {
                            const ImGuiTableColumnSortSpecs& spec = sort_specs->Specs[spec_idx];
                            int cmp = 0;
                            switch ((CriticalPointColumn)spec.ColumnUserID) {
                            case CriticalPointColumn_Idx:
                                cmp = (lhs_idx < rhs_idx) ? -1 : (lhs_idx > rhs_idx ? 1 : 0);
                                break;
                            case CriticalPointColumn_Type: {
                                const int a = critical_point_type_rank(graph.types[lhs_idx]), b = critical_point_type_rank(graph.types[rhs_idx]);
                                cmp = (a < b) ? -1 : (a > b ? 1 : 0);
                                break;
                            }
                            case CriticalPointColumn_Value:
                                cmp = (lhs.value < rhs.value) ? -1 : (lhs.value > rhs.value ? 1 : 0);
                                break;
                            case CriticalPointColumn_Laplacian: {
                                const double a = laplacian(lhs_idx), b = laplacian(rhs_idx);
                                cmp = (a < b) ? -1 : (a > b ? 1 : 0);
                                break;
                            }
                            case CriticalPointColumn_X:
                                cmp = (lhs.x < rhs.x) ? -1 : (lhs.x > rhs.x ? 1 : 0);
                                break;
                            case CriticalPointColumn_Y:
                                cmp = (lhs.y < rhs.y) ? -1 : (lhs.y > rhs.y ? 1 : 0);
                                break;
                            case CriticalPointColumn_Z:
                                cmp = (lhs.z < rhs.z) ? -1 : (lhs.z > rhs.z ? 1 : 0);
                                break;
                            default:
                                break;
                            }

                            if (cmp != 0) {
                                return spec.SortDirection == ImGuiSortDirection_Ascending ? cmp < 0 : cmp > 0;
                            }
                        }

                        return lhs_idx < rhs_idx;
                    };

                    std::sort(row_indices, row_indices + graph.num_vertices, compare_rows);
                    sort_specs->SpecsDirty = false;
                }
            }

            ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
            ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);

            ImGuiListClipper clipper;
            clipper.Begin((int)graph.num_vertices);
            while (clipper.Step()) {
                for (int row = clipper.DisplayStart; row < clipper.DisplayEnd; ++row) {
                    const uint32_t idx = row_indices[row];
                    const md_topo_vert_t& v = graph.vertices[idx];
                    const int type = graph.types[idx];

                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();

                    ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                    bool is_sel = md_bitmask_test_bit(&critical_points.selection_mask, idx);
                    bool is_hov = md_bitmask_test_bit(&critical_points.highlight_mask, idx);

                    if (is_hov) {
                        ImGui::PushStyleColor(ImGuiCol_Header, IM_YELLOW);
                    }
                    else {
                        ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                    }

                    char lbl[16];
                    snprintf(lbl, sizeof(lbl), "%u", idx);
                    ImGui::Selectable(lbl, is_sel || is_hov, selectable_flags);

                    // The selectable spans the row: hovered anywhere on it (a clipped table's row numbers are not the rows')
                    if (ImGui::IsItemHovered()) {
                        md_bitmask_clear(&critical_points.highlight_mask);
                        md_bitmask_or(&critical_points.highlight_mask, md_bitset_range(idx, idx + 1));

                        //Selection
                        if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
                            if (ImGui::IsKeyDown(ImGuiKey_MouseLeft)) {
                                md_bitmask_or(&critical_points.selection_mask, md_bitset_range(idx, idx + 1));
                            } else if (ImGui::IsKeyDown(ImGuiKey_MouseRight)) {
                                md_bitmask_andnot(&critical_points.selection_mask, md_bitset_range(idx, idx + 1));
                            }
                        }
                    }
                    ImGui::PopStyleColor(1);

                    ImGui::TableNextColumn();
                    {
                        const CriticalPointStyle& style = critical_point_style(type);
                        critical_point_swatch(type);
                        ImGui::SameLine();
                        ImGui::TextUnformatted(style.label);
                        ImGui::SetItemTooltip("%s %s" TOOLTIP_SEP "%s", style.name, style.signature, style.topology);
                    }
                    ImGui::TableNextColumn();
                    ImGui::Text("%.6g", f64 ? f64[idx].value : (double)v.value);
                    ImGui::TableNextColumn();
                    if (f64) {
                        const double* ev = f64[idx].hessian_eigenvalues;
                        ImGui::Text("%.6g", ev[0] + ev[1] + ev[2]);
                        ImGui::SetItemTooltip("The Laplacian of the density, the sum of the Hessian's eigenvalues:\n%.6g, %.6g, %.6g (au)", ev[0], ev[1], ev[2]);
                    }
                    // Positions in \xC3\x85, as everywhere else in the application (the graph's are in bohr). Shown to
                    // four decimals; the export has them in full.
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4f", (f64 ? f64[idx].x : (double)v.x) * BOHR_TO_ANGSTROM);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4f", (f64 ? f64[idx].y : (double)v.y) * BOHR_TO_ANGSTROM);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4f", (f64 ? f64[idx].z : (double)v.z) * BOHR_TO_ANGSTROM);
                }
            }

            ImGui::PopStyleColor(2);
            ImGui::EndTable();
        }
    }

    // Writes the critical points and the paths between them to a file the user picks (critical_points_write).
    void critical_points_export(const ApplicationState& state) {
        char path_buf[2048];
        if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, MD_STR_LIT("json,vtp"))) {
            critical_points_write(state, md_str_from_cstr(path_buf));
        }
    }

    // The critical points and the paths between them, in the file at path: .vtp for VTK PolyData as the Topology
    // ToolKit writes a Morse-Smale complex (ParaView, TTK), anything else as JSON. No format for this is standard:
    // QTAIM programs have their own (AIMAll's .sum), and VTK is what topological data analysis exchanges.
    bool critical_points_write(const ApplicationState& state, md_str_t path) {
        const CriticalPoints& cp = critical_points;
        md_str_t ext = {};
        const bool vtp = md_extract_ext(&ext, path) && md_str_eq_cstr_ignore_case(ext, "vtp");

        md_file_t file = {};
        if (!md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
            MD_LOG_ERROR("Critical points: could not open '" MD_STR_FMT "' for writing", MD_STR_ARG(path));
            return false;
        }
        MD_DEFER { md_file_close(&file); };

        CriticalPointsExport ex = {};
        ex.graph     = &cp.simp_graph;
        ex.analytic  = cp.method == CriticalPoints::Method_Analytic;
        ex.rho_min   = cp.gto_rho_min;
        ex.complete  = cp.gto_info.complete;
        ex.cancelled = cp.gto_info.cancelled;
        ex.field     = critical_points_density_path(state.mold.sys);
        if (vtp) {
            critical_points_write_vtp(file, ex);
        } else {
            critical_points_write_json(file, ex);
        }
        MD_LOG_INFO("Critical points: %u points and %u paths written to '" MD_STR_FMT "'", cp.simp_graph.num_vertices, cp.simp_graph.num_edges, MD_STR_ARG(path));
        return true;
    }

    void draw_summary_window(ApplicationState& state) {
        const md_system_t& sys = state.mold.sys;

        // The QM atom domain, gathered per frame rather than cached. The System Information and
        // Geometry sections below both index it.
        QmAtoms qm = {};
        es_qm_atoms(&qm, sys);

        if (!summary.show_window) {
            // Parked on the converged geometry while the window is closed, so re-opening it does not
            // step the molecule back to wherever the last session left the slider.
            opt.selected = (int)attribute_series_count(sys, MD_STR_LIT("vlx/opt/energy")) - 1;
            return;
        }

        // Same as the Response window: keep a margin between the data and the plot frame instead
        // of letting the curves run into the edges. ImPlot applies this only to axes it FITS, so
        // it covers the optimization plot and the SCF plot's linear energy axis; the SCF plot's
        // log gradient axis cannot use it and pads itself, see below.
        ImPlot::PushStyleVar(ImPlotStyleVar_FitPadding, ImVec2(0.05f, 0.1f));
        MD_DEFER { ImPlot::PopStyleVar(); };

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        const RspType rsp_type = response_type(sys);

        // The actual plot
        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Summary", &summary.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::TreeNode("Level of Calculation")) {
                const md_str_t basis_set = attribute_str(sys, MD_STR_LIT("vlx/basis_set"));
                const md_str_t dft_func   = attribute_str(sys, MD_STR_LIT("vlx/dft_functional"));

                // Named from the same two questions the orbital grid branches on, so the label
                // cannot disagree with the layout it describes. A reader with no word of its own for
                // the method still gets an accurate one.
                const bool unrestricted = es_has_distinct_beta_orbitals(state.mold.sys);
                const bool open_shell   = !unrestricted && es_has_distinct_beta_occupations(state.mold.sys);
                const char* scf_type = unrestricted ? "Unrestricted" : open_shell ? "Restricted Open-Shell" : "Restricted";

                ImGui::Text("Method: %s", md_str_ptr(dft_func));
                ImGui::Text("SCF Type: %s", scf_type);
                ImGui::Text("Basis Set: %s", md_str_ptr(basis_set));
                ImGui::Spacing();
                
                ImGui::TreePop();
            }
            if (ImGui::TreeNode("System Information")) {
                ImGui::Text("Num Atoms:           %-6zu", qm.count);
                // Summed from the occupation columns. The multiplicity follows from the same two
                // numbers - 2S+1 with S the half difference - so it is arithmetic here rather than a
                // third quantity to publish and keep in agreement with the other two.
                const double n_alpha = es_electron_count(state.mold.sys, es_path::alpha_occupation);
                const double n_beta  = es_electron_count(state.mold.sys, es_path::beta_occupation);
                ImGui::Text("Num Alpha Electrons: %-6.0f", n_alpha);
                ImGui::Text("Num Beta Electrons:  %-6.0f", n_beta);
                ImGui::Text("Molecular Charge:    %-6f",  attribute_scalar_f64(sys, MD_STR_LIT("vlx/molecular_charge")));
                ImGui::Text("Spin Multiplicity:   %-6.0f", fabs(n_alpha - n_beta) + 1.0);
                if (rsp_type == RSP_C6) {
                    ImGui::Text("C6 Value:            %-12.6f (au)", attribute_scalar_f64(sys, MD_STR_LIT("vlx/rsp/c6")));
                }
                ImGui::Spacing();
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("SCF")) {
                // The convergence history comes out of the attribute table, put there when the
                // file was loaded, so this panel is one consumer among others and a file which
                // carried no history simply has no such path - there is no separate "was it
                // parsed" question to ask alongside it.
                size_t num_iter = 0;
                size_t num_grad = 0;
                const double* energy    = attribute_series_f64(&num_iter, sys, MD_STR_LIT("vlx/scf/history/energy"));
                const double* grad_norm = attribute_series_f64(&num_grad, sys, MD_STR_LIT("vlx/scf/history/gradient_norm"));

                // Siblings in one group share a shape; a file which somehow disagrees is not a
                // history anyone can plot.
                if (!energy || !grad_norm || num_grad != num_iter) {
                    energy = nullptr;
                    grad_norm = nullptr;
                    num_iter = 0;
                }

                double* energy_offsets = nullptr;
                double ref_energy = 0.0;

                // We set up iterations as doubles for easier use
                if (num_iter > 0) {
                    energy_offsets = md_temp_alloc_array(temp, double, num_iter);
                    ref_energy = energy[num_iter - 1];
                    for (size_t i = 0; i < num_iter; i++) {
                        energy_offsets[i] = fabs(energy[i] - ref_energy);
                    }
                }

                if (num_iter > 1) {
                    // The gradient axis is log10, and FitPadding cannot pad it: ImPlotAxis::ApplyFit
                    // pads in VALUE space, while a log axis is constrained to (DBL_MIN, inf). A fit
                    // spanning 1e-9 to 1e-1 has an extent of ~0.1, so padding it by 10% drives the
                    // minimum negative and it clamps three hundred decades below the data. Pad in
                    // decades instead and set the limits directly.
                    double grad_lo = DBL_MAX;
                    double grad_hi = -DBL_MAX;
                    for (size_t i = 0; i < num_iter; ++i) {
                        // Non-positive values have no place on a log axis and would only drag the
                        // fit down to the clamp.
                        if (grad_norm[i] > 0.0) {
                            grad_lo = MD_MIN(grad_lo, grad_norm[i]);
                            grad_hi = MD_MAX(grad_hi, grad_norm[i]);
                        }
                    }
                    const bool has_grad_range = (grad_lo <= grad_hi);
                    if (has_grad_range) {
                        const double log_lo = log10(grad_lo);
                        const double log_hi = log10(grad_hi);
                        // Ten percent of the span at each end, matching the linear axes - with a
                        // floor, because a run that barely moved has no span to take a share of.
                        const double log_pad = MD_MAX((log_hi - log_lo) * 0.1, 0.1);
                        grad_lo = pow(10.0, log_lo - log_pad);
                        grad_hi = pow(10.0, log_hi + log_pad);
                    }

                    // Refit when a different history arrives - the same idiom the Optimization plot
                    // below uses. In between, ImPlotCond_Once leaves the axes the user's to zoom
                    // and pan, which ImPlotCond_Always every frame would take away.
                    uint64_t plot_hash = md_hash64((const uint8_t*)grad_norm,      sizeof(double) * num_iter, 0);
                    plot_hash          = md_hash64((const uint8_t*)energy_offsets, sizeof(double) * num_iter, plot_hash);
                    const ImPlotCond axis_cond = (plot_hash != scf.plot_hash) ? ImPlotCond_Always : ImPlotCond_Once;
                    scf.plot_hash = plot_hash;

                    if (ImPlot::BeginPlot("SCF")) {
                        ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, (double)num_iter, axis_cond);
                        ImPlot::SetupLegend(ImPlotLocation_NorthEast);
                        ImPlot::SetupAxes("Iteration", "Gradient Norm (au)");
                        // We draw 2 y axis as "Energy total" has values in a different range then the rest of the data
                        ImPlot::SetupAxis(ImAxis_Y2, "Energy (au)", ImPlotAxisFlags_AuxDefault);
                        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                        if (has_grad_range) {
                            // After SetupAxes, which is what enables Y1.
                            ImPlot::SetupAxisLimits(ImAxis_Y1, grad_lo, grad_hi, axis_cond);
                        }

                        ImPlot::PlotLine("Gradient", grad_norm, (int)num_iter, 1.0, 1.0);
                        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                        ImPlot::PlotLine("Energy", energy_offsets, (int)num_iter, 1.0, 1.0);
                        //ImPlot::PlotLine("Density Change", iter, vlx.scf.iter.density_change, (int)vlx.scf.iter.count);
                        //ImPlot::PlotLine("Energy Change", iter, vlx.scf.iter.energy_change, (int)vlx.scf.iter.count);
                        //ImPlot::PlotLine("Max Gradient", iter, vlx.scf.iter.max_gradient, (int)vlx.scf.iter.count);
                        ImPlot::EndPlot();
                    }
                } else {
                    ImGui::Text("This dataset carries no SCF convergence history");
                }
                ImGui::Spacing();
                if (num_iter > 0) {
                    ImGui::Text("Total energy:              %16.10f (au)", energy[num_iter - 1]);
                    ImGui::Text("Gradient norm:             %16.10f (au)", grad_norm[num_iter - 1]);
                }
                ImGui::Text("Nuclear repulsion energy:  %16.10f (au)", attribute_scalar_f64(sys, MD_STR_LIT("vlx/nuclear_repulsion_energy")));
                ImGui::Spacing();
                ImGui::TreePop();
            }

            {
                size_t num_steps = 0;
                const double* opt_energies = attribute_series_f64(&num_steps, sys, MD_STR_LIT("vlx/opt/energy"));
                if (opt_energies && num_steps > 0) {
                    if (ImGui::TreeNode("Optimization")) {
                        // FitPadding is pushed once for the whole window, at the top of this function.
                        const OptType opt_type = optimization_type(sys);
                        // Published for an IRC only, where it names the step the path was walked out
                        // from. num_steps is out of range on purpose everywhere else, so the branch
                        // below falls through to the minimum without a second test for the kind.
                        const size_t ts_index = (size_t)attribute_scalar_f64(sys, MD_STR_LIT("vlx/opt/irc_ts_index"), (double)num_steps);
                        const double* energies = opt_energies;
                        const char* y_axis_label = "Relative Energy [kJ/mol]";
                        const char* plot_label = "OPT";

                        uint64_t energy_hash = md_hash64((const uint8_t*)energies, sizeof(double) * num_steps, 0);

                        if (energy_hash != opt.energy_hash) {
                            opt.energy_hash = energy_hash;
                            ImPlot::SetNextAxesToFit();
                        }

                        // Find minima as the reference energy
                        double ref_energy = energies[0];
                        if (opt_type == OPT_IRC) {
                            if (ts_index < num_steps) {
                                ref_energy = energies[ts_index];
                            }
                            plot_label = "Intrinsic Reaction Coordinate (IRC)";
                        } else {
                            for (size_t i = 1; i < num_steps; i++) {
                                ref_energy = MD_MIN(ref_energy, energies[i]);
                            }
                        }

                        double* x_vals = md_temp_alloc_array(temp, double, num_steps);
                        double* y_vals = md_temp_alloc_array(temp, double, num_steps);

                        for (size_t i = 0; i < num_steps; ++i) {
                            y_vals[i] = (energies[i] - ref_energy) * HARTREE_TO_KJ_PER_MOL;
                            x_vals[i] = (double)(i + 1);
                        }

                        const ImPlotFlags plot_flags = ImPlotFlags_NoMouseText;
                        if (ImPlot::BeginPlot(plot_label, ImVec2(-1, 0), plot_flags)) {
                            
                            ImPlot::SetupAxes("Step", y_axis_label);
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast);

                            if (opt_type == OPT_CONSTRAINED) {
                                ImPlot::PlotLineSpline("Energy", x_vals, y_vals, (int)num_steps);
                            } else {
                                ImPlot::PlotLine("Energy", x_vals, y_vals, (int)num_steps);
                            }

                            if (opt.selected >= 0 && opt.selected < (int)num_steps) {
                                double val = x_vals[opt.selected];
                                ImPlot::PlotInfLines("##Selected", &val, 1);
                            }

                            plot_peaks("##opt_peaks", x_vals, y_vals, num_steps, opt.selected, opt.hovered, PlotPeaksFlags_Points);

                            if (opt_type == OPT_IRC) {
                                // Plot a scatter point for the ts_index
                                const double x = x_vals[ts_index];
                                const double y = y_vals[ts_index];

                                ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 6.0f, ImVec4(0,0,0,-1), 2.0f);
                                ImPlot::SetNextFillStyle(ImVec4(0, 0, 0, -1), 0.5f);

                                ImPlot::PlotScatter("Transition State", &x, &y, 1);
                            }

                            if (ImPlot::IsPlotHovered() && opt.hovered != -1) {
                                ImGui::SetTooltip("Step %d\nEnergy: %.8f kJ/mol", opt.hovered + 1, y_vals[opt.hovered]);
                            }

                            ImPlot::EndPlot();
                        }
                        ImGui::Spacing();
                        int value = opt.selected + 1;
                        //ImGui::SetNextItemWidth(-FLT_MIN);
						const char* step_label = "Step";
						const ImVec2 label_size = ImGui::CalcTextSize(step_label);
						const float width = label_size.x + ImGui::GetStyle().ItemInnerSpacing.x;
                        ImGui::SetNextItemWidth(-width);
                        ImGui::SliderInt("Step", &value, 0, (int)num_steps);
                        opt.selected = value - 1;

                        ImGui::TreePop();
                    }
                    else {
                        opt.selected = (int)num_steps - 1;
                    }

                    static int prev_idx = -1;
                    if (opt.selected != prev_idx) {
                        prev_idx = opt.selected;
                        set_atom_coordinates(state, attribute_vec3_row(sys, MD_STR_LIT("vlx/opt/coordinate"), (size_t)opt.selected, qm.count));
                    }
                }
            }

            size_t num_atoms = qm.count;
            if (num_atoms > 0 && ImGui::TreeNode("Geometry")) {
                static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollX |
                                                        ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;

                static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_NoSort;

                if (ImGui::BeginTable("Geometry Table", 5, flags, ImVec2(500, -1), 0)) {
                    const dvec3_t* opt_coord  = attribute_vec3_row(sys, MD_STR_LIT("vlx/opt/coordinate"), (size_t)opt.selected, num_atoms);
                    const dvec3_t* atom_coord = opt_coord ? opt_coord : qm.coordinate;
                    const uint8_t* atom_nr    = qm.atomic_number;

                    ImGui::TableSetupColumn("Atom", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Symbol", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Coord X", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Coord Y", columns_base_flags, 0.0f);
                    ImGui::TableSetupColumn("Coord Z", columns_base_flags | ImGuiTableColumnFlags_WidthFixed, 0.0f);
                    ImGui::TableSetupScrollFreeze(0, 1);
                    ImGui::TableHeadersRow();

                    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                    ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                    bool item_hovered = false;
                    for (size_t row_n = 0; row_n < num_atoms; row_n++) {

                        ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        bool is_sel = md_bitmask_test_bit(&state.selection.selection_mask, row_n); //If atom is selected, mark it as such
                        bool is_hov = md_bitmask_test_bit(&state.selection.highlight_mask, row_n); //If atom is hovered,  mark it as such
                        ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);
                        ImGui::TableNextColumn();

                        if (is_hov) {
                            ImGui::PushStyleColor(ImGuiCol_Header, IM_YELLOW);
                        }
                        else {
                            ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                        }

                        char label[24];
                        snprintf(label, sizeof(label), "%zu", row_n + 1);
                        ImGui::Selectable(label, is_sel || is_hov, selectable_flags);
                        if (ImGui::TableGetHoveredRow() == (int)(row_n + 1)) {
                            if (state.mold.sys.atom.count > row_n) {
                                md_bitmask_clear(&state.selection.highlight_mask);
                                md_bitmask_or(&state.selection.highlight_mask, md_bitset_range((uint32_t)(row_n), (uint32_t)(row_n) + 1));
                                item_hovered = true;

                                //Selection
                                if (ImGui::IsKeyDown(ImGuiKey_MouseLeft) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                    md_bitmask_or(&state.selection.selection_mask, md_bitset_range((uint32_t)(row_n), (uint32_t)(row_n) + 1));
                                }
                                //Deselect
                                else if (ImGui::IsKeyDown(ImGuiKey_MouseRight) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                    md_bitmask_andnot(&state.selection.selection_mask, md_bitset_range((uint32_t)(row_n), (uint32_t)(row_n) + 1));
                                }
                            }
                        }

                        // Both columns are optional members of the QM atom domain (see es_qm_atoms).
                        ImGui::TableNextColumn();
                        if (atom_nr) {
                            md_str_t sym = md_util_element_symbol(atom_nr[row_n]);
                            ImGui::Text(MD_STR_FMT, MD_STR_ARG(sym));
                        } else {
                            ImGui::TextUnformatted("-");
                        }
                        for (int c = 0; c < 3; ++c) {
                            ImGui::TableNextColumn();
                            if (atom_coord) {
                                ImGui::Text("%12.6f", atom_coord[row_n].elem[c]);
                            } else {
                                ImGui::TextUnformatted("-");
                            }
                        }

                        ImGui::PopStyleColor(1);
                                
                    }
                    if (!item_hovered && ImGui::IsWindowHovered()) {
                        //Makes sure that we clear the highlight if we are in this window, but don't hover an item
                        md_bitmask_clear(&state.selection.highlight_mask);
                    }

                    ImGui::PopStyleColor(2);
                    ImGui::EndTable();
                }
                ImGui::TreePop();
            }
        }
        ImGui::End();
    }

    typedef struct {
        double har_freq;
        double redu_mass;
        double force_const;
        double ir_intens;
        double* x;
        double* y;
        double* z;
    } vibration_mode;

    void draw_rsp_spectra_export_window(ApplicationState& state) {
#if 0
        (void)state;
        MD_ASSERT(&state);

        struct ExportFormat {
            md_str_t lbl;
            md_str_t ext;
        };

        ExportFormat table_formats[]{
            {MD_STR_INIT("XVG"), MD_STR_INIT("xvg")},
            {MD_STR_INIT("CSV"), MD_STR_INIT("csv")}
        };

        struct ExportProperty {
            double* x = 0;
            double* y = 0;
            md_str_t lable;
            md_str_t y_unit;
        };

        ExportProperty properties[]{
            {rsp.x_spectra, rsp.y_spectra_abs,   MD_STR_INIT("Absorption"), MD_STR_INIT((const char*)u8"ε (L mol⁻¹ cm⁻¹)")},
            {rsp.x_spectra, rsp.y_spectra_ecd,   MD_STR_INIT("ECD"),        MD_STR_INIT((const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)")},
            {rsp.vib_x,          rsp.vib_y, MD_STR_INIT("Vibration"),  MD_STR_INIT("IR Intensity (km/mol)")}
        };
        
        if (ImGui::Begin("Spectra Export", &rsp.show_export_window)) {
            static int table_format = 0;
            static int property_idx = 0;
            const char* x_unit = "DEBUG";
            
            //TODO: Add sanity checks
            
            ImGui::PushItemWidth(200);

            md_str_t file_extension = {};
            if (ImGui::BeginCombo("File Format", table_formats[table_format].lbl.ptr)) {
                for (int i = 0; i < (int)MD_COUNTOF(table_formats); ++i) {
                    if (ImGui::Selectable(table_formats[i].lbl.ptr, table_format == i)) {
                        table_format = i;
                    }
                }
                ImGui::EndCombo();
            }
            file_extension = table_formats[table_format].ext;

            if (ImGui::BeginCombo("Property", properties[property_idx].lable.ptr)) {
                for (int i = 0; i < (int)MD_COUNTOF(properties); ++i) {
                    if (ImGui::Selectable(properties[i].lable.ptr, property_idx == i)) {
                        property_idx = i;
                    }
                }
                ImGui::EndCombo();
            }

            if (property_idx == 0 || property_idx == 1) {
                x_unit = x_unit_full_str[rsp.x_unit];
            }
            else if (property_idx == 2) {
                x_unit = (const char*)u8"Harmonic Frequency (cm⁻¹)";
            }

            static bool export_valid = true;
            bool export_clicked = ImGui::Button("Export");
            if (export_clicked) {
                export_valid = true;

                if (export_valid) {
                    char path_buf[1024];
                    md_array(const float*)  column_data = 0;
                    md_array(md_str_t)   column_labels = 0;
                    md_array(md_str_t)         legends = 0;

                    md_array(float) x_values = md_array_create(float, 1024, arena);
                    md_array(float) y_values = md_array_create(float, 1024, arena);

                    for (size_t i = 0; i < 1024; i++) {
                        x_values[i] = (float)properties[property_idx].x[i];
                        y_values[i] = (float)properties[property_idx].y[i];
                    }

                    //md_str_t x_label = md_str_from_cstr(x_unit);
                    md_array_push(column_labels, md_str_from_cstr(x_unit), arena);
                    //md_str_t y_label = md_str_from_cstr("Y");
                    md_array_push(column_labels, properties[property_idx].y_unit, arena);

                    md_array_push(column_data, x_values, arena);
                    md_array_push(column_data, y_values, arena);

                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, file_extension)) {
                        md_str_t path = { path_buf, strnlen(path_buf, sizeof(path_buf)) };
                        if (table_format == 0) {
                            //md_xvg_write
                            //TODO: Implement md_xvg_write_to_file
                            md_str_t header = md_xvg_format_header(properties[property_idx].lable, md_str_from_cstr(x_unit), properties[property_idx].y_unit, 0, legends, arena);
                            md_str_t xvg = md_xvg_format(header, 2, 1024, column_data, arena);
                            md_file_t file = { 0 };
                            if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                                const size_t written_bytes = md_file_write(file, xvg.ptr, xvg.len);
                                if (written_bytes != xvg.len) {
                                    MD_LOG_ERROR("CSV: Unexpected error, some bytes were not written");
                                }
                                md_file_close(&file);
                            }
                            else {
                                MD_LOG_ERROR("CSV: File could not be opened for writing: '%.*s'", (int)path.len, path.ptr);
                            }

                        }
                        else if (table_format == 1) {
                            //export_csv(column_data, column_labels, 0, 0, path);
                            md_csv_write_to_file(column_data, column_labels, 2, 1024, path);
                        }
                    }
                }
            }
            ImGui::PopItemWidth();
            if (!export_valid) { ImGui::TextWrapped("Values are not valid, make sure that you have opened plot windows once to trigger calculations"); }
        }
        ImGui::End();
#endif
    }

    template<typename T>
    static inline void find_min_max(int* min_idx, int* max_idx, const T* arr, size_t n) {
        if (min_idx) *min_idx = 0;
        if (max_idx) *max_idx = 0;
        for (size_t i = 1; i < n; i++) {
            if (min_idx) {
                if (arr[i] < arr[*min_idx]) {
                    *min_idx = (int)i;
                }
            }
            if (max_idx) {
                if (arr[i] > arr[*max_idx]) {
                    *max_idx = (int)i;
                }
            }
        }
    }

    void draw_rsp_window(ApplicationState& state) {
        md_temp_scope_t temp = md_temp_begin_in(state.allocator.frame);
        MD_DEFER { md_temp_end(temp); };

        if (!rsp.show_window) return;

        const md_system_t& sys = state.mold.sys;

		const RspType rsp_type = response_type(sys);
        // The frequencies are given in a.u. but can either be peaks to broaden if LINEAR, or spectrum samples if CPP.
        // Its own length is the frequency count - there is no separate number to disagree with it.
        size_t num_frequencies = 0;
        const double* x_freq_au = attribute_series_f64(&num_frequencies, sys, MD_STR_LIT("vlx/rsp/frequency"));
        if (!x_freq_au) num_frequencies = 0;

        const size_t num_normal_modes = attribute_series_count(sys, MD_STR_LIT("vlx/vib/frequency"));
        // If the spetrum is broadened from peaks, this is the number of samples in the visual space (x_min -> x_max) used to sample that signal.
        constexpr int num_broadened_samples = 2048;

        // XPS is a delta-SCF property and can be the only thing in the file, so it has to be able to
        // keep this window alive on its own. Its ionization energies are the column everything else
        // in that section is indexed against, so their presence is the question.
        const bool has_xps = attribute_exists(sys, MD_STR_LIT("vlx/xps/ionization_energy"));

        if (num_frequencies == 0 && num_normal_modes == 0 && !has_xps) return;

        double x_freq_au_min = DBL_MAX;
        double x_freq_au_max = -DBL_MAX;
        for (size_t i = 0; i < num_frequencies; i++) {
            x_freq_au_min = MD_MIN(x_freq_au_min, x_freq_au[i]);
            x_freq_au_max = MD_MAX(x_freq_au_max, x_freq_au[i]);
        }

        ImPlot::PushStyleVar(ImPlotStyleVar_FitPadding, ImVec2(0.1f, 0.1f));
        MD_DEFER { ImPlot::PopStyleVar(); };

        // Hover hand-off between the spectra and their state tables. Each pairing runs through a
        // single index (rsp.hovered for the response spectra, vib.hovered for the vibrational ones);
        // whichever widget the cursor is over writes it as the window is drawn, and everything else
        // just reads it. Peak -> row lands in the same frame, since the plots are drawn before the
        // table; row -> peak lands on the next one.
        //
        // A table that has the cursor CLAIMS its index: the tables scroll, so they are child windows,
        // and the reset at the end of this window would otherwise see the outer window as un-hovered
        // and drop a hover the table only just resolved - which is what kept row -> peak from ever
        // working. Claiming per index rather than per window also keeps one table's cursor from
        // freezing a stale hover in the other one.
        bool rsp_hover_claimed = false;
        bool vib_hover_claimed = false;

        ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Response", &rsp.show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    if (ImGui::MenuItem("Export")) {
                        rsp.show_export_window = true;
                    }
                    ImGui::EndMenu();
                }
                if (num_normal_modes > 0) {
                    if (ImGui::BeginMenu("Settings")) {
                        ImGui::Checkbox("VIB: Displace Atomic Orbitals", &vib.displace_aos);
                        ImGui::EndMenu();
                    }
                }
                ImGui::EndMenuBar();
            }

            struct SpectrumGetterData {
                const double* x_peaks;
                const double* y_peaks;
                size_t num_peaks;
                double x_min;
                double x_max;
                double broadening_fwhm;
                size_t num_samples;
                x_unit_t x_unit = X_UNIT_EV;
            };

            const ImGuiTreeNodeFlags tree_flags = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_NoTreePushOnOpen;

            // RIXS is its own thing: the 2D map has no meaningful representation in the shared
            // 1D spectrum plots below, so it gets its own section.
            if (rsp_type == RSP_RIXS) {
                RixsMapInput rixs_input = {};
                if (rixs_map_input_gather(rixs_input, state.mold.sys)) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                    if (ImGui::TreeNodeEx("RIXS", tree_flags)) {
                        draw_rixs_section(rixs_input, rixs);
                    }
                }
            }

            // XPS is independent of the response type, so it sits outside the electronic spectroscopy
            // section below (which is driven by the rsp eigenvalues).
            if (has_xps) {
                ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                if (ImGui::TreeNodeEx("XPS", tree_flags)) {
                    draw_xps_plot(state);
                }
            }

            // Hide the ECD plot if the RSP is a VIB calculation, since ECD is not relevant for VIB, but can be involved as an intermediate result holder.
            if (num_normal_modes == 0 && num_frequencies > 0 && ImGui::TreeNode("Electronic Spectroscopy")) {
                // Explicit x axis limits for linking
                static double x_min = 0.0;
                static double x_max = 10.0;

                const float avail_width = ImGui::GetContentRegionAvail().x;
                ImGui::PushItemWidth(MD_MIN(avail_width, 200));
                // Broadening min max values
                static const double broadening_min = 0.05;
                static const double broadening_max = 1.0;

                ImGui::SliderScalar((const char*)u8"Broadening FWHM (eV)", ImGuiDataType_Double, &rsp.broadening_fwhm, &broadening_min, &broadening_max);
                ImGui::Combo("Broadening mode", (int*)(&rsp.broadening_mode), broadening_mode_str, BROADENING_MODE_COUNT);
                ImGui::Combo("X unit", (int*)(&rsp.x_unit), x_unit_full_str, X_UNIT_COUNT);
                ImGui::PopItemWidth();

                // Every y column below is indexed over the frequency axis, so one shorter than it is
                // treated as absent rather than read past its end.
                // Samples are only used in CPP RSP, these hold intermediate results, i.e. converted to the selected x unit.
                size_t num_samples = 0;
                double* x_samples = NULL;
                const double* y_samples_sigma = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/rsp/cpp/sigma"), num_frequencies);
                const double* y_samples_delta_epsilons = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/rsp/cpp/delta_epsilon"), num_frequencies);
                const double* y_samples_ord = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/rsp/cpp/optical_rotation"), num_frequencies);
                const double* y_samples_cs = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/rsp/tpa/cross_section"), num_frequencies);

                // Peaks are only used in LINEAR RSP, not in CPP
				size_t num_peaks = 0;
				double* x_peaks = NULL;
                const double* y_peaks_osc = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/rsp/oscillator_strength"), num_frequencies);
                const double* y_peaks_cgs = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/rsp/rotatory_strength"), num_frequencies);
                const double* y_peaks_tpa_trans_linear = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/rsp/tpa/linear"), num_frequencies);
                const double* y_peaks_tpa_trans_circular = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/rsp/tpa/circular"), num_frequencies);

                bool refit = false;

                // Construct a hash to compare against if we need to refit.
                uint64_t hash = md_hash64(&rsp.broadening_fwhm, sizeof(rsp.broadening_fwhm), rsp.broadening_mode);
                hash = md_hash64(x_freq_au, num_frequencies * sizeof(double), hash);
                hash = md_hash64(&rsp.x_unit, sizeof(rsp.x_unit), hash);

                if (hash != rsp.hash) {
                    rsp.hash = hash;
                    refit = true;
                }

				if (rsp_type == RSP_LINEAR || rsp_type == RSP_TPA_TRANSITION) {
					num_peaks = num_frequencies;
                    x_peaks = md_temp_alloc_array(temp, double, num_peaks);
                    convert_values_from_au(x_peaks, x_freq_au, num_peaks, rsp.x_unit);
				} else {
                    num_samples = num_frequencies;
                    x_samples = md_temp_alloc_array(temp, double, num_samples);
                    convert_values_from_au(x_samples, x_freq_au, num_samples, rsp.x_unit);
                }

                if ((num_samples > 0 && y_samples_sigma) || (num_peaks > 0 && y_peaks_osc)) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                    if (ImGui::TreeNodeEx("Absorption Plot", tree_flags)) {
                        ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                        if (refit) {
                            ImPlot::SetNextAxesToFit();
                        }

                        auto getter_lorentzian = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            // Convert x to a.u. for the broadening function
                            double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            double y = lorentzian_abs(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                            double eps = x_au * y * OSCILLATOR_STRENGTH_TO_EPSILON;

                            return ImPlotPoint{x, eps};
                        };

                        auto getter_gaussian = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            // Convert x to a.u. for the broadening function
                            double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            double y = gaussian_abs(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                            double eps = x_au * y * OSCILLATOR_STRENGTH_TO_EPSILON;

                            return ImPlotPoint{ x, eps };
                        };

                        ImPlotGetter getter = (rsp.broadening_mode == BROADENING_MODE_LORENTZIAN) ? getter_lorentzian : getter_gaussian;

                        double* y_samples = md_temp_alloc_zero_array(temp, double, num_samples);

                        if (ImPlot::BeginPlot("Absorption")) {
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                            ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                            ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"ε (L mol⁻¹ cm⁻¹)");
                            if (num_peaks > 0) {
                                ImPlot::SetupAxis(ImAxis_Y2, "f", ImPlotAxisFlags_AuxDefault);
                            }
                            ImPlot::SetupFinish();

                            ImPlot::SetAxis(ImAxis_Y1);
                            if (rsp_type == RSP_CPP) {
                                // Convert y_samples to eps
                                constexpr double sigma_to_epsilon = 7323.816924863764;
                                for (size_t i = 0; i < num_samples; ++i) {
                                    y_samples[i] = y_samples_sigma[i] * sigma_to_epsilon;
                                }
                                ImPlot::PlotLineSpline("Interpolated Spectrum", x_samples, y_samples, (int)num_samples);
                                int selected = -1;
                                int hovered = -1;
                                plot_peaks("Samples", x_samples, y_samples, num_samples, selected, hovered, PlotPeaksFlags_Points);

                                if (ImPlot::IsPlotHovered() && hovered != -1) {
                                    double x_val = x_samples[hovered];
                                    double y_val = y_samples[hovered];
                                    const char* x_unit = x_unit_short_str[rsp.x_unit];
                                    const char* y_unit = (const char*)u8"ε (L mol⁻¹ cm⁻¹)";
                                    if (ImGui::BeginTooltip()) {
                                        ImGui::Text("%f %s, %f %s", x_val, x_unit, y_val, y_unit);
                                        ImGui::EndTooltip();
                                    }
                                }
                            } else {
                                SpectrumGetterData data = { x_freq_au, y_peaks_osc, num_peaks, x_min, x_max, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                ImPlot::PlotLineG("Broadened Spectrum", getter, &data, num_broadened_samples);
                            }

                            if (num_peaks > 0 && x_peaks && y_peaks_osc) {
                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("Oscillator Strength", x_peaks, y_peaks_osc, num_peaks, rsp.selected, rsp.hovered);

                                if (ImPlot::IsPlotHovered() && rsp.hovered != -1) {
                                    double x_val = x_peaks[rsp.hovered];
                                    double y_val = y_peaks_osc[rsp.hovered];
                                    int state_idx = rsp.hovered + 1;
                                    const char* x_label = x_unit_label_str[rsp.x_unit];
                                    const char* x_unit = x_unit_short_str[rsp.x_unit];
                                    if (ImGui::BeginTooltip()) {
                                        ImGui::Text("State %i: %s = %f (%s), f = %f", state_idx, x_label, x_val, x_unit, y_val);
                                        ImGui::EndTooltip();
                                    }
                                }
                            }

                            if (ImPlot::FitThisFrame()) {
                                ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                if (plot) {
                                    const double min_x = convert_value_from_au(x_freq_au_min, rsp.x_unit);
                                    const double max_x = convert_value_from_au(x_freq_au_max, rsp.x_unit);
                                    plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(MD_MIN(min_x, max_x), MD_MAX(min_x, max_x));

                                    // FIND MAX Y AXES
                                    if (rsp_type == RSP_CPP) {
                                        int max_idx;
                                        find_min_max(NULL, &max_idx, y_samples, num_samples);
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, y_samples[max_idx]);
                                    } else {
                                        // Find the maximum peak and calculate the maximum value over the spectrum.
                                        int max_peak_idx = 0;
                                        find_min_max(NULL, &max_peak_idx, y_peaks_osc, num_peaks);
                                        double max_y = y_peaks_osc[max_peak_idx];

                                        double max_eps = 0.0;
                                        SpectrumGetterData data = { x_freq_au, y_peaks_osc, num_peaks, min_x, max_x, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                        for (int i = 0; i < num_broadened_samples; ++i) {
                                            ImPlotPoint point = getter(i, &data);
                                            max_eps = MD_MAX(max_eps, point.y);
                                        }
                                        
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, max_eps);
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(0.0, max_y);
                                    }
                                }
                            }
                            ImPlot::EndPlot();
                        }
                    }
                }
                
                if ((num_samples > 0 && y_samples_delta_epsilons) || (num_peaks > 0 && y_peaks_cgs)) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);
                    if (ImGui::TreeNodeEx("ECD Plot", tree_flags)) {
                        ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                        if (refit) {
                            ImPlot::SetNextAxesToFit();
                        }

                        auto getter_lorentzian = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            // Convert x to a.u. for the broadening function
                            double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            double y = lorentzian_ecd(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                            double delta_eps = x_au * y * ROTATORY_STRENGTH_TO_DELTA_EPSILON;

                            return ImPlotPoint{x, delta_eps};
                        };

                        auto getter_gaussian = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            // Convert x to a.u. for the broadening function
                            double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            double y = gaussian_ecd(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                            double delta_eps = x_au * y * ROTATORY_STRENGTH_TO_DELTA_EPSILON;

                            return ImPlotPoint{ x, delta_eps };
                        };

                        ImPlotGetter getter = (rsp.broadening_mode == BROADENING_MODE_LORENTZIAN) ? getter_lorentzian : getter_gaussian;

                        if (ImPlot::BeginPlot("ECD")) {
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                            ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                            ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)");
                            if (num_peaks > 0) {
                                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"R (10⁻⁴⁰ cgs)", ImPlotAxisFlags_AuxDefault);
                            }
                            ImPlot::SetupFinish();

                            ImPlot::SetAxis(ImAxis_Y1);
                            if (rsp_type == RSP_CPP) {
                                ImPlot::PlotLineSpline("Interpolated Spectrum", x_samples, y_samples_delta_epsilons, (int)num_samples);

                                int selected = -1;
                                int hovered = -1;
                                plot_peaks("Samples", x_samples, y_samples_delta_epsilons, num_samples, selected, hovered, PlotPeaksFlags_Points);

                                if (ImPlot::IsPlotHovered() && hovered != -1) {
                                    double x_val = x_samples[hovered];
                                    double y_val = y_samples_delta_epsilons[hovered];
                                    const char* x_unit = x_unit_short_str[rsp.x_unit];
                                    const char* y_unit = (const char*)u8"Δε(ω) (L mol⁻¹ cm⁻¹)";
                                    if (ImGui::BeginTooltip()) {
                                        ImGui::Text("%f %s, %f %s", x_val, x_unit, y_val, y_unit);
                                        ImGui::EndTooltip();
                                    }
                                }
                            } else {
                                SpectrumGetterData data = { x_freq_au, y_peaks_cgs, num_peaks, x_min, x_max, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                ImPlot::PlotLineG("Broadened Spectrum", getter, &data, num_broadened_samples);
                            }

                            if (num_peaks > 0 && x_peaks && y_peaks_cgs) {
                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("Rotatory Strength", x_peaks, y_peaks_cgs, num_peaks, rsp.selected, rsp.hovered);

                                if (ImPlot::IsPlotHovered() && rsp.hovered != -1) {
                                    double x_val = x_peaks[rsp.hovered];
                                    double y_val = y_peaks_cgs[rsp.hovered];
                                    int state_idx = rsp.hovered + 1;
                                    const char* x_label = x_unit_label_str[rsp.x_unit];
                                    const char* x_unit  = x_unit_short_str[rsp.x_unit];
                                    const char* y_unit  = (const char*)u8"10⁻⁴⁰ cgs";
                                    if (ImGui::BeginTooltip()) {
                                        ImGui::Text((const char*)u8"State %i: %s = %f (%s), R = %f (%s)", state_idx, x_label, x_val, x_unit, y_val, y_unit);
                                        ImGui::EndTooltip();
                                    }
                                }
                            }

                            if (ImPlot::FitThisFrame()) {
                                ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                if (plot) {
                                    const double min_x = convert_value_from_au(x_freq_au_min, rsp.x_unit);
                                    const double max_x = convert_value_from_au(x_freq_au_max, rsp.x_unit);
                                    plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(MD_MIN(min_x, max_x), MD_MAX(min_x, max_x));

                                    // FIND MAX Y AXES
                                    if (rsp_type == RSP_CPP) {
                                        int min_idx, max_idx;
                                        find_min_max(&min_idx, &max_idx, y_samples_delta_epsilons, num_samples);
                                        double max_delta_eps = MD_MAX(fabs(y_samples_delta_epsilons[min_idx]), fabs(y_samples_delta_epsilons[max_idx]));
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(-max_delta_eps, max_delta_eps);
                                    } else {
                                        // Find the maximum peak and calculate the maximum value over the spectrum.
                                        int min_peak_idx, max_peak_idx;
                                        find_min_max(&min_peak_idx, &max_peak_idx, y_peaks_cgs, num_peaks);
                                        int peak_idx = (fabs(y_peaks_cgs[min_peak_idx]) > fabs(y_peaks_cgs[max_peak_idx])) ? min_peak_idx : max_peak_idx;
                                        double max_y = fabs(y_peaks_cgs[peak_idx]);

                                        double max_delta_eps = 0.0;
                                        SpectrumGetterData data = { x_freq_au, y_peaks_cgs, num_peaks, min_x, max_x, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                        for (int i = 0; i < num_broadened_samples; ++i) {
                                            ImPlotPoint point = getter(i, &data);
                                            max_delta_eps = MD_MAX(max_delta_eps, fabs(point.y));
                                        }
                                        
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(-max_delta_eps, max_delta_eps);
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(-max_y, max_y);
                                    }
                                }
                            }
                            ImPlot::EndPlot();
                        }
                    }
                }

                if (num_samples > 0 && y_samples_ord) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Appearing);

                    if (ImGui::TreeNodeEx("ORD Plot", tree_flags)) {
                        ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                        if (refit) {
                            ImPlot::SetNextAxesToFit();
                        }

                        if (ImPlot::BeginPlot("Optical Rotatory Dispersion")) {
                            const char* y_unit = (const char*)u8"ORD (10³ ° cm^3 g⁻¹ dm⁻¹)";
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                            ImPlot::SetupAxis(ImAxis_X1, x_unit_full_str[rsp.x_unit]);
                            ImPlot::SetupAxis(ImAxis_Y1, y_unit);
                            ImPlot::SetupFinish();

                            ImPlot::SetAxis(ImAxis_Y1);
                            ImPlot::PlotLineSpline("Interpolated Spectrum", x_samples, y_samples_ord, (int)num_samples);
                            int selected = -1;
                            int hovered = -1;
                            plot_peaks("Samples", x_samples, y_samples_ord, num_samples, selected, hovered, PlotPeaksFlags_Points);

                            if (ImPlot::IsPlotHovered() && hovered != -1) {
                                double x_val = x_samples[hovered];
                                double y_val = y_samples_ord[hovered];
                                const char* x_unit = x_unit_short_str[rsp.x_unit];                                    
                                if (ImGui::BeginTooltip()) {
                                    ImGui::Text("%f %s, %f %s", x_val, x_unit, y_val, y_unit);
                                    ImGui::EndTooltip();
                                }
                            }

                            if (ImPlot::FitThisFrame()) {
                                ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                if (plot) {
                                    const double min_x = convert_value_from_au(x_freq_au_min, rsp.x_unit);
                                    const double max_x = convert_value_from_au(x_freq_au_max, rsp.x_unit);
                                    plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(MD_MIN(min_x, max_x), MD_MAX(min_x, max_x));

                                    int min_idx, max_idx;
                                    find_min_max(&min_idx, &max_idx, y_samples_ord, num_samples);
                                    double max_ord = MD_MAX(fabs(y_samples_ord[min_idx]), fabs(y_samples_ord[max_idx]));
                                    plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(-max_ord, max_ord);
                                }
                            }
                            ImPlot::EndPlot();
                        }
                    }
                }

                if ((num_samples > 0 && y_samples_cs) || (num_peaks > 0 && y_peaks_tpa_trans_linear)) {
                    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
                    if (ImGui::TreeNodeEx("Two Photon Absorption", tree_flags)) {

                        char x_label[64];
                        snprintf(x_label, sizeof(x_label), "Photon Energy [%s]", x_unit_short_str[rsp.x_unit]);
                        const char* title = (rsp_type == RSP_TPA_TRANSITION) ? "TPA Transition" : "TPA";

                        auto getter = [](int idx, void* user_data) -> ImPlotPoint {
                            const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                            // Construct an x value for idx as a linear space between x_min and x_max.
                            double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                            const double x_au = convert_value_to_au(x, data->x_unit);

                            // Compute the corresponding y value by broadening the peaks at this x value.
                            // Perform broadening in a.u.
                            // Lorentzian broadening
                            const double factor = AU_TO_GM * x_au * x_au;
                            double y = factor * lorentzian_tpa(x_au, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);

                            return ImPlotPoint{x, y};
                        };

                        ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                        if (refit) {
                            ImPlot::SetNextAxesToFit();
                        }

                        if (ImPlot::BeginPlot(title)) {
                            ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                            ImPlot::SetupAxis(ImAxis_X1, x_label);
                            ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"TPA Cross-Section [GM]");
                            if (rsp_type == RSP_TPA_TRANSITION) {
                                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"TPA Strengths (linear) [a.u.]", ImPlotAxisFlags_AuxDefault);
                            }
                            ImPlot::SetupFinish();

                            ImPlot::SetAxis(ImAxis_Y1);
                            if (rsp_type == RSP_TPA_TRANSITION) {
                                SpectrumGetterData getter_data = { x_freq_au, y_peaks_tpa_trans_linear, num_peaks, x_min, x_max, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                ImPlot::PlotLineG("Broadened Spectrum", getter, &getter_data, num_broadened_samples);

                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("Linear TPA Strength", x_peaks, y_peaks_tpa_trans_linear, num_peaks, rsp.selected, rsp.hovered);
                            } else {
                                ImPlot::PlotLineSpline("Interpolated Spectrum", x_samples, y_samples_cs, (int)num_samples);
                                plot_peaks("Computed Cross-Sections", x_samples, y_samples_cs, num_samples, rsp.selected, rsp.hovered, PlotPeaksFlags_Points);
                            }

                            if (ImPlot::FitThisFrame() || refit) {
                                ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                if (plot) {
                                    const double min_x = convert_value_from_au(x_freq_au_min, rsp.x_unit);
                                    const double max_x = convert_value_from_au(x_freq_au_max, rsp.x_unit);
                                    plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(MD_MIN(min_x, max_x), MD_MAX(min_x, max_x));

                                    // FIND MAX Y AXES
                                    if (rsp_type == RSP_TPA) {
                                        int max_idx;
                                        find_min_max(NULL, &max_idx, y_samples_cs, num_samples);
                                        double max_y = y_samples_cs[max_idx];
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, max_y);
                                    } else {
                                        int max_idx;
                                        find_min_max(NULL, &max_idx, y_peaks_tpa_trans_linear, num_peaks);
                                        double max_y = y_peaks_tpa_trans_linear[max_idx];
                                        double max_cs = 0.0;
                                        SpectrumGetterData data = { x_freq_au, y_peaks_tpa_trans_linear, num_peaks, min_x, max_x, rsp.broadening_fwhm * EV_TO_HARTREE, num_broadened_samples, rsp.x_unit};
                                        for (int i = 0; i < num_broadened_samples; ++i) {
                                            ImPlotPoint point = getter(i, &data);
                                            max_cs = MD_MAX(max_cs, point.y);
                                        }
                                        
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, max_cs);
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(0.0, max_y);
                                    }
                                }
                            }

                            ImPlot::EndPlot();
                        }
                    }
                }

                if (num_peaks > 0 && ImGui::TreeNodeEx("State Table", tree_flags)) {
                    // Table
                    static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY |
                        ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Sortable;

                    static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_DefaultSort;

                    ImVec2 table_size = { 0, 0 };

                    struct Column {
                        const char* label;
                        ImGuiTableColumnFlags flags;
                        const double* data;
                    };

                    Column cols[8];
                    int num_cols = 0;

                    // Populate colums with existing data
                    if (x_peaks) {
                        cols[num_cols++] = { x_unit_full_str[rsp.x_unit], columns_base_flags, x_peaks };
                    }
                    if (y_peaks_osc) {
                        cols[num_cols++] = { "Osc. Str.", columns_base_flags, y_peaks_osc };
                    }
                    if (y_peaks_cgs) {
                        cols[num_cols++] = { (const char*)u8"R (10⁻⁴⁰ cgs)", columns_base_flags, y_peaks_cgs };
                    }
                    if (y_peaks_tpa_trans_linear) {
                        cols[num_cols++] = { (const char*)u8"TPA Str. (lin) [a.u.]", columns_base_flags, y_peaks_tpa_trans_linear };
                    }
                    if (y_peaks_tpa_trans_circular) {
                        cols[num_cols++] = { (const char*)u8"TPA Str. (cir) [a.u.]", columns_base_flags, y_peaks_tpa_trans_circular };
                    }

                    if (ImGui::BeginTable("rsp table", num_cols + 1, flags, table_size, 0)) {
                        ImGui::TableSetupColumn("Index", columns_base_flags);
                        for (int col = 0; col < num_cols; ++col) {
                            ImGui::TableSetupColumn(cols[col].label, cols[col].flags);
                        }
                        ImGui::TableSetupScrollFreeze(0, 1);
                        ImGui::TableHeadersRow();

                        // TODO: Add sorting to the table once the vibs data structure is properly defined

                        // The row Selectable spans all columns, so it resolves the row hover on its own -
                        // no post-pass over TableGetHoveredRow needed. All this has to do is drop the
                        // previous value while the cursor is in the table and let the row under it set the
                        // new one; landing on no row (header, gap below the rows) then means no hover.
                        if (ImGui::IsWindowHovered()) {
                            rsp_hover_claimed = true;
                            rsp.hovered = -1;
                        }

                        const ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        for (int row_n = 0; row_n < (int)num_frequencies; row_n++) {
                            bool is_sel = row_n == rsp.selected;
                            bool is_hov = row_n == rsp.hovered;

                            // The Selectable is passed 'is_sel || is_hov', so it fills with ImGuiCol_Header
                            // for a hover that came from the plot and with ImGuiCol_HeaderHovered for one the
                            // cursor made here. Both are set to the same colour for a given row, so a row lights
                            // up identically whichever end the hover came from. Hover beats selection, matching
                            // plot_peaks and the atom table.
                            ImGui::PushStyleColor(ImGuiCol_Header, is_hov ? IM_YELLOW : IM_BLUE);
                            ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);

                            ImGui::TableNextRow();
                            ImGui::TableSetColumnIndex(0);

                            char label[16];
                            snprintf(label, sizeof(label), "%i", row_n + 1);
                            if (ImGui::Selectable(label, is_sel || is_hov, selectable_flags)) {
                                rsp.selected = (rsp.selected == row_n) ? -1 : row_n;
                            }
                            if (ImGui::IsItemHovered()) {
                                rsp.hovered = row_n;
                            }
                            ImGui::TableNextColumn();
                            for (int col = 0; col < num_cols; ++col) {
                                ImGui::Text("%12.6f", *(cols[col].data + row_n));
                                ImGui::TableNextColumn();
                            }

                            ImGui::PopStyleColor(2);
                        }

                        ImGui::EndTable();
                    }
                }
                ImGui::TreePop();
            }

            if (num_normal_modes > 0) {
                ImGui::SetNextItemOpen(true, ImGuiCond_Once);
                if (ImGui::TreeNode("Vibrational Spectroscopy")) {
                    // The displacement animation below moves the QM atoms, so it needs the domain
                    // and the map across to the system's own atoms.
                    QmAtoms qm = {};
                    es_qm_atoms(&qm, sys);

                    size_t num_atoms = qm.count;

                    // Broadening gamma limits
                    static const double gamma_min = 1.0;
                    static const double gamma_max = 100.0;

                    // Frequency scaling factor limits
                    static const double scale_min = 0.85;
                    static const double scale_max = 1.05;

                    // Explicit x axis limits for linking
                    static double x_min = 0.0;
                    static double x_max = 10.0;

                    const float avail_width = ImGui::GetContentRegionAvail().x;
                    ImGui::PushItemWidth(MD_MIN(avail_width, 200));

                    bool refit = false;

                    refit |= ImGui::SliderScalar((const char*)u8"Broadening FWHM (cm⁻¹)", ImGuiDataType_Double, &vib.broadening_fwhm, &gamma_min, &gamma_max, "%.3f", ImGuiSliderFlags_Logarithmic);
                    refit |= ImGui::Combo("Broadening mode", (int*)(&vib.broadening_mode), broadening_mode_str, IM_ARRAYSIZE(broadening_mode_str));
                    refit |= ImGui::SliderScalar("Frequency scaling factor", ImGuiDataType_Double, &vib.freq_scaling_factor, &scale_min, &scale_max, "%.4f");
                    ImGui::PopItemWidth();

                    ImGui::SetItemTooltip("Frequency scaling factor");

                    auto spectrum_getter_lorentzian = [](int idx, void* user_data) -> ImPlotPoint {
                        const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                        // Construct an x value for idx as a linear space between x_min and x_max.
                        double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                        // Compute the corresponding y value by broadening the peaks at this x value.
                        double y = lorentzian_vib(x, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                        return ImPlotPoint{x, y};
                    };

                    auto spectrum_getter_gaussian = [](int idx, void* user_data) -> ImPlotPoint {
                        const SpectrumGetterData* data = (const SpectrumGetterData*)user_data;

                        // Construct an x value for idx as a linear space between x_min and x_max.
                        double x = data->x_min + (data->x_max - data->x_min) * idx / (data->num_samples - 1);

                        // Compute the corresponding y value by broadening the peaks at this x value.
                        double y = gaussian_vib(x, data->x_peaks, data->y_peaks, data->num_peaks, data->broadening_fwhm);
                        return ImPlotPoint{ x, y };
                    };

					ImPlotGetter getter = nullptr;
					if (vib.broadening_mode == BROADENING_MODE_LORENTZIAN) {
						getter = spectrum_getter_lorentzian;
					} else {
						getter = spectrum_getter_gaussian;
					}

                    const double* x_peaks_raw = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/vib/frequency"), num_normal_modes);
                    // Optional, and indexed over the normal modes - dropped when it does not cover them.
                    const double* y_peaks_ir = attribute_series_f64_covering(sys, MD_STR_LIT("vlx/vib/ir_intensity"), num_normal_modes);

                    size_t num_external_frequencies = 0;
                    const double* external_frequencies = attribute_series_f64(&num_external_frequencies, sys, MD_STR_LIT("vlx/vib/external_frequency"));
                    if (!external_frequencies) num_external_frequencies = 0;

                    // {E,M}: one row of per mode activities per external frequency, so a row IS the
                    // spectrum at that frequency and the pointers below are into the table's own
                    // storage rather than a gathered copy. The row count is checked against the
                    // frequency axis it is indexed by - two siblings that disagree are not a
                    // spectrum, and dropping both is better than plotting one against the other.
                    size_t num_raman_rows = 0;
                    size_t num_raman_cols = 0;
                    const double* raman_activities = attribute_matrix_f64(&num_raman_rows, &num_raman_cols, sys, MD_STR_LIT("vlx/vib/raman_activity"));
                    const double** y_raman_activities = NULL;
                    if (raman_activities && num_raman_rows == num_external_frequencies && num_raman_cols == num_normal_modes && num_external_frequencies > 0) {
                        y_raman_activities = md_temp_alloc_array(temp, const double*, num_external_frequencies);
                        for (size_t i = 0; i < num_external_frequencies; ++i) {
                            y_raman_activities[i] = raman_activities + i * num_raman_cols;
                        }
                    } else {
                        num_external_frequencies = 0;
                        external_frequencies = NULL;
                    }

                    double* x_peaks = md_temp_alloc_array(temp, double, num_normal_modes);
                    for (size_t i = 0; i < num_normal_modes; ++i) {
                        x_peaks[i] = x_peaks_raw[i] * vib.freq_scaling_factor;
                    }

                    uint64_t x_hash = md_hash64(x_peaks_raw, sizeof(double) * num_normal_modes, 0);
                    if (vib.hash != x_hash) {
                        vib.hash = x_hash;
                        refit = true;
                    }

				    const size_t num_samples = 512;

                    /*
                    // These do not need to be exposed as the user can right click the plot and invert the axes
                    ImGui::Checkbox("Invert X", &vib.invert_x);
                    ImGui::SameLine();
                    ImGui::Checkbox("Invert Y", &vib.invert_y);
                    */

                    ImPlotAxisFlags x_flags = vib.invert_x ? ImPlotAxisFlags_Invert : 0;
                    ImPlotAxisFlags y_flags = vib.invert_y ? ImPlotAxisFlags_Invert : 0;

                    if (x_peaks && y_peaks_ir && num_normal_modes) {
                        if (ImGui::TreeNodeEx("IR", tree_flags)) {

                            ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                            if (refit) {
                                ImPlot::SetNextAxesToFit();
                            }

                            if (ImPlot::BeginPlot("IR Spectrum")) {
                                ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);

                                ImPlot::SetupAxis(ImAxis_X1, (const char*)u8"Harmonic Frequency (cm⁻¹)", x_flags);
                                ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"Absorption Cross Section", y_flags);
                                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"IR Intensity (km/mol)", y_flags | ImPlotAxisFlags_AuxDefault);

                                ImPlot::SetupFinish();

                                ImPlot::SetAxis(ImAxis_Y1);
                                SpectrumGetterData getter_data = { x_peaks, y_peaks_ir, num_normal_modes, x_min, x_max, vib.broadening_fwhm, num_samples };
                                ImPlot::PlotLineG("Broadening", getter, &getter_data, num_samples);

                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("IR Intensity", x_peaks, y_peaks_ir, num_normal_modes, vib.selected, vib.hovered);

                                if (ImPlot::FitThisFrame()) {
                                    ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                    if (plot) {
                                        double min_x = DBL_MAX, max_x = -DBL_MAX, max_y = -DBL_MAX;
                                        for (size_t i = 0; i < num_normal_modes; ++i) {
                                            min_x = MD_MIN(min_x, x_peaks[i]);
                                            max_x = MD_MAX(max_x, x_peaks[i]);
                                            max_y = MD_MAX(max_y, y_peaks_ir[i]);
                                        }

                                        double normalization_factor = 1.0;
                                        if (vib.broadening_mode == BROADENING_MODE_LORENTZIAN) {
                                            normalization_factor = 1.0 / MD_PI;
                                        }
                                        else {
                                            normalization_factor = sqrt(2.0) / (sigma_from_fwhm(vib.broadening_fwhm) * sqrt(MD_PI));
                                        }

                                        plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(min_x, max_x);
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, 1.2 * max_y * normalization_factor);
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(0.0, 1.1 * max_y);
                                    }
                                }

                                ImPlot::EndPlot();
                            }
                        }
                    }

                    if (num_external_frequencies > 0) {
                        if (ImGui::TreeNodeEx("Raman Activity", tree_flags)) {
                            
                            if (num_external_frequencies > 1) {
                                ImGui::PushItemWidth(MD_MIN(avail_width, 200));

                                char label[32];
                                snprintf(label, sizeof(label), "%.4f", external_frequencies[vib.external_frequency_index]);
                                if (ImGui::BeginCombo("External Frequency", label)) {
                                    for (size_t i = 0; i < num_external_frequencies; ++i) {
                                        bool is_selected = ((size_t)vib.external_frequency_index == i);

                                        snprintf(label, sizeof(label), "%.4f", external_frequencies[i]);
                                        if (ImGui::Selectable(label, is_selected)) {
                                            vib.external_frequency_index = (int)i;
                                        }
                                        if (is_selected) {
                                            ImGui::SetItemDefaultFocus();
                                        }
                                    }
                                    ImGui::EndCombo();
                                }

                                ImGui::PopItemWidth();
                            }

                            ImPlot::SetNextAxisLinks(ImAxis_X1, &x_min, &x_max);

                            if (refit) {
                                ImPlot::SetNextAxesToFit();
                            }

                            const double* y_peaks_raman = y_raman_activities[vib.external_frequency_index];

                            char label[32];
                            double freq = external_frequencies[vib.external_frequency_index];
                            if (num_external_frequencies == 1 && freq == 0.0) {
                                snprintf(label, sizeof(label), "Raman Spectrum (static)");
                            } else {
                                snprintf(label, sizeof(label), "Raman Spectrum (%.4f)", freq);
                            }

                            if (ImPlot::BeginPlot(label)) {
                                ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                                ImPlot::SetupAxis(ImAxis_X1, (const char*)u8"Harmonic Frequency (cm⁻¹)", x_flags);
                                ImPlot::SetupAxis(ImAxis_Y1, (const char*)u8"Scattering Cross Section", y_flags);
                                ImPlot::SetupAxis(ImAxis_Y2, (const char*)u8"Raman Activity (Å⁴/amu)", y_flags | ImPlotAxisFlags_AuxDefault);
                                ImPlot::SetupFinish();

                                ImPlot::SetAxis(ImAxis_Y1);
								SpectrumGetterData data = { x_peaks, y_peaks_raman, num_normal_modes, x_min, x_max, vib.broadening_fwhm, num_broadened_samples };
                                ImPlot::PlotLineG("Broadening", getter, &data, num_broadened_samples);

                                ImPlot::SetAxis(ImAxis_Y2);
                                plot_peaks("Raman Activity", x_peaks, y_peaks_raman, num_normal_modes, vib.selected, vib.hovered);

                                if (ImPlot::FitThisFrame() || refit) {
                                    ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                                    if (plot) {
                                        double min_x = DBL_MAX, max_x = -DBL_MAX, max_y = -DBL_MAX;
                                        for (size_t i = 0; i < num_normal_modes; ++i) {
                                            min_x = MD_MIN(min_x, x_peaks[i]);
                                            max_x = MD_MAX(max_x, x_peaks[i]);
                                            max_y = MD_MAX(max_y, y_peaks_raman[i]);
                                        }

                                        plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(min_x, max_x);
                                        plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(0.0, 1.2 * max_y / (MD_PI * vib.broadening_fwhm));
                                        plot->Axes[ImAxis_Y2].FitExtents = ImPlotRange(0.0, 1.1 * max_y);
                                    }
                                }

                                ImPlot::EndPlot();
                            }
                        }
                    }

                    ImGui::PushItemWidth(MD_MIN(avail_width, 200));

                    int idx = vib.selected + 1;
                    ImGui::SliderInt("Index", &idx, 0, (int)num_normal_modes);
                    vib.selected = MD_CLAMP(idx - 1, -1, (int)num_normal_modes - 1);
                    ImGui::SliderFloat((const char*)"Amplitude", &vib.displacement_amp_scl, 0.25f, 5.0f);
                    ImGui::SetItemTooltip("Displacement Amplitude Scale");
                    ImGui::SliderFloat((const char*)"Frequency", &vib.displacement_freq_scl, 0.25f, 5.0f);
                    ImGui::SetItemTooltip("Displacement Frequency Scale");
                    ImGui::PopItemWidth();

                    // Table
                    static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY |
                                                         ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_Sortable;

                    static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_DefaultSort;

                    ImVec2 table_size = {0, 0};

                    int num_cols = 1 + (int)num_external_frequencies;
					if (x_peaks) {
						num_cols++;
					}
					if (y_peaks_ir) {
						num_cols++;
					}
                    if (ImGui::BeginTable("table", num_cols, flags, table_size, 0)) {
                        ImGui::TableSetupColumn("Vibration mode", columns_base_flags, 0.0f);
                        if (x_peaks) {
                            ImGui::TableSetupColumn("Harmonic Frequency", columns_base_flags, 0.0f);
                        }
                        if (y_peaks_ir) {
                            ImGui::TableSetupColumn("IR Intensity", columns_base_flags, 0.0f);
                        }
                        for (size_t i = 0; i < num_external_frequencies; ++i) {
                            char label[32];
                            double freq = external_frequencies[i];
                            snprintf(label, sizeof(label), "Raman Activity (%.4f)", freq);
                            ImGui::TableSetupColumn(label, columns_base_flags, 0.0f);
                        }
                        ImGui::TableSetupScrollFreeze(0, 1);
                        ImGui::TableHeadersRow();

                        // TODO: Add sorting to the table once the vibs data structure is properly defined

                        // Same hover contract as the response state table above.
                        if (ImGui::IsWindowHovered()) {
                            vib_hover_claimed = true;
                            vib.hovered = -1;
                        }

                        const ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        for (int row_n = 0; row_n < (int)num_normal_modes; row_n++) {
                            bool is_sel = row_n == vib.selected;
                            bool is_hov = row_n == vib.hovered;

                            // The Selectable is passed 'is_sel || is_hov', so it fills with ImGuiCol_Header
                            // for a hover that came from the plot and with ImGuiCol_HeaderHovered for one the
                            // cursor made here. Both are set to the same colour for a given row, so a row lights
                            // up identically whichever end the hover came from. Hover beats selection, matching
                            // plot_peaks and the atom table.
                            ImGui::PushStyleColor(ImGuiCol_Header, is_hov ? IM_YELLOW : IM_BLUE);
                            ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);

                            ImGui::TableNextRow();
                            ImGui::TableSetColumnIndex(0);

                            char label[16];
                            snprintf(label, sizeof(label), "%i", row_n + 1);
                            if (ImGui::Selectable(label, is_sel || is_hov, selectable_flags)) {
                                vib.selected = (vib.selected == row_n) ? -1 : row_n;
                            }
                            if (ImGui::IsItemHovered()) {
                                vib.hovered = row_n;
                            }
                            if (x_peaks) {
                                ImGui::TableNextColumn();
                                ImGui::Text("%12.6f", x_peaks[row_n]);
                            }
							if (y_peaks_ir) {
                                ImGui::TableNextColumn();
                                ImGui::Text("%12.6f", y_peaks_ir[row_n]);
                            }
                            for (size_t i = 0; i < num_external_frequencies; ++i) {
                                ImGui::TableNextColumn();
                                ImGui::Text("%12.6f", y_raman_activities[i][row_n]);
                            }

                            ImGui::PopStyleColor(2);
                        }

                        ImGui::EndTable();
                    }

                    const dvec3_t* atom_coord = qm.coordinate;

                    // Check selected state
                    if (vib.selected != -1) {
                        // Animate
                        vib.t += state.app.timing.delta_s * vib.displacement_freq_scl * 8.0;
                        const dvec3_t* norm_modes = attribute_vec3_row(sys, MD_STR_LIT("qm/atom/normal_mode"), (size_t)vib.selected, num_atoms);

                        // Both columns are needed: the displacement is added to the reference
                        // geometry, and qm/atom/coordinate is not guaranteed to be published.
                        if (norm_modes && atom_coord) {

#if 0
                            md_temp_scope_t temp = md_temp_begin();
                            MD_DEFER { md_temp_end(temp); };

                            const double scl = vib.displacement_amp_scl * 1.0;

                            float* vx = (float*)md_temp_alloc(temp, sizeof(float) * num_atoms);
                            float* vy = (float*)md_temp_alloc(temp, sizeof(float) * num_atoms);
                            float* vz = (float*)md_temp_alloc(temp, sizeof(float) * num_atoms);

                            for (size_t i = 0; i < num_atoms; i++) {
                                /*
                                state.mold.sys.atom.x[i] = (float)(atom_coord[i].x + norm_modes[i].x * scl);
                                state.mold.sys.atom.y[i] = (float)(atom_coord[i].y + norm_modes[i].y * scl);
                                state.mold.sys.atom.z[i] = (float)(atom_coord[i].z + norm_modes[i].z * scl);
                                */
                                vx[i] = (float)(norm_modes[i].x * scl);
                                vy[i] = (float)(norm_modes[i].y * scl);
                                vz[i] = (float)(norm_modes[i].z * scl);
                            }

                            //md_gl_mol_set_atom_position(state.mold.gl_mol, 0, num_atoms, state.mold.sys.atom.x, state.mold.sys.atom.y, state.mold.sys.atom.z, 0);
                            md_gl_mol_set_atom_velocity(state.mold.gl_mol, 0, num_atoms, vx, vy, vz, 0);
                            

#else
                            const double scl = vib.displacement_amp_scl * 0.25 * sin(vib.t);
                            // Through the QM atom map and bounded by the system state, exactly as
                            // set_atom_coordinates writes: QM atom i is system atom i only when the
                            // calculation covers the whole system.
                            for (size_t i = 0; i < num_atoms; i++) {
                                const size_t idx = qm_to_system_atom(qm, i);
                                if (idx >= state.mold.state.num_atoms) continue;
                                state.mold.state.xyz[idx].x = (float)(atom_coord[i].x + norm_modes[i].x * scl);
                                state.mold.state.xyz[idx].y = (float)(atom_coord[i].y + norm_modes[i].y * scl);
                                state.mold.state.xyz[idx].z = (float)(atom_coord[i].z + norm_modes[i].z * scl);
                            }
                            if (vib.displace_aos) {
                            }
                            state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
#endif
                            vib.coord_modified = true;

                            if (vib.displace_aos) {
                                viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemStateChanged, viamd::EventPayloadType_ApplicationState, &state);
                                size_t num_reps = md_array_size(state.representation.reps);
                                for (size_t i = 0; i < num_reps; ++i) {
                                    if (state.representation.reps[i].type == RepresentationType::ElectronicStructure) {
                                        flag_representation_as_dirty(&state.representation.reps[i]);
                                    }
                                }
                            }
                        }
                    }
                    // If all is deselected, reset coords once
                    else if (vib.coord_modified) {
                        // coord_modified is only ever set with atom_coord present, but the column
                        // is re-read every frame and the guard costs nothing.
                        for (size_t i = 0; atom_coord && i < num_atoms; i++) {
                            const size_t idx = qm_to_system_atom(qm, i);
                            if (idx >= state.mold.state.num_atoms) continue;
                            state.mold.state.xyz[idx].x = (float)atom_coord[i].x;
                            state.mold.state.xyz[idx].y = (float)atom_coord[i].y;
                            state.mold.state.xyz[idx].z = (float)atom_coord[i].z;
                        }
                        viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemStateChanged, viamd::EventPayloadType_ApplicationState, &state);
                        state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
                        vib.coord_modified = false;
                    }

                    ImGui::TreePop();
                }
            }
        }
        // An index no table claimed this frame loses its hover, so nothing can stay lit after the
        // cursor has moved away. Dropping one the plots resolved is harmless: plot_peaks re-resolves
        // it as the plots draw, ahead of every reader. Same for the NTO index listbox, which drives
        // rsp.hovered from draw_nto_window - that runs after this window, so its value survives to
        // the plots on the next frame.
        if (!rsp_hover_claimed) {
            rsp.hovered = -1;
        }
        if (!vib_hover_claimed) {
            vib.hovered = -1;
        }
        ImGui::End();

        if (rsp.show_export_window) { draw_rsp_spectra_export_window(state); }
    }

    void draw_orb_window(ApplicationState& state) {
        // The XPS plot hands over the MO of the peak under the cursor. draw_rsp_window runs after
        // this one in the frame tick, so what is read here is one frame old - invisible for a hover
        // highlight, and cheaper than reordering the windows. Consumed unconditionally, ahead of the
        // early returns, so a hover that has gone away - or a Response window that stopped drawing
        // altogether - cannot leave a row lit up forever.
        const int32_t xps_highlight_mo = xps.highlight_mo_idx;
        xps.highlight_mo_idx = -1;

        if (!orb.show_window) return;
        if (num_molecular_orbitals() == 0) return;
        ImGui::SetNextWindowSize({600, 300}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Orbital Grid", &orb.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {

            // The neutral names, which mdlib already aliases onto whatever the reader called them.
            // Every one of these is optional - a TREXIO file, for one, need not carry orbital energies -
            // and each is indexed by mo index over the coefficient matrix's row count, so a column is
            // only used when it covers all of those rows. Everything below reads a null as "show no
            // value", never as a column to index.
            const size_t num_mo_rows = num_molecular_orbitals();
            const double* occ_alpha = attribute_series_f64_covering(state.mold.sys, es_path::alpha_occupation, num_mo_rows);
            const double* occ_beta  = attribute_series_f64_covering(state.mold.sys, es_path::beta_occupation,  num_mo_rows);
            const double* ene_alpha = attribute_series_f64_covering(state.mold.sys, es_path::alpha_energy,     num_mo_rows);
            const double* ene_beta  = attribute_series_f64_covering(state.mold.sys, es_path::beta_energy,      num_mo_rows);

            const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();

            // The SCF type as the DATA states it - see es_has_distinct_beta_* in viamd.h. Two spin
            // channels of coefficients is what an unrestricted grid is actually about, and one set
            // of orbitals with two occupations is what makes a SOMO; neither needs a reader's enum.
            const bool unrestricted = es_has_distinct_beta_orbitals(state.mold.sys);
            const bool open_shell   = !unrestricted && es_has_distinct_beta_occupations(state.mold.sys);

            int num_x = (unrestricted) ? 2 : orb.num_x;
            int num_y = orb.num_y;

            int orb_homo_idx = MD_MAX(homo_idx[0], homo_idx[1]);
            int orb_lumo_idx = MD_MIN(lumo_idx[0], lumo_idx[1]);
            if (open_shell) {
                orb_lumo_idx = MD_MAX(lumo_idx[0], lumo_idx[1]);
            }

            const int num_tot_mos = (int)num_molecular_orbitals();

            int num_mos = num_x * num_y;
            // Number of distinct orbitals on screen: an unrestricted grid spends its two columns
            // on the alpha and the beta half of the same orbital, so a grid row is one orbital.
            const int window_size = (unrestricted) ? num_y : num_mos;

            // The table is sortable and the grid shows a span of it, so the span is taken in ROW
            // space, not in mo index space. Under a sort by energy or occupancy the orbitals on
            // screen are no longer contiguous by index, and a span taken by index would show a
            // different set than the rows it highlights.
            md_temp_scope_t temp = md_temp_begin();
            MD_DEFER { md_temp_end(temp); };
            int* row_to_mo = md_temp_alloc_array(temp, int, num_tot_mos);
            int* mo_to_row = md_temp_alloc_array(temp, int, num_tot_mos);

            // A column the current scf type does not display must not stay selected, or a
            // restricted file opened after an unrestricted one would sort on a hidden key.
            if (!unrestricted &&
                (orb.sort_column == OrbitalColumn_OccBeta || orb.sort_column == OrbitalColumn_EneBeta)) {
                orb.sort_column = OrbitalColumn_Idx;
                orb.sort_ascending = false;
            }

            auto rebuild_order = [&]() {
                orbital_table_build_order(row_to_mo, mo_to_row, num_tot_mos, orb.sort_column, orb.sort_ascending,
                                          occ_alpha, occ_beta, ene_alpha, ene_beta, !unrestricted);
            };
            rebuild_order();

            // Centered on the selected orbital's row. Deliberately not clamped to the ends: the
            // window keeps its size and runs off the list, and panels past the end draw empty,
            // exactly as they did when the span was taken in index space.
            int beg_row = 0;
            auto recompute_window = [&]() {
                const int cur_row = (0 <= orb.mo_idx && orb.mo_idx < num_tot_mos) ? mo_to_row[orb.mo_idx] : 0;
                beg_row = cur_row - window_size / 2;
            };
            recompute_window();

            // Slot k of the window counted from the BOTTOM, which is the order the grid lays its
            // slots out in - slot 0 is the bottom right cell. Rows past either end return a
            // negative, distinct per slot: every mo_idx < 0 is skipped when drawing and when
            // evaluating, and keeping them distinct stops two empty slots from matching each
            // other in the volume reuse search below and swapping back and forth forever.
            auto window_mo_idx = [&](int k) -> int {
                const int row = beg_row + window_size - 1 - k;
                return (0 <= row && row < num_tot_mos) ? row_to_mo[row] : (-1 - k);
            };

            // LEFT PANE
            {
                ImGui::BeginChild("left pane", ImVec2(300, 0), ImGuiChildFlags_Borders | ImGuiChildFlags_ResizeX);
                ImGui::SameLine();

                ImGui::PushItemWidth(-1);
                ImGui::BeginGroup();

                ImGui::SliderInt("##Rows", &orb.num_y, 1, 4);

                // An unrestricted grid is always two columns wide, alpha next to beta. Show the
                // slider at the count actually in effect rather than at the stored preference,
                // but keep writing back only when it is live, so the preference survives for the
                // next restricted dataset.
                const bool cols_fixed = (unrestricted);
                int num_x_slider = cols_fixed ? num_x : orb.num_x;
                if (cols_fixed) ImGui::PushDisabled();
                if (ImGui::SliderInt("##Cols", &num_x_slider, 1, 4) && !cols_fixed) {
                    orb.num_x = num_x_slider;
                }
                if (cols_fixed) ImGui::PopDisabled();

                const double iso_min = 1.0e-4;
                const double iso_max = 5.0;
                double iso_val = orb.iso.values[0];
                ImGui::TextUnformatted(electronic_structure_iso_value_label());
                ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
                ImGui::SetItemTooltip("Visual surface threshold. Orbital grids render both +isovalue and -isovalue surfaces.");

                orb.iso.values[0] = (float)iso_val;
                orb.iso.values[1] = -(float)iso_val;
                orb.iso.count = 2;

                ImGui::ColorEdit4("##Color Positive", orb.iso.colors[0].elem);
                ImGui::SetItemTooltip("Color Positive");
                ImGui::ColorEdit4("##Color Negative", orb.iso.colors[1].elem);
                ImGui::SetItemTooltip("Color Negative");

                const char* btn_text = "Goto HOMO";
                if (open_shell && (homo_idx[0] != homo_idx[1])) {
                    btn_text = "Goto SOMO";
                }

                if (ImGui::IsWindowAppearing()) {
                    orb.scroll_to_idx = orb.mo_idx;
                }
                if (ImGui::Button(btn_text, ImVec2(-1, 0))) {
                    orb.scroll_to_idx = homo_idx[0];
                }

                int num_cols = (unrestricted) ? 5 : 3;

                const ImGuiTableFlags flags =
                    ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_RowBg |
                    ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_NoBordersInBody | ImGuiTableFlags_ScrollY |
                    ImGuiTableFlags_Sortable;
                if (ImGui::BeginTable("Molecular Orbitals", num_cols, flags))//, ImVec2(0.0f, TEXT_BASE_HEIGHT * 15), 0.0f))
                {
                    // Declare columns
                    // We use the "user_id" parameter of TableSetupColumn() to specify a user id that will be stored in the sort specifications.
                    // This is so our sort function can identify a column given our own identifier. We could also identify them based on their index!
                    // Demonstrate using a mixture of flags among available sort-related flags:
                    // - ImGuiTableColumnFlags_DefaultSort
                    // - ImGuiTableColumnFlags_NoSort / ImGuiTableColumnFlags_NoSortAscending / ImGuiTableColumnFlags_NoSortDescending
                    // - ImGuiTableColumnFlags_PreferSortAscending / ImGuiTableColumnFlags_PreferSortDescending
                    // MO carries PreferSortDescending as well as DefaultSort so the table opens on
                    // highest index first - the order it had before it was sortable.
                    const ImGuiTableColumnFlags mo_col_flags =
                        ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed;
                    if (unrestricted) {
                        ImGui::TableSetupColumn("MO", mo_col_flags, 0.0f, OrbitalColumn_Idx);
                        ImGui::TableSetupColumn((const char*)u8"Occ. α", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, OrbitalColumn_OccAlpha);
                        ImGui::TableSetupColumn((const char*)u8"Occ. β", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, OrbitalColumn_OccBeta);
                        ImGui::TableSetupColumn((const char*)u8"Ene. α", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, OrbitalColumn_EneAlpha);
                        ImGui::TableSetupColumn((const char*)u8"Ene. β", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, OrbitalColumn_EneBeta);
                    }
                    else {
                        ImGui::TableSetupColumn("MO", mo_col_flags, 0.0f, OrbitalColumn_Idx);
                        ImGui::TableSetupColumn("Occupancy", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, OrbitalColumn_OccAlpha);
                        ImGui::TableSetupColumn("Energy", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, OrbitalColumn_EneAlpha);
                    }
                    ImGui::TableSetupScrollFreeze(0, 1); // Make row always visible

                    // Re-sorted here rather than on the next frame: the rows below and the grid
                    // slots further down both read this order within this same frame, and a lag
                    // would show the grid a span of the previous order.
                    if (ImGuiTableSortSpecs* sort_specs = ImGui::TableGetSortSpecs()) {
                        if (sort_specs->SpecsDirty) {
                            const int  col = (sort_specs->SpecsCount > 0) ? (int)sort_specs->Specs[0].ColumnUserID : (int)OrbitalColumn_Idx;
                            const bool asc = (sort_specs->SpecsCount > 0) && (sort_specs->Specs[0].SortDirection == ImGuiSortDirection_Ascending);
                            if (col != orb.sort_column || asc != orb.sort_ascending) {
                                orb.sort_column = col;
                                orb.sort_ascending = asc;
                                rebuild_order();
                                recompute_window();
                            }
                            sort_specs->SpecsDirty = false;
                        }
                    }

                    ImGui::TableHeadersRow();

                    for (int row = 0; row < num_tot_mos; ++row) {
                        const int n = row_to_mo[row];
                        ImGui::PushID(n + 1);
                        ImGui::TableNextRow();
                        if (n == xps_highlight_mo) {
                            ImGui::TableSetBgColor(ImGuiTableBgTarget_RowBg1, ImGui::GetColorU32(COLOR_ROW_HIGHLIGHT));
                        }
                        // Contiguous in rows, so the highlighted block follows the sort and always
                        // names the same orbitals the grid is showing.
                        bool is_selected = (beg_row <= row && row < beg_row + window_size);
                        ImGui::TableNextColumn();
                        if (orb.scroll_to_idx != -1 && n == orb.scroll_to_idx) {
                            orb.scroll_to_idx = -1;
                            ImGui::SetScrollHereY();
                        }

                        char buf[32];
                        snprintf(buf, sizeof(buf), "%i", n + 1);

                        ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                        if (ImGui::Selectable(buf, is_selected, selectable_flags)) {
                            if (orb.mo_idx != n) {
                                orb.mo_idx = n;
                            }
                        }
                        ImGui::TableNextColumn();
                        if (occ_alpha) {
                            const char* lbl = "";
                            double occ = occ_alpha[n];
                            if (unrestricted) {
                                // Compared against alpha's own homo/lumo only. The window-wide
                                // orb_homo_idx/orb_lumo_idx are extremes over both spins, so
                                // requiring them here left whichever spin did not hold the
                                // extreme permanently unlabelled.
                                if (n == homo_idx[0]) {
                                    lbl = "HOMO";
                                }
                                else if (n == lumo_idx[0]) {
                                    lbl = "LUMO";
                                }
                            }
                            else {
                                if (occ_beta) {
                                    occ += occ_beta[n];
                                }
                                if (n == orb_homo_idx) {
                                    lbl = (homo_idx[0] == homo_idx[1]) ? "HOMO" : "SOMO";
                                }
                                else if (n == orb_lumo_idx) {
                                    lbl = "LUMO";
                                }
                            }
                            ImGui::Text("%.1f %s", occ, lbl);
                        }
                        else {
                            ImGui::Text("-");
                        }
                        if (unrestricted) {
                            ImGui::TableNextColumn();
                            if (occ_beta) {
                                const char* lbl = "";
                                if (n == homo_idx[1]) {
                                    lbl = "HOMO";
                                }
                                else if (n == lumo_idx[1]) {
                                    lbl = "LUMO";
                                }
                                ImGui::Text("%.1f %s", occ_beta[n], lbl);
                            }
                            else {
                                ImGui::Text("-");
                            }
                        }
                        ImGui::TableNextColumn();
                        if (ene_alpha) {
                            ImGui::Text("%.4f", ene_alpha[n]);
                        }
                        else {
                            ImGui::Text("-");
                        }
                        if (unrestricted) {
                            ImGui::TableNextColumn();
                            if (ene_beta) {
                                ImGui::Text("%.4f", ene_beta[n]);
                            }
                            else {
                                ImGui::Text("-");
                            }
                        }
                        ImGui::PopID();
                    }
                    ImGui::EndTable();
                }

                ImGui::EndGroup();
                ImGui::PopItemWidth();

                ImGui::EndChild();
            }

            ImGui::SameLine();

            // These represent the new mo_idx we want to have in each slot
            int vol_mo_idx[16] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
            SpinChannel vol_mo_type[16] = {};
            for (int i = 0; i < num_mos; ++i) {
                if (unrestricted) {
                    vol_mo_idx[i] = window_mo_idx(i / 2);
                    vol_mo_type[i] = (i & 1) ? SPIN_ALPHA : SPIN_BETA;
                }
                else {
                    vol_mo_idx[i] = window_mo_idx(i);
                    vol_mo_type[i] = SPIN_ALPHA;
                }
            }

            int job_queue[16];
            int num_jobs = 0;
            // Find and reuse volume data from existing slots (if applicable)
            // If there is no existing volume, we queue up a new job
            for (int i = 0; i < num_mos; ++i) {
                // Check if we already have that entry in the correct slot
                if (orb.vol_mo_idx[i] == vol_mo_idx[i] && orb.vol_mo_type[i] == vol_mo_type[i]) continue;

                // Try to find the entry in the existing list
                bool found = false;
                for (int j = 0; j < num_mos; ++j) {
                    if (i == j) continue;
                    // The spin has to match as well, not just the index. In the unrestricted
                    // layout two slots share every mo index and differ only in spin, so an
                    // index-only comparison lets a slot claim its own opposite-spin twin: it
                    // renders the wrong spin, and the swap drags an already-resolved slot's
                    // volume away with it. With the spin compared, no resolved slot can match,
                    // because (mo_idx, spin) is unique across slots.
                    if (vol_mo_idx[i] == orb.vol_mo_idx[j] && vol_mo_type[i] == orb.vol_mo_type[j]) {
                        // Swap to correct location
                        ImSwap(orb.vol[i], orb.vol[j]);
                        ImSwap(orb.vol_mo_idx[i], orb.vol_mo_idx[j]);
                        ImSwap(orb.vol_mo_type[i], orb.vol_mo_type[j]);
                        found = true;
                        break;
                    }
                }

                // If not found, put in job queue to compute the volume
                if (!found) {
                    job_queue[num_jobs++] = i;
                }
            }

            if (num_jobs > 0) {
                const float samples_per_unit_length = DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM;
                md_grid_t grid { 0 };
#if USE_AABB_GRID
                init_grid(&grid, mat3_ident(), aabb.min_ext, aabb.max_ext, samples_per_unit_length);
#else
                init_grid(&grid, oabb.orientation, oabb.min_ext, oabb.max_ext, samples_per_unit_length);
#endif
                for (int i = 0; i < num_jobs; ++i) {
                    int slot_idx = job_queue[i];
                    int mo_idx = vol_mo_idx[slot_idx];
                    SpinChannel mo_type = vol_mo_type[slot_idx];
                    orb.vol_mo_idx[slot_idx] = mo_idx;
                    orb.vol_mo_type[slot_idx] = mo_type;

                    if (-1 < mo_idx && mo_idx < num_tot_mos) {
                        init_volume(&orb.vol[slot_idx], grid);
                        evaluate_mo(state, orb.vol[slot_idx].tex_id, grid, mo_type, mo_idx, MD_GTO_OP_SET);
                    }
                }
            }

            ImVec2 canvas_sz = ImMax(ImVec2(50.0f, 50.0f), ImGui::GetContentRegionAvail());   // Resize canvas to what's available

            // Animate camera towards targets
            const double dt = state.app.timing.delta_s;
            camera_animate(&orb.camera, orb.target, dt);

            ImGui::Dummy(canvas_sz);

            // Draw border and background color
            ImGuiIO& io = ImGui::GetIO();

            ImVec2 canvas_min = ImGui::GetItemRectMin();
            ImVec2 canvas_max = ImGui::GetItemRectMax();
            ImVec2 orb_win_sz = ImFloor(canvas_sz / ImVec2((float)num_x, (float)num_y));

            const float aspect = orb_win_sz.x / orb_win_sz.y;

            //const mat4_t world_to_view = camera_world_to_view_matrix(orb.camera);
            //const mat4_t view_to_clip  = camera_view_to_clip_matrix_persp(orb.camera, aspect);

            const mat4_t view_to_world = camera_view_to_world_matrix(orb.camera);
            const mat4_t clip_to_view  = camera_clip_to_view_matrix_persp(orb.camera, aspect);
            const mat4_t clip_to_world = mat4_mul(view_to_world, clip_to_view);

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddRectFilled(canvas_min, canvas_max, IM_COL32(255, 255, 255, 255));

            // The render targets, sized before the panels go into the draw list: it records the texture
            // names, and a texture reallocated later in the frame would be drawn under a dead name
            int width  = MD_MAX(1, (int)(orb_win_sz.x * io.DisplayFramebufferScale.x));
            int height = MD_MAX(1, (int)(orb_win_sz.y * io.DisplayFramebufferScale.y));

            auto& gbuf = orb.gbuf;
            if ((int)gbuf.width != width || (int)gbuf.height != height) {
                gbuffer_init(&gbuf, width, height);
                for (int i = 0; i < (int)MD_COUNTOF(orb.iso_tex); ++i) {
                    if (orb.iso_tex[i]) gl::free_texture(&orb.iso_tex[i]);
                    if (orb.out_tex[i]) gl::free_texture(&orb.out_tex[i]);
                }
            }
            for (int i = 0; i < num_mos; ++i) {
                if (!orb.iso_tex[i]) gl::init_texture_2D(&orb.iso_tex[i], width, height, GL_RGBA16F);
                if (!orb.out_tex[i]) gl::init_texture_2D(&orb.out_tex[i], width, height, GL_RGBA8);
            }

            for (int i = 0; i < num_mos; ++i) {
                ImGui::PushID(i);
                MD_DEFER { ImGui::PopID(); };

                const int mo_idx = window_mo_idx((unrestricted) ? i / 2 : i);
                int x = num_x - i % num_x - 1;
                int y = num_y - i / num_x - 1;
                if (-1 < mo_idx && mo_idx < num_tot_mos) {
                    const ImVec2 p0 = canvas_min + orb_win_sz * ImVec2((float)(x + 0), (float)(y + 0));
                    const ImVec2 p1 = canvas_min + orb_win_sz * ImVec2((float)(x + 1), (float)(y + 1));
                    const ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    const ImVec2 text_pos_tl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p0.y + TEXT_BASE_HEIGHT * 0.25f);
                    const ImVec2 text_pos_br = ImVec2(p1.x - TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);

                    const ImVec2 sz = p1 - p0;

                    ImGui::SetCursorScreenPos(p0);
                    InteractionSurfaceState surface_state = interaction_surface(interaction_surface_orb, vec_cast(sz), InteractionSurfaceFlags_NoRegionSelect);

                    PickingHit hit = {};
                    if (surface_state.hovered) {
                        InteractionSurfaceHitArgs hit_args = {
                            .picking_surface = &orb.picking_surface,
                            .picking_handler = state.picking_handler,
                            .fbo = orb.gbuf.fbo,
                            .width = orb.gbuf.width,
                            .height = orb.gbuf.height,
                            .clip_to_world = clip_to_world,
                        };
                        interaction_surface_hit_extract(&hit, surface_state, hit_args);

                        InteractionSurfaceEvent event = {};
                        interaction_surface_event_extract(&event, surface_state, hit);

                        viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
                    }

                    InteractionSurfaceViewTransformArgs view_args = {
                        .camera = orb.camera,
                        .trackball_param = state.view.trackball_param,
                    };

                    InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&orb.target, surface_state, view_args);
                    if (view_result.reset_requested) {
                        ViewTransform reset_transform = default_view;
                        if (hit.depth < 1.0f) {
                            reset_transform.distance = orb.target.distance;
                            reset_transform.orientation = orb.camera.orientation;
                            reset_transform.position = hit.world_pos + orb.camera.orientation * vec3_set(0, 0, orb.target.distance);
                        }
                        orb.target = reset_transform;
                    }


                    const char* lbl = "";
                    if (unrestricted) {
                        // 0 = alpha, 1 = beta, matching vol_mo_type above. A panel shows one
                        // spin, so it is labelled against that spin's own homo/lumo; gating on
                        // the window-wide extremes left the other spin unlabelled.
                        int j = (i & 1) ? 0 : 1;
                        if (mo_idx == homo_idx[j]) {
                            lbl = "(HOMO)";
                        }
                        else if (mo_idx == lumo_idx[j]) {
                            lbl = "(LUMO)";
                        }
                    }
                    else {
                        if (mo_idx == orb_homo_idx) {
                            lbl = (homo_idx[0] == homo_idx[1]) ? "(HOMO)" : "(SOMO)";
                        }
                        else if (mo_idx == orb_lumo_idx) {
                            lbl = "(LUMO)";
                        }
                    }

                    char buf[32];
                    snprintf(buf, sizeof(buf), "%i %s", mo_idx + 1, lbl);
                    draw_list->AddImage((ImTextureID)(intptr_t)orb.out_tex[i], p0, p1, { 0,1 }, { 1,0 });
                    draw_list->AddText(text_pos_bl, ImColor(0, 0, 0), buf);

                    // The energy label is dropped, not faked, when the file carries no energies
                    // for this panel's spin.
                    const double* ene = ene_alpha;
                    if (unrestricted) {
                        draw_list->AddText(text_pos_tl, ImColor(0, 0, 0), (i & 1) ? (const char*)u8"α" : (const char*)u8"β");
                        ene = (i & 1) ? ene_alpha : ene_beta;
                    }
                    if (ene) {
                        snprintf(buf, sizeof(buf), "%.4f", ene[mo_idx]);
                        float width = ImGui::CalcTextSize(buf).x;
                        draw_list->AddText(text_pos_br - ImVec2(width, 0), ImColor(0, 0, 0), buf);
                    }
                }
            }

            // num_x, not orb.num_x: an unrestricted grid is forced to two columns, and the stored
            // preference drew a stray line on the canvas edge.
            for (int x = 1; x < num_x; ++x) {
                ImVec2 p0 = { canvas_min.x + orb_win_sz.x * x, canvas_min.y };
                ImVec2 p1 = { canvas_min.x + orb_win_sz.x * x, canvas_max.y };
                draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
            }

            for (int y = 1; y < num_y; ++y) {
                ImVec2 p0 = { canvas_min.x, canvas_min.y + orb_win_sz.y * y };
                ImVec2 p1 = { canvas_max.x, canvas_min.y + orb_win_sz.y * y };
                draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
            }

            const ImVec2 origin(canvas_min.x, canvas_min.y);  // Lock scrolled origin
            const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

            if (orb.show_coordinate_system_widget) {
                float ext = MD_MIN(orb_win_sz.x, orb_win_sz.y) * 0.4f;
                float pad = 20.0f;
                ImVec2 size = ImVec2(ext, ext);
                ImGui::SetCursorScreenPos(ImVec2(canvas_min.x + pad, canvas_max.y - ext - pad));

                quat_t out_orientation = orb.target.orientation;
                if (ImGui::CoordinateSystemWidget(&out_orientation, orb.camera.orientation, size)) {
                    const vec3_t look_at = camera_get_look_at(orb.target);
                    orb.target.orientation = quat_normalize(out_orientation);
                    orb.target.position = camera_position_from_look_at(look_at, orb.target.orientation, orb.target.distance);
                }
            }

            const float aspect_ratio = orb_win_sz.x / orb_win_sz.y;
            mat4_t view_mat = camera_world_to_view_matrix(orb.camera);
            mat4_t proj_mat = camera_view_to_clip_matrix_persp(orb.camera, aspect_ratio);
            mat4_t inv_proj_mat = camera_clip_to_view_matrix_persp(orb.camera, aspect_ratio);

            // The mesh rendering mode: each panel's orbital as a field, the meshes of all of them extracted together
            const bool mesh_mode = state.settings.isosurface_rendering == IsosurfaceRendering_Mesh;
            if (mesh_mode) {
                const uint64_t pos_hash = panel_positions_hash(state);
                for (int i = 0; i < PanelMeshes::MAX; ++i) {
                    const int mo_idx = (i < num_mos) ? orb.vol_mo_idx[i] : -1;
                    const SpinChannel mo_type = orb.vol_mo_type[i];
                    uint64_t key = 0;
                    if (-1 < mo_idx && mo_idx < num_tot_mos) {
                        const int id[2] = {mo_idx, (int)mo_type};
                        key = md_hash64(id, sizeof(id), pos_hash);
                        key = key ? key : 1;
                    }
                    panel_meshes_field(&orb.mesh, &state, i, key, [&]() {
                        const md_str_t path = (mo_type == SPIN_BETA) ? es_path::beta_coefficient : es_path::alpha_coefficient;
                        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)mo_idx);
                        return iso_field_create_orbital(&state, path, &slice);
                    });
                }
                panel_meshes_update(&orb.mesh, &state, num_mos, fabs((double)orb.iso.values[0]), MD_GTO_ISO_SHEET_BOTH);
            } else if (orb.mesh.meshes.key) {
                panel_meshes_free(&orb.mesh, &state);
            }

            // The panels are only redrawn when something they show changed, or while their meshes are refined
            // over the frames: a static grid costs nothing
            uint64_t render_hash = md_hash64(&view_mat, sizeof(view_mat), (uint64_t)gl_rep.id);
            {
                render_hash = md_hash64(&state.settings.isosurface_rendering, sizeof(state.settings.isosurface_rendering), render_hash);
                if (mesh_mode) {
                    render_hash = md_hash64(&orb.mesh.meshes.mesh.version, sizeof(orb.mesh.meshes.mesh.version), render_hash);
                }
                const int layout[3] = { width, height, num_mos };
                render_hash = md_hash64(layout, sizeof(layout), render_hash);
                for (int i = 0; i < num_mos; ++i) {
                    const uint64_t vol_version = volume::data_version(orb.vol[i].tex_id);
                    render_hash = md_hash64(&orb.vol[i].tex_id, sizeof(orb.vol[i].tex_id), render_hash);
                    render_hash = md_hash64(&vol_version, sizeof(vol_version), render_hash);
                    render_hash = md_hash64(&orb.vol[i].texture_to_world, sizeof(orb.vol[i].texture_to_world), render_hash);
                }
                render_hash = md_hash64(&proj_mat, sizeof(proj_mat), render_hash);
                render_hash = md_hash64(&orb.iso, sizeof(orb.iso), render_hash);
                render_hash = md_hash64(&state.settings.exact_isosurfaces, sizeof(state.settings.exact_isosurfaces), render_hash);
                render_hash = md_hash64(&state.visuals.tonemapping, sizeof(state.visuals.tonemapping), render_hash);
                render_hash = md_hash64(&state.visuals.shading, sizeof(state.visuals.shading), render_hash);
                render_hash = md_hash64(&state.mold.gpu_buffers_version, sizeof(state.mold.gpu_buffers_version), render_hash);
                // Baked into the picking buffer, which hovering reads
                const md_gl_picking_offset_t picking_offset = gl_picking_offset(state);
                render_hash = md_hash64(&picking_offset, sizeof(picking_offset), render_hash);
                render_hash = render_hash ? render_hash : 1;
            }

            if (gl_rep.id && (render_hash != orb.render_hash || (mesh_mode && !orb.mesh.converged))) {
                orb.render_hash = render_hash;

                gbuffer_clear(&gbuf);

                const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                    GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

                glEnable(GL_CULL_FACE);
                glCullFace(GL_BACK);

                glEnable(GL_DEPTH_TEST);
                glDepthFunc(GL_LESS);
                glDepthMask(GL_TRUE);
                glEnable(GL_SCISSOR_TEST);

                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
                glDrawBuffers((int)MD_COUNTOF(draw_buffers), draw_buffers);
                glViewport(0, 0, gbuf.width, gbuf.height);
                glScissor(0, 0, gbuf.width, gbuf.height);

                md_gl_draw_op_t draw_op = {};
                draw_op.type = MD_GL_REP_BALL_AND_STICK;
                draw_op.args.ball_and_stick.ball_scale   = 1.0f;
                draw_op.args.ball_and_stick.stick_radius = 1.0f;
                draw_op.args.ball_and_stick.bond_orders  = true;
                draw_op.rep = gl_rep;

                md_gl_draw_args_t draw_args = {
                    .shaders = state.gl.shaders,
                    .draw_operations = {
                        .count = 1,
                        .ops = &draw_op
                    },
                    .view_transform = {
                        .view_matrix = (const float*)view_mat.elem,
                        .proj_matrix = (const float*)proj_mat.elem,
                    },
                    .picking_offset = gl_picking_offset(state),
                };

                md_gl_draw(&draw_args);

                glDisable(GL_SCISSOR_TEST);

                postprocess_pipeline::Settings postprocess_settings = {};
                postprocess_settings.background_color = {1.f, 1.f, 1.f};
                postprocess_settings.shading.mode = state.visuals.shading.mode;
                postprocess_settings.shading.ambient = state.visuals.shading.ambient;
                postprocess_settings.shading.specular = state.visuals.shading.specular;
                postprocess_settings.shading.colors = state.visuals.shading.colors;
                postprocess_settings.shading.light = state.visuals.shading.light;
                postprocess_settings.shading.environment = state.visuals.shading.environment;
                postprocess_settings.tonemap.enabled = state.visuals.tonemapping.enabled;
                postprocess_settings.tonemap.mode = state.visuals.tonemapping.tonemapper;
                postprocess_settings.tonemap.exposure = state.visuals.tonemapping.exposure;
                postprocess_settings.tonemap.gamma = state.visuals.tonemapping.gamma;
                postprocess_settings.shadows.enabled = false;
                postprocess_settings.ssao.enabled = false;
                postprocess_settings.dof.enabled = false;
                postprocess_settings.fxaa.enabled = true;
                postprocess_settings.taa.enabled = false;
                postprocess_settings.sharpen.enabled = false;

                ViewParam view_param = {
                    .matrix = {
                        .curr = {
                        .view = view_mat,
                        .proj = proj_mat,
                        .norm = view_mat,
                    },
                    .inv = {
                        .proj = inv_proj_mat,
                    }
                    },
                    .clip_planes = {
                        .near = orb.camera.near_plane,
                        .far  = orb.camera.far_plane,
                    },
                    .resolution = {orb_win_sz.x, orb_win_sz.y},
                    .fov_y = orb.camera.fov_y,
                };

                if (!orb.out_fbo) glGenFramebuffers(1, &orb.out_fbo);

                // Every panel is its own orbital's isosurfaces over the one opaque image they share, put
                // through post-processing on its own: the surfaces are composited in HDR, before tone mapping
                PUSH_GPU_SECTION("ORB GRID PANELS")
                bool converged = true;
                for (int i = 0; i < num_mos; ++i) {
                    bool iso_written = false;
                    if (mesh_mode) {
                        // The sheet of each sign in the colour of the isovalue of that sign
                        const bool pos_first = orb.iso.values[0] >= 0.0f;
                        const vec4_t colors[2] = {orb.iso.colors[pos_first ? 0 : 1], orb.iso.colors[pos_first ? 1 : 0]};
                        isomesh::Surface surfaces[2];
                        const size_t num_surfaces = panel_meshes_surfaces(surfaces, orb.mesh, i, colors);
                        const isomesh::RenderDesc mesh_desc = {
                            .render_target = {
                                .depth  = orb.gbuf.tex.depth,
                                .color  = orb.iso_tex[i],
                                .width  = orb.gbuf.width,
                                .height = orb.gbuf.height,
                                .clear_color = true,
                            },
                            .matrix = {
                                .view = view_mat,
                                .proj = proj_mat,
                                .proj_no_jitter = proj_mat,
                            },
                            .surfaces = surfaces,
                            .num_surfaces = num_surfaces,
                            .shading = {
                                .roughness = 0.3f,
                                .ior = 1.5f,
                                .mode = state.visuals.shading.mode,
                                .ambient = state.visuals.shading.ambient,
                                .specular = state.visuals.shading.specular,
                                .colors = state.visuals.shading.colors,
                                .light = state.visuals.shading.light,
                                .environment = state.visuals.shading.environment,
                            },
                            .scene_hash = render_hash,
                        };
                        bool panel_converged = true;
                        iso_written = isomesh::render(&orb.mesh.target[i], mesh_desc, &panel_converged);
                        converged &= panel_converged;
                    } else {
                        volume::IsoRenderDesc vol_desc = {
                            .render_target = {
                                .depth  = orb.gbuf.tex.depth,
                                .color  = orb.iso_tex[i],
                                .width  = orb.gbuf.width,
                                .height = orb.gbuf.height,
                                .clear_color = true,
                            },
                            .texture = {
                                .density_volume = orb.vol[i].tex_id,
                            },
                            .matrix = {
                                .model = orb.vol[i].texture_to_world,
                                .view  = view_mat,
                                .proj  = proj_mat,
                                .inv_proj = inv_proj_mat,
                            },
                            .iso = {
                                .count  = (size_t)orb.iso.count,
                                .values = orb.iso.values,
                                .colors = orb.iso.colors,
                                .exact  = state.settings.exact_isosurfaces,
                            },
                            .shading = {
                                .roughness = 0.3f,
                                .ior = 1.5f,
                                .mode = state.visuals.shading.mode,
                                .ambient = state.visuals.shading.ambient,
                                .specular = state.visuals.shading.specular,
                                .colors = state.visuals.shading.colors,
                                .light = state.visuals.shading.light,
                                .environment = state.visuals.shading.environment,
                            },
                        };
                        iso_written = volume::render_isosurfaces(vol_desc);
                    }

                    postprocess_pipeline::Inputs postprocess_inputs = {};
                    postprocess_inputs.depth = orb.gbuf.tex.depth;
                    postprocess_inputs.color = orb.gbuf.tex.color;
                    postprocess_inputs.normal = orb.gbuf.tex.normal;
                    postprocess_inputs.velocity = orb.gbuf.tex.velocity;
                    postprocess_inputs.transparency_hdr = iso_written ? orb.iso_tex[i] : 0;

                    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, orb.out_fbo);
                    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, orb.out_tex[i], 0);
                    glDrawBuffer(GL_COLOR_ATTACHMENT0);
                    glViewport(0, 0, orb.gbuf.width, orb.gbuf.height);
                    postprocess_pipeline::execute(postprocess_inputs, postprocess_settings, view_param);
                }
                orb.mesh.converged = converged;
                POP_GPU_SECTION()

                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
                glDrawBuffer(GL_BACK);
            }
        }
        ImGui::End();
    }

    void draw_export_window(ApplicationState& state) {
        if (!export_state.show_window) return;
        if (ImGui::Begin("Quantum Chemistry Export", &export_state.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginCombo("Source", electronic_structure_source_str[(int)export_state.source])) {
                for (int i = 0; i < (int)ElectronicStructureSource::Count; i++) {
                    ElectronicStructureSource source = (ElectronicStructureSource)i;
                    bool is_selected = (export_state.source == source);
                    bool disabled = source == ElectronicStructureSource::DensityProperty ||
                                    ((source == ElectronicStructureSource::NaturalTransitionOrbital ||
                                     source == ElectronicStructureSource::TransitionDensity) &&
                                    (num_excited == 0));

                    if (disabled) ImGui::PushDisabled();
                    if (ImGui::Selectable(electronic_structure_source_str[i], is_selected)) {
                        export_state.source = source;
                        export_state.use_magnitude = false;
                        if (export_state.source == ElectronicStructureSource::ElectronDensity) {
                            export_state.spin = ElectronicStructureSpin::Total;
                        }
                    }
                    if (disabled) ImGui::PopDisabled();

                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            if (export_state.source == ElectronicStructureSource::MolecularOrbital ||
                export_state.source == ElectronicStructureSource::NaturalTransitionOrbital ||
                (export_state.source == ElectronicStructureSource::ElectronDensity && export_state.spin == ElectronicStructureSpin::Difference)) {
                const char* magnitude_label = "Magnitude";
                ImGui::Checkbox(magnitude_label, &export_state.use_magnitude);
            }

            if (export_state.source == ElectronicStructureSource::NaturalTransitionOrbital) {
                ImGui::Combo("Component", (int*)&export_state.nto_component, electronic_structure_nto_component_str, IM_ARRAYSIZE(electronic_structure_nto_component_str));
            }

            if (export_state.source == ElectronicStructureSource::TransitionDensity) {
                ImGui::Combo("Component", (int*)&export_state.transition_density_component, electronic_structure_transition_density_component_str, IM_ARRAYSIZE(electronic_structure_transition_density_component_str));
            }

            if (export_state.source == ElectronicStructureSource::ElectronDensity) {
                // Whether there is a beta density to offer is whether one is published. Restricted
                // publishes it as an alias of alpha's, so the selector is offered there too - which
                // is right: Total and Difference are still meaningful choices over the pair.
                const bool has_beta_density = attribute_exists(state.mold.sys, es_path::beta_density);
                if (has_beta_density) {
                    const ElectronicStructureSpin spin_options[] = {
                        ElectronicStructureSpin::Total,
                        ElectronicStructureSpin::Alpha,
                        ElectronicStructureSpin::Beta,
                        ElectronicStructureSpin::Difference,
                    };
                    const char* spin_labels[] = { "Total", "Alpha", "Beta", "Difference" };
                    int spin_idx = 0;
                    for (int i = 0; i < (int)IM_ARRAYSIZE(spin_options); ++i) {
                        if (export_state.spin == spin_options[i]) {
                            spin_idx = i;
                            break;
                        }
                    }
                    if (export_state.spin != spin_options[spin_idx]) {
                        export_state.spin = spin_options[spin_idx];
                    }
                    if (ImGui::Combo("Spin", &spin_idx, spin_labels, IM_ARRAYSIZE(spin_labels))) {
                        export_state.spin = spin_options[spin_idx];
                        if (export_state.spin != ElectronicStructureSpin::Difference) {
                            export_state.use_magnitude = false;
                        }
                    }
                } else {
                    export_state.spin = ElectronicStructureSpin::Total;
                    export_state.use_magnitude = false;
                }
            }

            int orb_idx = 0;
            bool show_mo_combo = false;
            bool show_lambda_combo = false;
            bool show_excited_state_combo = false;

            switch (export_state.source) {
            case ElectronicStructureSource::MolecularOrbital:
                show_mo_combo = true;
                break;
            case ElectronicStructureSource::NaturalTransitionOrbital:
                show_excited_state_combo = true;
                show_lambda_combo = true;
                break;
            case ElectronicStructureSource::TransitionDensity:
                show_excited_state_combo = true;
                break;
            case ElectronicStructureSource::ElectronDensity:
                break;
            default:
                break;
            }

            if (show_mo_combo) {
                const bool unrestricted = es_has_distinct_beta_orbitals(state.mold.sys);
                // Only an unrestricted calculation has two distinct sets of coefficients. For
                // restricted and restricted-openshell, scf.beta is a shallow copy of alpha that
                // shares the very same coefficient buffer, and the beta attribute is published
                // as an alias of alpha's, so a selector there would offer two names
                // for identical data. Forcing alpha in that case also keeps a Beta selection
                // left over from a previously loaded file from picking beta's homo/lumo labels.
                if (unrestricted) {
                    const char* options[2] = {"Alpha", "Beta"};
                    if (ImGui::BeginCombo("##MO_TYPE", options[export_state.mo.type])) {
                        for (size_t i = 0; i < MD_COUNTOF(options); ++i) {
                            const bool is_selected = (export_state.mo.type == (SpinChannel)i);
                            if (ImGui::Selectable(options[i], is_selected)) {
                                export_state.mo.type = (SpinChannel)i;
                            }
                            if (is_selected) {
                                ImGui::SetItemDefaultFocus();
                            }
                        }
                        ImGui::EndCombo();
                    }
                } else {
                    export_state.mo.type = SPIN_ALPHA;
                }

                int orb_homo_idx = homo_idx[export_state.mo.type];
                int orb_lumo_idx = lumo_idx[export_state.mo.type];

                auto write_lbl = [orb_homo_idx, orb_lumo_idx](char* buf, size_t cap, int idx) {
                    const char* suffix = "";
                    if (idx == orb_homo_idx) {
                        suffix = "(HOMO)";
                    } else if (idx == orb_lumo_idx) {
                        suffix = "(LUMO)";
                    }
                    snprintf(buf, cap, "%i %s", idx + 1, suffix);
                };

                char lbl[32];
                write_lbl(lbl, sizeof(lbl), export_state.mo.idx);
                int num_orbs = (int)num_molecular_orbitals();
                if (ImGui::BeginCombo("MO Idx", lbl)) {
                    for (int i = 0; i < num_orbs; i++) {                        
                        const bool is_selected = (export_state.mo.idx == i);

                        write_lbl(lbl, sizeof(lbl), i);
                        if (ImGui::Selectable(lbl, is_selected)) {
                            export_state.mo.idx = i;
                        }

                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
                orb_idx = export_state.mo.idx;
            }

            if (show_excited_state_combo) {
                int num_orbs = (int)num_natural_transition_orbitals();
                char lbl[32];
                snprintf(lbl, sizeof(lbl), "%i", export_state.nto.idx + 1);
                if (ImGui::BeginCombo("Excited State Idx", lbl)) {
                    for (int i = 0; i < num_orbs; i++) {
                        const bool is_selected = (export_state.nto.idx == i);

                        snprintf(lbl, sizeof(lbl), "%i", i + 1);
                        if (ImGui::Selectable(lbl, is_selected)) {
                            export_state.nto.idx = i;
                        }

                        if (is_selected) {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
                orb_idx = export_state.nto.idx;
            }

            if (show_lambda_combo) {
                // The same 1e-3 the rest of the tree stops at, and the same place the zero padded
                // tail of the ragged lambda axis stops it.
                const double LAMBDA_CUTOFF = 1.0e-3;
                double lambdas[4] = {0};
                const bool nto_ready = es_nto_ready(&state, (size_t)export_state.nto.idx);
                size_t lambda_count = nto_ready ? es_nto_lambdas(lambdas, MD_COUNTOF(lambdas), state.mold.sys, (size_t)export_state.nto.idx, LAMBDA_CUTOFF) : 0;
                if (!nto_ready) ImGui::TextDisabled("Lambda Idx: computing...");
                if (lambda_count > 0) {
                    if ((size_t)export_state.nto.lambda_idx >= lambda_count) {
                        export_state.nto.lambda_idx = (int)lambda_count - 1;
                    }
                    char lbl[32];
                    int i = export_state.nto.lambda_idx;
                    snprintf(lbl, sizeof(lbl), "%i (%.4f)", i + 1, lambdas[i]);
                    if (ImGui::BeginCombo("Lambda Idx", lbl)) {
                        for (i = 0; i < (int)lambda_count; i++) {
                            const bool is_selected = (export_state.nto.lambda_idx == i);
                            snprintf(lbl, sizeof(lbl), "%i (%.4f)", i + 1, lambdas[i]);
                            if (ImGui::Selectable(lbl, is_selected)) {
                                export_state.nto.lambda_idx = i;
                            }

                            if (is_selected) {
                                ImGui::SetItemDefaultFocus();
                            }
                        }
                        ImGui::EndCombo();
                    }
                }
            }

            vec3_t extent;
            if (export_state.use_obb) {
                extent = oabb.max_ext - oabb.min_ext;
            } else {
                extent = aabb.max_ext - aabb.min_ext;
            }

            int dim[3];
            compute_dim(dim, extent, volume_resolution_samples_per_angstrom[(int)export_state.resolution] * BOHR_TO_ANGSTROM);

            char lbl[32];
            static const char* fmt = (const char*)u8"%s (%i×%i×%i)";
            snprintf(lbl, sizeof(lbl), fmt, volume_resolution_str[(int)export_state.resolution], dim[0], dim[1], dim[2]);
            if (ImGui::BeginCombo("Volume Resolution (XYZ)", lbl)) {
                for (int i = 0; i < IM_ARRAYSIZE(volume_resolution_str); ++i) {
                    const bool is_selected = ((int)export_state.resolution == i);
                    compute_dim(dim, extent, volume_resolution_samples_per_angstrom[i] * BOHR_TO_ANGSTROM);
                    snprintf(lbl, sizeof(lbl), fmt, volume_resolution_str[i], dim[0], dim[1], dim[2]);
                    if (ImGui::Selectable(lbl, is_selected)) {
                        export_state.resolution = (VolumeResolution)i;
                    }

                    if (is_selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            enum {
                EXPORT_FILE_FORMAT_CUBE = 0,
                EXPORT_FILE_FORMAT_RAW_MHD_XYZ,
                EXPORT_FILE_FORMAT_COUNT,
            };
            const char* file_format_str[] = {"cube", "raw + mhd + xyz"};
            md_str_t file_format_ext[] = {MD_STR_INIT("cube"), MD_STR_INIT("raw")};

            MD_STATIC_ASSERT(MD_COUNTOF(file_format_str) == EXPORT_FILE_FORMAT_COUNT);
            MD_STATIC_ASSERT(MD_COUNTOF(file_format_ext) == EXPORT_FILE_FORMAT_COUNT);

            static int file_format = 0;
            ImGui::Combo("File format", &file_format, file_format_str, IM_ARRAYSIZE(file_format_str));

            ImGui::Checkbox("Use OBB", &export_state.use_obb);
            ImGui::SetItemTooltip("Use Oriented Bounding Box for volume (rotate and fit the volume to the data)");

            if (ImGui::Button("Export")) {
                md_str_t ext = file_format_ext[file_format];
                char path_buf[2048];
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ext)) {
                    md_str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};

                    md_file_t file = { 0 };
                    if (!md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                        MD_LOG_ERROR("Failed to open file for writing: '" MD_STR_FMT "'", MD_STR_ARG(path));
                        return;
                    }

                    const double samples_per_unit_length = volume_resolution_samples_per_angstrom[(int)export_state.resolution] * BOHR_TO_ANGSTROM;
                    md_grid_t grid = {};
                    if (export_state.use_obb) {
                        init_grid(&grid, oabb.orientation, oabb.min_ext, oabb.max_ext, samples_per_unit_length);
                    } else {
                        init_grid(&grid, mat3_ident(), aabb.min_ext, aabb.max_ext, samples_per_unit_length);
                    }

                    md_temp_scope_t temp = md_temp_begin();

                    MD_DEFER {
                        md_file_close(&file); 
                        md_temp_end(temp);
                    };

                    Volume vol = {};
                    init_volume(&vol, grid, GL_R32F);

                    switch (export_state.source) {
                    case ElectronicStructureSource::MolecularOrbital:
                    {
                        md_gto_op_t op = es_gto_op(export_state.use_magnitude);
                        evaluate_mo(state, vol.tex_id, grid, export_state.mo.type, export_state.mo.idx, op);
                        break;
                    }
                    case ElectronicStructureSource::NaturalTransitionOrbital:
                    {
                        md_gto_op_t op = es_gto_op(export_state.use_magnitude);

                        // The export's own component selection IS what evaluate_nto takes; there is
                        // nothing left to translate it through.
                        evaluate_nto(state, vol.tex_id, grid, export_state.nto.idx, export_state.nto.lambda_idx, export_state.nto_component, op);
                        break;
                    }
                    case ElectronicStructureSource::TransitionDensity:
                    {
                        evaluate_transition_density(state, vol.tex_id, grid, export_state.nto.idx, export_state.transition_density_component);
                        break;
                    }
                    case ElectronicStructureSource::ElectronDensity:
                        evaluate_electron_density(state, vol.tex_id, grid, export_state.spin, es_gto_op(export_state.spin == ElectronicStructureSpin::Difference && export_state.use_magnitude));
                        break;
                    default:
                        MD_LOG_ERROR("Unsupported export type");
                        goto done;
                    }

                    int num_samples = vol.dim[0] * vol.dim[1] * vol.dim[2];

                    // The geometry written into the cube and xyz headers is the one the calculation
                    // ran at, which is the QM domain's own - not wherever the system's atoms are now.
                    QmAtoms qm = {};
                    es_qm_atoms(&qm, state.mold.sys);

                    size_t              natoms = qm.count;
                    const dvec3_t* qm_coords  = qm.coordinate;
                    const uint8_t* qm_numbers = qm.atomic_number;
                    float* data = md_temp_alloc_array(temp, float, num_samples);

                    vec3_t origin = grid.origin;
                    vec3_t step_x = grid.orientation[0] * grid.spacing[0];
                    vec3_t step_y = grid.orientation[1] * grid.spacing[1];
                    vec3_t step_z = grid.orientation[2] * grid.spacing[2];

                    // Readbacks are queued, not completed, by the evaluate_* calls above.
                    // Exporting reads the texture directly, so it is one of the few places that
                    // genuinely needs the data now rather than next frame.
                    gpu_volume_jobs_drain(&state);

                    // Extract data from OpenGL Texture
                    glBindTexture(GL_TEXTURE_3D, vol.tex_id);
                    glGetTexImage(GL_TEXTURE_3D, 0, GL_RED, GL_FLOAT, data);
                    glBindTexture(GL_TEXTURE_3D, 0);

                    switch (file_format) {
                    case EXPORT_FILE_FORMAT_CUBE: {
                        md_file_printf(file, "Cube file generated by VIAMD\n");
                        if (export_state.source == ElectronicStructureSource::ElectronDensity) {
                            md_file_printf(file, "%s %s %s\n", electronic_structure_source_str[(int)export_state.source], electronic_structure_value_mode_str(export_state.source, export_state.spin, export_state.use_magnitude), electronic_structure_spin_str[(int)export_state.spin]);
                        } else {
                            md_file_printf(file, "%s %s\n", electronic_structure_source_str[(int)export_state.source], electronic_structure_value_mode_str(export_state.source, export_state.spin, export_state.use_magnitude));
                        }

                        md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", -(int)natoms, origin.x, origin.y, origin.z);
                        md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol.dim[0], step_x.x, step_x.y, step_x.z);
                        md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol.dim[1], step_y.x, step_y.y, step_y.z);
                        md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol.dim[2], step_z.x, step_z.y, step_z.z);

                        for (size_t i = 0; i < natoms; ++i) {
                            md_file_printf(file, "%5i %12.6f %12.6f %12.6f %12.6f\n", (int)qm_numbers[i], (float)qm_numbers[i], qm_coords[i].x * ANGSTROM_TO_BOHR, qm_coords[i].y * ANGSTROM_TO_BOHR, qm_coords[i].z * ANGSTROM_TO_BOHR);
                        }

                        // This entry somehow relates to the number of densities
                        md_file_printf(file, "%5i %5i\n", 1, orb_idx + 1);

                        // Write density data
                        int count = 0;
                        for (int x = 0; x < vol.dim[0]; ++x) {
                            for (int y = 0; y < vol.dim[1]; ++y) {
                                for (int z = 0; z < vol.dim[2]; ++z) {
                                    int idx = z * vol.dim[0] * vol.dim[1] + y * vol.dim[0] + x;
                                    float val = data[idx];
                                    md_file_printf(file, " %12.6E", val);
                                    if (++count % 6 == 0) md_file_printf(file, "\n");
                                }
                            }
                        }
                        MD_LOG_INFO("Successfully exported electronic structure to '" MD_STR_FMT "'", MD_STR_ARG(path));
                        break;
                    }
                    case EXPORT_FILE_FORMAT_RAW_MHD_XYZ: {
                        /*
                        * # EXAMPLE
                        ObjectType = Image
                        NDims = 3
                        DimSize = 256 256 64
                        ElementType = MET_USHORT
                        HeaderSize = -1
                        ElementSpacing = 1 1 1
                        ElementByteOrderMSB = False
                        ElementDataFile = image.raw
                        */

                        md_str_t basepath;
                        if (!md_extract_file_path_without_ext(&basepath, path)) {
                            MD_LOG_ERROR("Failed to extract base path");
                            goto done;
                        }

                        md_str_t filename_raw;
                        if (!md_extract_file(&filename_raw, path)) {
                            MD_LOG_ERROR("Failed to extract filename raw");
                            goto done;
                        }

                        md_str_t filename;
                        if (!md_extract_file_path_without_ext(&filename, filename_raw)) {
                            MD_LOG_ERROR("Failed to extract filename");
                            goto done;
                        }

                        md_file_write(file, data, num_samples * sizeof(float));

                        md_allocator_i* temp_arena = md_temp_allocator(temp);
                        md_str_t mhd_path = md_str_printf(temp_arena, MD_STR_FMT ".mhd", MD_STR_ARG(basepath));
                        md_str_t xyz_path = md_str_printf(temp_arena, MD_STR_FMT ".xyz", MD_STR_ARG(basepath));

                        md_file_t mhd_file = { 0 };
                        if (!md_file_open(&mhd_file, mhd_path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                            MD_LOG_ERROR("Failed to open .mhd file");
                            goto done;
                        }

                        vec3_t pos   = origin * BOHR_TO_ANGSTROM;
                        vec3_t step  = vol.voxel_size * BOHR_TO_ANGSTROM;

                        md_file_printf(mhd_file, "ObjectType = Image\n");
                        md_file_printf(mhd_file, "NDims = 3\n");
                        md_file_printf(mhd_file, "DimSize = %i %i %i\n", vol.dim[0], vol.dim[1], vol.dim[2]);
                        md_file_printf(mhd_file, "ElementType = MET_FLOAT\n");
                        md_file_printf(mhd_file, "HeaderSize = -1\n");
                        md_file_printf(mhd_file, "ElementSpacing = %12.6f %12.6f %12.6f\n", step.x, step.y, step.z);
                        md_file_printf(mhd_file, "ElementByteOrderMSB = False\n");
                        md_file_printf(mhd_file, "ElementDataFile = " MD_STR_FMT "\n", MD_STR_ARG(filename_raw));
                        md_file_printf(mhd_file, "Position = %12.6f %12.6f %12.6f\n", pos.x, pos.y, pos.z);
                        md_file_printf(mhd_file, "Orientation =\n");
                        md_file_printf(mhd_file, "%12.6f %12.6f %12.6f\n", grid.orientation[0].x, grid.orientation[0].y, grid.orientation[0].z);
                        md_file_printf(mhd_file, "%12.6f %12.6f %12.6f\n", grid.orientation[1].x, grid.orientation[1].y, grid.orientation[1].z);
                        md_file_printf(mhd_file, "%12.6f %12.6f %12.6f\n", grid.orientation[2].x, grid.orientation[2].y, grid.orientation[2].z);

                        md_file_close(&mhd_file);

                        md_file_t xyz_file = { 0 };
                        if (!md_file_open(&xyz_file, xyz_path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                            MD_LOG_ERROR("Failed to open .xyz file");
                            goto done;
                        }
                        MD_DEFER { md_file_close(&xyz_file); };
                        
                        // XYZ
                        md_file_printf(xyz_file, "%zu\n", natoms);
                        const FileEntry* system_file = file_list_system_file(&state.files.list);
                        md_file_printf(xyz_file, "Geometry extracted from: '%s'\n", system_file ? system_file->path : "");
                        for (size_t i = 0; i < natoms; ++i) {
                            md_str_t sym = md_util_element_symbol(qm_numbers[i]);
                            md_file_printf(xyz_file, "%-2s %12.6f %12.6f %12.6f\n", sym.ptr, qm_coords[i].x, qm_coords[i].y, qm_coords[i].z);
                        }

                        MD_LOG_INFO("Successfully exported electronic structure to '" MD_STR_FMT "'", MD_STR_ARG(path));

                        break;
                    }
                    default:
                        MD_ASSERT(false);
                        break;
                    }
                }
            }
        }
        done:
        ImGui::End();
    }

    // Calculates the transition matrix heuristic
    static inline void compute_transition_matrix(float* out_matrix, const size_t num_groups, const float* hole_charges, const float* part_charges) {
        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        if (num_groups == 0) return;

        int* donors         = md_temp_alloc_zero_array(temp, int, num_groups);
        int* acceptors      = md_temp_alloc_zero_array(temp, int, num_groups);
        double* charge_diff = md_temp_alloc_zero_array(temp, double, num_groups);

        int num_donors = 0;
        int num_acceptors = 0;

        // Sanitize inputs: clamp small negative numerical noise to zero before normalization.
        double hole_sum = 0;
        double part_sum = 0;

        for (size_t i = 0; i < num_groups; i++) {
            const double h = MD_MAX(0.0, (double)hole_charges[i]);
            const double p = MD_MAX(0.0, (double)part_charges[i]);
            hole_sum += h;
            part_sum += p;
        }

        // Guard against degenerate inputs (no charge to attribute).
        if (hole_sum <= 0.0 || part_sum <= 0.0) {
            md_memset(out_matrix, 0, sizeof(float) * num_groups * num_groups);
            return;
        }

        for (size_t i = 0; i < num_groups; i++) {
            // percentages (clamped, normalized so each sums to 1)
            double gsCharge = MD_MAX(0.0, (double)hole_charges[i]) / hole_sum;
            double esCharge = MD_MAX(0.0, (double)part_charges[i]) / part_sum;

            if (gsCharge > esCharge) {
                donors[num_donors++] = (int)i;
            } else {
                acceptors[num_acceptors++] = (int)i;
            }

            out_matrix[i * num_groups + i] = (float)MD_MIN(gsCharge, esCharge);
            charge_diff[i] = esCharge - gsCharge;
        }

        double total_acceptor_charge = 0;
        for (int i = 0; i < num_acceptors; i++) {
            total_acceptor_charge += MD_MAX(0.0, charge_diff[acceptors[i]]);
        }

        if (total_acceptor_charge <= 0.0) {
            return;
        }

        for (int i = 0; i < num_donors; i++) {
            int donor = donors[i];
            double charge_deficit = MD_MAX(0.0, -charge_diff[donor]);
            for (int j = 0; j < num_acceptors; j++) {
                int acceptor = acceptors[j];
                double acc = MD_MAX(0.0, charge_diff[acceptor]);
                double contrib = charge_deficit * acc / total_acceptor_charge;
                if (contrib < 0.0) contrib = 0.0;
                out_matrix[acceptor * num_groups + donor] = (float)contrib;
            }
        }
    }

    void delete_group(int index) {
        vec4_t deleted_color = nto.group.color[index];
        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
            if (nto.atom_group_idx[i] == (uint32_t)index) {
                nto.atom_group_idx[i] = 0;
            }
            else if (nto.atom_group_idx[i] > (uint32_t)index) {
                nto.atom_group_idx[i] -= 1;
            }
        }
        for (int i = index; i < (int)nto.group.count - 1; ++i) {
            nto.group.color[i] = nto.group.color[i + 1];
            md_memcpy(nto.group.label[i], nto.group.label[i + 1], sizeof(nto.group.label[i]));
        }
        nto.group.count--;
        snprintf(nto.group.label[nto.group.count], sizeof(nto.group.label[nto.group.count]), "Group %i", (int)nto.group.count);
        nto.group.color[nto.group.count] = deleted_color;
    }

    void copy_group(int from_idx, int to_idx) {
        nto.group.color[to_idx] = nto.group.color[from_idx];
        md_memcpy(nto.group.label[to_idx], nto.group.label[from_idx], sizeof(nto.group.label[from_idx]));
        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
            if ((int)nto.atom_group_idx[i] == to_idx) {
                nto.atom_group_idx[i] = (uint32_t)from_idx;
            }
        }
    }

    void swap_groups(int a, int b) {
        vec4_t color = nto.group.color[a];
        nto.group.color[a] = nto.group.color[b];
        nto.group.color[b] = color;

        char label_buf[sizeof(nto.group.label[a])];
        //sprintf(label_buf, nto.group.label[row_n]);
        md_memcpy(label_buf, nto.group.label[a], sizeof(nto.group.label[a])); //Current to buf
        md_memcpy(nto.group.label[a], nto.group.label[b], sizeof(nto.group.label[b])); //Next to current
        md_memcpy(nto.group.label[b], label_buf, sizeof(label_buf)); //Buf to next

        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
            if ((int)nto.atom_group_idx[i] == a) {
                nto.atom_group_idx[i] = (uint32_t)b;
            }
            else if ((int)nto.atom_group_idx[i] == b) {
                nto.atom_group_idx[i] = (uint32_t)a;
            }
        }
    }

    static inline void imgui_delayed_hover_tooltip(const char* text) {
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal)) {
            ImGui::SetTooltip("%s", text);
        }
    }

    void draw_nto_window(ApplicationState& state) {
        if (!nto.show_window) return;

        // Excited states carrying NTOs, from the weight table's leading axis. Cached at load
        // because the menu handler has no system to read it from.
        const size_t num_excited_states = num_excited;
        if (num_excited_states == 0) return;

        bool open_context_menu = false;
        bool export_transition_diagram_requested = false;
        static bool edit_mode = false;

        ImGui::SetNextWindowSize(ImVec2(500, 300), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Transition Analysis", &nto.show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            nto.group.hovered_index = -1;
            bool viewport_hovered = false;

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    export_transition_diagram_requested = ImGui::MenuItem("Export Transition Diagram");
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Settings")) {
                    ImGui::Text("Iso Colors");
                    ImGui::Checkbox("Link Attachment/Detachment Density Colors", &nto.link_attachment_detachment_density);
                    if (nto.link_attachment_detachment_density) {
                        ImGui::ColorEdit4("##Color Density", nto.col_den.elem);
                        ImGui::SetItemTooltip("Color Density");
                    } else {
                        ImGui::ColorEdit4("##Color Attachment", nto.col_att.elem);
                        ImGui::SetItemTooltip("Color Attachment");
                        ImGui::ColorEdit4("##Color Detachment", nto.col_det.elem);
                        ImGui::SetItemTooltip("Color Detachment");
                    }

                    ImGui::ColorEdit4("##Color Positive", nto.col_pos.elem);
                    ImGui::SetItemTooltip("Color Positive");
                    ImGui::ColorEdit4("##Color Negative", nto.col_neg.elem);
                    ImGui::SetItemTooltip("Color Negative");

                    ImGui::Text("Transition Dipole Moments");
                    ImGui::Checkbox("Enabled", &nto.dipole.enabled);
                    ImGui::SliderFloat("Scale", &nto.dipole.vector_scale, 1.0f, 10.0f);
                    ImGui::Checkbox("Show Angle", &nto.dipole.show_angle);
                    ImGui::ColorEdit4("Color Electric", nto.dipole.colors[0].elem);
                    ImGui::SetItemTooltip("Color Electric");
                    ImGui::ColorEdit4("Color Magnetic", nto.dipole.colors[1].elem);
                    ImGui::SetItemTooltip("Color Magnetic");
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }

            const ImVec2 outer_size = {300.f + edit_mode * 80, 0.f};
            ImGui::PushItemWidth(outer_size.x);
            ImGui::BeginGroup();

            ImGui::Text("Transition State Index");
            ImGui::Spacing();
            if (ImGui::BeginListBox("##NTO Index", outer_size)) {
                if (ImGui::IsWindowHovered()) {
                    rsp.hovered = -1;
                }
                for (int i = 0; i < (int)num_excited_states; ++i) {
                    bool is_selected = rsp.selected == i;
                    bool is_hovered  = rsp.hovered  == i;
                    char buf[32];
                    snprintf(buf, sizeof(buf), "%i", i + 1);
                    if (is_hovered) {
                        ImGui::PushStyleColor(ImGuiCol_Header, ImGui::GetColorU32(ImGuiCol_HeaderHovered));
                    }
                    if (ImGui::Selectable(buf, is_selected || is_hovered)) {
                        rsp.selected = i;
                    }
                    if (is_hovered) {
                        ImGui::PopStyleColor();
                    }
                    if (ImGui::IsItemHovered()) {
                        rsp.hovered = i;
                    }
                }
                ImGui::EndListBox();
            }
            ImGui::Spacing();
            const double iso_min = 1.0e-8;
            const double iso_max = 5.0;
            
            ImGui::Spacing();
            ImGui::TextUnformatted(electronic_structure_iso_value_label()); 
            ImGui::SliderScalar("##Iso Value", ImGuiDataType_Double, &nto.iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
            ImGui::SetItemTooltip("Visual surface threshold. Attachment and detachment density views render at isovalue^2 for perceptual consistency with NTO amplitudes.");

            ImGui::Spacing();

            static const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersH | ImGuiTableFlags_SizingFixedFit;
            static const ImGuiTableColumnFlags columns_base_flags = ImGuiTableColumnFlags_NoSort;
            
            static bool hide_overlap_text = true;
            ImVec2 button_size(ImGui::GetFontSize() + ImGui::GetStyle().FramePadding.x, 0.f);

            ImGui::Checkbox("Edit mode", &edit_mode);
            ImGui::SameLine();
            ImGui::Checkbox("Hide overlapping text", &hide_overlap_text);

            int group_counts[MAX_NTO_GROUPS] = {0};
            for (size_t i = 0; i < md_array_size(nto.atom_group_idx); ++i) {
                int group_idx = nto.atom_group_idx[i] < nto.group.count ? nto.atom_group_idx[i] : 0;
                group_counts[group_idx] += 1;
            }

            if (ImGui::BeginTable("Group Table", 3 + edit_mode * 4, flags, outer_size, 0)) {
                ImGui::TableSetupColumn("Group", columns_base_flags, 150.f);
                ImGui::TableSetupColumn("Count", columns_base_flags, 0.0f);
                ImGui::TableSetupColumn("Color", columns_base_flags, 0.0f);
                if (edit_mode) {
                    ImGui::TableSetupColumn("##Delete", columns_base_flags, button_size.x);
                    ImGui::TableSetupColumn("##Move down", columns_base_flags, button_size.x);
                    ImGui::TableSetupColumn("##Move up", columns_base_flags, button_size.x);
                    ImGui::TableSetupColumn("##Merge up", columns_base_flags, button_size.x);
                }

                //ImGui::TableSetupColumn("Coord Y", columns_base_flags, 0.0f);
                //ImGui::TableSetupColumn("Coord Z", columns_base_flags | ImGuiTableColumnFlags_WidthFixed, 0.0f);
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableHeadersRow();

                ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_YELLOW);
                ImGui::PushStyleColor(ImGuiCol_Header, IM_BLUE);
                bool show_unassigned = group_counts[0] > 0;
                int row_n = show_unassigned ? 0 : 1;
                for (; row_n < (int)nto.group.count; row_n++) {
                    ImGui::PushID((int)row_n);
                    ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
                    ImGui::TableNextRow(ImGuiTableRowFlags_None, 0);

                    ImGui::TableNextColumn();
                    ImGui::AlignTextToFramePadding();
                    if (edit_mode) {
                        ImGui::PushItemWidth(150.f);
                        if (row_n > 0) {
                            ImGui::InputText("##label", nto.group.label[row_n], sizeof(nto.group.label[row_n]));
                        } else {
                            ImGui::Text("%s", nto.group.label[0]);
                        }
                        ImGui::PopItemWidth();
                    }
                    else {
                        ImGui::Selectable(nto.group.label[row_n], false, selectable_flags);
                        if (ImGui::TableGetHoveredRow() == (int)(row_n + show_unassigned)) {
                            nto.group.hovered_index = (int8_t)row_n;
                            {
                                md_bitmask_t* hb = &state.selection.highlight_mask;
                                md_bitmask_clear(hb);
                                for (size_t j = 0; j < md_array_size(nto.atom_group_idx); j++) {
                                    if (nto.atom_group_idx[j] == (uint32_t)row_n) {
                                        //Add j to highlight mask
                                        md_bitmask_set_bit(hb, (uint32_t)j); //TODO: Hover only works if you hold leftMouseBtn
                                    }
                                }
                            }

                            //Selection
                            if (ImGui::IsKeyDown(ImGuiKey_MouseLeft) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                md_bitmask_or(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
                            }
                            //Deselect
                            else if (ImGui::IsKeyDown(ImGuiKey_MouseRight) && ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
                                md_bitmask_andnot(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
                            }
                        }
                    }


                    ImGui::TableNextColumn();
                    ImGui::Text("%i", group_counts[row_n]);

                    ImGui::TableNextColumn();

                    //Center the color picker
                    ImVec2 cell_size = ImGui::GetContentRegionAvail();
                    float color_size = ImGui::GetFrameHeight();
                    ImVec2 padding((cell_size.x - color_size) * 0.5f, 0.0f);
                    ImGui::SetCursorPos(ImGui::GetCursorPos() + padding);

                    if (edit_mode && row_n > 0) { //You cannot edit the Unassigned color
                        ImGui::ColorEdit4Minimal("##color", nto.group.color[row_n].elem);
                    } else {
                        ImGui::ColorButton("##color", vec_cast(nto.group.color[row_n]));
                    }

                    if (edit_mode) {
                        if (row_n > 0) {
                            ImGui::TableNextColumn();
                            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x - button_size.x);
                            if (ImGui::DeleteButton(ICON_FA_XMARK, button_size)) {
                                delete_group(row_n);
                            }
                            imgui_delayed_hover_tooltip("Delete this group (related atoms will be unassigned)");

                            //Down Arrow
                            ImGui::TableNextColumn();
                            if (row_n == (int)nto.group.count - 1) {
                                ImGui::Dummy(button_size);
                            }
                            else {
                                if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                                    if (ImGui::Button(ICON_FA_ANGLES_DOWN, button_size)) {
                                        for (int i = row_n; i < (int)nto.group.count - 1; i++) {
                                            swap_groups(i, i + 1);
                                        }
                                    }
                                    imgui_delayed_hover_tooltip("Move to bottom");

                                }
                                else {
                                    if (ImGui::Button(ICON_FA_ANGLE_DOWN, button_size)) {
                                        swap_groups(row_n, row_n + 1);
                                    }
                                    imgui_delayed_hover_tooltip("Move down (Ctrl-click to move to bottom)");
                                }
                            }

                            //Up Arrow
                            ImGui::TableNextColumn();
                            if (row_n == 1) {
                                ImGui::Dummy(button_size);
                            }
                            else {
                                if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                                    if (ImGui::Button(ICON_FA_ANGLES_UP, button_size)) {
                                        int last = (int)nto.group.count - 1;
                                        for (int i = last; i >= (int)row_n - 1; i--) {
                                            swap_groups(i, i - 1);
                                        }
                                    }
                                    imgui_delayed_hover_tooltip("Move to top");
                                }
                                else {
                                    if (ImGui::Button(ICON_FA_ANGLE_UP, button_size)) {
                                        swap_groups(row_n, row_n - 1);
                                    }
                                    imgui_delayed_hover_tooltip("Move up (Ctrl-click to move to top)");
                                }
                            }

                            //Merge with above arrow
                            ImGui::TableNextColumn();
                            if (row_n == 1) {
                                ImGui::Dummy(button_size);
                            }
                            else {
                                if (ImGui::Button(ICON_FA_TRASH_ARROW_UP, button_size)) {
                                    for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
                                        if (nto.atom_group_idx[i] == (uint32_t)row_n) {
                                            nto.atom_group_idx[i] = (uint32_t)(row_n - 1);
                                        }
                                    }
                                    delete_group(row_n);
                                }
                                imgui_delayed_hover_tooltip("Merge with above");
                            }
                        }
                    }
                    ImGui::PopID();
                }

                ImGui::PopStyleColor(2);
                ImGui::EndTable();
            }

            if (ImGui::Button("Add new group")) {
                nto.group.count++;
            }

            ImGui::EndGroup();
            ImGui::PopItemWidth();

            ImGui::SameLine();

            // Animate camera towards targets
            const double dt = state.app.timing.delta_s;
            camera_animate(&nto.camera, nto.target, dt);

            ImVec2 canvas_sz = ImGui::GetContentRegionAvail();   // Resize canvas to what's available
            canvas_sz.x = MD_MAX(canvas_sz.x, 50.0f);
            canvas_sz.y = MD_MAX(canvas_sz.y, 50.0f);

            if (export_transition_diagram_requested) {
                ImVec2 export_sz = { MD_MAX(canvas_sz.x * 0.5f, 1.0f), MD_MAX(canvas_sz.y, 1.0f) };
                export_sz *= 2.0f;
                export_transition_diagram(export_sz, &nto, hide_overlap_text);
            }

            ImGui::Dummy(canvas_sz);

            // Draw border and background color
            ImGuiIO& io = ImGui::GetIO();

            ImVec2 canvas_p0 = ImGui::GetItemRectMin();
            ImVec2 canvas_p1 = ImGui::GetItemRectMax();

            bool reset_view = false;

            if (rsp.selected != -1) {
                if (nto.sel_nto_idx != rsp.selected) {
                    nto.sel_nto_idx = rsp.selected;
                    size_t nto_idx = (size_t)rsp.selected;

                    init_volume(&nto.vol[NTO_Attachment], nto.grid, GL_R32F);
                    init_volume(&nto.vol[NTO_Detachment], nto.grid, GL_R32F);

                    evaluate_transition_density(state, nto.vol[NTO_Attachment].tex_id, nto.grid, nto_idx, ElectronicStructureTransitionDensityComponent::Attachment);
                    evaluate_transition_density(state, nto.vol[NTO_Detachment].tex_id, nto.grid, nto_idx, ElectronicStructureTransitionDensityComponent::Detachment);
                }
            }

            const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();

            ImVec2 grid_p0 = canvas_p0;
            ImVec2 grid_p1 = canvas_p0 + canvas_sz * ImVec2(0.5f, 1.0f);
            ImVec2 win_sz = (grid_p1 - grid_p0) / ImVec2(1.0f, 2.0f);
            win_sz.x = floorf(win_sz.x);
            win_sz.y = floorf(win_sz.y);

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->ChannelsSplit(2);
            draw_list->ChannelsSetCurrent(0);
            draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));

            const float aspect_ratio = win_sz.x / win_sz.y;

            mat4_t view_mat     = camera_world_to_view_matrix(nto.camera);
            mat4_t proj_mat     = camera_view_to_clip_matrix_persp(nto.camera, aspect_ratio);
            mat4_t inv_view_mat = camera_view_to_world_matrix(nto.camera);
            mat4_t inv_proj_mat = camera_clip_to_view_matrix_persp(nto.camera, aspect_ratio);

            mat4_t world_to_clip = proj_mat * view_mat;
            mat4_t clip_to_world = inv_view_mat * inv_proj_mat;

            if (rsp.selected != -1) {
                const int nto_target_idx[2] = {
                    NTO_Attachment,
                    NTO_Detachment,
                };

                // The render targets, sized before the panels go into the draw list: it records the texture
                // names, and a texture reallocated later in the frame would be drawn under a dead name
                int width  = MD_MAX(1, (int)(win_sz.x * io.DisplayFramebufferScale.x));
                int height = MD_MAX(1, (int)(win_sz.y * io.DisplayFramebufferScale.y));

                auto& gbuf = nto.gbuf;
                if ((int)gbuf.width != width || (int)gbuf.height != height) {
                    gbuffer_init(&gbuf, width, height);
                    for (int i = 0; i < (int)MD_COUNTOF(nto.iso_tex); ++i) {
                        if (nto.iso_tex[i]) gl::free_texture(&nto.iso_tex[i]);
                        if (nto.out_tex[i]) gl::free_texture(&nto.out_tex[i]);
                    }
                }
                for (int i : nto_target_idx) {
                    if (!nto.iso_tex[i]) gl::init_texture_2D(&nto.iso_tex[i], width, height, GL_RGBA16F);
                    if (!nto.out_tex[i]) gl::init_texture_2D(&nto.out_tex[i], width, height, GL_RGBA8);
                }

                // Draw Attachment / Detachment
                for (int i = 0; i < 2; ++i) {
					ImGui::PushID(i);
					MD_DEFER { ImGui::PopID(); };

                    int idx = nto_target_idx[i];

                    ImVec2 p0 = grid_p0 + win_sz * ImVec2(0.0f, (float)(i+0));
                    ImVec2 p1 = grid_p0 + win_sz * ImVec2(1.0f, (float)(i+1));
                    ImRect rect = {p0, p1};
                    ImGui::SetCursorScreenPos(p0);

                    draw_list->PushClipRect(p0, p1);
                    draw_list->ChannelsSetCurrent(1);
                    ImGui::PushID(i);
                    MD_DEFER {
                        draw_list->PopClipRect();
                        ImGui::PopID();
                    };

                    InteractionSurfaceState surface = interaction_surface(interaction_surface_nto, vec_cast(rect.GetSize()));

                    InteractionSurfaceViewTransformArgs view_args = {
                        .camera = nto.camera,
                        .trackball_param = state.view.trackball_param,
                    };
                    InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&nto.target, surface, view_args);
                    if (view_result.reset_requested) {
                        nto.target = default_view;
                    }

                    if (surface.hovered) {
                        InteractionSurfaceHitArgs hit_args = {
                            .picking_surface = &nto.picking_surface,
                            .picking_handler = state.picking_handler,
                            .fbo = nto.gbuf.fbo,
                            .width = nto.gbuf.width,
                            .height = nto.gbuf.height,
                            .clip_to_world = clip_to_world,
                        };

                        PickingHit hit = {};
                        interaction_surface_hit_extract(&hit, surface, hit_args);

                        InteractionSurfaceEvent event = {};
                        interaction_surface_event_extract(&event, surface, hit);

                        event.clip_to_world = clip_to_world;
                        event.world_to_clip = world_to_clip;

                        if (event.kind == InteractionSurfaceEventKind::ContextMenu) {
                            open_context_menu = true;
                        } else if (event.kind == InteractionSurfaceEventKind::RegionSelect) {
                            const md_bitset_t selected = md_bitmask_view(&state.selection.selection_mask);
                            const md_bitset_t* candidate_mask = nullptr;
                            if (event.selection_mode == InteractionSelectionMode::Remove) {
                                // When removing, only consider currently selected atoms as candidates for region selection
                                candidate_mask = &selected;
                            }
                            point_set_region_mask_compute(&state.selection.highlight_mask,
                                state.mold.state.xyz, state.mold.state.num_atoms,
                                candidate_mask, event.world_to_clip, event.region_min, event.region_max, event.surface_size);

                            grow_mask_by_selection_granularity(&state.selection.highlight_mask, state.selection.granularity, state.mold.sys);
                            if (event.region_phase == InteractionSurfaceEventPhase::Commit) {
                                // Merge highlight into selection
                                if (event.selection_mode == InteractionSelectionMode::Append) {
                                    md_bitmask_or(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
                                }
                                else if (event.selection_mode == InteractionSelectionMode::Remove) {
                                    md_bitmask_andnot(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
                                }
                                md_bitmask_clear(&state.selection.highlight_mask);
                            }
                        } else {
                            viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
                        }
                    }

                    draw_list->ChannelsSetCurrent(0);

                    //ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    ImVec2 text_pos_tl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p0.y + TEXT_BASE_HEIGHT * 0.5f);
                    const char* lbl = ((i & 1) == 0) ? "Attachment" : "Detachment";

                    //char buf[32];
                    //snprintf(buf, sizeof(buf), (const char*)u8"λ: %.3f", nto_lambda[i / 2]);
                    draw_list->AddImage((ImTextureID)(intptr_t)nto.out_tex[idx], p0, p1, { 0,1 }, { 1,0 });
                    //draw_list->AddText(text_pos_bl, ImColor(0,0,0), buf);
                    draw_list->AddText(text_pos_tl, ImColor(0,0,0), lbl);

                    if (nto.dipole.enabled) {
                        const vec3_t mid = vec3_lerp(aabb.min_ext, aabb.max_ext, 0.5f);

                        // Out of the dipole groups in the table, not out of a reader. Each is one
                        // vector per excited state, so the selection indexes them directly - and the
                        // count is checked, because a selection can outlive the file it was made in.
                        size_t num_magn = 0;
                        size_t num_elec = 0;
                        const dvec3_t* magnetic_dp = attribute_vec3_series(&num_magn, state.mold.sys, MD_STR_LIT("dipole/magnetic_transition/vector"));
                        const dvec3_t* electric_dp = attribute_vec3_series(&num_elec, state.mold.sys, MD_STR_LIT("dipole/electric_transition/vector"));

                        const bool has_dipoles = magnetic_dp && electric_dp && rsp.selected >= 0 &&
                                                 (size_t)rsp.selected < num_magn && (size_t)rsp.selected < num_elec;
                        const dvec3_t magn_raw = has_dipoles ? magnetic_dp[rsp.selected] : dvec3_t{0, 0, 0};
                        const dvec3_t elec_raw = has_dipoles ? electric_dp[rsp.selected] : dvec3_t{0, 0, 0};

                        vec3_t magn_vec = {(float)magn_raw.x, (float)magn_raw.y, (float)magn_raw.z};
                        vec3_t elec_vec = {(float)elec_raw.x, (float)elec_raw.y, (float)elec_raw.z};

                        magn_vec *= (float)(nto.dipole.vector_scale * BOHR_TO_ANGSTROM);
                        elec_vec *= (float)(nto.dipole.vector_scale * BOHR_TO_ANGSTROM);

                        auto proj_point = [world_to_clip, win_sz, p0](vec3_t point) -> ImVec2 {
                            const vec4_t p = mat4_mul_vec4(world_to_clip, vec4_from_vec3(point, 1.0));
                            const ImVec2 c = {
                                p0.x + ( p.x / p.w * 0.5f + 0.5f) * win_sz.x,
                                p0.y + (-p.y / p.w * 0.5f + 0.5f) * win_sz.y,
                            };
                            return c;
                        };

                        auto draw_arrow = [draw_list](ImVec2 beg, ImVec2 end, ImU32 color, float thickness) {
                            float len2 = ImLengthSqr(end-beg);
                            if (len2 > 1.0e-3f) {
                                float len = sqrtf(len2);
                                ImVec2 d = (end - beg) / len;
                                ImVec2 o = {-d.y, d.x};
                                float d_scl = ImMin(thickness * 3, len);
                                float o_scl = 0.57735026918962576451f * d_scl;
                                draw_list->AddLine(beg, end - d * thickness, color, thickness);
                                draw_list->AddTriangleFilled(end, end - d * d_scl + o * o_scl, end - d * d_scl - o * o_scl, color);
                            }
                        };

                        const ImU32 elec_color = convert_color(nto.dipole.colors[0]);
                        const ImU32 magn_color = convert_color(nto.dipole.colors[1]);
                        const ImU32 text_color = IM_COL32_BLACK;
                        const float arrow_thickness = 3.0f;

                        ImVec2 c      = proj_point(mid);
                        ImVec2 c_elec = proj_point(mid + elec_vec);
                        ImVec2 c_magn = proj_point(mid + magn_vec);

                        draw_arrow(c, c_elec, elec_color, arrow_thickness);
                        draw_list->AddText(c_elec, text_color, (const char*)u8"μe");

                        draw_arrow(c, c_magn, magn_color, arrow_thickness);
                        draw_list->AddText(c_magn, text_color, (const char*)u8"μm");

                        if (nto.dipole.show_angle) {
                            const vec3_t magn_dir = vec3_normalize(magn_vec);
                            const vec3_t elec_dir = vec3_normalize(elec_vec);
                            const float angle = acosf(vec3_dot(magn_dir, elec_dir));
                            const ImU32 angle_color = IM_COL32(0, 0, 0, 128);
                            const float angle_thickness = 1.0f;
                            const float angle_vec_scale = 0.1f;

                            const size_t num_seg = 20;
                            const vec3_t axis = vec3_normalize(vec3_cross(elec_dir, magn_dir));

                            for (size_t seg = 0; seg <= num_seg; ++seg) {
                                float  t = (float)(seg) / (float)num_seg;
                                vec3_t v = quat_mul_vec3(quat_axis_angle(axis, angle * t), elec_dir);
                                ImVec2 p = proj_point(mid + v * angle_vec_scale);
                                draw_list->PathLineTo(p);
                            }
                            draw_list->PathStroke(angle_color, 0, angle_thickness);

                            char buf[32];
                            snprintf(buf, sizeof(buf), (const char*)u8"%.2f°", MD_RAD_TO_DEG(angle));
                            draw_list->AddText(c, text_color, buf);
                        }
                    }
                }
                // Draw Sankey Diagram of Transition Matrix
                {
                    //compute_transition_matrix(nto.transition_matrix, nto.group.count, nto.transition_density_hole, nto.transition_density_part);
                    ImVec2 p0 = canvas_p0 + canvas_sz * ImVec2(0.5f, 0.0f);
                    ImVec2 p1 = canvas_p1;
                    im_sankey_diagram(&state, {p0.x, p0.y, p1.x, p1.y}, &nto, hide_overlap_text);
                    ImVec2 text_pos_bl = ImVec2(p0.x + TEXT_BASE_HEIGHT * 0.5f, p1.y - TEXT_BASE_HEIGHT);
                    draw_list->AddText(text_pos_bl, ImColor(0, 0, 0, 255), "Transition Diagram");
                }
                // Draw grid
                {
                    ImVec2 p0 = {floorf(canvas_p0.x + canvas_sz.x * 0.5f), canvas_p0.y};
                    ImVec2 p1 = {floorf(canvas_p0.x + canvas_sz.x * 0.5f), canvas_p1.y};
                    draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255));
                }

                int num_subwindows = 2;
                for (int i = 1; i < num_subwindows; ++i) {
                    float y = floorf(canvas_p0.y + canvas_sz.y / ((float)num_subwindows) * i);
                    float x0 = canvas_p0.x;
                    float x1 = floorf(canvas_p0.x + canvas_sz.x * (i & 1 ? 0.5f : 1.0f));
                    draw_list->AddLine({x0, y}, {x1, y}, IM_COL32(0, 0, 0, 255));
                }

                // Draw stuff

                const bool is_hovered = ImGui::IsItemHovered();
                const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
                const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x, io.MousePos.y - origin.y);

                if (is_hovered) {
                    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                        reset_view = true;
                    }
                }

                if (reset_view) {
					vec3_t center = vec3_lerp(aabb.min_ext, aabb.max_ext, 0.5f) * BOHR_TO_ANGSTROM;
					vec3_t half_ext = (aabb.max_ext - aabb.min_ext) * 0.5f * BOHR_TO_ANGSTROM;
                    nto.target = compute_optimal_view(center, half_ext, oabb.orientation, nto.distance_scale);
                }

                if (nto.show_coordinate_system_widget) {
                    float ext = MD_MIN(win_sz.x, win_sz.y) * 0.3f;
                    float pad = 20.0f;
                    ImVec2 size = ImVec2(ext, ext);
                    ImGui::SetCursorScreenPos(canvas_p0 + ImVec2(pad, 2.0f * win_sz.y - ext - pad));

                    quat_t out_orientation = nto.target.orientation;
                    if (ImGui::CoordinateSystemWidget(&out_orientation, nto.camera.orientation, size)) {
                        const vec3_t look_at = camera_get_look_at(nto.target);
                        nto.target.orientation = quat_normalize(out_orientation);
                        nto.target.position = camera_position_from_look_at(look_at, nto.target.orientation, nto.target.distance);
                    }
                }

                // The mesh rendering mode: the attachment and detachment densities as fields, extracted together
                const bool mesh_mode = state.settings.isosurface_rendering == IsosurfaceRendering_Mesh;
                if (mesh_mode) {
                    const uint64_t pos_hash = panel_positions_hash(state);
                    for (int i : nto_target_idx) {
                        const int id[2] = {nto.sel_nto_idx, i};
                        uint64_t key = md_hash64(id, sizeof(id), pos_hash);
                        key = key ? key : 1;
                        panel_meshes_field(&nto.mesh, &state, i, key, [&]() -> md_gto_iso_field_t* {
                            if (!td_cache_ensure(state.mold.sys, nto.sel_nto_idx)) return nullptr;
                            const double* D = (i == NTO_Attachment) ? td_cache.attach : td_cache.detach;
                            return iso_field_create_density_matrix(&state, D, td_cache.dim);
                        });
                    }
                    panel_meshes_update(&nto.mesh, &state, 2, nto.iso_val * nto.iso_val, MD_GTO_ISO_SHEET_POSITIVE);
                } else if (nto.mesh.meshes.key) {
                    panel_meshes_free(&nto.mesh, &state);
                }

                // The panels are only redrawn when something they show changed, or while their meshes are
                // refined over the frames: a static view costs nothing
                uint64_t render_hash = md_hash64(&view_mat, sizeof(view_mat), (uint64_t)nto.gl_rep.id);
                {
                    render_hash = md_hash64(&state.settings.isosurface_rendering, sizeof(state.settings.isosurface_rendering), render_hash);
                    if (mesh_mode) {
                        render_hash = md_hash64(&nto.mesh.meshes.mesh.version, sizeof(nto.mesh.meshes.mesh.version), render_hash);
                    }
                    const int layout[2] = { width, height };
                    render_hash = md_hash64(layout, sizeof(layout), render_hash);
                    for (int i : nto_target_idx) {
                        const uint64_t vol_version = volume::data_version(nto.vol[i].tex_id);
                        render_hash = md_hash64(&i, sizeof(i), render_hash);
                        render_hash = md_hash64(&nto.vol[i].tex_id, sizeof(nto.vol[i].tex_id), render_hash);
                        render_hash = md_hash64(&vol_version, sizeof(vol_version), render_hash);
                        render_hash = md_hash64(&nto.vol[i].texture_to_world, sizeof(nto.vol[i].texture_to_world), render_hash);
                    }
                    render_hash = md_hash64(&proj_mat, sizeof(proj_mat), render_hash);
                    render_hash = md_hash64(&nto.iso_val, sizeof(nto.iso_val), render_hash);
                    render_hash = md_hash64(&state.settings.exact_isosurfaces, sizeof(state.settings.exact_isosurfaces), render_hash);
                    const vec4_t iso_colors[5] = { nto.col_pos, nto.col_neg, nto.col_den, nto.col_att, nto.col_det };
                    render_hash = md_hash64(iso_colors, sizeof(iso_colors), render_hash);
                    render_hash = md_hash64(&nto.link_attachment_detachment_density, sizeof(nto.link_attachment_detachment_density), render_hash);
                    render_hash = md_hash64(&state.visuals.background, sizeof(state.visuals.background), render_hash);
                    render_hash = md_hash64(&state.visuals.tonemapping, sizeof(state.visuals.tonemapping), render_hash);
                    render_hash = md_hash64(&state.visuals.shading, sizeof(state.visuals.shading), render_hash);
                    render_hash = md_hash64(&state.selection.color, sizeof(state.selection.color), render_hash);
                    render_hash = md_hash64(&state.mold.gpu_buffers_version, sizeof(state.mold.gpu_buffers_version), render_hash);
                    // Baked into the picking buffer, which hovering reads
                    const md_gl_picking_offset_t picking_offset = gl_picking_offset(state);
                    render_hash = md_hash64(&picking_offset, sizeof(picking_offset), render_hash);
                    render_hash = render_hash ? render_hash : 1;
                }

                if (render_hash != nto.render_hash || (mesh_mode && !nto.mesh.converged)) {
                    nto.render_hash = render_hash;

                    gbuffer_clear(&gbuf);

                    const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                        GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

                    glEnable(GL_CULL_FACE);
                    glCullFace(GL_BACK);

                    glEnable(GL_DEPTH_TEST);
                    glDepthFunc(GL_LESS);
                    glDepthMask(GL_TRUE);
                    glEnable(GL_SCISSOR_TEST);

                    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
                    glDrawBuffers((int)MD_COUNTOF(draw_buffers), draw_buffers);
                    glViewport(0, 0, gbuf.width, gbuf.height);
                    glScissor(0, 0, gbuf.width, gbuf.height);

                    auto draw_rep = [&state](md_gl_rep_t& rep, md_gl_shaders_t& shaders, mat4_t& view_mat, mat4_t& proj_mat, uint32_t atom_mask = 0) {
                        md_gl_draw_op_t draw_op = {};
                        draw_op.type = MD_GL_REP_BALL_AND_STICK;
                        draw_op.args.ball_and_stick.ball_scale   = 1.0f;
                        draw_op.args.ball_and_stick.stick_radius = 1.0f;
                        draw_op.args.ball_and_stick.bond_orders  = true;
                        draw_op.rep = rep;

                        md_gl_draw_args_t draw_args = {
                            .shaders = shaders,
                            .draw_operations = {
                                .count = 1,
                                .ops = &draw_op
                            },
                            .view_transform = {
                                .view_matrix = (const float*)view_mat.elem,
                                .proj_matrix = (const float*)proj_mat.elem,
                            },

                            .picking_offset = gl_picking_offset(state),

                            .atom_mask = atom_mask,
                        };

                        md_gl_draw(&draw_args);
                    };

                    draw_rep(nto.gl_rep, state.gl.shaders, view_mat, proj_mat);

                    // The marks of the selection and the highlight, put over each panel as it is presented
                    const bool draw_marks = !md_bitmask_empty(&state.selection.selection_mask) || !md_bitmask_empty(&state.selection.highlight_mask);
                    if (draw_marks) {
                        PUSH_GPU_SECTION("Marks")
                        gbuffer_marks_begin(&gbuf);
                        draw_rep(nto.gl_rep, state.gl.shaders_marks, view_mat, proj_mat, AtomBit_Marked);
                        gbuffer_marks_end(&gbuf);
                        POP_GPU_SECTION()
                    }

                    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
                    glViewport(0, 0, gbuf.width, gbuf.height);
                    glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
                    glClearColor(0, 0, 0, 0);   // premultiplied: nothing
                    glClear(GL_COLOR_BUFFER_BIT);
                    glDepthFunc(GL_LESS);
                    glDepthMask(0);
                    glColorMask(1,1,1,1);

                    PUSH_GPU_SECTION("Postprocessing")
                    postprocess_pipeline::Settings postprocess_settings = {};
                    postprocess_pipeline::Inputs postprocess_inputs = {};

                    postprocess_settings.background_color = state.visuals.background.color;
                    postprocess_settings.shading.mode = state.visuals.shading.mode;
                    postprocess_settings.shading.ambient = state.visuals.shading.ambient;
                    postprocess_settings.shading.specular = state.visuals.shading.specular;
                    postprocess_settings.shading.colors = state.visuals.shading.colors;
                    postprocess_settings.shading.light = state.visuals.shading.light;
                    postprocess_settings.shading.environment = state.visuals.shading.environment;
                    postprocess_settings.tonemap.enabled = state.visuals.tonemapping.enabled;
                    postprocess_settings.tonemap.mode = state.visuals.tonemapping.tonemapper;
                    postprocess_settings.tonemap.exposure = state.visuals.tonemapping.exposure;
                    postprocess_settings.tonemap.gamma = state.visuals.tonemapping.gamma;
                    postprocess_settings.shadows.enabled = false;
                    postprocess_settings.ssao.enabled = false;
                    postprocess_settings.dof.enabled = false;
                    postprocess_settings.fxaa.enabled = true;
                    postprocess_settings.taa.enabled = false;
                    postprocess_settings.sharpen.enabled = false;

                    postprocess_inputs.depth = nto.gbuf.tex.depth;
                    postprocess_inputs.color = nto.gbuf.tex.color;
                    postprocess_inputs.normal = nto.gbuf.tex.normal;
                    postprocess_inputs.velocity = nto.gbuf.tex.velocity;
                    postprocess_inputs.transparency = nto.gbuf.tex.transparency;
                    if (draw_marks) {
                        postprocess_inputs.marks = nto.gbuf.tex.marks;
                        postprocess_inputs.marks_group = nto.gbuf.tex.marks_group;
                        marks_settings(&postprocess_settings, state, 1.0f, MD_MAX(1, (int)roundf(ImGui::GetIO().DisplayFramebufferScale.x)));
                    }

                    ViewParam view_param = {
                        .matrix = {
                            .curr = {
                            .view = view_mat,
                            .proj = proj_mat,
                            .norm = view_mat,
                        },
                        .inv = {
                            .proj = inv_proj_mat,
                        }
                        },
                        .clip_planes = {
                            .near = nto.camera.near_plane,
                            .far  = nto.camera.far_plane,
                        },
                        .resolution = {win_sz.x, win_sz.y},
                        .fov_y = nto.camera.fov_y,
                    };
                    POP_GPU_SECTION()

                    glDisable(GL_SCISSOR_TEST);
                    if (!nto.out_fbo) glGenFramebuffers(1, &nto.out_fbo);

                    // Every panel is its own density's isosurfaces over the opaque image (and selection) the two
                    // share, put through post-processing on its own: the surfaces are composited in HDR
                    PUSH_GPU_SECTION("NTO PANELS")
                    bool converged = true;
                    for (int i : nto_target_idx) {
                        bool is_density = (i == NTO_Attachment || i == NTO_Detachment);
                    
                        bool enabled = true;
                        bool iso_written = false;
                        size_t count = is_density ? 1 : 2;
                        float  values[2];
                        vec4_t colors[2];
                    
                        if (is_density) {
                            if (nto.link_attachment_detachment_density) {
                                colors[0] = nto.col_den;
                            } else {
                                if (i == NTO_Attachment) {
                                    colors[0] = nto.col_att;
                                } else {
                                    colors[0] = nto.col_det;
                                }
                            }
                            values[0] = (float)(nto.iso_val * nto.iso_val);
                        } else {
                            colors[0] = nto.col_pos;
                            colors[1] = nto.col_neg;
                            values[0] =  (float)nto.iso_val;
                            values[1] = -(float)nto.iso_val;
                        }

                        if (enabled && mesh_mode) {
                            const vec4_t sheet_colors[2] = {colors[0], colors[0]};
                            isomesh::Surface surfaces[2];
                            const size_t num_surfaces = panel_meshes_surfaces(surfaces, nto.mesh, i, sheet_colors);
                            const isomesh::RenderDesc mesh_desc = {
                                .render_target = {
                                    .depth  = nto.gbuf.tex.depth,
                                    .color  = nto.iso_tex[i],
                                    .width  = nto.gbuf.width,
                                    .height = nto.gbuf.height,
                                    .clear_color = true,
                                },
                                .matrix = {
                                    .view = view_mat,
                                    .proj = proj_mat,
                                    .proj_no_jitter = proj_mat,
                                },
                                .surfaces = surfaces,
                                .num_surfaces = num_surfaces,
                                .shading = {
                                    .roughness = 0.3f,
                                    .ior = 1.5f,
                                    .mode = postprocess_settings.shading.mode,
                                    .ambient = postprocess_settings.shading.ambient,
                                    .specular = postprocess_settings.shading.specular,
                                    .colors = postprocess_settings.shading.colors,
                                    .light = postprocess_settings.shading.light,
                                    .environment = postprocess_settings.shading.environment,
                                },
                                .scene_hash = render_hash,
                            };
                            bool panel_converged = true;
                            iso_written = isomesh::render(&nto.mesh.target[i], mesh_desc, &panel_converged);
                            converged &= panel_converged;
                        } else if (enabled) {
                            volume::IsoRenderDesc vol_desc = {
                                .render_target = {
                                    .depth  = nto.gbuf.tex.depth,
                                    .color  = nto.iso_tex[i],
                                    .width  = nto.gbuf.width,
                                    .height = nto.gbuf.height,
                                    .clear_color = true,
                                },
                                .texture = {
                                    .density_volume = nto.vol[i].tex_id,
                                },
                                .matrix = {
                                    .model = nto.vol[i].texture_to_world,
                                    .view  = view_mat,
                                    .proj  = proj_mat,
                                    .inv_proj = inv_proj_mat,
                                },
                                .iso = {
                                    .count   = count,
                                    .values  = values,
                                    .colors  = colors,
                                    .exact   = state.settings.exact_isosurfaces,
                                },
                                .shading = {
                                    .roughness = 0.3f,
                                    .ior = 1.5f,
                                    .mode = postprocess_settings.shading.mode,
                                    .ambient = postprocess_settings.shading.ambient,
                                    .specular = postprocess_settings.shading.specular,
                                    .colors = postprocess_settings.shading.colors,
                                    .light = postprocess_settings.shading.light,
                                    .environment = postprocess_settings.shading.environment,
                                },
                            };
                            iso_written = volume::render_isosurfaces(vol_desc);
                        }

                        postprocess_inputs.transparency_hdr = iso_written ? nto.iso_tex[i] : 0;
                        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, nto.out_fbo);
                        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, nto.out_tex[i], 0);
                        glDrawBuffer(GL_COLOR_ATTACHMENT0);
                        glViewport(0, 0, nto.gbuf.width, nto.gbuf.height);
                        postprocess_pipeline::execute(postprocess_inputs, postprocess_settings, view_param);
                    }
                    nto.mesh.converged = converged;
                    POP_GPU_SECTION();

                    // Reset state
                    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
                    glDrawBuffer(GL_BACK);
                    glDisable(GL_SCISSOR_TEST);
                }
            }
            if (ImGui::IsWindowHovered() && nto.group.hovered_index == -1 && !viewport_hovered) {
                md_bitmask_clear(&state.selection.highlight_mask);
            }
        }
        ImGui::End();

        if (open_context_menu) {
            ImGui::OpenPopup("Context Menu");
        }

        if (ImGui::BeginPopup("Context Menu")) {
			size_t count = md_bitmask_count(&state.selection.selection_mask);
			if (count == 0) {
                ImGui::Text("No atoms selected");
            } else {
                if (ImGui::BeginMenu("Assign to group")) {
                    for (size_t i = 0; i < nto.group.count; i++) {
                        if (ImGui::MenuItem(nto.group.label[i])) {
                            for (size_t j = 0; j < md_array_size(nto.atom_group_idx); j++) {
                                //Is the atom selected?
                                if (md_bitmask_test_bit(&state.selection.selection_mask, j)) {
                                    //If so, set its group to the selected group index
                                    nto.atom_group_idx[j] = (uint32_t)i;
                                }
                            }
                        }
                    }
                    ImGui::EndMenu();
                }
                if (nto.group.count < MAX_NTO_GROUPS) {
                    if (ImGui::MenuItem("Assign to new group")) {
                        //Create a new group and add an item to it
                        uint32_t group_idx = (uint32_t)(nto.group.count++);
                        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); i++) {
                            if (md_bitmask_test_bit(&state.selection.selection_mask, i)) {
                                //If so, set its group to the selected group index
                                nto.atom_group_idx[i] = group_idx;
                            }
                        }
                    }
                }
            }
            ImGui::EndPopup();
        }

    }
    // Everything derived from the group assignment: the representation colours and the transition
    // matrix. Lifted out of draw_nto_window so it runs whether or not that window is open - the
    // Charge Transfer window needs the same numbers, and a window should not have to be visible
    // for another one to have data.
    //
    // Driven by rsp.selected rather than nto.sel_nto_idx: the latter tracks which state's VOLUMES
    // are resident, which is a rendering concern and is only updated while the NTO window draws.
    // ---------------------------------------------------------------------------
    // Workspace persistence
    // ---------------------------------------------------------------------------

    void ensure_pending_group_storage() {
        if (pending_groups.atoms_initialized) return;
        for (size_t i = 0; i < MAX_NTO_GROUPS; ++i) {
            pending_groups.atoms[i] = md_bitset_empty_set();    // payloads on the heap
        }
        pending_groups.atoms_initialized = true;
    }

    void serialize_workspace(viamd::serialization_state_t& ser) {
        // Nothing to write without a QM system loaded, and the group assignment is what this
        // section is FOR: num_qm_atoms is "did init_qm_from_system run", which is the same question
        // the reader pointer used to answer and one this component can still ask for itself.
        if (num_qm_atoms == 0 || !nto.atom_group_idx) return;

        // Renamed from [VeloxChem] along with the component. Deliberately NOT accepting the old
        // spelling on read: this section holds nothing but window flags and an NTO group
        // assignment, both of which a user re-makes in a minute, and one name in the file beats a
        // compatibility branch that has to be carried forever.
        viamd::write_section_header(ser, MD_STR_LIT("QuantumChemistry"));

        viamd::write_bool(ser, MD_STR_LIT("ShowTransitionAnalysis"), nto.show_window);
        viamd::write_bool(ser, MD_STR_LIT("ShowChargeTransfer"),     flow.show_window);
        viamd::write_int (ser, MD_STR_LIT("FlowRoute"),              (int64_t)flow.mode);
        viamd::write_int (ser, MD_STR_LIT("FlowRows"),               (int64_t)flow.max_rows);
        viamd::write_flt (ser, MD_STR_LIT("FlowThreshold"),          flow.cut.threshold);
        viamd::write_flt (ser, MD_STR_LIT("FlowOthersCap"),          flow.cut.other_max);
        viamd::write_flt (ser, MD_STR_LIT("FlowSizeNormalize"),      flow.size_normalize);
        viamd::write_int (ser, MD_STR_LIT("FlowOthersExpanded"),     (int64_t)flow.cut.other_expanded);
        viamd::write_flt (ser, MD_STR_LIT("FlowBarWidth"),           flow.params.node_thickness);
        viamd::write_flt (ser, MD_STR_LIT("FlowNodeGap"),            flow.params.node_gap);
        viamd::write_flt (ser, MD_STR_LIT("FlowRibbonOpacity"),      flow.style.ribbon_alpha);
        viamd::write_flt (ser, MD_STR_LIT("FlowRibbonCurvature"),    flow.style.ribbon_curvature);
        viamd::write_int (ser, MD_STR_LIT("FlowOrder"),              (int64_t)flow.params.order);
        viamd::write_int (ser, MD_STR_LIT("FlowSweeps"),             (int64_t)flow.params.barycentre_sweeps);
        viamd::write_flt (ser, MD_STR_LIT("FlowRibbonMinPx"),        flow.style.ribbon_min_px);
        viamd::write_flt (ser, MD_STR_LIT("FlowNodeTint"),           flow.style.node_fill_tint);
        viamd::write_flt (ser, MD_STR_LIT("FlowNodeBorder"),         flow.style.border_width);
        viamd::write_bool(ser, MD_STR_LIT("FlowColorLabels"),        flow.style.color_labels);
        viamd::write_bool(ser, MD_STR_LIT("FlowColorBySide"),        flow.style.color_by_side);
        viamd::write_bool(ser, MD_STR_LIT("FlowRibbonLabels"),       flow.style.show_ribbon_labels);
        viamd::write_bool(ser, MD_STR_LIT("FlowLegend"),             flow.show_legend);

        viamd::write_int(ser, MD_STR_LIT("GroupCount"), (int64_t)nto.group.count);

        // One bitfield per group rather than one index per atom: the same information, but it
        // survives being looked at in a text editor.
        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        md_allocator_i* temp_alloc = md_temp_allocator(temp);

        for (size_t g = 0; g < nto.group.count; ++g) {
            md_bitmask_t b = md_bitmask_create(0, MD_BITSET_MAX_INDEX, temp_alloc);
            for (size_t i = 0; i < md_array_size(nto.atom_group_idx); ++i) {
                if (nto.atom_group_idx[i] == (uint32_t)g) {
                    md_bitmask_set_bit(&b, (uint32_t)i);
                }
            }
            const md_bitset_t mask = md_bitmask_view(&b);
            viamd::write_str(ser, MD_STR_LIT("GroupLabel"), md_str_from_cstr(nto.group.label[g]));
            viamd::write_vec4(ser, MD_STR_LIT("GroupColor"), nto.group.color[g]);
            viamd::write_bitset(ser, MD_STR_LIT("GroupAtoms"), &mask);
        }
    }

    void deserialize_workspace(viamd::deserialization_state_t& deser) {
        ensure_pending_group_storage();

        pending_groups.pending = true;
        pending_groups.count = 0;
        for (size_t i = 0; i < MAX_NTO_GROUPS; ++i) {
            md_bitset_free(&pending_groups.atoms[i], md_get_heap_allocator());
            pending_groups.label[i][0] = '\0';
            pending_groups.color[i] = vec4_t{0.5f, 0.5f, 0.5f, 1.0f};
        }

        // 'GroupLabel' opens a group; the color and atom entries that follow belong to it. The
        // index only advances on a label, so a malformed section cannot walk off the end.
        int64_t group_idx = -1;

        md_str_t ident = {};
        md_str_t arg   = {};
        while (viamd::next_entry(ident, arg, deser)) {
            if (md_str_eq_cstr(ident, "ShowTransitionAnalysis")) {
                viamd::extract_bool(nto.show_window, arg);
            } else if (md_str_eq_cstr(ident, "ShowChargeTransfer")) {
                viamd::extract_bool(flow.show_window, arg);
            } else if (md_str_eq_cstr(ident, "FlowRoute")) {
                int mode = FLOW_MODE_FRAGMENTS;
                if (viamd::extract_int(mode, arg)) {
                    flow.mode = (flow_mode_t)MD_CLAMP(mode, 0, FLOW_MODE_COUNT - 1);
                }
            } else if (md_str_eq_cstr(ident, "FlowRows")) {
                int rows = flow.max_rows;
                if (viamd::extract_int(rows, arg)) flow.max_rows = MD_CLAMP(rows, 2, 40);
            } else if (md_str_eq_cstr(ident, "FlowMode")) {
                // The prototype's numbering (groups, atoms via NTOs, orbitals): only its orbitals survive
                int mode = 0;
                if (viamd::extract_int(mode, arg)) {
                    flow.mode = (mode == 2) ? FLOW_MODE_ORBITALS : FLOW_MODE_FRAGMENTS;
                }
            } else if (md_str_eq_cstr(ident, "FlowThreshold")) {
                viamd::extract_flt(flow.cut.threshold, arg);
            } else if (md_str_eq_cstr(ident, "FlowOthersCap")) {
                viamd::extract_flt(flow.cut.other_max, arg);
            } else if (md_str_eq_cstr(ident, "FlowSizeNormalize")) {
                viamd::extract_flt(flow.size_normalize, arg);
            } else if (md_str_eq_cstr(ident, "FlowOthersExpanded")) {
                int mask = 0;
                if (viamd::extract_int(mask, arg)) flow.cut.other_expanded = (uint32_t)MD_MAX(0, mask);
            } else if (md_str_eq_cstr(ident, "FlowBarWidth")) {
                if (viamd::extract_flt(flow.bar_width, arg)) flow.params.node_thickness = flow.bar_width;
            } else if (md_str_eq_cstr(ident, "FlowNodeThickness") || md_str_eq_cstr(ident, "FlowRibbonAlpha") ||
                       md_str_eq_cstr(ident, "FlowNodeFillTint") || md_str_eq_cstr(ident, "FlowBorderWidth")) {
                // the boxed style's settings: they do not carry over to the bars
            } else if (md_str_eq_cstr(ident, "FlowNodeGap")) {
                viamd::extract_flt(flow.params.node_gap, arg);
            } else if (md_str_eq_cstr(ident, "FlowRibbonOpacity")) {
                viamd::extract_flt(flow.style.ribbon_alpha, arg);
            } else if (md_str_eq_cstr(ident, "FlowRibbonCurvature")) {
                viamd::extract_flt(flow.style.ribbon_curvature, arg);
            } else if (md_str_eq_cstr(ident, "FlowOrder")) {
                int order = MD_FLOW_ORDER_WEIGHT;
                if (viamd::extract_int(order, arg)) {
                    flow.params.order = (md_flow_order_t)MD_CLAMP(order, 0, (int)MD_FLOW_ORDER_GRAPH);
                }
            } else if (md_str_eq_cstr(ident, "FlowRibbonMinPx")) {
                viamd::extract_flt(flow.style.ribbon_min_px, arg);
            } else if (md_str_eq_cstr(ident, "FlowNodeTint")) {
                viamd::extract_flt(flow.style.node_fill_tint, arg);
            } else if (md_str_eq_cstr(ident, "FlowNodeBorder")) {
                viamd::extract_flt(flow.style.border_width, arg);
            } else if (md_str_eq_cstr(ident, "FlowColorLabels")) {
                viamd::extract_bool(flow.style.color_labels, arg);
            } else if (md_str_eq_cstr(ident, "FlowColorBySide")) {
                viamd::extract_bool(flow.style.color_by_side, arg);
            } else if (md_str_eq_cstr(ident, "FlowRibbonLabels")) {
                viamd::extract_bool(flow.style.show_ribbon_labels, arg);
            } else if (md_str_eq_cstr(ident, "FlowLegend")) {
                viamd::extract_bool(flow.show_legend, arg);
            } else if (md_str_eq_cstr(ident, "FlowSweeps")) {
                int sweeps = 2;
                if (viamd::extract_int(sweeps, arg)) {
                    flow.params.barycentre_sweeps = (uint32_t)MD_CLAMP(sweeps, 0, 8);
                }
            } else if (md_str_eq_cstr(ident, "GroupCount")) {
                int count = 0;
                if (viamd::extract_int(count, arg)) {
                    pending_groups.count = (size_t)MD_CLAMP(count, 0, MAX_NTO_GROUPS);
                }
            } else if (md_str_eq_cstr(ident, "GroupLabel")) {
                group_idx += 1;
                if (group_idx < MAX_NTO_GROUPS) {
                    viamd::extract_to_char_buf(pending_groups.label[group_idx], sizeof(pending_groups.label[0]), arg);
                }
            } else if (md_str_eq_cstr(ident, "GroupColor")) {
                if (0 <= group_idx && group_idx < MAX_NTO_GROUPS) {
                    viamd::extract_vec4(pending_groups.color[group_idx], arg);
                }
            } else if (md_str_eq_cstr(ident, "GroupAtoms")) {
                if (0 <= group_idx && group_idx < MAX_NTO_GROUPS) {
                    md_bitset_free(&pending_groups.atoms[group_idx], md_get_heap_allocator());
                    viamd::extract_bitset(&pending_groups.atoms[group_idx], arg, md_get_heap_allocator());
                }
            }
        }

        // Trust the labels actually present over a GroupCount that disagrees with them.
        pending_groups.count = MD_MAX(pending_groups.count, (size_t)(group_idx + 1));
        pending_groups.count = MD_MIN(pending_groups.count, (size_t)MAX_NTO_GROUPS);
    }

    // Consume-once: a buffered grouping belongs to the workspace that was just opened, and must
    // not be stamped onto whatever file is loaded next.
    void apply_pending_groups(const md_system_t& sys) {
        if (!pending_groups.pending) return;
        pending_groups.pending = false;

        if (!nto.atom_group_idx || pending_groups.count == 0) return;

        const size_t num_atoms = md_array_size(nto.atom_group_idx);
        md_memset(nto.atom_group_idx, 0, sizeof(uint32_t) * num_atoms);

        nto.group.count = pending_groups.count;
        for (size_t g = 0; g < pending_groups.count; ++g) {
            if (pending_groups.label[g][0]) {
                snprintf(nto.group.label[g], sizeof(nto.group.label[g]), "%s", pending_groups.label[g]);
            }
            nto.group.color[g] = pending_groups.color[g];

            for (size_t i = 0; i < num_atoms; ++i) {
                if (md_bitset_test(pending_groups.atoms[g], i)) {
                    nto.atom_group_idx[i] = (uint32_t)g;
                }
            }
        }

        update_nto_group_colors(sys);
        MD_LOG_INFO("Restored %zu atom groups from the workspace", pending_groups.count);
    }

    // ---------------------------------------------------------------------------
    // Charge Transfer Analysis
    // ---------------------------------------------------------------------------

    // Node keys are (column, level, index) so that expansion, selection and hover survive a
    // rebuild. Anything keyed by array index would reset every time the state or the grouping
    // changes, which is exactly when the user least wants to lose what they had open.
    enum flow_level_t {
        FLOW_LEVEL_GROUP = 0,
        FLOW_LEVEL_ATOM  = 1,
        FLOW_LEVEL_MO    = 2,   // a canonical molecular orbital; owns no atoms
        FLOW_LEVEL_NTO   = 3,   // an NTO pair; owns no atoms either
        FLOW_LEVEL_FRAGMENT = 4,   // a node of the fragment tree (fragment_tree.h), index into flow.frags
    };

    // Whether a node stands for a set of atoms. Orbital nodes - NTOs and MOs alike - do not, so
    // they take no part in atom selection, region sweeps or group assignment.
    static inline bool flow_level_has_atoms(uint32_t level) {
        return level == FLOW_LEVEL_GROUP || level == FLOW_LEVEL_ATOM || level == FLOW_LEVEL_FRAGMENT;
    }

    static inline uint64_t flow_node_key(uint32_t column, uint32_t level, uint32_t index) {
        return ((uint64_t)column << 40) | ((uint64_t)level << 32) | (uint64_t)index;
    }
    static inline uint32_t flow_key_index(uint64_t key) { return (uint32_t)(key & 0xFFFFFFFFu); }
    static inline uint32_t flow_key_level(uint64_t key) { return (uint32_t)((key >> 32) & 0xFFu); }
    // What a node stands for regardless of its column (md_flow_node_t.entity): never 0
    static inline uint64_t flow_entity(uint32_t level, uint32_t index) {
        return ((uint64_t)(level + 1) << 32) | (uint64_t)index;
    }

    static const size_t FLOW_ATOM_LABEL_STRIDE = 12;   // "Xx" + up to 8 digits + nul

    const char* flow_atom_label(uint32_t atom_idx) const {
        if (!flow.atom_labels || atom_idx >= flow.num_atom_labels) return "?";
        return flow.atom_labels + (size_t)atom_idx * FLOW_ATOM_LABEL_STRIDE;
    }

    // Colour says what an atom belongs to. Inside a group that is the group; outside one there is
    // no group to name, so it falls back to the atom's ROLE - donating or accepting.
    //
    // Not CPK here, deliberately. Now that colour lives on the border rather than the fill, a
    // white element (hydrogen) would draw a white outline on a near-white box and vanish, and a
    // near-black one (carbon) would say nothing about the transition. Role colour is both always
    // visible and the thing the diagram is actually about; the 3D viewport still uses CPK, because
    // there a molecule should look like a molecule.
    vec4_t flow_atom_color(uint32_t atom_idx, uint32_t column) const {
        const uint32_t group_idx = (atom_idx < md_array_size(nto.atom_group_idx)) ? nto.atom_group_idx[atom_idx] : 0;
        if (group_idx > 0 && group_idx < nto.group.count) {
            return nto.group.color[group_idx];
        }
        return (column == 0) ? vec4_t{0.78f, 0.24f, 0.22f, 1.0f}    // donating
                             : vec4_t{0.13f, 0.56f, 0.35f, 1.0f};   // accepting
    }

    char* flow_mo_label_storage(uint32_t mo_idx) {
        if (!flow.mo_labels || mo_idx >= flow.num_mo_labels) return nullptr;
        return flow.mo_labels + (size_t)mo_idx * FLOW_ATOM_LABEL_STRIDE;
    }

    void flow_build_label_storage(const md_system_t& sys) {
        const size_t num_mo = num_mos;
        if (num_mo > 0 && !flow.mo_labels) {
            flow.mo_labels = (char*)md_alloc(arena, num_mo * FLOW_ATOM_LABEL_STRIDE);
            md_memset(flow.mo_labels, 0, num_mo * FLOW_ATOM_LABEL_STRIDE);
            flow.num_mo_labels = num_mo;
        }
        flow_build_atom_labels(sys);
    }

    void flow_build_atom_labels(const md_system_t& sys) {
        QmAtoms qm = {};
        if (!es_qm_atoms(&qm, sys)) return;

        const size_t num_atoms = qm.count;
        if (num_atoms == 0 || flow.atom_labels) return;

        const md_element_t* z = qm.atomic_number;
        flow.atom_labels = (char*)md_alloc(arena, num_atoms * FLOW_ATOM_LABEL_STRIDE);
        flow.num_atom_labels = num_atoms;

        for (size_t i = 0; i < num_atoms; ++i) {
            char* dst = flow.atom_labels + i * FLOW_ATOM_LABEL_STRIDE;
            const md_str_t sym = z ? md_util_element_symbol(z[i]) : md_str_from_cstr("X");
            // One-based, matching how the rest of VIAMD numbers atoms to the user.
            snprintf(dst, FLOW_ATOM_LABEL_STRIDE, "%.*s%zu", (int)sym.len, sym.ptr, i + 1);
        }
    }

    // Fills the flow procedures' input out of the system's attribute table. Every field here is
    // either a read of a published path or a pure function of the basis - nothing reaches for
    // md_vlx - which is what lets a system published by another QM reader drive the same diagram.
    //
    // Allocates into 'temp', so that scope has to outlive the graph build which consumes it.
    // Returns false only when there is no calculation to speak of; a missing OPTIONAL part leaves
    // its field null and is the caller's to check, not an error.
    bool flow_qm_data_gather(FlowData* out, md_temp_scope_t temp, const md_system_t& sys) {
        MD_ASSERT(out);
        *out = {};

        if (!es_orbital_extent(sys, &out->num_mo, &out->num_ao)) return false;

        const size_t num_ao = out->num_ao;
        const size_t num_mo = out->num_mo;
        if (num_ao == 0 || num_mo == 0) return false;

        out->homo_idx = homo_idx[0];

        if (lowdin_job_ready(lowdin) && lowdin->num_mo == num_mo && lowdin->num_ao == num_ao) {
            out->C_lowdin = lowdin->C_lowdin;
        }

        // Extracted rather than pointed at, so a producer that publishes either of these through a
        // provider instead of storing them works here unchanged.
        if (const md_attribute_t* attr = md_attributes_find(&sys.attributes, es_path::overlap)) {
            if (attr->format.rank == 2 && attr->format.components == 1 &&
                attr->format.shape[0] == num_ao && attr->format.shape[1] == num_ao) {
                double* S = md_temp_alloc_array(temp, double, num_ao * num_ao);
                if (md_attribute_extract_f64(S, num_ao * num_ao, attr, md_attribute_slice_all(), md_unit_none()) == num_ao * num_ao) {
                    out->S = S;
                }
            }
        }

        if (const md_attribute_t* attr = md_attributes_find(&sys.attributes, es_path::alpha_coefficient)) {
            double* C = md_temp_alloc_array(temp, double, num_mo * num_ao);
            if (md_attribute_extract_f64(C, num_mo * num_ao, attr, md_attribute_slice_all(), md_unit_none()) == num_mo * num_ao) {
                out->C_mo = C;
            }
        }

        // A pure function of the shell list, and md_vlx derives its own copy by calling this same
        // procedure - so there is nothing to publish, and nothing that can drift out of step with
        // the AO ordering the evaluator walks.
        if (basis.num_shells > 0 && md_gto_basis_num_ao(&basis) == num_ao) {
            uint32_t* shell_atom = md_temp_alloc_array(temp, uint32_t, num_ao);
            if (md_gto_basis_ao_to_atom(shell_atom, &basis) == num_ao) {
                int* ao_to_atom = md_temp_alloc_array(temp, int, num_ao);
                for (size_t i = 0; i < num_ao; ++i) {
                    ao_to_atom[i] = (int)shell_atom[i];
                }
                out->ao_to_atom = ao_to_atom;
            }
        }

        // One excited state's NTO pairs. The attributes are {S,Lmax,A} and {S,Lmax}, so narrowing
        // by the state gives exactly the {K,A} block and the [K] weights the builders want, with no
        // reshaping on either side.
        if (rsp.selected >= 0) {
            const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)rsp.selected);

            const md_attribute_t* lambda_attr = md_attributes_find(&sys.attributes, es_path::nto_lambda);
            const md_attribute_t* hole_attr   = md_attributes_find(&sys.attributes, es_path::nto_hole);
            const md_attribute_t* part_attr   = md_attributes_find(&sys.attributes, es_path::nto_particle);

            if (lambda_attr && hole_attr && part_attr) {
                // Bounded by the label storage the builders write into, which is the real
                // constraint - not by any reader's idea of how many pairs a state may have.
                const size_t num_nto = MD_MIN(md_attribute_slice_count(lambda_attr, slice), MD_COUNTOF(flow.nto_label));
                const size_t block   = num_nto * num_ao;

                if (num_nto > 0 &&
                    md_attribute_slice_count(hole_attr, slice) >= block &&
                    md_attribute_slice_count(part_attr, slice) >= block) {
                    double* lambda = md_temp_alloc_array(temp, double, num_nto);
                    double* hole   = md_temp_alloc_array(temp, double, block);
                    double* part   = md_temp_alloc_array(temp, double, block);

                    if (md_attribute_extract_f64(lambda, num_nto, lambda_attr, slice, md_unit_none()) >= num_nto &&
                        md_attribute_extract_f64(hole,   block,  hole_attr,   slice, md_unit_none()) >= block &&
                        md_attribute_extract_f64(part,   block,  part_attr,   slice, md_unit_none()) >= block) {
                        out->num_nto    = num_nto;
                        out->nto_lambda = lambda;
                        out->C_nto_hole = hole;
                        out->C_nto_part = part;
                    }
                }
            }
        }

        return true;
    }

    // The automatic fragments and the QM <-> system atom maps, once per system
    void flow_build_fragments(const md_system_t& sys) {
        if (flow.frags_built && flow.frags_version == sys.topology_version) return;
        const bool again = flow.frags_built;
        flow.frags_built = true;
        flow.frags_version = sys.topology_version;
        if (again) {
            // what was made of the old tree: the strip's units
            strip_job_destroy(strip_job);
            strip_job = nullptr;
            flow.ao_unit = nullptr;
            md_array_shrink(flow.strip, 0);
        }

        QmAtoms qm = {};
        if (!es_qm_atoms(&qm, sys) || qm.count == 0) return;
        const size_t ns = md_atom_count(&sys.atom);
        flow.num_qm_atoms  = qm.count;
        flow.num_sys_atoms = ns;
        flow.qm_to_sys = (uint32_t*)md_alloc(arena, sizeof(uint32_t) * qm.count);
        flow.sys_to_qm = (int32_t*) md_alloc(arena, sizeof(int32_t) * (ns ? ns : 1));
        for (size_t i = 0; i < ns; ++i) flow.sys_to_qm[i] = -1;
        for (size_t a = 0; a < qm.count; ++a) {
            const uint32_t s = qm.system_index ? qm.system_index[a] : (uint32_t)a;
            flow.qm_to_sys[a] = s;
            if (s < ns) flow.sys_to_qm[s] = (int32_t)a;
        }
        // the molecules with an atom in the calculation
        bool* include = (bool*)md_alloc(arena, sizeof(bool) * (ns ? ns : 1));
        for (size_t i = 0; i < ns; ++i) include[i] = flow.sys_to_qm[i] >= 0;
        if (!frag_tree_build(&flow.frags, &sys, include, arena)) {
            MD_LOG_ERROR("Charge transfer: could not work out the fragments of the system");
        }
    }

    uint32_t flow_qm_to_sys(uint32_t qm_atom) const {
        return (flow.qm_to_sys && qm_atom < flow.num_qm_atoms) ? flow.qm_to_sys[qm_atom] : qm_atom;
    }

    // The leaf fragment of a QM atom, MD_FRAGMENT_NONE for none
    uint32_t flow_fragment_of(uint32_t qm_atom) const {
        const uint32_t s = flow_qm_to_sys(qm_atom);
        return (s < flow.frags.num_atoms) ? flow.frags.atom_fragment[s] : MD_FRAGMENT_NONE;
    }

    uint32_t flow_group_of(uint32_t qm_atom) const {
        return (qm_atom < md_array_size(nto.atom_group_idx)) ? nto.atom_group_idx[qm_atom] : 0u;
    }

    // A fragment's colour. Leaves (what a molecule is made of) cycle a palette; the levels above them are
    // a quiet grey-blue, so the colour on screen says which unit, not which tree level.
    vec4_t flow_fragment_color(uint32_t f) const {
        if (!frag_tree_is_leaf(&flow.frags, f)) return vec4_t{0.42f, 0.47f, 0.58f, 1.0f};
        uint32_t leaf_rank = 0;
        for (uint32_t i = 0; i < f; ++i) leaf_rank += frag_tree_is_leaf(&flow.frags, i) ? 1 : 0;
        return vec_cast(ImPlot::GetColormapColor((int)(leaf_rank % 10), ImPlotColormap_Deep));
    }

    // hole fragment -> particle fragment: the exact transfer matrix (md_qm_ct.h), from the selected state's
    // job (CtStateJob), which holds it per ATOM pair.
    //
    // The atoms are the leaves, and every coarser node is a sum of its atoms - exact, because the blocks
    // add - so collapsing and expanding cannot disagree. Above the atoms: the user's groups at the top
    // level, then the automatic fragments (md_fragment.h) of whatever is in none of them.
    //
    // Both columns list the same nodes in the same order, and a node is the same ENTITY on both sides
    // (md_flow_node_t.entity): a flow from a fragment to itself is local excitation, drawn inside the
    // node, and every band that crosses the diagram is a genuine transfer. Between atoms most of any
    // excitation crosses - neighbours in a pi system always exchange - which is why fragments open first.
    bool build_flow_graph_fragments() {
        md_flow_graph_init(&flow.graph, 2, arena);
        flow.nto_dropped = 0.0f;
        flow.link_dropped = 0.0f;

        const CtStateJob* job = ct_job;
        if (!ct_state_job_ready(job) || job->state != rsp.selected) return false;
        const size_t N = job->num_atoms;
        if (N == 0) return false;
        flow.nto_dropped  = job->nto_dropped;
        flow.link_dropped = job->link_dropped;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        const FragmentTree& F = flow.frags;
        const size_t num_frags = F.count;
        const size_t num_groups = MD_MIN((size_t)nto.group.count, (size_t)MAX_NTO_GROUPS);

        // Fragment weights per column, over the atoms in no user group
        double* frag_weight[2] = {
            md_temp_alloc_zero_array(temp, double, num_frags ? num_frags : 1),
            md_temp_alloc_zero_array(temp, double, num_frags ? num_frags : 1),
        };
        for (uint32_t a = 0; a < N; ++a) {
            if (flow_group_of(a) != 0) continue;
            const uint32_t f = flow_fragment_of(a);
            if (f == MD_FRAGMENT_NONE) continue;
            for (int c = 0; c < 2; ++c) {
                for (uint32_t g = f; g != MD_FRAGMENT_NONE; g = F.parent[g]) frag_weight[c][g] += job->weight[c][a];
            }
        }

        uint32_t* atom_node[2] = {
            md_temp_alloc_array(temp, uint32_t, N),
            md_temp_alloc_array(temp, uint32_t, N),
        };
        for (size_t a = 0; a < N; ++a) atom_node[0][a] = atom_node[1][a] = MD_FLOW_INVALID_INDEX;
        uint32_t* frag_node = md_temp_alloc_array(temp, uint32_t, num_frags ? num_frags : 1);

        for (uint32_t column = 0; column < 2; ++column) {
            const double* weight = job->weight[column];
            uint32_t* nodes = atom_node[column];

            auto add_atom = [&](uint32_t a, uint32_t parent, vec4_t color) {
                if (weight[a] <= 0.0) return;
                md_flow_node_t node = {
                    .column = column, .parent = parent, .level = FLOW_LEVEL_ATOM,
                    .weight = (float)weight[a], .color = color,
                    .label  = md_str_from_cstr(flow_atom_label(a)),
                    .key    = flow_node_key(column, FLOW_LEVEL_ATOM, a),
                    .entity = flow_entity(FLOW_LEVEL_ATOM, a),
                };
                nodes[a] = md_flow_graph_add_node(&flow.graph, &node);
            };

            // The user's groups, at the top
            for (uint32_t g = 1; g < num_groups; ++g) {
                double group_weight = 0.0;
                for (uint32_t a = 0; a < N; ++a) if (flow_group_of(a) == g) group_weight += weight[a];
                if (group_weight <= 0.0) continue;
                md_flow_node_t group_node = {
                    .column = column, .parent = MD_FLOW_INVALID_INDEX, .level = FLOW_LEVEL_GROUP,
                    .weight = (float)group_weight, .color = nto.group.color[g],
                    .label  = md_str_from_cstr(nto.group.label[g]),
                    .key    = flow_node_key(column, FLOW_LEVEL_GROUP, g),
                    .entity = flow_entity(FLOW_LEVEL_GROUP, g),
                };
                const uint32_t parent = md_flow_graph_add_node(&flow.graph, &group_node);
                if (parent == MD_FLOW_INVALID_INDEX) continue;
                for (uint32_t a = 0; a < N; ++a) {
                    if (flow_group_of(a) == g) add_atom(a, parent, nto.group.color[g]);
                }
            }

            // The automatic fragments of the rest, in their own (pre-)order
            for (size_t f = 0; f < num_frags; ++f) {
                frag_node[f] = MD_FLOW_INVALID_INDEX;
                if (frag_weight[column][f] <= 0.0) continue;
                uint32_t parent = MD_FLOW_INVALID_INDEX;
                for (uint32_t p = F.parent[f]; p != MD_FRAGMENT_NONE; p = F.parent[p]) {
                    if (frag_node[p] != MD_FLOW_INVALID_INDEX) { parent = frag_node[p]; break; }
                }
                const vec4_t color = flow_fragment_color((uint32_t)f);
                md_flow_node_t node = {
                    .column = column, .parent = parent, .level = FLOW_LEVEL_FRAGMENT,
                    .weight = (float)frag_weight[column][f], .color = color,
                    .label  = frag_tree_name(&F, f),
                    .key    = flow_node_key(column, FLOW_LEVEL_FRAGMENT, (uint32_t)f),
                    .entity = flow_entity(FLOW_LEVEL_FRAGMENT, (uint32_t)f),
                };
                frag_node[f] = md_flow_graph_add_node(&flow.graph, &node);
                if (frag_node[f] == MD_FLOW_INVALID_INDEX || !frag_tree_is_leaf(&F, f)) continue;
                for (uint32_t k = F.atom_beg[f]; k < F.atom_end[f]; ++k) {
                    const int32_t s = F.atom_idx[k];
                    const int32_t a = (flow.sys_to_qm && s >= 0 && (size_t)s < flow.num_sys_atoms) ? flow.sys_to_qm[s] : s;
                    if (a < 0 || (size_t)a >= N || flow_group_of((uint32_t)a) != 0) continue;
                    add_atom((uint32_t)a, frag_node[f], color);
                }
            }

            // QM atoms in no fragment (none, when the system has its atoms): roots of their own
            for (uint32_t a = 0; a < N; ++a) {
                if (flow_group_of(a) != 0 || flow_fragment_of(a) != MD_FRAGMENT_NONE) continue;
                add_atom(a, MD_FLOW_INVALID_INDEX, flow_atom_color(a, column));
            }
        }

        for (size_t a = 0; a < N; ++a) {
            if (atom_node[0][a] == MD_FLOW_INVALID_INDEX) continue;
            const double* row = job->omega + a * N;
            for (size_t b = 0; b < N; ++b) {
                if (atom_node[1][b] == MD_FLOW_INVALID_INDEX) continue;
                const float w = (float)row[b];
                if (w > 0.0f) md_flow_graph_add_link(&flow.graph, atom_node[0][a], atom_node[1][b], w);
            }
        }

        return true;
    }

    // "HOMO-2", "LUMO+3" - how a chemist names an orbital, which is the entire point of showing
    // this column. Falls back to a bare index when the HOMO is unknown.
    void flow_mo_label(char* buf, size_t cap, size_t mo_idx, size_t homo_idx) {
        if (mo_idx <= homo_idx) {
            const size_t below = homo_idx - mo_idx;
            if (below == 0) snprintf(buf, cap, "HOMO");
            else            snprintf(buf, cap, "HOMO-%zu", below);
        } else {
            const size_t above = mo_idx - homo_idx - 1;
            if (above == 0) snprintf(buf, cap, "LUMO");
            else            snprintf(buf, cap, "LUMO+%zu", above);
        }
    }

    // The MO pair an NTO is mostly made of, as a label like "HOMO-1 -> LUMO".
    //
    // This is what the sketch's "pi* / sigma* / Rydberg" is reaching for, but those are symmetry
    // assignments and we have no classifier - inventing one would put a confident chemical claim on
    // screen that nothing computed. The dominant MO pair is the honest version: it is measured, it
    // is what people write in papers, and it says something the NTO index alone never can.
    //
    // Left blank unless one MO leads on each side by a clear margin. A 40/35 split has no dominant
    // character and saying otherwise would be worse than saying nothing.
    void flow_nto_character(char* buf, size_t cap, const FlowData& qm, const double* C_hole_k, const double* C_part_k, float min_share) {
        buf[0] = '\0';

        const size_t num_ao  = qm.num_ao;
        const size_t num_mo  = qm.num_mo;
        const double* S      = qm.S;
        const double* C_all  = qm.C_mo;
        if (!S || !C_all || qm.homo_idx < 0) return;

        const size_t homo = (size_t)qm.homo_idx;
        if (num_ao == 0 || num_mo == 0 || homo + 1 >= num_mo) return;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        double* sc = md_temp_alloc_array(temp, double, num_ao);

        size_t best[2] = { num_mo, num_mo };
        double share[2] = { 0.0, 0.0 };

        for (int side = 0; side < 2; ++side) {
            const double* C_nto = (side == 0) ? C_hole_k : C_part_k;
            const size_t beg = (side == 0) ? 0    : homo + 1;
            const size_t end = (side == 0) ? homo + 1 : num_mo;

            for (size_t mu = 0; mu < num_ao; ++mu) {
                double sum = 0.0;
                for (size_t nu = 0; nu < num_ao; ++nu) sum += S[mu * num_ao + nu] * C_nto[nu];
                sc[mu] = sum;
            }

            double norm = 0.0, top = 0.0;
            size_t top_mo = num_mo;
            for (size_t mo = beg; mo < end; ++mo) {
                const double* C_mo = C_all + mo * num_ao;
                double ov = 0.0;
                for (size_t mu = 0; mu < num_ao; ++mu) ov += C_mo[mu] * sc[mu];
                const double w = ov * ov;
                norm += w;
                if (w > top) { top = w; top_mo = mo; }
            }
            if (norm <= 1.0e-9 || top_mo >= num_mo) return;
            best[side] = top_mo;
            share[side] = top / norm;
        }

        if (share[0] < (double)min_share || share[1] < (double)min_share) return;

        char lo[24], hi[24];
        flow_mo_label(lo, sizeof(lo), best[0], homo);
        flow_mo_label(hi, sizeof(hi), best[1], homo);
        snprintf(buf, cap, "%s " "\xe2\x86\x92" " %s", lo, hi);
    }

    // occupied MO -> NTO pair -> virtual MO.
    //
    // NTOs ARE a rotation of the MOs: hole NTO k = sum_i U_ik phi_i over occupied MOs, particle
    // NTO k = sum_a V_ak phi_a over virtual ones. U and V are orthogonal, so |U_ik|^2 is an exact
    // normalized share and every band here is exact - no Mulliken partitioning anywhere.
    //
    // This is deliberately a route WITHOUT atoms rather than an extra column inside the atom
    // route. Chaining atom <- MO <- NTO would keep only the diagonal of
    //     q_k(A) = sum_ij U_ik U_jk P_ij(A)
    // and drop the MO-MO interference. Measured on test_data/vlx/amide.h5, that moves 33.6% of the
    // population on the dominant transition: a diagram that adds up perfectly and is wrong by a
    // third about which atoms donated. See docs/transition_flow_design.md.
    bool build_flow_graph_orbitals(const FlowData& qm) {
        md_flow_graph_init(&flow.graph, 3, arena);

        if (rsp.selected < 0) return false;

        const size_t num_ao  = qm.num_ao;
        const size_t num_mo  = qm.num_mo;
        const double* S      = qm.S;
        const double* C_all  = qm.C_mo;

        if (num_ao == 0 || num_mo == 0 || !S || !C_all || qm.homo_idx < 0) return false;
        const size_t homo = (size_t)qm.homo_idx;
        const size_t nocc = homo + 1;
        if (nocc == 0 || nocc >= num_mo) return false;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        const double* lambdas = qm.nto_lambda;
        const double* C_hole  = qm.C_nto_hole;
        const double* C_part  = qm.C_nto_part;
        const size_t num_nto  = qm.num_nto;
        if (num_nto == 0 || !lambdas || !C_hole || !C_part) return false;

        // Projection onto the MO basis: U_ik = <phi_i | NTO_k> = C_i^T S C_k. S*C_k is formed once
        // per NTO and then dotted with each MO, which keeps this at O(k * nao^2) rather than
        // O(k * nmo * nao^2).
        double* sc = md_temp_alloc_array(temp, double, num_ao);
        double* w_occ = md_temp_alloc_zero_array(temp, double, num_mo);
        double* w_vir = md_temp_alloc_zero_array(temp, double, num_mo);
        double* proj_occ = md_temp_alloc_zero_array(temp, double, num_nto * num_mo);
        double* proj_vir = md_temp_alloc_zero_array(temp, double, num_nto * num_mo);

        double* lambda_keep = md_temp_alloc_zero_array(temp, double, num_nto);
        double lambda_sum = 0.0;
        size_t num_keep = 0;

        for (size_t k = 0; k < num_nto; ++k) {
            if (lambdas[k] < (double)flow.lambda_cutoff) continue;
            if (num_keep >= MAX_NTO_LABELS) break;  // the labels' room; the rest weigh least (descending)

            for (int side = 0; side < 2; ++side) {
                const double* C_nto = (side == 0 ? C_hole : C_part) + k * num_ao;
                double* proj = (side == 0 ? proj_occ : proj_vir) + num_keep * num_mo;
                const size_t mo_beg = (side == 0) ? 0    : nocc;
                const size_t mo_end = (side == 0) ? nocc : num_mo;

                for (size_t mu = 0; mu < num_ao; ++mu) {
                    double sum = 0.0;
                    for (size_t nu = 0; nu < num_ao; ++nu) sum += S[mu * num_ao + nu] * C_nto[nu];
                    sc[mu] = sum;
                }

                double norm = 0.0;
                for (size_t mo = mo_beg; mo < mo_end; ++mo) {
                    const double* C_mo = C_all + mo * num_ao;
                    double overlap = 0.0;
                    for (size_t mu = 0; mu < num_ao; ++mu) overlap += C_mo[mu] * sc[mu];
                    proj[mo] = overlap * overlap;   // |U_ik|^2
                    norm += proj[mo];
                }
                // Should already be 1 by orthogonality; renormalizing absorbs the numerical drift
                // rather than letting it leak into the conservation check.
                if (norm > 1.0e-9) {
                    for (size_t mo = mo_beg; mo < mo_end; ++mo) proj[mo] /= norm;
                } else {
                    continue;
                }
            }

            snprintf(flow.nto_label[num_keep], sizeof(flow.nto_label[0]), "NTO %zu", k + 1);
            flow.nto_index[num_keep] = (uint32_t)k;
            flow_nto_character(flow.nto_character[num_keep], sizeof(flow.nto_character[0]), qm,
                               C_hole + k * num_ao, C_part + k * num_ao, 0.35f);

            lambda_keep[num_keep] = lambdas[k];
            lambda_sum += lambdas[k];
            num_keep++;
        }

        if (num_keep == 0 || lambda_sum <= 0.0) return false;
        for (size_t k = 0; k < num_keep; ++k) lambda_keep[k] /= lambda_sum;

        for (size_t k = 0; k < num_keep; ++k) {
            for (size_t mo = 0; mo < num_mo; ++mo) {
                w_occ[mo] += lambda_keep[k] * proj_occ[k * num_mo + mo];
                w_vir[mo] += lambda_keep[k] * proj_vir[k * num_mo + mo];
            }
        }

        uint32_t* occ_node = md_temp_alloc_array(temp, uint32_t, num_mo);
        uint32_t* vir_node = md_temp_alloc_array(temp, uint32_t, num_mo);
        uint32_t* mid_node = md_temp_alloc_array(temp, uint32_t, num_keep);
        for (size_t mo = 0; mo < num_mo; ++mo) { occ_node[mo] = vir_node[mo] = MD_FLOW_INVALID_INDEX; }

        // Every orbital with any weight becomes a node. There is no truncation here on purpose:
        // dropping the tail would break conservation, and the cut's threshold already folds the
        // small ones into "Others" at view time, where it can be undone.
        char label[32];
        for (size_t mo = 0; mo < num_mo; ++mo) {
            const bool occupied = mo <= homo;
            const double w = occupied ? w_occ[mo] : w_vir[mo];
            if (w <= 0.0) continue;

            flow_mo_label(label, sizeof(label), mo, homo);
            char* dst = flow_mo_label_storage((uint32_t)mo);
            if (!dst) continue;
            snprintf(dst, FLOW_ATOM_LABEL_STRIDE, "%s", label);

            md_flow_node_t node = {
                .column = occupied ? 0u : 2u,
                .parent = MD_FLOW_INVALID_INDEX,
                .level  = FLOW_LEVEL_MO,
                .weight = (float)w,
                .color  = occupied ? vec4_t{0.78f, 0.24f, 0.22f, 1.0f} : vec4_t{0.13f, 0.56f, 0.35f, 1.0f},
                .label  = md_str_from_cstr(dst),
                .key    = flow_node_key(occupied ? 0u : 2u, FLOW_LEVEL_MO, (uint32_t)mo),
            };
            (occupied ? occ_node : vir_node)[mo] = md_flow_graph_add_node(&flow.graph, &node);
        }

        for (size_t k = 0; k < num_keep; ++k) {
            const ImVec4 c = ImPlot::GetColormapColor((int)k, ImPlotColormap_Deep);
            md_flow_node_t node = {
                .column = 1, .parent = MD_FLOW_INVALID_INDEX, .level = FLOW_LEVEL_NTO,
                .weight = (float)lambda_keep[k], .color = vec_cast(c),
                .label = md_str_from_cstr(flow.nto_label[k]),
                .sublabel = md_str_from_cstr(flow.nto_character[k]),
                .key = flow_node_key(1, FLOW_LEVEL_NTO, (uint32_t)k),
            };
            mid_node[k] = md_flow_graph_add_node(&flow.graph, &node);
        }

        for (size_t k = 0; k < num_keep; ++k) {
            for (size_t mo = 0; mo < num_mo; ++mo) {
                if (occ_node[mo] != MD_FLOW_INVALID_INDEX) {
                    const float w = (float)(lambda_keep[k] * proj_occ[k * num_mo + mo]);
                    if (w > 0.0f) md_flow_graph_add_link(&flow.graph, occ_node[mo], mid_node[k], w);
                }
                if (vir_node[mo] != MD_FLOW_INVALID_INDEX) {
                    const float w = (float)(lambda_keep[k] * proj_vir[k * num_mo + mo]);
                    if (w > 0.0f) md_flow_graph_add_link(&flow.graph, mid_node[k], vir_node[mo], w);
                }
            }
        }

        return true;
    }

    static int flow_key_cmp(const void* a, const void* b) {
        const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
        return x < y ? -1 : x > y ? 1 : 0;
    }
    static bool flow_key_in(const md_array(uint64_t) sorted, uint64_t key) {
        return sorted && bsearch(&key, sorted, md_array_size(sorted), sizeof(uint64_t), flow_key_cmp) != nullptr;
    }

    // qm: what the Orbitals route is built from; the Fragments route reads the state's job instead
    void rebuild_flow_graph(const FlowData* qm) {
        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        md_allocator_i* talloc = md_temp_allocator(temp);

        // The cut stores what the user opened, closed and unfolded by node INDEX, and a rebuild renumbers
        // every node - a new excited state or a regroup would otherwise leave the pins on whatever now
        // happens to sit at those indices. Carried across by KEY, which is what the key field exists for.
        md_array(uint64_t) pinned_open = 0;
        md_array(uint64_t) pinned_closed = 0;
        md_array(uint64_t) others_open = 0;
        md_array(uint64_t) known_keys = 0;
        for (size_t i = 0; i < md_array_size(flow.graph.nodes); ++i) {
            const uint64_t key = flow.graph.nodes[i].key;
            md_array_push(known_keys, key, talloc);
            if (md_bitmask_test_bit(&flow.cut.expanded,   (uint32_t)i)) md_array_push(pinned_open,   key, talloc);
            if (md_bitmask_test_bit(&flow.cut.collapsed,  (uint32_t)i)) md_array_push(pinned_closed, key, talloc);
            if (md_bitmask_test_bit(&flow.cut.other_open, (uint32_t)i)) md_array_push(others_open,   key, talloc);
        }
        md_array(uint64_t) lists[4] = { pinned_open, pinned_closed, others_open, known_keys };
        for (int l = 0; l < 4; ++l) {
            if (lists[l]) qsort(lists[l], md_array_size(lists[l]), sizeof(uint64_t), flow_key_cmp);
        }

        md_flow_graph_free(&flow.graph);
        flow.built = false;
        flow.hover = FlowHit{};     // a cut index, and stale now (selection is by key, so it lives)

        const bool ok = (flow.mode == FLOW_MODE_ORBITALS) ? (qm && build_flow_graph_orbitals(*qm))
                                                          : build_flow_graph_fragments();
        if (!ok) {
            md_flow_graph_free(&flow.graph);
            return;
        }

        // Validating every build is cheap and it is the only thing standing between a arithmetic
        // slip upstream and a diagram that looks plausible while being wrong.
        const md_flow_result_t result = md_flow_graph_validate(&flow.graph, 1.0e-3f);
        if (result != MD_FLOW_OK) {
            MD_LOG_ERROR("Charge transfer graph is inconsistent: %s", md_flow_result_str(result));
            md_flow_graph_clear(&flow.graph);
            return;
        }

        // Fragments: what opens is the row budget's call unless the user pinned it, and the pins live.
        // Orbitals has no hierarchy; everything new opens.
        const bool fragments = flow.mode == FLOW_MODE_FRAGMENTS;
        flow.cut.mirror   = fragments;
        flow.cut.max_rows = fragments ? (uint32_t)MD_MAX(flow.max_rows, 0) : 0u;
        md_flow_cut_clear_pins(&flow.cut);
        md_bitmask_clear(&flow.cut.other_open);
        for (size_t i = 0; i < md_array_size(flow.graph.nodes); ++i) {
            const uint64_t key = flow.graph.nodes[i].key;
            if (flow_key_in(lists[0], key) || (!fragments && !flow_key_in(lists[3], key))) md_flow_cut_set_expanded(&flow.cut, (uint32_t)i, true);
            if (flow_key_in(lists[1], key)) md_bitmask_set_bit(&flow.cut.collapsed,  (uint32_t)i);
            if (flow_key_in(lists[2], key)) md_bitmask_set_bit(&flow.cut.other_open, (uint32_t)i);
        }

        flow.built = true;
        flow.fit_requested = true;
        md_flow_cut_resolve(&flow.cut, &flow.graph);
    }

    // Frame the diagram with room to breathe.
    //
    // Resetting to pan=0, zoom=1 makes the unit square exactly fill the canvas, which sounds right
    // and looks wrong: the layout's unit square bounds the NODES, while bands hang off the sides
    // and captions sit above the top row. Those are sized in pixels by the drawing backend, so the
    // fit has to reserve pixels for them and only then convert to a zoom.
    // What the diagram's columns mean, for the backend to draw: captions, the state with its split into
    // charge transfer and local excitation at the level open now, and a legend.
    void flow_annotations(FlowAnnotations* notes, char* subtitle, size_t cap) const {
        *notes = {};
        if (flow.mode == FLOW_MODE_ORBITALS) {
            notes->caption[0] = "Occupied";
            notes->caption[1] = "NTO pairs";
            notes->caption[2] = "Virtual";
            snprintf(subtitle, cap, "S%d", rsp.selected + 1);
            notes->subtitle = subtitle;
            notes->legend_ribbon = "Ribbon = hole to particle through an orbital pair";
            return;
        }
        notes->caption[0] = "Hole";
        notes->caption[2] = "Particle";
        float local = 0.0f, total = 0.0f;
        for (size_t i = 0; i < md_array_size(flow.cut.links); ++i) {
            const md_flow_link_t* l = flow.cut.links + i;
            const md_flow_node_t* a = md_flow_cut_node(&flow.cut, &flow.graph, l->src);
            const md_flow_node_t* b = md_flow_cut_node(&flow.cut, &flow.graph, l->dst);
            if (a && b && a->entity != 0 && a->entity == b->entity) local += l->weight;
            total += l->weight;
        }
        if (total > 0.0f) {
            snprintf(subtitle, cap, "S%d \xC2\xB7 CT %.0f%% \xC2\xB7 LE %.0f%%", rsp.selected + 1,
                     100.0f * (1.0f - local / total), 100.0f * local / total);
        } else {
            snprintf(subtitle, cap, "S%d", rsp.selected + 1);
        }
        notes->subtitle = subtitle;
        if (flow.show_legend) {
            notes->legend_ribbon = "Ribbon = transfer between fragments";
            notes->legend_local  = "Solid = local excitation";
        }
    }

    // Pan and zoom that put the nodes' box, with the room its labels and annotations take, inside a
    // canvas of 'canvas_size' pixels.
    void flow_fit(vec2_t* pan, float* zoom, const md_flow_layout_t& layout, ImVec2 canvas_size, const FlowMargins& m) const {
        *pan = { 0, 0 };
        *zoom = 1.0f;
        if (md_array_size(layout.nodes) == 0) return;

        vec2_t bb_min = { 1.0e30f,  1.0e30f };
        vec2_t bb_max = {-1.0e30f, -1.0e30f };
        for (size_t i = 0; i < md_array_size(layout.nodes); ++i) {
            const md_flow_layout_node_t* n = layout.nodes + i;
            bb_min.x = MD_MIN(bb_min.x, n->min.x);  bb_min.y = MD_MIN(bb_min.y, n->min.y);
            bb_max.x = MD_MAX(bb_max.x, n->max.x);  bb_max.y = MD_MAX(bb_max.y, n->max.y);
        }
        const float bb_w = MD_MAX(bb_max.x - bb_min.x, 1.0e-4f);
        const float bb_h = MD_MAX(bb_max.y - bb_min.y, 1.0e-4f);

        // a plain margin as well, so the diagram is not welded to the window edge
        const float margin = 12.0f;
        const float pad_l = margin + m.left,  pad_r = margin + m.right;
        const float pad_t = margin + m.top,   pad_b = margin + m.bottom;
        const float avail_w = MD_MAX(canvas_size.x - pad_l - pad_r, 16.0f);
        const float avail_h = MD_MAX(canvas_size.y - pad_t - pad_b, 16.0f);

        // Graph units are fractions of the canvas, so a graph-space extent of bb_w occupies
        // bb_w * zoom * canvas_size.x pixels.
        *zoom = MD_CLAMP(MD_MIN(avail_w / (bb_w * canvas_size.x), avail_h / (bb_h * canvas_size.y)), 0.1f, 12.0f);

        // Centred in what the pads leave
        const vec2_t bb_center = { (bb_min.x + bb_max.x) * 0.5f, (bb_min.y + bb_max.y) * 0.5f };
        const float target_x = 0.5f + (pad_l - pad_r) * 0.5f / canvas_size.x;
        const float target_y = 0.5f + (pad_t - pad_b) * 0.5f / canvas_size.y;
        pan->x = target_x - bb_center.x * *zoom;
        pan->y = target_y - bb_center.y * *zoom;
    }

    void flow_fit_view(ImVec2 canvas_size) {
        if (!flow.built) {
            flow.pan = {0, 0};
            flow.zoom = 1.0f;
            return;
        }
        FlowAnnotations notes = {};
        char subtitle[96];
        flow_annotations(&notes, subtitle, sizeof(subtitle));
        const FlowMargins m = flow_measure_margins(&flow.graph, &flow.cut, &flow.layout, flow.style, &notes);
        flow_fit(&flow.pan, &flow.zoom, flow.layout, canvas_size, m);
    }

    // Visible nodes whose on-screen rect overlaps a region, in cut index space.
    //
    // Orbital nodes are skipped: they own no atoms, so sweeping over one and calling that
    // "selecting atoms" would either do nothing or, worse, quietly select the whole molecule.
    void flow_nodes_in_region(md_array(uint32_t)* out_nodes, const FlowView& view, ImVec2 rect_min, ImVec2 rect_max,
                              md_allocator_i* alloc)
    {
        for (size_t i = 0; i < md_array_size(flow.layout.nodes); ++i) {
            const md_flow_layout_node_t* ln = flow.layout.nodes + i;
            const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, ln->cut_idx);
            if (!node || !flow_level_has_atoms(node->level)) continue;

            const ImVec2 min = flow_view_to_screen(view, ln->min);
            const ImVec2 max = flow_view_to_screen(view, ln->max);

            // Overlap, not containment: a sweep across a column should catch the rows it crosses,
            // not only the ones it happens to swallow whole.
            if (max.x < rect_min.x || min.x > rect_max.x) continue;
            if (max.y < rect_min.y || min.y > rect_max.y) continue;

            md_array_push(*out_nodes, ln->cut_idx, alloc);
        }
    }

    // Group membership from whatever is currently selected - wherever that selection was made.
    void flow_assign_selection_to_group(const ApplicationState& state, uint32_t group_idx) {
        if (group_idx >= nto.group.count) return;
        for (size_t i = 0; i < md_array_size(nto.atom_group_idx); ++i) {
            if (md_bitmask_test_bit(&state.selection.selection_mask, flow_qm_to_sys((uint32_t)i))) {
                nto.atom_group_idx[i] = group_idx;
            }
        }
    }

    // Returns the new group index, or 0 (the ungrouped bucket) when there is no room.
    uint32_t flow_create_group_from_selection(const ApplicationState& state) {
        if (nto.group.count >= MAX_NTO_GROUPS) {
            MD_LOG_ERROR("Cannot create another group: the limit of %d is reached", (int)MAX_NTO_GROUPS);
            return 0;
        }
        const uint32_t group_idx = (uint32_t)nto.group.count++;
        flow_assign_selection_to_group(state, group_idx);
        return group_idx;
    }

    bool flow_orbital_is_selected(uint64_t key) const {
        for (size_t i = 0; i < md_array_size(flow.selected_orbital_keys); ++i) {
            if (flow.selected_orbital_keys[i] == key) return true;
        }
        return false;
    }

    void flow_orbital_set_selected(uint64_t key, bool selected) {
        for (size_t i = 0; i < md_array_size(flow.selected_orbital_keys); ++i) {
            if (flow.selected_orbital_keys[i] != key) continue;
            if (!selected) {
                flow.selected_orbital_keys[i] = flow.selected_orbital_keys[md_array_size(flow.selected_orbital_keys) - 1];
                md_array_shrink(flow.selected_orbital_keys, md_array_size(flow.selected_orbital_keys) - 1);
            }
            return;
        }
        if (selected) md_array_push(flow.selected_orbital_keys, key, arena);
    }

    // One node's contribution to a selection change. Atom-bearing nodes edit the application mask
    // directly - there is no diagram-local copy to keep in step - while NTO nodes are a separate
    // axis and never touch it.
    void flow_set_selected(ApplicationState& state, uint32_t cut_idx, bool selected) {
        const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, cut_idx);
        if (!node) return;

        if (!flow_level_has_atoms(node->level)) {
            flow_orbital_set_selected(node->key, selected);
            return;
        }

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        const md_bitset_t mask = flow_node_atom_mask(cut_idx, md_temp_allocator(temp));
        if (selected) {
            md_bitmask_or(&state.selection.selection_mask, mask);
        } else {
            md_bitmask_andnot(&state.selection.selection_mask, mask);
        }
    }

    void flow_clear_selection(ApplicationState& state) {
        md_bitmask_clear(&state.selection.selection_mask);
        md_array_shrink(flow.selected_orbital_keys, 0);
    }

    // Atoms belonging to a visible node, for the two-way link with the viewport.
    // Allocated with alloc
    md_bitset_t flow_node_atom_mask(uint32_t cut_idx, md_allocator_i* alloc) {
        const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, cut_idx);
        if (!node || (node->flags & MD_FLOW_NODE_FLAG_OTHER)) return md_bitset_empty_set();

        // Orbitals own no atoms; selecting "all the atoms an NTO touches" would be the whole
        // molecule and mean nothing. Gated on the LEVEL rather than the column, because in the
        // Orbitals route the outer columns are orbitals too.
        if (!flow_level_has_atoms(node->level)) return md_bitset_empty_set();

        // Atoms are QM atoms here and system atoms in the mask
        const uint32_t index = flow_key_index(node->key);
        if (node->level == FLOW_LEVEL_ATOM) {
            const uint32_t s = flow_qm_to_sys(index);
            return md_bitset_range(s, s + 1);
        }
        md_temp_scope_t temp = md_temp_begin_avoid(alloc);
        md_bitmask_t b = md_bitmask_create(0, MD_BITSET_MAX_INDEX, md_temp_allocator(temp));
        for (size_t a = 0; a < md_array_size(nto.atom_group_idx); ++a) {
            bool in = false;
            if (node->level == FLOW_LEVEL_GROUP) {
                in = nto.atom_group_idx[a] == index;
            } else if (node->level == FLOW_LEVEL_FRAGMENT && nto.atom_group_idx[a] == 0) {
                for (uint32_t f = flow_fragment_of((uint32_t)a); f != MD_FRAGMENT_NONE && !in; f = flow.frags.parent[f]) in = f == index;
            }
            if (in) md_bitmask_set_bit(&b, flow_qm_to_sys((uint32_t)a));
        }
        const md_bitset_t result = md_bitmask_extract(&b, alloc);
        md_temp_end(temp);
        return result;
    }

    // The part of a visible node's weight that flows to the same entity on the other side: local
    // excitation, in the fragment route.
    float flow_node_self_weight(uint32_t cut_idx) const {
        float sum = 0.0f;
        for (size_t i = 0; i < md_array_size(flow.cut.links); ++i) {
            const md_flow_link_t* l = flow.cut.links + i;
            if (l->src != cut_idx && l->dst != cut_idx) continue;
            const md_flow_node_t* a = md_flow_cut_node(&flow.cut, &flow.graph, l->src);
            const md_flow_node_t* b = md_flow_cut_node(&flow.cut, &flow.graph, l->dst);
            if (a && b && a->entity != 0 && a->entity == b->entity) sum += l->weight;
        }
        return sum;
    }

    // The QM atoms a visible node stands for; false for an "Others" and for orbitals
    bool flow_node_qm_atoms(md_array(uint32_t)* out, uint32_t cut_idx, md_allocator_i* alloc) const {
        const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, cut_idx);
        if (!node || (node->flags & MD_FLOW_NODE_FLAG_OTHER) || !flow_level_has_atoms(node->level)) return false;
        const uint32_t index = flow_key_index(node->key);
        const size_t N = md_array_size(nto.atom_group_idx);
        if (node->level == FLOW_LEVEL_ATOM) {
            md_array_push(*out, index, alloc);
        } else if (node->level == FLOW_LEVEL_GROUP) {
            for (uint32_t a = 0; a < N; ++a) if (nto.atom_group_idx[a] == index) md_array_push(*out, a, alloc);
        } else if (index < flow.frags.count) {
            for (uint32_t k = flow.frags.atom_beg[index]; k < flow.frags.atom_end[index]; ++k) {
                const int32_t s = flow.frags.atom_idx[k];
                const int32_t a = (flow.sys_to_qm && s >= 0 && (size_t)s < flow.num_sys_atoms) ? flow.sys_to_qm[s] : s;
                if (a >= 0 && (size_t)a < N && nto.atom_group_idx[a] == 0) md_array_push(*out, (uint32_t)a, alloc);
            }
        }
        return true;
    }

    // What each visible ribbon is made of: the singular value decomposition of its block of the Lowdin
    // transition density (md_qm_ct_block_nto), from the K x K blocks of its two ends, which are sums of
    // their atoms' blocks. Exact, never negative, and adding up to the ribbon. Recomputed when the
    // visible links change.
    void flow_update_stripes() {
        const bool active = flow.mode == FLOW_MODE_FRAGMENTS && flow.built && ct_state_job_ready(ct_job) && ct_job->state == rsp.selected;
        const size_t L = md_array_size(flow.cut.links);
        const uint64_t hash = active ? md_hash64(flow.cut.links, sizeof(md_flow_link_t) * L, flow.data_hash ^ 0x57121e5ull) : 0;
        if (hash == flow.stripe_hash && (active || md_array_size(flow.stripe_offset) == 0)) return;
        flow.stripe_hash = hash;
        md_array_shrink(flow.stripe_offset, 0);
        md_array_shrink(flow.stripe_frac, 0);
        if (!active) return;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        md_allocator_i* talloc = md_temp_allocator(temp);

        const size_t K = ct_job->K, KK = K * K;
        const size_t cut_space = md_array_size(flow.graph.nodes) + md_array_size(flow.cut.other);
        double** block = md_temp_alloc_zero_array(temp, double*, cut_space ? cut_space : 1);
        auto node_block = [&](uint32_t cut_idx, const double* atom_blocks) -> const double* {
            if (cut_idx >= cut_space) return nullptr;
            if (block[cut_idx]) return block[cut_idx];
            md_array(uint32_t) atoms = 0;
            if (!flow_node_qm_atoms(&atoms, cut_idx, talloc) || md_array_size(atoms) == 0) return nullptr;
            double* b = md_temp_alloc_zero_array(temp, double, KK);
            for (size_t i = 0; i < md_array_size(atoms); ++i) {
                if (atoms[i] >= ct_job->num_atoms) continue;
                const double* ab = atom_blocks + (size_t)atoms[i] * KK;
                for (size_t k = 0; k < KK; ++k) b[k] += ab[k];
            }
            return block[cut_idx] = b;
        };

        double* w = md_temp_alloc_array(temp, double, K);
        for (size_t l = 0; l < L; ++l) {
            md_array_push(flow.stripe_offset, (uint32_t)md_array_size(flow.stripe_frac), arena);
            const md_flow_link_t* link = flow.cut.links + l;
            const md_flow_node_t* a = md_flow_cut_node(&flow.cut, &flow.graph, link->src);
            const md_flow_node_t* b = md_flow_cut_node(&flow.cut, &flow.graph, link->dst);
            if (!a || !b || (a->entity != 0 && a->entity == b->entity)) continue;    // local: drawn in the node
            const double* H = node_block(link->src, ct_job->hole_blocks);
            const double* P = node_block(link->dst, ct_job->part_blocks);
            if (!H || !P) continue;
            const size_t n = md_qm_ct_block_nto(w, nullptr, nullptr, ct_job->lambda, K, H, P);
            double sum = 0.0;
            for (size_t j = 0; j < n; ++j) sum += w[j];
            if (n < 2 || sum <= 0.0) continue;
            // the pairs that matter, and the rest as one
            double cum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                const double f = w[j] / sum;
                if (f < 0.01 || cum > 0.995) {
                    md_array_push(flow.stripe_frac, (float)(1.0 - cum), arena);
                    break;
                }
                md_array_push(flow.stripe_frac, (float)f, arena);
                cum += f;
            }
        }
        md_array_push(flow.stripe_offset, (uint32_t)md_array_size(flow.stripe_frac), arena);
    }

    // The strip's inputs, once the fragments and the Lowdin basis are there: the unit of every AO (the leaf
    // fragments, numbered compactly) and the QM atoms' positions
    void flow_build_strip_inputs(const md_system_t& sys) {
        if (flow.ao_unit || !lowdin_job_ready(lowdin)) return;
        const size_t A = lowdin->num_ao, N = lowdin->num_atoms;
        const FragmentTree& F = flow.frags;

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        uint32_t* leaf_id = md_temp_alloc_array(temp, uint32_t, F.count ? F.count : 1);
        uint32_t U = 0;
        for (size_t f = 0; f < F.count; ++f) leaf_id[f] = frag_tree_is_leaf(&F, f) ? U++ : MD_FRAGMENT_NONE;

        flow.ao_unit = (uint32_t*)md_alloc(arena, sizeof(uint32_t) * A);
        for (size_t mu = 0; mu < A; ++mu) {
            const uint32_t atom = lowdin->ao_atom[mu];
            const uint32_t f = atom < N ? flow_fragment_of(atom) : MD_FRAGMENT_NONE;
            flow.ao_unit[mu] = (f != MD_FRAGMENT_NONE) ? leaf_id[f] : U;    // U and above: in no unit, left out
        }
        flow.num_units = U;

        QmAtoms qm = {};
        if (es_qm_atoms(&qm, sys) && qm.coordinate && qm.count == N) {
            flow.atom_xyz = (double*)md_alloc(arena, sizeof(double) * 3 * N);
            for (size_t a = 0; a < N; ++a) {
                flow.atom_xyz[a * 3 + 0] = qm.coordinate[a].x;
                flow.atom_xyz[a * 3 + 1] = qm.coordinate[a].y;
                flow.atom_xyz[a * 3 + 2] = qm.coordinate[a].z;
            }
        }
        md_array_resize(flow.strip, num_excited, arena);
        for (size_t s = 0; s < num_excited; ++s) flow.strip[s] = {};
    }

    // One state at a time: collect what finished, start the next (the selected one first). Asking for a
    // state's NTOs starts reading them on a worker (es_nto_ready); the strip job starts once they are in.
    void flow_step_strip(ApplicationState& state) {
        flow_build_strip_inputs(state.mold.sys);
        if (!flow.ao_unit || md_array_size(flow.strip) == 0) return;

        if (strip_job && strip_job->done.load(std::memory_order_acquire)) {
            if ((size_t)strip_job->state < md_array_size(flow.strip)) {
                QuantumChemistry::Flow::StripEntry& e = flow.strip[strip_job->state];
                e.done = true;
                e.ok = strip_job->ok;
                e.ct = strip_job->ct;
                e.distance = strip_job->distance;
                e.pr = strip_job->pr;
            }
            strip_job_destroy(strip_job);
            strip_job = nullptr;
        }
        if (strip_job) return;

        int32_t next = -1;
        if (rsp.selected >= 0 && (size_t)rsp.selected < md_array_size(flow.strip) && !flow.strip[rsp.selected].done) {
            next = rsp.selected;
        } else {
            for (size_t s = 0; s < md_array_size(flow.strip); ++s) {
                if (!flow.strip[s].done) { next = (int32_t)s; break; }
            }
        }
        if (next < 0 || !es_nto_ready(&state, (size_t)next)) return;
        strip_job = strip_job_start(state.mold.sys, lowdin, next, flow.ao_unit, flow.num_units, flow.atom_xyz);
        if (!strip_job) flow.strip[next].done = true;   // nothing to compute it from: not to be retried
    }

    // Every excited state at a glance, above the diagram: bar height the oscillator strength, colour the
    // charge transfer between units (grey-blue local, orange transfer), hollow while not yet computed.
    // Click to select a state.
    void draw_flow_strip(ApplicationState& state) {
        const size_t S = md_array_size(flow.strip);
        if (S == 0) return;
        size_t num_f = 0, num_e = 0;
        const double* osc = attribute_series_f64(&num_f, state.mold.sys, MD_STR_LIT("vlx/rsp/oscillator_strength"));
        const double* ene = attribute_series_f64(&num_e, state.mold.sys, MD_STR_LIT("vlx/rsp/frequency"));
        double f_max = 0.0;
        for (size_t s = 0; s < MD_MIN(S, num_f); ++s) f_max = MD_MAX(f_max, osc[s]);

        const float h = 46.0f;
        const ImVec2 p0 = ImGui::GetCursorScreenPos();
        const float w = MD_MAX(ImGui::GetContentRegionAvail().x, 64.0f);
        ImGui::InvisibleButton("##state_strip", ImVec2(w, h));
        const bool hovered = ImGui::IsItemHovered();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        dl->AddRectFilled(p0, p0 + ImVec2(w, h), IM_COL32(250, 250, 250, 255));

        const float slot = w / (float)S;
        const float bar_w = MD_CLAMP(slot * 0.6f, 2.0f, 14.0f);
        const float base_y = p0.y + h - 12.0f;
        const float max_h = h - 16.0f;
        const vec4_t local_col = { 0.55f, 0.60f, 0.70f, 1.0f };
        const vec4_t ct_col    = { 0.88f, 0.42f, 0.12f, 1.0f };
        const ImVec2 mouse = ImGui::GetMousePos();
        int32_t hover = -1;

        for (size_t s = 0; s < S; ++s) {
            const float cx = p0.x + slot * ((float)s + 0.5f);
            const double f = (s < num_f && f_max > 0.0) ? osc[s] / f_max : 0.0;
            // square root: a dark state is still a visible stub, a bright one does not flatten the rest
            const float bh = MD_MAX(3.0f, (float)sqrt(MD_MAX(f, 0.0)) * max_h);
            const ImVec2 b0 = { cx - bar_w * 0.5f, base_y - bh };
            const ImVec2 b1 = { cx + bar_w * 0.5f, base_y };
            const QuantumChemistry::Flow::StripEntry& e = flow.strip[s];
            if (e.done && e.ok) {
                const vec4_t c = vec4_lerp(local_col, ct_col, MD_CLAMP(e.ct, 0.0f, 1.0f));
                dl->AddRectFilled(b0, b1, ImGui::ColorConvertFloat4ToU32(vec_cast(c)), 1.5f);
            } else {
                dl->AddRect(b0, b1, IM_COL32(150, 150, 150, 255), 1.5f);
            }
            if ((int32_t)s == rsp.selected) {
                dl->AddTriangleFilled({ cx, base_y + 2.0f }, { cx - 4.0f, base_y + 9.0f }, { cx + 4.0f, base_y + 9.0f }, IM_COL32(30, 30, 30, 255));
            }
            if (hovered && mouse.x >= p0.x + slot * (float)s && mouse.x < p0.x + slot * (float)(s + 1)) hover = (int32_t)s;
        }

        if (hover >= 0) {
            const QuantumChemistry::Flow::StripEntry& e = flow.strip[hover];
            char details[160] = "computing...";
            if (e.done && e.ok) {
                snprintf(details, sizeof(details), "Charge transfer %.0f%% (between units)\nHole to particle %.2f \xc3\x85\nNTO pairs %.1f",
                         100.0f * e.ct, e.distance, e.pr);
            } else if (e.done) {
                snprintf(details, sizeof(details), "no NTO pairs");
            }
            const double ev = (size_t)hover < num_e ? ene[hover] * HARTREE_TO_EV : 0.0;
            const double fo = (size_t)hover < num_f ? osc[hover] : 0.0;
            ImGui::SetTooltip("S%d   %.3f eV   f = %.4f\n%s\n\nClick to select", hover + 1, ev, fo, details);
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) rsp.selected = hover;
        }
    }

    // ---------------------------------------------------------------------------
    // The side views of the charge transfer window (CtViews)
    // ---------------------------------------------------------------------------

    size_t ctv_entry_bytes() const {
        // two volumes of half floats
        return (size_t)ctv.grid.dim[0] * (size_t)ctv.grid.dim[1] * (size_t)ctv.grid.dim[2] * 2 * 2;
    }

    static uint64_t ctv_density_key(int32_t s) { return (uint64_t)(uint32_t)s + 1; }    // never 0, top bit clear
    static uint64_t ctv_orbital_key(uint64_t h) { return h | (1ull << 63); }

    CtViews::Entry* ctv_find_key(uint64_t key) {
        if (!key) return nullptr;
        for (size_t i = 0; i < md_array_size(ctv.cache); ++i) {
            if (ctv.cache[i].key == key) return &ctv.cache[i];
        }
        return nullptr;
    }
    CtViews::Entry* ctv_find(int32_t s) { return ctv_find_key(ctv_density_key(s)); }

    // Both volumes filled: a GPU evaluation lands in its texture frames after it was asked for, and every landing
    // moves the texture's data version on from the one it had when it was made
    static bool ctv_entry_ready(const CtViews::Entry& e) {
        if (e.failed) return false;
        for (int c = 0; c < 2; ++c) {
            if (!e.vol[c].tex_id || volume::data_version(e.vol[c].tex_id) == e.v0[c]) return false;
        }
        return true;
    }

    static void ctv_entry_free(CtViews::Entry* e) {
        for (int c = 0; c < 2; ++c) if (e->vol[c].tex_id) gl::free_texture(&e->vol[c].tex_id);
    }

    // Room for one more entry: the least recently used ready entries out, never the one asked for or one on
    // screen, never one whose evaluation is still under way (its readback would land in a texture gone). With
    // evict false nothing goes, which is what a prefetch asks: it is never worth what is already there.
    bool ctv_make_room(uint64_t keep, bool evict) {
        const size_t budget = (size_t)MD_MAX(ctv.budget_mb, 16) << 20;
        const size_t entry = ctv_entry_bytes();
        for (;;) {
            size_t used = 0;
            for (size_t i = 0; i < md_array_size(ctv.cache); ++i) if (!ctv.cache[i].failed) used += entry;
            if (used + entry <= budget || md_array_size(ctv.cache) == 0) return true;
            if (!evict) return false;
            int64_t victim = -1;
            for (size_t i = 0; i < md_array_size(ctv.cache); ++i) {
                const CtViews::Entry& e = ctv.cache[i];
                if (e.failed || e.key == keep || e.key == ctv_density_key(ctv.shown) || e.key == ctv.orb_shown || !ctv_entry_ready(e)) continue;
                if (victim < 0 || e.used < ctv.cache[victim].used) victim = (int64_t)i;
            }
            if (victim < 0) return true;    // over the budget by what is in use: one more state, for a while
            ctv_entry_free(&ctv.cache[victim]);
            md_array_swap_back_and_pop(ctv.cache, (size_t)victim);
        }
    }

    // Once a frame while the views are on. Collects the density job that finished and puts its matrices on the
    // GPU, then starts the next: the selected state when it is not cached, otherwise, when nothing is asked and
    // the evaluation runs off this thread, the uncached state nearest the selected one that fits in the budget.
    void ctv_update(ApplicationState& state) {
        if (!ctv.initialized) return;
        md_allocator_i* heap = md_get_heap_allocator();

        if (ctv.job && ctv.job->done.load(std::memory_order_acquire)) {
            CtDensityJob* job = ctv.job;
            ctv.job = nullptr;
            if (!ctv_find(job->state) && ctv_make_room(ctv_density_key(job->state), job->state == rsp.selected)) {
                CtViews::Entry e = {};
                e.key = ctv_density_key(job->state);
                e.state = job->state;
                e.used = ++ctv.tick;
                e.failed = !job->ok;
                if (job->ok) {
                    for (int c = 0; c < 2; ++c) {
                        init_volume(&e.vol[c], ctv.grid, GL_R16F);
                        e.v0[c] = volume::data_version(e.vol[c].tex_id);
                        if (!density_matrix_evaluate(&state, e.vol[c].tex_id, ctv.grid, job->D[c], job->A, MD_GTO_OP_SET)) e.failed = true;
                    }
                    if (e.failed) ctv_entry_free(&e);
                }
                md_array_push(ctv.cache, e, heap);
            }
            ct_density_job_destroy(job);
        }
        if (ctv.job) return;

        const int32_t S = (int32_t)num_excited;
        int32_t next = -1;
        if (rsp.selected >= 0 && rsp.selected < S && !ctv_find(rsp.selected)) {
            next = rsp.selected;
        } else if (rsp.selected >= 0 && gpu_volume_evaluation_available(&state)) {
            // nothing waiting on its evaluation: one state at a time, and only with room to spare
            for (size_t i = 0; i < md_array_size(ctv.cache); ++i) {
                if (!ctv.cache[i].failed && !ctv_entry_ready(ctv.cache[i])) return;
            }
            if (!ctv_make_room(0, false)) return;
            for (int32_t d = 1; d < S && next < 0; ++d) {
                const int32_t c[2] = { rsp.selected + d, rsp.selected - d };
                for (int32_t s : c) if (s >= 0 && s < S && !ctv_find(s)) { next = s; break; }
            }
        }
        if (next < 0 || !es_nto_ready(&state, (size_t)next)) return;
        ctv.job = ct_density_job_start(state.mold.sys, next);
        if (!ctv.job) {
            // nothing to make the densities of: remembered, so that it is not asked again every frame
            CtViews::Entry e = {};
            e.key = ctv_density_key(next);
            e.state = next;
            e.failed = true;
            md_array_push(ctv.cache, e, heap);
        }
    }

    // ---- ribbon orbitals ----

    // An orbital the views can show: a ribbon's block NTO pair (fragment route) or an NTO pair of the state
    // (orbital route). key 0 for none.
    struct CtOrbReq {
        uint64_t key = 0;
        int32_t  link = -1;     // cut link (fragment route)
        uint32_t pair = 0;      // block NTO pair of the ribbon, descending
        int32_t  nto = -1;      // NTO pair index of the state (orbital route)
    };

    bool ctv_block_ready() const {
        return flow.mode == FLOW_MODE_FRAGMENTS && flow.built && ct_state_job_ready(ct_job) && ct_job->state == rsp.selected && lowdin_job_ready(lowdin);
    }

    // The NTO pair a node of the orbital route's middle column stands for, -1 for any other node
    int32_t ctv_nto_of_cut_node(uint32_t cut_idx) const {
        const md_flow_node_t* n = md_flow_cut_node(&flow.cut, &flow.graph, cut_idx);
        if (!n || (n->flags & MD_FLOW_NODE_FLAG_OTHER) || n->level != FLOW_LEVEL_NTO) return -1;
        const uint32_t k = flow_key_index(n->key);
        return k < MAX_NTO_LABELS ? (int32_t)flow.nto_index[k] : -1;
    }

    CtOrbReq ctv_request(int32_t link, int32_t stripe, int32_t node) const {
        CtOrbReq r = {};
        if (rsp.selected < 0) return r;
        if (flow.mode == FLOW_MODE_FRAGMENTS) {
            if (link < 0 || (size_t)link >= md_array_size(flow.cut.links) || !ctv_block_ready()) return r;
            const md_flow_link_t* l = flow.cut.links + link;
            const md_flow_node_t* a = md_flow_cut_node(&flow.cut, &flow.graph, l->src);
            const md_flow_node_t* b = md_flow_cut_node(&flow.cut, &flow.graph, l->dst);
            if (!a || !b || (a->flags & MD_FLOW_NODE_FLAG_OTHER) || (b->flags & MD_FLOW_NODE_FLAG_OTHER)) return r;
            if (!flow_level_has_atoms(a->level) || !flow_level_has_atoms(b->level)) return r;
            r.link = link;
            r.pair = stripe > 0 ? (uint32_t)stripe : 0;
            const uint64_t h[4] = { (uint64_t)rsp.selected, a->key, b->key, r.pair };
            r.key = ctv_orbital_key(md_hash64(h, sizeof(h), 0xb10c4a0ull));
        } else if (flow.mode == FLOW_MODE_ORBITALS) {
            int32_t k = node >= 0 ? ctv_nto_of_cut_node((uint32_t)node) : -1;
            if (k < 0 && link >= 0 && (size_t)link < md_array_size(flow.cut.links)) {
                k = ctv_nto_of_cut_node(flow.cut.links[link].src);
                if (k < 0) k = ctv_nto_of_cut_node(flow.cut.links[link].dst);
            }
            if (k < 0) return r;
            r.nto = k;
            const uint64_t h[2] = { (uint64_t)rsp.selected, (uint64_t)k };
            r.key = ctv_orbital_key(md_hash64(h, sizeof(h), 0x7e0ull));
        }
        return r;
    }

    // The pinned orbital as a request, its ribbon found again among the visible ones by its ends
    CtOrbReq ctv_pin_request() const {
        const CtViews::Pin& p = ctv.pin;
        if (!p.on || p.state != rsp.selected) return {};
        if (flow.mode == FLOW_MODE_ORBITALS) {
            if (p.nto < 0) return {};
            for (size_t i = 0; i < md_array_size(flow.graph.nodes); ++i) {
                if (flow.graph.nodes[i].level == FLOW_LEVEL_NTO && ctv_nto_of_cut_node((uint32_t)i) == p.nto) return ctv_request(-1, -1, (int32_t)i);
            }
            return {};
        }
        const int32_t l = ctv_pin_link();
        return l >= 0 ? ctv_request(l, (int32_t)p.pair, -1) : CtOrbReq{};
    }

    int32_t ctv_pin_link() const {
        const CtViews::Pin& p = ctv.pin;
        if (!p.on || p.state != rsp.selected || flow.mode != FLOW_MODE_FRAGMENTS) return -1;
        for (size_t i = 0; i < md_array_size(flow.cut.links); ++i) {
            const md_flow_node_t* a = md_flow_cut_node(&flow.cut, &flow.graph, flow.cut.links[i].src);
            const md_flow_node_t* b = md_flow_cut_node(&flow.cut, &flow.graph, flow.cut.links[i].dst);
            if (a && b && a->key == p.src_key && b->key == p.dst_key) return (int32_t)i;
        }
        return -1;
    }

    void ctv_pin_toggle(const CtOrbReq& r) {
        const CtOrbReq cur = ctv_pin_request();
        if (cur.key && cur.key == r.key) { ctv.pin = {}; return; }
        ctv.pin = {};
        if (!r.key) return;
        ctv.pin.on = true;
        ctv.pin.state = rsp.selected;
        ctv.pin.nto = r.nto;
        ctv.pin.pair = r.pair;
        if (r.link >= 0) {
            const md_flow_node_t* a = md_flow_cut_node(&flow.cut, &flow.graph, flow.cut.links[r.link].src);
            const md_flow_node_t* b = md_flow_cut_node(&flow.cut, &flow.graph, flow.cut.links[r.link].dst);
            ctv.pin.src_key = a ? a->key : 0;
            ctv.pin.dst_key = b ? b->key : 0;
        }
    }

    // A ribbon's block NTO pair in the Lowdin basis: the pair'th singular pair of the ribbon's block of the
    // transition density, each confined to its end's atoms. x_hole and x_part are A doubles; false when the
    // ribbon has no such pair.
    bool ctv_block_pair(double* x_hole, double* x_part, char info[2][64], char* sub, size_t sub_cap, int32_t link, uint32_t pair) {
        const size_t K = ct_job->K, KK = K * K, A = lowdin->num_ao, N = ct_job->num_atoms;
        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };
        md_allocator_i* talloc = md_temp_allocator(temp);

        const md_flow_link_t* l = flow.cut.links + link;
        const uint32_t end[2] = { l->src, l->dst };
        const double* x_lowdin[2] = { ct_job->hole, ct_job->part };
        const double* blocks[2]   = { ct_job->hole_blocks, ct_job->part_blocks };
        double* B[2];
        bool* in[2];
        for (int side = 0; side < 2; ++side) {
            md_array(uint32_t) atoms = 0;
            if (!flow_node_qm_atoms(&atoms, end[side], talloc) || md_array_size(atoms) == 0) return false;
            B[side]  = md_temp_alloc_zero_array(temp, double, KK);
            in[side] = md_temp_alloc_zero_array(temp, bool, N + 1);
            for (size_t i = 0; i < md_array_size(atoms); ++i) {
                if (atoms[i] >= N) continue;
                in[side][atoms[i]] = true;
                const double* ab = blocks[side] + (size_t)atoms[i] * KK;
                for (size_t k = 0; k < KK; ++k) B[side][k] += ab[k];
            }
        }
        double* w  = md_temp_alloc_array(temp, double, K);
        double* hc = md_temp_alloc_array(temp, double, KK);
        double* pc = md_temp_alloc_array(temp, double, KK);
        const size_t n = md_qm_ct_block_nto(w, hc, pc, ct_job->lambda, K, B[0], B[1]);
        if (n == 0) return false;
        const size_t j = MD_MIN((size_t)pair, n - 1);
        double sum = 0.0;
        for (size_t i = 0; i < n; ++i) sum += w[i];

        const double* coef[2] = { hc + j * K, pc + j * K };
        double* out[2] = { x_hole, x_part };
        for (int side = 0; side < 2; ++side) {
            for (size_t mu = 0; mu < A; ++mu) {
                const uint32_t atom = lowdin->ao_atom[mu];
                double v = 0.0;
                if (atom < N && in[side][atom]) {
                    for (size_t k = 0; k < K; ++k) v += coef[side][k] * x_lowdin[side][k * A + mu];
                }
                out[side][mu] = v;
            }
            const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, end[side]);
            snprintf(info[side], 64, "%.*s", node ? (int)node->label.len : 0, node && node->label.ptr ? node->label.ptr : "");
        }
        if (n > 1) {
            snprintf(sub, sub_cap, "pair %zu of %zu, %.0f%% of the ribbon (%.1f%%)", j + 1, n, sum > 0.0 ? 100.0 * w[j] / sum : 0.0, 100.0 * l->weight);
        } else {
            snprintf(sub, sub_cap, "the ribbon's one pair (%.1f%%)", 100.0 * l->weight);
        }
        return true;
    }

    // Once a frame while the views are on, after ctv_update: collects the orbital job that finished, and starts
    // the one asked for when it is neither cached nor under way. The orbital route's NTO pairs need no job: they
    // are published as AO coefficients and go straight to the evaluation.
    void ctv_orbital_update(ApplicationState& state, const CtOrbReq& req) {
        if (!ctv.initialized) return;
        md_allocator_i* heap = md_get_heap_allocator();

        auto evaluate = [&](CtViews::Entry* e, const double* c0, const double* c1, size_t A) {
            const double* c[2] = { c0, c1 };
            for (int i = 0; i < 2 && !e->failed; ++i) {
                init_volume(&e->vol[i], ctv.grid, GL_R16F);
                e->v0[i] = volume::data_version(e->vol[i].tex_id);
                if (!orbital_coefficients_evaluate(&state, e->vol[i].tex_id, ctv.grid, c[i], A, MD_GTO_EVAL_MODE_PSI, MD_GTO_OP_SET, DEFAULT_GTO_CUTOFF_VALUE)) e->failed = true;
            }
            if (e->failed) ctv_entry_free(e);
        };

        if (ctv.orb_job && ctv.orb_job->done.load(std::memory_order_acquire)) {
            CtOrbitalJob* job = ctv.orb_job;
            ctv.orb_job = nullptr;
            if (!ctv_find_key(job->key) && ctv_make_room(job->key, true)) {
                CtViews::Entry e = {};
                e.key = job->key;
                e.state = rsp.selected;
                e.orbital = true;
                e.used = ++ctv.tick;
                e.failed = !job->ok;
                md_memcpy(e.info, job->info, sizeof(e.info));
                md_memcpy(e.sub, job->sub, sizeof(e.sub));
                if (job->ok) evaluate(&e, job->c[0], job->c[1], job->lowdin->num_ao);
                md_array_push(ctv.cache, e, heap);
            }
            ct_orbital_job_destroy(job);
        }

        if (!req.key || ctv_find_key(req.key)) return;
        if (req.nto >= 0) {
            // the state's own NTO pair: published, so evaluated from the attribute table
            if (!es_nto_ready(&state, (size_t)rsp.selected) || !ctv_make_room(req.key, true)) return;
            CtViews::Entry e = {};
            e.key = req.key;
            e.state = rsp.selected;
            e.orbital = true;
            e.used = ++ctv.tick;
            snprintf(e.info[0], sizeof(e.info[0]), "NTO %d", req.nto + 1);
            snprintf(e.info[1], sizeof(e.info[1]), "NTO %d", req.nto + 1);
            const md_attribute_t* lam = md_attributes_find(&state.mold.sys.attributes, es_path::nto_lambda);
            double lambda = 0.0;
            if (lam && md_attribute_extract_f64(&lambda, 1, lam, md_attribute_slice_2((uint32_t)rsp.selected, (uint32_t)req.nto), md_unit_none()) == 1) {
                snprintf(e.sub, sizeof(e.sub), "\xCE\xBB = %.3f", lambda);
            }
            for (int i = 0; i < 2 && !e.failed; ++i) {
                init_volume(&e.vol[i], ctv.grid, GL_R16F);
                e.v0[i] = volume::data_version(e.vol[i].tex_id);
                const ElectronicStructureNtoComponent comp = i == 0 ? ElectronicStructureNtoComponent::Hole : ElectronicStructureNtoComponent::Particle;
                if (!evaluate_nto(state, e.vol[i].tex_id, ctv.grid, (size_t)rsp.selected, (size_t)req.nto, comp, MD_GTO_OP_SET)) e.failed = true;
            }
            if (e.failed) ctv_entry_free(&e);
            md_array_push(ctv.cache, e, heap);
            return;
        }
        if (ctv.orb_job || req.link < 0 || !ctv_block_ready()) return;

        const size_t A = lowdin->num_ao;
        double* x_hole = (double*)md_alloc(heap, sizeof(double) * A);
        double* x_part = (double*)md_alloc(heap, sizeof(double) * A);
        char info[2][64] = {};
        char sub[96] = {};
        if (!ctv_block_pair(x_hole, x_part, info, sub, sizeof(sub), req.link, req.pair)) {
            md_free(heap, x_hole, sizeof(double) * A);
            md_free(heap, x_part, sizeof(double) * A);
            CtViews::Entry e = {};
            e.key = req.key;
            e.state = rsp.selected;
            e.orbital = true;
            e.failed = true;
            md_array_push(ctv.cache, e, heap);
            return;
        }
        ctv.orb_job = ct_orbital_job_start(lowdin, req.key, x_hole, x_part);
        md_memcpy(ctv.orb_job->info, info, sizeof(info));
        md_memcpy(ctv.orb_job->sub, sub, sizeof(sub));
    }

    // The two views, hole and particle, in r[0] and r[1]: the molecule, the selection and highlight marks over it,
    // and each its own density - or, while an orbital is asked for (orb_key), its hole and particle orbitals.
    // keep_orbital: the cursor is on something with orbitals, so the last orbital shown stays until the one asked
    // for is in, rather than the densities flashing between them. One camera for both views, so they turn
    // together and read side by side.
    void draw_ct_views(ApplicationState& state, const ImRect r[2], uint64_t orb_key, bool keep_orbital) {
        ImGuiIO& io = ImGui::GetIO();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImVec2 sz = r[0].GetSize();
        if (sz.x < 16.0f || sz.y < 16.0f || !ctv.initialized) return;

        camera_animate(&ctv.camera, ctv.target, state.app.timing.delta_s);

        // What to show: the selected state when its volumes are in, the last state shown until then
        CtViews::Entry* sel = rsp.selected >= 0 ? ctv_find(rsp.selected) : nullptr;
        if (sel && ctv_entry_ready(*sel)) {
            ctv.shown = rsp.selected;
        }
        CtViews::Entry* shown = ctv.shown >= 0 ? ctv_find(ctv.shown) : nullptr;
        if (shown && !ctv_entry_ready(*shown)) shown = nullptr;
        if (sel && sel->failed) shown = nullptr;    // the last state's densities would read as this one's
        bool updating = ctv.shown != rsp.selected && !(sel && sel->failed);
        const char* failed_status = (sel && sel->failed) ? "no NTO pairs" : nullptr;

        // An orbital over the densities, when one is asked for
        if (orb_key || keep_orbital) {
            CtViews::Entry* want = ctv_find_key(orb_key);
            if (want && ctv_entry_ready(*want)) {
                ctv.orb_shown = orb_key;
                shown = want;
                updating = false;
                failed_status = nullptr;
            } else if (want && want->failed) {
                shown = nullptr;
                updating = false;
                failed_status = "no orbital pair";
            } else {
                CtViews::Entry* last = keep_orbital ? ctv_find_key(ctv.orb_shown) : nullptr;
                if (last && ctv_entry_ready(*last) && last->state == rsp.selected) shown = last;
                updating = orb_key != 0;
            }
        }
        if (shown) shown->used = ++ctv.tick;
        const bool orbital = shown && shown->orbital;

        const float aspect = sz.x / sz.y;
        const mat4_t view_mat     = camera_world_to_view_matrix(ctv.camera);
        const mat4_t proj_mat     = camera_view_to_clip_matrix_persp(ctv.camera, aspect);
        const mat4_t inv_view_mat = camera_view_to_world_matrix(ctv.camera);
        const mat4_t inv_proj_mat = camera_clip_to_view_matrix_persp(ctv.camera, aspect);
        const mat4_t world_to_clip = proj_mat * view_mat;
        const mat4_t clip_to_world = inv_view_mat * inv_proj_mat;

        const int width  = MD_MAX(1, (int)(sz.x * io.DisplayFramebufferScale.x));
        const int height = MD_MAX(1, (int)(sz.y * io.DisplayFramebufferScale.y));
        // Sized before the images go into the draw list, which records the texture names
        if ((int)ctv.gbuf.width != width || (int)ctv.gbuf.height != height) {
            gbuffer_init(&ctv.gbuf, width, height);
            for (int i = 0; i < 2; ++i) {
                if (ctv.iso_tex[i]) gl::free_texture(&ctv.iso_tex[i]);
                if (ctv.out_tex[i]) gl::free_texture(&ctv.out_tex[i]);
            }
        }
        for (int i = 0; i < 2; ++i) {
            if (!ctv.iso_tex[i]) gl::init_texture_2D(&ctv.iso_tex[i], width, height, GL_RGBA16F);
            if (!ctv.out_tex[i]) gl::init_texture_2D(&ctv.out_tex[i], width, height, GL_RGBA8);
        }

        vec4_t colors[2] = { flow.style.side_color[0], flow.style.side_color[2] };
        colors[0].w = colors[1].w = 0.75f;
        const float iso_value = (float)(ctv.iso_val * ctv.iso_val);
        // orbitals in the orbital window's lobe colours, at +-iso
        const float  orb_values[2] = { (float)ctv.iso_val, -(float)ctv.iso_val };
        const vec4_t orb_colors[2] = { orb.iso.colors[0], orb.iso.colors[1] };

        // Camera control in either view; the picking reads the shared buffer, the same in both
        for (int i = 0; i < 2; ++i) {
            ImGui::PushID(i + 0x5c7);
            ImGui::SetCursorScreenPos(r[i].Min);
            InteractionSurfaceState surface = interaction_surface(interaction_surface_flow, vec_cast(sz));
            InteractionSurfaceViewTransformArgs view_args = {
                .camera = ctv.camera,
                .trackball_param = state.view.trackball_param,
            };
            InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&ctv.target, surface, view_args);
            if (view_result.reset_requested) ctv.target = default_view;
            if (surface.hovered) {
                InteractionSurfaceHitArgs hit_args = {
                    .picking_surface = &ctv.picking_surface,
                    .picking_handler = state.picking_handler,
                    .fbo = ctv.gbuf.fbo,
                    .width = ctv.gbuf.width,
                    .height = ctv.gbuf.height,
                    .clip_to_world = clip_to_world,
                };
                PickingHit hit = {};
                interaction_surface_hit_extract(&hit, surface, hit_args);
                InteractionSurfaceEvent event = {};
                interaction_surface_event_extract(&event, surface, hit);
                event.clip_to_world = clip_to_world;
                event.world_to_clip = world_to_clip;
                if (event.kind != InteractionSurfaceEventKind::ContextMenu && event.kind != InteractionSurfaceEventKind::RegionSelect) {
                    viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
                }
            }
            ImGui::PopID();
        }

        uint64_t render_hash = md_hash64(&view_mat, sizeof(view_mat), (uint64_t)ctv.gl_rep.id);
        {
            const int layout[2] = { width, height };
            render_hash = md_hash64(layout, sizeof(layout), render_hash);
            render_hash = md_hash64(&proj_mat, sizeof(proj_mat), render_hash);
            const uint64_t s = shown ? shown->key : 0;
            render_hash = md_hash64(&s, sizeof(s), render_hash);
            if (shown) {
                for (int c = 0; c < 2; ++c) {
                    const uint64_t v = volume::data_version(shown->vol[c].tex_id);
                    render_hash = md_hash64(&shown->vol[c].tex_id, sizeof(shown->vol[c].tex_id), render_hash);
                    render_hash = md_hash64(&v, sizeof(v), render_hash);
                }
            }
            render_hash = md_hash64(colors, sizeof(colors), render_hash);
            render_hash = md_hash64(&iso_value, sizeof(iso_value), render_hash);
            render_hash = md_hash64(orb_colors, sizeof(orb_colors), render_hash);
            render_hash = md_hash64(&state.settings.exact_isosurfaces, sizeof(state.settings.exact_isosurfaces), render_hash);
            render_hash = md_hash64(&state.visuals.background, sizeof(state.visuals.background), render_hash);
            render_hash = md_hash64(&state.visuals.tonemapping, sizeof(state.visuals.tonemapping), render_hash);
            render_hash = md_hash64(&state.visuals.shading, sizeof(state.visuals.shading), render_hash);
            render_hash = md_hash64(&state.selection.color, sizeof(state.selection.color), render_hash);
            render_hash = md_hash64(&state.mold.gpu_buffers_version, sizeof(state.mold.gpu_buffers_version), render_hash);
            const md_gl_picking_offset_t picking_offset = gl_picking_offset(state);
            render_hash = md_hash64(&picking_offset, sizeof(picking_offset), render_hash);
            render_hash = render_hash ? render_hash : 1;
        }

        if (render_hash != ctv.render_hash) {
            ctv.render_hash = render_hash;
            GBuffer& gbuf = ctv.gbuf;
            gbuffer_clear(&gbuf);

            const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                                            GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LESS);
            glDepthMask(GL_TRUE);
            glEnable(GL_SCISSOR_TEST);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
            glDrawBuffers((int)MD_COUNTOF(draw_buffers), draw_buffers);
            glViewport(0, 0, gbuf.width, gbuf.height);
            glScissor(0, 0, gbuf.width, gbuf.height);

            auto draw_rep = [&state, &view_mat, &proj_mat](md_gl_rep_t rep, md_gl_shaders_t& shaders, uint32_t atom_mask) {
                md_gl_draw_op_t draw_op = {};
                draw_op.type = MD_GL_REP_BALL_AND_STICK;
                draw_op.args.ball_and_stick.ball_scale   = 1.0f;
                draw_op.args.ball_and_stick.stick_radius = 1.0f;
                draw_op.args.ball_and_stick.bond_orders  = true;
                draw_op.rep = rep;
                md_gl_draw_args_t draw_args = {
                    .shaders = shaders,
                    .draw_operations = { .count = 1, .ops = &draw_op },
                    .view_transform = {
                        .view_matrix = (const float*)view_mat.elem,
                        .proj_matrix = (const float*)proj_mat.elem,
                    },
                    .picking_offset = gl_picking_offset(state),
                    .atom_mask = atom_mask,
                };
                md_gl_draw(&draw_args);
            };
            draw_rep(ctv.gl_rep, state.gl.shaders, 0);

            const bool draw_marks = !md_bitmask_empty(&state.selection.selection_mask) || !md_bitmask_empty(&state.selection.highlight_mask);
            if (draw_marks) {
                gbuffer_marks_begin(&gbuf);
                draw_rep(ctv.gl_rep, state.gl.shaders_marks, AtomBit_Marked);
                gbuffer_marks_end(&gbuf);
            }

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
            glViewport(0, 0, gbuf.width, gbuf.height);
            glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
            glClearColor(0, 0, 0, 0);
            glClear(GL_COLOR_BUFFER_BIT);
            glDepthFunc(GL_LESS);
            glDepthMask(0);
            glColorMask(1, 1, 1, 1);

            postprocess_pipeline::Settings pp = {};
            postprocess_pipeline::Inputs in = {};
            pp.background_color   = state.visuals.background.color;
            pp.shading.mode       = state.visuals.shading.mode;
            pp.shading.ambient    = state.visuals.shading.ambient;
            pp.shading.specular   = state.visuals.shading.specular;
            pp.shading.colors     = state.visuals.shading.colors;
            pp.shading.light      = state.visuals.shading.light;
            pp.shading.environment = state.visuals.shading.environment;
            pp.tonemap.enabled    = state.visuals.tonemapping.enabled;
            pp.tonemap.mode       = state.visuals.tonemapping.tonemapper;
            pp.tonemap.exposure   = state.visuals.tonemapping.exposure;
            pp.tonemap.gamma      = state.visuals.tonemapping.gamma;
            pp.shadows.enabled = false;
            pp.ssao.enabled    = false;
            pp.dof.enabled     = false;
            pp.fxaa.enabled    = true;
            pp.taa.enabled     = false;
            pp.sharpen.enabled = false;
            in.depth        = gbuf.tex.depth;
            in.color        = gbuf.tex.color;
            in.normal       = gbuf.tex.normal;
            in.velocity     = gbuf.tex.velocity;
            in.transparency = gbuf.tex.transparency;
            if (draw_marks) {
                in.marks       = gbuf.tex.marks;
                in.marks_group = gbuf.tex.marks_group;
                marks_settings(&pp, state, 1.0f, MD_MAX(1, (int)roundf(io.DisplayFramebufferScale.x)));
            }
            ViewParam view_param = {
                .matrix = {
                    .curr = { .view = view_mat, .proj = proj_mat, .norm = view_mat },
                    .inv  = { .proj = inv_proj_mat },
                },
                .clip_planes = { .near = ctv.camera.near_plane, .far = ctv.camera.far_plane },
                .resolution = { sz.x, sz.y },
                .fov_y = ctv.camera.fov_y,
            };

            glDisable(GL_SCISSOR_TEST);
            if (!ctv.out_fbo) glGenFramebuffers(1, &ctv.out_fbo);
            for (int i = 0; i < 2; ++i) {
                bool iso_written = false;
                if (shown) {
                    volume::IsoRenderDesc desc = {
                        .render_target = {
                            .depth  = gbuf.tex.depth,
                            .color  = ctv.iso_tex[i],
                            .width  = gbuf.width,
                            .height = gbuf.height,
                            .clear_color = true,
                        },
                        .texture = { .density_volume = shown->vol[i].tex_id },
                        .matrix = {
                            .model = shown->vol[i].texture_to_world,
                            .view  = view_mat,
                            .proj  = proj_mat,
                            .inv_proj = inv_proj_mat,
                        },
                        .iso = {
                            .count  = orbital ? 2u : 1u,
                            .values = orbital ? orb_values : &iso_value,
                            .colors = orbital ? orb_colors : &colors[i],
                            .exact  = state.settings.exact_isosurfaces,
                        },
                        .shading = {
                            .roughness = 0.3f,
                            .ior = 1.5f,
                            .mode = pp.shading.mode,
                            .ambient = pp.shading.ambient,
                            .specular = pp.shading.specular,
                            .colors = pp.shading.colors,
                            .light = pp.shading.light,
                            .environment = pp.shading.environment,
                        },
                    };
                    iso_written = volume::render_isosurfaces(desc);
                }
                in.transparency_hdr = iso_written ? ctv.iso_tex[i] : 0;
                glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ctv.out_fbo);
                glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ctv.out_tex[i], 0);
                glDrawBuffer(GL_COLOR_ATTACHMENT0);
                glViewport(0, 0, gbuf.width, gbuf.height);
                postprocess_pipeline::execute(in, pp, view_param);
            }
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glDrawBuffer(GL_BACK);
            glDisable(GL_SCISSOR_TEST);
        }

        static const char* lbl[2] = { "Hole", "Particle" };
        for (int i = 0; i < 2; ++i) {
            dl->AddImage((ImTextureID)(intptr_t)ctv.out_tex[i], r[i].Min, r[i].Max, { 0, 1 }, { 1, 0 });
            dl->AddRect(r[i].Min, r[i].Max, IM_COL32(0, 0, 0, 22));
            const ImVec2 tp = r[i].Min + ImVec2(6.0f, 4.0f);
            const ImU32 caption_col = ImGui::ColorConvertFloat4ToU32(vec_cast(vec4_t{ colors[i].x * 0.8f, colors[i].y * 0.8f, colors[i].z * 0.8f, 1.0f }));
            dl->AddText(tp, caption_col, lbl[i]);
            float y = tp.y + ImGui::GetTextLineHeight();
            if (orbital) {
                // what it is the orbital of, under the caption, clipped to the view
                dl->PushClipRect(r[i].Min, r[i].Max, true);
                dl->AddText(ImVec2(tp.x, y), IM_COL32(0, 0, 0, 200), shown->info[i]);
                y += ImGui::GetTextLineHeight();
                if (shown->sub[0]) dl->AddText(ImVec2(tp.x, y), IM_COL32(0, 0, 0, 130), shown->sub);
                dl->PopClipRect();
            }
            const char* status = failed_status ? failed_status : updating ? "updating..." : nullptr;
            if (status) {
                const ImVec2 ts = ImGui::CalcTextSize(status);
                dl->AddText(ImVec2(r[i].Max.x - ts.x - 6.0f, r[i].Min.y + 4.0f), IM_COL32(0, 0, 0, 140), status);
            }
            if (i == 0 && ctv.pin.on && orbital && shown->key == ctv_pin_request().key) {
                const char* pinned = "pinned (click it again, or Esc)";
                const ImVec2 ts = ImGui::CalcTextSize(pinned);
                dl->AddText(ImVec2(r[i].Min.x + 6.0f, r[i].Max.y - ts.y - 4.0f), IM_COL32(0, 0, 0, 120), pinned);
            }
        }
    }

    void draw_flow_window(ApplicationState& state) {
        if (!flow.show_window) return;

        const size_t num_excited_states = num_excited;
        if (num_excited_states == 0) return;

        if (!flow.cut.alloc) {
            md_flow_cut_init(&flow.cut, arena);
            md_flow_layout_init(&flow.layout, arena);
            flow.params = md_flow_layout_params_default();
            // The same fragments in the same order on both sides: a straight band is local, a
            // slanted one is transfer.
            flow.params.order = MD_FLOW_ORDER_GRAPH;
            // Thin bars: the labels sit outside them, and the width between belongs to the ribbons
            flow.params.node_thickness = flow.bar_width;
            flow.cut.threshold = flow.initial_threshold;
            // Uncapped on purpose. A cap makes "Show above 1%" stop being true - folding halts
            // early and the LARGEST sub-threshold rows stay on screen, with nothing to tell the
            // reader why. A big "Others" is fine now that it opens.
            flow.cut.other_max = 1.0f;
            flow_build_label_storage(state.mold.sys);
        }

        // Rebuild only when the numbers behind the picture actually changed. Hashing the matrix
        // itself rather than a bookkeeping counter keeps this window independent of whatever the
        // NTO window happens to be doing.
        uint64_t hash = md_hash64(nto.atom_group_idx, md_array_bytes(nto.atom_group_idx), 0);
        hash = md_hash64(nto.group.color, sizeof(nto.group.color), hash);
        hash ^= (uint64_t)nto.group.count;
        hash ^= (uint64_t)(rsp.selected + 1) << 8;
        hash ^= (uint64_t)flow.mode << 24;
        hash ^= flow.frags_version * 0x9E3779B97F4A7C15ull;

        // The fragment route: the fragments and the Lowdin basis once per system, then one job per selected
        // state. Both on workers; the last diagram stays on screen while they run.
        const bool nto_pending = rsp.selected >= 0 && !es_nto_ready(&state, (size_t)rsp.selected);
        bool lowdin_pending = false;
        bool ct_pending = false;
        if (flow.mode == FLOW_MODE_FRAGMENTS || flow.show_strip) {
            flow_build_fragments(state.mold.sys);
            if (!lowdin && flow.num_qm_atoms > 0) {
                lowdin = lowdin_job_start(state.mold.sys, basis, flow.num_qm_atoms);
            }
        }
        if (flow.show_strip) {
            flow_step_strip(state);
        }
        if (flow.mode == FLOW_MODE_FRAGMENTS) {
            lowdin_pending = lowdin && !lowdin->done.load(std::memory_order_acquire);
            if (lowdin_job_ready(lowdin) && rsp.selected >= 0 && !nto_pending && (!ct_job || ct_job->state != rsp.selected) &&
                flow.ct_failed_state != rsp.selected) {
                ct_state_job_destroy(ct_job);
                ct_job = ct_state_job_start(state.mold.sys, lowdin, rsp.selected);
                if (!ct_job) flow.ct_failed_state = rsp.selected;   // no NTO pairs to read: do not ask every frame
            }
            ct_pending = ct_job && (!ct_job->done.load(std::memory_order_acquire) || ct_job->state != rsp.selected);
            hash ^= (uint64_t)lowdin_job_ready(lowdin) << 32;
            if (ct_job) hash ^= ((uint64_t)(ct_job->state + 1) << 40) ^ ((uint64_t)ct_state_job_ready(ct_job) << 56);
        }

        const bool waiting = nto_pending || lowdin_pending || ct_pending;
        if (hash != flow.data_hash && !waiting) {
            flow.data_hash = hash;
            if (flow.mode == FLOW_MODE_ORBITALS) {
                // Gathered HERE and not inside: the flow graph procedures work on arrays and counts and
                // know nothing about a system, an attribute table or a reader.
                md_temp_scope_t temp = md_temp_begin();
                MD_DEFER { md_temp_end(temp); };
                FlowData qm = {};
                flow_qm_data_gather(&qm, temp, state.mold.sys);
                rebuild_flow_graph(&qm);
            } else {
                rebuild_flow_graph(nullptr);
            }
        }

        bool export_requested = false;

        ImGui::SetNextWindowSize(ImVec2(720, 420), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Charge Transfer Analysis", &flow.show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    export_requested = ImGui::MenuItem("Export Diagram (SVG)", nullptr, false, flow.built);
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Settings")) {
                    ImGui::Checkbox("Colour by side", &flow.style.color_by_side);
                    ImGui::SetItemTooltip("Hole nodes in one colour, particle nodes in another, ribbons in a third.\n"
                                          "Off: every node in its own colour, ribbons blending between them.");
                    ImGui::BeginDisabled(flow.style.color_by_side);
                    ImGui::Checkbox("Colour the labels", &flow.style.color_labels);
                    ImGui::EndDisabled();
                    ImGui::SliderFloat("Node tint", &flow.style.node_fill_tint, 0.0f, 0.6f);
                    ImGui::SetItemTooltip("How much of a node's colour fills its interior. 0 is white.");
                    ImGui::SliderFloat("Border width", &flow.style.border_width, 0.5f, 4.0f, "%.1f px");
                    ImGui::Checkbox("Ribbon labels", &flow.style.show_ribbon_labels);
                    ImGui::SetItemTooltip("The share of every ribbon thick enough to carry it, at its middle.");
                    ImGui::Checkbox("Legend", &flow.show_legend);
                    ImGui::Separator();
                    ImGui::SliderFloat("Ribbon opacity", &flow.style.ribbon_alpha, 0.05f, 1.0f);
                    ImGui::SliderFloat("Ribbon curvature", &flow.style.ribbon_curvature, 0.0f, 1.0f);
                    ImGui::SliderInt("Ribbon segments", &flow.style.ribbon_segments, 4, 64);
                    ImGui::Separator();
                    ImGui::SliderFloat("Bar width", &flow.params.node_thickness, 0.01f, 0.25f);
                    ImGui::SliderFloat("Node gap", &flow.params.node_gap, 0.0f, 0.1f);
                    ImGui::BeginDisabled(flow.params.order != MD_FLOW_ORDER_BARYCENTRE);
                    int sweeps = (int)flow.params.barycentre_sweeps;
                    if (ImGui::SliderInt("Crossing sweeps", &sweeps, 0, 4)) {
                        flow.params.barycentre_sweeps = (uint32_t)MD_MAX(0, sweeps);
                    }
                    ImGui::EndDisabled();
                    ImGui::SetItemTooltip("Only used by the 'Fewest crossings' row order.");
                    ImGui::Separator();
                    ImGui::Checkbox("State strip", &flow.show_strip);
                    ImGui::SetItemTooltip("Every excited state above the diagram: bar height the oscillator strength,\n"
                                          "colour the charge transfer between the automatic units.");
                    ImGui::SliderFloat("Local fill opacity", &flow.style.self_alpha, 0.0f, 1.0f);
                    ImGui::SetItemTooltip("How strongly the local excitation inside a node is filled in.");
                    {
                        float other_pct = flow.cut.other_max * 100.0f;
                        if (ImGui::SliderFloat("Others cap", &other_pct, 1.0f, 100.0f, "%.0f%%")) {
                            flow.cut.other_max = other_pct * 0.01f;
                            md_flow_cut_resolve(&flow.cut, &flow.graph);
                        }
                        ImGui::SetItemTooltip("The most a single 'Others' may hold. Below 100%% folding stops early,\n"
                                              "which means rows under 'Fold below' STAY on screen. Leave it at 100%%\n"
                                              "unless you want that.");
                    }
                    ImGui::SliderFloat("Hide ribbons under", &flow.style.ribbon_min_px, 0.0f, 6.0f, "%.1f px");
                    ImGui::SetItemTooltip("A ribbon thinner than this is not drawn - it reads as noise rather\n"
                                          "than as a quantity. The flow is still in the model and still counted.");
                    ImGui::Separator();
                    ImGui::Checkbox("3D views", &ctv.show);
                    ImGui::BeginDisabled(!ctv.show);
                    ImGui::SliderFloat("View width", &ctv.width_frac, 0.1f, 0.4f, "%.2f");
                    ImGui::SetItemTooltip("Of the window's width, each side.");
                    float iso = (float)ctv.iso_val;
                    if (ImGui::SliderFloat("Iso value", &iso, 0.005f, 0.2f, "%.3f", ImGuiSliderFlags_Logarithmic)) ctv.iso_val = iso;
                    ImGui::SetItemTooltip("Orbitals are drawn at plus and minus it, densities at its square, as the\n"
                                          "transition analysis draws them.");
                    ImGui::SliderInt("Cache", &ctv.budget_mb, 32, 2048, "%d MB");
                    ImGui::SetItemTooltip("GPU memory for the states' densities, kept so that switching state shows them at once.\n"
                                          "%zu states fit.", ctv.initialized && ctv_entry_bytes() ? ((size_t)ctv.budget_mb << 20) / ctv_entry_bytes() : (size_t)0);
                    ImGui::EndDisabled();
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }

            ImGui::PushItemWidth(160.0f);
            int state_idx = rsp.selected + 1;
            if (ImGui::SliderInt("State", &state_idx, 1, (int)num_excited_states)) {
                rsp.selected = state_idx - 1;
            }
            ImGui::PopItemWidth();

            ImGui::SameLine();
            ImGui::PushItemWidth(130.0f);
            int mode = (int)flow.mode;
            if (ImGui::Combo("Route", &mode, flow_mode_str, FLOW_MODE_COUNT)) {
                flow.mode = (flow_mode_t)mode;
            }
            ImGui::SetItemTooltip("'Fragments': where the excitation takes its electron from (left) and\n"
                                  "where it puts it (right), as the exact charge transfer numbers in the\n"
                                  "Lowdin basis. Local excitation fills the nodes; bands are transfer.\n"
                                  "'Orbitals': which canonical MOs the NTO pairs are built from.\n"
                                  "Every band there is exact too.");
            ImGui::PopItemWidth();

            const bool fragments = flow.mode == FLOW_MODE_FRAGMENTS;
            if (fragments) {
                ImGui::SameLine();
                ImGui::PushItemWidth(120.0f);
                if (ImGui::SliderInt("Rows", &flow.max_rows, 2, 40)) {
                    flow.cut.max_rows = (uint32_t)flow.max_rows;
                    md_flow_cut_resolve(&flow.cut, &flow.graph);
                }
                ImGui::SetItemTooltip("How many rows a column may show. Fragments open on their own, the\n"
                                      "heaviest first, until the next would not fit - so the detail goes where\n"
                                      "the hole and the particle are. What you open or close yourself stays so.");
                ImGui::PopItemWidth();
            }

            ImGui::SameLine();
            ImGui::PushItemWidth(120.0f);
            float threshold_pct = flow.cut.threshold * 100.0f;
            if (ImGui::SliderFloat("Fold below", &threshold_pct, 0.0f, 25.0f, "%.1f%%")) {
                flow.cut.threshold = threshold_pct * 0.01f;
                md_flow_cut_resolve(&flow.cut, &flow.graph);
            }
            ImGui::SetItemTooltip("Share of a COLUMN, not of the whole diagram. Siblings below it fold into\n"
                                  "an 'Others' inside their parent; double click that to draw them.");
            ImGui::PopItemWidth();

            if (!fragments) {
                ImGui::SameLine();
                ImGui::PushItemWidth(120.0f);
                ImGui::SliderFloat("Even sizes", &flow.size_normalize, 0.0f, 1.0f, "%.2f");
                ImGui::SetItemTooltip("Flattens the orbital columns' heights so a 1%% channel is still big\n"
                                      "enough to read and click. Ribbons keep their true widths and simply\n"
                                      "taper into the column.");
                ImGui::PopItemWidth();
            }

            if (fragments) {
                ImGui::SameLine();
                if (ImGui::Button("Expand all")) {
                    md_flow_cut_expand_all(&flow.cut, &flow.graph);
                    md_bitmask_clear(&flow.cut.collapsed);
                    md_flow_cut_resolve(&flow.cut, &flow.graph);
                }
                ImGui::SameLine();
                if (ImGui::Button("Collapse all")) {
                    md_bitmask_clear(&flow.cut.expanded);
                    for (size_t i = 0; i < md_array_size(flow.graph.nodes); ++i) {
                        if (!md_flow_node_is_leaf(&flow.graph, (uint32_t)i)) md_bitmask_set_bit(&flow.cut.collapsed, (uint32_t)i);
                    }
                    md_flow_cut_resolve(&flow.cut, &flow.graph);
                }
                ImGui::SameLine();
                if (ImGui::Button("Auto")) {
                    md_flow_cut_clear_pins(&flow.cut);
                    md_bitmask_clear(&flow.cut.other_open);
                    flow.cut.other_expanded = 0;
                    md_flow_cut_resolve(&flow.cut, &flow.graph);
                }
                ImGui::SetItemTooltip("Forget what was opened and closed by hand: the row budget decides again.");
            }

            ImGui::SameLine();
            ImGui::PushItemWidth(150.0f);
            int order = (int)flow.params.order;
            // In md_flow_order_t's order
            static const char* order_str[] = { "By contribution", "Fewest crossings", "As listed" };
            if (ImGui::Combo("Order", &order, order_str, IM_ARRAYSIZE(order_str))) {
                flow.params.order = (md_flow_order_t)order;
            }
            ImGui::SetItemTooltip("Rows within a column. 'As listed' keeps the fragments in the same order\n"
                                  "on both sides, so a straight band is local and a slanted one is transfer.\n"
                                  "'By contribution' puts the largest at the top. Either way a group's\n"
                                  "atoms stay together inside the group's slot.");
            ImGui::PopItemWidth();

            ImGui::SameLine();
            ImGui::Checkbox("3D", &ctv.show);
            ImGui::SetItemTooltip("The molecule on either side of the diagram: the hole density left, the particle\n"
                                  "density right, from the state's NTO pairs. The states are prepared in the\n"
                                  "background, so switching between them does not wait. Sizes under Settings.\n\n"
                                  "Rest the cursor on a ribbon for the orbital pair that carries it (the stripe\n"
                                  "under the cursor, where there are several); click it to keep it there, click\n"
                                  "it again, empty space or Esc to let go. In the orbital route, an NTO pair.");

            ImGui::SameLine();
            ImGui::TextDisabled("(?)");
            ImGui::SetItemTooltip("Left drag pans, wheel zooms at the cursor.\n"
                                  "Shift+left selects, shift+right deselects - click one item,\n"
                                  "or drag a box over several. Shift+right on empty space clears.\n"
                                  "Double click opens a node; double click on empty space refits.\n"
                                  "Right click for group actions on the selection.");


            if (flow.show_strip) {
                draw_flow_strip(state);
            }

            // The side views, when on, take a share of the width on either side: the hole density beside the hole
            // column, the particle density beside the particle column. Too narrow to read, they are left out.
            const ImVec2 area_p0 = ImGui::GetCursorScreenPos();
            const ImVec2 area_sz = ImGui::GetContentRegionAvail();
            const float side_gap = 4.0f;
            float side_w = 0.0f;
            if (ctv.show && ctv.initialized) {
                side_w = floorf(area_sz.x * MD_CLAMP(ctv.width_frac, 0.1f, 0.4f));
                if (side_w < 96.0f || area_sz.y < 64.0f) side_w = 0.0f;
            }
            const float side_take = side_w > 0.0f ? side_w + side_gap : 0.0f;
            const ImVec2 canvas_p0 = area_p0 + ImVec2(side_take, 0.0f);
            ImVec2 canvas_sz = area_sz - ImVec2(2.0f * side_take, 0.0f);
            canvas_sz.x = MD_MAX(canvas_sz.x, 64.0f);
            canvas_sz.y = MD_MAX(canvas_sz.y, 64.0f);
            const ImVec2 canvas_p1 = canvas_p0 + canvas_sz;

            if (side_w > 0.0f) {
                const ImRect side_rect[2] = {
                    { area_p0, area_p0 + ImVec2(side_w, canvas_sz.y) },
                    { ImVec2(canvas_p1.x + side_gap, area_p0.y), ImVec2(canvas_p1.x + side_gap + side_w, area_p0.y + canvas_sz.y) },
                };
                // The orbital asked for: the one under the cursor once it has rested there, else the pinned one.
                // From the last frame's hit, the diagram being drawn after the views.
                if (ctv.pin.on && ctv.pin.state != rsp.selected) ctv.pin = {};
                const CtOrbReq cand = ctv_request(flow.hover.link, flow.hover.stripe, flow.hover.node);
                const double now = ImGui::GetTime();
                if (cand.key != ctv.hover_cand) {
                    ctv.hover_cand = cand.key;
                    ctv.hover_t0 = now;
                }
                const CtOrbReq pinned = ctv_pin_request();
                CtOrbReq req = {};
                if (cand.key && now - ctv.hover_t0 >= (double)ctv.hover_delay) req = cand;
                else if (!cand.key) req = pinned;

                ctv_update(state);
                ctv_orbital_update(state, req);
                draw_ct_views(state, side_rect, req.key, cand.key != 0 || pinned.key != 0);
                ImGui::SetCursorScreenPos(canvas_p0);
            }

            ImGui::InvisibleButton("##flow_canvas", canvas_sz,
                ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle);
            const bool canvas_hovered = ImGui::IsItemHovered();
            const bool canvas_active  = ImGui::IsItemActive();

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddRectFilled(canvas_p0, canvas_p1, IM_COL32(255, 255, 255, 255));
            draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(0, 0, 0, 22));

            if (!flow.built) {
                const char* msg = lowdin_pending ? "Computing the Lowdin basis (once per system)..."
                                : nto_pending    ? "Computing the natural transition orbitals..."
                                : ct_pending     ? "Computing the charge transfer numbers..."
                                : (flow.mode == FLOW_MODE_FRAGMENTS && lowdin && !lowdin_job_ready(lowdin))
                                                 ? "The Lowdin basis could not be computed (see the log)"
                                                 : "No transition data for this state";
                draw_list->AddText(canvas_p0 + (canvas_sz - ImGui::CalcTextSize(msg)) * 0.5f, IM_COL32(0, 0, 0, 160), msg);
                ImGui::End();
                return;
            }

            ImGuiIO& io = ImGui::GetIO();

            // Zoom about the cursor. Screen = area + (p * zoom + pan) * size, so the graph point
            // under the pointer is p = (rel - pan) / zoom, and holding it fixed across a zoom
            // change means pan' = rel - p * zoom'. That is an ASSIGNMENT, not an increment - the
            // increment form leaves an extra pan term in and the view creeps toward the corner.
            if (canvas_hovered && io.MouseWheel != 0.0f) {
                const ImVec2 rel = (io.MousePos - canvas_p0) / canvas_sz;
                const vec2_t p = { (rel.x - flow.pan.x) / flow.zoom, (rel.y - flow.pan.y) / flow.zoom };
                flow.zoom = MD_CLAMP(flow.zoom * (1.0f + io.MouseWheel * 0.12f), 0.25f, 12.0f);
                flow.pan.x = rel.x - p.x * flow.zoom;
                flow.pan.y = rel.y - p.y * flow.zoom;
            }

            // Left drag pans. Selection moved onto shift, so the plain drag is free for the thing
            // people reach for first - and a shift drag is a region, so panning has to stand down
            // for it or the diagram would slide out from under the rubber band.
            if (canvas_active && !flow.region_active && !io.KeyShift &&
                                 (ImGui::IsMouseDragging(ImGuiMouseButton_Left) ||
                                  ImGui::IsMouseDragging(ImGuiMouseButton_Middle) ||
                                  ImGui::IsMouseDragging(ImGuiMouseButton_Right))) {
                flow.pan.x += io.MouseDelta.x / canvas_sz.x;
                flow.pan.y += io.MouseDelta.y / canvas_sz.y;
            }

            // Map the slider onto a per-column exponent. Orbital columns get it; fragment columns stay
            // proportional, because their proportions ARE the charge transfer.
            {
                const float e = 1.0f - flow.size_normalize * 0.85f;   // 1.0 .. 0.15
                for (uint32_t c = 0; c < MD_FLOW_MAX_COLUMNS; ++c) {
                    flow.params.weight_exponent[c] = (flow.mode == FLOW_MODE_ORBITALS) ? e : 1.0f;
                }
            }

            md_flow_layout_compute(&flow.layout, &flow.graph, &flow.cut, &flow.params);

            // The layout lives in graph space and knows nothing about pan or zoom, so fitting
            // after it is computed costs nothing and needs no second pass.
            if (flow.fit_requested) {
                flow.fit_requested = false;
                flow_fit_view(canvas_sz);
            }

            FlowView view = {};
            view.area_min = canvas_p0;
            view.area_max = canvas_p1;
            view.pan  = { flow.pan.x, flow.pan.y };
            view.zoom = flow.zoom;

            // Node state is DERIVED from the application masks every frame rather than stored.
            // That is what makes a selection made in the 3D view, or by a script, light up here
            // without any event plumbing - and it is why there is no local copy to fall out of
            // step with the rest of VIAMD.
            md_temp_scope_t draw_temp = md_temp_begin();
            MD_DEFER { md_temp_end(draw_temp); };

            const size_t cut_space = md_array_size(flow.graph.nodes) + md_array_size(flow.cut.other);
            uint8_t* node_state = md_temp_alloc_zero_array(draw_temp, uint8_t, cut_space ? cut_space : 1);

            const size_t num_atoms = md_array_size(nto.atom_group_idx);

            // One pass over the atoms gives the per-atom answer and the tallies per group and per
            // fragment (summed up the fragment tree), so a group or fragment node is a lookup.
            // Atoms are QM atoms; the masks are over system atoms.
            uint32_t* group_total = md_temp_alloc_zero_array(draw_temp, uint32_t, MAX_NTO_GROUPS);
            uint32_t* group_sel   = md_temp_alloc_zero_array(draw_temp, uint32_t, MAX_NTO_GROUPS);
            uint32_t* group_hi    = md_temp_alloc_zero_array(draw_temp, uint32_t, MAX_NTO_GROUPS);
            const size_t num_frags = flow.frags.count;
            uint32_t* frag_total = md_temp_alloc_zero_array(draw_temp, uint32_t, num_frags ? num_frags : 1);
            uint32_t* frag_sel   = md_temp_alloc_zero_array(draw_temp, uint32_t, num_frags ? num_frags : 1);
            uint32_t* frag_hi    = md_temp_alloc_zero_array(draw_temp, uint32_t, num_frags ? num_frags : 1);
            bool* atom_sel = md_temp_alloc_zero_array(draw_temp, bool, num_atoms ? num_atoms : 1);
            bool* atom_hi  = md_temp_alloc_zero_array(draw_temp, bool, num_atoms ? num_atoms : 1);
            bool selection_touches_atoms = false;

            for (size_t a = 0; a < num_atoms; ++a) {
                const uint32_t g = nto.atom_group_idx[a] < MAX_NTO_GROUPS ? nto.atom_group_idx[a] : 0;
                const uint32_t s = flow_qm_to_sys((uint32_t)a);
                const bool sel = md_bitmask_test_bit(&state.selection.selection_mask, s);
                const bool hi  = md_bitmask_test_bit(&state.selection.highlight_mask, s);
                atom_sel[a] = sel;
                atom_hi[a]  = hi;
                group_total[g] += 1;
                group_sel[g]   += sel ? 1 : 0;
                group_hi[g]    += hi  ? 1 : 0;
                selection_touches_atoms |= sel;
                if (g == 0) {
                    for (uint32_t f = flow_fragment_of((uint32_t)a); f != MD_FRAGMENT_NONE; f = flow.frags.parent[f]) {
                        frag_total[f] += 1;
                        frag_sel[f]   += sel ? 1 : 0;
                        frag_hi[f]    += hi  ? 1 : 0;
                    }
                }
            }

            // A set is only "selected" when all of it is. Anything less is partial - otherwise one
            // selected atom would claim the whole group, and there would be no way to see the difference.
            auto set_state = [](uint32_t total, uint32_t sel, uint32_t hi) -> uint8_t {
                uint8_t st = 0;
                if (total == 0) return st;
                if (sel == total) st |= FlowNodeState_Selected;
                else if (sel > 0) st |= FlowNodeState_Partial;
                if (hi > 0)       st |= FlowNodeState_Highlight;
                return st;
            };

            for (size_t i = 0; i < md_array_size(flow.graph.nodes); ++i) {
                const md_flow_node_t* n = flow.graph.nodes + i;
                uint8_t st = 0;
                const uint32_t index = flow_key_index(n->key);
                switch (n->level) {
                case FLOW_LEVEL_ATOM:
                    if (index < num_atoms) {
                        if (atom_sel[index]) st |= FlowNodeState_Selected;
                        if (atom_hi[index])  st |= FlowNodeState_Highlight;
                    }
                    break;
                case FLOW_LEVEL_GROUP:
                    if (index < MAX_NTO_GROUPS) st = set_state(group_total[index], group_sel[index], group_hi[index]);
                    break;
                case FLOW_LEVEL_FRAGMENT:
                    if (index < num_frags) st = set_state(frag_total[index], frag_sel[index], frag_hi[index]);
                    break;
                default:    // orbitals
                    if (flow_orbital_is_selected(n->key)) st |= FlowNodeState_Selected;
                    break;
                }
                node_state[i] = st;
            }

            FlowEmphasis emphasis = {};
            emphasis.node_state = node_state;
            emphasis.num_state_entries = cut_space;
            emphasis.any_selected = selection_touches_atoms || md_array_size(flow.selected_orbital_keys) > 0;
            emphasis.hover_node = flow.hover.node;
            emphasis.hover_link = flow.hover.link;
            if (side_w > 0.0f && flow.hover.node < 0 && flow.hover.link < 0) {
                emphasis.hover_link = ctv_pin_link();     // the pinned ribbon stands out while nothing is hovered
            }

            flow_update_stripes();
            FlowStripes stripes = {};
            stripes.offset    = flow.stripe_offset;
            stripes.frac      = flow.stripe_frac;
            stripes.num_links = md_array_size(flow.stripe_offset) ? md_array_size(flow.stripe_offset) - 1 : 0;

            FlowAnnotations notes = {};
            char subtitle[96];
            flow_annotations(&notes, subtitle, sizeof(subtitle));

            draw_list->PushClipRect(canvas_p0, canvas_p1, true);
            const FlowHit hit = flow_draw_imgui(draw_list, view, &flow.graph, &flow.cut, &flow.layout,
                                                flow.style, emphasis, &stripes, &notes);
            draw_list->PopClipRect();

            flow.hover = canvas_hovered ? hit : FlowHit{};

            // With the side views on: a click on a ribbon (or, in the orbital route, an NTO pair) pins its orbitals
            // in them, again unpins; a click on empty space or Esc lets go. A press and release in place.
            if (side_w > 0.0f) {
                if (canvas_hovered && hit.toggle < 0 && !io.KeyShift && !flow.region_active &&
                    ImGui::IsMouseReleased(ImGuiMouseButton_Left) && io.MouseDragMaxDistanceSqr[ImGuiMouseButton_Left] < 16.0f) {
                    const CtOrbReq r = ctv_request(hit.link, hit.stripe, hit.node);
                    if (r.key) ctv_pin_toggle(r);
                    else if (hit.node < 0 && hit.link < 0) ctv.pin = {};
                }
                if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) && ImGui::IsKeyPressed(ImGuiKey_Escape, false)) ctv.pin = {};
            }

            // A group's arrow: one click opens or closes it, pinned (the row budget leaves it so) and on both
            // sides. A click is a press and release in place, as everywhere here, since a drag pans.
            if (canvas_hovered && hit.toggle >= 0 && !io.KeyShift && !flow.region_active) {
                ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left) && io.MouseDragMaxDistanceSqr[ImGuiMouseButton_Left] < 16.0f) {
                    const uint32_t idx = (uint32_t)hit.toggle;
                    md_flow_cut_pin(&flow.cut, &flow.graph, idx, !md_flow_cut_is_open(&flow.cut, idx));
                    md_flow_cut_resolve(&flow.cut, &flow.graph);
                }
            }
            if (canvas_hovered && hit.subtitle && flow.mode == FLOW_MODE_FRAGMENTS) {
                ImGui::SetTooltip("Between the nodes as they are open now: a fragment's share that stays on\n"
                                  "it is local excitation (LE), the rest is charge transfer (CT). Opening a\n"
                                  "fragment into atoms turns most of its local part into transfer between its atoms.%s%s",
                                  flow.nto_dropped  > 1.0e-3f ? "\n\nNTO pairs below round off are left out." : "",
                                  flow.link_dropped > 1.0e-3f ? "\n\nAtom pairs below 1e-6 of the total are not drawn." : "");
            }

            // ---- selection: rubber band, or a click when the band never opened ----
            //
            // Press, drag and release are one gesture, so they are handled in one place rather than
            // split across the "is something under the cursor" branches - a drag that starts on a
            // node routinely ends over empty space, and the release is the half that matters.
            bool region_preview = false;
            {
                // Every button here also pans, so a "click" means pressed and released without
                // travelling. Deciding that at RELEASE, from the distance actually covered, is what
                // lets one gesture be either a click or a band without the user declaring which.
                const float slop = 4.0f;

                if (canvas_hovered && io.KeyShift && !flow.region_active &&
                    (ImGui::IsMouseClicked(ImGuiMouseButton_Left) || ImGui::IsMouseClicked(ImGuiMouseButton_Right)))
                {
                    flow.region_active = true;
                    flow.region_add    = ImGui::IsMouseClicked(ImGuiMouseButton_Left);
                    flow.region_start  = io.MousePos;
                }

                if (flow.region_active) {
                    const ImGuiMouseButton button = flow.region_add ? ImGuiMouseButton_Left : ImGuiMouseButton_Right;
                    const ImVec2 rect_min = { ImMin(flow.region_start.x, io.MousePos.x), ImMin(flow.region_start.y, io.MousePos.y) };
                    const ImVec2 rect_max = { ImMax(flow.region_start.x, io.MousePos.x), ImMax(flow.region_start.y, io.MousePos.y) };
                    const bool dragged = ImLengthSqr(io.MousePos - flow.region_start) >= slop * slop;

                    md_temp_scope_t rtemp = md_temp_begin();
                    MD_DEFER { md_temp_end(rtemp); };
                    md_allocator_i* ralloc = md_temp_allocator(rtemp);

                    if (ImGui::IsMouseDown(button)) {
                        if (dragged) {
                            region_preview = true;
                            const ImU32 fill = flow.region_add ? IM_COL32(60, 130, 220, 48) : IM_COL32(210, 70, 70, 48);
                            const ImU32 line = flow.region_add ? IM_COL32(60, 130, 220, 200) : IM_COL32(210, 70, 70, 200);
                            draw_list->PushClipRect(canvas_p0, canvas_p1, true);
                            draw_list->AddRectFilled(rect_min, rect_max, fill);
                            draw_list->AddRect(rect_min, rect_max, line);
                            draw_list->PopClipRect();

                            // Preview into the highlight mask so the 3D view shows what the release
                            // is about to do, for both add and remove.
                            md_array(uint32_t) nodes = 0;
                            flow_nodes_in_region(&nodes, view, rect_min, rect_max, ralloc);

                            md_bitmask_clear(&state.selection.highlight_mask);
                            for (size_t i = 0; i < md_array_size(nodes); ++i) {
                                md_bitmask_or(&state.selection.highlight_mask, flow_node_atom_mask(nodes[i], ralloc));
                            }
                        }
                    } else {
                        if (dragged) {
                            md_array(uint32_t) nodes = 0;
                            flow_nodes_in_region(&nodes, view, rect_min, rect_max, ralloc);
                            for (size_t i = 0; i < md_array_size(nodes); ++i) {
                                flow_set_selected(state, nodes[i], flow.region_add);
                            }
                        } else if (hit.node >= 0) {
                            // Never opened into a band: the same gesture on a single item.
                            flow_set_selected(state, (uint32_t)hit.node, flow.region_add);
                        } else if (!flow.region_add) {
                            // Shift+right on nothing clears everything - the gesture that removes
                            // one item, aimed at no item.
                            flow_clear_selection(state);
                        }
                        flow.region_active = false;
                    }
                }
            }

            // Plain right click (no shift, no drag) opens the group menu. Shift+right is the
            // deselect gesture and a right DRAG pans, so this is the one right-button meaning left
            // over - and it is the one that needs somewhere to live, since the new window had no
            // way to make a group at all.
            if (canvas_hovered && !io.KeyShift && !flow.region_active &&
                ImGui::IsMouseReleased(ImGuiMouseButton_Right) &&
                io.MouseDragMaxDistanceSqr[ImGuiMouseButton_Right] < 16.0f)
            {
                ImGui::OpenPopup("##flow_context");
            }

            if (ImGui::BeginPopup("##flow_context")) {
                const size_t num_selected = md_bitmask_count(&state.selection.selection_mask);
                if (num_selected == 0) {
                    ImGui::TextDisabled("No atoms selected");
                    ImGui::Separator();
                    ImGui::TextDisabled("Shift+click or shift+drag to select");
                } else {
                    ImGui::Text("%zu atom%s selected", num_selected, num_selected == 1 ? "" : "s");
                    ImGui::Separator();

                    if (ImGui::MenuItem("Create new group", nullptr, false, nto.group.count < MAX_NTO_GROUPS)) {
                        flow_create_group_from_selection(state);
                    }
                    if (ImGui::BeginMenu("Add to group", nto.group.count > 1)) {
                        // Slot 0 is the ungrouped bucket, not a destination anyone means to pick.
                        for (size_t g = 1; g < nto.group.count; ++g) {
                            if (ImGui::MenuItem(nto.group.label[g])) {
                                flow_assign_selection_to_group(state, (uint32_t)g);
                            }
                        }
                        ImGui::EndMenu();
                    }
                    if (ImGui::MenuItem("Remove from groups")) {
                        flow_assign_selection_to_group(state, 0);
                    }
                    ImGui::Separator();
                    if (ImGui::MenuItem("Clear selection")) {
                        flow_clear_selection(state);
                    }
                }
                ImGui::EndPopup();
            }

            if (canvas_hovered && hit.node >= 0) {
                const md_flow_node_t* node = md_flow_cut_node(&flow.cut, &flow.graph, (uint32_t)hit.node);
                if (node && (node->flags & MD_FLOW_NODE_FLAG_OTHER)) {
                    // An "Others" is a sum, so its tooltip has to say what it is a sum OF - an
                    // unexplained grey box holding a chunk of the transfer is worse than no box.
                    const uint32_t other_idx = (uint32_t)hit.node - (uint32_t)md_array_size(flow.graph.nodes);
                    const uint32_t count = other_idx < md_array_size(flow.cut.other_count)
                                         ? flow.cut.other_count[other_idx] : 0;
                    ImGui::SetTooltip("Others\n%.1f%% across %u entries, each below %.1f%%\n\n"
                                      "Double click to draw them individually",
                                      node->weight * 100.0f, count, flow.cut.threshold * 100.0f);
                } else if (node) {
                    const bool two_columns = flow.graph.num_columns == 2;
                    const char* what = (node->column == 0) ? (two_columns ? "hole" : "detachment")
                                     : (node->column == 1 && !two_columns) ? "transition"
                                     : (two_columns ? "particle" : "attachment");
                    const bool has_children = (uint32_t)hit.node < md_array_size(flow.graph.nodes)
                                              && !md_flow_node_is_leaf(&flow.graph, (uint32_t)hit.node);
                    // Of a fragment's share, what stays on it (drawn inside the node) and what moves
                    char local[96] = "";
                    if (flow.mode == FLOW_MODE_FRAGMENTS && node->weight > 0.0f) {
                        const float self = flow_node_self_weight((uint32_t)hit.node);
                        snprintf(local, sizeof(local), "\n%.1f%% of it local, %.1f%% %s",
                                 100.0f * self / node->weight, 100.0f * (1.0f - self / node->weight),
                                 node->column == 0 ? "transferred away" : "transferred in");
                    }
                    if (has_children) {
                        const bool open = md_flow_cut_is_open(&flow.cut, (uint32_t)hit.node);
                        // a group with an arrow (its band) says so; a fragment of atoms has none
                        bool has_arrow = false;
                        for (size_t b = 0; b < md_array_size(flow.layout.bands); ++b) has_arrow |= flow.layout.bands[b].node == (uint32_t)hit.node;
                        ImGui::SetTooltip("%.*s\n%.1f%% of the %s%s\n\n%s to %s",
                            (int)node->label.len, node->label.ptr, node->weight * 100.0f, what, local,
                            has_arrow ? "Click the arrow, or double click," : "Double click",
                            open ? "collapse" : "open it");
                    } else {
                        ImGui::SetTooltip("%.*s\n%.1f%% of the %s%s",
                            (int)node->label.len, node->label.ptr, node->weight * 100.0f, what, local);
                    }
                }
                // Hovering previews in the viewport; committing to a selection takes shift, so a
                // plain drag can pan without the selection changing under the cursor. While a band
                // is open it owns the preview, or the node under the cursor would keep overwriting
                // what the band is about to do.
                if (!region_preview) {
                    md_bitmask_assign(&state.selection.highlight_mask, flow_node_atom_mask((uint32_t)hit.node, state.allocator.frame));
                }

                // Double click opens or closes a node - including an "Others", which opens by drawing
                // its members individually. That is what lets the threshold hide things honestly:
                // everything below it goes, and you can still look. What the user opens or closes is
                // PINNED, so the row budget leaves it as it is; a fragment is pinned on both sides.
                // Not on an arrow: its two clicks have already closed and opened it again
                if (!io.KeyShift && hit.toggle < 0 && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    const uint32_t idx = (uint32_t)hit.node;
                    const md_flow_node_t* n = md_flow_cut_node(&flow.cut, &flow.graph, idx);
                    if (n && (n->flags & MD_FLOW_NODE_FLAG_OTHER)) {
                        md_flow_cut_open_other(&flow.cut, &flow.graph, n, true);
                        md_flow_cut_resolve(&flow.cut, &flow.graph);
                    } else if (idx < md_array_size(flow.graph.nodes) && !md_flow_node_is_leaf(&flow.graph, idx)) {
                        md_flow_cut_pin(&flow.cut, &flow.graph, idx, !md_flow_cut_is_open(&flow.cut, idx));
                        md_flow_cut_resolve(&flow.cut, &flow.graph);
                    } else if (n && md_flow_cut_is_other_open_for(&flow.cut, n)) {
                        // Double clicking a member of an opened "Others" folds it again - the pair of
                        // gestures is symmetric, since the "Others" itself is gone while open.
                        const md_flow_node_t other = { .column = n->column, .parent = n->parent };
                        md_flow_cut_open_other(&flow.cut, &flow.graph, &other, false);
                        md_flow_cut_set_other_expanded(&flow.cut, n->column, false);
                        md_flow_cut_resolve(&flow.cut, &flow.graph);
                    }
                }
            } else if (canvas_hovered && hit.node < 0 && hit.link < 0) {
                // Cursor is in this canvas and over nothing, so whatever this window last put in
                // the highlight mask is stale. Cleared only under those two conditions: clearing it
                // unconditionally would stamp on the 3D viewport's own hover.
                if (!region_preview) {
                    md_bitmask_clear(&state.selection.highlight_mask);
                }

                // Empty space: double click reframes. A button for this would occupy toolbar room
                // permanently to undo something that only occasionally needs undoing.
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    flow_fit_view(canvas_sz);
                }
            }
            if (canvas_hovered && hit.link >= 0) {
                const md_flow_link_t* link = flow.cut.links + hit.link;
                const md_flow_node_t* src = md_flow_cut_node(&flow.cut, &flow.graph, link->src);
                const md_flow_node_t* dst = md_flow_cut_node(&flow.cut, &flow.graph, link->dst);
                if (src && dst) {
                    // What it is made of: the orbital pairs of this transfer, when there are several
                    char pairs[160] = "";
                    const size_t num_links = md_array_size(flow.stripe_offset) ? md_array_size(flow.stripe_offset) - 1 : 0;
                    if ((size_t)hit.link < num_links) {
                        const uint32_t beg = flow.stripe_offset[hit.link], end = flow.stripe_offset[hit.link + 1];
                        if (end - beg >= 2) {
                            int len = snprintf(pairs, sizeof(pairs), "\nCarried by %u orbital pairs:", end - beg);
                            for (uint32_t j = beg; j < end && len > 0 && (size_t)len < sizeof(pairs); ++j) {
                                len += snprintf(pairs + len, sizeof(pairs) - (size_t)len, " %.0f%%", 100.0f * flow.stripe_frac[j]);
                            }
                            if (hit.stripe >= 0 && (uint32_t)hit.stripe < end - beg && len > 0 && (size_t)len < sizeof(pairs)) {
                                snprintf(pairs + len, sizeof(pairs) - (size_t)len, "\n(pair %d under the cursor)", hit.stripe + 1);
                            }
                        }
                    }
                    ImGui::SetTooltip("%.*s  ->  %.*s\n%.2f%%%s",
                        (int)src->label.len, src->label.ptr,
                        (int)dst->label.len, dst->label.ptr,
                        link->weight * 100.0f, pairs);
                }
            }
        }
        ImGui::End();

        if (export_requested) {
            export_flow_diagram(ImVec2(1024, 640));
        }
    }

    bool export_flow_diagram(ImVec2 size) {
        if (!flow.built) {
            MD_LOG_ERROR("Charge transfer export failed: no diagram to export");
            return false;
        }

        char path_buf[1024] = {};
        if (!application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, MD_STR_LIT("svg"))) {
            return false;
        }

        md_temp_scope_t temp = md_temp_begin();
        MD_DEFER { md_temp_end(temp); };

        md_vg_scene_t scene = {};
        md_vg_scene_init(&scene, { size.x, size.y }, md_temp_allocator(temp));

        // The export draws the SAME layout the screen does, at a different size. That is the
        // point of keeping layout separate from drawing: there is no second implementation to
        // drift out of agreement with the first.
        md_flow_layout_t layout = {};
        md_flow_layout_init(&layout, md_temp_allocator(temp));
        md_flow_layout_compute(&layout, &flow.graph, &flow.cut, &flow.params);

        FlowAnnotations notes = {};
        char subtitle[96];
        flow_annotations(&notes, subtitle, sizeof(subtitle));

        FlowView view = {};
        view.area_min = { 0.0f, 0.0f };
        view.area_max = size;
        vec2_t pan = {};
        flow_fit(&pan, &view.zoom, layout, size, flow_measure_margins_vg(&flow.graph, &flow.cut, &layout, flow.style, &notes));
        view.pan = { pan.x, pan.y };

        // the page under it, so the file reads the same in any viewer's theme
        md_vg_add_rect_filled(&scene, { 0.0f, 0.0f }, { size.x, size.y }, MD_VG_COLOR_RGB(255, 255, 255), 0.0f);

        FlowStripes stripes = {};
        stripes.offset    = flow.stripe_offset;
        stripes.frac      = flow.stripe_frac;
        stripes.num_links = md_array_size(flow.stripe_offset) ? md_array_size(flow.stripe_offset) - 1 : 0;
        flow_draw_vg(&scene, view, &flow.graph, &flow.cut, &layout, flow.style, &stripes, &notes);

        const bool success = md_vg_scene_write_svg_file(&scene, md_str_from_cstr(path_buf));
        if (success) {
            MD_LOG_INFO("Exported charge transfer diagram to '%s'", path_buf);
        } else {
            MD_LOG_ERROR("Failed to export charge transfer diagram to '%s'", path_buf);
        }
        return success;
    }

    // Bookkeeping rather than numerics: it owns the hashes and the transition matrix storage, and
    // the two procedures it hands the numbers to - attribute_charge_density and
    // compute_transition_matrix - already take nothing but arrays. So this is where the system is
    // allowed to be, and the fetches stay INSIDE the hash guard: a transition density is rebuilt
    // from the solution vector on every read, and pulling two A x A of them per frame would be a
    // real cost for a value that only changes when the selection does.
    void update_nto_derived_data(const md_system_t& sys) {
        if (!nto.atom_group_idx) return;

        // What the table carries, not what a reader parsed. A system with no transition densities
        // has no transition to attribute, whichever program produced it.
        if (!md_attributes_find(&sys.attributes, es_path::attachment_density)) return;

        // Create hash to check for changes to trigger update of colors in gl representation
        uint64_t atom_idx_hash = md_hash64(nto.atom_group_idx, md_array_bytes(nto.atom_group_idx), 0);
        uint64_t color_hash    = md_hash64(nto.group.color, sizeof(nto.group.color), atom_idx_hash);

        static uint64_t cur_color_hash = 0;
        if (color_hash != cur_color_hash) {
            cur_color_hash = color_hash;
            update_nto_group_colors(sys);
        }

        // Create hash to check for changes to trigger recomputation of transition matrix
        uint64_t matrix_hash = atom_idx_hash ^ (uint64_t)rsp.selected ^ nto.group.count;
        static uint64_t cur_matrix_hash = 0;

        if (matrix_hash != cur_matrix_hash) {
            cur_matrix_hash = matrix_hash;

            // Resize transition matrix to the correct size
            if (nto.transition_matrix_dim != nto.group.count) {
                if (nto.transition_matrix) {
                    // The allocated size contains matrix N*N + 2*N for hole/part arrays
                    const size_t cur_mem = sizeof(float) * nto.transition_matrix_dim * (nto.transition_matrix_dim + 2);
                    md_free(arena, nto.transition_matrix, cur_mem);
                }

                nto.transition_matrix_dim = nto.group.count;
                const size_t new_mem = sizeof(float) * nto.transition_matrix_dim * (nto.transition_matrix_dim + 2);
                nto.transition_matrix = (float*)md_alloc(arena, new_mem);
                nto.transition_density_hole = nto.transition_matrix + (nto.transition_matrix_dim * nto.transition_matrix_dim);
                nto.transition_density_part = nto.transition_density_hole + nto.transition_matrix_dim;
            }

            // Clear all of the values to zero (matrix, group values)
            md_memset(nto.transition_matrix, 0, sizeof(float) * nto.transition_matrix_dim * (nto.transition_matrix_dim + 2));

			if (rsp.selected != -1) {
				const size_t nto_idx = (size_t)rsp.selected;
#if 1
				{
					md_temp_scope_t temp = md_temp_begin();
                    MD_DEFER { md_temp_end(temp); };

					// Full-density Mulliken attribution from attachment / detachment AO matrices.
					// hole charges  <- detachment density (D-)
					// part charges  <- attachment density (D+)
					//
					// Both are VIRTUAL attributes indexed by excited state, and the Transition
					// Analysis window wants the same two matrices in this same frame - so they
					// come from the component's cache, which rebuilds them once per state.
					const bool have_td = td_cache_ensure(sys, (int)nto_idx);

					size_t s_dim = 0;
					const double* D_attach = have_td ? td_cache.attach : nullptr;
					const double* D_detach = have_td ? td_cache.detach : nullptr;
					// The overlap is STORED, so extracting it is a copy and not a reconstruction.
					const double* S        = density_matrix_extract(&s_dim, temp, sys, es_path::overlap, nullptr);

					const size_t num_aos = s_dim;

					if (num_aos > 0 && S && D_attach && D_detach) {
						if (td_cache.dim == num_aos) {
							double group_density_part[MAX_NTO_GROUPS] = {0};
							double group_density_hole[MAX_NTO_GROUPS] = {0};

							// A pure function of the shell list - the same call md_vlx uses to build
							// the copy this used to read - so there is no map to publish and none
							// that can drift from the AO ordering.
							uint32_t* ao_to_atom = md_temp_alloc_array(temp, uint32_t, num_aos);
							if (md_gto_basis_num_ao(&basis) != num_aos ||
								md_gto_basis_ao_to_atom(ao_to_atom, &basis) != num_aos) {
								return;
							}

							// Bounded by the array actually being indexed rather than by any atom
							// count reported elsewhere, which is the only bound that makes the read
							// below safe.
							const size_t num_group_atoms = md_array_size(nto.atom_group_idx);
                            int* ao_to_group = md_temp_alloc_array(temp, int, num_aos);

							for (size_t i = 0; i < num_aos; i++) {
								const uint32_t atom_idx = ao_to_atom[i];
								MD_ASSERT(atom_idx < num_group_atoms);
								int g = (atom_idx < num_group_atoms) ? (int)nto.atom_group_idx[atom_idx] : 0;
								if (g < 0 || g >= MAX_NTO_GROUPS) g = 0;
								ao_to_group[i] = g;
							}

							attribute_charge_density(group_density_part, ao_to_group, D_attach, S, num_aos);
							attribute_charge_density(group_density_hole, ao_to_group, D_detach, S, num_aos);

							for (size_t g1 = 0; g1 < nto.group.count; g1++) {
								// Mulliken populations of the attachment / detachment density are non-negative
								// in exact arithmetic. Clamp tiny negative numerical noise to zero so that
								// the downstream transition matrix cannot acquire spurious negative entries.
								nto.transition_density_part[g1] = (float)MD_MAX(0.0, group_density_part[g1]);
								nto.transition_density_hole[g1] = (float)MD_MAX(0.0, group_density_hole[g1]);
							}
							compute_transition_matrix(nto.transition_matrix, nto.group.count, nto.transition_density_hole, nto.transition_density_part);
						}
					}
				}

#else
                const float samples_per_unit_length = DEFAULT_SAMPLES_PER_ANGSTROM * BOHR_TO_ANGSTROM;

                if (use_gpu_path) {
                    md_gto_segment_and_attribute_to_groups_GPU(nto.transition_density_part, nto.group.count, nto.vol[NTO_Attachment].tex_id, &nto.grid, (const float*)nto.atom_xyzr, nto.atom_group_idx, nto.num_atoms);
                    md_gto_segment_and_attribute_to_groups_GPU(nto.transition_density_hole, nto.group.count, nto.vol[NTO_Detachment].tex_id, &nto.grid, (const float*)nto.atom_xyzr, nto.atom_group_idx, nto.num_atoms);
                    compute_transition_matrix(nto.transition_matrix, nto.group.count, nto.transition_density_hole, nto.transition_density_part);
                } else {
                    task_system::ID eval_attach = 0;
                    task_system::ID seg_attach  = 0;
                    task_system::ID eval_detach = 0;
                    task_system::ID seg_detach  = 0;

                    if (compute_transition_group_values_async(&eval_attach, &seg_attach, nto.transition_density_part, nto.group.count, nto.grid, nto.atom_group_idx, nto.atom_xyzr, nto.num_atoms, nto_idx, ElectronicStructureNtoComponent::Particle, MD_GTO_EVAL_MODE_PSI, samples_per_unit_length) &&
                        compute_transition_group_values_async(&eval_detach, &seg_detach, nto.transition_density_hole, nto.group.count, nto.grid, nto.atom_group_idx, nto.atom_xyzr, nto.num_atoms, nto_idx, ElectronicStructureNtoComponent::Hole,     MD_GTO_EVAL_MODE_PSI, samples_per_unit_length))
                    {
                        task_system::ID compute_matrix_task = task_system::create_main_task(MD_STR_LIT("##Compute Transition Matrix"), [nto = &nto]() {
                            compute_transition_matrix(nto->transition_matrix, nto->group.count, nto->transition_density_hole, nto->transition_density_part);
                            // Print out the groups and the values for debugging
#if 1
                            printf("Computed Transition Matrix for NTO index %d\n", nto->sel_nto_idx);
                            printf("Number of Groups: %zu\n", nto->group.count);
                            printf("Atom to Group Mapping:\n");
                            for (size_t i = 0; i < nto->num_atoms; i++) {
                                printf("    Atom %zu: Group %u\n", i, nto->atom_group_idx[i]);
                            }
                            printf("Transition Matrix:\n");
                            for (size_t i = 0; i < nto->group.count; i++) {
                                printf("  From Group %zu:\n", i);
                                for (size_t j = 0; j < nto->group.count; j++) {
                                    printf("    To Group %zu: %f\n", j, nto->transition_matrix[i * nto->group.count + j]);
                                }
                            }
                            printf("Group Densities (Particle):\n");
                            for (size_t i = 0; i < nto->group.count; i++) {
                                printf("    Group %zu: %f\n", i, nto->transition_density_part[i]);
                            }
                            printf("Group Densities (Hole):\n");
                            for (size_t i = 0; i < nto->group.count; i++) {
                                printf("    Group %zu: %f\n", i, nto->transition_density_hole[i]);
                            }
#endif          
                        });

                        task_system::set_task_dependency(compute_matrix_task, seg_attach);
                        task_system::set_task_dependency(compute_matrix_task, seg_detach);

                        task_system::enqueue_task(eval_attach);
                        task_system::enqueue_task(eval_detach);

                        nto.seg_task[0] = seg_attach;
                        nto.seg_task[1] = seg_detach;
                    } else {
                        MD_LOG_DEBUG("An error occured when computing nto group values");
                    }
                }
#endif
            }
        }
    }

};

static QuantumChemistry instance = {};
