#pragma once

#include <stddef.h>
#include <stdint.h>

#include <functional>
#include <vector>

// Text and a scale bar that are put on the frames of a movie, shown for a stretch of its timeline
enum class MovieOverlayType : int {
    Text,        // A title or a caption
    Timestamp,   // The time of the trajectory frame that is shown
    ScaleBar,    // A bar of a given length in the structure, which follows the camera
    Logo,        // The VIAMD logo
    Image,       // An image file (png or jpg)
    TimeBar,     // How far the trajectory has gone, counted forward even when it is played backward
    Timeline,    // Subplots of the Timelines window, stacked, drawn as the movie plays: wide and low, at the bottom
    Distribution,// Subplots of the Distributions window, stacked, filled as the movie plays: narrow and tall, at the side
    PropertyVis, // The visualization of a script property in the viewport, with its labels
    Figure,      // Only in the first workspaces with plots: both kinds in one. It is split into a Timeline and a Distribution when it is read.
    Count,
};

enum class MovieOverlayAnchor : int {
    TopLeft, TopCenter, TopRight,
    MiddleLeft, Center, MiddleRight,
    BottomLeft, BottomCenter, BottomRight,
    Count,
};

// What 'size' of an overlay is. A point is a pixel of a frame that is 1080 pixels high, scaled with the frame, so
// the same number looks the same at any resolution.
enum class MovieOverlaySizeUnit : int {
    FrameHeight,   // A part of the height of the frame
    Points,
    Count,
};

constexpr float MOVIE_OVERLAY_POINT_REFERENCE_HEIGHT = 1080.0f;

// Which window a panel of a figure takes its subplot from
enum class MoviePlotView : int {
    Timeline,
    Distribution,
    Count,
};

// One subplot of a figure, by the id the subplot has (not its position, which changes)
struct MoviePlotPanel {
    MoviePlotView view = MoviePlotView::Timeline;
    uint32_t      subplot = 0;
    double        begin = 0.0;   // Movie time (seconds) it appears at. Inside the overlay's own range: before it, the overlay is not shown.
    double        end = 0.0;     // ... and goes away at; not after 'begin' means it stays to the end of the overlay
    char          title[48] = "";   // Written above it with the overlay's Titles; empty takes the name of the subplot
};

// A note on the timeline of the movie, drawn on the timelines of a figure where the movie gets to it
struct MovieMarker {
    double time = 0.0;      // Movie time, in seconds
    char   label[48] = "";
    uint32_t subplot = 0;   // Stable timeline subplot id; 0 draws on all timeline subplots
    float  color[4] = {0.0f, 0.0f, 0.0f, 0.0f}; // RGBA; alpha 0 picks one from the palette by the marker's place in time
};

bool movie_marker_matches_subplot(const MovieMarker& marker, uint32_t subplot);

// The color of marker i: its own, or else the palette's, by its place in time among the markers (so that markers close in time
// differ and keep their color as the movie plays)
void movie_marker_color(const MovieMarker* markers, size_t n, size_t i, float out[4]);

// What the horizontal axis of a timeline overlay is
enum class MoviePlotAxis : int {
    Elapsed,           // Trajectory time that the movie has covered, like the time bar: the curve always grows to the right
    TrajectoryTime,    // Trajectory time itself, turned around when the movie plays the trajectory backward
    Count,
};

struct MovieOverlay {
    MovieOverlayType   type = MovieOverlayType::Text;
    bool               enabled = true;
    double             begin = 0.0;           // On the movie's timeline, in seconds
    double             end = 5.0;
    float              fade_in = 0.5f;        // Seconds to appear and to go away, within begin..end
    float              fade_out = 0.5f;
    MovieOverlayAnchor anchor = MovieOverlayAnchor::BottomLeft;
    float              size = 0.05f;          // The height of the text (of the logo), in size_unit
    MovieOverlaySizeUnit size_unit = MovieOverlaySizeUnit::FrameHeight;
    float              color[4] = {1, 1, 1, 1};
    float              background[4] = {0, 0, 0, 0};  // A plate behind it, none while its alpha is 0
    char               text[128] = "";
    float              length = 0.0f;         // Scale bar: its length in Angstrom, 0 chooses one
    char               path[512] = "";         // Image: the file
    float              width = 0.4f;          // Time bar: its width, as a part of the width of the frame
    bool               show_elapsed = true;   // Time bar: the time that has gone, over the whole
    bool               show_speed = false;    // Time bar: how fast the trajectory plays, relative to the Animation panel
    std::vector<MoviePlotPanel> panels;       // Timeline, Distribution: the subplots that are drawn, stacked in this order
    int                legacy_subplot_mask = 0;   // Old workspaces: bit i is the subplot at position i
    MoviePlotAxis      plot_axis = MoviePlotAxis::Elapsed;   // Timeline, Distribution: what the horizontal axis of its timelines is
    float              font_points = 0.0f;    // Timeline, Distribution: the text, in points (see MovieOverlaySizeUnit), 0 follows the height
    float              line_points = 0.0f;    // Timeline, Distribution: the width of the lines, in points, 0 follows the text
    int                palette = 0;           // Timeline, Distribution: 0 the colours of the plots, else a colour set of its own (movie_plot_palette_name)
    int                num_bins = 0;          // Distribution: 0 keeps each source series' bin count
    bool               show_markers = true;   // Timeline, Distribution: the markers of the movie, on its timelines
    bool               show_titles = false;   // Timeline, Distribution: the name of each subplot above it (when it has one)
    bool               reveal = true;         // Timeline, Distribution: only what the movie has played so far
    bool               show_value = true;     // Timeline, Distribution: the value at the frame that is shown, in the legend
};

constexpr int MOVIE_DISTRIBUTION_MIN_BINS = 2;
constexpr int MOVIE_DISTRIBUTION_MAX_BINS = 4096;
int movie_distribution_bins(const MovieOverlay& overlay, int source_bins);

// Names of the colour sets a figure can have of its own (index 0 is the colours of the plots themselves)
constexpr int MOVIE_PLOT_PALETTE_COUNT = 5;
const char* movie_plot_palette_name(int palette);

// Whether the horizontal axis is labelled under panel 'i' of a stack of 'n': under the last, under a distribution (which has
// an axis of its own) and under a timeline that has a distribution below it. Timelines above timelines share the axis.
bool movie_figure_x_labels(const MoviePlotView* views, size_t n, size_t i);

// What a new Timeline or Distribution overlay has for a place and a size: a timeline is wide and low at the bottom centre of the
// frame, a distribution narrow and tall at the middle of the right side. Other overlays are left.
void movie_overlay_plot_defaults(MovieOverlay* o);

// For workspaces from before the two kinds were separate: an overlay that holds a mask of the positions of subplots gets the
// ids of those subplots (given by position for each window), and a Figure is split into a Timeline and a Distribution overlay,
// each of the panels of its kind and with the place and size of its kind (everything else is kept).
void movie_overlays_migrate(std::vector<MovieOverlay>* overlays, const uint32_t* timeline_ids, int num_timeline, const uint32_t* distribution_ids, int num_distribution);

// The rectangle an overlay takes on the frame (with its plate), at its anchor, and the movie times it is shown
struct MovieOverlayBox {
    int    anchor = 0;   // MovieOverlayAnchor
    float  x0 = 0.0f, y0 = 0.0f, x1 = 0.0f, y1 = 0.0f;
    double begin = 0.0, end = 0.0;
};

// How far down (positive) or up each box is moved so that boxes shown at the same time do not overlap. Boxes in the top and
// bottom rows keep their place before those in the middle row, and in each the centre column before the others, then list
// order. A box that hits a placed one moves past it, 'gap' pixels clear: away from its edge of the frame in the top and
// bottom rows, the shorter way in the middle row. A box is never moved out of the frame [frame_y0, frame_y1]: it tries the
// other way, and stays where it is (overlapping) when neither fits.
void movie_overlay_avoid(const MovieOverlayBox* boxes, size_t n, float gap, float frame_y0, float frame_y1, float* dy_out);

// A label written on a line of text, from x0 to x1
struct MovieLabelSpan {
    float x0, x1;
    int row;
};

// Where the label (width w) of a mark at x goes, between lo and hi: right of the mark in the first half, left of it
// (ending at the mark) in the second half. It takes the first line where it does not come within 'gap' of a label
// placed before; then the other side of the mark; if every line is taken, where it overlaps least.
// Returns the line (0 is the top one) and writes where the label starts.
int movie_label_place(const MovieLabelSpan* placed, size_t n, float x, float w, float lo, float hi, float offset, float gap, int max_rows, float* x0_out);

// The logo in the top left corner for the whole movie, which a movie starts with
MovieOverlay movie_overlay_default_logo();

// The height of an overlay in pixels of a frame 'frame_height_px' high
float movie_overlay_size_px(const MovieOverlay& o, float frame_height_px);

// The same height written in another unit, e.g. to keep it as it is when the unit is changed
float movie_overlay_convert_size(float size, MovieOverlaySizeUnit from, MovieOverlaySizeUnit to);

// The range a size can be set in, in its unit
void movie_overlay_size_range(MovieOverlaySizeUnit unit, float* lo, float* hi);

// How visible an overlay is at a time, 0..1. Zero outside begin..end, rising over fade_in and falling over
// fade_out; where they do not both fit, they share the time in proportion.
float movie_overlay_alpha(const MovieOverlay& o, double time);

// A length in the structure that suits a bar a fraction 'target' of a frame 'span_px' wide, of the form 1, 2
// or 5 times a power of ten. units_per_pixel is how much of the structure one pixel covers.
float movie_scale_bar_length(double units_per_pixel, double span_px, double target);

// The rectangle that a frame of fw x fh takes in a viewport of vw x vh: as large as fits, in its proportions, centered, using
// the part 'fill' (0..1) of the room
void movie_frame_fit(float vw, float vh, float fw, float fh, float fill, float* x, float* y, float* w, float* h);

// The vertical field of view that makes a rectangle 'guide_h' high, in a viewport 'view_h' high, show what 'fov_y' shows in a whole frame
// (never wider than 170 degrees)
float movie_guide_fov_y(float fov_y, float view_h, float guide_h);

// How much of the structure one pixel of the frame covers, at the distance the camera looks at
double movie_units_per_pixel(float distance, float fov_y, float frame_height_px);

// How far a movie has taken the trajectory by each time, for the time bar. What is counted is the movement along the
// trajectory, whatever its direction: a trajectory that is played backward still fills the bar forward, and a stretch
// that is played fast fills it fast, a slow one slowly, a hold not at all.
struct MovieTimeBarProfile {
    std::vector<double> distance;   // Moved by each sample, in the unit of the quantity
    std::vector<double> q;          // The quantity at each sample
    std::vector<double> lo, hi;     // The least and the most the quantity has been at by each sample (the stretch of it that has been visited)
    double duration = 0.0;          // The movie's length, the samples are even over it
    double total() const { return distance.empty() ? 0.0 : distance.back(); }
};

// 'quantity_at(t)' is the trajectory time (or frame) at a movie time. Sampled 'samples' times over 0 .. duration.
void movie_time_bar_profile(MovieTimeBarProfile* out, double duration, int samples, const std::function<double(double)>& quantity_at);

// 0 (nothing moved yet) .. 1 (all moved) at a movie time. A trajectory that does not move at all gives 0.
double movie_time_bar_progress(const MovieTimeBarProfile& p, double time);

// The stretch of the quantity that the movie has visited by a movie time, from its least to its most value so far. At
// time 0 it is the single value the movie starts at. False while the profile is empty.
bool movie_time_bar_visited(const MovieTimeBarProfile& p, double time, double* lo, double* hi);

// Round tick positions in [lo, hi], 1, 2 or 5 times a power of ten apart, at most 'max_ticks' of them. Returns how many
// were written to 'out' (capacity 'cap') and the step between them.
int movie_nice_ticks(double lo, double hi, int max_ticks, double* out, int cap, double* step);

// Counts the values y[i * stride] * y_scale of the samples whose x lies in [x_lo, x_hi] into 'num_bins' bins evenly spread
// over [v_min, v_max]. Values outside it are left out. 'counts' is resized.
void movie_histogram_counts(std::vector<float>* counts, int num_bins, double v_min, double v_max, const float* x, const float* y,
    int stride, double y_scale, int num_samples, double x_lo, double x_hi);

// A series along the path of the movie: how much has moved (s) and the value there (v)
struct MovieCurvePoint { double s, v; };

// The points of a series along the elapsed axis from the start of the movie to a movie time, in order. Where the movie plays
// the trajectory backward the series is read backward, so s only ever grows. 'xs' are the ascending positions of the series'
// samples in the unit of the quantity, 'value_of_sample' the value of sample i and 'value_at' the value at any position.
void movie_elapsed_curve(std::vector<MovieCurvePoint>* out, const MovieTimeBarProfile& p, double time, const float* xs, int num_samples,
    const std::function<double(int)>& value_of_sample, const std::function<double(double)>& value_at);

// The movie times a panel is there from and to: its own range inside the range of the overlay
void movie_panel_span(const MoviePlotPanel& p, double overlay_begin, double overlay_end, double* begin, double* end);

// How visible a panel is at a time, 0..1, fading in and out over the overlay's fades at the ends of its span (zero outside it)
float movie_panel_alpha(const MoviePlotPanel& p, const MovieOverlay& overlay, double time);

// The stretch of the quantity visited from one movie time to another, from its least to its most value. False while the profile is empty.
bool movie_time_bar_visited_between(const MovieTimeBarProfile& p, double t0, double t1, double* lo, double* hi);

// The same as movie_elapsed_curve for the movie from t0 to t1: s counts from 0 at t0
void movie_elapsed_curve_between(std::vector<MovieCurvePoint>* out, const MovieTimeBarProfile& p, double t0, double t1, const float* xs, int num_samples,
    const std::function<double(int)>& value_of_sample, const std::function<double(double)>& value_at);

// How much has moved by a movie time, in the unit of the quantity
double movie_time_bar_moved(const MovieTimeBarProfile& p, double time);

// How fast the quantity changes at a movie time, per second of the movie, whatever the direction
double movie_quantity_speed(const std::function<double(double)>& quantity_at, double time, double step);
