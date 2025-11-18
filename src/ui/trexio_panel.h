#pragma once

#ifdef MD_TREXIO

namespace trexio_ui {

// Initialize the TREXIO UI panel
void initialize();

// Shutdown and cleanup resources
void shutdown();

// Update and render the TREXIO UI panel
void update();

// Show/hide the main tools window
void set_main_window_visible(bool visible);

} // namespace trexio_ui

#endif // MD_TREXIO
