#pragma once

// Filters as viamd keeps and edits them: the text of one and what came of it when it was last evaluated, and the input
// it is typed in.
//
// The input completes what is typed. What can be written at the cursor comes from mdlib (md_filter_complete), which
// knows that a filter picks atoms: outside of a call it offers only what does (and 'not'), after an operand only the
// operators that join it to another, and in a call everything, with what the parameter takes first. Suggestions are
// matched and ranked as in the script editor (script_editor::completion_match / completion_before). They pop up while
// a name or a string is typed, and Ctrl+Space asks for them; Tab or Enter takes the selected one, Up and Down move
// through them, Escape closes them, and a click takes one.
//
// A filter that does not evaluate turns its input red, the part of the text the error is about is underlined, and
// hovering it says why.

#include <md_filter.h>
#include <core/md_array.h>

#include <imgui.h>

struct md_allocator_i;
struct md_bitset_t;

struct AtomFilter {
    char text[256] = "";
    // Of the last evaluation: whether it picked atoms, why not and where, and whether what it picks depends on the
    // coordinates (and so has to be evaluated again in another frame)
    md_filter_result_t result = {};
};

// Evaluates the filter in ctx. The atoms it picks are written to atoms, whose payload is released and allocated with
// alloc (see md_filter). Returns whether it picked atoms.
bool filter_evaluate(AtomFilter* filter, md_bitset_t* atoms, const md_filter_context_t& ctx, md_allocator_i* alloc);
// The same, one set for each structure it gives (see md_filter_sets)
bool filter_evaluate_sets(AtomFilter* filter, md_array(md_bitset_t)* sets, const md_filter_context_t& ctx, md_allocator_i* alloc);

// The input of a filter. ctx is what it is completed against: its system, the script whose names it may use, and the
// atoms it is evaluated within, if any. Returns true when the text changes. With confirmed given, it is set when Enter
// is pressed to confirm the filter (and not to take a suggestion); ImGuiInputTextFlags_EnterReturnsTrue in flags is
// ignored either way. hint, if given, is shown while it is empty.
bool filter_input(const char* label, AtomFilter* filter, const md_filter_context_t& ctx, ImGuiInputTextFlags flags = 0, bool* confirmed = nullptr, const char* hint = nullptr);

// Whether the suggestions of a filter input are open: Enter and Escape are theirs then, not the window's
bool filter_completions_open();
