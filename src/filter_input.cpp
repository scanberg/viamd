#include <filter_input.h>

#include <script_editor.h>
#include <imgui_widgets.h>

#include <core/md_allocator.h>
#include <core/md_bitset.h>
#include <core/md_str.h>

#include <imgui_internal.h>

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

bool filter_evaluate(AtomFilter* filter, md_bitset_t* atoms, const md_filter_context_t& ctx, md_allocator_i* alloc) {
    filter->result = md_filter(atoms, str_from_cstr(filter->text), &ctx, alloc);
    return filter->result.success;
}

bool filter_evaluate_sets(AtomFilter* filter, md_array(md_bitset_t)* sets, const md_filter_context_t& ctx, md_allocator_i* alloc) {
    filter->result = md_filter_sets(sets, str_from_cstr(filter->text), &ctx, alloc);
    return filter->result.success;
}

namespace {

struct Suggestion {
    std::string label;
    std::string text;       // what replaces the range of the text
    md_script_completion_kind_t kind;
    bool fits;              // of the type expected where it goes
};

// The suggestions of the one input they are open for: only one is typed in at a time
struct Completions {
    ImGuiID id = 0;                     // the input they belong to
    bool open = false;
    std::vector<Suggestion> items;
    int  selected = 0;
    bool scroll = false;                // bring the selected one into view
    int  range_beg = 0, range_end = 0;  // the bytes of the text a suggestion replaces

    // What the next callback of the input does
    int  take = -1;                     // takes this suggestion (Enter)
    bool ask = false;                   // looks for suggestions, also with nothing typed (Ctrl+Space)
    bool refocus = false;               // puts the cursor at 'cursor', once the input has the keyboard back (a click)
    int  cursor = 0;

    bool list_hovered = false;          // the list, the last time it was drawn
};

Completions comp;

const int MAX_SUGGESTIONS = 50;
const int MAX_ROWS = 10;

const char* kind_name(md_script_completion_kind_t kind) {
    switch (kind) {
    case MD_SCRIPT_COMPLETION_VALUE:     return "value";
    case MD_SCRIPT_COMPLETION_PARAMETER: return "parameter";
    case MD_SCRIPT_COMPLETION_VARIABLE:  return "variable";
    case MD_SCRIPT_COMPLETION_PROCEDURE: return "procedure";
    case MD_SCRIPT_COMPLETION_CONSTANT:  return "constant";
    case MD_SCRIPT_COMPLETION_KEYWORD:   return "keyword";
    default:                             return "";
    }
}

// What can be written at the cursor of text, matched against what has been typed of it and ranked. Shown when asked
// for, or while a name or a string is being typed; not when nothing matches or what has been typed is complete.
void suggest(const char* text, int len, int cursor, const md_filter_context_t& ctx, bool asked) {
    comp.items.clear();
    md_temp_scope_t temp = md_temp_begin();
    const md_script_completions_t c = md_filter_complete(str_t{text, (size_t)len}, cursor, &ctx, md_temp_allocator(temp));
    const std::string_view typed(c.prefix.ptr, c.prefix.len);

    struct Match {
        const md_script_completion_t* item;
        int score;
    };
    std::vector<Match> matches;
    for (size_t i = 0; i < c.count; ++i) {
        int score;
        if (script_editor::completion_match(typed, c.items[i], &score)) matches.push_back({&c.items[i], score});
    }
    std::sort(matches.begin(), matches.end(), [&](const Match& a, const Match& b) {
        return script_editor::completion_before(*a.item, a.score, *b.item, b.score, typed);
    });
    for (size_t i = 0; i < matches.size() && i < (size_t)MAX_SUGGESTIONS; ++i) {
        const md_script_completion_t& item = *matches[i].item;
        comp.items.push_back({std::string(item.label.ptr, item.label.len), std::string(item.text.ptr, item.text.len), item.kind, item.fits});
    }
    comp.range_beg = c.range.beg;
    comp.range_end = c.range.end;

    const std::string_view current(text + c.range.beg, (size_t)(c.range.end - c.range.beg));
    const bool complete = comp.items.size() == 1 && comp.items[0].text == current;
    const bool in_string = c.range.beg > 0 && (text[c.range.beg - 1] == '"' || text[c.range.beg - 1] == '\'');
    comp.open = !comp.items.empty() && !complete && (asked || !typed.empty() || in_string);
    comp.selected = 0;
    comp.scroll = true;
    md_temp_end(temp);
}

// Puts a suggestion in place of the range in the input's buffer
void take(ImGuiInputTextCallbackData* data, int index) {
    const Suggestion& s = comp.items[index];
    const int beg = ImClamp(comp.range_beg, 0, data->BufTextLen);
    const int end = ImClamp(comp.range_end, beg, data->BufTextLen);
    data->DeleteChars(beg, end - beg);
    data->InsertChars(beg, s.text.c_str());   // as much as there is room for
    data->CursorPos = data->SelectionStart = data->SelectionEnd = ImMin(beg + (int)s.text.size(), data->BufTextLen);
    comp.open = false;
}

struct CallbackData {
    const md_filter_context_t* ctx;
};

// One event per frame: Tab, Up or Down, an edit, or else the one of every frame
int callback(ImGuiInputTextCallbackData* data) {
    const CallbackData* cb = (const CallbackData*)data->UserData;
    switch (data->EventFlag) {
    case ImGuiInputTextFlags_CallbackAlways:
        if (comp.refocus) {
            comp.refocus = false;
            data->CursorPos = data->SelectionStart = data->SelectionEnd = ImClamp(comp.cursor, 0, data->BufTextLen);
        }
        if (comp.take >= 0) {
            if (comp.open && comp.take < (int)comp.items.size()) take(data, comp.take);
            comp.take = -1;
        }
        if (comp.ask) {
            comp.ask = false;
            suggest(data->Buf, data->BufTextLen, data->CursorPos, *cb->ctx, true);
        }
        // The cursor left what the suggestions are for
        if (comp.open && (data->CursorPos < comp.range_beg || data->CursorPos > comp.range_end)) {
            comp.open = false;
        }
        break;
    case ImGuiInputTextFlags_CallbackEdit:
        suggest(data->Buf, data->BufTextLen, data->CursorPos, *cb->ctx, false);
        break;
    case ImGuiInputTextFlags_CallbackCompletion:
        if (comp.open && !comp.items.empty()) {
            take(data, comp.selected);
        } else {
            suggest(data->Buf, data->BufTextLen, data->CursorPos, *cb->ctx, true);
        }
        break;
    case ImGuiInputTextFlags_CallbackHistory:
        if (comp.open && !comp.items.empty()) {
            const int n = (int)comp.items.size();
            comp.selected = (comp.selected + (data->EventKey == ImGuiKey_DownArrow ? 1 : n - 1)) % n;
            comp.scroll = true;
        }
        break;
    default:
        break;
    }
    return 0;
}

// Where byte i of text is drawn in the input that was just submitted
float text_x(ImGuiID id, const char* text, int i) {
    const ImGuiContext& g = *ImGui::GetCurrentContext();
    float scroll = 0.0f;
    if (g.ActiveId == id) {
        if (const ImGuiInputTextState* state = ImGui::GetInputTextState(id)) scroll = state->Scroll.x;
    }
    return ImGui::GetItemRectMin().x + g.Style.FramePadding.x + ImGui::CalcTextSize(text, text + i).x - scroll;
}

// The list under the input, from where the range begins. Returns the suggestion clicked, -1 if none.
int draw_list(ImGuiID id, const char* text) {
    const ImGuiStyle& style = ImGui::GetStyle();
    const ImVec2 item_min = ImGui::GetItemRectMin();
    const ImVec2 item_max = ImGui::GetItemRectMax();
    const float x = ImClamp(text_x(id, text, comp.range_beg) - style.WindowPadding.x, item_min.x, item_max.x);

    float label_w = 0.0f, kind_w = 0.0f;
    for (const Suggestion& s : comp.items) {
        label_w = ImMax(label_w, ImGui::CalcTextSize(s.label.c_str()).x);
        kind_w  = ImMax(kind_w,  ImGui::CalcTextSize(kind_name(s.kind)).x);
    }
    const int   num  = (int)comp.items.size();
    const float gap  = style.ItemSpacing.x * 3.0f;
    const float row_w = label_w + gap + kind_w;
    const float rows_h = ImMin(num, MAX_ROWS) * ImGui::GetTextLineHeightWithSpacing();
    const float bar_w = num > MAX_ROWS ? style.ScrollbarSize : 0.0f;

    int clicked = -1;
    ImGui::SetNextWindowPos(ImVec2(x, item_max.y));
    const ImGuiWindowFlags flags = ImGuiWindowFlags_Tooltip | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav;
    if (ImGui::Begin("##filter_suggestions", nullptr, flags)) {
        ImGui::BeginChild("##rows", ImVec2(row_w + bar_w, rows_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoNav);
        for (int i = 0; i < num; ++i) {
            const Suggestion& s = comp.items[i];
            ImGui::PushID(i);
            const float y = ImGui::GetCursorPosY();
            // What does not fit where it goes (in an argument) is there, but faint
            if (!s.fits) ImGui::PushStyleColor(ImGuiCol_Text, style.Colors[ImGuiCol_TextDisabled]);
            if (ImGui::Selectable(s.label.c_str(), i == comp.selected, ImGuiSelectableFlags_None, ImVec2(row_w, 0))) {
                clicked = i;
            }
            if (!s.fits) ImGui::PopStyleColor();
            if (i == comp.selected && comp.scroll) {
                ImGui::SetScrollHereY();
                comp.scroll = false;
            }
            ImGui::SetCursorPos(ImVec2(label_w + gap, y));
            ImGui::TextDisabled("%s", kind_name(s.kind));
            ImGui::PopID();
        }
        ImGui::EndChild();
        // Also while a suggestion is held down, which is an active item
        comp.list_hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows | ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
    }
    ImGui::End();
    return clicked;
}

// The part of the text an error is about, underlined in the input that was just submitted
void underline(ImGuiID id, const char* text, md_script_range_marker_t range) {
    const int len = (int)strlen(text);
    const int beg = ImClamp(range.beg, 0, len);
    const int end = ImClamp(range.end, beg, len);
    if (end <= beg) return;
    const ImGuiStyle& style = ImGui::GetStyle();
    const ImVec2 item_min = ImGui::GetItemRectMin();
    const ImVec2 item_max = ImGui::GetItemRectMax();
    const float x0 = ImMax(text_x(id, text, beg), item_min.x + style.FramePadding.x);
    const float x1 = ImMin(text_x(id, text, end), item_max.x - style.FramePadding.x);
    if (x1 <= x0) return;
    const float y = item_min.y + style.FramePadding.y + ImGui::GetFontSize() + 1.0f;
    ImGui::GetWindowDrawList()->AddLine(ImVec2(x0, y), ImVec2(x1, y), IM_COL32(255, 80, 80, 255), 2.0f);
}

}  // namespace

bool filter_completions_open() {
    // Those of an input that is not typed in any more (gone with its window) are not
    return comp.open && ImGui::GetActiveID() == comp.id;
}

bool filter_input(const char* label, AtomFilter* filter, const md_filter_context_t& ctx, ImGuiInputTextFlags flags, bool* confirmed, const char* hint) {
    const ImGuiID id = ImGui::GetID(label);
    const bool active = ImGui::GetActiveID() == id;
    if (active && comp.id != id) {
        // This input has the keyboard now, the suggestions of another go
        comp = Completions{};
        comp.id = id;
    }
    const bool mine = comp.id == id;

    // The suggestions of an input which no longer has the keyboard go, unless the mouse is on their way to them
    if (mine && comp.open && !active && !comp.list_hovered && !comp.refocus) {
        comp.open = false;
    }

    // Keys the list takes before the input sees them: Enter takes the selected suggestion, Escape closes the list
    if (mine && comp.open && active) {
        const ImGuiID owner = id + 1;
        if (ImGui::IsKeyPressed(ImGuiKey_Enter, false) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter, false)) {
            comp.take = comp.selected;
            ImGui::SetKeyOwner(ImGuiKey_Enter, owner, ImGuiInputFlags_LockThisFrame);
            ImGui::SetKeyOwner(ImGuiKey_KeypadEnter, owner, ImGuiInputFlags_LockThisFrame);
        }
        if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
            comp.open = false;
            ImGui::SetKeyOwner(ImGuiKey_Escape, owner, ImGuiInputFlags_LockThisFrame);
        }
    }
    if (active && ImGui::IsKeyChordPressed(ImGuiMod_Ctrl | ImGuiKey_Space)) {
        comp.ask = true;
    }
    if (mine && comp.refocus && !active) {
        ImGui::SetKeyboardFocusHere();
    }

    // Red until it evaluates, but not while it is empty
    const bool invalid = !filter->result.success && filter->text[0] != '\0';
    CallbackData data = {&ctx};
    const ImGuiInputTextFlags input_flags = (flags & ~ImGuiInputTextFlags_EnterReturnsTrue) |
        (confirmed ? ImGuiInputTextFlags_EnterReturnsTrue : 0) |
        ImGuiInputTextFlags_CallbackAlways | ImGuiInputTextFlags_CallbackEdit |
        ImGuiInputTextFlags_CallbackCompletion | ImGuiInputTextFlags_CallbackHistory;
    if (invalid) ImGui::PushInvalid();
    const bool entered = hint ? ImGui::InputTextWithHint(label, hint, filter->text, sizeof(filter->text), input_flags, callback, &data)
                              : ImGui::InputText(label, filter->text, sizeof(filter->text), input_flags, callback, &data);
    if (invalid) ImGui::PopInvalid();
    bool changed = ImGui::IsItemEdited();
    if (confirmed) {
        *confirmed = entered;
    } else {
        changed |= entered;
    }

    if (invalid) {
        underline(id, filter->text, filter->result.error_range);
        if (ImGui::IsItemHovered() && !(comp.open && comp.id == id) && filter->result.error[0]) {
            ImGui::SetTooltip("%s", filter->result.error);
        }
    }

    if (comp.id == id && comp.open) {
        const int clicked = draw_list(id, filter->text);
        if (clicked >= 0) {
            // The input lost the keyboard to the click: the suggestion goes into the text, and the keyboard back to it
            const Suggestion& s = comp.items[clicked];
            std::string text = filter->text;
            const int beg = ImClamp(comp.range_beg, 0, (int)text.size());
            const int end = ImClamp(comp.range_end, beg, (int)text.size());
            text.replace((size_t)beg, (size_t)(end - beg), s.text);
            str_copy_to_char_buf(filter->text, sizeof(filter->text), str_t{text.data(), text.size()});
            comp.open = false;
            comp.refocus = true;
            comp.cursor = ImMin(beg + (int)s.text.size(), (int)strlen(filter->text));
            changed = true;
        }
    } else if (comp.id == id) {
        comp.list_hovered = false;
    }

    return changed;
}
