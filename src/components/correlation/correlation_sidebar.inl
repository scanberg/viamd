// CORRELATION: THE SIDEBAR
//
// The list of layers (grouped by their axes), the inspector of the selected one, and the picker of series for its axes.
// Members of the Correlation component: included within its definition in correlation.cpp, of which it sees everything.

    // Statistics of a layer as one line: its correlation, or how many pairs
    void layer_stats(char* buf, size_t cap, const Layer& l) const {
        if (!l.data_valid) {
            snprintf(buf, cap, "%s", l.data_error[0] ? l.data_error : "no samples");
            return;
        }
        const size_t n = md_array_size(l.sx);
        const char* what = l.binned ? "bins" : "samples";
        if (l.num_pairs == 1) {
            char st[96];
            stats_text(st, sizeof(st), l.pairs[0]);
            snprintf(buf, cap, "%s \xc2\xb7 %zu %s", st, n, what);
        } else {
            snprintf(buf, cap, "%d pairs \xc2\xb7 %zu %s", l.num_pairs, n, what);
        }
    }

    // The pairs of a layer, with their statistics
    void draw_pair_table(const Layer& l) const {
        constexpr int COLS = 6;
        if (!ImGui::BeginTable("##pairs", COLS, ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_RowBg)) return;
        ImGui::TableSetupColumn("pair", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("n");
        ImGui::TableSetupColumn("r");
        ImGui::TableSetupColumn("95 %");
        ImGui::TableSetupColumn("\xcf\x81");
        ImGui::TableSetupColumn("slope");
        ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
        static const char* tips[COLS] = {
            nullptr,
            "Samples: frames (times members), or places along the bins of distributions",
            "Pearson's correlation. Of periodic values (angles), the circular one: Jammalamadaka-SenGupta's of two (-1..1), Mardia's of one against a linear value (0..1, no sign; hover for the signed parts)",
            "The 95 % confidence interval of Pearson's r, by Fisher's z over the effective number of samples: successive frames are alike, and members move together, so they are worth fewer independent ones (Bretherton et al. 1999). Dimmed when it includes 0",
            "Spearman's rank correlation: how monotonic, whatever the shape",
            "Of the least squares line, Y on X",
        };
        for (int c = 0; c < COLS; ++c) {
            ImGui::TableSetColumnIndex(c);
            ImGui::TableHeader(ImGui::TableGetColumnName(c));
            if (tips[c]) ImGui::SetItemTooltip("%s", tips[c]);
        }
        auto row = [&](const PairStats& s, const char* name, ImVec4 c, bool indent) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            if (indent) ImGui::Indent();
            ImGui::ColorButton("##c", c, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder, ImVec2(ImGui::GetTextLineHeight() * 0.7f, ImGui::GetTextLineHeight() * 0.7f));
            ImGui::SameLine();
            ImGui::TextUnformatted(name);
            if (indent) ImGui::Unindent();
            if (s.members > 1) {
                ImGui::SameLine();
                ImGui::TextDisabled("(%d)", s.members);
            }
            ImGui::TableNextColumn();
            ImGui::Text("%llu", (unsigned long long)s.n);
            if (s.n_eff > 0.0) {
                ImGui::SetItemTooltip("Worth about %.0f independent samples\nlag-1 autocorrelation of X %.2f, of Y %.2f", s.n_eff, s.ax, s.ay);
            }
            ImGui::TableNextColumn();
            if (!s.has_corr()) ImGui::TextDisabled("-");
            else if (s.significant()) ImGui::Text("%.3f", s.corr());
            else ImGui::TextDisabled("%.3f", s.corr());
            if (s.circular == 2) ImGui::SetItemTooltip("Circular (both periodic)");
            else if (s.circular == 1 && s.parts_valid) {
                ImGui::SetItemTooltip("Circular-linear: how much the one predicts the other, no sign\nAgainst the cosine of the angle %+.3f, the sine %+.3f", s.r_cos, s.r_sin);
            } else if (s.circular == 1) ImGui::SetItemTooltip("Circular-linear: how much the one predicts the other, no sign");
            ImGui::TableNextColumn();
            if (s.ci_valid) {
                if (s.significant()) ImGui::Text("%.2f..%.2f", s.ci_lo, s.ci_hi);
                else ImGui::TextDisabled("%.2f..%.2f", s.ci_lo, s.ci_hi);
            } else ImGui::TextDisabled("-");
            ImGui::TableNextColumn();
            if (s.rho_valid) ImGui::Text("%.3f", s.rho);
            else ImGui::TextDisabled("-");
            if (s.rho_valid && s.rho_n < s.n) ImGui::SetItemTooltip("Of %llu samples, evenly spaced", (unsigned long long)s.rho_n);
            ImGui::TableNextColumn();
            if (s.has_fit()) ImGui::Text("%.4g", s.slope());
            else ImGui::TextDisabled("-");
        };
        for (int p = 0; p < l.num_pairs; ++p) {
            char name[192];
            pair_label(name, sizeof(name), l, p);
            row(l.pairs[p], name, l.num_pairs > 1 ? pair_color(l, p, 0) : vec_cast(l.color), false);
        }
        // Of a population, each of its members under all of them
        if (layer_by_member(l)) {
            for (int g = 0; g < layer_num_groups(l); ++g) {
                char name[192];
                group_label(name, sizeof(name), l, g);
                ImGui::PushID(g);
                row(group_stats(l, g), name, group_color(l, g), true);
                ImGui::PopID();
            }
        }
        ImGui::EndTable();
    }

    void draw_layer_tooltip(const Layer& l) const {
        char label[128];
        layer_label(label, sizeof(label), l);
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(label);
        if (!l.data_valid) {
            ImGui::TextColored(error_color, "%s", l.data_error[0] ? l.data_error : "no samples");
            ImGui::EndTooltip();
            return;
        }
        draw_pair_table(l);
        if (layer_draws_points(l) && l.point_stride > 1) {
            ImGui::TextDisabled("Every %u samples are shown", l.point_stride);
        }
        ImGui::EndTooltip();
    }

    // What the rows asked for, done once the list is drawn
    struct SidebarActions {
        uint32_t remove = 0;
        uint32_t duplicate = 0;
        bool     move = false;
        LayerNodeRef move_src, move_target;
        RepresentationDrop move_drop = RepresentationDrop_Before;
    };

    // Selects l. Of other axes than the plot's: the plot turns to them, l shown and the layers of the axes it had hidden,
    // so that what is shown is what is on the plot (more of these axes are then shown with their eyes).
    void select_layer(Layer& l) {
        selected_layer = l.id;
        if (l.axes_sig == 0 || l.axes_sig == active_sig) return;
        for (int j = 0; j < num_layers; ++j) {
            if (layers[j].axes_sig == active_sig) layers[j].visible = false;
        }
        l.visible = true;
        active_sig = l.axes_sig;
    }

    // Over the layers whose axes are of the same: what they are of, and the ones the plot shows marked. Clicked, those
    // shown, by selecting the first of them.
    void draw_axes_header(uint64_t sig) {
        Layer* first = nullptr;
        int count = 0;
        for (int i = 0; i < num_layers; ++i) {
            if (layers[i].axes_sig != sig) continue;
            if (!first) first = &layers[i];
            count += 1;
        }
        if (!first) return;
        const bool active = sig != 0 && sig == active_sig;
        char text[128];
        char tip[256];
        if (sig == 0) {
            snprintf(text, sizeof(text), "Incomplete");
            snprintf(tip, sizeof(tip), "Layers without a series on each axis (available): not shown");
        } else {
            char x[48], y[48];
            slot_axis_label(x, sizeof(x), first->x_info[0]);
            slot_axis_label(y, sizeof(y), first->y_info[0]);
            snprintf(text, sizeof(text), "%s  \xc3\x97  %s", x, y);
            snprintf(tip, sizeof(tip), "X in %s%s, Y in %s%s: the layers below share the axes.%s",
                     first->x_info[0].unit[0] ? first->x_info[0].unit : "no unit", slot_periodic(first->x_info[0]) ? " (periodic)" : "",
                     first->y_info[0].unit[0] ? first->y_info[0].unit : "no unit", slot_periodic(first->y_info[0]) ? " (periodic)" : "",
                     active ? " Shown." : " Click (or one of them, or its eye) to turn the plot to them: those shown now are hidden.");
        }

        const float h = ImGui::GetTextLineHeight() + 4.0f;
        const ImVec2 p0 = ImGui::GetCursorScreenPos();
        const float w = ImGui::GetContentRegionAvail().x;
        ImGui::Dummy(ImVec2(0, 2));
        ImGui::SetCursorScreenPos(ImVec2(p0.x, p0.y + 2));
        if (ImGui::InvisibleButton("##axes", ImVec2(w, h)) && sig != 0 && !active) select_layer(*first);
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) ImGui::SetTooltip("%s", tip);
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImVec2 r0 = ImGui::GetItemRectMin(), r1 = ImGui::GetItemRectMax();
        // Shown: marked along its edge, as the plot's
        if (active) dl->AddRectFilled(r0, ImVec2(r0.x + 3, r1.y), ImGui::GetColorU32(ImGuiCol_CheckMark), 1.0f);
        const ImU32 col = ImGui::GetColorU32(active ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        const float fs = ImGui::GetFontSize() * 0.85f;
        const float ty = r0.y + (h - fs) * 0.5f;
        dl->AddText(ImGui::GetFont(), fs, ImVec2(r0.x + 8, ty), col, text);
        char n[16];
        snprintf(n, sizeof(n), "%d", count);
        const ImVec2 ns = ImGui::GetFont()->CalcTextSizeA(fs, FLT_MAX, 0.0f, n);
        dl->AddText(ImGui::GetFont(), fs, ImVec2(r1.x - ns.x - 4, ty), ImGui::GetColorU32(ImGuiCol_TextDisabled), n);
        dl->AddLine(ImVec2(r0.x, r1.y), ImVec2(r1.x, r1.y), ImGui::GetColorU32(ImGuiCol_Separator, 0.6f));
        ImGui::SetCursorScreenPos(ImVec2(p0.x, r1.y + 2));
    }

    void draw_layer_row(Layer& l, RepresentationNode* const nodes[], size_t num_nodes, SidebarActions* act) {
        ImGui::PushID("layer");
        ImGui::PushID((int)l.id);

        const float line_h = ImGui::GetTextLineHeight();
        const float row_h = line_h * 2.0f + 8.0f;
        const ImVec2 row_min = ImGui::GetCursorScreenPos();
        const float row_w = ImGui::GetContentRegionAvail().x;
        const bool selected = selected_layer == l.id;
        const ImRect row_rect(row_min, ImVec2(row_min.x + row_w, row_min.y + row_h));

        char label[128];
        layer_label(label, sizeof(label), l);

        // The row, a selectable under the eye, the swatch and the text
        ImGui::SetNextItemAllowOverlap();
        if (ImGui::Selectable("##row", selected, ImGuiSelectableFlags_AllowOverlap, ImVec2(row_w, row_h))) {
            select_layer(l);
        }
        const bool row_hovered = ImGui::IsItemHovered();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNormal | ImGuiHoveredFlags_NoSharedDelay) && !ImGui::IsMouseDragging(0)) {
            draw_layer_tooltip(l);
        }
        layer_tree_drag_source({LIST_ID, l.id, false}, label);
        if (ImGui::BeginPopupContextItem("##menu")) {
            if (ImGui::IsWindowAppearing()) select_layer(l);
            for (int a = 0; a < 2; ++a) {
                if (ImGui::BeginMenu(a == 0 ? ICON_FA_PLUS " Series on X" : ICON_FA_PLUS " Series on Y")) {
                    draw_series_picker(l, a);
                    ImGui::EndMenu();
                }
            }
            ImGui::Separator();
            if (ImGui::MenuItem(ICON_FA_COPY " Duplicate", nullptr, false, num_layers < MAX_LAYERS)) act->duplicate = l.id;
            if (ImGui::MenuItem(ICON_FA_RIGHT_LEFT " Swap X and Y", nullptr, false, l.num_x + l.num_y > 0)) swap_axes(l);
            ImGui::Separator();
            if (ImGui::MenuItem(ICON_FA_XMARK " Remove")) act->remove = l.id;
            ImGui::EndPopup();
        }
        // A period from the timeline: this layer over its frames; a series: onto the axis it fits
        if (ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                l.frames = LayerFrames_Period;
                l.period = *(const uint32_t*)payload->Data;
                selected_layer = l.id;
            }
            ImGui::EndDragDropTarget();
        }
        // A layer of the list: moved here
        LayerNodeRef src;
        RepresentationDrop drop;
        if (layer_tree_drop_target(row_rect, row_min.x, {LIST_ID, l.id, false}, false, nodes, num_nodes, &src, &drop)) {
            act->move = true;
            act->move_src = src;
            act->move_target = {LIST_ID, l.id, false};
            act->move_drop = drop;
        }
        if (row_hovered && l.frames == LayerFrames_Period) {
            timeline_periods_highlight(l.period);
        }

        ImDrawList* dl = ImGui::GetWindowDrawList();
        const float pad = 4.0f;
        const ImVec2 text_col_min = ImVec2(row_min.x + pad + line_h + 6 + line_h + 6, row_min.y + 4);

        // The eye
        ImGui::SetCursorScreenPos(ImVec2(row_min.x + pad, row_min.y + (row_h - ImGui::GetFrameHeight()) * 0.5f));
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, ImGui::GetStyle().FramePadding.y));
        // Of other axes than the plot's: dim, whether shown or not
        const bool on_plot = l.axes_sig != 0 && l.axes_sig == active_sig;
        if (!on_plot) ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
        // Of other axes: shown, and selected, so that the plot turns to its axes (as a click on its card does)
        const bool switches = !on_plot && l.axes_sig != 0;
        if (ImGui::Button(l.visible ? ICON_FA_EYE "##vis" : ICON_FA_EYE_SLASH "##vis")) {
            if (switches) {
                select_layer(l);
            } else {
                l.visible = !l.visible;
            }
        }
        if (!on_plot) ImGui::PopStyleColor();
        ImGui::PopStyleVar();
        ImGui::PopStyleColor();
        if (switches) ImGui::SetItemTooltip("Show, on its axes: the plot turns to them, and the layers of the axes shown are hidden");
        else ImGui::SetItemTooltip("%s", l.visible ? "Hide" : "Show");

        // The swatch: what it looks like
        const ImVec2 sw_min = ImVec2(row_min.x + pad + line_h + 8, row_min.y + (row_h - line_h) * 0.5f);
        const ImVec2 sw_max = ImVec2(sw_min.x + line_h, sw_min.y + line_h);
        const float vis_alpha = l.visible && on_plot ? 1.0f : 0.35f;
        const DensityMapSwatch kind = l.style == LayerStyle_Points ? DensityMapSwatch_Points : l.style == LayerStyle_Contours ? DensityMapSwatch_Contours : DensityMapSwatch_Heatmap;
        density_map_swatch(dl, sw_min, sw_max, kind, vec_cast(l.color), l.fill, l.colormap, vis_alpha);

        // The name, and below it what it holds
        char stats[160];
        layer_stats(stats, sizeof(stats), l);
        const ImU32 text_col = ImGui::GetColorU32(l.visible && on_plot ? ImGuiCol_Text : ImGuiCol_TextDisabled);
        // On the selection's colour the disabled one is too faint to read
        const ImU32 dim_col = selected ? ImGui::GetColorU32(ImGuiCol_Text, 0.7f) : ImGui::GetColorU32(ImGuiCol_TextDisabled);
        const ImVec4 clip = ImVec4(text_col_min.x, row_min.y, row_min.x + row_w - pad, row_min.y + row_h);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(), text_col_min, text_col, label, nullptr, 0.0f, &clip);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 0.85f, ImVec2(text_col_min.x, text_col_min.y + line_h + 1),
                    !l.data_valid ? ImGui::ColorConvertFloat4ToU32(error_color) : dim_col, stats, nullptr, 0.0f, &clip);

        // Working on its density
        if (layer_busy(l)) density_map_spinner(dl, ImVec2(row_min.x + row_w - pad - 6, row_min.y + 9), dim_col);

        ImGui::SetCursorScreenPos(ImVec2(row_min.x, row_min.y + row_h + 2));
        ImGui::Dummy(ImVec2(0, 0));
        ImGui::PopID();
        ImGui::PopID();
    }

    // A new layer, first in the list, of the selected layer's series
    Layer* add_layer_like_selected(LayerFrames frames, LayerStyle style) {
        const Layer* src = find_layer(selected_layer);
        const uint32_t src_id = src ? src->id : 0;
        Layer* l = add_layer(frames, style, 0);
        if (!l) return nullptr;
        if (const Layer* s = find_layer(src_id)) copy_series(*l, *s);
        selected_layer = l->id;
        return l;
    }

    void draw_add_menu() {
        const bool room = num_layers < MAX_LAYERS;
        if (ImGui::MenuItem("Current frame", nullptr, false, room)) {
            if (Layer* l = add_layer_like_selected(LayerFrames_Current, LayerStyle_Points)) l->point_size = 8.0f;
        }
        if (ImGui::MenuItem("All frames", nullptr, false, room)) {
            add_layer_like_selected(LayerFrames_All, LayerStyle_Points);
        }
        const TimelinePeriods& ps = app->timeline.periods;
        if (ImGui::BeginMenu("Period", ps.count > 0 && room)) {
            for (int i = 0; i < ps.count; ++i) {
                char label[TIMELINE_PERIOD_NAME_CAP + 16];
                snprintf(label, sizeof(label), "%s###p%u", ps.items[i].name, ps.items[i].id);
                if (ImGui::MenuItem(label)) add_period_layer(ps.items[i].id);
                if (ImGui::IsItemHovered()) timeline_periods_highlight(ps.items[i].id);
            }
            ImGui::EndMenu();
        }
        if (ps.count == 0) ImGui::SetItemTooltip("Mark periods in the timeline (shift-drag) to correlate the series over them");
        ImGui::TextDisabled("New layers start with the selected layer's series");
    }

    void draw_levels(Layer& l) {
        prop::row("levels", "The contours enclose these fractions of the samples (highest density regions)");
        density_map_levels_editor(l.levels, &l.num_levels);
    }

    // ### Picking series

    // What there is to put on an axis, as the export window lists it: the script's properties over time, its
    // distributions, and the series loaded along the run (energies, ...), each under a heading
    struct PickSection {
        char label[96];
        size_t beg, end;            // Of the keys
    };
    std::vector<SeriesKey>   pick_keys;
    std::vector<PickSection> pick_sections;
    ImGuiTextFilter pick_filter;
    char pick_refusal[192] = "";    // Why the last one picked was not put on the axis

    void gather_pick_candidates() {
        pick_keys.clear();
        pick_sections.clear();
        auto section = [&](const char* label) {
            if (!pick_sections.empty()) pick_sections.back().end = pick_keys.size();
            PickSection ps = {};
            snprintf(ps.label, sizeof(ps.label), "%s", label);
            ps.beg = ps.end = pick_keys.size();
            pick_sections.push_back(ps);
        };
        section("Script, over time");
        series_for_each_script_property(app, SeriesSource_Script, MD_SCRIPT_PROPERTY_FLAG_TEMPORAL, [&](const SeriesKey& k) { pick_keys.push_back(k); });
        section("Script, distributions");
        series_for_each_script_property(app, SeriesSource_Script, MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION, [&](const SeriesKey& k) {
            if (is_script_distribution(app, k)) pick_keys.push_back(k);
        });
        str_t groups[64];
        const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), app), ARRAY_SIZE(groups));
        std::vector<str_t> paths;
        for (size_t g = 0; g < num_groups; ++g) {
            char label[96];
            system_series_group_label(label, sizeof(label), app, groups[g]);
            section(label);
            paths.resize(system_series_members(nullptr, 0, app, groups[g]) + 1);
            const size_t num = MIN(system_series_members(paths.data(), paths.size(), app, groups[g]), paths.size());
            for (size_t i = 0; i < num; ++i) pick_keys.push_back(series_key(SeriesSource_System, paths[i]));
        }
        if (!pick_sections.empty()) pick_sections.back().end = pick_keys.size();
    }

    // The candidates as items of the popup or menu open: picked, a series goes on the axis (the popup staying open, for
    // more). Those on it already are checked; one unlike them is refused, saying why.
    void draw_series_picker(Layer& l, int axis) {
        if (ImGui::IsWindowAppearing()) {
            gather_pick_candidates();
            pick_filter.Clear();
            pick_refusal[0] = '\0';
            ImGui::SetKeyboardFocusHere();
        }
        pick_filter.Draw("##filter", ImGui::GetFontSize() * 16);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Filter by name");
        const SeriesKey* on_axis = axis == 0 ? l.xs : l.ys;
        const int n_on = axis == 0 ? l.num_x : l.num_y;
        const bool full = n_on >= MAX_AXIS_SERIES;
        bool any = false;
        ImGui::BeginChild("##candidates", ImVec2(ImGui::GetFontSize() * 22, ImGui::GetTextLineHeightWithSpacing() * 14));
        for (const PickSection& ps : pick_sections) {
            bool header = false;
            for (size_t i = ps.beg; i < ps.end; ++i) {
                const SeriesKey& key = pick_keys[i];
                char label[128];
                series_label(label, sizeof(label), app, key);
                if (!pick_filter.PassFilter(label)) continue;
                if (!header) {
                    ImGui::SeparatorText(ps.label);
                    header = true;
                }
                any = true;
                bool present = false;
                for (int j = 0; j < n_on; ++j) present |= series_key_equal(on_axis[j], plotted_key(key));
                ImGui::PushID((int)i);
                ImGui::BeginDisabled(present || full);
                if (ImGui::Selectable(label, present, ImGuiSelectableFlags_DontClosePopups)) {
                    const char* reason = nullptr;
                    if (!series_droppable(key, &reason) || (reason = axis_refusal(l, axis, key)) != nullptr) {
                        snprintf(pick_refusal, sizeof(pick_refusal), "%s: %s", label, reason ? reason : "cannot be plotted");
                    } else {
                        pick_refusal[0] = '\0';
                        drop_series_on(l, axis, key);
                    }
                }
                ImGui::EndDisabled();
                ImGui::PopID();
            }
        }
        if (!any) ImGui::TextDisabled(pick_keys.empty() ? "No properties: write some in the script" : "None of that name");
        ImGui::EndChild();
        if (full) ImGui::TextDisabled("As many series on %s as there may be", axis == 0 ? "X" : "Y");
        if (pick_refusal[0]) {
            ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + ImGui::GetFontSize() * 22);
            ImGui::TextColored(error_color, "%s", pick_refusal);
            ImGui::PopTextWrapPos();
        }
    }

    // The series of an axis of a layer: listed, removed, picked (+) and dropped on from the timeline
    void draw_series_list(Layer& l, int axis) {
        ImGui::PushID(axis);
        const int n = axis == 0 ? l.num_x : l.num_y;
        const SlotInfo* infos = axis == 0 ? l.x_info : l.y_info;
        const float line_h = ImGui::GetFrameHeight();
        const float w = prop::control_width();
        const float h = MAX(1, n) * (line_h + 2.0f) + 4.0f;
        const ImVec2 p0 = ImGui::GetCursorScreenPos();

        ImGui::BeginGroup();
        int remove = -1;
        if (n == 0) {
            ImGui::Dummy(ImVec2(w, 2));
            ImGui::TextDisabled("  drop series here, or +");
        }
        for (int i = 0; i < n; ++i) {
            ImGui::PushID(i);
            const SlotInfo& info = infos[i];
            char name[128];
            slot_label(name, sizeof(name), l, axis, i);
            // On the axis the pairs vary along, the colour they are drawn in
            const bool varies = l.num_pairs > 1 && i < l.num_pairs && (axis == 0 ? l.num_x > 1 && l.num_y == 1 : l.num_y > 1);
            ImGui::ColorButton("##c", varies ? pair_color(l, i, 0) : info.color, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder, ImVec2(line_h * 0.5f, line_h * 0.5f));
            ImGui::SameLine();
            const float bw = line_h;
            const ImVec2 tp = ImGui::GetCursorScreenPos();
            ImGui::PushClipRect(tp, ImVec2(p0.x + w - bw - 4, tp.y + line_h), true);
            if (!info.ok) ImGui::TextColored(error_color, "%s", name);
            else if (info.unit[0]) ImGui::Text("%s  \xc2\xb7  %s", name, info.unit);
            else ImGui::TextUnformatted(name);
            ImGui::PopClipRect();
            if (ImGui::IsItemHovered()) {
                ImGui::BeginTooltip();
                ImGui::TextUnformatted(name);
                if (!info.ok) ImGui::TextColored(error_color, "Not available");
                else if (info.dim > 1) ImGui::TextDisabled("%d members", info.dim);
                ImGui::EndTooltip();
            }
            ImGui::SameLine();
            ImGui::SetCursorScreenPos(ImVec2(p0.x + w - bw, ImGui::GetCursorScreenPos().y));
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
            if (ImGui::Button(ICON_FA_XMARK "##rm", ImVec2(bw, 0))) remove = i;
            ImGui::PopStyleColor();
            ImGui::SetItemTooltip("Remove");
            ImGui::PopID();
        }
        // Picked from all there is
        {
            const float bw = line_h;
            const ImVec2 cur = ImGui::GetCursorScreenPos();
            ImGui::SetCursorScreenPos(ImVec2(p0.x + w - bw, n == 0 ? p0.y + 2.0f : cur.y));
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
            if (ImGui::Button(ICON_FA_PLUS "##pick", ImVec2(bw, 0))) ImGui::OpenPopup("##picker");
            ImGui::PopStyleColor();
            ImGui::SetItemTooltip("Add a series on %s", axis == 0 ? "X" : "Y");
            if (ImGui::BeginPopup("##picker")) {
                draw_series_picker(l, axis);
                ImGui::EndPopup();
            }
            if (n == 0) {
                ImGui::SetCursorScreenPos(cur);
                ImGui::Dummy(ImVec2(0, 0));
            }
        }
        ImGui::EndGroup();

        const ImRect r(p0, ImVec2(p0.x + w, MAX(p0.y + h, ImGui::GetItemRectMax().y)));
        ImGui::GetWindowDrawList()->AddRect(r.Min, r.Max, ImGui::GetColorU32(ImGuiCol_Border), 3.0f);
        if (ImGui::BeginDragDropTargetCustom(r, ImGui::GetID("##drop"))) {
            const ImGuiPayload* p = ImGui::AcceptDragDropPayload(TIMELINE_SERIES_DND, ImGuiDragDropFlags_AcceptBeforeDelivery);
            if (!p) p = ImGui::AcceptDragDropPayload(DISTRIBUTION_SERIES_DND, ImGuiDragDropFlags_AcceptBeforeDelivery);
            if (const SeriesDragPayload* sp = series_payload(p)) {
                const SeriesKey key = sp->key;
                const char* reason = nullptr;
                bool ok = series_droppable(key, &reason);
                if (ok && (reason = axis_refusal(l, axis, key)) != nullptr) ok = false;
                if (!ok && reason) ImGui::SetTooltip("%s", reason);
                if (p->IsDelivery() && ok) drop_series_on(l, axis, key);
            }
            ImGui::EndDragDropTarget();
        }
        if (remove >= 0) remove_series(l, axis, remove);
        ImGui::PopID();
    }

    void draw_inspector(Layer& l) {
        prop::begin("correlation");
        ImGui::PushID((int)l.id);

        char auto_name[128];
        {
            char saved[sizeof(l.name)];
            MEMCPY(saved, l.name, sizeof(saved));
            l.name[0] = '\0';
            layer_label(auto_name, sizeof(auto_name), l);
            MEMCPY(l.name, saved, sizeof(saved));
        }
        prop::row("name");
        ImGui::InputTextWithHint("##name", auto_name, l.name, sizeof(l.name));

        prop::row("X", "Series along X: one, or as many as Y to pair them in order");
        draw_series_list(l, 0);
        prop::row("Y", "Series along Y: one, or as many as X to pair them in order");
        draw_series_list(l, 1);
        prop::row("");
        if (ImGui::Button(ICON_FA_RIGHT_LEFT " Swap X and Y")) swap_axes(l);
        if (layer_shape(l) == CorrShape_Mismatch) {
            prop::hint("One series on an axis against several on the other, or as many on both (paired in order)");
        }
        if (l.members_truncated) {
            prop::hint("Populations of different sizes: paired up to the smaller");
        }
        if (l.binned) {
            prop::hint("Distributions: correlated bin by bin over their common range, of the frames below");
        } else {
            bool periodic = false;
            for (int i = 0; i < l.num_x; ++i) periodic |= series_periodic_valid(l.x_info[i].periodic);
            for (int i = 0; i < l.num_y; ++i) periodic |= series_periodic_valid(l.y_info[i].periodic);
            if (periodic) prop::hint("Periodic values: their correlation is the circular one. An axis of one period throughout pans around the circle (right-click it to cut it elsewhere); otherwise they are cut at their seam");
        }

        prop::row("frames");
        {
            char preview[TIMELINE_PERIOD_NAME_CAP + 16];
            const TimelinePeriods& ps = app->timeline.periods;
            const TimelinePeriod* cur = timeline_periods_find(&ps, l.period);
            if (l.frames == LayerFrames_Current)  snprintf(preview, sizeof(preview), "Current frame");
            else if (l.frames == LayerFrames_All) snprintf(preview, sizeof(preview), "All frames");
            else snprintf(preview, sizeof(preview), "%s", cur ? cur->name : "Period (missing)");
            if (ImGui::BeginCombo("##frames", preview)) {
                if (ImGui::Selectable("Current frame", l.frames == LayerFrames_Current)) l.frames = LayerFrames_Current;
                if (ImGui::Selectable("All frames", l.frames == LayerFrames_All)) l.frames = LayerFrames_All;
                if (ps.count > 0) ImGui::SeparatorText("Periods");
                for (int i = 0; i < ps.count; ++i) {
                    char label[TIMELINE_PERIOD_NAME_CAP + 16];
                    snprintf(label, sizeof(label), "%s###p%u", ps.items[i].name, ps.items[i].id);
                    if (ImGui::Selectable(label, l.frames == LayerFrames_Period && l.period == ps.items[i].id)) {
                        l.frames = LayerFrames_Period;
                        l.period = ps.items[i].id;
                    }
                    if (ImGui::IsItemHovered()) timeline_periods_highlight(ps.items[i].id);
                }
                if (ps.count == 0) ImGui::TextDisabled("Mark periods in the timeline to use them here");
                ImGui::EndCombo();
            }
            if (l.frames != LayerFrames_Current && run_num_frames(app) < 2 && l.style != LayerStyle_Points) {
                prop::hint("No trajectory to find a density over: drawn as points");
            }
            if (l.frames == LayerFrames_Current && l.style != LayerStyle_Points) {
                prop::hint("One frame: drawn as points");
            }
        }

        prop::row("style");
        {
            int s = (int)l.style;
            static const int values[] = { LayerStyle_Points, LayerStyle_Contours, LayerStyle_Heatmap };
            density_map_segmented("##style", style_names, values, LayerStyle_Count, &s, nullptr, true);
            l.style = (LayerStyle)s;
        }

        // A heat map is coloured by its colour map, points by time too
        if (l.style != LayerStyle_Heatmap && !(layer_draws_points(l) && l.point_coloring == PointColoring_Time)) {
            prop::row("color");
            prop::color("##color", l.color.elem, true);
        }
        prop::row("opacity");
        ImGui::SliderFloat("##opacity", &l.opacity, 0.0f, 1.0f, "%.2f");

        if (l.style == LayerStyle_Contours) {
            prop::row("draw");
            ImGui::Checkbox("fill", &l.fill);
            ImGui::SameLine();
            ImGui::Checkbox("lines", &l.lines);
            if (l.num_pairs > 1) {
                prop::row("split", "One density per pair, each of its own samples, in the colour of its series: rather than one of them all, in the layer's");
                ImGui::Checkbox("by pair", &l.split_pairs);
            } else if (layer_by_member(l)) {
                prop::row("split", "One density per member of the population, each of its own samples, in its colour: rather than one of them all, in the layer's");
                ImGui::Checkbox("by member", &l.split_pairs);
            } else if (l.num_pairs == 1 && l.pairs[0].members > MAX_GROUPS) {
                prop::row("split");
                ImGui::TextDisabled("%d members: too many to split (%d at most)", l.pairs[0].members, MAX_GROUPS);
            }
            if (l.lines) {
                prop::row("line width");
                ImGui::SliderFloat("##line_width", &l.line_width, 0.5f, 6.0f, "%.1f px");
            }
            draw_levels(l);
        }
        if (l.style == LayerStyle_Heatmap) {
            prop::row("colormap");
            ImPlot::ColormapSelection("##cmap", &l.colormap);
            prop::row("log scale", "Four decades below the maximum");
            ImGui::Checkbox("##log", &l.log_scale);
        }
        if (l.style != LayerStyle_Points) {
            prop::row("bandwidth", "The standard deviation of the gaussian the density is smoothed with, in % of the extent of the samples");
            ImGui::SliderFloat("##bandwidth", &l.bandwidth, 0.5f, 20.0f, "%.1f %%", ImGuiSliderFlags_AlwaysClamp);
        }
        if (layer_draws_points(l)) {
            prop::row("point size");
            ImGui::SliderFloat("##point_size", &l.point_size, 1.0f, 12.0f, "%.0f px");
            prop::row("color by", "Layer: its colour. Pair: the series that varies between pairs, or the members of a population. Time: by the colour map");
            int c = (int)l.point_coloring;
            if (ImGui::Combo("##coloring", &c, point_coloring_names, PointColoring_Count)) l.point_coloring = (PointColoring)c;
            const bool members = l.num_pairs == 1 && l.pairs[0].members > 1;
            if (l.point_coloring == PointColoring_Time || (l.point_coloring == PointColoring_Pair && members)) {
                prop::row("colormap");
                ImPlot::ColormapSelection("##cmap_points", &l.colormap);
            }
            if (l.point_stride > 1) {
                prop::hint("Every %u samples, to keep to %u points", l.point_stride, (uint32_t)MAX_POINTS);
            }
        }
        prop::row("fit", "The least squares line of each pair");
        ImGui::Checkbox("##fit", &l.show_fit);

        ImGui::PopID();
        prop::end();

        if (l.data_valid && l.num_pairs > 0) {
            ImGui::Spacing();
            draw_pair_table(l);
        } else if (l.data_error[0]) {
            ImGui::Spacing();
            ImGui::TextColored(error_color, "%s", l.data_error);
        }
    }

    void draw_sidebar() {
        // The header: the title and adding
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("LAYERS");
        ImGui::SameLine();
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(0.0f, ImGui::GetContentRegionAvail().x - ImGui::GetFrameHeight()));
        if (ImGui::Button(ICON_FA_PLUS "##add", ImVec2(ImGui::GetFrameHeight(), 0))) ImGui::OpenPopup("##add_menu");
        ImGui::SetItemTooltip("Add a layer");
        if (ImGui::BeginPopup("##add_menu")) {
            draw_add_menu();
            ImGui::EndPopup();
        }

        SidebarActions act;
        {
            RepresentationNode* nodes[MAX_LAYERS];
            const size_t n = tree_nodes(nodes);
            // Under what their axes are of, in the order of the first layer of each; those still without, last
            uint64_t sigs[MAX_LAYERS + 1];
            int num_sigs = 0;
            bool incomplete = false;
            for (int i = 0; i < num_layers; ++i) {
                const uint64_t sig = layers[i].axes_sig;
                if (sig == 0) {
                    incomplete = true;
                    continue;
                }
                bool seen = false;
                for (int k = 0; k < num_sigs && !seen; ++k) seen = sigs[k] == sig;
                if (!seen) sigs[num_sigs++] = sig;
            }
            if (incomplete) sigs[num_sigs++] = 0;
            for (int k = 0; k < num_sigs; ++k) {
                ImGui::PushID(k);
                draw_axes_header(sigs[k]);
                for (int i = 0; i < num_layers; ++i) {
                    if (layers[i].axes_sig == sigs[k]) draw_layer_row(layers[i], nodes, n, &act);
                }
                ImGui::PopID();
            }
        }
        if (num_layers == 0) {
            ImGui::TextDisabled("No layers: add one with +");
        }
        // Below the list: the end of it
        {
            RepresentationNode* nodes[MAX_LAYERS];
            const size_t n = tree_nodes(nodes);
            const ImVec2 p = ImGui::GetCursorScreenPos();
            const ImRect r(p, ImVec2(p.x + ImGui::GetContentRegionAvail().x, p.y + ImGui::GetFrameHeight()));
            ImGui::Dummy(ImVec2(0, ImGui::GetFrameHeight() * 0.25f));
            LayerNodeRef src;
            if (layer_tree_drop_target_end(r, r.Min.x, LIST_ID, 0, nodes, n, &src)) {
                act.move = true;
                act.move_src = src;
                act.move_target = {};
                act.move_drop = RepresentationDrop_IntoLast;
            }
        }

        if (act.move) move_node(act.move_src, act.move_target, act.move_drop);
        if (act.duplicate) duplicate_layer(act.duplicate);
        if (act.remove) remove_layer(act.remove);

        if (Layer* l = find_layer(selected_layer)) {
            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();
            draw_inspector(*l);
        }

        // Delete takes the selected layer away
        if (ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows) && !ImGui::GetIO().WantTextInput && ImGui::IsKeyPressed(ImGuiKey_Delete)) {
            if (selected_layer) remove_layer(selected_layer);
        }
        if (tree_dirty) apply_tree();
    }
