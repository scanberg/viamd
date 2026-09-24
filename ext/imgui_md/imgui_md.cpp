/*
 * imgui_md: Markdown for Dear ImGui using MD4C
 * (https://github.com/mekhontsev/imgui_md)
 *
 * Copyright (c) 2021 Dmitry Mekhontsev
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 * Modified for viamd, see imgui_md.h for a summary of the changes.
 */

#include "imgui_md.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

void append_utf8(std::string& s, unsigned cp) {
	if (cp < 0x80) {
		s += (char)cp;
	} else if (cp < 0x800) {
		s += (char)(0xC0 | (cp >> 6));
		s += (char)(0x80 | (cp & 0x3F));
	} else if (cp < 0x10000) {
		s += (char)(0xE0 | (cp >> 12));
		s += (char)(0x80 | ((cp >> 6) & 0x3F));
		s += (char)(0x80 | (cp & 0x3F));
	} else {
		s += (char)(0xF0 | (cp >> 18));
		s += (char)(0x80 | ((cp >> 12) & 0x3F));
		s += (char)(0x80 | ((cp >> 6) & 0x3F));
		s += (char)(0x80 | (cp & 0x3F));
	}
}

// Decodes "&name;" / "&#123;" / "&#x1F;". Returns an empty string for entities we do not know.
std::string decode_entity(const char* str, const char* end) {
	std::string out;
	if (end - str < 3 || str[0] != '&' || end[-1] != ';') return out;
	const std::string name(str + 1, end - 1);
	if (name[0] == '#') {
		const bool hex = name.size() > 1 && (name[1] == 'x' || name[1] == 'X');
		const unsigned cp = (unsigned)strtoul(name.c_str() + (hex ? 2 : 1), nullptr, hex ? 16 : 10);
		if (cp) append_utf8(out, cp);
		return out;
	}
	static const struct { const char* name; unsigned cp; } table[] = {
		{"amp", '&'}, {"lt", '<'}, {"gt", '>'}, {"quot", '"'}, {"apos", '\''}, {"nbsp", ' '},
		{"copy", 0xA9}, {"deg", 0xB0}, {"plusmn", 0xB1}, {"times", 0xD7}, {"middot", 0xB7},
		{"ndash", 0x2013}, {"mdash", 0x2014}, {"hellip", 0x2026}, {"rarr", 0x2192}, {"larr", 0x2190},
		{"le", 0x2264}, {"ge", 0x2265}, {"pi", 0x3C0}, {"Aring", 0xC5},
	};
	for (const auto& e : table) {
		if (name == e.name) { append_utf8(out, e.cp); return out; }
	}
	return out;
}

ImVec4 mix(const ImVec4& a, const ImVec4& b, float t) {
	return ImVec4(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t, a.w + (b.w - a.w) * t);
}

// Pushes the text colour towards white on dark themes and towards black on light ones.
ImVec4 emphasize(const ImVec4& c, float t) {
	const float lum = (c.x + c.y + c.z) / 3.0f;
	return mix(c, lum > 0.5f ? ImVec4(1, 1, 1, c.w) : ImVec4(0, 0, 0, c.w), t);
}

////////////////////////////////////////////////////////////////////////////////
// Table pre-pass: estimate how wide every column wants to be.
////////////////////////////////////////////////////////////////////////////////

struct TablePrepass {
	std::vector<std::vector<float>>* out = nullptr;
	int   col     = 0;
	float len     = 0.0f;
	bool  in_cell = false;
};

constexpr float kMinColumnWeight = 6.0f;
constexpr float kMaxColumnWeight = 64.0f;

} // namespace

imgui_md::imgui_md()
{
	m_md.abi_version = 0;
	m_md.flags = MD_FLAG_TABLES | MD_FLAG_STRIKETHROUGH;

	m_md.enter_block = [](MD_BLOCKTYPE t, void* d, void* u) { return ((imgui_md*)u)->block(t, d, true); };
	m_md.leave_block = [](MD_BLOCKTYPE t, void* d, void* u) { return ((imgui_md*)u)->block(t, d, false); };
	m_md.enter_span  = [](MD_SPANTYPE t, void* d, void* u) { return ((imgui_md*)u)->span(t, d, true); };
	m_md.leave_span  = [](MD_SPANTYPE t, void* d, void* u) { return ((imgui_md*)u)->span(t, d, false); };
	m_md.text = [](MD_TEXTTYPE t, const MD_CHAR* text, MD_SIZE size, void* u) {
		return ((imgui_md*)u)->text(t, text, text + size);
	};
	m_md.debug_log = nullptr;
	m_md.syntax = nullptr;
}

void imgui_md::prepass_tables(const char* str, const char* str_end)
{
	m_table_weights.clear();

	TablePrepass state;
	state.out = &m_table_weights;

	MD_PARSER p = {};
	p.abi_version = 0;
	p.flags = m_md.flags;
	p.enter_block = [](MD_BLOCKTYPE t, void* d, void* u) {
		TablePrepass& s = *(TablePrepass*)u;
		if (t == MD_BLOCK_TABLE) {
			s.out->emplace_back((size_t)((MD_BLOCK_TABLE_DETAIL*)d)->col_count, kMinColumnWeight);
		} else if (t == MD_BLOCK_TR) {
			s.col = 0;
		} else if (t == MD_BLOCK_TH || t == MD_BLOCK_TD) {
			s.in_cell = true;
			s.len = 0.0f;
		}
		return 0;
	};
	p.leave_block = [](MD_BLOCKTYPE t, void*, void* u) {
		TablePrepass& s = *(TablePrepass*)u;
		if (t == MD_BLOCK_TH || t == MD_BLOCK_TD) {
			auto& w = s.out->back();
			if ((size_t)s.col < w.size()) {
				w[s.col] = std::max(w[s.col], std::min(s.len, kMaxColumnWeight));
			}
			++s.col;
			s.in_cell = false;
		}
		return 0;
	};
	p.enter_span = [](MD_SPANTYPE, void*, void*) { return 0; };
	p.leave_span = [](MD_SPANTYPE, void*, void*) { return 0; };
	p.text = [](MD_TEXTTYPE, const MD_CHAR* text, MD_SIZE size, void* u) {
		TablePrepass& s = *(TablePrepass*)u;
		if (s.in_cell) {
			for (MD_SIZE i = 0; i < size; ++i) {
				if (((unsigned char)text[i] & 0xC0) != 0x80) s.len += 1.0f; // count code points
			}
		}
		return 0;
	};

	md_parse(str, (MD_SIZE)(str_end - str), &p, &state);
}

int imgui_md::print(const char* str, const char* str_end)
{
	if (str >= str_end) return 0;

	m_href.clear();
	m_is_underline = m_is_strikethrough = m_is_em = m_is_strong = m_is_code = false;
	m_is_table_header = m_is_image = false;
	m_hlevel = 0;
	m_code_index = 0;
	m_inline = false;
	m_first_block = true;
	m_skip = false;
	m_in_code_block = false;
	m_list_stack.clear();
	m_quote_stack.clear();
	m_in_table = false;
	m_table_index = -1;

	prepass_tables(str, str_end);

	return md_parse(str, (MD_SIZE)(str_end - str), &m_md, this);
}

////////////////////////////////////////////////////////////////////////////////
// Style hooks

float imgui_md::get_text_scale() const
{
	switch (m_hlevel) {
	case 1:  return 1.55f;
	case 2:  return 1.30f;
	case 3:  return 1.15f;
	default: return 1.0f;
	}
}

ImVec4 imgui_md::get_text_color() const
{
	const ImGuiStyle& style = ImGui::GetStyle();
	const ImVec4 text = style.Colors[ImGuiCol_Text];

	if (!m_href.empty())  return style.Colors[ImGuiCol_TextLink];
	if (m_hlevel)         return emphasize(text, m_hlevel <= 2 ? 0.6f : 0.35f);
	if (m_is_code)        return mix(text, ImVec4(1.00f, 0.66f, 0.36f, text.w), 0.85f);
	if (m_is_strong || m_is_table_header) return emphasize(text, 0.85f);
	if (m_is_em)          return mix(text, ImVec4(0.55f, 0.78f, 1.00f, text.w), 0.6f);
	return text;
}

float imgui_md::get_block_gap() const
{
	return ImGui::GetFontSize() * 0.3f;
}

void imgui_md::open_url(const std::string&) {}

void imgui_md::on_heading(int, const std::string&, float) {}

////////////////////////////////////////////////////////////////////////////////
// Layout helpers

void imgui_md::block_gap()
{
	end_line();
	if (m_first_block) {
		m_first_block = false;
		return;
	}
	if (!m_list_stack.empty() || m_in_table || m_skip) return;
	ImGui::Dummy(ImVec2(0.0f, get_block_gap()));
}

void imgui_md::render_text(const char* str, const char* end)
{
	if (m_skip || m_is_image) return;

	// Spaces at the start of a line are not shown (they are what a wrap leaves behind).
	if (!m_inline) {
		while (str < end && *str == ' ') ++str;
	}
	if (str >= end) return;

	const float scale = get_text_scale();
	ImFont* font = get_font();
	const bool push_font = (font != nullptr) || (scale != 1.0f);
	if (push_font) ImGui::PushFont(font, ImGui::GetStyle().FontSizeBase * scale);
	ImGui::PushStyleColor(ImGuiCol_Text, get_text_color());

	const ImGuiStyle& style = ImGui::GetStyle();
	const float font_size = ImGui::GetFontSize();
	ImDrawList* dl = ImGui::GetWindowDrawList();

	while (str < end) {
		const ImVec2 cursor = ImGui::GetCursorScreenPos();
		const float line_start = cursor.x;
		const float right = cursor.x + ImGui::GetContentRegionAvail().x;
		const float line_x = m_inline ? ImGui::GetItemRectMax().x : line_start;

		const char* te = ImGui::GetFont()->CalcWordWrapPosition(font_size, str, end, right - line_x);
		if (te == str) {
			if (m_inline && line_x > line_start + 0.5f) {
				// Not even the first word fits on the rest of this line: start a new one.
				m_inline = false;
				while (str < end && *str == ' ') ++str;
				continue;
			}
			// A single glyph wider than the whole line. Emit it anyway to guarantee progress.
			te = str + 1;
			while (te < end && ((unsigned char)*te & 0xC0) == 0x80) ++te;
		}

		if (m_inline && te < end && *te != ' ' && !memchr(str, ' ', (size_t)(te - str)) && line_x > line_start + 0.5f) {
			// The line ran out in the middle of a word. Words are not split unless they are longer than a whole
			// line, so move the word to the next line instead.
			m_inline = false;
			continue;
		}

		if (m_inline) ImGui::SameLine(0.0f, 0.0f);

		if (m_is_code) {
			const ImVec2 p = ImGui::GetCursorScreenPos();
			const ImVec2 sz = ImGui::CalcTextSize(str, te);
			ImVec4 bg = style.Colors[ImGuiCol_FrameBg];
			bg.w *= 0.9f;
			dl->AddRectFilled(ImVec2(p.x - 1.0f, p.y), ImVec2(p.x + sz.x + 1.0f, p.y + sz.y), ImGui::GetColorU32(bg), 2.0f);
		}

		ImGui::TextUnformatted(str, te);

		const ImVec2 mi = ImGui::GetItemRectMin();
		const ImVec2 ma = ImGui::GetItemRectMax();

		if (!m_href.empty()) {
			const bool hovered = ImGui::IsItemHovered();
			ImVec4 c = style.Colors[ImGuiCol_TextLink];
			if (hovered) {
				c = emphasize(c, 0.35f);
				ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
				if (m_href[0] != '#') ImGui::SetTooltip("%s", m_href.c_str());
				if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) open_url(m_href);
			} else {
				c.w *= 0.55f;
			}
			dl->AddLine(ImVec2(mi.x, ma.y - 1.0f), ImVec2(ma.x, ma.y - 1.0f), ImGui::GetColorU32(c), 1.0f);
		}
		if (m_is_underline) {
			dl->AddLine(ImVec2(mi.x, ma.y - 1.0f), ImVec2(ma.x, ma.y - 1.0f), ImGui::GetColorU32(ImGuiCol_Text), 1.0f);
		}
		if (m_is_strikethrough) {
			const float y = (mi.y + ma.y) * 0.5f;
			dl->AddLine(ImVec2(mi.x, y), ImVec2(ma.x, y), ImGui::GetColorU32(ImGuiCol_Text), 1.0f);
		}

		str = te;
		if (str < end) {
			// Wrapped: the rest goes on a new line and leading spaces are dropped.
			m_inline = false;
			while (str < end && *str == ' ') ++str;
		} else {
			m_inline = true;
		}
	}

	ImGui::PopStyleColor();
	if (push_font) ImGui::PopFont();
}

bool imgui_md::render_entity(const char* str, const char* end)
{
	const std::string s = decode_entity(str, end);
	if (s.empty()) return false;
	if (m_hlevel) m_heading_text += s;
	render_text(s.data(), s.data() + s.size());
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Blocks

void imgui_md::BLOCK_DOC(bool) {}

void imgui_md::BLOCK_QUOTE(bool e)
{
	if (e) {
		block_gap();
		m_quote_stack.push_back(ImGui::GetCursorScreenPos().y);
		ImGui::Indent(ImGui::GetFontSize() * 0.9f);
	} else {
		end_line();
		ImGui::Unindent(ImGui::GetFontSize() * 0.9f);
		const float y0 = m_quote_stack.back();
		m_quote_stack.pop_back();
		const float y1 = ImGui::GetCursorScreenPos().y - ImGui::GetStyle().ItemSpacing.y;
		const float x  = ImGui::GetCursorScreenPos().x + 2.0f;
		ImVec4 c = ImGui::GetStyle().Colors[ImGuiCol_TextDisabled];
		c.w *= 0.8f;
		ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(x, y0), ImVec2(x + 3.0f, std::max(y0 + 1.0f, y1)), ImGui::GetColorU32(c));
	}
}

void imgui_md::BLOCK_UL(const MD_BLOCK_UL_DETAIL* d, bool e)
{
	if (e) {
		if (m_list_stack.empty()) block_gap();
		end_line();
		m_list_stack.push_back(list_info{0, d->mark, false});
	} else {
		m_list_stack.pop_back();
		end_line();
	}
}

void imgui_md::BLOCK_OL(const MD_BLOCK_OL_DETAIL* d, bool e)
{
	if (e) {
		if (m_list_stack.empty()) block_gap();
		end_line();
		m_list_stack.push_back(list_info{d->start, d->mark_delimiter, true});
	} else {
		m_list_stack.pop_back();
		end_line();
	}
}

void imgui_md::BLOCK_LI(const MD_BLOCK_LI_DETAIL*, bool e)
{
	const float indent = ImGui::GetFontSize() * 1.6f;
	if (e) {
		end_line();
		m_first_block = false;
		list_info& info = m_list_stack.back();
		const float fs = ImGui::GetFontSize();
		const ImVec2 p = ImGui::GetCursorScreenPos();
		ImDrawList* dl = ImGui::GetWindowDrawList();
		const ImU32 col = ImGui::GetColorU32(ImGuiCol_Text);
		if (info.is_ol) {
			char buf[24];
			snprintf(buf, sizeof(buf), "%u%c", info.cur_ol++, info.delim);
			const ImVec2 sz = ImGui::CalcTextSize(buf);
			dl->AddText(ImVec2(p.x + indent - fs * 0.35f - sz.x, p.y), col, buf);
		} else {
			const ImVec2 c(p.x + indent * 0.5f - fs * 0.15f, p.y + fs * 0.5f);
			const float r = fs * 0.14f;
			const size_t depth = m_list_stack.size();
			if (depth <= 1) {
				dl->AddCircleFilled(c, r, col);
			} else if (depth == 2) {
				dl->AddCircle(c, r, col, 0, 1.2f);
			} else {
				dl->AddRectFilled(ImVec2(c.x - r * 0.8f, c.y - r * 0.8f), ImVec2(c.x + r * 0.8f, c.y + r * 0.8f), col);
			}
		}
		ImGui::Indent(indent);
	} else {
		end_line();
		ImGui::Unindent(indent);
	}
}

void imgui_md::BLOCK_HR(bool e)
{
	if (e) {
		block_gap();
		ImGui::Separator();
	}
}

void imgui_md::BLOCK_H(const MD_BLOCK_H_DETAIL* d, bool e)
{
	if (e) {
		const bool first = m_first_block;
		block_gap();
		if (!first && !m_in_table && m_list_stack.empty()) {
			ImGui::Dummy(ImVec2(0.0f, get_block_gap() * (d->level <= 3 ? 1.5f : 0.5f)));
		}
		m_hlevel = d->level;
		m_heading_text.clear();
		m_heading_y = ImGui::GetCursorScreenPos().y;
	} else {
		on_heading((int)d->level, m_heading_text, m_heading_y);
		m_hlevel = 0;
		end_line();
		if (d->level <= 2) ImGui::Separator();
	}
}

void imgui_md::BLOCK_CODE(const MD_BLOCK_CODE_DETAIL* d, bool e)
{
	if (e) {
		m_in_code_block = true;
		m_code_buf.clear();
		m_code_lang.assign(d->lang.text ? d->lang.text : "", d->lang.text ? d->lang.size : 0);
	} else {
		m_in_code_block = false;
		while (!m_code_buf.empty() && (m_code_buf.back() == '\n' || m_code_buf.back() == '\r')) m_code_buf.pop_back();
		if (!m_skip) {
			block_gap();
			ImGui::PushID(m_code_index);
			render_code_block(m_code_lang, m_code_buf);
			ImGui::PopID();
		}
		++m_code_index;
	}
}

void imgui_md::render_code_block(const std::string&, const std::string& code)
{
	const ImGuiStyle& style = ImGui::GetStyle();
	const float pad = style.FramePadding.y + 2.0f;
	const ImVec2 sz = ImGui::CalcTextSize(code.c_str(), code.c_str() + code.size());
	const ImVec2 p = ImGui::GetCursorScreenPos();
	const float w = ImGui::GetContentRegionAvail().x;
	ImDrawList* dl = ImGui::GetWindowDrawList();
	dl->AddRectFilled(p, ImVec2(p.x + w, p.y + sz.y + pad * 2.0f), ImGui::GetColorU32(ImGuiCol_FrameBg), 3.0f);
	dl->AddText(ImVec2(p.x + pad, p.y + pad), ImGui::GetColorU32(ImGuiCol_Text), code.c_str(), code.c_str() + code.size());
	ImGui::Dummy(ImVec2(w, sz.y + pad * 2.0f));
}

void imgui_md::BLOCK_HTML(bool) {}

void imgui_md::BLOCK_P(bool e)
{
	if (e) {
		block_gap();
	} else {
		end_line();
	}
}

void imgui_md::BLOCK_TABLE(const MD_BLOCK_TABLE_DETAIL* d, bool e)
{
	if (e) {
		block_gap();
		++m_table_index;
		m_table_ncols = (int)d->col_count;
		char id[32];
		snprintf(id, sizeof(id), "##mdtable%d", m_table_counter++);
		const ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_NoSavedSettings;
		if (ImGui::BeginTable(id, m_table_ncols, flags)) {
			m_in_table = true;
			for (int c = 0; c < m_table_ncols; ++c) {
				float w = kMinColumnWeight;
				if ((size_t)m_table_index < m_table_weights.size() && (size_t)c < m_table_weights[m_table_index].size()) {
					w = m_table_weights[m_table_index][c];
				}
				ImGui::TableSetupColumn(nullptr, ImGuiTableColumnFlags_WidthStretch, w);
			}
		} else {
			m_skip = true;
		}
	} else {
		end_line();
		if (m_in_table) ImGui::EndTable();
		m_in_table = false;
		m_skip = false;
	}
}

void imgui_md::BLOCK_THEAD(bool e) { m_is_table_header = e; }
void imgui_md::BLOCK_TBODY(bool) {}

void imgui_md::BLOCK_TR(bool e)
{
	if (e && m_in_table) {
		ImGui::TableNextRow(m_is_table_header ? ImGuiTableRowFlags_Headers : 0);
	}
	end_line();
}

void imgui_md::BLOCK_TH(const MD_BLOCK_TD_DETAIL* d, bool e) { BLOCK_TD(d, e); }

void imgui_md::BLOCK_TD(const MD_BLOCK_TD_DETAIL*, bool e)
{
	end_line();
	if (e && m_in_table) {
		ImGui::TableNextColumn();
	}
}

////////////////////////////////////////////////////////////////////////////////
// Spans

void imgui_md::SPAN_EM(bool e)     { m_is_em = e; }
void imgui_md::SPAN_STRONG(bool e) { m_is_strong = e; }
void imgui_md::SPAN_U(bool e)      { m_is_underline = e; }
void imgui_md::SPAN_DEL(bool e)    { m_is_strikethrough = e; }
void imgui_md::SPAN_CODE(bool e)   { m_is_code = e; }

void imgui_md::SPAN_A(const MD_SPAN_A_DETAIL* d, bool e)
{
	if (e) {
		m_href.assign(d->href.text ? d->href.text : "", d->href.text ? d->href.size : 0);
	} else {
		m_href.clear();
	}
}

void imgui_md::SPAN_IMG(const MD_SPAN_IMG_DETAIL*, bool e)
{
	m_is_image = e; // alt text is suppressed, images are not drawn
}

////////////////////////////////////////////////////////////////////////////////
// md4c glue

int imgui_md::text(MD_TEXTTYPE type, const char* str, const char* str_end)
{
	if (m_in_code_block) {
		if (type != MD_TEXT_NULLCHAR) m_code_buf.append(str, str_end);
		return 0;
	}

	switch (type) {
	case MD_TEXT_NORMAL:
	case MD_TEXT_CODE:
		if (m_hlevel) m_heading_text.append(str, str_end);
		render_text(str, str_end);
		break;
	case MD_TEXT_ENTITY:
		if (!render_entity(str, str_end)) {
			if (m_hlevel) m_heading_text.append(str, str_end);
			render_text(str, str_end);
		}
		break;
	case MD_TEXT_SOFTBR:
		if (m_hlevel) m_heading_text += ' ';
		render_text(" ", " " + 1);
		break;
	case MD_TEXT_BR:
		end_line();
		break;
	case MD_TEXT_HTML:
		if (str_end - str >= 3 && strncmp(str, "<br", 3) == 0) end_line();
		break;
	default:
		break;
	}
	return 0;
}

int imgui_md::block(MD_BLOCKTYPE type, void* d, bool e)
{
	switch (type) {
	case MD_BLOCK_DOC:   BLOCK_DOC(e); break;
	case MD_BLOCK_QUOTE: BLOCK_QUOTE(e); break;
	case MD_BLOCK_UL:    BLOCK_UL((MD_BLOCK_UL_DETAIL*)d, e); break;
	case MD_BLOCK_OL:    BLOCK_OL((MD_BLOCK_OL_DETAIL*)d, e); break;
	case MD_BLOCK_LI:    BLOCK_LI((MD_BLOCK_LI_DETAIL*)d, e); break;
	case MD_BLOCK_HR:    BLOCK_HR(e); break;
	case MD_BLOCK_H:     BLOCK_H((MD_BLOCK_H_DETAIL*)d, e); break;
	case MD_BLOCK_CODE:  BLOCK_CODE((MD_BLOCK_CODE_DETAIL*)d, e); break;
	case MD_BLOCK_HTML:  BLOCK_HTML(e); break;
	case MD_BLOCK_P:     BLOCK_P(e); break;
	case MD_BLOCK_TABLE: BLOCK_TABLE((MD_BLOCK_TABLE_DETAIL*)d, e); break;
	case MD_BLOCK_THEAD: BLOCK_THEAD(e); break;
	case MD_BLOCK_TBODY: BLOCK_TBODY(e); break;
	case MD_BLOCK_TR:    BLOCK_TR(e); break;
	case MD_BLOCK_TH:    BLOCK_TH((MD_BLOCK_TD_DETAIL*)d, e); break;
	case MD_BLOCK_TD:    BLOCK_TD((MD_BLOCK_TD_DETAIL*)d, e); break;
	default: break;
	}
	return 0;
}

int imgui_md::span(MD_SPANTYPE type, void* d, bool e)
{
	switch (type) {
	case MD_SPAN_EM:     SPAN_EM(e); break;
	case MD_SPAN_STRONG: SPAN_STRONG(e); break;
	case MD_SPAN_A:      SPAN_A((MD_SPAN_A_DETAIL*)d, e); break;
	case MD_SPAN_IMG:    SPAN_IMG((MD_SPAN_IMG_DETAIL*)d, e); break;
	case MD_SPAN_CODE:   SPAN_CODE(e); break;
	case MD_SPAN_DEL:    SPAN_DEL(e); break;
	case MD_SPAN_U:      SPAN_U(e); break;
	default: break;
	}
	return 0;
}
