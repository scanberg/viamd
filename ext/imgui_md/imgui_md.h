/*
 * imgui_md: Markdown for Dear ImGui using MD4C
 * (https://github.com/mekhontsev/imgui_md)
 *
 * Copyright (c) 2021 Dmitry Mekhontsev (MIT, see LICENSE)
 *
 * This copy is adapted for viamd (Dear ImGui >= 1.92, dynamic fonts). Changes with respect to upstream:
 *   - Fonts: no font swapping. Headings are drawn with a scaled size (PushFont(NULL, size)), emphasis / strong /
 *     inline code are distinguished by colour, since the application only loads a single (monospace) font.
 *   - Tables use ImGui::BeginTable, so cells wrap and columns are laid out by ImGui. Column widths are seeded from
 *     a pre-pass over the source so that wide text columns get more room than narrow ones.
 *   - Inline layout: text fragments are flowed onto the current line with correct word wrapping between spans.
 *   - Code blocks are collected and handed to render_code_block(), which derived classes can override to add
 *     syntax colouring, copy buttons, etc.
 *   - Headings report themselves through on_heading() so a host can scroll to them.
 *   - Lists, block quotes, entities, and a few inline HTML tags (<br>) are handled; other HTML is ignored.
 *   - Images are not rendered.
 */

#ifndef IMGUI_MD_H
#define IMGUI_MD_H

#include "md4c.h"
#include "imgui.h"
#include <string>
#include <vector>

struct imgui_md
{
	imgui_md();
	virtual ~imgui_md() {}

	// Parse and draw. Returns 0 on success.
	int print(const char* str, const char* str_end);

protected:

	virtual void BLOCK_DOC(bool);
	virtual void BLOCK_QUOTE(bool);
	virtual void BLOCK_UL(const MD_BLOCK_UL_DETAIL*, bool);
	virtual void BLOCK_OL(const MD_BLOCK_OL_DETAIL*, bool);
	virtual void BLOCK_LI(const MD_BLOCK_LI_DETAIL*, bool);
	virtual void BLOCK_HR(bool);
	virtual void BLOCK_H(const MD_BLOCK_H_DETAIL*, bool);
	virtual void BLOCK_CODE(const MD_BLOCK_CODE_DETAIL*, bool);
	virtual void BLOCK_HTML(bool);
	virtual void BLOCK_P(bool);
	virtual void BLOCK_TABLE(const MD_BLOCK_TABLE_DETAIL*, bool);
	virtual void BLOCK_THEAD(bool);
	virtual void BLOCK_TBODY(bool);
	virtual void BLOCK_TR(bool);
	virtual void BLOCK_TH(const MD_BLOCK_TD_DETAIL*, bool);
	virtual void BLOCK_TD(const MD_BLOCK_TD_DETAIL*, bool);

	virtual void SPAN_EM(bool);
	virtual void SPAN_STRONG(bool);
	virtual void SPAN_A(const MD_SPAN_A_DETAIL*, bool);
	virtual void SPAN_IMG(const MD_SPAN_IMG_DETAIL*, bool);
	virtual void SPAN_CODE(bool);
	virtual void SPAN_DEL(bool);
	virtual void SPAN_U(bool);

	////////////////////////////////////////////////////////////////////////////
	// Hooks for derived classes

	// Draw a fenced / indented code block. 'lang' is the info string (may be empty).
	// Default: a framed, non-scrolling monospace box.
	virtual void render_code_block(const std::string& lang, const std::string& code);

	// Called when a heading has been laid out. 'screen_y' is the y position (screen space) of its top.
	virtual void on_heading(int level, const std::string& text, float screen_y);

	// A link was clicked. 'url' is the raw destination ("#anchor", "https://...").
	virtual void open_url(const std::string& url);

	// Font size multiplier and text colour for the *current* state (m_hlevel, m_is_strong, m_is_code, m_href, ...).
	virtual float  get_text_scale() const;
	virtual ImVec4 get_text_color() const;

	// Optional font override (nullptr -> current font).
	virtual ImFont* get_font() const { return nullptr; }

	// Vertical gap that separates block level elements, in pixels.
	virtual float get_block_gap() const;

	////////////////////////////////////////////////////////////////////////////
	// Current state, readable from the hooks above

	std::string m_href;             // destination of the link we are inside, empty if none
	bool     m_is_underline     = false;
	bool     m_is_strikethrough = false;
	bool     m_is_em            = false;
	bool     m_is_strong        = false;
	bool     m_is_code          = false; // inside an inline code span
	bool     m_is_table_header  = false;
	bool     m_is_image         = false;
	unsigned m_hlevel           = 0;     // 0 - no heading

	int m_code_index = 0; // running index of code blocks in the current print() call, useful for ImGui ids

private:

	int text(MD_TEXTTYPE type, const char* str, const char* str_end);
	int block(MD_BLOCKTYPE type, void* d, bool e);
	int span(MD_SPANTYPE type, void* d, bool e);

	void render_text(const char* str, const char* str_end);
	bool render_entity(const char* str, const char* str_end);
	void block_gap();
	void end_line() { m_inline = false; }

	// Layout state
	bool  m_inline      = false; // the last thing emitted was a text fragment that the next one continues
	bool  m_first_block = true;
	bool  m_skip        = false; // inside a table that is completely clipped

	// Headings
	std::string m_heading_text;
	float       m_heading_y = 0.0f;

	// Code blocks
	bool        m_in_code_block = false;
	std::string m_code_lang;
	std::string m_code_buf;

	// Lists
	struct list_info
	{
		unsigned cur_ol;
		char     delim;
		bool     is_ol;
	};
	std::vector<list_info> m_list_stack;

	// Block quotes: y position where each open quote started
	std::vector<float> m_quote_stack;

	// Tables
	int                                m_table_counter = 0;
	int                                m_table_index   = 0;  // index of the table being drawn (document order)
	bool                               m_in_table      = false;
	int                                m_table_ncols   = 0;
	std::vector<std::vector<float>>    m_table_weights;      // per table, per column

	void prepass_tables(const char* str, const char* str_end);

	MD_PARSER m_md;
};

#endif
