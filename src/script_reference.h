#pragma once

// In-app viewer for the mdlib script reference (ext/mdlib/docs/script_reference.md).
//
// The reference is plain Markdown that lives next to the implementation. This viewer parses its outline
// (headings, procedure metadata comments, alias anchors), shows one section ("page") at a time with a navigable
// sidebar, history, and search, and renders the Markdown with imgui_md.
//
// The document is expected to follow the conventions listed at the top of script_reference.md:
//   - `##` headings are the pages, `###` headings are entries within a page.
//   - Every procedure has a `### name` heading followed by `<!-- proc name=... aliases=a,b category=... -->`
//     and optionally `<a id="alias"></a>` anchors.
//   - Code blocks tagged `mdscript` are examples that can be copied or inserted into the script editor.

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

struct ScriptReference {
    struct Section {
        std::string title;                  // plain text of the heading
        std::string slug;                   // GitHub style anchor of the heading
        std::string proc;                   // procedure name from <!-- proc name= -->, empty if this is not a procedure
        std::string category;               // e.g. "selector.atom"
        std::vector<std::string> aliases;   // other names that lead here
        std::string body_lower;             // lower-cased text of the section, for searching
        int level = 0;                      // heading level, 1..6
        int page = 0;
        int heading_index = 0;              // index among the headings of its page, in document order
    };

    struct Page {
        std::string title;
        size_t begin = 0;                   // byte range in the document
        size_t end = 0;
        std::vector<int> sections;          // indices into sections, in document order
    };

    ScriptReference();

    // Parses the outline. The text is copied.
    void set_document(const char* markdown, size_t size);
    bool has_document() const { return !m_doc.empty(); }

    // Draw the viewer into the current window (toolbar, sidebar, content).
    void draw();

    // Convenience: draw as its own window. Returns false if the window is closed.
    void draw_window(bool* p_open, const char* title = "Script Reference");

    // Navigation. show() returns false if there is no entry with that name or anchor, in which case nothing changes.
    bool show(const char* topic);
    void search(const char* query);     // fills the search box, and shows the results
    void focus_search() { m_focus_search = true; }

    // Host callbacks. Both are optional.
    std::function<void(const std::string& code)> on_insert_code;   // "Insert" button of example blocks
    std::function<void(const std::string& url)>  on_open_url;      // external links

    // Read access to the parsed outline (used by tools and tests).
    const std::vector<Section>& sections() const { return m_sections; }
    const std::vector<Page>&    pages()    const { return m_pages; }
    int find_anchor(const std::string& name) const;   // -1 if not found
    int current_section() const { return m_current; }
    int num_headings_rendered() const { return m_headings_rendered; }
    const std::string& document() const { return m_doc; }

private:
    friend struct ScriptReferenceRenderer;

    void navigate(int section, bool push_history = true);
    void follow_link(const std::string& url);
    void history_step(int dir);
    void update_search_results();

    void draw_toolbar();
    void draw_sidebar();
    void draw_content();

    std::string m_doc;
    std::vector<Page> m_pages;
    std::vector<Section> m_sections;
    std::unordered_map<std::string, int> m_anchors;        // anchor / procedure name / alias -> section
    std::unordered_map<std::string, int> m_anchors_lower;

    // View state
    int  m_page    = 0;
    int  m_current = 0;                 // section that was navigated to last
    int  m_target_heading = -1;         // heading (on the current page) to scroll to, -1 for none
    bool m_scroll_top   = false;
    int  m_scroll_settle = 0;           // frames the scroll target stays armed
    bool m_reveal_nav   = false;        // scroll the sidebar to the current entry
    bool m_focus_search = false;
    int  m_headings_rendered = 0;       // headings seen while drawing the current page (for scroll targets)

    std::vector<int> m_history;
    int m_history_pos = -1;

    char  m_query[128] = "";
    std::string m_query_applied;
    struct Hit { int section; int score; };
    std::vector<Hit> m_hits;

    float m_nav_width = 0.0f;           // sidebar width in pixels (0 -> derive from font size)

    // Feedback for the copy / insert buttons of code blocks
    int    m_feedback_id = -1;
    double m_feedback_until = 0.0;
    const char* m_feedback_text = "";
};

// The reference that ships with the application (embedded in the executable). Created on first use.
ScriptReference& script_reference();
