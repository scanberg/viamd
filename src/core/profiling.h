#pragma once

#include <core/types.h>
#include <core/string_utils.h>

namespace profiling {
	void initialize();
	void shutdown();

	void push_section(CString sec);
	void pop_section();

	void finish();
	void print_log();

	struct ScopedSection {
		ScopedSection(CString sec) { push_section(sec); }
		~ScopedSection() { pop_section(); }
	};
}  // namespace math