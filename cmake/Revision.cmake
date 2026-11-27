# Writes OUT, a header naming the revisions of viamd (SRC) and of its mdlib (SRC/ext/mdlib) as git has them: the short
# hash, "-dirty" with uncommitted changes to tracked files, "unknown" without git. Run at every build (the target
# viamd_revision); the file is rewritten only when they change, so nothing is recompiled for nothing.
find_package(Git QUIET)

function(revision_of dir out)
    set(r "unknown")
    if(GIT_FOUND)
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${dir}" rev-parse --short=10 HEAD
                        OUTPUT_VARIABLE hash OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE rc ERROR_QUIET)
        if(rc EQUAL 0)
            set(r "${hash}")
            execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${dir}" status --porcelain --untracked-files=no
                            OUTPUT_VARIABLE changes RESULT_VARIABLE rc ERROR_QUIET)
            if(rc EQUAL 0 AND NOT changes STREQUAL "")
                set(r "${r}-dirty")
            endif()
        endif()
    endif()
    set(${out} "${r}" PARENT_SCOPE)
endfunction()

revision_of("${SRC}" viamd_rev)
revision_of("${SRC}/ext/mdlib" mdlib_rev)
set(text "// Generated at build time (cmake/Revision.cmake): the revisions this build is made from\n#pragma once\n#define VIAMD_REVISION \"${viamd_rev}\"\n#define MDLIB_REVISION \"${mdlib_rev}\"\n")
set(old "")
if(EXISTS "${OUT}")
    file(READ "${OUT}" old)
endif()
if(NOT old STREQUAL text)
    file(WRITE "${OUT}" "${text}")
endif()
