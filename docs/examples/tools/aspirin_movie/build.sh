#!/bin/sh
# Builds the two small C++ helpers into ASPIRIN_WORK. They link VIAMD's own camera code and mdlib, so the movie's
# camera path and script are evaluated as the application does. Needs the main build first (build/lib/libmdlib.a).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../../.." && pwd)"
WORK="${ASPIRIN_WORK:-$HERE/work}"
mkdir -p "$WORK"
INC="-I $ROOT/src -I $ROOT/ext/mdlib/src -I $ROOT/ext/mdlib/ext/simde -I $ROOT/ext/mdlib/ext"
LIBS="-L$ROOT/build/lib -lmdlib -lm -lpthread -ldl"
g++ -std=c++17 -O1 -w -mavx2 $INC "$HERE/camtool.cpp" "$ROOT/src/gfx/camera_utils.cpp" -o "$WORK/camtool" $LIBS
g++ -std=c++17 -O1 -w -mavx2 $INC "$HERE/scripttest.cpp" -o "$WORK/scripttest" $LIBS
echo "built $WORK/camtool and $WORK/scripttest"
echo "check the script of a workspace: $WORK/scripttest <data folder> <workspace.via>"
