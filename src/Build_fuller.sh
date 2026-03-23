#!/bin/bash
# Build_fuller.sh — Fortran 95版フラーレンMDビルドスクリプト
#
# 使用方法:
#   src/Build_fuller.sh          # Serial をビルド (デフォルト)
#   src/Build_fuller.sh serial   # Serial のみ
#   src/Build_fuller.sh clean    # ビルド成果物を削除

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC_DIR="$ROOT_DIR/src"
BIN_DIR="$ROOT_DIR/bin"

# Fortran コンパイラ検出
if [ -n "$FC" ]; then
    FCOMP="$FC"
elif command -v gfortran-15 &>/dev/null; then
    FCOMP=gfortran-15
elif command -v gfortran-14 &>/dev/null; then
    FCOMP=gfortran-14
elif command -v gfortran-13 &>/dev/null; then
    FCOMP=gfortran-13
elif command -v gfortran-12 &>/dev/null; then
    FCOMP=gfortran-12
elif command -v gfortran &>/dev/null; then
    FCOMP=gfortran
else
    echo "ERROR: gfortran が見つかりません。"
    echo "  macOS: brew install gcc"
    echo "  Linux: sudo apt install gfortran"
    exit 1
fi

echo "Fortran compiler: $FCOMP"
$FCOMP --version 2>/dev/null | head -1

MODE=${1:-serial}

build_serial() {
    echo ""
    echo "=== Building Serial ==="
    mkdir -p "$BIN_DIR"

    # [1] LJ コア版 Serial
    if [ -f "$SRC_DIR/fuller_LJ_npt_md_core_serial.f90" ]; then
        echo "  [1] fuller_LJ_core_serial"
        $FCOMP -O3 \
            -o "$BIN_DIR/fuller_LJ_core_serial" \
            "$SRC_DIR/fuller_LJ_npt_md_core_serial.f90"
    fi

    echo ""
    echo "=== Build complete ==="
    echo "Executables in: $BIN_DIR/"
    ls -la "$BIN_DIR/"
}

do_clean() {
    echo "Cleaning..."
    rm -rf "$BIN_DIR"
    echo "Done."
}

case "$MODE" in
    serial)  build_serial ;;
    clean)   do_clean ;;
    *)
        echo "Usage: $0 [serial|clean]"
        exit 1
        ;;
esac
