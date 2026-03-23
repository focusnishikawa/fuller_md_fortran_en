#!/bin/bash
# Build_fuller.sh — Fortran 95版フラーレンMDビルドスクリプト
#
# 使用方法:
#   src/Build_fuller.sh          # Serial + OpenMP をビルド (デフォルト)
#   src/Build_fuller.sh serial   # Serial のみ
#   src/Build_fuller.sh omp      # OpenMP のみ
#   src/Build_fuller.sh all      # Serial + OpenMP
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

MODE=${1:-all}

build_serial() {
    echo ""
    echo "=== Building Serial ==="
    mkdir -p "$BIN_DIR"

    # [1] LJ コア版 Serial (pure)
    if [ -f "$SRC_DIR/fuller_LJ_npt_md_core_serial.f90" ]; then
        echo "  [1] fuller_LJ_core_serial_pure"
        $FCOMP -O3 \
            -o "$BIN_DIR/fuller_LJ_core_serial_pure" \
            "$SRC_DIR/fuller_LJ_npt_md_core_serial.f90"
    fi

    # [2] LJ コア版 Serial (OMP/ACC統合ソースのSerial)
    if [ -f "$SRC_DIR/fuller_LJ_npt_md_core_serial_omp_acc.f90" ]; then
        echo "  [2] fuller_LJ_core_serial"
        $FCOMP -O3 \
            -o "$BIN_DIR/fuller_LJ_core_serial" \
            "$SRC_DIR/fuller_LJ_npt_md_core_serial_omp_acc.f90"
    fi

    # [3] LJ フル版 Serial
    if [ -f "$SRC_DIR/fuller_LJ_npt_md_serial_omp_acc.f90" ]; then
        echo "  [3] fuller_LJ_npt_md_serial"
        $FCOMP -O3 \
            -o "$BIN_DIR/fuller_LJ_npt_md_serial" \
            "$SRC_DIR/fuller_LJ_npt_md_serial_omp_acc.f90"
    fi

    # [4] 分子力学版 Serial
    if [ -f "$SRC_DIR/fuller_LJ_npt_mmmd_serial_omp_acc.f90" ]; then
        echo "  [4] fuller_LJ_npt_mmmd_serial"
        $FCOMP -O3 \
            -o "$BIN_DIR/fuller_LJ_npt_mmmd_serial" \
            "$SRC_DIR/fuller_LJ_npt_mmmd_serial_omp_acc.f90"
    fi

    # [5] AIREBO版 Serial
    if [ -f "$SRC_DIR/fuller_airebo_npt_md_serial_omp_acc.f90" ]; then
        echo "  [5] fuller_airebo_npt_md_serial"
        $FCOMP -O3 \
            -o "$BIN_DIR/fuller_airebo_npt_md_serial" \
            "$SRC_DIR/fuller_airebo_npt_md_serial_omp_acc.f90"
    fi

    echo ""
    echo "=== Serial build complete ==="
}

build_omp() {
    echo ""
    echo "=== Building OpenMP ==="
    mkdir -p "$BIN_DIR"

    # [2] LJ コア版 OpenMP
    if [ -f "$SRC_DIR/fuller_LJ_npt_md_core_serial_omp_acc.f90" ]; then
        echo "  [2] fuller_LJ_core_omp"
        $FCOMP -O3 -fopenmp \
            -o "$BIN_DIR/fuller_LJ_core_omp" \
            "$SRC_DIR/fuller_LJ_npt_md_core_serial_omp_acc.f90"
    fi

    # [3] LJ フル版 OpenMP
    if [ -f "$SRC_DIR/fuller_LJ_npt_md_serial_omp_acc.f90" ]; then
        echo "  [3] fuller_LJ_npt_md_omp"
        $FCOMP -O3 -fopenmp \
            -o "$BIN_DIR/fuller_LJ_npt_md_omp" \
            "$SRC_DIR/fuller_LJ_npt_md_serial_omp_acc.f90"
    fi

    # [4] 分子力学版 OpenMP
    if [ -f "$SRC_DIR/fuller_LJ_npt_mmmd_serial_omp_acc.f90" ]; then
        echo "  [4] fuller_LJ_npt_mmmd_omp"
        $FCOMP -O3 -fopenmp \
            -o "$BIN_DIR/fuller_LJ_npt_mmmd_omp" \
            "$SRC_DIR/fuller_LJ_npt_mmmd_serial_omp_acc.f90"
    fi

    # [5] AIREBO版 OpenMP
    if [ -f "$SRC_DIR/fuller_airebo_npt_md_serial_omp_acc.f90" ]; then
        echo "  [5] fuller_airebo_npt_md_omp"
        $FCOMP -O3 -fopenmp \
            -o "$BIN_DIR/fuller_airebo_npt_md_omp" \
            "$SRC_DIR/fuller_airebo_npt_md_serial_omp_acc.f90"
    fi

    echo ""
    echo "=== OpenMP build complete ==="
}

do_clean() {
    echo "Cleaning..."
    rm -rf "$BIN_DIR"
    rm -f "$SRC_DIR"/*.mod "$SRC_DIR"/*.o
    echo "Done."
}

case "$MODE" in
    serial)  build_serial ;;
    omp)     build_omp ;;
    all)     build_serial; build_omp ;;
    clean)   do_clean ;;
    *)
        echo "Usage: $0 [serial|omp|all|clean]"
        exit 1
        ;;
esac

echo ""
echo "Executables in: $BIN_DIR/"
ls -la "$BIN_DIR/" 2>/dev/null || echo "(empty)"
