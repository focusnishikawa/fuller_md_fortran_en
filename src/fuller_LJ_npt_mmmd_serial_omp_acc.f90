! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2025, Takeshi Nishikawa
!===========================================================================
!  fuller_LJ_npt_mmmd_serial_omp_acc.f90
!  Fullerene Crystal NPT-MD (Molecular Mechanics Force Field + LJ, Fortran 95)
!
!  V_total = V_bond + V_angle + V_dihedral + V_improper + V_LJ + V_Coulomb
!
!  コンパイル:
!    gfortran -O3 -o fuller_LJ_npt_mmmd_serial fuller_LJ_npt_mmmd_serial_omp_acc.f90
!    gfortran -O3 -fopenmp -o fuller_LJ_npt_mmmd_omp fuller_LJ_npt_mmmd_serial_omp_acc.f90
!
!  実行時オプション: --cell=nc --temp=K --pres=GPa --step=N --dt=fs --seed=n
!    --coldstart=N --warmup=N --from=step --to=step --mon=N
!    --ovito=N --restart=N --resfile=path
!    --ff_kb=kcal/mol --ff_kth=kcal/mol --ff_v2=kcal/mol --ff_kimp=kcal/mol
!  単位系: A, amu, eV, fs, K, GPa
!===========================================================================
module fuller_mmmd_mod
  implicit none
  double precision, parameter :: CONV=9.64853321d-3, kB=8.617333262d-5
  double precision, parameter :: eV2GPa=160.21766208d0, eV2kcalmol=23.06054783d0
  double precision, parameter :: kcal2eV=1.0d0/eV2kcalmol
  double precision, parameter :: PI_=3.14159265358979323846d0
  double precision, parameter :: sigma_LJ=3.431d0, eps_LJ=2.635d-3
  double precision, parameter :: RCUT=3.0d0*sigma_LJ, RCUT2=RCUT*RCUT
  double precision, parameter :: sig2_LJ=sigma_LJ*sigma_LJ, mC=12.011d0
  double precision, parameter :: sr_v=1.0d0/3.0d0, sr2_v=sr_v*sr_v
  double precision, parameter :: sr6_v=sr2_v*sr2_v*sr2_v
  double precision, parameter :: VSHFT=4.0d0*eps_LJ*(sr6_v*sr6_v-sr6_v)
  double precision, parameter :: COULOMB_K=14.3996d0
  double precision, parameter :: COUL_RCUT=RCUT, COUL_RCUT2=COUL_RCUT*COUL_RCUT
  integer, parameter :: MAX_NATOM=84, MAX_LJ_NEIGH=400

  type :: NPTState
    double precision :: xi, Q, Vg(9), W, Pe, Tt
    integer :: Nf
  end type

  ! ═══ トポロジー構造体 (GPU-friendly flat arrays) ═══
  type :: FlatTopology
    integer :: Nb, Nang, Ndih, Nimp
    integer, allocatable :: b_i0(:), b_i1(:)
    double precision, allocatable :: b_kb(:), b_r0(:)
    integer, allocatable :: ang_i0(:), ang_i1(:), ang_i2(:)
    double precision, allocatable :: ang_kth(:), ang_th0(:)
    integer, allocatable :: dih_i0(:), dih_i1(:), dih_i2(:), dih_i3(:)
    double precision, allocatable :: dih_Vn(:), dih_gamma(:)
    integer, allocatable :: dih_mult(:)
    integer, allocatable :: imp_i0(:), imp_i1(:), imp_i2(:), imp_i3(:)
    double precision, allocatable :: imp_ki(:), imp_gamma(:)
    integer, allocatable :: mol_id(:)
  end type

  ! ═══ 力の結果構造体 ═══
  type :: ForceResult
    double precision :: Eb, Ea, Ed, Ei, Elj, Ecoul, Etot
  end type

contains

  ! ═══ 行列式 ═══
  double precision function mat_det9(h)
    double precision, intent(in) :: h(9)
    mat_det9 = h(1)*(h(5)*h(9)-h(6)*h(8)) &
             - h(2)*(h(4)*h(9)-h(6)*h(7)) &
             + h(3)*(h(4)*h(8)-h(5)*h(7))
  end function

  ! ═══ 逆行列 ═══
  subroutine mat_inv9(h, hi)
    double precision, intent(in) :: h(9)
    double precision, intent(out) :: hi(9)
    double precision :: d, id
    d = mat_det9(h); id = 1.0d0/d
    hi(1) = id*(h(5)*h(9)-h(6)*h(8)); hi(2) = id*(h(3)*h(8)-h(2)*h(9))
    hi(3) = id*(h(2)*h(6)-h(3)*h(5)); hi(4) = id*(h(6)*h(7)-h(4)*h(9))
    hi(5) = id*(h(1)*h(9)-h(3)*h(7)); hi(6) = id*(h(3)*h(4)-h(1)*h(6))
    hi(7) = id*(h(4)*h(8)-h(5)*h(7)); hi(8) = id*(h(2)*h(7)-h(1)*h(8))
    hi(9) = id*(h(1)*h(5)-h(2)*h(4))
  end subroutine

  ! ═══ 最小イメージ変換 ═══
  subroutine mimg_flat(dx, dy, dz, hi, h)
    double precision, intent(inout) :: dx, dy, dz
    double precision, intent(in) :: hi(9), h(9)
    double precision :: s0, s1, s2
    s0 = hi(1)*dx + hi(2)*dy + hi(3)*dz
    s1 = hi(4)*dx + hi(5)*dy + hi(6)*dz
    s2 = hi(7)*dx + hi(8)*dy + hi(9)*dz
    s0 = s0 - anint(s0); s1 = s1 - anint(s1); s2 = s2 - anint(s2)
    dx = h(1)*s0 + h(2)*s1 + h(3)*s2
    dy = h(4)*s0 + h(5)*s1 + h(6)*s2
    dz = h(7)*s0 + h(8)*s1 + h(9)*s2
  end subroutine

  ! ═══ PBC適用 ═══
  subroutine apply_pbc(pos, h, hi, Na)
    double precision, intent(inout) :: pos(:)
    double precision, intent(in) :: h(9), hi(9)
    integer, intent(in) :: Na
    double precision :: px, py, pz, s0, s1, s2
    integer :: i, idx
    !$OMP PARALLEL DO PRIVATE(i,idx,px,py,pz,s0,s1,s2)
    do i = 1, Na
      idx = (i-1)*3
      px = pos(idx+1); py = pos(idx+2); pz = pos(idx+3)
      s0 = hi(1)*px + hi(2)*py + hi(3)*pz
      s1 = hi(4)*px + hi(5)*py + hi(6)*pz
      s2 = hi(7)*px + hi(8)*py + hi(9)*pz
      s0 = s0 - floor(s0); s1 = s1 - floor(s1); s2 = s2 - floor(s2)
      pos(idx+1) = h(1)*s0 + h(2)*s1 + h(3)*s2
      pos(idx+2) = h(4)*s0 + h(5)*s1 + h(6)*s2
      pos(idx+3) = h(7)*s0 + h(8)*s1 + h(9)*s2
    end do
    !$OMP END PARALLEL DO
  end subroutine

  ! ═══ 文字列分割ユーティリティ ═══
  subroutine split_line(line, words, nw)
    character(len=*), intent(in) :: line
    character(len=32), intent(out) :: words(100)
    integer, intent(out) :: nw
    integer :: i, ln, wstart
    logical :: in_word
    ln = len_trim(line)
    nw = 0; in_word = .false.
    do i = 1, ln
      if (line(i:i) /= ' ' .and. line(i:i) /= char(9)) then
        if (.not. in_word) then
          in_word = .true.; wstart = i
        end if
      else
        if (in_word) then
          in_word = .false.; nw = nw + 1
          if (nw <= 100) words(nw) = line(wstart:i-1)
        end if
      end if
    end do
    if (in_word) then
      nw = nw + 1
      if (nw <= 100) words(nw) = line(wstart:ln)
    end if
  end subroutine

  ! ═══ cc1ファイル読み込み (ボンド接続情報付き) ═══
  subroutine load_cc1_mmmd(path, coords, natom, Rmol, Dmol, &
                           bond_i, bond_j, nbonds)
    character(len=*), intent(in) :: path
    double precision, intent(out) :: coords(MAX_NATOM, 3), Rmol, Dmol
    integer, intent(out) :: natom, nbonds
    integer, intent(out) :: bond_i(4096), bond_j(4096)
    integer :: adj(MAX_NATOM, 20), nadj(MAX_NATOM)
    double precision :: cm(3), r2, r, dx, dy, dz, d, x, y, z
    integer :: i, j, k, ios, idx_val, fl, bval, nw, ja
    character(len=1024) :: line
    character(len=32) :: words(100)
    logical :: found

    open(unit=20, file=path, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: cannot open ', trim(path); stop 1
    end if
    read(20, *) natom
    if (natom > MAX_NATOM) then
      write(*,*) 'Error: natom > MAX_NATOM'; stop 1
    end if

    nadj = 0
    do i = 1, natom
      read(20, '(A)', iostat=ios) line
      do while (len_trim(line) == 0 .and. ios == 0)
        read(20, '(A)', iostat=ios) line
      end do
      ! フィールド分割: elem idx x y z fl [bond1 bond2 ...]
      call split_line(line, words, nw)
      if (nw < 6) cycle
      read(words(3), *, iostat=ios) x
      read(words(4), *, iostat=ios) y
      read(words(5), *, iostat=ios) z
      coords(i, 1) = x; coords(i, 2) = y; coords(i, 3) = z
      ! 7番目以降はボンド接続情報 (1-indexed)
      do k = 7, nw
        read(words(k), *, iostat=ios) bval
        if (ios == 0 .and. bval >= 1 .and. bval <= natom) then
          nadj(i) = nadj(i) + 1
          if (nadj(i) <= 20) adj(i, nadj(i)) = bval
        end if
      end do
    end do
    close(20)

    ! 重心補正
    cm = 0.0d0
    do i = 1, natom
      cm(1) = cm(1) + coords(i, 1)
      cm(2) = cm(2) + coords(i, 2)
      cm(3) = cm(3) + coords(i, 3)
    end do
    cm = cm / dble(natom)
    do i = 1, natom
      coords(i, 1) = coords(i, 1) - cm(1)
      coords(i, 2) = coords(i, 2) - cm(2)
      coords(i, 3) = coords(i, 3) - cm(3)
    end do

    ! Rmol, Dmol計算
    Rmol = 0.0d0; Dmol = 0.0d0
    do i = 1, natom
      r2 = coords(i,1)**2 + coords(i,2)**2 + coords(i,3)**2; r = sqrt(r2)
      if (r > Rmol) Rmol = r
      do j = i+1, natom
        dx = coords(i,1) - coords(j,1); dy = coords(i,2) - coords(j,2)
        dz = coords(i,3) - coords(j,3)
        d = sqrt(dx*dx + dy*dy + dz*dz); if (d > Dmol) Dmol = d
      end do
    end do

    ! ボンドリスト構築 (重複排除: i < j のみ)
    nbonds = 0
    do i = 1, natom
      do k = 1, nadj(i)
        ja = adj(i, k)
        if (i < ja) then
          found = .false.
          do j = 1, nbonds
            if (bond_i(j) == i .and. bond_j(j) == ja) then
              found = .true.; exit
            end if
          end do
          if (.not. found) then
            nbonds = nbonds + 1
            bond_i(nbonds) = i; bond_j(nbonds) = ja
          end if
        end if
      end do
    end do
  end subroutine

  ! ═══ C60座標の生成 (cc1ファイルなしの場合のフォールバック) ═══
  subroutine generate_c60_mmmd(coords, natom, Rmol, Dmol, &
                               bond_i, bond_j, nbonds)
    double precision, intent(out) :: coords(MAX_NATOM, 3), Rmol, Dmol
    integer, intent(out) :: natom, nbonds
    integer, intent(out) :: bond_i(4096), bond_j(4096)
    double precision :: phi, tmp(60, 3), cm(3), r2, r, dx, dy, dz, d
    integer :: n, p, s1, s2, s3, cyc(3, 3), signs(2), i, j

    natom = 60; phi = (1.0d0 + sqrt(5.0d0)) / 2.0d0
    signs(1) = -1; signs(2) = 1
    cyc(1,1)=1; cyc(1,2)=2; cyc(1,3)=3
    cyc(2,1)=2; cyc(2,2)=3; cyc(2,3)=1
    cyc(3,1)=3; cyc(3,2)=1; cyc(3,3)=2
    tmp = 0.0d0; n = 0
    do p = 1, 3; do s2 = 1, 2; do s3 = 1, 2
      n = n + 1; tmp(n, cyc(p,1)) = 0.0d0
      tmp(n, cyc(p,2)) = dble(signs(s2))
      tmp(n, cyc(p,3)) = dble(signs(s3)) * 3.0d0 * phi
    end do; end do; end do
    do p = 1, 3; do s1 = 1, 2; do s2 = 1, 2; do s3 = 1, 2
      n = n + 1
      tmp(n, cyc(p,1)) = dble(signs(s1)) * 2.0d0
      tmp(n, cyc(p,2)) = dble(signs(s2)) * (1 + 2*phi)
      tmp(n, cyc(p,3)) = dble(signs(s3)) * phi
    end do; end do; end do; end do
    do p = 1, 3; do s1 = 1, 2; do s2 = 1, 2; do s3 = 1, 2
      n = n + 1
      tmp(n, cyc(p,1)) = dble(signs(s1))
      tmp(n, cyc(p,2)) = dble(signs(s2)) * (2 + phi)
      tmp(n, cyc(p,3)) = dble(signs(s3)) * 2 * phi
    end do; end do; end do; end do
    cm = 0
    do n = 1, 60; cm = cm + tmp(n, :); end do
    cm = cm / 60
    do n = 1, 60; tmp(n, :) = (tmp(n, :) - cm) * 0.72d0; end do
    Rmol = 0; Dmol = 0
    do i = 1, 60
      r2 = tmp(i,1)**2 + tmp(i,2)**2 + tmp(i,3)**2; r = sqrt(r2)
      if (r > Rmol) Rmol = r
      do j = i+1, 60
        dx = tmp(i,1)-tmp(j,1); dy = tmp(i,2)-tmp(j,2); dz = tmp(i,3)-tmp(j,3)
        d = sqrt(dx*dx + dy*dy + dz*dz); if (d > Dmol) Dmol = d
      end do
    end do
    coords(1:60, :) = tmp
    ! C60のボンドリスト: 距離ベースで生成 (1.8A以内)
    nbonds = 0
    do i = 1, 60
      do j = i+1, 60
        dx = coords(i,1)-coords(j,1); dy = coords(i,2)-coords(j,2)
        dz = coords(i,3)-coords(j,3)
        d = sqrt(dx*dx + dy*dy + dz*dz)
        if (d < 1.8d0) then
          nbonds = nbonds + 1
          bond_i(nbonds) = i; bond_j(nbonds) = j
        end if
      end do
    end do
  end subroutine

  ! ═══ FCC結晶構築 ═══
  function make_fcc(a, nc, pos, h) result(Nmol)
    double precision, intent(in) :: a
    integer, intent(in) :: nc
    double precision, intent(out) :: pos(:), h(9)
    integer :: Nmol, ix, iy, iz, b, idx
    double precision :: bas(4, 3)
    bas(1,:) = (/0d0, 0d0, 0d0/)
    bas(2,:) = (/.5d0*a, .5d0*a, 0d0/)
    bas(3,:) = (/.5d0*a, 0d0, .5d0*a/)
    bas(4,:) = (/0d0, .5d0*a, .5d0*a/)
    Nmol = 0
    do ix = 0, nc-1; do iy = 0, nc-1; do iz = 0, nc-1; do b = 1, 4
      Nmol = Nmol + 1; idx = (Nmol-1)*3
      pos(idx+1) = a*dble(ix) + bas(b,1)
      pos(idx+2) = a*dble(iy) + bas(b,2)
      pos(idx+3) = a*dble(iz) + bas(b,3)
    end do; end do; end do; end do
    h = 0; h(1) = dble(nc)*a; h(5) = dble(nc)*a; h(9) = dble(nc)*a
  end function

  ! ═══ 二面角計算 (トポロジー構築用) ═══
  double precision function compute_phi0(cc, natom_in, i, j, k, l)
    double precision, intent(in) :: cc(:)
    integer, intent(in) :: natom_in, i, j, k, l
    double precision :: b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z
    double precision :: mx, my, mz, nx, ny, nz, mm, nn, b2len
    double precision :: cosphi, sinphi
    ! i,j,k,l は1-indexed
    b1x = cc((j-1)*3+1) - cc((i-1)*3+1)
    b1y = cc((j-1)*3+2) - cc((i-1)*3+2)
    b1z = cc((j-1)*3+3) - cc((i-1)*3+3)
    b2x = cc((k-1)*3+1) - cc((j-1)*3+1)
    b2y = cc((k-1)*3+2) - cc((j-1)*3+2)
    b2z = cc((k-1)*3+3) - cc((j-1)*3+3)
    b3x = cc((l-1)*3+1) - cc((k-1)*3+1)
    b3y = cc((l-1)*3+2) - cc((k-1)*3+2)
    b3z = cc((l-1)*3+3) - cc((k-1)*3+3)
    mx = b1y*b2z - b1z*b2y; my = b1z*b2x - b1x*b2z; mz = b1x*b2y - b1y*b2x
    nx = b2y*b3z - b2z*b3y; ny = b2z*b3x - b2x*b3z; nz = b2x*b3y - b2y*b3x
    mm = mx*mx + my*my + mz*mz; nn = nx*nx + ny*ny + nz*nz
    if (mm < 1d-20 .or. nn < 1d-20) then
      compute_phi0 = 0.0d0; return
    end if
    b2len = sqrt(b2x*b2x + b2y*b2y + b2z*b2z)
    cosphi = (mx*nx + my*ny + mz*nz) / sqrt(mm*nn)
    if (cosphi > 1.0d0) cosphi = 1.0d0
    if (cosphi < -1.0d0) cosphi = -1.0d0
    sinphi = (mx*b3x + my*b3y + mz*b3z) * b2len / sqrt(mm*nn)
    compute_phi0 = atan2(sinphi, cosphi)
  end function

  ! ═══ フラットトポロジー構築 ═══
  subroutine build_flat_topology(ft, mol_coords, natom, &
       bond_i_mol, bond_j_mol, nb_mol_in, Nmol, &
       kb, kth, v2_dih, k_imp)
    type(FlatTopology), intent(out) :: ft
    double precision, intent(in) :: mol_coords(MAX_NATOM, 3)
    integer, intent(in) :: natom, nb_mol_in, Nmol
    integer, intent(in) :: bond_i_mol(4096), bond_j_mol(4096)
    double precision, intent(in) :: kb, kth, v2_dih, k_imp

    ! 隣接リスト
    integer :: adj(MAX_NATOM, 20), nadj(MAX_NATOM)
    ! per-moleculeトポロジー
    integer :: bnd_a0(4096), bnd_a1(4096), nb_mol
    double precision :: bnd_r0(4096)
    integer :: ang_a0(8192), ang_a1(8192), ang_a2(8192), nang_mol
    double precision :: ang_th0_tmp(8192)
    integer :: dih_a0(16384), dih_a1(16384), dih_a2(16384), dih_a3(16384), ndih_mol
    double precision :: dih_gamma_tmp(16384)
    integer :: imp_a0(4096), imp_a1(4096), imp_a2(4096), imp_a3(4096), nimp_mol
    double precision :: imp_gamma_tmp(4096)
    double precision, allocatable :: cc(:)
    double precision :: dx, dy, dz
    double precision :: rji(3), rjk(3), dji, djk, costh, phi0, psi0
    integer :: i, j, k, l, bi, bj, bk, m, b, a, d, p, ii, kk
    integer :: off, ob, oa, od, oi

    ! 隣接リスト構築
    nadj = 0
    do b = 1, nb_mol_in
      bi = bond_i_mol(b); bj = bond_j_mol(b)
      nadj(bi) = nadj(bi) + 1; adj(bi, nadj(bi)) = bj
      nadj(bj) = nadj(bj) + 1; adj(bj, nadj(bj)) = bi
    end do

    ! Flat coords for phi0
    allocate(cc(natom*3))
    do i = 1, natom
      cc((i-1)*3+1) = mol_coords(i,1)
      cc((i-1)*3+2) = mol_coords(i,2)
      cc((i-1)*3+3) = mol_coords(i,3)
    end do

    ! ── ボンド ──
    nb_mol = nb_mol_in
    do b = 1, nb_mol
      bnd_a0(b) = bond_i_mol(b); bnd_a1(b) = bond_j_mol(b)
      dx = mol_coords(bond_i_mol(b),1) - mol_coords(bond_j_mol(b),1)
      dy = mol_coords(bond_i_mol(b),2) - mol_coords(bond_j_mol(b),2)
      dz = mol_coords(bond_i_mol(b),3) - mol_coords(bond_j_mol(b),3)
      bnd_r0(b) = sqrt(dx*dx + dy*dy + dz*dz)
    end do

    ! ── 角度 ──
    nang_mol = 0
    do j = 1, natom
      do ii = 1, nadj(j)
        do kk = ii+1, nadj(j)
          i = adj(j, ii); k = adj(j, kk)
          rji(1) = mol_coords(i,1) - mol_coords(j,1)
          rji(2) = mol_coords(i,2) - mol_coords(j,2)
          rji(3) = mol_coords(i,3) - mol_coords(j,3)
          rjk(1) = mol_coords(k,1) - mol_coords(j,1)
          rjk(2) = mol_coords(k,2) - mol_coords(j,2)
          rjk(3) = mol_coords(k,3) - mol_coords(j,3)
          dji = sqrt(rji(1)**2 + rji(2)**2 + rji(3)**2)
          djk = sqrt(rjk(1)**2 + rjk(2)**2 + rjk(3)**2)
          costh = (rji(1)*rjk(1) + rji(2)*rjk(2) + rji(3)*rjk(3)) / (dji*djk)
          if (costh > 1.0d0) costh = 1.0d0
          if (costh < -1.0d0) costh = -1.0d0
          nang_mol = nang_mol + 1
          ang_a0(nang_mol) = i; ang_a1(nang_mol) = j; ang_a2(nang_mol) = k
          ang_th0_tmp(nang_mol) = acos(costh)
        end do
      end do
    end do

    ! ── 二面角 ──
    ndih_mol = 0
    do b = 1, nb_mol
      bj = bond_i_mol(b); bk = bond_j_mol(b)
      do ii = 1, nadj(bj)
        i = adj(bj, ii)
        if (i == bk) cycle
        do kk = 1, nadj(bk)
          l = adj(bk, kk)
          if (l == bj .or. l == i) cycle
          phi0 = compute_phi0(cc, natom, i, bj, bk, l)
          ndih_mol = ndih_mol + 1
          dih_a0(ndih_mol) = i; dih_a1(ndih_mol) = bj
          dih_a2(ndih_mol) = bk; dih_a3(ndih_mol) = l
          dih_gamma_tmp(ndih_mol) = 2.0d0 * phi0 + PI_
        end do
      end do
    end do

    ! ── 不適切二面角 (3つの隣接原子を持つ原子) ──
    nimp_mol = 0
    do i = 1, natom
      if (nadj(i) == 3) then
        j = adj(i, 1); k = adj(i, 2); l = adj(i, 3)
        psi0 = compute_phi0(cc, natom, j, i, k, l)
        nimp_mol = nimp_mol + 1
        imp_a0(nimp_mol) = i; imp_a1(nimp_mol) = j
        imp_a2(nimp_mol) = k; imp_a3(nimp_mol) = l
        imp_gamma_tmp(nimp_mol) = 2.0d0 * psi0 + PI_
      end if
    end do

    deallocate(cc)

    ! ── グローバル配列確保 ──
    ft%Nb = nb_mol * Nmol
    ft%Nang = nang_mol * Nmol
    ft%Ndih = ndih_mol * Nmol
    ft%Nimp = nimp_mol * Nmol

    allocate(ft%b_i0(ft%Nb), ft%b_i1(ft%Nb))
    allocate(ft%b_kb(ft%Nb), ft%b_r0(ft%Nb))
    allocate(ft%ang_i0(ft%Nang), ft%ang_i1(ft%Nang), ft%ang_i2(ft%Nang))
    allocate(ft%ang_kth(ft%Nang), ft%ang_th0(ft%Nang))
    allocate(ft%dih_i0(ft%Ndih), ft%dih_i1(ft%Ndih))
    allocate(ft%dih_i2(ft%Ndih), ft%dih_i3(ft%Ndih))
    allocate(ft%dih_Vn(ft%Ndih), ft%dih_mult(ft%Ndih), ft%dih_gamma(ft%Ndih))
    allocate(ft%imp_i0(ft%Nimp), ft%imp_i1(ft%Nimp))
    allocate(ft%imp_i2(ft%Nimp), ft%imp_i3(ft%Nimp))
    allocate(ft%imp_ki(ft%Nimp), ft%imp_gamma(ft%Nimp))
    allocate(ft%mol_id(Nmol * natom))

    ! ── 全分子に複製 (1-indexed) ──
    do m = 1, Nmol
      off = (m-1) * natom
      ob = (m-1) * nb_mol
      oa = (m-1) * nang_mol
      od = (m-1) * ndih_mol
      oi = (m-1) * nimp_mol
      do a = 1, natom
        ft%mol_id(off + a) = m
      end do
      do b = 1, nb_mol
        ft%b_i0(ob+b) = bnd_a0(b) + off
        ft%b_i1(ob+b) = bnd_a1(b) + off
        ft%b_kb(ob+b) = kb; ft%b_r0(ob+b) = bnd_r0(b)
      end do
      do a = 1, nang_mol
        ft%ang_i0(oa+a) = ang_a0(a) + off
        ft%ang_i1(oa+a) = ang_a1(a) + off
        ft%ang_i2(oa+a) = ang_a2(a) + off
        ft%ang_kth(oa+a) = kth
        ft%ang_th0(oa+a) = ang_th0_tmp(a)
      end do
      do d = 1, ndih_mol
        ft%dih_i0(od+d) = dih_a0(d) + off
        ft%dih_i1(od+d) = dih_a1(d) + off
        ft%dih_i2(od+d) = dih_a2(d) + off
        ft%dih_i3(od+d) = dih_a3(d) + off
        ft%dih_Vn(od+d) = v2_dih
        ft%dih_mult(od+d) = 2
        ft%dih_gamma(od+d) = dih_gamma_tmp(d)
      end do
      do p = 1, nimp_mol
        ft%imp_i0(oi+p) = imp_a0(p) + off
        ft%imp_i1(oi+p) = imp_a1(p) + off
        ft%imp_i2(oi+p) = imp_a2(p) + off
        ft%imp_i3(oi+p) = imp_a3(p) + off
        ft%imp_ki(oi+p) = k_imp
        ft%imp_gamma(oi+p) = imp_gamma_tmp(p)
      end do
    end do

    write(*,'(A,I0,A,I0,A,I0,A,I0)') &
      '  Topology/mol: ', nb_mol, ' bonds, ', nang_mol, &
      ' angles, ', ndih_mol, ' dihedrals, ', nimp_mol
    write(*,'(A,I0,A,I0,A,I0,A,I0)') &
      '  Total:        ', ft%Nb, ' bonds, ', ft%Nang, &
      ' angles, ', ft%Ndih, ' dihedrals, ', ft%Nimp
  end subroutine

  ! ═══ LJ近接リスト構築 (half-list, 分子間のみ) ═══
  subroutine build_nlist_lj(pos, h, hi, Na, mol_id, nlc, nll)
    double precision, intent(in) :: pos(:), h(9), hi(9)
    integer, intent(in) :: Na, mol_id(:)
    integer, intent(out) :: nlc(:), nll(:)
    double precision :: rc2, dx, dy, dz
    integer :: i, j
    rc2 = (RCUT + 2.0d0)**2
    nlc(1:Na) = 0
    do i = 1, Na-1
      do j = i+1, Na
        if (mol_id(j) == mol_id(i)) cycle
        dx = pos((j-1)*3+1) - pos((i-1)*3+1)
        dy = pos((j-1)*3+2) - pos((i-1)*3+2)
        dz = pos((j-1)*3+3) - pos((i-1)*3+3)
        call mimg_flat(dx, dy, dz, hi, h)
        if (dx*dx + dy*dy + dz*dz < rc2) then
          if (nlc(i) < MAX_LJ_NEIGH) then
            nlc(i) = nlc(i) + 1
            nll((i-1)*MAX_LJ_NEIGH + nlc(i)) = j
          end if
        end if
      end do
    end do
  end subroutine

  ! ═══ 力計算カーネル ═══
  function compute_forces(Fv, vir9, pos, h, hi, ft, nlc, nll, Na) result(fr)
    double precision, intent(out) :: Fv(:), vir9(9)
    double precision, intent(in) :: pos(:), h(9), hi(9)
    type(FlatTopology), intent(in) :: ft
    integer, intent(in) :: nlc(:), nll(:), Na
    type(ForceResult) :: fr

    double precision :: Eb, Ea, Ed, Ei, Elj
    double precision :: dx, dy, dz, r, dr, fm, fx, fy, fz, r2, ri2, sr2, sr6, sr12
    double precision :: rji0, rji1, rji2, rjk0, rjk1, rjk2
    double precision :: dji, djk, costh, th, dth, sinth, dV
    double precision :: rji_arr(3), rjk_arr(3), fi_c, fk_c
    double precision :: b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z
    double precision :: mx, my, mz, nx, ny, nz, mm, nn, imm, inn, b2len
    double precision :: cosphi, sinphi, phi, dphi, Vn, gamma
    double precision :: b1db2, b2db3, coef2i, coef2k
    double precision :: f1x, f1y, f1z, f4x, f4y, f4z
    double precision :: f2x, f2y, f2z, f3x, f3y, f3z
    double precision :: ki, gm
    integer :: b, a, d, p, i, j, k, l, c, nni, jn, mult
    integer :: i0, i1, i2, i3, c0, c1, c2, c3

    ! 力とビリアルをゼロ初期化
    Fv(1:Na*3) = 0.0d0
    vir9 = 0.0d0
    Eb = 0.0d0; Ea = 0.0d0; Ed = 0.0d0; Ei = 0.0d0; Elj = 0.0d0

    ! ── 1. ボンド伸縮  V_bond = 0.5*kb*(r-r0)^2 ──
    !$OMP PARALLEL DO PRIVATE(b,i,j,dx,dy,dz,r,dr,fm,fx,fy,fz) &
    !$OMP REDUCTION(+:Eb) SCHEDULE(STATIC)
    do b = 1, ft%Nb
      i = ft%b_i0(b); j = ft%b_i1(b)
      dx = pos((j-1)*3+1) - pos((i-1)*3+1)
      dy = pos((j-1)*3+2) - pos((i-1)*3+2)
      dz = pos((j-1)*3+3) - pos((i-1)*3+3)
      call mimg_flat(dx, dy, dz, hi, h)
      r = sqrt(dx*dx + dy*dy + dz*dz)
      if (r < 1d-10) cycle
      dr = r - ft%b_r0(b)
      Eb = Eb + 0.5d0 * ft%b_kb(b) * dr * dr
      fm = -ft%b_kb(b) * dr / r
      fx = fm * dx; fy = fm * dy; fz = fm * dz
      !$OMP ATOMIC
      Fv((i-1)*3+1) = Fv((i-1)*3+1) + fx
      !$OMP ATOMIC
      Fv((i-1)*3+2) = Fv((i-1)*3+2) + fy
      !$OMP ATOMIC
      Fv((i-1)*3+3) = Fv((i-1)*3+3) + fz
      !$OMP ATOMIC
      Fv((j-1)*3+1) = Fv((j-1)*3+1) - fx
      !$OMP ATOMIC
      Fv((j-1)*3+2) = Fv((j-1)*3+2) - fy
      !$OMP ATOMIC
      Fv((j-1)*3+3) = Fv((j-1)*3+3) - fz
      !$OMP ATOMIC
      vir9(1) = vir9(1) + dx * fx
      !$OMP ATOMIC
      vir9(5) = vir9(5) + dy * fy
      !$OMP ATOMIC
      vir9(9) = vir9(9) + dz * fz
    end do
    !$OMP END PARALLEL DO

    ! ── 2. 角度曲げ  V_angle = 0.5*kth*(th-th0)^2 ──
    !$OMP PARALLEL DO PRIVATE(a,i,j,k,rji0,rji1,rji2,rjk0,rjk1,rjk2, &
    !$OMP   dji,djk,costh,th,dth,sinth,dV,rji_arr,rjk_arr,fi_c,fk_c,c) &
    !$OMP REDUCTION(+:Ea) SCHEDULE(STATIC)
    do a = 1, ft%Nang
      i = ft%ang_i0(a); j = ft%ang_i1(a); k = ft%ang_i2(a)
      rji0 = pos((i-1)*3+1) - pos((j-1)*3+1)
      rji1 = pos((i-1)*3+2) - pos((j-1)*3+2)
      rji2 = pos((i-1)*3+3) - pos((j-1)*3+3)
      rjk0 = pos((k-1)*3+1) - pos((j-1)*3+1)
      rjk1 = pos((k-1)*3+2) - pos((j-1)*3+2)
      rjk2 = pos((k-1)*3+3) - pos((j-1)*3+3)
      call mimg_flat(rji0, rji1, rji2, hi, h)
      call mimg_flat(rjk0, rjk1, rjk2, hi, h)
      dji = sqrt(rji0*rji0 + rji1*rji1 + rji2*rji2)
      djk = sqrt(rjk0*rjk0 + rjk1*rjk1 + rjk2*rjk2)
      if (dji < 1d-10 .or. djk < 1d-10) cycle
      costh = (rji0*rjk0 + rji1*rjk1 + rji2*rjk2) / (dji * djk)
      if (costh > 0.999999d0) costh = 0.999999d0
      if (costh < -0.999999d0) costh = -0.999999d0
      th = acos(costh); dth = th - ft%ang_th0(a)
      Ea = Ea + 0.5d0 * ft%ang_kth(a) * dth * dth
      sinth = sqrt(1.0d0 - costh*costh + 1d-30)
      dV = -ft%ang_kth(a) * dth / sinth
      rji_arr(1) = rji0; rji_arr(2) = rji1; rji_arr(3) = rji2
      rjk_arr(1) = rjk0; rjk_arr(2) = rjk1; rjk_arr(3) = rjk2
      do c = 1, 3
        fi_c = dV * (rjk_arr(c)/(dji*djk) - costh*rji_arr(c)/(dji*dji))
        fk_c = dV * (rji_arr(c)/(dji*djk) - costh*rjk_arr(c)/(djk*djk))
        !$OMP ATOMIC
        Fv((i-1)*3+c) = Fv((i-1)*3+c) + fi_c
        !$OMP ATOMIC
        Fv((k-1)*3+c) = Fv((k-1)*3+c) + fk_c
        !$OMP ATOMIC
        Fv((j-1)*3+c) = Fv((j-1)*3+c) - fi_c - fk_c
      end do
    end do
    !$OMP END PARALLEL DO

    ! ── 3. 二面角ねじれ (Bekker解析的力) ──
    !$OMP PARALLEL DO PRIVATE(d,i0,i1,i2,i3,b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z, &
    !$OMP   mx,my,mz,nx,ny,nz,mm,nn,imm,inn,b2len,cosphi,sinphi,phi,mult,gamma,Vn, &
    !$OMP   dphi,b1db2,b2db3,coef2i,coef2k, &
    !$OMP   f1x,f1y,f1z,f4x,f4y,f4z,f2x,f2y,f2z,f3x,f3y,f3z) &
    !$OMP REDUCTION(+:Ed) SCHEDULE(DYNAMIC)
    do d = 1, ft%Ndih
      i0 = ft%dih_i0(d); i1 = ft%dih_i1(d)
      i2 = ft%dih_i2(d); i3 = ft%dih_i3(d)
      b1x = pos((i1-1)*3+1) - pos((i0-1)*3+1)
      b1y = pos((i1-1)*3+2) - pos((i0-1)*3+2)
      b1z = pos((i1-1)*3+3) - pos((i0-1)*3+3)
      b2x = pos((i2-1)*3+1) - pos((i1-1)*3+1)
      b2y = pos((i2-1)*3+2) - pos((i1-1)*3+2)
      b2z = pos((i2-1)*3+3) - pos((i1-1)*3+3)
      b3x = pos((i3-1)*3+1) - pos((i2-1)*3+1)
      b3y = pos((i3-1)*3+2) - pos((i2-1)*3+2)
      b3z = pos((i3-1)*3+3) - pos((i2-1)*3+3)
      call mimg_flat(b1x, b1y, b1z, hi, h)
      call mimg_flat(b2x, b2y, b2z, hi, h)
      call mimg_flat(b3x, b3y, b3z, hi, h)
      mx = b1y*b2z - b1z*b2y; my = b1z*b2x - b1x*b2z; mz = b1x*b2y - b1y*b2x
      nx = b2y*b3z - b2z*b3y; ny = b2z*b3x - b2x*b3z; nz = b2x*b3y - b2y*b3x
      mm = mx*mx + my*my + mz*mz; nn = nx*nx + ny*ny + nz*nz
      if (mm < 1d-20 .or. nn < 1d-20) cycle
      imm = 1.0d0/mm; inn = 1.0d0/nn
      b2len = sqrt(b2x*b2x + b2y*b2y + b2z*b2z)
      cosphi = (mx*nx + my*ny + mz*nz) / sqrt(mm*nn)
      if (cosphi > 1.0d0) cosphi = 1.0d0
      if (cosphi < -1.0d0) cosphi = -1.0d0
      b1db2 = b1x*b2x + b1y*b2y + b1z*b2z
      b2db3 = b2x*b3x + b2y*b3y + b2z*b3z
      sinphi = (mx*b3x + my*b3y + mz*b3z) * b2len / sqrt(mm*nn)
      phi = atan2(sinphi, cosphi)
      mult = ft%dih_mult(d); gamma = ft%dih_gamma(d); Vn = ft%dih_Vn(d)
      Ed = Ed + 0.5d0*Vn*(1.0d0 + cos(dble(mult)*phi - gamma))
      dphi = -0.5d0*Vn*dble(mult)*sin(dble(mult)*phi - gamma)
      f1x = dphi*b2len*imm*mx; f1y = dphi*b2len*imm*my; f1z = dphi*b2len*imm*mz
      f4x = -dphi*b2len*inn*nx; f4y = -dphi*b2len*inn*ny; f4z = -dphi*b2len*inn*nz
      coef2i = b1db2 / (b2len*b2len); coef2k = b2db3 / (b2len*b2len)
      f2x = -f1x + coef2i*f1x - coef2k*f4x
      f2y = -f1y + coef2i*f1y - coef2k*f4y
      f2z = -f1z + coef2i*f1z - coef2k*f4z
      f3x = -f4x - coef2i*f1x + coef2k*f4x
      f3y = -f4y - coef2i*f1y + coef2k*f4y
      f3z = -f4z - coef2i*f1z + coef2k*f4z
      !$OMP ATOMIC
      Fv((i0-1)*3+1) = Fv((i0-1)*3+1) + f1x
      !$OMP ATOMIC
      Fv((i0-1)*3+2) = Fv((i0-1)*3+2) + f1y
      !$OMP ATOMIC
      Fv((i0-1)*3+3) = Fv((i0-1)*3+3) + f1z
      !$OMP ATOMIC
      Fv((i1-1)*3+1) = Fv((i1-1)*3+1) + f2x
      !$OMP ATOMIC
      Fv((i1-1)*3+2) = Fv((i1-1)*3+2) + f2y
      !$OMP ATOMIC
      Fv((i1-1)*3+3) = Fv((i1-1)*3+3) + f2z
      !$OMP ATOMIC
      Fv((i2-1)*3+1) = Fv((i2-1)*3+1) + f3x
      !$OMP ATOMIC
      Fv((i2-1)*3+2) = Fv((i2-1)*3+2) + f3y
      !$OMP ATOMIC
      Fv((i2-1)*3+3) = Fv((i2-1)*3+3) + f3z
      !$OMP ATOMIC
      Fv((i3-1)*3+1) = Fv((i3-1)*3+1) + f4x
      !$OMP ATOMIC
      Fv((i3-1)*3+2) = Fv((i3-1)*3+2) + f4y
      !$OMP ATOMIC
      Fv((i3-1)*3+3) = Fv((i3-1)*3+3) + f4z
    end do
    !$OMP END PARALLEL DO

    ! ── 4. 不適切二面角 (同じ解析式, center=i0) ──
    !$OMP PARALLEL DO PRIVATE(p,c0,c1,c2,c3,b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z, &
    !$OMP   mx,my,mz,nx,ny,nz,mm,nn,imm,inn,b2len,cosphi,sinphi,phi,ki,gm,dphi, &
    !$OMP   b1db2,b2db3,coef2i,coef2k, &
    !$OMP   f1x,f1y,f1z,f4x,f4y,f4z,f2x,f2y,f2z,f3x,f3y,f3z) &
    !$OMP REDUCTION(+:Ei) SCHEDULE(STATIC)
    do p = 1, ft%Nimp
      ! center=c0, neighbors=c1,c2,c3  ->  二面角順序: c1-c0-c2-c3
      c0 = ft%imp_i0(p); c1 = ft%imp_i1(p)
      c2 = ft%imp_i2(p); c3 = ft%imp_i3(p)
      b1x = pos((c0-1)*3+1) - pos((c1-1)*3+1)
      b1y = pos((c0-1)*3+2) - pos((c1-1)*3+2)
      b1z = pos((c0-1)*3+3) - pos((c1-1)*3+3)
      b2x = pos((c2-1)*3+1) - pos((c0-1)*3+1)
      b2y = pos((c2-1)*3+2) - pos((c0-1)*3+2)
      b2z = pos((c2-1)*3+3) - pos((c0-1)*3+3)
      b3x = pos((c3-1)*3+1) - pos((c2-1)*3+1)
      b3y = pos((c3-1)*3+2) - pos((c2-1)*3+2)
      b3z = pos((c3-1)*3+3) - pos((c2-1)*3+3)
      call mimg_flat(b1x, b1y, b1z, hi, h)
      call mimg_flat(b2x, b2y, b2z, hi, h)
      call mimg_flat(b3x, b3y, b3z, hi, h)
      mx = b1y*b2z - b1z*b2y; my = b1z*b2x - b1x*b2z; mz = b1x*b2y - b1y*b2x
      nx = b2y*b3z - b2z*b3y; ny = b2z*b3x - b2x*b3z; nz = b2x*b3y - b2y*b3x
      mm = mx*mx + my*my + mz*mz; nn = nx*nx + ny*ny + nz*nz
      if (mm < 1d-20 .or. nn < 1d-20) cycle
      imm = 1.0d0/mm; inn = 1.0d0/nn
      b2len = sqrt(b2x*b2x + b2y*b2y + b2z*b2z)
      sinphi = (mx*b3x + my*b3y + mz*b3z) * b2len / sqrt(mm*nn)
      cosphi = (mx*nx + my*ny + mz*nz) / sqrt(mm*nn)
      if (cosphi > 1.0d0) cosphi = 1.0d0
      if (cosphi < -1.0d0) cosphi = -1.0d0
      phi = atan2(sinphi, cosphi)
      ki = ft%imp_ki(p); gm = ft%imp_gamma(p)
      Ei = Ei + 0.5d0*ki*(1.0d0 + cos(2.0d0*phi - gm))
      dphi = -0.5d0*ki*2.0d0*sin(2.0d0*phi - gm)
      b1db2 = b1x*b2x + b1y*b2y + b1z*b2z
      b2db3 = b2x*b3x + b2y*b3y + b2z*b3z
      f1x = dphi*b2len*imm*mx; f1y = dphi*b2len*imm*my; f1z = dphi*b2len*imm*mz
      f4x = -dphi*b2len*inn*nx; f4y = -dphi*b2len*inn*ny; f4z = -dphi*b2len*inn*nz
      coef2i = b1db2 / (b2len*b2len); coef2k = b2db3 / (b2len*b2len)
      ! atom1->c1, atom2->c0, atom3->c2, atom4->c3
      !$OMP ATOMIC
      Fv((c1-1)*3+1) = Fv((c1-1)*3+1) + f1x
      !$OMP ATOMIC
      Fv((c1-1)*3+2) = Fv((c1-1)*3+2) + f1y
      !$OMP ATOMIC
      Fv((c1-1)*3+3) = Fv((c1-1)*3+3) + f1z
      f2x = -f1x + coef2i*f1x - coef2k*f4x
      f2y = -f1y + coef2i*f1y - coef2k*f4y
      f2z = -f1z + coef2i*f1z - coef2k*f4z
      !$OMP ATOMIC
      Fv((c0-1)*3+1) = Fv((c0-1)*3+1) + f2x
      !$OMP ATOMIC
      Fv((c0-1)*3+2) = Fv((c0-1)*3+2) + f2y
      !$OMP ATOMIC
      Fv((c0-1)*3+3) = Fv((c0-1)*3+3) + f2z
      f3x = -f4x - coef2i*f1x + coef2k*f4x
      f3y = -f4y - coef2i*f1y + coef2k*f4y
      f3z = -f4z - coef2i*f1z + coef2k*f4z
      !$OMP ATOMIC
      Fv((c2-1)*3+1) = Fv((c2-1)*3+1) + f3x
      !$OMP ATOMIC
      Fv((c2-1)*3+2) = Fv((c2-1)*3+2) + f3y
      !$OMP ATOMIC
      Fv((c2-1)*3+3) = Fv((c2-1)*3+3) + f3z
      !$OMP ATOMIC
      Fv((c3-1)*3+1) = Fv((c3-1)*3+1) + f4x
      !$OMP ATOMIC
      Fv((c3-1)*3+2) = Fv((c3-1)*3+2) + f4y
      !$OMP ATOMIC
      Fv((c3-1)*3+3) = Fv((c3-1)*3+3) + f4z
    end do
    !$OMP END PARALLEL DO

    ! ── 5. LJ分子間相互作用 (half-list) ──
    !$OMP PARALLEL DO PRIVATE(i,nni,jn,j,dx,dy,dz,r2,ri2,sr2,sr6,sr12,fm) &
    !$OMP REDUCTION(+:Elj) SCHEDULE(DYNAMIC,4)
    do i = 1, Na
      nni = nlc(i)
      do jn = 1, nni
        j = nll((i-1)*MAX_LJ_NEIGH + jn)
        dx = pos((j-1)*3+1) - pos((i-1)*3+1)
        dy = pos((j-1)*3+2) - pos((i-1)*3+2)
        dz = pos((j-1)*3+3) - pos((i-1)*3+3)
        call mimg_flat(dx, dy, dz, hi, h)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2 > RCUT2) cycle
        if (r2 < 0.25d0) r2 = 0.25d0
        ri2 = 1.0d0/r2; sr2 = sig2_LJ*ri2; sr6 = sr2*sr2*sr2; sr12 = sr6*sr6
        fm = 24.0d0*eps_LJ*(2.0d0*sr12 - sr6)*ri2
        Elj = Elj + 4.0d0*eps_LJ*(sr12 - sr6) - VSHFT
        Fv((i-1)*3+1) = Fv((i-1)*3+1) - fm*dx
        Fv((i-1)*3+2) = Fv((i-1)*3+2) - fm*dy
        Fv((i-1)*3+3) = Fv((i-1)*3+3) - fm*dz
        !$OMP ATOMIC
        Fv((j-1)*3+1) = Fv((j-1)*3+1) + fm*dx
        !$OMP ATOMIC
        Fv((j-1)*3+2) = Fv((j-1)*3+2) + fm*dy
        !$OMP ATOMIC
        Fv((j-1)*3+3) = Fv((j-1)*3+3) + fm*dz
        !$OMP ATOMIC
        vir9(1) = vir9(1) + dx*fm*dx
        !$OMP ATOMIC
        vir9(5) = vir9(5) + dy*fm*dy
        !$OMP ATOMIC
        vir9(9) = vir9(9) + dz*fm*dz
      end do
    end do
    !$OMP END PARALLEL DO

    fr%Eb = Eb; fr%Ea = Ea; fr%Ed = Ed; fr%Ei = Ei
    fr%Elj = Elj; fr%Ecoul = 0.0d0
    fr%Etot = Eb + Ea + Ed + Ei + Elj
  end function

  ! ═══ 運動エネルギー (全原子) ═══
  double precision function ke_total(vel, Na)
    double precision, intent(in) :: vel(:)
    integer, intent(in) :: Na
    double precision :: s
    integer :: i
    s = 0.0d0
    !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:s)
    do i = 1, Na
      s = s + mC*(vel((i-1)*3+1)**2 + vel((i-1)*3+2)**2 + vel((i-1)*3+3)**2)
    end do
    !$OMP END PARALLEL DO
    ke_total = 0.5d0 * s / CONV
  end function

  ! ═══ 瞬間温度 ═══
  double precision function inst_T(KE, Nf)
    double precision, intent(in) :: KE
    integer, intent(in) :: Nf
    inst_T = 2.0d0 * KE / (dble(Nf) * kB)
  end function

  ! ═══ 瞬間圧力 ═══
  double precision function inst_P(Wm9, KE, V)
    double precision, intent(in) :: Wm9(9), KE, V
    inst_P = (2.0d0*KE + Wm9(1) + Wm9(5) + Wm9(9)) / (3.0d0*V) * eV2GPa
  end function

  ! ═══ NPT状態初期化 (自由度 = 3*Na - 3, 全原子マイナスCM) ═══
  subroutine make_npt(npt, T, Pe, Na)
    type(NPTState), intent(out) :: npt
    double precision, intent(in) :: T, Pe
    integer, intent(in) :: Na
    npt%Nf = 3*Na - 3
    npt%xi = 0.0d0
    npt%Q = max(dble(npt%Nf)*kB*T*1d4, 1d-20)
    npt%Vg = 0.0d0
    npt%W = max(dble(npt%Nf+9)*kB*T*1d6, 1d-20)
    npt%Pe = Pe; npt%Tt = T
  end subroutine

  ! ═══ クランプ ═══
  double precision function clamp_val(x, lo, hi)
    double precision, intent(in) :: x, lo, hi
    clamp_val = max(lo, min(x, hi))
  end function

  ! ═══ NPTステップ ═══
  subroutine step_npt_mmmd(pos, vel, Fv, vir9, h, hi, &
       Na, dt, npt, ft, nlc, nll, fr_out, KE_out)
    double precision, intent(inout) :: pos(:), vel(:), Fv(:)
    double precision, intent(inout) :: vir9(9), h(9), hi(9)
    integer, intent(in) :: Na
    double precision, intent(in) :: dt
    type(NPTState), intent(inout) :: npt
    type(FlatTopology), intent(in) :: ft
    integer, intent(in) :: nlc(:), nll(:)
    type(ForceResult), intent(out) :: fr_out
    double precision, intent(out) :: KE_out

    double precision :: hdt, V, KE, dP, eps_tr, sc_nh, sc_pr, sc_v, mi_inv
    double precision :: px, py, pz, vx, vy, vz, sx, sy, sz, vsx, vsy, vsz
    double precision :: eps2, sc_v2, V2
    integer :: i, a, idx

    hdt = 0.5d0 * dt
    call mat_inv9(h, hi)
    V = abs(mat_det9(h))
    KE = ke_total(vel, Na)

    ! ── Nose-Hooverサーモスタット (半ステップ) ──
    npt%xi = npt%xi + hdt*(2.0d0*KE - dble(npt%Nf)*kB*npt%Tt) / npt%Q
    npt%xi = clamp_val(npt%xi, -0.05d0, 0.05d0)

    ! ── バロスタット (半ステップ) ──
    dP = inst_P(vir9, KE, V) - npt%Pe
    do a = 0, 2
      npt%Vg(a*4+1) = npt%Vg(a*4+1) + hdt*V*dP/(npt%W*eV2GPa)
      npt%Vg(a*4+1) = clamp_val(npt%Vg(a*4+1), -0.005d0, 0.005d0)
    end do

    eps_tr = npt%Vg(1)*hi(1) + npt%Vg(5)*hi(5) + npt%Vg(9)*hi(9)
    sc_nh = exp(-hdt*npt%xi)
    sc_pr = exp(-hdt*eps_tr/3.0d0)
    sc_v = sc_nh * sc_pr
    mi_inv = CONV / mC

    ! ── 速度半ステップ更新 ──
    !$OMP PARALLEL DO PRIVATE(i,idx,a)
    do i = 1, Na
      idx = (i-1)*3
      do a = 1, 3
        vel(idx+a) = vel(idx+a)*sc_v + hdt*Fv(idx+a)*mi_inv
      end do
    end do
    !$OMP END PARALLEL DO

    ! ── 位置更新 (分率座標経由) ──
    !$OMP PARALLEL DO PRIVATE(i,idx,px,py,pz,vx,vy,vz,sx,sy,sz,vsx,vsy,vsz)
    do i = 1, Na
      idx = (i-1)*3
      px = pos(idx+1); py = pos(idx+2); pz = pos(idx+3)
      vx = vel(idx+1); vy = vel(idx+2); vz = vel(idx+3)
      sx = hi(1)*px + hi(2)*py + hi(3)*pz
      sy = hi(4)*px + hi(5)*py + hi(6)*pz
      sz = hi(7)*px + hi(8)*py + hi(9)*pz
      vsx = hi(1)*vx + hi(2)*vy + hi(3)*vz
      vsy = hi(4)*vx + hi(5)*vy + hi(6)*vz
      vsz = hi(7)*vx + hi(8)*vy + hi(9)*vz
      sx = sx + dt*vsx; sy = sy + dt*vsy; sz = sz + dt*vsz
      sx = sx - floor(sx); sy = sy - floor(sy); sz = sz - floor(sz)
      pos(idx+1) = sx; pos(idx+2) = sy; pos(idx+3) = sz
    end do
    !$OMP END PARALLEL DO

    ! ── セル行列更新 ──
    do a = 0, 2
      h(a*3+1) = h(a*3+1) + dt*npt%Vg(a*3+1)
      h(a*3+2) = h(a*3+2) + dt*npt%Vg(a*3+2)
      h(a*3+3) = h(a*3+3) + dt*npt%Vg(a*3+3)
    end do
    call mat_inv9(h, hi)

    ! ── 分率→カーテシアン変換 ──
    !$OMP PARALLEL DO PRIVATE(i,idx,sx,sy,sz)
    do i = 1, Na
      idx = (i-1)*3
      sx = pos(idx+1); sy = pos(idx+2); sz = pos(idx+3)
      pos(idx+1) = h(1)*sx + h(2)*sy + h(3)*sz
      pos(idx+2) = h(4)*sx + h(5)*sy + h(6)*sz
      pos(idx+3) = h(7)*sx + h(8)*sy + h(9)*sz
    end do
    !$OMP END PARALLEL DO

    ! ── 力再計算 ──
    fr_out = compute_forces(Fv, vir9, pos, h, hi, ft, nlc, nll, Na)

    ! ── 速度残り半ステップ ──
    eps2 = npt%Vg(1)*hi(1) + npt%Vg(5)*hi(5) + npt%Vg(9)*hi(9)
    sc_v2 = sc_nh * exp(-hdt*eps2/3.0d0)

    !$OMP PARALLEL DO PRIVATE(i,idx,a)
    do i = 1, Na
      idx = (i-1)*3
      do a = 1, 3
        vel(idx+a) = (vel(idx+a) + hdt*Fv(idx+a)*mi_inv) * sc_v2
      end do
    end do
    !$OMP END PARALLEL DO

    ! ── Nose-Hoover & バロスタット残り半ステップ ──
    KE = ke_total(vel, Na); KE_out = KE
    npt%xi = npt%xi + hdt*(2.0d0*KE - dble(npt%Nf)*kB*npt%Tt) / npt%Q
    npt%xi = clamp_val(npt%xi, -0.05d0, 0.05d0)
    V2 = abs(mat_det9(h))
    dP = inst_P(vir9, KE, V2) - npt%Pe
    do a = 0, 2
      npt%Vg(a*4+1) = npt%Vg(a*4+1) + hdt*V2*dP/(npt%W*eV2GPa)
      npt%Vg(a*4+1) = clamp_val(npt%Vg(a*4+1), -0.005d0, 0.005d0)
    end do
  end subroutine

  ! ═══ OVITO XYZ出力 ═══
  subroutine write_ovito_mmmd(iu, istep, dt, pos, vel, mol_id, h, Na)
    integer, intent(in) :: iu, istep, Na
    double precision, intent(in) :: dt, pos(:), vel(:), h(9)
    integer, intent(in) :: mol_id(:)
    integer :: i
    write(iu, '(I0)') Na
    write(iu, '(A,9(ES18.8),A,F10.4,A,I0,A)') 'Lattice="', &
      h(1), h(4), h(7), h(2), h(5), h(8), h(3), h(6), h(9), &
      '" Properties=species:S:1:pos:R:3:c_mol:I:1 Time=', istep*dt, &
      ' Step=', istep, ' pbc="T T T"'
    do i = 1, Na
      write(iu, '(A,3(ES16.8),I6)') 'C ', &
        pos((i-1)*3+1), pos((i-1)*3+2), pos((i-1)*3+3), mol_id(i)
    end do
  end subroutine

  ! ═══ リスタート保存 ═══
  subroutine write_restart_mmmd(fname, istep, h, npt, pos, vel, Na, Nmol, natom)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: istep, Na, Nmol, natom
    double precision, intent(in) :: h(9), pos(:), vel(:)
    type(NPTState), intent(in) :: npt
    integer :: i, ios
    open(unit=30, file=fname, status='replace', iostat=ios)
    if (ios /= 0) return
    write(30, '(A)') '# RESTART fuller_LJ_npt_mmmd_serial_omp_acc'
    write(30, '(A,I0)') 'STEP ', istep
    write(30, '(A,I0)') 'NMOL ', Nmol
    write(30, '(A,I0)') 'NATOM_MOL ', natom
    write(30, '(A,I0)') 'NATOM ', Na
    write(30, '(A,9(ES20.12))') 'H ', h
    write(30, '(A,5(ES20.12),I8)') 'NPT ', npt%xi, npt%Q, npt%W, npt%Pe, npt%Tt, npt%Nf
    write(30, '(A,9(ES20.12))') 'VG ', npt%Vg
    do i = 1, Na
      write(30, '(A,I6,6(ES20.12))') 'ATOM', i, &
        pos((i-1)*3+1), pos((i-1)*3+2), pos((i-1)*3+3), &
        vel((i-1)*3+1), vel((i-1)*3+2), vel((i-1)*3+3)
    end do
    write(30, '(A)') 'END'
    close(30)
  end subroutine

  ! ═══ ガウス乱数 (Box-Muller) ═══
  double precision function gauss_rand()
    double precision :: u1, u2
    call random_number(u1); call random_number(u2)
    if (u1 < 1d-30) u1 = 1d-30
    gauss_rand = sqrt(-2.0d0*log(u1)) * cos(2.0d0*PI_*u2)
  end function

  ! ═══ コマンドライン引数パーサー ═══
  subroutine get_opt_int(key, defval, result)
    character(len=*), intent(in) :: key
    integer, intent(in) :: defval
    integer, intent(out) :: result
    integer :: i, nargs, kl, ios
    character(len=256) :: arg
    result = defval; nargs = command_argument_count(); kl = len_trim(key)
    do i = 1, nargs
      call get_command_argument(i, arg)
      if (arg(1:kl+3) == '--'//key(1:kl)//'=') then
        read(arg(kl+4:), *, iostat=ios) result
      end if
    end do
  end subroutine

  subroutine get_opt_dbl(key, defval, result)
    character(len=*), intent(in) :: key
    double precision, intent(in) :: defval
    double precision, intent(out) :: result
    integer :: i, nargs, kl, ios
    character(len=256) :: arg
    result = defval; nargs = command_argument_count(); kl = len_trim(key)
    do i = 1, nargs
      call get_command_argument(i, arg)
      if (arg(1:kl+3) == '--'//key(1:kl)//'=') then
        read(arg(kl+4:), *, iostat=ios) result
      end if
    end do
  end subroutine

end module fuller_mmmd_mod

! ═══════════════════════════════════════════════════════════════════════
!  メインプログラム
! ═══════════════════════════════════════════════════════════════════════
program fuller_LJ_npt_mmmd
  use fuller_mmmd_mod
  implicit none

  integer :: nc, Nmol, natom, Na, nsteps, nlup, gstep
  integer :: coldstart, warmup, avg_from, avg_to, total_steps, gavg_from, gavg_to
  integer :: nrec_o, nrec_rst, start_step, mon_int, cur_prn, prn, prn_pre, nav
  double precision :: T, Pe, dt, a0, Dmol_val, Rmol_val
  double precision :: KE, V_val, Tn, Pn, an_val, scale_v, tgt
  double precision :: sT, sP, sa, sEb, sEa, sEd, sElj, sEt
  double precision :: T_cold, T_init, sv
  double precision :: t_start, t_now, elapsed
  double precision :: ff_kb_val, ff_kth_val, ff_v2_val, ff_kimp_val
  double precision :: h(9), hi(9), vir9(9), vcm(3)
  double precision :: mol_coords(MAX_NATOM, 3)
  double precision, allocatable :: pos(:), vel(:), Fv(:), mol_centers(:)
  integer, allocatable :: nlc(:), nll(:)
  integer :: bond_i_mol(4096), bond_j_mol(4096), nb_mol
  type(NPTState) :: npt
  type(FlatTopology) :: ft
  type(ForceResult) :: fr
  integer :: i, a, idx, m
  character(len=6) :: phase_str

  T_cold = 4.0d0

  ! ── コマンドライン引数の取得 ──
  call get_opt_int('cell', 3, nc)
  call get_opt_dbl('temp', 298.0d0, T)
  call get_opt_dbl('pres', 0.0d0, Pe)
  call get_opt_int('step', 10000, nsteps)
  call get_opt_dbl('dt', 0.1d0, dt)       ! MM版デフォルト = 0.1 fs
  call get_opt_int('coldstart', 0, coldstart)
  call get_opt_int('warmup', 0, warmup)
  call get_opt_int('from', 0, avg_from)
  call get_opt_int('to', 0, avg_to)
  call get_opt_int('ovito', 0, nrec_o)
  call get_opt_int('restart', 0, nrec_rst)
  call get_opt_int('mon', 0, mon_int)

  ! 力場パラメータ (kcal/mol入力 → eV変換)
  call get_opt_dbl('ff_kb', 469.0d0, ff_kb_val);   ff_kb_val = ff_kb_val * kcal2eV
  call get_opt_dbl('ff_kth', 63.0d0, ff_kth_val);  ff_kth_val = ff_kth_val * kcal2eV
  call get_opt_dbl('ff_v2', 14.5d0, ff_v2_val);    ff_v2_val = ff_v2_val * kcal2eV
  call get_opt_dbl('ff_kimp', 15.0d0, ff_kimp_val); ff_kimp_val = ff_kimp_val * kcal2eV

  if (avg_to <= 0) avg_to = nsteps
  if (avg_from <= 0) avg_from = max(1, nsteps - nsteps/4)
  total_steps = coldstart + warmup + nsteps
  gavg_from = coldstart + warmup + avg_from
  gavg_to = coldstart + warmup + avg_to
  start_step = 0; nlup = 20

  ! ── C60フォールバック生成 ──
  call generate_c60_mmmd(mol_coords, natom, Rmol_val, Dmol_val, &
                         bond_i_mol, bond_j_mol, nb_mol)
  a0 = Dmol_val * sqrt(2.0d0) * 1.4d0

  ! ── FCC結晶構築 ──
  Nmol = 4 * nc * nc * nc
  Na = Nmol * natom

  ! ── 配列確保 ──
  allocate(mol_centers(Nmol*3))
  allocate(pos(Na*3), vel(Na*3), Fv(Na*3))
  allocate(nlc(Na), nll(Na*MAX_LJ_NEIGH))

  pos = 0; vel = 0; Fv = 0; h = 0; hi = 0; vir9 = 0
  nlc = 0; nll = 0

  ! 分子中心をFCC格子に配置
  Nmol = make_fcc(a0, nc, mol_centers, h)
  Na = Nmol * natom
  call mat_inv9(h, hi)

  ! 各原子の座標を分子中心+テンプレートで初期化
  do m = 1, Nmol
    do a = 1, natom
      idx = ((m-1)*natom + a - 1)*3
      pos(idx+1) = mol_centers((m-1)*3+1) + mol_coords(a, 1)
      pos(idx+2) = mol_centers((m-1)*3+2) + mol_coords(a, 2)
      pos(idx+3) = mol_centers((m-1)*3+3) + mol_coords(a, 3)
    end do
  end do
  deallocate(mol_centers)

  ! ── トポロジー構築 ──
  call build_flat_topology(ft, mol_coords, natom, &
       bond_i_mol, bond_j_mol, nb_mol, Nmol, &
       ff_kb_val, ff_kth_val, ff_v2_val, ff_kimp_val)

  write(*,'(A)') '========================================================================'
  write(*,'(A)') '  Fullerene Crystal NPT-MD — MM Force Field (Serial, Fortran 95)'
  write(*,'(A)') '========================================================================'
  write(*,'(A,I0)') '  Atoms/molecule  : ', natom
  write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '  FCC ', nc, 'x', nc, 'x', nc, &
    '  Nmol=', Nmol, '  Natom=', Na
  write(*,'(A,F8.3,A,F6.1,A,F6.4,A,F5.3,A)') &
    '  a0=', a0, ' A  T=', T, ' K  P=', Pe, ' GPa  dt=', dt, ' fs'
  if (coldstart > 0) write(*,'(A,I0,A,F4.1,A)') &
    '  Coldstart       : ', coldstart, ' steps at ', T_cold, ' K'
  if (warmup > 0) write(*,'(A,I0,A)') &
    '  Warmup          : ', warmup, ' steps'
  write(*,'(A,I0,A,I0,A,I0)') &
    '  Production      : ', nsteps, ' steps  avg=', avg_from, '-', avg_to
  write(*,'(A,I0)') '  Total           : ', total_steps
  write(*,'(A)') '========================================================================'
  write(*,*)

  ! ── 初期速度 ──
  T_init = T; if (coldstart > 0 .or. warmup > 0) T_init = T_cold
  call random_seed()
  sv = sqrt(kB * T_init * CONV / mC)
  do i = 1, Na
    idx = (i-1)*3
    do a = 1, 3
      vel(idx+a) = sv * gauss_rand()
    end do
  end do
  ! 重心速度除去
  vcm = 0.0d0
  do i = 1, Na
    idx = (i-1)*3
    vcm(1) = vcm(1) + vel(idx+1)
    vcm(2) = vcm(2) + vel(idx+2)
    vcm(3) = vcm(3) + vel(idx+3)
  end do
  vcm = vcm / dble(Na)
  do i = 1, Na
    idx = (i-1)*3
    vel(idx+1) = vel(idx+1) - vcm(1)
    vel(idx+2) = vel(idx+2) - vcm(2)
    vel(idx+3) = vel(idx+3) - vcm(3)
  end do

  call make_npt(npt, T, Pe, Na); npt%Tt = T_init
  call mat_inv9(h, hi)
  call build_nlist_lj(pos, h, hi, Na, ft%mol_id, nlc, nll)
  call apply_pbc(pos, h, hi, Na)
  fr = compute_forces(Fv, vir9, pos, h, hi, ft, nlc, nll, Na)

  prn = mon_int; if (prn <= 0) prn = max(1, total_steps/50)
  prn_pre = prn; if (coldstart + warmup > 0) prn_pre = max(1, (coldstart+warmup)/100)
  sT = 0; sP = 0; sa = 0; sEb = 0; sEa = 0; sEd = 0; sElj = 0; sEt = 0; nav = 0
  call cpu_time(t_start)

  ! OVITO出力ファイル
  if (nrec_o > 0) open(unit=40, file='ovito_traj_mmmd_serial.xyz', status='replace')

  write(*,'(A8,A6,A8,A10,A9,A10,A10,A10,A10,A10,A8)') &
    'step', 'phase', 'T[K]', 'P[GPa]', 'a[A]', 'E_bond', 'E_angle', &
    'E_dih', 'E_LJ', 'E_total', 't[s]'

  ! ═══ MDメインループ ═══
  do gstep = start_step + 1, total_steps
    ! ── フェーズ判定 ──
    if (gstep <= coldstart) then
      npt%Tt = T_cold; phase_str = ' COLD'
    else if (gstep <= coldstart + warmup) then
      npt%Tt = T_cold + (T - T_cold) * dble(gstep - coldstart) / dble(warmup)
      phase_str = ' WARM'
    else
      npt%Tt = T; phase_str = ' PROD'
    end if
    if (coldstart > 0 .and. gstep == coldstart + 1) then
      npt%xi = 0; npt%Vg = 0
    end if
    if (gstep <= coldstart) npt%Vg = 0
    cur_prn = prn; if (gstep <= coldstart + warmup) cur_prn = prn_pre

    ! ── 近接リスト更新 ──
    if (mod(gstep, nlup) == 0) then
      call mat_inv9(h, hi)
      call build_nlist_lj(pos, h, hi, Na, ft%mol_id, nlc, nll)
    end if

    ! ── NPTステップ ──
    call step_npt_mmmd(pos, vel, Fv, vir9, h, hi, &
         Na, dt, npt, ft, nlc, nll, fr, KE)

    V_val = abs(mat_det9(h))
    Tn = inst_T(KE, npt%Nf)
    Pn = inst_P(vir9, KE, V_val)

    ! ── COLD/WARM velocity rescaling ──
    if ((gstep <= coldstart .or. gstep <= coldstart + warmup) .and. Tn > 0.1d0) then
      tgt = T_cold; if (gstep > coldstart) tgt = npt%Tt
      scale_v = sqrt(max(tgt, 0.1d0) / Tn)
      !$OMP PARALLEL DO PRIVATE(i,idx)
      do i = 1, Na
        idx = (i-1)*3
        vel(idx+1) = vel(idx+1) * scale_v
        vel(idx+2) = vel(idx+2) * scale_v
        vel(idx+3) = vel(idx+3) * scale_v
      end do
      !$OMP END PARALLEL DO
      KE = ke_total(vel, Na)
      Tn = inst_T(KE, npt%Nf); npt%xi = 0
      if (gstep <= coldstart) npt%Vg = 0
    end if

    an_val = h(1) / dble(nc)

    ! ── 平均値蓄積 ──
    if (gstep >= gavg_from .and. gstep <= gavg_to) then
      sT = sT + Tn; sP = sP + Pn; sa = sa + an_val
      sEb = sEb + fr%Eb / dble(Nmol)
      sEa = sEa + fr%Ea / dble(Nmol)
      sEd = sEd + (fr%Ed + fr%Ei) / dble(Nmol)
      sElj = sElj + fr%Elj / dble(Nmol)
      sEt = sEt + fr%Etot / dble(Nmol)
      nav = nav + 1
    end if

    ! ── OVITO出力 ──
    if (nrec_o > 0 .and. mod(gstep, nrec_o) == 0) then
      call write_ovito_mmmd(40, gstep, dt, pos, vel, ft%mol_id, h, Na)
      flush(40)
    end if

    ! ── リスタート保存 ──
    if (nrec_rst > 0 .and. (mod(gstep, nrec_rst) == 0 .or. gstep == total_steps)) then
      call write_restart_mmmd('restart_mmmd.rst', gstep, h, npt, pos, vel, Na, Nmol, natom)
    end if

    ! ── モニタリング出力 ──
    if (mod(gstep, cur_prn) == 0 .or. gstep == total_steps) then
      call cpu_time(t_now); elapsed = t_now - t_start
      write(*,'(I8,A6,F8.1,F10.3,F9.3,F10.4,F10.4,F10.4,F10.4,F10.4,F8.0)') &
        gstep, phase_str, Tn, Pn, an_val, &
        fr%Eb/dble(Nmol), fr%Ea/dble(Nmol), &
        (fr%Ed+fr%Ei)/dble(Nmol), fr%Elj/dble(Nmol), &
        fr%Etot/dble(Nmol), elapsed
    end if
  end do

  if (nrec_o > 0) close(40)

  ! ── 平均値出力 ──
  if (nav > 0) then
    write(*,*)
    write(*,'(A)') '========================================================================'
    write(*,'(A,I0,A,F7.2,A,F8.4,A,F8.4,A)') &
      '  Averages (', nav, '): T=', sT/dble(nav), ' K  P=', sP/dble(nav), &
      ' GPa  a=', sa/dble(nav), ' A'
    write(*,'(A,F9.4,A,F9.4,A,F9.4,A,F9.4,A,F9.4)') &
      '  bond=', sEb/dble(nav), '  ang=', sEa/dble(nav), &
      '  dih=', sEd/dble(nav), '  LJ=', sElj/dble(nav), &
      '  tot=', sEt/dble(nav)
    write(*,'(A)') '========================================================================'
  end if
  call cpu_time(t_now)
  write(*,'(A,F8.1,A)') '  Done (', t_now - t_start, ' sec)'

  ! ── 後始末 ──
  deallocate(pos, vel, Fv, nlc, nll)
  if (allocated(ft%b_i0)) deallocate(ft%b_i0)
  if (allocated(ft%b_i1)) deallocate(ft%b_i1)
  if (allocated(ft%b_kb)) deallocate(ft%b_kb)
  if (allocated(ft%b_r0)) deallocate(ft%b_r0)
  if (allocated(ft%ang_i0)) deallocate(ft%ang_i0)
  if (allocated(ft%ang_i1)) deallocate(ft%ang_i1)
  if (allocated(ft%ang_i2)) deallocate(ft%ang_i2)
  if (allocated(ft%ang_kth)) deallocate(ft%ang_kth)
  if (allocated(ft%ang_th0)) deallocate(ft%ang_th0)
  if (allocated(ft%dih_i0)) deallocate(ft%dih_i0)
  if (allocated(ft%dih_i1)) deallocate(ft%dih_i1)
  if (allocated(ft%dih_i2)) deallocate(ft%dih_i2)
  if (allocated(ft%dih_i3)) deallocate(ft%dih_i3)
  if (allocated(ft%dih_Vn)) deallocate(ft%dih_Vn)
  if (allocated(ft%dih_mult)) deallocate(ft%dih_mult)
  if (allocated(ft%dih_gamma)) deallocate(ft%dih_gamma)
  if (allocated(ft%imp_i0)) deallocate(ft%imp_i0)
  if (allocated(ft%imp_i1)) deallocate(ft%imp_i1)
  if (allocated(ft%imp_i2)) deallocate(ft%imp_i2)
  if (allocated(ft%imp_i3)) deallocate(ft%imp_i3)
  if (allocated(ft%imp_ki)) deallocate(ft%imp_ki)
  if (allocated(ft%imp_gamma)) deallocate(ft%imp_gamma)
  if (allocated(ft%mol_id)) deallocate(ft%mol_id)

end program fuller_LJ_npt_mmmd
