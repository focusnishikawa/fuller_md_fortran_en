! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2025, Takeshi Nishikawa
!===========================================================================
!  fuller_LJ_npt_md_serial_omp_acc.f90
!  Fullerene Crystal NPT-MD (Rigid-body LJ, Serial / OpenMP, Fortran 95)
!
!  コンパイル:
!    Serial:  gfortran -O3 -o fuller_LJ_npt_md_serial &
!               fuller_LJ_npt_md_serial_omp_acc.f90
!    OpenMP:  gfortran -O3 -fopenmp -o fuller_LJ_npt_md_omp &
!               fuller_LJ_npt_md_serial_omp_acc.f90
!
!  実行時オプション (全て --key=value 形式):
!    --help                  このヘルプを表示
!    --fullerene=<名前>      フラーレン種 (デフォルト: C60)
!    --crystal=<fcc|hcp|bcc> 結晶構造 (デフォルト: fcc)
!    --cell=<nc>             単位胞の繰り返し数 (デフォルト: 3)
!    --temp=<T_K>            目標温度 [K] (デフォルト: 298.0)
!    --pres=<P_GPa>          目標圧力 [GPa] (デフォルト: 0.0)
!    --step=<N>              本計算ステップ数 (デフォルト: 10000)
!    --dt=<fs>               時間刻み [fs] (デフォルト: 1.0)
!    --init_scale=<s>        格子定数スケール因子 (デフォルト: 1.0)
!    --seed=<n>              乱数シード (デフォルト: 42)
!    --coldstart=<N>         低温(4K)ステップ数 (デフォルト: 0)
!    --warmup=<N>            昇温ステップ数 4K→T (デフォルト: 0)
!    --from=<step>           平均開始ステップ (デフォルト: 本計算の3/4地点)
!    --to=<step>             平均終了ステップ (デフォルト: nsteps)
!    --mon=<N>               モニタリング出力間隔 (デフォルト: 自動)
!    --warmup_mon=<mode>     昇温中の出力頻度 norm|freq|some (デフォルト: norm)
!    --ovito=<N>             OVITO XYZ出力間隔 (0=無効, デフォルト: 0)
!    --ofile=<filename>      OVITO出力ファイル名 (デフォルト: 自動生成)
!    --restart=<N>           リスタート保存間隔 (0=無効, デフォルト: 0)
!    --resfile=<path>        リスタートファイルから再開
!    --libdir=<path>         フラーレンライブラリ (デフォルト: FullereneLib)
!
!  停止制御:
!    実行中にカレントディレクトリに以下を作成すると動作を制御できる:
!    - abort.md: 即座に停止 (リスタート有効時は保存してから終了)
!    - stop.md:  次のリスタートチェックポイントで停止
!
!  単位系: A, amu, eV, fs, K, GPa
!===========================================================================

module fuller_LJ_full_mod
  implicit none

  ! ═══ 物理定数・単位変換 ═══
  double precision, parameter :: CONV       = 9.64853321d-3   ! eV*fs^2/(amu*A^2)
  double precision, parameter :: kB         = 8.617333262d-5  ! ボルツマン定数 [eV/K]
  double precision, parameter :: eV2GPa     = 160.21766208d0  ! eV/A^3 → GPa
  double precision, parameter :: eV2kcalmol = 23.06054783d0   ! eV → kcal/mol
  double precision, parameter :: PI_        = 3.14159265358979323846d0

  ! ═══ LJポテンシャルパラメータ ═══
  double precision, parameter :: sigma_LJ = 3.431d0
  double precision, parameter :: eps_LJ   = 2.635d-3
  double precision, parameter :: RCUT     = 3.0d0 * sigma_LJ
  double precision, parameter :: RCUT2    = RCUT * RCUT
  double precision, parameter :: sig2_LJ  = sigma_LJ * sigma_LJ
  double precision, parameter :: mC       = 12.011d0

  ! VSHFT: V(RCUT) を差し引いてカットオフ不連続を除去
  double precision, parameter :: sr_v  = 1.0d0 / 3.0d0
  double precision, parameter :: sr2_v = sr_v * sr_v
  double precision, parameter :: sr6_v = sr2_v * sr2_v * sr2_v
  double precision, parameter :: VSHFT = 4.0d0*eps_LJ*(sr6_v*sr6_v - sr6_v)

  integer, parameter :: MAX_NATOM = 84   ! C84 support
  integer, parameter :: MAX_NEIGH = 80

  ! ═══ NPT状態 ═══
  type :: NPTState
    double precision :: xi, Q
    double precision :: Vg(9)
    double precision :: W, Pe, Tt
    integer :: Nf
  end type

contains

  ! ═══════════ H行列ユーティリティ ═══════════
  ! h(9) = {H00,H01,H02, H10,H11,H12, H20,H21,H22}
  ! アクセス: h(3*i + j + 1)  (i,j = 0,1,2)

  double precision function mat_det9(h)
    double precision, intent(in) :: h(9)
    mat_det9 = h(1)*(h(5)*h(9) - h(6)*h(8)) &
             - h(2)*(h(4)*h(9) - h(6)*h(7)) &
             + h(3)*(h(4)*h(8) - h(5)*h(7))
  end function

  double precision function mat_tr9(h)
    double precision, intent(in) :: h(9)
    mat_tr9 = h(1) + h(5) + h(9)
  end function

  subroutine mat_inv9(h, hi)
    double precision, intent(in)  :: h(9)
    double precision, intent(out) :: hi(9)
    double precision :: d, id
    d = mat_det9(h)
    id = 1.0d0 / d
    hi(1) = id * (h(5)*h(9) - h(6)*h(8))
    hi(2) = id * (h(3)*h(8) - h(2)*h(9))
    hi(3) = id * (h(2)*h(6) - h(3)*h(5))
    hi(4) = id * (h(6)*h(7) - h(4)*h(9))
    hi(5) = id * (h(1)*h(9) - h(3)*h(7))
    hi(6) = id * (h(3)*h(4) - h(1)*h(6))
    hi(7) = id * (h(4)*h(8) - h(5)*h(7))
    hi(8) = id * (h(2)*h(7) - h(1)*h(8))
    hi(9) = id * (h(1)*h(5) - h(2)*h(4))
  end subroutine

  ! ═══ 最小像規約 (周期境界条件) ═══
  subroutine mimg_flat(dx, dy, dz, hi, h)
    double precision, intent(inout) :: dx, dy, dz
    double precision, intent(in)    :: hi(9), h(9)
    double precision :: s0, s1, s2
    s0 = hi(1)*dx + hi(2)*dy + hi(3)*dz
    s1 = hi(4)*dx + hi(5)*dy + hi(6)*dz
    s2 = hi(7)*dx + hi(8)*dy + hi(9)*dz
    s0 = s0 - anint(s0)
    s1 = s1 - anint(s1)
    s2 = s2 - anint(s2)
    dx = h(1)*s0 + h(2)*s1 + h(3)*s2
    dy = h(4)*s0 + h(5)*s1 + h(6)*s2
    dz = h(7)*s0 + h(8)*s1 + h(9)*s2
  end subroutine

  ! ═══ 四元数操作 ═══
  subroutine q2R_flat(q, R)
    double precision, intent(in)  :: q(4)
    double precision, intent(out) :: R(9)
    double precision :: w, x, y, z
    w = q(1); x = q(2); y = q(3); z = q(4)
    R(1) = 1.0d0 - 2.0d0*(y*y + z*z)
    R(2) = 2.0d0*(x*y - w*z)
    R(3) = 2.0d0*(x*z + w*y)
    R(4) = 2.0d0*(x*y + w*z)
    R(5) = 1.0d0 - 2.0d0*(x*x + z*z)
    R(6) = 2.0d0*(y*z - w*x)
    R(7) = 2.0d0*(x*z - w*y)
    R(8) = 2.0d0*(y*z + w*x)
    R(9) = 1.0d0 - 2.0d0*(x*x + y*y)
  end subroutine

  subroutine qmul_flat(a, b, o)
    double precision, intent(in)  :: a(4), b(4)
    double precision, intent(out) :: o(4)
    o(1) = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4)
    o(2) = a(1)*b(2) + a(2)*b(1) + a(3)*b(4) - a(4)*b(3)
    o(3) = a(1)*b(3) - a(2)*b(4) + a(3)*b(1) + a(4)*b(2)
    o(4) = a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + a(4)*b(1)
  end subroutine

  subroutine qnorm_flat(q)
    double precision, intent(inout) :: q(4)
    double precision :: n, inv
    n = sqrt(q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4))
    inv = 1.0d0 / n
    q(1) = q(1)*inv; q(2) = q(2)*inv
    q(3) = q(3)*inv; q(4) = q(4)*inv
  end subroutine

  subroutine omega2dq_flat(wx, wy, wz, dt, dq)
    double precision, intent(in)  :: wx, wy, wz, dt
    double precision, intent(out) :: dq(4)
    double precision :: wm, th, s
    wm = sqrt(wx*wx + wy*wy + wz*wz)
    th = wm * dt * 0.5d0
    if (th < 1.0d-14) then
      dq(1) = 1.0d0
      dq(2) = 0.5d0*dt*wx
      dq(3) = 0.5d0*dt*wy
      dq(4) = 0.5d0*dt*wz
    else
      s = sin(th) / wm
      dq(1) = cos(th)
      dq(2) = s*wx; dq(3) = s*wy; dq(4) = s*wz
    end if
  end subroutine

  ! ═══ cc1ファイル読み込み ═══
  subroutine load_cc1(path, coords, natom, Rmol, Dmol, I0, Mmol)
    character(len=*), intent(in) :: path
    double precision, intent(out) :: coords(MAX_NATOM, 3)
    double precision, intent(out) :: Rmol, Dmol, I0, Mmol
    integer, intent(out) :: natom
    double precision :: cm(3), r2, r, dx, dy, dz, d, Ixx, Iyy, Izz
    double precision :: x, y, z
    integer :: i, j, ios, idx_val, fl
    character(len=512) :: line
    character(len=8) :: elem

    open(unit=20, file=path, status='old', iostat=ios)
    if (ios /= 0) then
      write(*, '(A,A)') 'Error: cannot open ', trim(path)
      stop 1
    end if
    read(20, *) natom
    if (natom > MAX_NATOM) then
      write(*, '(A,I0,A,I0)') 'Error: natom=', natom, ' > MAX_NATOM=', MAX_NATOM
      stop 1
    end if
    coords = 0.0d0
    do i = 1, natom
      read(20, '(A)', iostat=ios) line
      do while (len_trim(line) == 0 .and. ios == 0)
        read(20, '(A)', iostat=ios) line
      end do
      read(line, *, iostat=ios) elem, idx_val, x, y, z, fl
      coords(i, 1) = x; coords(i, 2) = y; coords(i, 3) = z
    end do
    close(20)

    ! 重心除去
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

    ! Rmol, Dmol の計算
    Rmol = 0.0d0; Dmol = 0.0d0
    do i = 1, natom
      r2 = coords(i,1)**2 + coords(i,2)**2 + coords(i,3)**2
      r = sqrt(r2)
      if (r > Rmol) Rmol = r
      do j = i+1, natom
        dx = coords(i,1) - coords(j,1)
        dy = coords(i,2) - coords(j,2)
        dz = coords(i,3) - coords(j,3)
        d = sqrt(dx*dx + dy*dy + dz*dz)
        if (d > Dmol) Dmol = d
      end do
    end do

    ! 慣性モーメント
    Ixx = 0.0d0; Iyy = 0.0d0; Izz = 0.0d0
    do i = 1, natom
      r2 = coords(i,1)**2 + coords(i,2)**2 + coords(i,3)**2
      Ixx = Ixx + mC * (r2 - coords(i,1)**2)
      Iyy = Iyy + mC * (r2 - coords(i,2)**2)
      Izz = Izz + mC * (r2 - coords(i,3)**2)
    end do
    I0 = (Ixx + Iyy + Izz) / 3.0d0
    Mmol = dble(natom) * mC
  end subroutine

  ! ═══ フラーレン名前解決 ═══
  subroutine resolve_fullerene(spec, libdir, fpath, label)
    character(len=*), intent(in) :: spec, libdir
    character(len=512), intent(out) :: fpath
    character(len=128), intent(out) :: label
    character(len=256) :: sl
    integer :: i, slen

    ! 小文字化
    sl = spec
    slen = len_trim(sl)
    do i = 1, slen
      if (sl(i:i) >= 'A' .and. sl(i:i) <= 'Z') then
        sl(i:i) = achar(iachar(sl(i:i)) + 32)
      end if
    end do

    if (trim(sl) == 'buckyball' .or. trim(sl) == 'c60' .or. &
        trim(sl) == 'c60:ih') then
      fpath = trim(libdir) // '/C60-76/C60-Ih.cc1'
      label = 'C60(Ih)'
      return
    end if
    if (trim(sl) == 'c70' .or. trim(sl) == 'c70:d5h') then
      fpath = trim(libdir) // '/C60-76/C70-D5h.cc1'
      label = 'C70(D5h)'
      return
    end if
    if (trim(sl) == 'c72' .or. trim(sl) == 'c72:d6d') then
      fpath = trim(libdir) // '/C60-76/C72-D6d.cc1'
      label = 'C72(D6d)'
      return
    end if
    if (trim(sl) == 'c74' .or. trim(sl) == 'c74:d3h') then
      fpath = trim(libdir) // '/C60-76/C74-D3h.cc1'
      label = 'C74(D3h)'
      return
    end if
    ! C76:sym 形式
    if (slen >= 5 .and. sl(1:4) == 'c76:') then
      fpath = trim(libdir) // '/C60-76/C76-' // trim(spec(5:)) // '.cc1'
      label = 'C76(' // trim(spec(5:)) // ')'
      return
    end if
    ! C84:num:sym or C84:num 形式 (簡易実装)
    if (slen >= 5 .and. sl(1:4) == 'c84:') then
      ! C84:num の場合、番号のみで .cc1 を探す
      call resolve_c84(spec(5:), libdir, fpath, label)
      return
    end if

    write(*, '(A,A)') 'Unknown fullerene: ', trim(spec)
    stop 1
  end subroutine

  ! ═══ C84解決ヘルパー ═══
  subroutine resolve_c84(rest, libdir, fpath, label)
    character(len=*), intent(in) :: rest, libdir
    character(len=512), intent(out) :: fpath
    character(len=128), intent(out) :: label
    integer :: colon_pos, num, i
    character(len=64) :: sym, numstr
    logical :: fexist

    ! コロン位置を探す
    colon_pos = 0
    do i = 1, len_trim(rest)
      if (rest(i:i) == ':') then
        colon_pos = i
        exit
      end if
    end do

    if (colon_pos > 0) then
      ! C84:num:sym 形式
      read(rest(1:colon_pos-1), *, err=900) num
      sym = rest(colon_pos+1:)
      write(fpath, '(A,A,I2.2,A,A,A)') trim(libdir), '/C84/C84-No.', num, '-', trim(sym), '.cc1'
      write(label, '(A,I0)') 'C84 No.', num
      return
    else
      ! C84:num 形式 — 番号のみ
      read(rest, *, err=900) num
      write(numstr, '(A,I2.2,A)') 'C84-No.', num, '-'
      ! 既知の対称性を試す
      write(fpath, '(A,A,I2.2,A)') trim(libdir), '/C84/C84-No.', num, '-D2.cc1'
      inquire(file=trim(fpath), exist=fexist)
      if (fexist) then
        write(label, '(A,I0)') 'C84 No.', num
        return
      end if
      write(fpath, '(A,A,I2.2,A)') trim(libdir), '/C84/C84-No.', num, '-D2d.cc1'
      inquire(file=trim(fpath), exist=fexist)
      if (fexist) then
        write(label, '(A,I0)') 'C84 No.', num
        return
      end if
      write(fpath, '(A,A,I2.2,A)') trim(libdir), '/C84/C84-No.', num, '-Td.cc1'
      inquire(file=trim(fpath), exist=fexist)
      if (fexist) then
        write(label, '(A,I0)') 'C84 No.', num
        return
      end if
      write(fpath, '(A,A,I2.2,A)') trim(libdir), '/C84/C84-No.', num, '-D6h.cc1'
      inquire(file=trim(fpath), exist=fexist)
      if (fexist) then
        write(label, '(A,I0)') 'C84 No.', num
        return
      end if
    end if
900 continue
    write(*, '(A,A)') 'Unknown C84 isomer: ', trim(rest)
    stop 1
  end subroutine

  ! ═══ C60座標の生成 (cc1ファイルなしの場合のフォールバック) ═══
  subroutine generate_c60(coords, natom, I0, Mmol, Rmol, Dmol)
    double precision, intent(out) :: coords(MAX_NATOM, 3)
    double precision, intent(out) :: I0, Mmol, Rmol, Dmol
    integer, intent(out) :: natom
    double precision :: phi, tmp(60, 3), cm(3), r2, r, Isum, dx, dy, dz, d
    integer :: n, p, s1, s2, s3, cyc(3, 3), signs(2), i, j

    natom = 60
    Mmol = 60.0d0 * mC
    phi = (1.0d0 + sqrt(5.0d0)) / 2.0d0
    signs(1) = -1; signs(2) = 1
    cyc(1,1) = 1; cyc(1,2) = 2; cyc(1,3) = 3
    cyc(2,1) = 2; cyc(2,2) = 3; cyc(2,3) = 1
    cyc(3,1) = 3; cyc(3,2) = 1; cyc(3,3) = 2
    tmp = 0.0d0; n = 0

    ! 群1: (0, +-1, +-3*phi) の巡回置換 -> 12頂点
    do p = 1, 3
      do s2 = 1, 2
        do s3 = 1, 2
          n = n + 1
          tmp(n, cyc(p,1)) = 0.0d0
          tmp(n, cyc(p,2)) = dble(signs(s2))
          tmp(n, cyc(p,3)) = dble(signs(s3)) * 3.0d0 * phi
        end do
      end do
    end do

    ! 群2: (+-2, +-(1+2*phi), +-phi) の巡回置換 -> 24頂点
    do p = 1, 3
      do s1 = 1, 2
        do s2 = 1, 2
          do s3 = 1, 2
            n = n + 1
            tmp(n, cyc(p,1)) = dble(signs(s1)) * 2.0d0
            tmp(n, cyc(p,2)) = dble(signs(s2)) * (1.0d0 + 2.0d0*phi)
            tmp(n, cyc(p,3)) = dble(signs(s3)) * phi
          end do
        end do
      end do
    end do

    ! 群3: (+-1, +-(2+phi), +-2*phi) の巡回置換 -> 24頂点
    do p = 1, 3
      do s1 = 1, 2
        do s2 = 1, 2
          do s3 = 1, 2
            n = n + 1
            tmp(n, cyc(p,1)) = dble(signs(s1))
            tmp(n, cyc(p,2)) = dble(signs(s2)) * (2.0d0 + phi)
            tmp(n, cyc(p,3)) = dble(signs(s3)) * 2.0d0 * phi
          end do
        end do
      end do
    end do

    ! 重心除去 + スケーリング
    cm = 0.0d0
    do n = 1, 60
      cm(1) = cm(1) + tmp(n, 1)
      cm(2) = cm(2) + tmp(n, 2)
      cm(3) = cm(3) + tmp(n, 3)
    end do
    cm = cm / 60.0d0
    do n = 1, 60
      tmp(n, 1) = (tmp(n, 1) - cm(1)) * 0.72d0
      tmp(n, 2) = (tmp(n, 2) - cm(2)) * 0.72d0
      tmp(n, 3) = (tmp(n, 3) - cm(3)) * 0.72d0
    end do

    Rmol = 0.0d0; Isum = 0.0d0; Dmol = 0.0d0
    do i = 1, 60
      r2 = tmp(i,1)**2 + tmp(i,2)**2 + tmp(i,3)**2
      r = sqrt(r2)
      if (r > Rmol) Rmol = r
      Isum = Isum + mC * r2
      do j = i+1, 60
        dx = tmp(i,1) - tmp(j,1)
        dy = tmp(i,2) - tmp(j,2)
        dz = tmp(i,3) - tmp(j,3)
        d = sqrt(dx*dx + dy*dy + dz*dz)
        if (d > Dmol) Dmol = d
      end do
    end do
    I0 = Isum * 2.0d0 / 3.0d0
    coords = 0.0d0
    coords(1:60, :) = tmp
  end subroutine

  ! ═══ 結晶構造: FCC ═══
  function make_fcc(a, nc, pos, h) result(Nmol)
    double precision, intent(in)  :: a
    integer, intent(in)           :: nc
    double precision, intent(out) :: pos(:), h(9)
    integer :: Nmol
    double precision :: bas(4, 3)
    integer :: ix, iy, iz, b, idx

    bas(1,:) = (/ 0.0d0,   0.0d0,   0.0d0   /)
    bas(2,:) = (/ 0.5d0*a, 0.5d0*a, 0.0d0   /)
    bas(3,:) = (/ 0.5d0*a, 0.0d0,   0.5d0*a /)
    bas(4,:) = (/ 0.0d0,   0.5d0*a, 0.5d0*a /)

    Nmol = 0
    do ix = 0, nc-1
      do iy = 0, nc-1
        do iz = 0, nc-1
          do b = 1, 4
            Nmol = Nmol + 1
            idx = (Nmol - 1) * 3
            pos(idx+1) = a*dble(ix) + bas(b, 1)
            pos(idx+2) = a*dble(iy) + bas(b, 2)
            pos(idx+3) = a*dble(iz) + bas(b, 3)
          end do
        end do
      end do
    end do
    h = 0.0d0
    h(1) = dble(nc) * a; h(5) = dble(nc) * a; h(9) = dble(nc) * a
  end function

  ! ═══ 結晶構造: HCP ═══
  function make_hcp(a, nc, pos, h) result(Nmol)
    double precision, intent(in)  :: a
    integer, intent(in)           :: nc
    double precision, intent(out) :: pos(:), h(9)
    integer :: Nmol
    double precision :: c_lat, a1(3), a2(3), a3(3), bas(2, 3)
    double precision :: fx, fy, fz
    integer :: ix, iy, iz, b, idx

    c_lat = a * sqrt(8.0d0 / 3.0d0)
    a1(1) = a;        a1(2) = 0.0d0;               a1(3) = 0.0d0
    a2(1) = a/2.0d0;  a2(2) = a*sqrt(3.0d0)/2.0d0; a2(3) = 0.0d0
    a3(1) = 0.0d0;    a3(2) = 0.0d0;                a3(3) = c_lat

    bas(1, :) = (/ 0.0d0,       0.0d0,       0.0d0 /)
    bas(2, :) = (/ 1.0d0/3.0d0, 2.0d0/3.0d0, 0.5d0 /)

    Nmol = 0
    do ix = 0, nc-1
      do iy = 0, nc-1
        do iz = 0, nc-1
          do b = 1, 2
            Nmol = Nmol + 1
            idx = (Nmol - 1) * 3
            fx = dble(ix) + bas(b, 1)
            fy = dble(iy) + bas(b, 2)
            fz = dble(iz) + bas(b, 3)
            pos(idx+1) = fx*a1(1) + fy*a2(1) + fz*a3(1)
            pos(idx+2) = fx*a1(2) + fy*a2(2) + fz*a3(2)
            pos(idx+3) = fx*a1(3) + fy*a2(3) + fz*a3(3)
          end do
        end do
      end do
    end do

    h = 0.0d0
    ! H行列: h(1)=nc*a1x, h(4)=nc*a1y, h(2)=nc*a2x, h(5)=nc*a2y, h(9)=nc*a3z
    h(1) = dble(nc) * a1(1)   ! H(0,0)
    h(4) = dble(nc) * a1(2)   ! H(1,0)
    h(2) = dble(nc) * a2(1)   ! H(0,1)
    h(5) = dble(nc) * a2(2)   ! H(1,1)
    h(9) = dble(nc) * a3(3)   ! H(2,2)
  end function

  ! ═══ 結晶構造: BCC ═══
  function make_bcc(a, nc, pos, h) result(Nmol)
    double precision, intent(in)  :: a
    integer, intent(in)           :: nc
    double precision, intent(out) :: pos(:), h(9)
    integer :: Nmol
    double precision :: bas(2, 3)
    integer :: ix, iy, iz, b, idx

    bas(1,:) = (/ 0.0d0,   0.0d0,   0.0d0   /)
    bas(2,:) = (/ 0.5d0*a, 0.5d0*a, 0.5d0*a /)

    Nmol = 0
    do ix = 0, nc-1
      do iy = 0, nc-1
        do iz = 0, nc-1
          do b = 1, 2
            Nmol = Nmol + 1
            idx = (Nmol - 1) * 3
            pos(idx+1) = a*dble(ix) + bas(b, 1)
            pos(idx+2) = a*dble(iy) + bas(b, 2)
            pos(idx+3) = a*dble(iz) + bas(b, 3)
          end do
        end do
      end do
    end do
    h = 0.0d0
    h(1) = dble(nc) * a; h(5) = dble(nc) * a; h(9) = dble(nc) * a
  end function

  ! ═══ デフォルト格子定数 ═══
  double precision function default_a0(dmax, st, s)
    double precision, intent(in) :: dmax, s
    character(len=*), intent(in) :: st
    double precision :: m_factor
    m_factor = 1.4d0
    if (trim(st) == 'FCC') then
      default_a0 = dmax * sqrt(2.0d0) * m_factor * s
    else if (trim(st) == 'HCP') then
      default_a0 = dmax * m_factor * s
    else
      default_a0 = dmax * 2.0d0 / sqrt(3.0d0) * m_factor * s
    end if
  end function

  ! ═══ 近傍リスト構築 (対称フルリスト) ═══
  subroutine nlist_build_sym(pos, h, hi, N, rmcut, nl_count, nl_list)
    double precision, intent(in) :: pos(:), h(9), hi(9), rmcut
    integer, intent(in)          :: N
    integer, intent(out)         :: nl_count(:), nl_list(:)
    double precision :: rc2, dx, dy, dz, r2
    integer :: i, j, ci, cj

    rc2 = (rmcut + 3.0d0)**2
    nl_count(1:N) = 0

    do i = 1, N
      do j = i+1, N
        dx = pos((j-1)*3+1) - pos((i-1)*3+1)
        dy = pos((j-1)*3+2) - pos((i-1)*3+2)
        dz = pos((j-1)*3+3) - pos((i-1)*3+3)
        call mimg_flat(dx, dy, dz, hi, h)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2 < rc2) then
          ci = nl_count(i) + 1
          cj = nl_count(j) + 1
          if (ci <= MAX_NEIGH) then
            nl_list((i-1)*MAX_NEIGH + ci) = j
            nl_count(i) = ci
          end if
          if (cj <= MAX_NEIGH) then
            nl_list((j-1)*MAX_NEIGH + cj) = i
            nl_count(j) = cj
          end if
        end if
      end do
    end do
  end subroutine

  ! ═══ PBC適用 ═══
  subroutine apply_pbc(pos, h, hi, N)
    double precision, intent(inout) :: pos(:)
    double precision, intent(in)    :: h(9), hi(9)
    integer, intent(in)             :: N
    double precision :: px, py, pz, s0, s1, s2
    integer :: i, idx

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,px,py,pz,s0,s1,s2)
    do i = 1, N
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

  ! ═══ 力・トルク・ビリアル計算 (メインカーネル) ═══
  double precision function calc_forces(Fv, Tv, Wm9, pos, qv, body, &
       h, hi, nl_count, nl_list, N, natom, rmcut2, lab)
    double precision, intent(out)   :: Fv(:), Tv(:), Wm9(9)
    double precision, intent(in)    :: pos(:), qv(:), body(:,:)
    double precision, intent(in)    :: h(9), hi(9)
    integer, intent(in)             :: nl_count(:), nl_list(:)
    integer, intent(in)             :: N, natom
    double precision, intent(in)    :: rmcut2
    double precision, intent(inout) :: lab(:)

    double precision :: R(9), bx, by, bz
    double precision :: fi0, fi1, fi2, ti0, ti1, ti2, my_Ep
    double precision :: w00,w01,w02, w10,w11,w12, w20,w21,w22
    double precision :: dmx, dmy, dmz, rax, ray, raz
    double precision :: ddx, ddy, ddz, r2, ri2, sr2, sr6, sr12, fm
    double precision :: fx, fy, fz, Ep
    integer :: i, j, k, ai, bj, nni, ia, jb, idx

    ! ステップ1: lab座標計算
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(R,ai,bx,by,bz,idx)
    do i = 1, N
      call q2R_flat(qv((i-1)*4+1:(i-1)*4+4), R)
      do ai = 1, natom
        bx = body(ai, 1); by = body(ai, 2); bz = body(ai, 3)
        idx = (i-1)*natom*3 + (ai-1)*3
        lab(idx+1) = R(1)*bx + R(2)*by + R(3)*bz
        lab(idx+2) = R(4)*bx + R(5)*by + R(6)*bz
        lab(idx+3) = R(7)*bx + R(8)*by + R(9)*bz
      end do
    end do
    !$OMP END PARALLEL DO

    ! ゼロ初期化
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    do i = 1, N*3
      Fv(i) = 0.0d0; Tv(i) = 0.0d0
    end do
    !$OMP END PARALLEL DO
    Wm9 = 0.0d0; Ep = 0.0d0

    ! ステップ2: LJ力計算メインカーネル
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) REDUCTION(+:Ep) &
    !$OMP   PRIVATE(fi0,fi1,fi2,ti0,ti1,ti2,my_Ep, &
    !$OMP           w00,w01,w02,w10,w11,w12,w20,w21,w22, &
    !$OMP           nni,k,j,dmx,dmy,dmz,ai,ia,rax,ray,raz, &
    !$OMP           bj,jb,ddx,ddy,ddz, &
    !$OMP           r2,ri2,sr2,sr6,sr12,fm,fx,fy,fz)
    do i = 1, N
      fi0=0.0d0; fi1=0.0d0; fi2=0.0d0
      ti0=0.0d0; ti1=0.0d0; ti2=0.0d0
      my_Ep = 0.0d0
      w00=0.0d0;w01=0.0d0;w02=0.0d0
      w10=0.0d0;w11=0.0d0;w12=0.0d0
      w20=0.0d0;w21=0.0d0;w22=0.0d0

      nni = nl_count(i)
      do k = 1, nni
        j = nl_list((i-1)*MAX_NEIGH + k)
        dmx = pos((j-1)*3+1) - pos((i-1)*3+1)
        dmy = pos((j-1)*3+2) - pos((i-1)*3+2)
        dmz = pos((j-1)*3+3) - pos((i-1)*3+3)
        call mimg_flat(dmx, dmy, dmz, hi, h)
        if (dmx*dmx + dmy*dmy + dmz*dmz > rmcut2) cycle

        do ai = 1, natom
          ia = (i-1)*natom*3 + (ai-1)*3
          rax = lab(ia+1); ray = lab(ia+2); raz = lab(ia+3)

          do bj = 1, natom
            jb = (j-1)*natom*3 + (bj-1)*3
            ddx = dmx + lab(jb+1) - rax
            ddy = dmy + lab(jb+2) - ray
            ddz = dmz + lab(jb+3) - raz
            r2 = ddx*ddx + ddy*ddy + ddz*ddz

            if (r2 < RCUT2) then
              if (r2 < 0.25d0) r2 = 0.25d0
              ri2 = 1.0d0 / r2
              sr2 = sig2_LJ * ri2
              sr6 = sr2*sr2*sr2
              sr12 = sr6*sr6
              fm = 24.0d0*eps_LJ*(2.0d0*sr12 - sr6)*ri2
              fx = fm*ddx; fy = fm*ddy; fz = fm*ddz

              fi0 = fi0 - fx; fi1 = fi1 - fy; fi2 = fi2 - fz
              ti0 = ti0 - (ray*fz - raz*fy)
              ti1 = ti1 - (raz*fx - rax*fz)
              ti2 = ti2 - (rax*fy - ray*fx)
              my_Ep = my_Ep + 0.5d0*(4.0d0*eps_LJ*(sr12-sr6) - VSHFT)
              w00=w00+0.5d0*ddx*fx; w01=w01+0.5d0*ddx*fy; w02=w02+0.5d0*ddx*fz
              w10=w10+0.5d0*ddy*fx; w11=w11+0.5d0*ddy*fy; w12=w12+0.5d0*ddy*fz
              w20=w20+0.5d0*ddz*fx; w21=w21+0.5d0*ddz*fy; w22=w22+0.5d0*ddz*fz
            end if
          end do
        end do
      end do

      Fv((i-1)*3+1) = fi0; Fv((i-1)*3+2) = fi1; Fv((i-1)*3+3) = fi2
      Tv((i-1)*3+1) = ti0; Tv((i-1)*3+2) = ti1; Tv((i-1)*3+3) = ti2

      ! Wm9 virial update: use ATOMIC for thread safety
      !$OMP ATOMIC
      Wm9(1) = Wm9(1) + w00
      !$OMP ATOMIC
      Wm9(2) = Wm9(2) + w01
      !$OMP ATOMIC
      Wm9(3) = Wm9(3) + w02
      !$OMP ATOMIC
      Wm9(4) = Wm9(4) + w10
      !$OMP ATOMIC
      Wm9(5) = Wm9(5) + w11
      !$OMP ATOMIC
      Wm9(6) = Wm9(6) + w12
      !$OMP ATOMIC
      Wm9(7) = Wm9(7) + w20
      !$OMP ATOMIC
      Wm9(8) = Wm9(8) + w21
      !$OMP ATOMIC
      Wm9(9) = Wm9(9) + w22

      Ep = Ep + my_Ep
    end do
    !$OMP END PARALLEL DO

    calc_forces = Ep
  end function

  ! ═══ 運動エネルギー ═══
  double precision function ke_trans(vel, N, Mmol)
    double precision, intent(in) :: vel(:), Mmol
    integer, intent(in) :: N
    double precision :: s
    integer :: i, idx
    s = 0.0d0
    !$OMP PARALLEL DO REDUCTION(+:s) PRIVATE(idx)
    do i = 1, N
      idx = (i-1)*3
      s = s + vel(idx+1)**2 + vel(idx+2)**2 + vel(idx+3)**2
    end do
    !$OMP END PARALLEL DO
    ke_trans = 0.5d0 * Mmol * s / CONV
  end function

  double precision function ke_rot(omg, N, I0)
    double precision, intent(in) :: omg(:), I0
    integer, intent(in) :: N
    double precision :: s
    integer :: i, idx
    s = 0.0d0
    !$OMP PARALLEL DO REDUCTION(+:s) PRIVATE(idx)
    do i = 1, N
      idx = (i-1)*3
      s = s + omg(idx+1)**2 + omg(idx+2)**2 + omg(idx+3)**2
    end do
    !$OMP END PARALLEL DO
    ke_rot = 0.5d0 * I0 * s / CONV
  end function

  double precision function inst_T(KE, Nf)
    double precision, intent(in) :: KE
    integer, intent(in) :: Nf
    inst_T = 2.0d0 * KE / (dble(Nf) * kB)
  end function

  double precision function inst_P(Wm9, KEt, V)
    double precision, intent(in) :: Wm9(9), KEt, V
    inst_P = (2.0d0*KEt + Wm9(1) + Wm9(5) + Wm9(9)) / (3.0d0*V) * eV2GPa
  end function

  ! ═══ NPT状態変数の初期化 ═══
  subroutine make_npt(npt, T, Pe, N)
    type(NPTState), intent(out) :: npt
    double precision, intent(in) :: T, Pe
    integer, intent(in) :: N
    npt%Nf = 6*N - 3
    npt%xi = 0.0d0
    npt%Q  = max(dble(npt%Nf)*kB*T*100.0d0*100.0d0, 1.0d-20)
    npt%Vg = 0.0d0
    npt%W  = max(dble(npt%Nf+9)*kB*T*1000.0d0*1000.0d0, 1.0d-20)
    npt%Pe = Pe
    npt%Tt = T
  end subroutine

  ! ═══ クランプ関数 ═══
  double precision function clamp_val(x, lo, hi)
    double precision, intent(in) :: x, lo, hi
    clamp_val = max(lo, min(x, hi))
  end function

  ! ═══ NPT Velocity-Verlet 1ステップ ═══
  subroutine step_npt(pos, vel, qv, omg, Fv, Tv, Wm9, h, hi, &
       body, I0, Mmol, N, natom, rmcut2, dt, npt, &
       nl_count, nl_list, lab, Ep_out, KE_out)
    double precision, intent(inout) :: pos(:), vel(:), qv(:), omg(:)
    double precision, intent(inout) :: Fv(:), Tv(:), Wm9(9)
    double precision, intent(inout) :: h(9), hi(9)
    double precision, intent(in)    :: body(:,:), I0, Mmol, rmcut2, dt
    integer, intent(in)             :: N, natom
    type(NPTState), intent(inout)   :: npt
    integer, intent(in)             :: nl_count(:), nl_list(:)
    double precision, intent(inout) :: lab(:)
    double precision, intent(out)   :: Ep_out, KE_out

    double precision :: hdt, V, kt, kr, KE, dP, eps_tr, sc_nh, sc_pr, sc_v
    double precision :: cF, cT, px, py, pz, vx, vy, vz
    double precision :: sx, sy, sz, vsx, vsy, vsz
    double precision :: dq(4), tmp_q(4), eps_tr2, sc_v2, V2
    integer :: i, a, idx

    hdt = 0.5d0 * dt
    call mat_inv9(h, hi)
    V = abs(mat_det9(h))
    kt = ke_trans(vel, N, Mmol)
    kr = ke_rot(omg, N, I0)
    KE = kt + kr

    ! (A) サーモスタット前半
    npt%xi = npt%xi + hdt*(2.0d0*KE - dble(npt%Nf)*kB*npt%Tt) / npt%Q
    npt%xi = clamp_val(npt%xi, -0.1d0, 0.1d0)

    ! (B) バロスタット前半
    dP = inst_P(Wm9, kt, V) - npt%Pe
    do a = 0, 2
      npt%Vg(a*4+1) = npt%Vg(a*4+1) + hdt*V*dP/(npt%W*eV2GPa)
      npt%Vg(a*4+1) = clamp_val(npt%Vg(a*4+1), -0.01d0, 0.01d0)
    end do

    eps_tr = npt%Vg(1)*hi(1) + npt%Vg(5)*hi(5) + npt%Vg(9)*hi(9)
    sc_nh = exp(-hdt*npt%xi)
    sc_pr = exp(-hdt*eps_tr/3.0d0)
    sc_v = sc_nh * sc_pr
    cF = CONV / Mmol
    cT = CONV / I0

    ! (C) 速度の前半更新
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,a)
    do i = 1, N
      idx = (i-1)*3
      do a = 1, 3
        vel(idx+a) = vel(idx+a)*sc_v + hdt*Fv(idx+a)*cF
        omg(idx+a) = omg(idx+a)*sc_nh + hdt*Tv(idx+a)*cT
      end do
    end do
    !$OMP END PARALLEL DO

    ! (D) 座標更新 (分率座標)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,px,py,pz,vx,vy,vz, &
    !$OMP   sx,sy,sz,vsx,vsy,vsz)
    do i = 1, N
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

    ! (E) セルH行列の更新
    do a = 0, 2
      h(a*3+1) = h(a*3+1) + dt*npt%Vg(a*3+1)
      h(a*3+2) = h(a*3+2) + dt*npt%Vg(a*3+2)
      h(a*3+3) = h(a*3+3) + dt*npt%Vg(a*3+3)
    end do

    ! (F) 分率->実座標
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,sx,sy,sz)
    do i = 1, N
      idx = (i-1)*3
      sx = pos(idx+1); sy = pos(idx+2); sz = pos(idx+3)
      pos(idx+1) = h(1)*sx + h(2)*sy + h(3)*sz
      pos(idx+2) = h(4)*sx + h(5)*sy + h(6)*sz
      pos(idx+3) = h(7)*sx + h(8)*sy + h(9)*sz
    end do
    !$OMP END PARALLEL DO

    ! (G) 四元数更新
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,dq,tmp_q)
    do i = 1, N
      idx = (i-1)*3
      call omega2dq_flat(omg(idx+1), omg(idx+2), omg(idx+3), dt, dq)
      call qmul_flat(qv((i-1)*4+1:(i-1)*4+4), dq, tmp_q)
      qv((i-1)*4+1) = tmp_q(1); qv((i-1)*4+2) = tmp_q(2)
      qv((i-1)*4+3) = tmp_q(3); qv((i-1)*4+4) = tmp_q(4)
      call qnorm_flat(qv((i-1)*4+1:(i-1)*4+4))
    end do
    !$OMP END PARALLEL DO

    ! (H) 力の再計算
    call mat_inv9(h, hi)
    Ep_out = calc_forces(Fv, Tv, Wm9, pos, qv, body, h, hi, &
                         nl_count, nl_list, N, natom, rmcut2, lab)

    ! (I) 速度の後半更新
    eps_tr2 = npt%Vg(1)*hi(1) + npt%Vg(5)*hi(5) + npt%Vg(9)*hi(9)
    sc_v2 = sc_nh * exp(-hdt*eps_tr2/3.0d0)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,a)
    do i = 1, N
      idx = (i-1)*3
      do a = 1, 3
        vel(idx+a) = (vel(idx+a) + hdt*Fv(idx+a)*cF) * sc_v2
        omg(idx+a) = (omg(idx+a) + hdt*Tv(idx+a)*cT) * sc_nh
      end do
    end do
    !$OMP END PARALLEL DO

    ! (J)(K) サーモ/バロ後半
    kt = ke_trans(vel, N, Mmol)
    kr = ke_rot(omg, N, I0)
    KE = kt + kr
    KE_out = KE
    npt%xi = npt%xi + hdt*(2.0d0*KE - dble(npt%Nf)*kB*npt%Tt) / npt%Q
    npt%xi = clamp_val(npt%xi, -0.1d0, 0.1d0)
    V2 = abs(mat_det9(h))
    dP = inst_P(Wm9, kt, V2) - npt%Pe
    do a = 0, 2
      npt%Vg(a*4+1) = npt%Vg(a*4+1) + hdt*V2*dP/(npt%W*eV2GPa)
      npt%Vg(a*4+1) = clamp_val(npt%Vg(a*4+1), -0.01d0, 0.01d0)
    end do
  end subroutine

  ! ═══ OVITO XYZ出力 ═══
  subroutine write_ovito(iu, istep, dt, pos, vel, qv, body, h, N, natom)
    integer, intent(in) :: iu, istep, N, natom
    double precision, intent(in) :: dt, pos(:), vel(:), qv(:), body(:,:), h(9)
    double precision :: R(9), bx, by, bz, rx, ry, rz
    integer :: i, a

    write(iu, '(I0)') N*natom
    write(iu, '(A,9(ES18.8),A,F10.4,A,I0,A)') &
      'Lattice="', &
      h(1),h(4),h(7), h(2),h(5),h(8), h(3),h(6),h(9), &
      '" Properties=species:S:1:pos:R:3:c_mol:I:1:vx:R:1:vy:R:1:vz:R:1 Time=', &
      istep*dt, ' Step=', istep, ' pbc="T T T"'
    do i = 1, N
      call q2R_flat(qv((i-1)*4+1:(i-1)*4+4), R)
      do a = 1, natom
        bx = body(a,1); by = body(a,2); bz = body(a,3)
        rx = pos((i-1)*3+1) + R(1)*bx + R(2)*by + R(3)*bz
        ry = pos((i-1)*3+2) + R(4)*bx + R(5)*by + R(6)*bz
        rz = pos((i-1)*3+3) + R(7)*bx + R(8)*by + R(9)*bz
        write(iu, '(A,3(ES16.8),I6,3(ES16.8))') &
          'C ', rx, ry, rz, i, vel((i-1)*3+1), vel((i-1)*3+2), vel((i-1)*3+3)
      end do
    end do
  end subroutine

  ! ═══ リスタート保存 ═══
  subroutine write_restart_lj(fname, istep, st, nc, T, Pe, nsteps, dt, &
       seed, fspec, init_scale, h, npt, pos, qv, vel, omg, N, natom)
    character(len=*), intent(in) :: fname, st, fspec
    integer, intent(in) :: istep, nc, nsteps, seed, N, natom
    double precision, intent(in) :: T, Pe, dt, init_scale
    double precision, intent(in) :: h(9), pos(:), qv(:), vel(:), omg(:)
    type(NPTState), intent(in) :: npt
    integer :: i, ios

    open(unit=30, file=fname, status='replace', iostat=ios)
    if (ios /= 0) return
    write(30, '(A)') '# RESTART fuller_LJ_npt_md_serial_omp_acc'
    write(30, '(A,I0)') 'STEP ', istep
    write(30, '(A,I0)') 'NSTEPS ', nsteps
    write(30, '(A,ES20.12)') 'DT ', dt
    write(30, '(A,ES20.12)') 'TEMP ', T
    write(30, '(A,ES20.12)') 'PRES ', Pe
    write(30, '(A,A)') 'CRYSTAL ', trim(st)
    write(30, '(A,I0)') 'NC ', nc
    write(30, '(A,A)') 'FULLERENE ', trim(fspec)
    write(30, '(A,ES20.12)') 'INIT_SCALE ', init_scale
    write(30, '(A,I0)') 'SEED ', seed
    write(30, '(A,I0)') 'NMOL ', N
    write(30, '(A,I0)') 'NATOM_MOL ', natom
    write(30, '(A,9(ES20.12))') 'H ', h
    write(30, '(A,5(ES20.12),I8)') 'NPT ', npt%xi, npt%Q, npt%W, npt%Pe, npt%Tt, npt%Nf
    write(30, '(A,9(ES20.12))') 'VG ', npt%Vg
    do i = 1, N
      write(30, '(A,I6,13(ES20.12))') 'MOL ', i, &
        pos((i-1)*3+1), pos((i-1)*3+2), pos((i-1)*3+3), &
        qv((i-1)*4+1), qv((i-1)*4+2), qv((i-1)*4+3), qv((i-1)*4+4), &
        vel((i-1)*3+1), vel((i-1)*3+2), vel((i-1)*3+3), &
        omg((i-1)*3+1), omg((i-1)*3+2), omg((i-1)*3+3)
    end do
    write(30, '(A)') 'END'
    close(30)
  end subroutine

  ! ═══ リスタート読み込み ═══
  subroutine read_restart_lj(fname, rd_istep, rd_N, rd_h, rd_npt, &
       rd_pos, rd_qv, rd_vel, rd_omg, rd_ok)
    character(len=*), intent(in) :: fname
    integer, intent(out) :: rd_istep, rd_N
    double precision, intent(out) :: rd_h(9)
    type(NPTState), intent(out) :: rd_npt
    double precision, allocatable, intent(out) :: rd_pos(:), rd_qv(:), rd_vel(:), rd_omg(:)
    logical, intent(out) :: rd_ok

    character(len=1024) :: line, tag
    integer :: ios, idx_val, mol_count, i
    double precision :: v(13)

    rd_ok = .false.
    rd_istep = 0; rd_N = 0; rd_h = 0.0d0
    rd_npt%xi = 0.0d0; rd_npt%Q = 0.0d0; rd_npt%W = 0.0d0
    rd_npt%Pe = 0.0d0; rd_npt%Tt = 0.0d0; rd_npt%Nf = 0; rd_npt%Vg = 0.0d0

    open(unit=25, file=fname, status='old', iostat=ios)
    if (ios /= 0) return

    ! 最初のパスでNMOLを取得
    mol_count = 0
    do
      read(25, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) tag
      if (trim(tag) == 'NMOL') then
        read(line(5:), *, iostat=ios) rd_N
      else if (trim(tag) == 'MOL') then
        mol_count = mol_count + 1
      else if (trim(tag) == 'END') then
        exit
      end if
    end do

    if (rd_N <= 0) rd_N = mol_count
    if (rd_N <= 0) then
      close(25)
      return
    end if

    ! 配列確保
    allocate(rd_pos(rd_N*3), rd_qv(rd_N*4), rd_vel(rd_N*3), rd_omg(rd_N*3))
    rd_pos = 0.0d0; rd_qv = 0.0d0; rd_vel = 0.0d0; rd_omg = 0.0d0

    ! 2回目のパスでデータを読む
    rewind(25)
    mol_count = 0
    do
      read(25, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) tag

      if (trim(tag) == 'STEP') then
        read(line(5:), *, iostat=ios) rd_istep
      else if (trim(tag) == 'H') then
        read(line(3:), *, iostat=ios) rd_h
      else if (trim(tag) == 'NPT') then
        read(line(4:), *, iostat=ios) rd_npt%xi, rd_npt%Q, rd_npt%W, &
             rd_npt%Pe, rd_npt%Tt, rd_npt%Nf
      else if (trim(tag) == 'VG') then
        read(line(3:), *, iostat=ios) rd_npt%Vg
      else if (trim(tag) == 'MOL') then
        mol_count = mol_count + 1
        if (mol_count <= rd_N) then
          read(line(4:), *, iostat=ios) idx_val, v
          do i = 1, 3
            rd_pos((mol_count-1)*3+i) = v(i)
          end do
          do i = 1, 4
            rd_qv((mol_count-1)*4+i) = v(i+3)
          end do
          do i = 1, 3
            rd_vel((mol_count-1)*3+i) = v(i+7)
          end do
          do i = 1, 3
            rd_omg((mol_count-1)*3+i) = v(i+10)
          end do
        end if
      else if (trim(tag) == 'END') then
        exit
      end if
    end do
    close(25)

    rd_N = mol_count
    rd_ok = (rd_N > 0)
    if (rd_ok) then
      write(*, '(A,A,A,I0,A,I0,A)') '  Restart loaded: ', trim(fname), &
        ' (step ', rd_istep, ', ', rd_N, ' mols)'
    end if
  end subroutine

  ! ═══ リスタートファイル名生成 ═══
  subroutine restart_filename(oname, istep, nsteps, rfn)
    character(len=*), intent(in) :: oname
    integer, intent(in) :: istep, nsteps
    character(len=512), intent(out) :: rfn
    character(len=512) :: base
    integer :: dt_pos, p_pos, dg, x_val
    character(len=32) :: fmt_str, step_str

    ! 拡張子を除去
    base = oname
    dt_pos = index(base, '.', back=.true.)
    if (dt_pos > 0) base = base(1:dt_pos-1)

    ! "ovito_traj" を "restart" に置換
    p_pos = index(base, 'ovito_traj')
    if (p_pos > 0) then
      base = base(1:p_pos-1) // 'restart' // base(p_pos+10:)
    end if

    if (istep == nsteps) then
      rfn = trim(base) // '.rst'
    else
      ! 桁数計算
      dg = 1
      x_val = nsteps
      do while (x_val >= 10)
        x_val = x_val / 10
        dg = dg + 1
      end do
      write(fmt_str, '(A,I0,A,I0,A)') '(A,A,I', dg, '.', dg, ',A)'
      write(rfn, fmt_str) trim(base), '_', istep, '.rst'
    end if
  end subroutine

  ! ═══ 出力サフィックスの生成 ═══
  subroutine build_suffix(fspec, crystal, nc, T, Pe, nsteps, dt, &
       init_scale, seed, sfx)
    character(len=*), intent(in) :: fspec, crystal
    integer, intent(in) :: nc, nsteps, seed
    double precision, intent(in) :: T, Pe, dt, init_scale
    character(len=512), intent(out) :: sfx
    character(len=64) :: tmp_str

    sfx = ''
    ! fullerene
    if (trim(fspec) /= 'C60') then
      sfx = trim(sfx) // '_fullerene_' // trim(fspec)
    end if
    ! crystal
    if (trim(crystal) /= 'fcc') then
      sfx = trim(sfx) // '_crystal_' // trim(crystal)
    end if
    ! cell
    if (nc /= 3) then
      write(tmp_str, '(I0)') nc
      sfx = trim(sfx) // '_cell_' // trim(tmp_str)
    end if
    ! temp
    if (abs(T - 298.0d0) > 1.0d-10) then
      write(tmp_str, '(F12.1)') T
      sfx = trim(sfx) // '_temp_' // trim(adjustl(tmp_str))
    end if
    ! pres
    if (abs(Pe) > 1.0d-10) then
      write(tmp_str, '(F12.1)') Pe
      sfx = trim(sfx) // '_pres_' // trim(adjustl(tmp_str))
    end if
    ! step
    if (nsteps /= 10000) then
      write(tmp_str, '(I0)') nsteps
      sfx = trim(sfx) // '_step_' // trim(tmp_str)
    end if
    ! dt
    if (abs(dt - 1.0d0) > 1.0d-10) then
      write(tmp_str, '(F12.1)') dt
      sfx = trim(sfx) // '_dt_' // trim(adjustl(tmp_str))
    end if
    ! init_scale
    if (abs(init_scale - 1.0d0) > 1.0d-10) then
      write(tmp_str, '(F12.4)') init_scale
      sfx = trim(sfx) // '_init_scale_' // trim(adjustl(tmp_str))
    end if
    ! seed
    if (seed /= 42) then
      write(tmp_str, '(I0)') seed
      sfx = trim(sfx) // '_seed_' // trim(tmp_str)
    end if
  end subroutine

  ! ═══ ユニークファイル名生成 ═══
  subroutine unique_file(base, ext, result_path)
    character(len=*), intent(in) :: base, ext
    character(len=512), intent(out) :: result_path
    logical :: fexist
    integer :: n
    character(len=16) :: n_str

    result_path = trim(base) // trim(ext)
    inquire(file=trim(result_path), exist=fexist)
    if (.not. fexist) return

    do n = 1, 9999
      write(n_str, '(I0)') n
      result_path = trim(base) // '_' // trim(n_str) // trim(ext)
      inquire(file=trim(result_path), exist=fexist)
      if (.not. fexist) return
    end do
    result_path = trim(base) // '_9999' // trim(ext)
  end subroutine

  ! ═══ ガウス乱数 (Box-Muller法) ═══
  double precision function gauss_rand()
    double precision :: u1, u2
    call random_number(u1)
    call random_number(u2)
    if (u1 < 1.0d-30) u1 = 1.0d-30
    gauss_rand = sqrt(-2.0d0*log(u1)) * cos(2.0d0*PI_*u2)
  end function

  ! ═══ コマンドライン引数パーサー ═══
  subroutine get_opt_int(key, defval, result_val)
    character(len=*), intent(in) :: key
    integer, intent(in) :: defval
    integer, intent(out) :: result_val
    integer :: i, nargs, kl, ios
    character(len=512) :: arg
    result_val = defval
    nargs = command_argument_count()
    kl = len_trim(key)
    do i = 1, nargs
      call get_command_argument(i, arg)
      if (len_trim(arg) >= kl+3) then
        if (arg(1:kl+3) == '--' // key(1:kl) // '=') then
          read(arg(kl+4:), *, iostat=ios) result_val
        end if
      end if
    end do
  end subroutine

  subroutine get_opt_dbl(key, defval, result_val)
    character(len=*), intent(in) :: key
    double precision, intent(in) :: defval
    double precision, intent(out) :: result_val
    integer :: i, nargs, kl, ios
    character(len=512) :: arg
    result_val = defval
    nargs = command_argument_count()
    kl = len_trim(key)
    do i = 1, nargs
      call get_command_argument(i, arg)
      if (len_trim(arg) >= kl+3) then
        if (arg(1:kl+3) == '--' // key(1:kl) // '=') then
          read(arg(kl+4:), *, iostat=ios) result_val
        end if
      end if
    end do
  end subroutine

  subroutine get_opt_str(key, defval, result_val)
    character(len=*), intent(in) :: key, defval
    character(len=*), intent(out) :: result_val
    integer :: i, nargs, kl
    character(len=512) :: arg
    result_val = defval
    nargs = command_argument_count()
    kl = len_trim(key)
    do i = 1, nargs
      call get_command_argument(i, arg)
      if (len_trim(arg) >= kl+3) then
        if (arg(1:kl+3) == '--' // key(1:kl) // '=') then
          result_val = arg(kl+4:)
        end if
      end if
    end do
  end subroutine

  ! ═══ --help フラグ検出 ═══
  logical function has_help_flag()
    integer :: i, nargs
    character(len=512) :: arg
    has_help_flag = .false.
    nargs = command_argument_count()
    do i = 1, nargs
      call get_command_argument(i, arg)
      if (trim(arg) == '--help') then
        has_help_flag = .true.
        return
      end if
    end do
  end function

  ! ═══ ディレクトリ存在確認 (停止制御用) ═══
  ! INQUIRE(FILE=) はgfortranではディレクトリにも対応
  logical function dir_exists(path)
    character(len=*), intent(in) :: path
    inquire(file=trim(path), exist=dir_exists)
  end function

  ! ═══ 文字列の小文字化 ═══
  subroutine to_upper(str, out)
    character(len=*), intent(in) :: str
    character(len=*), intent(out) :: out
    integer :: i
    out = str
    do i = 1, len_trim(out)
      if (out(i:i) >= 'a' .and. out(i:i) <= 'z') then
        out(i:i) = achar(iachar(out(i:i)) - 32)
      end if
    end do
  end subroutine

end module fuller_LJ_full_mod


! ═══════════════════════════════════════════════════════════
!  メインプログラム
! ═══════════════════════════════════════════════════════════
program fuller_LJ_npt_md
  use fuller_LJ_full_mod
!$ use omp_lib
  implicit none

  ! ── 変数宣言 ──
  integer :: nc, N, natom, nsteps, nlup, Nmax, i, a, idx, nav, gstep
  integer :: coldstart, warmup, avg_from, avg_to, total_steps
  integer :: gavg_from, gavg_to
  integer :: nrec_o, nrec_rst, start_step, mon_int, cur_prn, prn, prn_pre
  integer :: seed, omp_threads, nn, div
  double precision :: T, Pe, dt, a0, I0, Mmol, Rmol, Dmol
  double precision :: init_scale
  double precision :: sv, sw, n_norm, Ep, KE, kt_val, V_val, Tn, Pn, Ec, an_val
  double precision :: sum_T, sum_P, sum_a, sum_Ep, T_cold, T_init, scale_v, tgt
  double precision :: t_start, t_now, elapsed
  double precision :: RMCUT_v, RMCUT2_v
  double precision :: h(9), hi(9), Wm9(9), vcm(3)
  double precision :: mol_coords(MAX_NATOM, 3)
  double precision, allocatable :: pos(:), vel(:), omg(:), qv(:)
  double precision, allocatable :: Fv(:), Tv(:), lab(:)
  double precision, allocatable :: body(:,:)
  integer, allocatable :: nl_count(:), nl_list(:)
  type(NPTState) :: npt
  logical :: stop_requested, use_cc1
  character(len=512) :: fspec, crystal, st, libdir, resfile
  character(len=512) :: warmup_mon_mode, ofile_opt
  character(len=512) :: fpath, ovito_file, rfn, sfx
  character(len=128) :: label
  character(len=8) :: phase

  ! リスタート用
  integer :: rd_istep, rd_N
  double precision :: rd_h(9)
  type(NPTState) :: rd_npt
  double precision, allocatable :: rd_pos(:), rd_qv(:), rd_vel(:), rd_omg(:)
  logical :: rd_ok

  T_cold = 4.0d0

  ! ── ヘルプ表示 ──
  if (has_help_flag()) then
    write(*, '(A)') 'fuller_LJ_npt_md_serial_omp_acc -- LJ rigid-body fullerene NPT-MD'
    write(*, '(A)') ''
    write(*, '(A)') 'Options:'
    write(*, '(A)') '  --help                  Show this help'
    write(*, '(A)') '  --fullerene=<name>      Fullerene species (default: C60)'
    write(*, '(A)') '  --crystal=<fcc|hcp|bcc> Crystal structure (default: fcc)'
    write(*, '(A)') '  --cell=<nc>             Unit cell repeats (default: 3)'
    write(*, '(A)') '  --temp=<K>              Target temperature [K] (default: 298.0)'
    write(*, '(A)') '  --pres=<GPa>            Target pressure [GPa] (default: 0.0)'
    write(*, '(A)') '  --step=<N>              Production steps (default: 10000)'
    write(*, '(A)') '  --dt=<fs>               Time step [fs] (default: 1.0)'
    write(*, '(A)') '  --init_scale=<s>        Lattice scale factor (default: 1.0)'
    write(*, '(A)') '  --seed=<n>              Random seed (default: 42)'
    write(*, '(A)') '  --coldstart=<N>         Cold-start steps at 4K (default: 0)'
    write(*, '(A)') '  --warmup=<N>            Warm-up ramp steps 4K->T (default: 0)'
    write(*, '(A)') '  --from=<step>           Averaging start step (default: auto)'
    write(*, '(A)') '  --to=<step>             Averaging end step (default: nsteps)'
    write(*, '(A)') '  --mon=<N>               Monitor print interval (default: auto)'
    write(*, '(A)') '  --warmup_mon=<mode>     Warmup monitor: norm|freq|some (default: norm)'
    write(*, '(A)') '  --ovito=<N>             OVITO XYZ output interval, 0=off (default: 0)'
    write(*, '(A)') '  --ofile=<filename>      OVITO output filename (default: auto)'
    write(*, '(A)') '  --restart=<N>           Restart save interval, 0=off (default: 0)'
    write(*, '(A)') '  --resfile=<path>        Resume from restart file'
    write(*, '(A)') '  --libdir=<path>         Fullerene library dir (default: FullereneLib)'
    write(*, '(A)') ''
    write(*, '(A)') 'Examples:'
    write(*, '(A)') '  ./prog --temp=500 --pres=1.0 --step=50000'
    write(*, '(A)') '  ./prog --coldstart=2000 --warmup=3000 --step=20000'
    write(*, '(A)') '  ./prog --step=10000 --ovito=100'
    write(*, '(A)') '  ./prog --step=50000 --restart=5000'
    write(*, '(A)') '  ./prog --resfile=restart_LJ_serial_00010000.rst'
    stop
  end if

  ! ── コマンドライン引数の取得 ──
  call get_opt_str('fullerene', 'C60', fspec)
  call get_opt_str('crystal', 'fcc', crystal)
  call get_opt_int('cell', 3, nc)
  call get_opt_dbl('temp', 298.0d0, T)
  call get_opt_dbl('pres', 0.0d0, Pe)
  call get_opt_int('step', 10000, nsteps)
  call get_opt_dbl('dt', 1.0d0, dt)
  call get_opt_dbl('init_scale', 1.0d0, init_scale)
  call get_opt_int('seed', 42, seed)
  call get_opt_int('coldstart', 0, coldstart)
  call get_opt_int('warmup', 0, warmup)
  call get_opt_int('from', 0, avg_from)
  call get_opt_int('to', 0, avg_to)
  call get_opt_int('mon', 0, mon_int)
  call get_opt_str('warmup_mon', 'norm', warmup_mon_mode)
  call get_opt_int('ovito', 0, nrec_o)
  call get_opt_str('ofile', '', ofile_opt)
  call get_opt_int('restart', 0, nrec_rst)
  call get_opt_str('resfile', '', resfile)
  call get_opt_str('libdir', 'FullereneLib', libdir)
  start_step = 0

  ! 結晶構造タイプを大文字化
  call to_upper(crystal, st)

  if (avg_to <= 0) avg_to = nsteps
  if (avg_from <= 0) avg_from = max(1, nsteps - nsteps/4)
  total_steps = coldstart + warmup + nsteps
  gavg_from = coldstart + warmup + avg_from
  gavg_to = coldstart + warmup + avg_to

  if (avg_from < 1 .or. avg_from >= avg_to .or. avg_to > nsteps) then
    write(*, '(A)') 'Error: invalid --from/--to range'
    stop 1
  end if

  nlup = 25

  ! ── フラーレン解決 + cc1読み込み ──
  ! cc1ファイルからの読み込みを試みる。失敗したらC60フォールバック
  use_cc1 = .false.
  call resolve_fullerene(trim(fspec), trim(libdir), fpath, label)
  ! fpath が存在するか確認
  block
    logical :: cc1_exists
    inquire(file=trim(fpath), exist=cc1_exists)
    if (cc1_exists) then
      call load_cc1(trim(fpath), mol_coords, natom, Rmol, Dmol, I0, Mmol)
      use_cc1 = .true.
    else
      ! cc1が見つからない場合: C60ならフォールバック生成
      if (trim(fspec) == 'C60' .or. trim(fspec) == 'c60' .or. &
          trim(fspec) == 'Buckyball' .or. trim(fspec) == 'buckyball') then
        call generate_c60(mol_coords, natom, I0, Mmol, Rmol, Dmol)
        label = 'C60(generated)'
      else
        write(*, '(A,A)') 'Error: cc1 file not found: ', trim(fpath)
        stop 1
      end if
    end if
  end block

  ! ── カットオフ・格子定数 ──
  RMCUT_v = RCUT + 2.0d0*Rmol + 1.0d0
  RMCUT2_v = RMCUT_v * RMCUT_v
  a0 = default_a0(Dmol, trim(st), init_scale)

  ! ── Nmax計算 ──
  if (trim(st) == 'FCC') then
    Nmax = 4 * nc * nc * nc
  else if (trim(st) == 'HCP') then
    Nmax = 2 * nc * nc * nc
  else
    Nmax = 2 * nc * nc * nc
  end if

  ! ── 配列確保 ──
  allocate(pos(Nmax*3), vel(Nmax*3), omg(Nmax*3), qv(Nmax*4))
  allocate(Fv(Nmax*3), Tv(Nmax*3), lab(Nmax*natom*3))
  allocate(body(natom, 3))
  allocate(nl_count(Nmax), nl_list(Nmax*MAX_NEIGH))
  pos = 0.0d0; vel = 0.0d0; omg = 0.0d0; qv = 0.0d0
  Fv = 0.0d0; Tv = 0.0d0; lab = 0.0d0
  h = 0.0d0; hi = 0.0d0; Wm9 = 0.0d0
  nl_count = 0; nl_list = 0

  body(1:natom, :) = mol_coords(1:natom, :)

  ! ── 結晶構造の構築 ──
  if (trim(st) == 'FCC') then
    N = make_fcc(a0, nc, pos, h)
  else if (trim(st) == 'HCP') then
    N = make_hcp(a0, nc, pos, h)
  else
    N = make_bcc(a0, nc, pos, h)
  end if
  call mat_inv9(h, hi)

  ! ── OVITO出力ファイル名の決定 ──
  call build_suffix(trim(fspec), trim(crystal), nc, T, Pe, nsteps, dt, &
       init_scale, seed, sfx)
  if (len_trim(ofile_opt) > 0) then
    ovito_file = ofile_opt
  else
    call unique_file('ovito_traj_LJ_serial' // trim(sfx), '.xyz', ovito_file)
  end if

  ! ── OpenMPスレッド数取得 ──
  omp_threads = 1
!$ omp_threads = omp_get_max_threads()

  ! ── バナー表示 ──
  write(*, '(A)') '========================================================================'
  if (omp_threads > 1) then
    write(*, '(A,I0,A)') &
      '  Fullerene Crystal NPT-MD -- LJ rigid-body (OpenMP, ', omp_threads, ' threads)'
  else
    write(*, '(A)') '  Fullerene Crystal NPT-MD -- LJ rigid-body (Serial)'
  end if
  write(*, '(A)') '========================================================================'
  write(*, '(A,A,A,I0,A)') '  Fullerene       : ', trim(label), &
    ' (', natom, ' atoms/mol)'
  write(*, '(A,A,1X,I0,A,I0,A,I0,A,I0)') '  Crystal         : ', &
    trim(st), nc, 'x', nc, 'x', nc, '  Nmol=', N
  write(*, '(A,F8.3,A,F6.1,A,F8.4,A,F5.2,A)') &
    '  a0=', a0, ' A  T=', T, ' K  P=', Pe, ' GPa  dt=', dt, ' fs'
  if (coldstart > 0) write(*, '(A,I0,A,F4.1,A)') &
    '  Coldstart       : ', coldstart, ' steps at ', T_cold, ' K'
  if (warmup > 0) write(*, '(A,I0,A,F4.1,A,F6.1,A)') &
    '  Warmup          : ', warmup, ' steps (', T_cold, 'K->', T, 'K)'
  write(*, '(A,I0,A,I0,A,I0)') &
    '  Production      : ', nsteps, ' steps  avg=', avg_from, '-', avg_to
  write(*, '(A,I0,A)') '  Total           : ', total_steps, ' steps'
  write(*, '(A)') '========================================================================'
  write(*, *)

  ! ── 初期速度 (Maxwell-Boltzmann分布) ──
  T_init = T
  if (coldstart > 0 .or. warmup > 0) T_init = T_cold

  ! 乱数シードの設定
  block
    integer :: seed_size, j
    integer, allocatable :: seed_array(:)
    call random_seed(size=seed_size)
    allocate(seed_array(seed_size))
    do j = 1, seed_size
      seed_array(j) = seed + (j - 1) * 37
    end do
    call random_seed(put=seed_array)
    deallocate(seed_array)
  end block

  sv = sqrt(kB * T_init * CONV / Mmol)
  sw = sqrt(kB * T_init * CONV / I0)
  do i = 1, N
    idx = (i-1)*3
    do a = 1, 3
      vel(idx+a) = sv * gauss_rand()
      omg(idx+a) = sw * gauss_rand()
    end do
    do a = 1, 4
      qv((i-1)*4+a) = gauss_rand()
    end do
    n_norm = sqrt(qv((i-1)*4+1)**2 + qv((i-1)*4+2)**2 + &
                  qv((i-1)*4+3)**2 + qv((i-1)*4+4)**2)
    do a = 1, 4
      qv((i-1)*4+a) = qv((i-1)*4+a) / n_norm
    end do
  end do

  ! 重心速度除去
  vcm = 0.0d0
  do i = 1, N
    idx = (i-1)*3
    vcm(1) = vcm(1) + vel(idx+1)
    vcm(2) = vcm(2) + vel(idx+2)
    vcm(3) = vcm(3) + vel(idx+3)
  end do
  vcm = vcm / dble(N)
  do i = 1, N
    idx = (i-1)*3
    vel(idx+1) = vel(idx+1) - vcm(1)
    vel(idx+2) = vel(idx+2) - vcm(2)
    vel(idx+3) = vel(idx+3) - vcm(3)
  end do

  ! ── NPT初期化 ──
  call make_npt(npt, T, Pe, N)
  npt%Tt = T_init

  ! ── リスタートファイルからの再開 ──
  if (len_trim(resfile) > 0) then
    call read_restart_lj(trim(resfile), rd_istep, rd_N, rd_h, rd_npt, &
         rd_pos, rd_qv, rd_vel, rd_omg, rd_ok)
    if (rd_ok) then
      start_step = rd_istep
      h = rd_h
      npt = rd_npt
      nn = min(N, rd_N)
      do i = 1, nn*3
        pos(i) = rd_pos(i); vel(i) = rd_vel(i); omg(i) = rd_omg(i)
      end do
      do i = 1, nn*4
        qv(i) = rd_qv(i)
      end do
      write(*, '(A,I0)') '  Restarting from global step ', start_step
    end if
    if (allocated(rd_pos)) deallocate(rd_pos)
    if (allocated(rd_qv)) deallocate(rd_qv)
    if (allocated(rd_vel)) deallocate(rd_vel)
    if (allocated(rd_omg)) deallocate(rd_omg)
  end if

  ! ── 初回近傍リスト + PBC + 力計算 ──
  call mat_inv9(h, hi)
  call nlist_build_sym(pos, h, hi, N, RMCUT_v, nl_count, nl_list)
  call apply_pbc(pos, h, hi, N)
  Ep = calc_forces(Fv, Tv, Wm9, pos, qv, body, h, hi, &
                   nl_count, nl_list, N, natom, RMCUT2_v, lab)

  ! ── モニタリング間隔 ──
  prn = mon_int
  if (prn <= 0) prn = max(1, total_steps / 50)
  prn_pre = prn
  if (coldstart + warmup > 0) then
    div = 100
    if (trim(warmup_mon_mode) == 'freq') then
      div = 10
    else if (trim(warmup_mon_mode) == 'some') then
      div = 1000
    end if
    prn_pre = max(1, (coldstart + warmup) / div)
  end if

  sum_T = 0.0d0; sum_P = 0.0d0; sum_a = 0.0d0; sum_Ep = 0.0d0; nav = 0
  call cpu_time(t_start)

  ! ── OVITO出力ファイルオープン ──
  if (nrec_o > 0) then
    open(unit=40, file=trim(ovito_file), status='replace')
  end if

  write(*, '(A8,A6,A8,A10,A9,A11,A14,A8)') &
    'step', 'phase', 'T[K]', 'P[GPa]', 'a[A]', 'Ecoh[eV]', 'Ecoh[kcal/m]', 't[s]'

  stop_requested = .false.

  ! ═══ MDメインループ ═══
  do gstep = start_step + 1, total_steps
    ! フェーズ判定
    if (gstep <= coldstart) then
      phase = 'COLD'
    else if (gstep <= coldstart + warmup) then
      phase = 'WARM'
    else
      phase = 'PROD'
    end if
    cur_prn = prn
    if (gstep <= coldstart + warmup) cur_prn = prn_pre

    ! 目標温度の設定
    if (gstep <= coldstart) then
      npt%Tt = T_cold
    else if (gstep <= coldstart + warmup) then
      npt%Tt = T_cold + (T - T_cold) * dble(gstep - coldstart) / dble(warmup)
    else
      npt%Tt = T
    end if

    ! コールドスタート終了時のリセット
    if (coldstart > 0 .and. gstep == coldstart + 1) then
      npt%xi = 0.0d0; npt%Vg = 0.0d0
    end if
    if (gstep <= coldstart) npt%Vg = 0.0d0

    ! 近傍リスト再構築
    if (mod(gstep, nlup) == 0) then
      call mat_inv9(h, hi)
      call nlist_build_sym(pos, h, hi, N, RMCUT_v, nl_count, nl_list)
    end if

    ! 1ステップ
    call step_npt(pos, vel, qv, omg, Fv, Tv, Wm9, h, hi, &
         body, I0, Mmol, N, natom, RMCUT2_v, dt, npt, &
         nl_count, nl_list, lab, Ep, KE)

    ! 瞬時物理量
    kt_val = ke_trans(vel, N, Mmol)
    V_val = abs(mat_det9(h))
    Tn = inst_T(KE, npt%Nf)
    Pn = inst_P(Wm9, kt_val, V_val)

    ! COLD/WARM velocity rescaling
    if ((gstep <= coldstart .or. gstep <= coldstart + warmup) .and. Tn > 0.1d0) then
      if (gstep <= coldstart) then
        tgt = T_cold
      else
        tgt = npt%Tt
      end if
      scale_v = sqrt(max(tgt, 0.1d0) / Tn)
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx)
      do i = 1, N
        idx = (i-1)*3
        vel(idx+1) = vel(idx+1) * scale_v
        vel(idx+2) = vel(idx+2) * scale_v
        vel(idx+3) = vel(idx+3) * scale_v
        omg(idx+1) = omg(idx+1) * scale_v
        omg(idx+2) = omg(idx+2) * scale_v
        omg(idx+3) = omg(idx+3) * scale_v
      end do
      !$OMP END PARALLEL DO
      kt_val = ke_trans(vel, N, Mmol)
      KE = kt_val + ke_rot(omg, N, I0)
      Tn = inst_T(KE, npt%Nf)
      npt%xi = 0.0d0
      if (gstep <= coldstart) npt%Vg = 0.0d0
    end if

    Ec = Ep / dble(N)
    an_val = h(1) / dble(nc)

    ! 平均値の蓄積
    if (gstep >= gavg_from .and. gstep <= gavg_to) then
      sum_T = sum_T + Tn; sum_P = sum_P + Pn; sum_a = sum_a + an_val; sum_Ep = sum_Ep + Ec
      nav = nav + 1
    end if

    ! OVITO出力
    if (nrec_o > 0 .and. mod(gstep, nrec_o) == 0) then
      call write_ovito(40, gstep, dt, pos, vel, qv, body, h, N, natom)
      flush(40)
    end if

    ! リスタート保存
    if (nrec_rst > 0 .and. &
        (mod(gstep, nrec_rst) == 0 .or. gstep == total_steps)) then
      call restart_filename(trim(ovito_file), gstep, total_steps, rfn)
      call write_restart_lj(trim(rfn), gstep, trim(st), nc, T, Pe, &
           nsteps, dt, seed, trim(fspec), init_scale, &
           h, npt, pos, qv, vel, omg, N, natom)
      if (stop_requested) then
        write(*, '(A,I0,A)') &
          '  *** Stopped at restart checkpoint (step ', gstep, ') ***'
        exit
      end if
    end if

    ! モニタリング出力 + 停止制御
    if (mod(gstep, cur_prn) == 0 .or. gstep == total_steps) then
      ! abort.md チェック
      if (dir_exists('abort.md')) then
        write(*, '(A,I0,A)') '  *** abort.md detected at step ', gstep, ' ***'
        if (nrec_rst > 0) then
          call restart_filename(trim(ovito_file), gstep, total_steps, rfn)
          call write_restart_lj(trim(rfn), gstep, trim(st), nc, T, Pe, &
               nsteps, dt, seed, trim(fspec), init_scale, &
               h, npt, pos, qv, vel, omg, N, natom)
        end if
        exit
      end if

      ! stop.md チェック
      if (.not. stop_requested .and. dir_exists('stop.md')) then
        stop_requested = .true.
        write(*, '(A,I0,A)') &
          '  *** stop.md detected at step ', gstep, &
          ' -- will stop at next checkpoint ***'
      end if

      call cpu_time(t_now)
      elapsed = t_now - t_start
      write(*, '(I8,A6,F8.1,F10.3,F9.3,F11.5,F14.4,F8.0)') &
        gstep, phase, Tn, Pn, an_val, Ec, Ec*eV2kcalmol, elapsed
    end if
  end do

  ! ── OVITO出力ファイルクローズ ──
  if (nrec_o > 0) close(40)

  ! ── 平均値出力 ──
  if (nav > 0) then
    write(*, *)
    write(*, '(A)') '========================================================================'
    write(*, '(A,I0,A,F7.2,A,F8.4,A,F8.4,A,F10.5,A)') &
      '  Averages (', nav, ' samples): T=', sum_T/dble(nav), &
      ' K  P=', sum_P/dble(nav), ' GPa  a=', sum_a/dble(nav), &
      ' A  Ecoh=', sum_Ep/dble(nav), ' eV'
    write(*, '(A)') '========================================================================'
  end if
  call cpu_time(t_now)
  write(*, '(A,F8.1,A)') '  Done (', t_now - t_start, ' sec)'

  ! ── 解放 ──
  deallocate(pos, vel, omg, qv, Fv, Tv, lab, body, nl_count, nl_list)

end program fuller_LJ_npt_md
