! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2025, Takeshi Nishikawa
!===========================================================================
!  fuller_airebo_npt_md_serial_omp_acc.f90
!  Fullerene Crystal NPT-MD
!  (AIREBO: REBO-II + LJ, all-atom, Serial / OpenMP, Fortran 95)
!
!  REBO-II (Brenner 2002) + LJ (Stuart 2000) / NH+PR / 3D PBC
!
!  コンパイル:
!    gfortran -O3 -o fuller_airebo_npt_md_serial fuller_airebo_npt_md_serial_omp_acc.f90
!    gfortran -O3 -fopenmp -o fuller_airebo_npt_md_omp fuller_airebo_npt_md_serial_omp_acc.f90
!
!  実行時オプション (全て --key=value 形式):
!    --fullerene=<名前>      フラーレン種 (デフォルト: C60)
!    --crystal=<fcc|hcp|bcc> 結晶構造 (デフォルト: fcc)
!    --cell=<nc>             単位胞の繰り返し数 (デフォルト: 3)
!    --temp=<T_K>            目標温度 [K] (デフォルト: 298.0)
!    --pres=<P_GPa>          目標圧力 [GPa] (デフォルト: 0.0)
!    --step=<N>              本計算ステップ数 (デフォルト: 10000)
!    --dt=<fs>               時間刻み [fs] (デフォルト: 0.5)
!    --init_scale=<s>        格子定数スケール因子 (デフォルト: 1.0)
!    --seed=<n>              乱数シード (デフォルト: 42)
!    --coldstart=<N>         低温(4K)ステップ数 (デフォルト: 0)
!    --warmup=<N>            昇温ステップ数 4K→T (デフォルト: 0)
!    --from=<step>           平均開始ステップ (デフォルト: 本計算の3/4地点)
!    --to=<step>             平均終了ステップ (デフォルト: nsteps)
!    --mon=<N>               モニタリング出力間隔 (デフォルト: 自動)
!    --ovito=<N>             OVITO XYZ出力間隔 (0=無効, デフォルト: 0)
!    --restart=<N>           リスタート保存間隔 (0=無効, デフォルト: 0)
!    --resfile=<path>        リスタートファイルから再開
!  単位系: A, amu, eV, fs, K, GPa
!===========================================================================
module fuller_airebo_mod
  implicit none

  ! ═══ 物理定数 ═══
  double precision, parameter :: CONV = 9.64853321d-3
  double precision, parameter :: kB_ = 8.617333262d-5
  double precision, parameter :: eV2GPa = 160.21766208d0
  double precision, parameter :: eV2kcalmol = 23.06054783d0
  double precision, parameter :: PI_ = 3.14159265358979323846d0
  double precision, parameter :: mC_ = 12.011d0

  ! ═══ REBO-II C-C パラメータ ═══
  double precision, parameter :: Q_CC = 0.3134602960833d0
  double precision, parameter :: A_CC = 10953.544162170d0
  double precision, parameter :: alpha_CC = 4.7465390606595d0
  double precision, parameter :: B1_CC = 12388.79197798d0
  double precision, parameter :: beta1_CC = 4.7204523127d0
  double precision, parameter :: B2_CC = 17.56740646509d0
  double precision, parameter :: beta2_CC = 1.4332132499d0
  double precision, parameter :: B3_CC = 30.71493208065d0
  double precision, parameter :: beta3_CC = 1.3826912506d0
  double precision, parameter :: Dmin_CC = 1.7d0
  double precision, parameter :: Dmax_CC = 2.0d0
  double precision, parameter :: BO_DELTA = 0.5d0
  double precision, parameter :: G_a0 = 0.00020813d0
  double precision, parameter :: G_c0 = 330.0d0
  double precision, parameter :: G_d0 = 3.5d0
  double precision, parameter :: REBO_RCUT = Dmax_CC + 0.3d0

  ! ═══ LJ 分子間パラメータ ═══
  double precision, parameter :: sig_LJ = 3.40d0
  double precision, parameter :: eps_LJ_ = 2.84d-3
  double precision, parameter :: LJ_RCUT = 3.0d0 * sig_LJ
  double precision, parameter :: LJ_RCUT2 = LJ_RCUT * LJ_RCUT
  double precision, parameter :: sig2_LJ = sig_LJ * sig_LJ
  double precision, parameter :: sr_v_ = sig_LJ / LJ_RCUT
  double precision, parameter :: sr2_ = sr_v_ * sr_v_
  double precision, parameter :: sr6_ = sr2_ * sr2_ * sr2_
  double precision, parameter :: LJ_VSHFT = 4.0d0*eps_LJ_*(sr6_*sr6_ - sr6_)

  ! ═══ 近傍リストサイズ ═══
  integer, parameter :: MAX_REBO_NEIGH = 12
  integer, parameter :: MAX_LJ_NEIGH = 400
  integer, parameter :: MAX_NATOM = 84

  ! ═══ NPT状態 ═══
  type :: NPTState
    double precision :: xi, Q, Vg(9), W_, Pe, Tt
    integer :: Nf
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
    hi(1)=id*(h(5)*h(9)-h(6)*h(8)); hi(2)=id*(h(3)*h(8)-h(2)*h(9))
    hi(3)=id*(h(2)*h(6)-h(3)*h(5)); hi(4)=id*(h(6)*h(7)-h(4)*h(9))
    hi(5)=id*(h(1)*h(9)-h(3)*h(7)); hi(6)=id*(h(3)*h(4)-h(1)*h(6))
    hi(7)=id*(h(4)*h(8)-h(5)*h(7)); hi(8)=id*(h(2)*h(7)-h(1)*h(8))
    hi(9)=id*(h(1)*h(5)-h(2)*h(4))
  end subroutine

  ! ═══ 最小イメージ ═══
  subroutine mimg9(dx, dy, dz, hi, h)
    double precision, intent(inout) :: dx, dy, dz
    double precision, intent(in) :: hi(9), h(9)
    double precision :: s0, s1, s2
    s0=hi(1)*dx+hi(2)*dy+hi(3)*dz
    s1=hi(4)*dx+hi(5)*dy+hi(6)*dz
    s2=hi(7)*dx+hi(8)*dy+hi(9)*dz
    s0=s0-anint(s0); s1=s1-anint(s1); s2=s2-anint(s2)
    dx=h(1)*s0+h(2)*s1+h(3)*s2
    dy=h(4)*s0+h(5)*s1+h(6)*s2
    dz=h(7)*s0+h(8)*s1+h(9)*s2
  end subroutine

  ! ═══ PBC適用 ═══
  subroutine apply_pbc(pos, h, hi, N)
    double precision, intent(inout) :: pos(:)
    double precision, intent(in) :: h(9), hi(9)
    integer, intent(in) :: N
    double precision :: px, py, pz, s0, s1, s2
    integer :: i, idx
    !$OMP PARALLEL DO PRIVATE(i,idx,px,py,pz,s0,s1,s2)
    do i = 1, N
      idx = (i-1)*3
      px=pos(idx+1); py=pos(idx+2); pz=pos(idx+3)
      s0=hi(1)*px+hi(2)*py+hi(3)*pz
      s1=hi(4)*px+hi(5)*py+hi(6)*pz
      s2=hi(7)*px+hi(8)*py+hi(9)*pz
      s0=s0-floor(s0); s1=s1-floor(s1); s2=s2-floor(s2)
      pos(idx+1)=h(1)*s0+h(2)*s1+h(3)*s2
      pos(idx+2)=h(4)*s0+h(5)*s1+h(6)*s2
      pos(idx+3)=h(7)*s0+h(8)*s1+h(9)*s2
    end do
    !$OMP END PARALLEL DO
  end subroutine

  ! ═══ REBO カットオフ関数 fc(r) ═══
  double precision function fc_d(r, dmin, dmax)
    double precision, intent(in) :: r, dmin, dmax
    if (r <= dmin) then; fc_d = 1.0d0
    else if (r >= dmax) then; fc_d = 0.0d0
    else; fc_d = 0.5d0*(1.0d0 + cos(PI_*(r-dmin)/(dmax-dmin))); end if
  end function

  ! ═══ REBO カットオフ関数の微分 dfc/dr ═══
  double precision function dfc_d(r, dmin, dmax)
    double precision, intent(in) :: r, dmin, dmax
    if (r <= dmin .or. r >= dmax) then; dfc_d = 0.0d0
    else; dfc_d = -0.5d0*PI_/(dmax-dmin)*sin(PI_*(r-dmin)/(dmax-dmin)); end if
  end function

  ! ═══ REBO 反発ポテンシャル VR(r) ═══
  double precision function VR_CC(r)
    double precision, intent(in) :: r
    VR_CC = (1.0d0 + Q_CC/r) * A_CC * exp(-alpha_CC*r)
  end function

  ! ═══ REBO dVR/dr ═══
  double precision function dVR_CC(r)
    double precision, intent(in) :: r
    double precision :: ex
    ex = A_CC*exp(-alpha_CC*r)
    dVR_CC = (-Q_CC/(r*r))*ex + (1.0d0 + Q_CC/r)*(-alpha_CC)*ex
  end function

  ! ═══ REBO 引力ポテンシャル VA(r) ═══
  double precision function VA_CC(r)
    double precision, intent(in) :: r
    VA_CC = B1_CC*exp(-beta1_CC*r) + B2_CC*exp(-beta2_CC*r) &
          + B3_CC*exp(-beta3_CC*r)
  end function

  ! ═══ REBO dVA/dr ═══
  double precision function dVA_CC(r)
    double precision, intent(in) :: r
    dVA_CC = -beta1_CC*B1_CC*exp(-beta1_CC*r) &
             -beta2_CC*B2_CC*exp(-beta2_CC*r) &
             -beta3_CC*B3_CC*exp(-beta3_CC*r)
  end function

  ! ═══ REBO 角度関数 G(cos_theta) ═══
  double precision function G_C(x)
    double precision, intent(in) :: x
    double precision :: c2, d2, hv
    c2=G_c0*G_c0; d2=G_d0*G_d0; hv=1.0d0+x
    G_C = G_a0*(1.0d0 + c2/d2 - c2/(d2 + hv*hv))
  end function

  ! ═══ REBO dG/d(cos_theta) ═══
  double precision function dG_C(x)
    double precision, intent(in) :: x
    double precision :: c2, d2, hv, dn
    c2=G_c0*G_c0; d2=G_d0*G_d0; hv=1.0d0+x; dn=d2+hv*hv
    dG_C = G_a0*2.0d0*c2*hv/(dn*dn)
  end function

  ! ═══ cc1ファイル読み込み ═══
  subroutine load_cc1(path, coords, natom, Rmol, Dmol)
    character(len=*), intent(in) :: path
    double precision, intent(out) :: coords(MAX_NATOM,3), Rmol, Dmol
    integer, intent(out) :: natom
    double precision :: cm(3), r, dx, dy, dz, d, x, y, z
    integer :: i, j, ios, idx_a, fl
    character(len=512) :: line
    character(len=8) :: elem
    open(unit=20, file=path, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: cannot open ', trim(path); stop 1
    end if
    read(20,*) natom
    if (natom > MAX_NATOM) then
      write(*,*) 'Error: natom > MAX_NATOM'; stop 1
    end if
    do i = 1, natom
      read(20,'(A)',iostat=ios) line
      do while (len_trim(line) == 0 .and. ios == 0)
        read(20,'(A)',iostat=ios) line
      end do
      read(line,*,iostat=ios) elem, idx_a, x, y, z, fl
      coords(i,1) = x; coords(i,2) = y; coords(i,3) = z
    end do
    close(20)
    cm = 0.0d0
    do i = 1, natom
      cm(1)=cm(1)+coords(i,1); cm(2)=cm(2)+coords(i,2); cm(3)=cm(3)+coords(i,3)
    end do
    cm = cm / dble(natom)
    do i = 1, natom
      coords(i,1)=coords(i,1)-cm(1); coords(i,2)=coords(i,2)-cm(2)
      coords(i,3)=coords(i,3)-cm(3)
    end do
    Rmol = 0.0d0; Dmol = 0.0d0
    do i = 1, natom
      r = sqrt(coords(i,1)**2 + coords(i,2)**2 + coords(i,3)**2)
      if (r > Rmol) Rmol = r
      do j = i+1, natom
        dx=coords(i,1)-coords(j,1); dy=coords(i,2)-coords(j,2)
        dz=coords(i,3)-coords(j,3)
        d = sqrt(dx*dx+dy*dy+dz*dz)
        if (d > Dmol) Dmol = d
      end do
    end do
  end subroutine

  ! ═══ C60座標の生成 (cc1ファイルなしの場合のフォールバック) ═══
  subroutine generate_c60_airebo(coords, natom, Rmol, Dmol)
    double precision, intent(out) :: coords(MAX_NATOM,3)
    integer, intent(out) :: natom
    double precision, intent(out) :: Rmol, Dmol
    double precision :: phi, tmp(60,3), cm(3), r2, r, dx, dy, dz, d
    integer :: n, p, s1, s2, s3, cyc(3,3), signs(2), i, j
    natom = 60; phi = (1.0d0+sqrt(5.0d0))/2.0d0
    signs(1) = -1; signs(2) = 1
    cyc(1,1)=1;cyc(1,2)=2;cyc(1,3)=3
    cyc(2,1)=2;cyc(2,2)=3;cyc(2,3)=1
    cyc(3,1)=3;cyc(3,2)=1;cyc(3,3)=2
    tmp = 0.0d0; n = 0
    do p=1,3; do s2=1,2; do s3=1,2; n=n+1
      tmp(n,cyc(p,1))=0.0d0
      tmp(n,cyc(p,2))=dble(signs(s2))
      tmp(n,cyc(p,3))=dble(signs(s3))*3.0d0*phi
    end do; end do; end do
    do p=1,3; do s1=1,2; do s2=1,2; do s3=1,2; n=n+1
      tmp(n,cyc(p,1))=dble(signs(s1))*2.0d0
      tmp(n,cyc(p,2))=dble(signs(s2))*(1.0d0+2.0d0*phi)
      tmp(n,cyc(p,3))=dble(signs(s3))*phi
    end do; end do; end do; end do
    do p=1,3; do s1=1,2; do s2=1,2; do s3=1,2; n=n+1
      tmp(n,cyc(p,1))=dble(signs(s1))
      tmp(n,cyc(p,2))=dble(signs(s2))*(2.0d0+phi)
      tmp(n,cyc(p,3))=dble(signs(s3))*2.0d0*phi
    end do; end do; end do; end do
    cm = 0.0d0
    do n=1,60; cm(1)=cm(1)+tmp(n,1); cm(2)=cm(2)+tmp(n,2); cm(3)=cm(3)+tmp(n,3); end do
    cm = cm/60.0d0
    do n=1,60; tmp(n,:)=(tmp(n,:)-cm)*0.72d0; end do
    Rmol = 0.0d0; Dmol = 0.0d0
    do i=1,60
      r2=tmp(i,1)**2+tmp(i,2)**2+tmp(i,3)**2; r=sqrt(r2)
      if (r > Rmol) Rmol = r
      do j=i+1,60
        dx=tmp(i,1)-tmp(j,1); dy=tmp(i,2)-tmp(j,2); dz=tmp(i,3)-tmp(j,3)
        d=sqrt(dx*dx+dy*dy+dz*dz); if (d > Dmol) Dmol = d
      end do
    end do
    coords(1:60,:) = tmp
  end subroutine

  ! ═══ FCC結晶構造生成 ═══
  function make_fcc(a, nc, pos, h) result(Nmol)
    double precision, intent(in) :: a
    integer, intent(in) :: nc
    double precision, intent(out) :: pos(:), h(9)
    integer :: Nmol, ix, iy, iz, b, idx
    double precision :: bas(4,3)
    bas(1,:)=(/0d0,0d0,0d0/); bas(2,:)=(/.5d0*a,.5d0*a,0d0/)
    bas(3,:)=(/.5d0*a,0d0,.5d0*a/); bas(4,:)=(/0d0,.5d0*a,.5d0*a/)
    Nmol = 0
    do ix=0,nc-1; do iy=0,nc-1; do iz=0,nc-1; do b=1,4
      Nmol=Nmol+1; idx=(Nmol-1)*3
      pos(idx+1)=a*dble(ix)+bas(b,1)
      pos(idx+2)=a*dble(iy)+bas(b,2)
      pos(idx+3)=a*dble(iz)+bas(b,3)
    end do; end do; end do; end do
    h=0; h(1)=dble(nc)*a; h(5)=dble(nc)*a; h(9)=dble(nc)*a
  end function

  ! ═══ HCP結晶構造生成 ═══
  function make_hcp(a, nc, pos, h) result(Nmol)
    double precision, intent(in) :: a
    integer, intent(in) :: nc
    double precision, intent(out) :: pos(:), h(9)
    integer :: Nmol, ix, iy, iz, b, idx
    double precision :: c, a1(3), a2(3), a3(3), bas(2,3), fx, fy, fz
    c = a*sqrt(8.0d0/3.0d0)
    a1 = (/a, 0.0d0, 0.0d0/)
    a2 = (/a/2.0d0, a*sqrt(3.0d0)/2.0d0, 0.0d0/)
    a3 = (/0.0d0, 0.0d0, c/)
    bas(1,:) = (/0.0d0, 0.0d0, 0.0d0/)
    bas(2,:) = (/1.0d0/3.0d0, 2.0d0/3.0d0, 0.5d0/)
    Nmol = 0
    do ix=0,nc-1; do iy=0,nc-1; do iz=0,nc-1; do b=1,2
      Nmol=Nmol+1; idx=(Nmol-1)*3
      fx=dble(ix)+bas(b,1); fy=dble(iy)+bas(b,2); fz=dble(iz)+bas(b,3)
      pos(idx+1)=fx*a1(1)+fy*a2(1)+fz*a3(1)
      pos(idx+2)=fx*a1(2)+fy*a2(2)+fz*a3(2)
      pos(idx+3)=fx*a1(3)+fy*a2(3)+fz*a3(3)
    end do; end do; end do; end do
    h=0; h(1)=dble(nc)*a1(1); h(4)=dble(nc)*a1(2)
    h(2)=dble(nc)*a2(1); h(5)=dble(nc)*a2(2); h(9)=dble(nc)*a3(3)
  end function

  ! ═══ BCC結晶構造生成 ═══
  function make_bcc(a, nc, pos, h) result(Nmol)
    double precision, intent(in) :: a
    integer, intent(in) :: nc
    double precision, intent(out) :: pos(:), h(9)
    integer :: Nmol, ix, iy, iz, b, idx
    double precision :: bas(2,3)
    bas(1,:)=(/0d0,0d0,0d0/); bas(2,:)=(/.5d0*a,.5d0*a,.5d0*a/)
    Nmol = 0
    do ix=0,nc-1; do iy=0,nc-1; do iz=0,nc-1; do b=1,2
      Nmol=Nmol+1; idx=(Nmol-1)*3
      pos(idx+1)=a*dble(ix)+bas(b,1)
      pos(idx+2)=a*dble(iy)+bas(b,2)
      pos(idx+3)=a*dble(iz)+bas(b,3)
    end do; end do; end do; end do
    h=0; h(1)=dble(nc)*a; h(5)=dble(nc)*a; h(9)=dble(nc)*a
  end function

  ! ═══ デフォルト格子定数 ═══
  double precision function default_a0(dmax, cst, s)
    double precision, intent(in) :: dmax, s
    character(len=*), intent(in) :: cst
    double precision :: m
    m = 1.4d0
    if (cst == "FCC") then
      default_a0 = dmax*sqrt(2.0d0)*m*s
    else if (cst == "HCP") then
      default_a0 = dmax*m*s
    else
      default_a0 = dmax*2.0d0/sqrt(3.0d0)*m*s
    end if
  end function

  ! ═══ REBO近傍リスト構築 (分子内のみ) ═══
  subroutine build_nlist_rebo(pos, h, hi, Na, mol_id, nlc, nll)
    double precision, intent(in) :: pos(:), h(9), hi(9)
    integer, intent(in) :: Na, mol_id(:)
    integer, intent(out) :: nlc(:), nll(:)
    double precision :: rc2, dx, dy, dz, r2
    integer :: i, j
    rc2 = REBO_RCUT*REBO_RCUT
    nlc(1:Na) = 0
    do i = 1, Na-1
      do j = i+1, Na
        if (mol_id(j) /= mol_id(i)) cycle
        dx=pos((j-1)*3+1)-pos((i-1)*3+1)
        dy=pos((j-1)*3+2)-pos((i-1)*3+2)
        dz=pos((j-1)*3+3)-pos((i-1)*3+3)
        call mimg9(dx,dy,dz,hi,h)
        r2 = dx*dx+dy*dy+dz*dz
        if (r2 < rc2) then
          if (nlc(i) < MAX_REBO_NEIGH) then
            nlc(i)=nlc(i)+1; nll((i-1)*MAX_REBO_NEIGH+nlc(i))=j
          end if
          if (nlc(j) < MAX_REBO_NEIGH) then
            nlc(j)=nlc(j)+1; nll((j-1)*MAX_REBO_NEIGH+nlc(j))=i
          end if
        end if
      end do
    end do
  end subroutine

  ! ═══ LJ近傍リスト構築 (分子間のみ, half-list) ═══
  subroutine build_nlist_lj(pos, h, hi, Na, mol_id, nlc, nll)
    double precision, intent(in) :: pos(:), h(9), hi(9)
    integer, intent(in) :: Na, mol_id(:)
    integer, intent(out) :: nlc(:), nll(:)
    double precision :: rc2, dx, dy, dz, r2
    integer :: i, j
    rc2 = (LJ_RCUT+2.0d0)**2
    nlc(1:Na) = 0
    do i = 1, Na-1
      do j = i+1, Na
        if (mol_id(j) == mol_id(i)) cycle
        dx=pos((j-1)*3+1)-pos((i-1)*3+1)
        dy=pos((j-1)*3+2)-pos((i-1)*3+2)
        dz=pos((j-1)*3+3)-pos((i-1)*3+3)
        call mimg9(dx,dy,dz,hi,h)
        r2 = dx*dx+dy*dy+dz*dz
        if (r2 < rc2 .and. nlc(i) < MAX_LJ_NEIGH) then
          nlc(i)=nlc(i)+1; nll((i-1)*MAX_LJ_NEIGH+nlc(i))=j
        end if
      end do
    end do
  end subroutine

  ! ═══ REBO-II 力計算 (3体力含む) ═══
  subroutine compute_rebo(F, vir9, pos, h, hi, nlc, nll, Na, Ep_rebo)
    double precision, intent(inout) :: F(:), vir9(9)
    double precision, intent(in) :: pos(:), h(9), hi(9)
    integer, intent(in) :: nlc(:), nll(:), Na
    double precision, intent(out) :: Ep_rebo
    double precision :: dx, dy, dz, rij, fcut_v, dfcut_v, vr, dvr, va, dva, rij_inv
    double precision :: rhat0, rhat1, rhat2
    double precision :: Gs_ij, Gs_ji, bij, bji, bbar, fpair, costh, fc_ik_v, rik
    double precision :: dkx, dky, dkz, dlx, dly, dlz, rjl, fc_jl_v
    double precision :: dbp, vh, dfc_ik_v, rik_inv
    double precision :: rhat_ik0, rhat_ik1, rhat_ik2
    double precision :: gv, dgv, coeff_v, fk1, dc_v, fk_v
    double precision :: rh(3), rk(3)
    double precision :: dfc_jl_v, rjl_inv
    double precision :: rhat_jl0, rhat_jl1, rhat_jl2
    double precision :: fl1, fl_v, rji(3), rl(3)
    double precision :: da, db
    integer :: i, j, k, l, jn, kn, ln, nni, nnj, aa, bb

    Ep_rebo = 0.0d0

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,4) REDUCTION(+:Ep_rebo) &
    !$OMP PRIVATE(i,nni,jn,j,dx,dy,dz,rij,fcut_v,dfcut_v,vr,dvr,va,dva, &
    !$OMP   rij_inv,rhat0,rhat1,rhat2,Gs_ij,kn,k,dkx,dky,dkz,rik, &
    !$OMP   fc_ik_v,costh,Gs_ji,nnj,ln,l,dlx,dly,dlz,rjl,fc_jl_v, &
    !$OMP   bij,bji,bbar,fpair,dbp,vh,dfc_ik_v,rik_inv, &
    !$OMP   rhat_ik0,rhat_ik1,rhat_ik2,gv,dgv,coeff_v,fk1,dc_v,fk_v, &
    !$OMP   rh,rk,dfc_jl_v,rjl_inv,rhat_jl0,rhat_jl1,rhat_jl2, &
    !$OMP   fl1,fl_v,rji,rl,da,db,aa,bb)
    do i = 1, Na
      nni = nlc(i)
      do jn = 1, nni
        j = nll((i-1)*MAX_REBO_NEIGH+jn)
        if (j <= i) cycle

        dx=pos((j-1)*3+1)-pos((i-1)*3+1)
        dy=pos((j-1)*3+2)-pos((i-1)*3+2)
        dz=pos((j-1)*3+3)-pos((i-1)*3+3)
        call mimg9(dx,dy,dz,hi,h)
        rij = sqrt(dx*dx+dy*dy+dz*dz)
        if (rij > Dmax_CC) cycle

        fcut_v = fc_d(rij, Dmin_CC, Dmax_CC)
        dfcut_v = dfc_d(rij, Dmin_CC, Dmax_CC)
        if (fcut_v < 1.0d-15 .and. dfcut_v == 0.0d0) cycle

        vr = VR_CC(rij); dvr = dVR_CC(rij)
        va = VA_CC(rij); dva = dVA_CC(rij)
        rij_inv = 1.0d0/rij
        rhat0 = dx*rij_inv; rhat1 = dy*rij_inv; rhat2 = dz*rij_inv

        ! ── b_ij の計算 ──
        Gs_ij = 0.0d0
        do kn = 1, nni
          k = nll((i-1)*MAX_REBO_NEIGH+kn); if (k == j) cycle
          dkx=pos((k-1)*3+1)-pos((i-1)*3+1)
          dky=pos((k-1)*3+2)-pos((i-1)*3+2)
          dkz=pos((k-1)*3+3)-pos((i-1)*3+3)
          call mimg9(dkx,dky,dkz,hi,h)
          rik = sqrt(dkx*dkx+dky*dky+dkz*dkz)
          if (rik > Dmax_CC) cycle
          fc_ik_v = fc_d(rik, Dmin_CC, Dmax_CC)
          if (fc_ik_v < 1.0d-15) cycle
          costh = (dx*dkx+dy*dky+dz*dkz)/(rij*rik)
          costh = max(-1.0d0, min(1.0d0, costh))
          Gs_ij = Gs_ij + fc_ik_v*G_C(costh)
        end do
        bij = (1.0d0+Gs_ij)**(-BO_DELTA)

        ! ── b_ji の計算 ──
        Gs_ji = 0.0d0; nnj = nlc(j)
        do ln = 1, nnj
          l = nll((j-1)*MAX_REBO_NEIGH+ln); if (l == i) cycle
          dlx=pos((l-1)*3+1)-pos((j-1)*3+1)
          dly=pos((l-1)*3+2)-pos((j-1)*3+2)
          dlz=pos((l-1)*3+3)-pos((j-1)*3+3)
          call mimg9(dlx,dly,dlz,hi,h)
          rjl = sqrt(dlx*dlx+dly*dly+dlz*dlz)
          if (rjl > Dmax_CC) cycle
          fc_jl_v = fc_d(rjl, Dmin_CC, Dmax_CC)
          if (fc_jl_v < 1.0d-15) cycle
          costh = (-dx*dlx-dy*dly-dz*dlz)/(rij*rjl)
          costh = max(-1.0d0, min(1.0d0, costh))
          Gs_ji = Gs_ji + fc_jl_v*G_C(costh)
        end do
        bji = (1.0d0+Gs_ji)**(-BO_DELTA)
        bbar = 0.5d0*(bij+bji)

        ! ── エネルギー・ペア力 ──
        Ep_rebo = Ep_rebo + fcut_v*(vr - bbar*va)
        fpair = (dfcut_v*(vr - bbar*va) + fcut_v*(dvr - bbar*dva))*rij_inv

        ! ── ペア力: F[i]は安全、F[j]はatomic ──
        F((i-1)*3+1) = F((i-1)*3+1) + fpair*dx
        F((i-1)*3+2) = F((i-1)*3+2) + fpair*dy
        F((i-1)*3+3) = F((i-1)*3+3) + fpair*dz
        !$OMP ATOMIC
        F((j-1)*3+1) = F((j-1)*3+1) - fpair*dx
        !$OMP ATOMIC
        F((j-1)*3+2) = F((j-1)*3+2) - fpair*dy
        !$OMP ATOMIC
        F((j-1)*3+3) = F((j-1)*3+3) - fpair*dz

        ! ── ペアビリアル (全9成分, atomic) ──
        do aa = 1, 3
          if (aa == 1) then; da = dx
          else if (aa == 2) then; da = dy
          else; da = dz; end if
          do bb = 1, 3
            if (bb == 1) then; db = dx
            else if (bb == 2) then; db = dy
            else; db = dz; end if
            !$OMP ATOMIC
            vir9((aa-1)*3+bb) = vir9((aa-1)*3+bb) - da*fpair*db
          end do
        end do

        ! ── 3体力: db_ij/dr_k ──
        if (abs(Gs_ij) > 1.0d-20 .and. va > 1.0d-20) then
          dbp = -BO_DELTA*(1.0d0+Gs_ij)**(-BO_DELTA-1.0d0)
          vh = 0.5d0*fcut_v*va
          do kn = 1, nni
            k = nll((i-1)*MAX_REBO_NEIGH+kn); if (k == j) cycle
            dkx=pos((k-1)*3+1)-pos((i-1)*3+1)
            dky=pos((k-1)*3+2)-pos((i-1)*3+2)
            dkz=pos((k-1)*3+3)-pos((i-1)*3+3)
            call mimg9(dkx,dky,dkz,hi,h)
            rik = sqrt(dkx*dkx+dky*dky+dkz*dkz)
            if (rik > Dmax_CC) cycle
            fc_ik_v = fc_d(rik, Dmin_CC, Dmax_CC)
            if (fc_ik_v < 1.0d-15) cycle
            dfc_ik_v = dfc_d(rik, Dmin_CC, Dmax_CC)
            rik_inv = 1.0d0/rik
            rhat_ik0 = dkx*rik_inv; rhat_ik1 = dky*rik_inv; rhat_ik2 = dkz*rik_inv
            costh = (dx*dkx+dy*dky+dz*dkz)/(rij*rik)
            costh = max(-1.0d0, min(1.0d0, costh))
            gv = G_C(costh); dgv = dG_C(costh)
            coeff_v = -vh*dbp
            rh(1)=rhat0; rh(2)=rhat1; rh(3)=rhat2
            rk(1)=rhat_ik0; rk(2)=rhat_ik1; rk(3)=rhat_ik2
            do aa = 1, 3
              fk1 = coeff_v*dfc_ik_v*gv*rk(aa)
              dc_v = (rh(aa) - costh*rk(aa))*rik_inv
              fk_v = fk1 + coeff_v*fc_ik_v*dgv*dc_v
              !$OMP ATOMIC
              F((k-1)*3+aa) = F((k-1)*3+aa) + fk_v
              F((i-1)*3+aa) = F((i-1)*3+aa) - fk_v
            end do
          end do
        end if

        ! ── 3体力: db_ji/dr_l ──
        if (abs(Gs_ji) > 1.0d-20 .and. va > 1.0d-20) then
          dbp = -BO_DELTA*(1.0d0+Gs_ji)**(-BO_DELTA-1.0d0)
          vh = 0.5d0*fcut_v*va
          do ln = 1, nnj
            l = nll((j-1)*MAX_REBO_NEIGH+ln); if (l == i) cycle
            dlx=pos((l-1)*3+1)-pos((j-1)*3+1)
            dly=pos((l-1)*3+2)-pos((j-1)*3+2)
            dlz=pos((l-1)*3+3)-pos((j-1)*3+3)
            call mimg9(dlx,dly,dlz,hi,h)
            rjl = sqrt(dlx*dlx+dly*dly+dlz*dlz)
            if (rjl > Dmax_CC) cycle
            fc_jl_v = fc_d(rjl, Dmin_CC, Dmax_CC)
            if (fc_jl_v < 1.0d-15) cycle
            dfc_jl_v = dfc_d(rjl, Dmin_CC, Dmax_CC)
            rjl_inv = 1.0d0/rjl
            rhat_jl0 = dlx*rjl_inv; rhat_jl1 = dly*rjl_inv; rhat_jl2 = dlz*rjl_inv
            costh = (-dx*dlx-dy*dly-dz*dlz)/(rij*rjl)
            costh = max(-1.0d0, min(1.0d0, costh))
            gv = G_C(costh); dgv = dG_C(costh)
            coeff_v = -vh*dbp
            rji(1)=-rhat0; rji(2)=-rhat1; rji(3)=-rhat2
            rl(1)=rhat_jl0; rl(2)=rhat_jl1; rl(3)=rhat_jl2
            do aa = 1, 3
              fl1 = coeff_v*dfc_jl_v*gv*rl(aa)
              dc_v = (rji(aa) - costh*rl(aa))*rjl_inv
              fl_v = fl1 + coeff_v*fc_jl_v*dgv*dc_v
              !$OMP ATOMIC
              F((l-1)*3+aa) = F((l-1)*3+aa) + fl_v
              !$OMP ATOMIC
              F((j-1)*3+aa) = F((j-1)*3+aa) - fl_v
            end do
          end do
        end if

      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine

  ! ═══ LJ分子間力計算 (half-list) ═══
  subroutine compute_lj(F, vir9, pos, h, hi, nlc, nll, Na, Ep_lj)
    double precision, intent(inout) :: F(:), vir9(9)
    double precision, intent(in) :: pos(:), h(9), hi(9)
    integer, intent(in) :: nlc(:), nll(:), Na
    double precision, intent(out) :: Ep_lj
    double precision :: dx, dy, dz, r2, ri2, sr2, sr6, sr12, fm, da, db
    integer :: i, jn, j, nni, aa, bb

    Ep_lj = 0.0d0
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,4) REDUCTION(+:Ep_lj) &
    !$OMP PRIVATE(i,nni,jn,j,dx,dy,dz,r2,ri2,sr2,sr6,sr12,fm,da,db,aa,bb)
    do i = 1, Na
      nni = nlc(i)
      do jn = 1, nni
        j = nll((i-1)*MAX_LJ_NEIGH+jn)
        dx=pos((j-1)*3+1)-pos((i-1)*3+1)
        dy=pos((j-1)*3+2)-pos((i-1)*3+2)
        dz=pos((j-1)*3+3)-pos((i-1)*3+3)
        call mimg9(dx,dy,dz,hi,h)
        r2 = dx*dx+dy*dy+dz*dz
        if (r2 > LJ_RCUT2) cycle
        if (r2 < 0.25d0) r2 = 0.25d0
        ri2=1.0d0/r2; sr2=sig2_LJ*ri2; sr6=sr2*sr2*sr2; sr12=sr6*sr6
        fm = 24.0d0*eps_LJ_*(2.0d0*sr12-sr6)*ri2
        Ep_lj = Ep_lj + 4.0d0*eps_LJ_*(sr12-sr6) - LJ_VSHFT

        F((i-1)*3+1)=F((i-1)*3+1)-fm*dx
        F((i-1)*3+2)=F((i-1)*3+2)-fm*dy
        F((i-1)*3+3)=F((i-1)*3+3)-fm*dz
        !$OMP ATOMIC
        F((j-1)*3+1)=F((j-1)*3+1)+fm*dx
        !$OMP ATOMIC
        F((j-1)*3+2)=F((j-1)*3+2)+fm*dy
        !$OMP ATOMIC
        F((j-1)*3+3)=F((j-1)*3+3)+fm*dz

        do aa = 1, 3
          if (aa == 1) then; da = dx
          else if (aa == 2) then; da = dy
          else; da = dz; end if
          do bb = 1, 3
            if (bb == 1) then; db = dx
            else if (bb == 2) then; db = dy
            else; db = dz; end if
            !$OMP ATOMIC
            vir9((aa-1)*3+bb) = vir9((aa-1)*3+bb) + da*fm*db
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine

  ! ═══ 運動エネルギー ═══
  double precision function ke_total(vel, mass, Na)
    double precision, intent(in) :: vel(:), mass(:)
    integer, intent(in) :: Na
    double precision :: s
    integer :: i, idx
    s = 0.0d0
    !$OMP PARALLEL DO PRIVATE(i,idx) REDUCTION(+:s)
    do i = 1, Na
      idx=(i-1)*3
      s = s + mass(i)*(vel(idx+1)**2 + vel(idx+2)**2 + vel(idx+3)**2)
    end do
    !$OMP END PARALLEL DO
    ke_total = 0.5d0*s/CONV
  end function

  ! ═══ 瞬間温度 ═══
  double precision function inst_T(KE, Nf)
    double precision, intent(in) :: KE
    integer, intent(in) :: Nf
    inst_T = 2.0d0*KE/(dble(Nf)*kB_)
  end function

  ! ═══ 瞬間圧力 ═══
  double precision function inst_P(W, KE, V)
    double precision, intent(in) :: W(9), KE, V
    inst_P = (2.0d0*KE + W(1) + W(5) + W(9))/(3.0d0*V)*eV2GPa
  end function

  ! ═══ NPT初期化 ═══
  subroutine make_npt(npt, T, Pe, Na)
    type(NPTState), intent(out) :: npt
    double precision, intent(in) :: T, Pe
    integer, intent(in) :: Na
    npt%Nf = 3*Na-3; npt%xi = 0.0d0
    npt%Q = max(dble(npt%Nf)*kB_*T*1.0d4, 1.0d-20)
    npt%Vg = 0.0d0
    npt%W_ = max(dble(npt%Nf+9)*kB_*T*1.0d6, 1.0d-20)
    npt%Pe = Pe; npt%Tt = T
  end subroutine

  ! ═══ clamp関数 ═══
  double precision function clamp_val(x, lo, hi)
    double precision, intent(in) :: x, lo, hi
    clamp_val = max(lo, min(x, hi))
  end function

  ! ═══ NPT 1ステップ ═══
  subroutine step_npt_airebo(pos, vel, F, vir9, h, hi, mass, Na, dt, npt, &
       nlc_r, nll_r, nlc_l, nll_l, Ep_rebo, Ep_lj, KE_out)
    double precision, intent(inout) :: pos(:), vel(:), F(:), vir9(9), h(9), hi(9)
    double precision, intent(in) :: mass(:), dt
    integer, intent(in) :: Na, nlc_r(:), nll_r(:), nlc_l(:), nll_l(:)
    type(NPTState), intent(inout) :: npt
    double precision, intent(out) :: Ep_rebo, Ep_lj, KE_out
    double precision :: hdt, V, KE, dP, eps_v, sc_nh, sc_pr, sc_v, mi
    double precision :: px, py, pz, vx, vy, vz
    double precision :: sx, sy, sz, vsx, vsy, vsz
    double precision :: eps2, sv2, V2
    integer :: i, a, idx

    hdt = 0.5d0*dt
    V = abs(mat_det9(h))
    KE = ke_total(vel, mass, Na)

    ! 半ステップ NH
    npt%xi = npt%xi + hdt*(2.0d0*KE - dble(npt%Nf)*kB_*npt%Tt)/npt%Q
    npt%xi = clamp_val(npt%xi, -0.05d0, 0.05d0)

    ! 半ステップ PR
    dP = inst_P(vir9, KE, V) - npt%Pe
    do a = 0, 2
      npt%Vg(a*4+1) = npt%Vg(a*4+1) + hdt*V*dP/(npt%W_*eV2GPa)
      npt%Vg(a*4+1) = clamp_val(npt%Vg(a*4+1), -0.005d0, 0.005d0)
    end do

    eps_v = npt%Vg(1)*hi(1) + npt%Vg(5)*hi(5) + npt%Vg(9)*hi(9)
    sc_nh = exp(-hdt*npt%xi)
    sc_pr = exp(-hdt*eps_v/3.0d0)
    sc_v = sc_nh*sc_pr

    ! 速度半ステップ
    !$OMP PARALLEL DO PRIVATE(i,idx,mi)
    do i = 1, Na
      idx = (i-1)*3; mi = CONV/mass(i)
      vel(idx+1) = vel(idx+1)*sc_v + hdt*F(idx+1)*mi
      vel(idx+2) = vel(idx+2)*sc_v + hdt*F(idx+2)*mi
      vel(idx+3) = vel(idx+3)*sc_v + hdt*F(idx+3)*mi
    end do
    !$OMP END PARALLEL DO

    ! 位置更新 (分率座標経由)
    !$OMP PARALLEL DO PRIVATE(i,idx,px,py,pz,vx,vy,vz,sx,sy,sz,vsx,vsy,vsz)
    do i = 1, Na
      idx = (i-1)*3
      px=pos(idx+1); py=pos(idx+2); pz=pos(idx+3)
      vx=vel(idx+1); vy=vel(idx+2); vz=vel(idx+3)
      sx=hi(1)*px+hi(2)*py+hi(3)*pz
      sy=hi(4)*px+hi(5)*py+hi(6)*pz
      sz=hi(7)*px+hi(8)*py+hi(9)*pz
      vsx=hi(1)*vx+hi(2)*vy+hi(3)*vz
      vsy=hi(4)*vx+hi(5)*vy+hi(6)*vz
      vsz=hi(7)*vx+hi(8)*vy+hi(9)*vz
      sx=sx+dt*vsx; sy=sy+dt*vsy; sz=sz+dt*vsz
      sx=sx-floor(sx); sy=sy-floor(sy); sz=sz-floor(sz)
      pos(idx+1)=sx; pos(idx+2)=sy; pos(idx+3)=sz
    end do
    !$OMP END PARALLEL DO

    ! セル行列更新
    do a = 0, 2
      h(a*3+1) = h(a*3+1) + dt*npt%Vg(a*3+1)
      h(a*3+2) = h(a*3+2) + dt*npt%Vg(a*3+2)
      h(a*3+3) = h(a*3+3) + dt*npt%Vg(a*3+3)
    end do
    call mat_inv9(h, hi)

    ! 分率→実空間
    !$OMP PARALLEL DO PRIVATE(i,idx,sx,sy,sz)
    do i = 1, Na
      idx = (i-1)*3
      sx=pos(idx+1); sy=pos(idx+2); sz=pos(idx+3)
      pos(idx+1)=h(1)*sx+h(2)*sy+h(3)*sz
      pos(idx+2)=h(4)*sx+h(5)*sy+h(6)*sz
      pos(idx+3)=h(7)*sx+h(8)*sy+h(9)*sz
    end do
    !$OMP END PARALLEL DO

    ! 力再計算
    F(1:Na*3) = 0.0d0; vir9 = 0.0d0
    call compute_rebo(F, vir9, pos, h, hi, nlc_r, nll_r, Na, Ep_rebo)
    call compute_lj(F, vir9, pos, h, hi, nlc_l, nll_l, Na, Ep_lj)

    ! 2回目の速度半ステップ
    eps2 = npt%Vg(1)*hi(1) + npt%Vg(5)*hi(5) + npt%Vg(9)*hi(9)
    sv2 = sc_nh*exp(-hdt*eps2/3.0d0)
    !$OMP PARALLEL DO PRIVATE(i,idx,mi)
    do i = 1, Na
      idx = (i-1)*3; mi = CONV/mass(i)
      vel(idx+1) = (vel(idx+1) + hdt*F(idx+1)*mi)*sv2
      vel(idx+2) = (vel(idx+2) + hdt*F(idx+2)*mi)*sv2
      vel(idx+3) = (vel(idx+3) + hdt*F(idx+3)*mi)*sv2
    end do
    !$OMP END PARALLEL DO

    ! 2回目のNH・PR
    KE = ke_total(vel, mass, Na); KE_out = KE
    npt%xi = npt%xi + hdt*(2.0d0*KE - dble(npt%Nf)*kB_*npt%Tt)/npt%Q
    npt%xi = clamp_val(npt%xi, -0.05d0, 0.05d0)
    V2 = abs(mat_det9(h))
    dP = inst_P(vir9, KE, V2) - npt%Pe
    do a = 0, 2
      npt%Vg(a*4+1) = npt%Vg(a*4+1) + hdt*V2*dP/(npt%W_*eV2GPa)
      npt%Vg(a*4+1) = clamp_val(npt%Vg(a*4+1), -0.005d0, 0.005d0)
    end do
  end subroutine

  ! ═══ OVITO XYZ出力 ═══
  subroutine write_ovito(iu, istep, dt, pos, vel, mol_id, h, Na)
    integer, intent(in) :: iu, istep, Na
    double precision, intent(in) :: dt, pos(:), vel(:), h(9)
    integer, intent(in) :: mol_id(:)
    integer :: i
    write(iu,'(I0)') Na
    write(iu,'(A,9(ES18.8),A,F10.4,A,I0,A)') 'Lattice="', &
      h(1),h(4),h(7),h(2),h(5),h(8),h(3),h(6),h(9), &
      '" Properties=species:S:1:pos:R:3:c_mol:I:1:vx:R:1:vy:R:1:vz:R:1 Time=', &
      istep*dt,' Step=',istep,' pbc="T T T"'
    do i = 1, Na
      write(iu,'(A,3(ES16.8),I6,3(ES16.8))') 'C ', &
        pos((i-1)*3+1), pos((i-1)*3+2), pos((i-1)*3+3), mol_id(i), &
        vel((i-1)*3+1), vel((i-1)*3+2), vel((i-1)*3+3)
    end do
  end subroutine

  ! ═══ リスタート保存 ═══
  subroutine write_restart_airebo(fname, istep, h, npt, pos, vel, &
       mol_id, mass, Na, Nmol, natom_mol)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: istep, Na, Nmol, natom_mol
    double precision, intent(in) :: h(9), pos(:), vel(:), mass(:)
    integer, intent(in) :: mol_id(:)
    type(NPTState), intent(in) :: npt
    integer :: i, ios
    open(unit=30, file=fname, status='replace', iostat=ios)
    if (ios /= 0) return
    write(30,'(A)') '# RESTART fuller_airebo_npt_md_serial_omp_acc'
    write(30,'(A,I0)') 'STEP ', istep
    write(30,'(A,I0)') 'NMOL ', Nmol
    write(30,'(A,I0)') 'NATOM_MOL ', natom_mol
    write(30,'(A,I0)') 'NATOM ', Na
    write(30,'(A,9(ES20.12))') 'H ', h
    write(30,'(A,5(ES20.12),I8)') 'NPT ', npt%xi, npt%Q, npt%W_, npt%Pe, npt%Tt, npt%Nf
    write(30,'(A,9(ES20.12))') 'VG ', npt%Vg
    do i = 1, Na
      write(30,'(A,I8,I8,7(ES20.12))') 'ATOM', i, mol_id(i), mass(i), &
        pos((i-1)*3+1), pos((i-1)*3+2), pos((i-1)*3+3), &
        vel((i-1)*3+1), vel((i-1)*3+2), vel((i-1)*3+3)
    end do
    write(30,'(A)') 'END'
    close(30)
  end subroutine

  ! ═══ リスタート読み込み ═══
  subroutine read_restart_airebo(fname, istep_out, h, npt, pos, vel, &
       mol_id, mass, Na_out, ok)
    character(len=*), intent(in) :: fname
    integer, intent(out) :: istep_out, Na_out
    double precision, intent(out) :: h(9)
    type(NPTState), intent(out) :: npt
    double precision, intent(inout) :: pos(:), vel(:), mass(:)
    integer, intent(inout) :: mol_id(:)
    logical, intent(out) :: ok
    character(len=1024) :: line
    character(len=32) :: tag
    integer :: ios, idx_a, mid_a, cnt
    double precision :: m_v, px, py, pz, vx, vy, vz
    ok = .false.; h = 0.0d0
    npt%xi=0; npt%Q=1; npt%W_=1; npt%Vg=0; npt%Pe=0; npt%Tt=298; npt%Nf=3
    istep_out = 0; Na_out = 0; cnt = 0
    open(unit=31, file=fname, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Cannot read: ', trim(fname); return
    end if
    do
      read(31,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line,*,iostat=ios) tag
      if (tag == 'STEP') then
        read(line(6:),*,iostat=ios) istep_out
      else if (tag == 'NATOM') then
        read(line(7:),*,iostat=ios) Na_out
      else if (tag == 'H') then
        read(line(3:),*,iostat=ios) h
      else if (tag == 'NPT') then
        read(line(5:),*,iostat=ios) npt%xi, npt%Q, npt%W_, npt%Pe, npt%Tt, npt%Nf
      else if (tag == 'VG') then
        read(line(4:),*,iostat=ios) npt%Vg
      else if (tag == 'ATOM') then
        cnt = cnt + 1
        read(line(5:),*,iostat=ios) idx_a, mid_a, m_v, px, py, pz, vx, vy, vz
        if (cnt <= size(pos)/3) then
          pos((cnt-1)*3+1)=px; pos((cnt-1)*3+2)=py; pos((cnt-1)*3+3)=pz
          vel((cnt-1)*3+1)=vx; vel((cnt-1)*3+2)=vy; vel((cnt-1)*3+3)=vz
          mol_id(cnt)=mid_a; mass(cnt)=m_v
        end if
      else if (tag == 'END') then
        exit
      end if
    end do
    close(31)
    if (cnt > 0) then
      Na_out = cnt; ok = .true.
      write(*,'(A,A,A,I0,A,I0,A)') '  Restart loaded: ', trim(fname), &
        ' (step ', istep_out, ', ', Na_out, ' atoms)'
    end if
  end subroutine

  ! ═══ ガウス乱数 ═══
  double precision function gauss_rand()
    double precision :: u1, u2
    call random_number(u1); call random_number(u2)
    if (u1 < 1.0d-30) u1 = 1.0d-30
    gauss_rand = sqrt(-2.0d0*log(u1))*cos(2.0d0*PI_*u2)
  end function

  ! ═══ コマンドライン引数パーサー (整数) ═══
  subroutine get_opt_int(key, defval, res)
    character(len=*), intent(in) :: key
    integer, intent(in) :: defval
    integer, intent(out) :: res
    integer :: i, n, kl, ios
    character(len=256) :: arg
    res = defval; n = command_argument_count(); kl = len_trim(key)
    do i = 1, n
      call get_command_argument(i, arg)
      if (arg(1:kl+3) == '--'//key(1:kl)//'=') read(arg(kl+4:),*,iostat=ios) res
    end do
  end subroutine

  ! ═══ コマンドライン引数パーサー (実数) ═══
  subroutine get_opt_dbl(key, defval, res)
    character(len=*), intent(in) :: key
    double precision, intent(in) :: defval
    double precision, intent(out) :: res
    integer :: i, n, kl, ios
    character(len=256) :: arg
    res = defval; n = command_argument_count(); kl = len_trim(key)
    do i = 1, n
      call get_command_argument(i, arg)
      if (arg(1:kl+3) == '--'//key(1:kl)//'=') read(arg(kl+4:),*,iostat=ios) res
    end do
  end subroutine

  ! ═══ コマンドライン引数パーサー (文字列) ═══
  subroutine get_opt_str(key, defval, res)
    character(len=*), intent(in) :: key, defval
    character(len=*), intent(out) :: res
    integer :: i, n, kl
    character(len=256) :: arg
    res = defval; n = command_argument_count(); kl = len_trim(key)
    do i = 1, n
      call get_command_argument(i, arg)
      if (arg(1:kl+3) == '--'//key(1:kl)//'=') res = arg(kl+4:)
    end do
  end subroutine

end module fuller_airebo_mod

! ═══════════════════════════════════════════════════════════════════════
!  メインプログラム
! ═══════════════════════════════════════════════════════════════════════
program fuller_airebo_npt_md
  use fuller_airebo_mod
  implicit none

  integer :: nc, Nmol, natom, Na, nsteps, nlup, Nmax, i, a, gstep, nav, j_idx
  integer :: coldstart, warmup, avg_from, avg_to, total_steps, gavg_from, gavg_to
  integer :: mon_int, prn, prn_pre, cur_prn, start_step
  integer :: nrec_o, nrec_rst
  double precision :: T, Pe, dt, a0, Rmol, Dmol, T_cold, T_init, sv, scale_v, tgt
  double precision :: init_scale
  double precision :: Ep_rebo, Ep_lj, KE, Tn, Pn, V_val, an_val
  double precision :: sT, sP, sa, sR, sL, sE, t_start, t_now, elapsed
  double precision :: h(9), hi(9), vir9(9), vcm(3)
  double precision :: mol_coords(MAX_NATOM,3)
  double precision, allocatable :: pos(:), vel(:), F(:), mass(:), mol_centers(:)
  integer, allocatable :: mol_id(:), nlc_r(:), nll_r(:), nlc_l(:), nll_l(:)
  type(NPTState) :: npt
  character(len=256) :: resfile, crystal_str, cryst_st
  character(len=4) :: phase_str
  character(len=256) :: rst_fname
  logical :: rst_ok, stop_requested
  integer :: Na_rst, istep_rst

  T_cold = 4.0d0

  ! ── コマンドライン引数の取得 ──
  call get_opt_int('cell', 3, nc)
  call get_opt_dbl('temp', 298.0d0, T)
  call get_opt_dbl('pres', 0.0d0, Pe)
  call get_opt_int('step', 10000, nsteps)
  call get_opt_dbl('dt', 0.5d0, dt)
  call get_opt_dbl('init_scale', 1.0d0, init_scale)
  call get_opt_int('coldstart', 0, coldstart)
  call get_opt_int('warmup', 0, warmup)
  call get_opt_int('from', 0, avg_from)
  call get_opt_int('to', 0, avg_to)
  call get_opt_int('ovito', 0, nrec_o)
  call get_opt_int('restart', 0, nrec_rst)
  call get_opt_int('mon', 0, mon_int)
  call get_opt_str('resfile', '', resfile)
  call get_opt_str('crystal', 'fcc', crystal_str)

  ! 結晶構造名を大文字化
  cryst_st = crystal_str
  do i = 1, len_trim(cryst_st)
    if (ichar(cryst_st(i:i)) >= ichar('a') .and. ichar(cryst_st(i:i)) <= ichar('z')) then
      cryst_st(i:i) = char(ichar(cryst_st(i:i)) - 32)
    end if
  end do

  if (avg_to <= 0) avg_to = nsteps
  if (avg_from <= 0) avg_from = max(1, nsteps-nsteps/4)
  total_steps = coldstart + warmup + nsteps
  gavg_from = coldstart + warmup + avg_from
  gavg_to = coldstart + warmup + avg_to
  start_step = 0; nlup = 20

  ! ── C60フォールバック生成 ──
  call generate_c60_airebo(mol_coords, natom, Rmol, Dmol)
  a0 = default_a0(Dmol, trim(cryst_st), init_scale)

  ! ── 最大分子数 ──
  if (trim(cryst_st) == 'FCC') then
    Nmax = 4*nc*nc*nc
  else if (trim(cryst_st) == 'HCP') then
    Nmax = 2*nc*nc*nc
  else
    Nmax = 2*nc*nc*nc
  end if

  ! ── 配列確保 ──
  allocate(mol_centers(Nmax*3))
  mol_centers = 0.0d0
  if (trim(cryst_st) == 'FCC') then
    Nmol = make_fcc(a0, nc, mol_centers, h)
  else if (trim(cryst_st) == 'HCP') then
    Nmol = make_hcp(a0, nc, mol_centers, h)
  else
    Nmol = make_bcc(a0, nc, mol_centers, h)
  end if
  Na = Nmol * natom

  allocate(pos(Na*3), vel(Na*3), F(Na*3), mass(Na), mol_id(Na))
  allocate(nlc_r(Na), nll_r(Na*MAX_REBO_NEIGH))
  allocate(nlc_l(Na), nll_l(Na*MAX_LJ_NEIGH))
  pos=0; vel=0; F=0; vir9=0; nlc_r=0; nll_r=0; nlc_l=0; nll_l=0

  ! ── 全原子位置・質量・分子IDの初期化 ──
  do i = 1, Nmol
    do a = 1, natom
      j_idx = (i-1)*natom + a
      pos((j_idx-1)*3+1) = mol_centers((i-1)*3+1) + mol_coords(a,1)
      pos((j_idx-1)*3+2) = mol_centers((i-1)*3+2) + mol_coords(a,2)
      pos((j_idx-1)*3+3) = mol_centers((i-1)*3+3) + mol_coords(a,3)
      mass(j_idx) = mC_; mol_id(j_idx) = i
    end do
  end do
  deallocate(mol_centers)

  ! ── バナー出力 ──
  write(*,'(A)') '========================================================================'
  write(*,'(A)') '  Fullerene Crystal NPT-MD — AIREBO (Serial, Fortran 95)'
  write(*,'(A)') '========================================================================'
  write(*,'(A,I0,A)') '  Fullerene       : C60 (', natom, ' atoms/mol)'
  write(*,'(A,A,1X,I0,A,I0,A,I0,A,I0,A,I0)') '  Crystal         : ', &
    trim(cryst_st), nc, 'x', nc, 'x', nc, '  Nmol=', Nmol, '  Natom=', Na
  write(*,'(A,F8.3,A,F6.1,A,F6.4,A,F5.2,A)') &
    '  a0=', a0, ' A  T=', T, ' K  P=', Pe, ' GPa  dt=', dt, ' fs'
  write(*,'(A,I0,A,I0,A,I0,A,I0)') &
    '  Production      : ', nsteps, ' steps  avg=', avg_from, '-', avg_to, &
    '  Total=', total_steps
  write(*,'(A)') '========================================================================'
  write(*,*)

  call mat_inv9(h, hi)

  ! ── 初期速度 ──
  T_init = T; if (coldstart > 0 .or. warmup > 0) T_init = T_cold
  call random_seed()
  do i = 1, Na
    sv = sqrt(kB_*T_init*CONV/mass(i))
    do a = 1, 3; vel((i-1)*3+a) = sv*gauss_rand(); end do
  end do
  vcm = 0.0d0
  do i = 1, Na
    vcm(1)=vcm(1)+vel((i-1)*3+1)
    vcm(2)=vcm(2)+vel((i-1)*3+2)
    vcm(3)=vcm(3)+vel((i-1)*3+3)
  end do
  vcm = vcm/dble(Na)
  do i = 1, Na
    vel((i-1)*3+1)=vel((i-1)*3+1)-vcm(1)
    vel((i-1)*3+2)=vel((i-1)*3+2)-vcm(2)
    vel((i-1)*3+3)=vel((i-1)*3+3)-vcm(3)
  end do

  ! ── NPT初期化 ──
  call make_npt(npt, T, Pe, Na); npt%Tt = T_init

  ! ── リスタートファイル読み込み ──
  if (len_trim(resfile) > 0) then
    call read_restart_airebo(trim(resfile), istep_rst, h, npt, &
         pos, vel, mol_id, mass, Na_rst, rst_ok)
    if (rst_ok) then
      start_step = istep_rst
      call mat_inv9(h, hi)
      Na = Na_rst
      write(*,'(A,I0)') '  Restarting from global step ', start_step
    end if
  end if

  ! ── 近傍リスト構築・初期力計算 ──
  call build_nlist_rebo(pos, h, hi, Na, mol_id, nlc_r, nll_r)
  call build_nlist_lj(pos, h, hi, Na, mol_id, nlc_l, nll_l)
  call apply_pbc(pos, h, hi, Na)
  F = 0.0d0; vir9 = 0.0d0
  call compute_rebo(F, vir9, pos, h, hi, nlc_r, nll_r, Na, Ep_rebo)
  call compute_lj(F, vir9, pos, h, hi, nlc_l, nll_l, Na, Ep_lj)

  ! ── 出力間隔設定 ──
  prn = mon_int; if (prn <= 0) prn = max(1, total_steps/50)
  prn_pre = prn; if (coldstart+warmup > 0) prn_pre = max(1, (coldstart+warmup)/100)
  sT=0; sP=0; sa=0; sR=0; sL=0; sE=0; nav=0
  call cpu_time(t_start)

  stop_requested = .false.

  ! ── OVITO出力ファイル ──
  if (nrec_o > 0) open(unit=40, file='ovito_traj_airebo_serial.xyz', status='replace')

  write(*,'(A8,A6,A8,A10,A9,A12,A12,A12,A8)') &
    'step', 'phase', 'T[K]', 'P[GPa]', 'a[A]', 'E_REBO', 'E_LJ', 'E_total', 't[s]'

  ! ═══ MDメインループ ═══
  do gstep = start_step+1, total_steps
    if (gstep <= coldstart) then
      npt%Tt = T_cold
    else if (gstep <= coldstart+warmup) then
      npt%Tt = T_cold + (T-T_cold)*dble(gstep-coldstart)/dble(max(warmup,1))
    else
      npt%Tt = T
    end if
    if (coldstart > 0 .and. gstep == coldstart+1) then
      npt%xi = 0.0d0; npt%Vg = 0.0d0
    end if
    if (gstep <= coldstart) npt%Vg = 0.0d0
    cur_prn = prn; if (gstep <= coldstart+warmup) cur_prn = prn_pre

    ! 近傍リスト再構築
    if (mod(gstep, nlup) == 0) then
      call mat_inv9(h, hi)
      call build_nlist_rebo(pos, h, hi, Na, mol_id, nlc_r, nll_r)
      call build_nlist_lj(pos, h, hi, Na, mol_id, nlc_l, nll_l)
    end if

    ! NPTステップ
    call step_npt_airebo(pos, vel, F, vir9, h, hi, mass, Na, dt, npt, &
         nlc_r, nll_r, nlc_l, nll_l, Ep_rebo, Ep_lj, KE)

    V_val = abs(mat_det9(h))
    Tn = inst_T(KE, npt%Nf)
    Pn = inst_P(vir9, KE, V_val)

    ! COLD/WARM 速度リスケーリング
    if ((gstep <= coldstart .or. gstep <= coldstart+warmup) .and. Tn > 0.1d0) then
      tgt = T_cold; if (gstep > coldstart) tgt = npt%Tt
      scale_v = sqrt(max(tgt, 0.1d0)/Tn)
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1, Na*3; vel(i) = vel(i)*scale_v; end do
      !$OMP END PARALLEL DO
      KE = ke_total(vel, mass, Na)
      Tn = inst_T(KE, npt%Nf)
      npt%xi = 0.0d0
      if (gstep <= coldstart) npt%Vg = 0.0d0
    end if

    an_val = h(1)/dble(nc)
    if (gstep >= gavg_from .and. gstep <= gavg_to) then
      sT=sT+Tn; sP=sP+Pn; sa=sa+an_val
      sR=sR+Ep_rebo/dble(Nmol); sL=sL+Ep_lj/dble(Nmol)
      sE=sE+(Ep_rebo+Ep_lj)/dble(Nmol); nav=nav+1
    end if

    ! OVITO出力
    if (nrec_o > 0 .and. mod(gstep, nrec_o) == 0) then
      call write_ovito(40, gstep, dt, pos, vel, mol_id, h, Na)
      flush(40)
    end if

    ! リスタート保存
    if (nrec_rst > 0 .and. (mod(gstep, nrec_rst) == 0 .or. gstep == total_steps)) then
      write(rst_fname,'(A,I0,A)') 'restart_airebo_serial_', gstep, '.rst'
      call write_restart_airebo(trim(rst_fname), gstep, h, npt, pos, vel, &
           mol_id, mass, Na, Nmol, natom)
      if (stop_requested) then
        write(*,'(A,I0,A)') &
          '  *** Stopped at restart checkpoint (global step ', gstep, ') ***'
        exit
      end if
    end if

    ! モニタリング出力
    if (mod(gstep, cur_prn) == 0 .or. gstep == total_steps) then
      call cpu_time(t_now); elapsed = t_now - t_start
      if (gstep <= coldstart) then
        phase_str = 'COLD'
      else if (gstep <= coldstart+warmup) then
        phase_str = 'WARM'
      else
        phase_str = 'PROD'
      end if
      write(*,'(I8,A6,F8.1,F10.3,F9.3,F12.4,F12.4,F12.4,F8.0)') &
        gstep, '  '//phase_str, Tn, Pn, an_val, &
        Ep_rebo/dble(Nmol), Ep_lj/dble(Nmol), &
        (Ep_rebo+Ep_lj)/dble(Nmol), elapsed
    end if
  end do

  if (nrec_o > 0) close(40)

  if (nav > 0) then
    write(*,*)
    write(*,'(A)') '========================================================================'
    write(*,'(A,I0,A,F7.2,A,F8.4,A,F8.4,A,F10.4,A,F10.4,A,F10.4)') &
      '  Avg(', nav, '): T=', sT/dble(nav), ' P=', sP/dble(nav), &
      ' a=', sa/dble(nav), ' REBO=', sR/dble(nav), &
      ' LJ=', sL/dble(nav), ' Total=', sE/dble(nav)
    write(*,'(A)') '========================================================================'
  end if
  call cpu_time(t_now)
  write(*,'(A,F8.1,A)') '  Done (', t_now-t_start, ' sec)'

  deallocate(pos, vel, F, mass, mol_id, nlc_r, nll_r, nlc_l, nll_l)
end program fuller_airebo_npt_md
