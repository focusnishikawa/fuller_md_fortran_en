! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2025, Takeshi Nishikawa
!===========================================================================
!  fuller_LJ_npt_md_core_serial_omp_acc.f90
!  C60フラーレン結晶 NPT分子動力学シミュレーション
!  (LJ剛体モデル・コア版 — Serial / OpenMP / OpenACC, Fortran 95)
!
!  コンパイル:
!    Serial (シングルスレッド):
!      gfortran -O3 -o fuller_LJ_core_serial &
!          fuller_LJ_npt_md_core_serial_omp_acc.f90
!
!    OpenMP (マルチスレッド):
!      gfortran -O3 -fopenmp -o fuller_LJ_core_omp &
!          fuller_LJ_npt_md_core_serial_omp_acc.f90
!
!    OpenACC GPU:
!      nvfortran -O3 -acc -gpu=cc80 -Minfo=accel -o fuller_LJ_core_gpu &
!          fuller_LJ_npt_md_core_serial_omp_acc.f90
!
!  実行:
!    ./fuller_LJ_core_serial [nc]
!    ./fuller_LJ_core_omp    [nc]
!    ./fuller_LJ_core_gpu    [nc]
!    ./fuller_LJ_core_serial --cell=5
!
!  固定パラメータ: T=300K, P=0GPa, dt=1fs, 1000ステップ
!  単位系: A (距離), amu (質量), eV (エネルギー), fs (時間), K (温度), GPa (圧力)
!
!  Note: !$OMP sentinels are ignored by non-OpenMP compilers (Fortran standard).
!        !$ACC directives are commented out for portability; uncomment for GPU.
!===========================================================================

module fuller_core_omp_mod
  implicit none

  ! ═══ 物理定数・単位変換 ═══
  double precision, parameter :: CONV       = 9.64853321d-3   ! eV*fs^2/(amu*A^2)
  double precision, parameter :: kB         = 8.617333262d-5  ! ボルツマン定数 [eV/K]
  double precision, parameter :: eV2GPa     = 160.21766208d0  ! eV/A^3 → GPa
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

  ! ═══ 分子パラメータ ═══
  integer, parameter :: C60_NATOM = 60
  integer, parameter :: MAX_NATOM = 84
  double precision, parameter :: MC60  = C60_NATOM * mC
  double precision, parameter :: RC60  = 3.55d0
  double precision, parameter :: RMCUT = RCUT + 2.0d0*RC60 + 1.0d0
  double precision, parameter :: RMCUT2 = RMCUT * RMCUT
  integer, parameter :: MAX_NEIGH = 80

  ! ═══ NPT状態 ═══
  type :: NPTState
    double precision :: xi, Q
    double precision :: Vg(9)
    double precision :: W, Pe, Tt
    integer :: Nf
  end type

contains

  ! ═══ H行列ユーティリティ ═══
  ! h(9) = {H00,H01,H02, H10,H11,H12, H20,H21,H22}
  ! アクセス: h(3*i + j + 1)  (i,j = 0,1,2)

  double precision function mat_det9(h)
    double precision, intent(in) :: h(9)
    mat_det9 = h(1)*(h(5)*h(9)-h(6)*h(8)) &
             - h(2)*(h(4)*h(9)-h(6)*h(7)) &
             + h(3)*(h(4)*h(8)-h(5)*h(7))
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

  subroutine qmul_flat(a, b, out)
    double precision, intent(in)  :: a(4), b(4)
    double precision, intent(out) :: out(4)
    out(1) = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4)
    out(2) = a(1)*b(2) + a(2)*b(1) + a(3)*b(4) - a(4)*b(3)
    out(3) = a(1)*b(3) - a(2)*b(4) + a(3)*b(1) + a(4)*b(2)
    out(4) = a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + a(4)*b(1)
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

  ! ═══ C60座標の生成 ═══
  subroutine generate_c60(coords, I0, Mmol, Rmol)
    double precision, intent(out) :: coords(60, 3)
    double precision, intent(out) :: I0, Mmol, Rmol
    double precision :: phi, tmp(60, 3), cm(3), r2, r, Isum
    integer :: n, p, s2, s3, s1
    integer :: cyc(3, 3)
    integer :: signs(2)

    Mmol = MC60
    phi = (1.0d0 + sqrt(5.0d0)) / 2.0d0
    signs(1) = -1; signs(2) = 1
    cyc(1,1) = 1; cyc(1,2) = 2; cyc(1,3) = 3
    cyc(2,1) = 2; cyc(2,2) = 3; cyc(2,3) = 1
    cyc(3,1) = 3; cyc(3,2) = 1; cyc(3,3) = 2

    tmp = 0.0d0
    n = 0

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

    Rmol = 0.0d0
    Isum = 0.0d0
    do n = 1, 60
      r2 = tmp(n,1)**2 + tmp(n,2)**2 + tmp(n,3)**2
      r = sqrt(r2)
      if (r > Rmol) Rmol = r
      Isum = Isum + mC * r2
    end do
    I0 = Isum * 2.0d0 / 3.0d0
    coords = tmp
  end subroutine

  ! ═══ FCC結晶生成 ═══
  function make_fcc(a, nc, pos, h) result(Nmol)
    double precision, intent(in)  :: a
    integer, intent(in)           :: nc
    double precision, intent(out) :: pos(:), h(9)
    integer :: Nmol
    double precision :: bas(4, 3)
    integer :: ix, iy, iz, b, idx

    bas(1,:) = (/ 0.0d0,    0.0d0,    0.0d0    /)
    bas(2,:) = (/ 0.5d0*a,  0.5d0*a,  0.0d0    /)
    bas(3,:) = (/ 0.5d0*a,  0.0d0,    0.5d0*a  /)
    bas(4,:) = (/ 0.0d0,    0.5d0*a,  0.5d0*a  /)

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
    h(1) = dble(nc) * a
    h(5) = dble(nc) * a
    h(9) = dble(nc) * a
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

    !!$ACC PARALLEL LOOP PRESENT(pos,h,hi)
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
    !!$ACC END PARALLEL LOOP
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
    double precision :: dmx, dmy, dmz, rax, ray, raz, rbx, rby, rbz
    double precision :: ddx, ddy, ddz, r2, ri2, sr2, sr6, sr12, fm
    double precision :: fx, fy, fz, Ep
    integer :: i, j, k, a, ai, bj, nni, ia, jb, idx

    ! ステップ1: lab座標計算
    !!$ACC PARALLEL LOOP PRESENT(qv,body,lab)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(R,a,bx,by,bz,idx)
    do i = 1, N
      call q2R_flat(qv((i-1)*4+1:(i-1)*4+4), R)
      do a = 1, natom
        bx = body(a, 1); by = body(a, 2); bz = body(a, 3)
        idx = (i-1)*natom*3 + (a-1)*3
        lab(idx+1) = R(1)*bx + R(2)*by + R(3)*bz
        lab(idx+2) = R(4)*bx + R(5)*by + R(6)*bz
        lab(idx+3) = R(7)*bx + R(8)*by + R(9)*bz
      end do
    end do
    !$OMP END PARALLEL DO
    !!$ACC END PARALLEL LOOP

    ! ゼロ初期化
    !!$ACC PARALLEL LOOP PRESENT(Fv,Tv)
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    do i = 1, N*3
      Fv(i) = 0.0d0
      Tv(i) = 0.0d0
    end do
    !$OMP END PARALLEL DO
    !!$ACC END PARALLEL LOOP

    !!$ACC PARALLEL LOOP PRESENT(Wm9)
    do i = 1, 9
      Wm9(i) = 0.0d0
    end do
    !!$ACC END PARALLEL LOOP
    Ep = 0.0d0

    ! ステップ2: LJ力計算メインカーネル
    !!$ACC PARALLEL LOOP PRESENT(pos,lab,Fv,Tv,Wm9,h,hi,nl_count,nl_list) &
    !!$ACC   REDUCTION(+:Ep)
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) REDUCTION(+:Ep) &
    !$OMP   PRIVATE(fi0,fi1,fi2,ti0,ti1,ti2,my_Ep, &
    !$OMP           w00,w01,w02,w10,w11,w12,w20,w21,w22, &
    !$OMP           nni,k,j,dmx,dmy,dmz,ai,ia,rax,ray,raz, &
    !$OMP           bj,jb,rbx,rby,rbz,ddx,ddy,ddz, &
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
            rbx = lab(jb+1); rby = lab(jb+2); rbz = lab(jb+3)
            ddx = dmx + rbx - rax
            ddy = dmy + rby - ray
            ddz = dmz + rbz - raz
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
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(1) = Wm9(1) + w00
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(2) = Wm9(2) + w01
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(3) = Wm9(3) + w02
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(4) = Wm9(4) + w10
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(5) = Wm9(5) + w11
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(6) = Wm9(6) + w12
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(7) = Wm9(7) + w20
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(8) = Wm9(8) + w21
      !!$ACC ATOMIC UPDATE
      !$OMP ATOMIC
      Wm9(9) = Wm9(9) + w22

      Ep = Ep + my_Ep
    end do
    !$OMP END PARALLEL DO
    !!$ACC END PARALLEL LOOP

    calc_forces = Ep
  end function

  ! ═══ 運動エネルギー ═══
  double precision function ke_trans(vel, N, Mmol)
    double precision, intent(in) :: vel(:), Mmol
    integer, intent(in) :: N
    double precision :: s
    integer :: i, idx
    s = 0.0d0
    !!$ACC PARALLEL LOOP PRESENT(vel) REDUCTION(+:s)
    !$OMP PARALLEL DO REDUCTION(+:s) PRIVATE(idx)
    do i = 1, N
      idx = (i-1)*3
      s = s + vel(idx+1)**2 + vel(idx+2)**2 + vel(idx+3)**2
    end do
    !$OMP END PARALLEL DO
    !!$ACC END PARALLEL LOOP
    ke_trans = 0.5d0 * Mmol * s / CONV
  end function

  double precision function ke_rot(omg, N, I0)
    double precision, intent(in) :: omg(:), I0
    integer, intent(in) :: N
    double precision :: s
    integer :: i, idx
    s = 0.0d0
    !!$ACC PARALLEL LOOP PRESENT(omg) REDUCTION(+:s)
    !$OMP PARALLEL DO REDUCTION(+:s) PRIVATE(idx)
    do i = 1, N
      idx = (i-1)*3
      s = s + omg(idx+1)**2 + omg(idx+2)**2 + omg(idx+3)**2
    end do
    !$OMP END PARALLEL DO
    !!$ACC END PARALLEL LOOP
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
    !!$ACC PARALLEL LOOP PRESENT(vel,omg,Fv,Tv)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,a)
    do i = 1, N
      idx = (i-1)*3
      do a = 1, 3
        vel(idx+a) = vel(idx+a)*sc_v + hdt*Fv(idx+a)*cF
        omg(idx+a) = omg(idx+a)*sc_nh + hdt*Tv(idx+a)*cT
      end do
    end do
    !$OMP END PARALLEL DO
    !!$ACC END PARALLEL LOOP

    ! (D) 座標更新 (分率座標)
    !!$ACC PARALLEL LOOP PRESENT(pos,vel,hi)
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
    !!$ACC END PARALLEL LOOP

    ! (E) セルH行列の更新
    do a = 0, 2
      h(a*3+1) = h(a*3+1) + dt*npt%Vg(a*3+1)
      h(a*3+2) = h(a*3+2) + dt*npt%Vg(a*3+2)
      h(a*3+3) = h(a*3+3) + dt*npt%Vg(a*3+3)
    end do

    ! (F) 分率->実座標
    !!$ACC PARALLEL LOOP PRESENT(pos,h)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,sx,sy,sz)
    do i = 1, N
      idx = (i-1)*3
      sx = pos(idx+1); sy = pos(idx+2); sz = pos(idx+3)
      pos(idx+1) = h(1)*sx + h(2)*sy + h(3)*sz
      pos(idx+2) = h(4)*sx + h(5)*sy + h(6)*sz
      pos(idx+3) = h(7)*sx + h(8)*sy + h(9)*sz
    end do
    !$OMP END PARALLEL DO
    !!$ACC END PARALLEL LOOP

    ! (G) 四元数更新
    !!$ACC PARALLEL LOOP PRESENT(qv,omg)
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
    !!$ACC END PARALLEL LOOP

    ! (H) 力の再計算
    call mat_inv9(h, hi)
    Ep_out = calc_forces(Fv, Tv, Wm9, pos, qv, body, h, hi, &
                         nl_count, nl_list, N, natom, rmcut2, lab)

    ! (I) 速度の後半更新
    eps_tr2 = npt%Vg(1)*hi(1) + npt%Vg(5)*hi(5) + npt%Vg(9)*hi(9)
    sc_v2 = sc_nh * exp(-hdt*eps_tr2/3.0d0)
    !!$ACC PARALLEL LOOP PRESENT(vel,omg,Fv,Tv)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(idx,a)
    do i = 1, N
      idx = (i-1)*3
      do a = 1, 3
        vel(idx+a) = (vel(idx+a) + hdt*Fv(idx+a)*cF) * sc_v2
        omg(idx+a) = (omg(idx+a) + hdt*Tv(idx+a)*cT) * sc_nh
      end do
    end do
    !$OMP END PARALLEL DO
    !!$ACC END PARALLEL LOOP

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

  ! ═══ ガウス乱数 (Box-Muller法) ═══
  double precision function gauss_rand()
    double precision :: u1, u2
    call random_number(u1)
    call random_number(u2)
    if (u1 < 1.0d-30) u1 = 1.0d-30
    gauss_rand = sqrt(-2.0d0*log(u1)) * cos(2.0d0*PI_*u2)
  end function

end module fuller_core_omp_mod


! ═══════════════════════════════════════════════════════════
!  メインプログラム
! ═══════════════════════════════════════════════════════════
program fuller_LJ_core_omp_acc
  use fuller_core_omp_mod
!$ use omp_lib
  implicit none

  integer :: nc, N, natom, nsteps, mon_int, nlup, avg_from, nav
  double precision :: T_target, Pe, dt, a0
  double precision :: c60_coords(60, 3), I0, Mmol, Rmol

  double precision, allocatable :: pos(:), vel(:), omg(:), qv(:)
  double precision, allocatable :: Fv(:), Tv(:), lab(:)
  double precision :: h(9), hi(9), Wm9(9)
  double precision, allocatable :: body(:,:)
  integer, allocatable :: nl_count(:), nl_list(:)

  double precision :: sv, sw, n_norm, vcm(3)
  double precision :: Ep, KE, kt_val, V_val, Tn, Pn, Ec, an_val
  double precision :: sT, sP, sa, sEp
  double precision :: mem_lab, mem_nl, mem_total
  double precision :: t_start, t_now, elapsed
  integer :: i, a, g, Nmax, idx
  integer :: omp_threads
  character(len=256) :: arg_str
  integer :: nargs

  ! -- パラメータ設定 --
  nc = 3
  nargs = command_argument_count()
  do i = 1, nargs
    call get_command_argument(i, arg_str)
    if (arg_str(1:7) == '--cell=') then
      read(arg_str(8:), *) nc
    else if (arg_str(1:1) /= '-') then
      read(arg_str, *) nc
    end if
  end do
  if (nc < 1 .or. nc > 8) then
    write(*,*) 'Error: nc must be 1-8 (got ', nc, ')'
    stop 1
  end if

  nsteps = 1000
  mon_int = 100
  nlup = 25
  T_target = 300.0d0
  Pe = 0.0d0
  dt = 1.0d0
  a0 = 14.17d0
  avg_from = nsteps - nsteps/4

  ! -- C60分子座標の生成 --
  call generate_c60(c60_coords, I0, Mmol, Rmol)
  natom = C60_NATOM

  ! -- 配列確保 --
  Nmax = 4 * nc * nc * nc
  allocate(pos(Nmax*3), vel(Nmax*3), omg(Nmax*3), qv(Nmax*4))
  allocate(Fv(Nmax*3), Tv(Nmax*3), lab(Nmax*natom*3))
  allocate(body(natom, 3))
  allocate(nl_count(Nmax), nl_list(Nmax*MAX_NEIGH))
  pos = 0.0d0; vel = 0.0d0; omg = 0.0d0; qv = 0.0d0
  Fv = 0.0d0; Tv = 0.0d0; lab = 0.0d0
  h = 0.0d0; hi = 0.0d0; Wm9 = 0.0d0
  nl_count = 0; nl_list = 0

  body = c60_coords

  N = make_fcc(a0, nc, pos, h)
  call mat_inv9(h, hi)

  ! -- OpenMPスレッド数取得 --
  omp_threads = 1
!$ omp_threads = omp_get_max_threads()

  ! -- バナー表示 --
  write(*,'(A)') '================================================================'
  if (omp_threads > 1) then
    write(*,'(A,I0,A)') &
         '  C60 LJ NPT-MD Core (Serial / OpenMP / OpenACC, ', &
         omp_threads, ' threads)'
  else
    write(*,'(A)') '  C60 LJ NPT-MD Core (Serial / OpenMP / OpenACC)'
  end if
  write(*,'(A)') '================================================================'
  write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '  FCC cell        : ', &
       nc, 'x', nc, 'x', nc, '  N=', N, ' molecules'
  write(*,'(A,I0)') '  Atoms/molecule  : ', natom
  write(*,'(A,F5.2,A,F5.0,A,F4.1,A,F4.1,A,I0)') &
       '  a0=', a0, ' A  T=', T_target, ' K  P=', Pe, &
       ' GPa  dt=', dt, ' fs  steps=', nsteps

  mem_lab = dble(N)*dble(natom)*3.0d0*8.0d0/1024.0d0/1024.0d0
  mem_nl  = dble(N)*dble(MAX_NEIGH)*4.0d0/1024.0d0/1024.0d0
  mem_total = dble(N*3*8*5 + N*4*8 + N*natom*3*8 + natom*3*8 &
              + N*4 + N*MAX_NEIGH*4 + 9*8*3) / 1024.0d0/1024.0d0
  write(*,'(A,F6.2,A,F6.2,A,F6.2,A)') &
       '  Memory          : lab=', mem_lab, ' MB  nl=', mem_nl, &
       ' MB  total=', mem_total, ' MB'
  write(*,'(A)') '================================================================'
  write(*,*)

  ! -- 初期速度 (Maxwell-Boltzmann分布) --
  call random_seed()   ! システム依存のシード
  sv = sqrt(kB*T_target*CONV/Mmol)
  sw = sqrt(kB*T_target*CONV/I0)
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

  ! -- NPT初期化 --
  block
    type(NPTState) :: npt_state
    call make_npt(npt_state, T_target, Pe, N)

    ! -- 初回近傍リスト + PBC + 力計算 --
    call nlist_build_sym(pos, h, hi, N, RMCUT, nl_count, nl_list)
    call apply_pbc(pos, h, hi, N)
    Ep = calc_forces(Fv, Tv, Wm9, pos, qv, body, h, hi, &
                     nl_count, nl_list, N, natom, RMCUT2, lab)

    sT = 0.0d0; sP = 0.0d0; sa = 0.0d0; sEp = 0.0d0; nav = 0
    call cpu_time(t_start)

    write(*,'(A8,A8,A10,A9,A11,A8)') &
         'step', 'T[K]', 'P[GPa]', 'a[A]', 'Ecoh[eV]', 't[s]'

    ! ═══ MDメインループ ═══
    do g = 1, nsteps

      ! 近傍リスト再構築
      if (mod(g, nlup) == 0) then
        call mat_inv9(h, hi)
        call nlist_build_sym(pos, h, hi, N, RMCUT, nl_count, nl_list)
      end if

      ! 1ステップ
      call step_npt(pos, vel, qv, omg, Fv, Tv, Wm9, h, hi, &
           body, I0, Mmol, N, natom, RMCUT2, dt, npt_state, &
           nl_count, nl_list, lab, Ep, KE)

      ! 瞬時物理量
      kt_val = ke_trans(vel, N, Mmol)
      V_val = abs(mat_det9(h))
      Tn = inst_T(KE, npt_state%Nf)
      Pn = inst_P(Wm9, kt_val, V_val)
      Ec = Ep / dble(N)
      an_val = h(1) / dble(nc)
      if (g >= avg_from) then
        sT = sT + Tn; sP = sP + Pn; sa = sa + an_val; sEp = sEp + Ec
        nav = nav + 1
      end if

      ! モニタリング出力
      if (mod(g, mon_int) == 0 .or. g == nsteps) then
        call cpu_time(t_now)
        elapsed = t_now - t_start
        write(*,'(I8,F8.1,F10.3,F9.3,F11.5,F8.0)') &
             g, Tn, Pn, an_val, Ec, elapsed
      end if
    end do

    if (nav > 0) then
      write(*,'(A,I0,A,F7.2,A,F8.4,A,F8.4,A,F10.5)') &
           'Avg(', nav, '): T=', sT/dble(nav), ' P=', sP/dble(nav), &
           ' a=', sa/dble(nav), ' Ecoh=', sEp/dble(nav)
    end if
    call cpu_time(t_now)
    write(*,'(A,F6.1,A)') 'Done ', t_now - t_start, 's'
  end block

  ! -- 解放 --
  deallocate(pos, vel, omg, qv, Fv, Tv, lab, body, nl_count, nl_list)

end program fuller_LJ_core_omp_acc
