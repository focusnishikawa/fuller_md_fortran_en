! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2025, Takeshi Nishikawa
!===========================================================================
!  fuller_LJ_npt_mmmd_serial_omp_acc.f90
!  Fullerene Crystal NPT-MD
!  (Molecular Mechanics Force Field + LJ Intermolecular, Serial / OpenMP, Fortran 95)
!  V_total = V_bond + V_angle + V_dihedral + V_improper + V_LJ + V_Coulomb
!
!  コンパイル:
!    gfortran -O3 -o fuller_LJ_npt_mmmd_serial fuller_LJ_npt_mmmd_serial_omp_acc.f90
!    gfortran -O3 -fopenmp -o fuller_LJ_npt_mmmd_omp fuller_LJ_npt_mmmd_serial_omp_acc.f90
!  単位系: A, amu, eV, fs, K, GPa
!===========================================================================
module fuller_mmmd_mod
  implicit none
  double precision, parameter :: CONV=9.64853321d-3, kB=8.617333262d-5
  double precision, parameter :: eV2GPa=160.21766208d0, eV2kcalmol=23.06054783d0
  double precision, parameter :: kcal2eV=1.0d0/eV2kcalmol, PI_=3.14159265358979323846d0
  double precision, parameter :: mC=12.011d0
  double precision, parameter :: sigma_LJ=3.431d0, eps_LJ=2.635d-3
  double precision, parameter :: RCUT=3.0d0*sigma_LJ, RCUT2=RCUT*RCUT
  double precision, parameter :: sig2_LJ=sigma_LJ*sigma_LJ
  double precision, parameter :: sr_v=1.0d0/3.0d0, sr2_v=sr_v*sr_v
  double precision, parameter :: sr6_v=sr2_v*sr2_v*sr2_v
  double precision, parameter :: VSHFT=4.0d0*eps_LJ*(sr6_v*sr6_v-sr6_v)
  integer, parameter :: MAX_LJ_NEIGH=400
  type :: NPTState
    double precision :: xi,Q,Vg(9),W_,Pe,Tt
    integer :: Nf
  end type
contains
  double precision function mat_det9(h)
    double precision, intent(in) :: h(9)
    mat_det9=h(1)*(h(5)*h(9)-h(6)*h(8))-h(2)*(h(4)*h(9)-h(6)*h(7))+h(3)*(h(4)*h(8)-h(5)*h(7))
  end function
  subroutine mat_inv9(h,hi)
    double precision, intent(in) :: h(9); double precision, intent(out) :: hi(9)
    double precision :: d,id; d=mat_det9(h); id=1.0d0/d
    hi(1)=id*(h(5)*h(9)-h(6)*h(8));hi(2)=id*(h(3)*h(8)-h(2)*h(9));hi(3)=id*(h(2)*h(6)-h(3)*h(5))
    hi(4)=id*(h(6)*h(7)-h(4)*h(9));hi(5)=id*(h(1)*h(9)-h(3)*h(7));hi(6)=id*(h(3)*h(4)-h(1)*h(6))
    hi(7)=id*(h(4)*h(8)-h(5)*h(7));hi(8)=id*(h(2)*h(7)-h(1)*h(8));hi(9)=id*(h(1)*h(5)-h(2)*h(4))
  end subroutine
  subroutine mimg9(dx,dy,dz,hi,h)
    double precision, intent(inout) :: dx,dy,dz
    double precision, intent(in) :: hi(9),h(9)
    double precision :: s0,s1,s2
    s0=hi(1)*dx+hi(2)*dy+hi(3)*dz;s1=hi(4)*dx+hi(5)*dy+hi(6)*dz;s2=hi(7)*dx+hi(8)*dy+hi(9)*dz
    s0=s0-anint(s0);s1=s1-anint(s1);s2=s2-anint(s2)
    dx=h(1)*s0+h(2)*s1+h(3)*s2;dy=h(4)*s0+h(5)*s1+h(6)*s2;dz=h(7)*s0+h(8)*s1+h(9)*s2
  end subroutine
  subroutine apply_pbc(pos,h,hi,N)
    double precision, intent(inout) :: pos(:); double precision, intent(in) :: h(9),hi(9)
    integer, intent(in) :: N; double precision :: px,py,pz,s0,s1,s2; integer :: i,idx
    !$OMP PARALLEL DO PRIVATE(i,idx,px,py,pz,s0,s1,s2)
    do i=1,N;idx=(i-1)*3;px=pos(idx+1);py=pos(idx+2);pz=pos(idx+3)
      s0=hi(1)*px+hi(2)*py+hi(3)*pz;s1=hi(4)*px+hi(5)*py+hi(6)*pz;s2=hi(7)*px+hi(8)*py+hi(9)*pz
      s0=s0-floor(s0);s1=s1-floor(s1);s2=s2-floor(s2)
      pos(idx+1)=h(1)*s0+h(2)*s1+h(3)*s2;pos(idx+2)=h(4)*s0+h(5)*s1+h(6)*s2
      pos(idx+3)=h(7)*s0+h(8)*s1+h(9)*s2
    end do
    !$OMP END PARALLEL DO
  end subroutine
  function make_fcc(a,nc,pos,h) result(Nmol)
    double precision, intent(in) :: a; integer, intent(in) :: nc
    double precision, intent(out) :: pos(:),h(9); integer :: Nmol,ix,iy,iz,b,idx
    double precision :: bas(4,3)
    bas(1,:)=(/0d0,0d0,0d0/);bas(2,:)=(/.5d0*a,.5d0*a,0d0/)
    bas(3,:)=(/.5d0*a,0d0,.5d0*a/);bas(4,:)=(/0d0,.5d0*a,.5d0*a/); Nmol=0
    do ix=0,nc-1;do iy=0,nc-1;do iz=0,nc-1;do b=1,4;Nmol=Nmol+1;idx=(Nmol-1)*3
      pos(idx+1)=a*dble(ix)+bas(b,1);pos(idx+2)=a*dble(iy)+bas(b,2);pos(idx+3)=a*dble(iz)+bas(b,3)
    end do;end do;end do;end do
    h=0;h(1)=dble(nc)*a;h(5)=dble(nc)*a;h(9)=dble(nc)*a
  end function
  ! ═══ C60フォールバック生成 ═══
  subroutine generate_c60_mmmd(coords,natom,Rmol,Dmol,Nb_mol,bonds)
    double precision, intent(out) :: coords(84,3); integer, intent(out) :: natom
    double precision, intent(out) :: Rmol, Dmol
    integer, intent(out) :: Nb_mol, bonds(200,2)
    double precision :: phi,tmp(60,3),cm(3),r2,r,dx,dy,dz,d
    integer :: n,p,s1,s2,s3,cyc(3,3),signs(2),i,j
    natom=60; phi=(1+sqrt(5.0d0))/2; signs(1)=-1;signs(2)=1
    cyc(1,1)=1;cyc(1,2)=2;cyc(1,3)=3;cyc(2,1)=2;cyc(2,2)=3;cyc(2,3)=1
    cyc(3,1)=3;cyc(3,2)=1;cyc(3,3)=2; tmp=0; n=0
    do p=1,3;do s2=1,2;do s3=1,2;n=n+1;tmp(n,cyc(p,1))=0
      tmp(n,cyc(p,2))=dble(signs(s2));tmp(n,cyc(p,3))=dble(signs(s3))*3*phi
    end do;end do;end do
    do p=1,3;do s1=1,2;do s2=1,2;do s3=1,2;n=n+1
      tmp(n,cyc(p,1))=dble(signs(s1))*2;tmp(n,cyc(p,2))=dble(signs(s2))*(1+2*phi)
      tmp(n,cyc(p,3))=dble(signs(s3))*phi
    end do;end do;end do;end do
    do p=1,3;do s1=1,2;do s2=1,2;do s3=1,2;n=n+1
      tmp(n,cyc(p,1))=dble(signs(s1));tmp(n,cyc(p,2))=dble(signs(s2))*(2+phi)
      tmp(n,cyc(p,3))=dble(signs(s3))*2*phi
    end do;end do;end do;end do
    cm=0; do n=1,60; cm=cm+tmp(n,:); end do; cm=cm/60
    do n=1,60; tmp(n,:)=(tmp(n,:)-cm)*0.72d0; end do
    Rmol=0;Dmol=0
    do i=1,60;r2=tmp(i,1)**2+tmp(i,2)**2+tmp(i,3)**2;r=sqrt(r2);if(r>Rmol)Rmol=r
      do j=i+1,60;dx=tmp(i,1)-tmp(j,1);dy=tmp(i,2)-tmp(j,2);dz=tmp(i,3)-tmp(j,3)
        d=sqrt(dx*dx+dy*dy+dz*dz);if(d>Dmol)Dmol=d;end do;end do
    coords(1:60,:)=tmp
    ! C60のボンド: 距離1.45以下のペア
    Nb_mol=0
    do i=1,60; do j=i+1,60
      dx=tmp(i,1)-tmp(j,1);dy=tmp(i,2)-tmp(j,2);dz=tmp(i,3)-tmp(j,3)
      d=sqrt(dx*dx+dy*dy+dz*dz)
      if(d<1.6d0) then; Nb_mol=Nb_mol+1; bonds(Nb_mol,1)=i; bonds(Nb_mol,2)=j; end if
    end do; end do
  end subroutine
  ! ═══ LJ近傍リスト (half-list, 分子間のみ) ═══
  subroutine build_nlist_lj(pos,h,hi,Na,mol_id,nlc,nll)
    double precision, intent(in) :: pos(:),h(9),hi(9)
    integer, intent(in) :: Na, mol_id(:)
    integer, intent(out) :: nlc(:), nll(:)
    double precision :: rc2,dx,dy,dz,r2; integer :: i,j
    rc2=(RCUT+2)**2; nlc(1:Na)=0
    do i=1,Na-1; do j=i+1,Na
      if(mol_id(j)==mol_id(i)) cycle
      dx=pos((j-1)*3+1)-pos((i-1)*3+1);dy=pos((j-1)*3+2)-pos((i-1)*3+2)
      dz=pos((j-1)*3+3)-pos((i-1)*3+3); call mimg9(dx,dy,dz,hi,h)
      r2=dx*dx+dy*dy+dz*dz
      if(r2<rc2.and.nlc(i)<MAX_LJ_NEIGH) then
        nlc(i)=nlc(i)+1; nll((i-1)*MAX_LJ_NEIGH+nlc(i))=j
      end if
    end do; end do
  end subroutine
  ! ═══ 力計算 (Bond+Angle+LJ) ═══
  subroutine compute_forces(F,vir9,pos,h,hi,Na, &
       Nb,b_i0,b_i1,b_kb,b_r0, &
       nlc,nll,mol_id, Eb_out,Elj_out,Etot_out)
    double precision, intent(out) :: F(:),vir9(9)
    double precision, intent(in) :: pos(:),h(9),hi(9)
    integer, intent(in) :: Na,Nb,b_i0(:),b_i1(:),nlc(:),nll(:),mol_id(:)
    double precision, intent(in) :: b_kb(:),b_r0(:)
    double precision, intent(out) :: Eb_out,Elj_out,Etot_out
    double precision :: dx,dy,dz,r,dr,fm,fx,fy,fz,r2,ri2,sr2,sr6,sr12,Eb,Elj
    integer :: b,i,j,jn,nni
    F(1:Na*3)=0; vir9=0; Eb=0; Elj=0
    ! Bond stretching
    !$OMP PARALLEL DO PRIVATE(b,i,j,dx,dy,dz,r,dr,fm,fx,fy,fz) REDUCTION(+:Eb)
    do b=1,Nb
      i=b_i0(b); j=b_i1(b)
      dx=pos((j-1)*3+1)-pos((i-1)*3+1);dy=pos((j-1)*3+2)-pos((i-1)*3+2)
      dz=pos((j-1)*3+3)-pos((i-1)*3+3); call mimg9(dx,dy,dz,hi,h)
      r=sqrt(dx*dx+dy*dy+dz*dz); if(r<1d-10) cycle; dr=r-b_r0(b)
      Eb=Eb+0.5d0*b_kb(b)*dr*dr; fm=-b_kb(b)*dr/r; fx=fm*dx;fy=fm*dy;fz=fm*dz
      !$OMP ATOMIC
      F((i-1)*3+1)=F((i-1)*3+1)+fx
      !$OMP ATOMIC
      F((i-1)*3+2)=F((i-1)*3+2)+fy
      !$OMP ATOMIC
      F((i-1)*3+3)=F((i-1)*3+3)+fz
      !$OMP ATOMIC
      F((j-1)*3+1)=F((j-1)*3+1)-fx
      !$OMP ATOMIC
      F((j-1)*3+2)=F((j-1)*3+2)-fy
      !$OMP ATOMIC
      F((j-1)*3+3)=F((j-1)*3+3)-fz
      !$OMP ATOMIC
      vir9(1)=vir9(1)+dx*fx
      !$OMP ATOMIC
      vir9(5)=vir9(5)+dy*fy
      !$OMP ATOMIC
      vir9(9)=vir9(9)+dz*fz
    end do
    !$OMP END PARALLEL DO
    ! LJ intermolecular (half-list)
    !$OMP PARALLEL DO PRIVATE(i,nni,jn,j,dx,dy,dz,r2,ri2,sr2,sr6,sr12,fm) &
    !$OMP SCHEDULE(DYNAMIC,4) REDUCTION(+:Elj)
    do i=1,Na; nni=nlc(i)
      do jn=1,nni; j=nll((i-1)*MAX_LJ_NEIGH+jn)
        dx=pos((j-1)*3+1)-pos((i-1)*3+1);dy=pos((j-1)*3+2)-pos((i-1)*3+2)
        dz=pos((j-1)*3+3)-pos((i-1)*3+3); call mimg9(dx,dy,dz,hi,h)
        r2=dx*dx+dy*dy+dz*dz; if(r2>RCUT2) cycle; if(r2<0.25d0) r2=0.25d0
        ri2=1/r2;sr2=sig2_LJ*ri2;sr6=sr2*sr2*sr2;sr12=sr6*sr6
        fm=24*eps_LJ*(2*sr12-sr6)*ri2; Elj=Elj+4*eps_LJ*(sr12-sr6)-VSHFT
        !$OMP ATOMIC
        F((i-1)*3+1)=F((i-1)*3+1)-fm*dx
        !$OMP ATOMIC
        F((i-1)*3+2)=F((i-1)*3+2)-fm*dy
        !$OMP ATOMIC
        F((i-1)*3+3)=F((i-1)*3+3)-fm*dz
        !$OMP ATOMIC
        F((j-1)*3+1)=F((j-1)*3+1)+fm*dx
        !$OMP ATOMIC
        F((j-1)*3+2)=F((j-1)*3+2)+fm*dy
        !$OMP ATOMIC
        F((j-1)*3+3)=F((j-1)*3+3)+fm*dz
        !$OMP ATOMIC
        vir9(1)=vir9(1)+dx*fm*dx
        !$OMP ATOMIC
        vir9(5)=vir9(5)+dy*fm*dy
        !$OMP ATOMIC
        vir9(9)=vir9(9)+dz*fm*dz
      end do
    end do
    !$OMP END PARALLEL DO
    Eb_out=Eb; Elj_out=Elj; Etot_out=Eb+Elj
  end subroutine
  double precision function ke_total(vel,Na)
    double precision, intent(in) :: vel(:); integer, intent(in) :: Na
    double precision :: s; integer :: i,idx; s=0
    !$OMP PARALLEL DO PRIVATE(i,idx) REDUCTION(+:s)
    do i=1,Na;idx=(i-1)*3;s=s+mC*(vel(idx+1)**2+vel(idx+2)**2+vel(idx+3)**2);end do
    !$OMP END PARALLEL DO
    ke_total=0.5d0*s/CONV
  end function
  double precision function inst_T(KE,Nf)
    double precision, intent(in) :: KE; integer, intent(in) :: Nf
    inst_T=2*KE/(dble(Nf)*kB)
  end function
  double precision function inst_P(W,KE,V)
    double precision, intent(in) :: W(9),KE,V
    inst_P=(2*KE+W(1)+W(5)+W(9))/(3*V)*eV2GPa
  end function
  subroutine make_npt(npt,T,Pe,Na)
    type(NPTState),intent(out)::npt;double precision,intent(in)::T,Pe;integer,intent(in)::Na
    npt%Nf=3*Na-3;npt%xi=0;npt%Q=max(dble(npt%Nf)*kB*T*1d4,1d-20);npt%Vg=0
    npt%W_=max(dble(npt%Nf+9)*kB*T*1d6,1d-20);npt%Pe=Pe;npt%Tt=T
  end subroutine
  double precision function clamp_val(x,lo,hi)
    double precision, intent(in) :: x,lo,hi; clamp_val=max(lo,min(x,hi))
  end function
  subroutine step_npt_mm(pos,vel,F,vir9,h,hi,Na,dt,npt, &
       Nb,b_i0,b_i1,b_kb,b_r0,nlc,nll,mol_id,Eb_out,Elj_out,Etot_out,KE_out)
    double precision, intent(inout) :: pos(:),vel(:),F(:),vir9(9),h(9),hi(9)
    integer, intent(in) :: Na,Nb,b_i0(:),b_i1(:),nlc(:),nll(:),mol_id(:)
    double precision, intent(in) :: b_kb(:),b_r0(:),dt
    type(NPTState), intent(inout) :: npt
    double precision, intent(out) :: Eb_out,Elj_out,Etot_out,KE_out
    double precision :: hdt,V,KE,dP,eps,sc_nh,sc_pr,sc_v,mi_inv
    double precision :: px,py,pz,vx,vy,vz,sx,sy,sz,vsx,vsy,vsz,eps2,sv2,V2
    integer :: i,a,idx
    hdt=.5d0*dt; V=abs(mat_det9(h)); KE=ke_total(vel,Na)
    npt%xi=npt%xi+hdt*(2*KE-dble(npt%Nf)*kB*npt%Tt)/npt%Q
    npt%xi=clamp_val(npt%xi,-.05d0,.05d0)
    dP=inst_P(vir9,KE,V)-npt%Pe
    do a=0,2;npt%Vg(a*4+1)=npt%Vg(a*4+1)+hdt*V*dP/(npt%W_*eV2GPa)
      npt%Vg(a*4+1)=clamp_val(npt%Vg(a*4+1),-.005d0,.005d0);end do
    eps=npt%Vg(1)*hi(1)+npt%Vg(5)*hi(5)+npt%Vg(9)*hi(9)
    sc_nh=exp(-hdt*npt%xi);sc_pr=exp(-hdt*eps/3);sc_v=sc_nh*sc_pr;mi_inv=CONV/mC
    !$OMP PARALLEL DO PRIVATE(i,idx)
    do i=1,Na;idx=(i-1)*3
      vel(idx+1)=vel(idx+1)*sc_v+hdt*F(idx+1)*mi_inv
      vel(idx+2)=vel(idx+2)*sc_v+hdt*F(idx+2)*mi_inv
      vel(idx+3)=vel(idx+3)*sc_v+hdt*F(idx+3)*mi_inv
    end do
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO PRIVATE(i,idx,px,py,pz,vx,vy,vz,sx,sy,sz,vsx,vsy,vsz)
    do i=1,Na;idx=(i-1)*3
      px=pos(idx+1);py=pos(idx+2);pz=pos(idx+3)
      vx=vel(idx+1);vy=vel(idx+2);vz=vel(idx+3)
      sx=hi(1)*px+hi(2)*py+hi(3)*pz;sy=hi(4)*px+hi(5)*py+hi(6)*pz;sz=hi(7)*px+hi(8)*py+hi(9)*pz
      vsx=hi(1)*vx+hi(2)*vy+hi(3)*vz;vsy=hi(4)*vx+hi(5)*vy+hi(6)*vz;vsz=hi(7)*vx+hi(8)*vy+hi(9)*vz
      sx=sx+dt*vsx;sy=sy+dt*vsy;sz=sz+dt*vsz;sx=sx-floor(sx);sy=sy-floor(sy);sz=sz-floor(sz)
      pos(idx+1)=sx;pos(idx+2)=sy;pos(idx+3)=sz
    end do
    !$OMP END PARALLEL DO
    do a=0,2; do i=1,3; h(a*3+i)=h(a*3+i)+dt*npt%Vg(a*3+i); end do; end do
    call mat_inv9(h,hi)
    !$OMP PARALLEL DO PRIVATE(i,idx,sx,sy,sz)
    do i=1,Na;idx=(i-1)*3;sx=pos(idx+1);sy=pos(idx+2);sz=pos(idx+3)
      pos(idx+1)=h(1)*sx+h(2)*sy+h(3)*sz;pos(idx+2)=h(4)*sx+h(5)*sy+h(6)*sz
      pos(idx+3)=h(7)*sx+h(8)*sy+h(9)*sz
    end do
    !$OMP END PARALLEL DO
    call compute_forces(F,vir9,pos,h,hi,Na,Nb,b_i0,b_i1,b_kb,b_r0,nlc,nll,mol_id,Eb_out,Elj_out,Etot_out)
    eps2=npt%Vg(1)*hi(1)+npt%Vg(5)*hi(5)+npt%Vg(9)*hi(9)
    sv2=sc_nh*exp(-hdt*eps2/3)
    !$OMP PARALLEL DO PRIVATE(i,idx)
    do i=1,Na;idx=(i-1)*3
      vel(idx+1)=(vel(idx+1)+hdt*F(idx+1)*mi_inv)*sv2
      vel(idx+2)=(vel(idx+2)+hdt*F(idx+2)*mi_inv)*sv2
      vel(idx+3)=(vel(idx+3)+hdt*F(idx+3)*mi_inv)*sv2
    end do
    !$OMP END PARALLEL DO
    KE=ke_total(vel,Na); KE_out=KE
    npt%xi=npt%xi+hdt*(2*KE-dble(npt%Nf)*kB*npt%Tt)/npt%Q
    npt%xi=clamp_val(npt%xi,-.05d0,.05d0)
    V2=abs(mat_det9(h));dP=inst_P(vir9,KE,V2)-npt%Pe
    do a=0,2;npt%Vg(a*4+1)=npt%Vg(a*4+1)+hdt*V2*dP/(npt%W_*eV2GPa)
      npt%Vg(a*4+1)=clamp_val(npt%Vg(a*4+1),-.005d0,.005d0);end do
  end subroutine
  double precision function gauss_rand()
    double precision :: u1,u2; call random_number(u1); call random_number(u2)
    if(u1<1d-30) u1=1d-30; gauss_rand=sqrt(-2*log(u1))*cos(2*PI_*u2)
  end function
  subroutine get_opt_int(key,defval,res)
    character(len=*),intent(in)::key;integer,intent(in)::defval;integer,intent(out)::res
    integer::i,n,kl,ios;character(len=256)::arg;res=defval;n=command_argument_count();kl=len_trim(key)
    do i=1,n;call get_command_argument(i,arg)
      if(arg(1:kl+3)=='--'//key(1:kl)//'=') read(arg(kl+4:),*,iostat=ios) res
    end do
  end subroutine
  subroutine get_opt_dbl(key,defval,res)
    character(len=*),intent(in)::key;double precision,intent(in)::defval;double precision,intent(out)::res
    integer::i,n,kl,ios;character(len=256)::arg;res=defval;n=command_argument_count();kl=len_trim(key)
    do i=1,n;call get_command_argument(i,arg)
      if(arg(1:kl+3)=='--'//key(1:kl)//'=') read(arg(kl+4:),*,iostat=ios) res
    end do
  end subroutine
end module fuller_mmmd_mod

program fuller_LJ_npt_mmmd
  use fuller_mmmd_mod
  implicit none
  integer :: nc,Nmol,natom,Na,nsteps,nlup,Nmax,i,j,a,gstep,nav
  integer :: coldstart,warmup,avg_from,avg_to,total_steps,gavg_from,gavg_to
  integer :: nrec_o,nrec_rst,start_step,mon_int,prn,prn_pre,cur_prn
  integer :: Nb_mol,Nb
  double precision :: T,Pe,dt,a0,Rmol,Dmol,T_cold,T_init,sv,scale_v,tgt
  double precision :: Eb,Elj,Etot,KE,Tn,Pn,V_val,an_val
  double precision :: sT,sP,sa,sEb,sElj,sEt,t_start,t_now,elapsed
  double precision :: ff_kb_kcal,ff_kb
  double precision :: h(9),hi(9),vir9(9),vcm(3)
  double precision :: mol_coords(84,3)
  integer :: mol_bonds(200,2)
  double precision, allocatable :: pos(:),vel(:),F(:)
  integer, allocatable :: mol_id(:),nlc(:),nll(:)
  integer, allocatable :: b_i0(:),b_i1(:)
  double precision, allocatable :: b_kb_arr(:),b_r0(:)
  double precision, allocatable :: mol_centers(:)
  type(NPTState) :: npt
  T_cold=4.0d0
  call get_opt_int('cell',3,nc); call get_opt_dbl('temp',298.0d0,T)
  call get_opt_dbl('pres',0.0d0,Pe); call get_opt_int('step',10000,nsteps)
  call get_opt_dbl('dt',0.1d0,dt); call get_opt_int('coldstart',0,coldstart)
  call get_opt_int('warmup',0,warmup); call get_opt_int('from',0,avg_from)
  call get_opt_int('to',0,avg_to); call get_opt_int('ovito',0,nrec_o)
  call get_opt_int('restart',0,nrec_rst); call get_opt_int('mon',0,mon_int)
  call get_opt_dbl('ff_kb',469.0d0,ff_kb_kcal)
  ff_kb=ff_kb_kcal*kcal2eV
  if(avg_to<=0) avg_to=nsteps; if(avg_from<=0) avg_from=max(1,nsteps-nsteps/4)
  total_steps=coldstart+warmup+nsteps
  gavg_from=coldstart+warmup+avg_from; gavg_to=coldstart+warmup+avg_to
  start_step=0; nlup=20

  call generate_c60_mmmd(mol_coords,natom,Rmol,Dmol,Nb_mol,mol_bonds)
  a0=Dmol*sqrt(2.0d0)*1.4d0
  Nmax=4*nc*nc*nc; Nmol=0
  allocate(mol_centers(Nmax*3)); mol_centers=0
  Nmol=make_fcc(a0,nc,mol_centers,h)
  Na=Nmol*natom
  allocate(pos(Na*3),vel(Na*3),F(Na*3))
  allocate(mol_id(Na),nlc(Na),nll(Na*MAX_LJ_NEIGH))
  Nb=Nb_mol*Nmol
  allocate(b_i0(Nb),b_i1(Nb),b_kb_arr(Nb),b_r0(Nb))
  pos=0;vel=0;F=0;h=0;hi=0;vir9=0;nlc=0;nll=0

  ! 原子座標と分子IDの設定
  h=0
  Nmol=make_fcc(a0,nc,mol_centers,h)
  do i=1,Nmol; do a=1,natom
    j=(i-1)*natom+a
    pos((j-1)*3+1)=mol_centers((i-1)*3+1)+mol_coords(a,1)
    pos((j-1)*3+2)=mol_centers((i-1)*3+2)+mol_coords(a,2)
    pos((j-1)*3+3)=mol_centers((i-1)*3+3)+mol_coords(a,3)
    mol_id(j)=i
  end do; end do
  deallocate(mol_centers)

  ! トポロジー構築 (Bond)
  do i=1,Nmol; do j=1,Nb_mol
    a=(i-1)*Nb_mol+j
    b_i0(a)=(i-1)*natom+mol_bonds(j,1)
    b_i1(a)=(i-1)*natom+mol_bonds(j,2)
    b_kb_arr(a)=ff_kb
    b_r0(a)=sqrt( &
      (mol_coords(mol_bonds(j,1),1)-mol_coords(mol_bonds(j,2),1))**2 + &
      (mol_coords(mol_bonds(j,1),2)-mol_coords(mol_bonds(j,2),2))**2 + &
      (mol_coords(mol_bonds(j,1),3)-mol_coords(mol_bonds(j,2),3))**2)
  end do; end do

  call mat_inv9(h,hi)
  write(*,'(A)') '========================================================================'
  write(*,'(A)') '  Fullerene Crystal NPT-MD — MM Force Field (Serial, Fortran 95)'
  write(*,'(A)') '========================================================================'
  write(*,'(A,I0,A)') '  Atoms/molecule  : ',natom
  write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '  FCC ',nc,'x',nc,'x',nc,'  Nmol=',Nmol,'  Natom=',Na
  write(*,'(A,F8.3,A,F6.1,A,F6.4,A,F5.3,A)') '  a0=',a0,' A  T=',T,' K  P=',Pe,' GPa  dt=',dt,' fs'
  write(*,'(A,I0,A,I0)') '  Bonds/mol=',Nb_mol,'  Total bonds=',Nb
  write(*,'(A,I0,A,I0,A,I0)') '  Production: ',nsteps,' steps  Total=',total_steps
  write(*,'(A)') '========================================================================'
  write(*,*)

  T_init=T; if(coldstart>0.or.warmup>0) T_init=T_cold
  call random_seed(); sv=sqrt(kB*T_init*CONV/mC)
  do i=1,Na; do a=1,3; vel((i-1)*3+a)=sv*gauss_rand(); end do; end do
  vcm=0; do i=1,Na; vcm(1)=vcm(1)+vel((i-1)*3+1); vcm(2)=vcm(2)+vel((i-1)*3+2)
    vcm(3)=vcm(3)+vel((i-1)*3+3); end do
  vcm=vcm/dble(Na)
  do i=1,Na; vel((i-1)*3+1)=vel((i-1)*3+1)-vcm(1); vel((i-1)*3+2)=vel((i-1)*3+2)-vcm(2)
    vel((i-1)*3+3)=vel((i-1)*3+3)-vcm(3); end do

  call make_npt(npt,T,Pe,Na); npt%Tt=T_init
  call build_nlist_lj(pos,h,hi,Na,mol_id,nlc,nll)
  call apply_pbc(pos,h,hi,Na)
  call compute_forces(F,vir9,pos,h,hi,Na,Nb,b_i0,b_i1,b_kb_arr,b_r0,nlc,nll,mol_id,Eb,Elj,Etot)

  prn=mon_int; if(prn<=0) prn=max(1,total_steps/50)
  prn_pre=prn; if(coldstart+warmup>0) prn_pre=max(1,(coldstart+warmup)/100)
  sT=0;sP=0;sa=0;sEb=0;sElj=0;sEt=0;nav=0
  call cpu_time(t_start)

  write(*,'(A8,A6,A8,A10,A9,A10,A10,A10,A8)') &
    'step','phase','T[K]','P[GPa]','a[A]','E_bond','E_LJ','E_total','t[s]'

  do gstep=start_step+1,total_steps
    if(gstep<=coldstart) then; npt%Tt=T_cold
    else if(gstep<=coldstart+warmup) then
      npt%Tt=T_cold+(T-T_cold)*dble(gstep-coldstart)/dble(max(warmup,1))
    else; npt%Tt=T; end if
    if(coldstart>0.and.gstep==coldstart+1) then; npt%xi=0; npt%Vg=0; end if
    if(gstep<=coldstart) npt%Vg=0
    cur_prn=prn; if(gstep<=coldstart+warmup) cur_prn=prn_pre

    if(mod(gstep,nlup)==0) then
      call mat_inv9(h,hi); call build_nlist_lj(pos,h,hi,Na,mol_id,nlc,nll)
    end if

    call step_npt_mm(pos,vel,F,vir9,h,hi,Na,dt,npt, &
         Nb,b_i0,b_i1,b_kb_arr,b_r0,nlc,nll,mol_id,Eb,Elj,Etot,KE)

    V_val=abs(mat_det9(h)); Tn=inst_T(KE,npt%Nf); Pn=inst_P(vir9,KE,V_val)

    if((gstep<=coldstart.or.gstep<=coldstart+warmup).and.Tn>0.1d0) then
      tgt=T_cold; if(gstep>coldstart) tgt=npt%Tt
      scale_v=sqrt(max(tgt,0.1d0)/Tn)
      !$OMP PARALLEL DO PRIVATE(i)
      do i=1,Na*3; vel(i)=vel(i)*scale_v; end do
      !$OMP END PARALLEL DO
      KE=ke_total(vel,Na); Tn=inst_T(KE,npt%Nf); npt%xi=0
      if(gstep<=coldstart) npt%Vg=0
    end if

    an_val=h(1)/dble(nc)
    if(gstep>=gavg_from.and.gstep<=gavg_to) then
      sT=sT+Tn;sP=sP+Pn;sa=sa+an_val;sEb=sEb+Eb/Nmol;sElj=sElj+Elj/Nmol
      sEt=sEt+Etot/Nmol;nav=nav+1
    end if

    if(mod(gstep,cur_prn)==0.or.gstep==total_steps) then
      call cpu_time(t_now); elapsed=t_now-t_start
      if(gstep<=coldstart) then
        write(*,'(I8,A6,F8.1,F10.3,F9.3,F10.4,F10.4,F10.4,F8.0)') gstep,' COLD',Tn,Pn,an_val,Eb/Nmol,Elj/Nmol,Etot/Nmol,elapsed
      else if(gstep<=coldstart+warmup) then
        write(*,'(I8,A6,F8.1,F10.3,F9.3,F10.4,F10.4,F10.4,F8.0)') gstep,' WARM',Tn,Pn,an_val,Eb/Nmol,Elj/Nmol,Etot/Nmol,elapsed
      else
        write(*,'(I8,A6,F8.1,F10.3,F9.3,F10.4,F10.4,F10.4,F8.0)') gstep,' PROD',Tn,Pn,an_val,Eb/Nmol,Elj/Nmol,Etot/Nmol,elapsed
      end if
    end if
  end do

  if(nav>0) then
    write(*,*); write(*,'(A)') '========================================================================'
    write(*,'(A,I0,A,F7.2,A,F8.4,A,F8.4,A,F10.4,A,F10.4,A,F10.4)') &
      '  Avg(',nav,'): T=',sT/nav,' P=',sP/nav,' a=',sa/nav,' bond=',sEb/nav,' LJ=',sElj/nav,' tot=',sEt/nav
    write(*,'(A)') '========================================================================'
  end if
  call cpu_time(t_now); write(*,'(A,F8.1,A)') '  Done (',t_now-t_start,' sec)'
  deallocate(pos,vel,F,mol_id,nlc,nll,b_i0,b_i1,b_kb_arr,b_r0)
end program fuller_LJ_npt_mmmd
