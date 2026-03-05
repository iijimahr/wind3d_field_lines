subroutine bbtobln(i_bln,j_bln,k_bln,lmin_bln,lmax_bln, &
     & lcen_bln,icen_bln,jcen_bln,kcen_bln, &
     & bx,by,bz,dx,dy,dz,nsubstepx,lx_bln,nx_bln,ix,jx,kx,margin)
  implicit none
  real(8),dimension(nx_bln,lx_bln),intent(INOUT) :: i_bln,j_bln,k_bln
  integer,dimension(nx_bln),intent(INOUT) :: lmin_bln,lmax_bln
  integer,intent(IN) :: lcen_bln
  real(8),dimension(nx_bln),intent(IN) :: icen_bln,jcen_bln,kcen_bln
  real(8),dimension(ix,jx,kx),intent(IN) :: bx,by,bz
  real(8),dimension(kx),intent(IN) :: dx,dy,dz
  integer,intent(IN) :: nsubstepx
  integer,intent(IN) :: lx_bln,nx_bln
  integer,intent(IN) :: ix,jx,kx,margin

  integer,parameter :: niterx=5
  real(8),parameter :: bb_eps=1.d-40

  integer :: niter,nsubstep
  integer :: i,j,k,l,n
  integer :: ks,kb0,kb1
  real(8) :: fkb0,fkb1
  real(8) :: i_s,j_s,k_s
  real(8) :: i_e,j_e,k_e
  real(8) :: bnx_s,bny_s,bnz_s
  real(8) :: bnx_e,bny_e,bnz_e
  real(8) :: dx_s,dy_s,dz_s,ds_s
  real(8) :: dx_e,dy_e,dz_e,ds_e
  real(8),dimension(ix,jx,kx) :: bnx,bny,bnz
  real(8) :: bbi0

  ks=min(1,kx-1)
  kb0=1+margin*ks
  kb1=kx-margin*ks
  fkb0=real(kb0,8)-0.5d0
  fkb1=real(kb1,8)+0.5d0

  do n=1,nx_bln
     i_bln(n,1:lx_bln)=icen_bln(n)
     j_bln(n,1:lx_bln)=jcen_bln(n)
     k_bln(n,1:lx_bln)=kcen_bln(n)
  end do

  do k=1,kx
     do j=1,jx
        do i=1,ix
           bbi0=1.d0/max(bb_eps,sqrt(bx(i,j,k)**2 +by(i,j,k)**2+bz(i,j,k)**2))
           bnx(i,j,k)=bbi0*bx(i,j,k)
           bny(i,j,k)=bbi0*by(i,j,k)
           bnz(i,j,k)=bbi0*bz(i,j,k)
        end do
     end do
  end do

  lmin_bln(1:nx_bln)=1
  lmax_bln(1:nx_bln)=lx_bln

  do n=1,nx_bln
     l_forward: do l=lcen_bln+1,lx_bln
        if (lmax_bln(n)<lx_bln) cycle

        i_s=i_bln(n,l-1)
        j_s=j_bln(n,l-1)
        k_s=k_bln(n,l-1)

        call interp_k(dx_s,k_s,dx,kx)
        call interp_k(dy_s,k_s,dy,kx)
        call interp_k(dz_s,k_s,dz,kx)

        call interp_ijk(bnx_s,i_s,j_s,k_s,bnx,ix,jx,kx,margin)
        call interp_ijk(bny_s,i_s,j_s,k_s,bny,ix,jx,kx,margin)
        call interp_ijk(bnz_s,i_s,j_s,k_s,bnz,ix,jx,kx,margin)

        ds_s=dz_s/real(nsubstepx,8)

        i_e=i_s+ds_s*bnx_s/dx_s
        j_e=j_s+ds_s*bny_s/dy_s
        k_e=k_s+ds_s*bnz_s/dz_s

        do nsubstep=1,nsubstepx
           if (k_e<fkb0.or.k_e>fkb1) then
              lmax_bln(n)=l-1
              cycle l_forward
           end if

           do niter=1,niterx
              call interp_k(dx_e,k_e,dx,kx)
              call interp_k(dy_e,k_e,dy,kx)
              call interp_k(dz_e,k_e,dz,kx)

              call interp_ijk(bnx_e,i_e,j_e,k_e,bnx,ix,jx,kx,margin)
              call interp_ijk(bny_e,i_e,j_e,k_e,bny,ix,jx,kx,margin)
              call interp_ijk(bnz_e,i_e,j_e,k_e,bnz,ix,jx,kx,margin)

              ds_e=dz_e/real(nsubstepx,8)

              i_e=i_s+0.5d0*(ds_s*bnx_s/dx_s+ds_e*bnx_e/dx_e)
              j_e=j_s+0.5d0*(ds_s*bny_s/dy_s+ds_e*bny_e/dy_e)
              k_e=k_s+0.5d0*(ds_s*bnz_s/dz_s+ds_e*bnz_e/dz_e)
           end do

           i_s=i_e
           j_s=j_e
           k_s=k_e

           call interp_k(dx_s,k_s,dx,kx)
           call interp_k(dy_s,k_s,dy,kx)
           call interp_k(dz_s,k_s,dz,kx)

           call interp_ijk(bnx_s,i_s,j_s,k_s,bnx,ix,jx,kx,margin)
           call interp_ijk(bny_s,i_s,j_s,k_s,bny,ix,jx,kx,margin)
           call interp_ijk(bnz_s,i_s,j_s,k_s,bnz,ix,jx,kx,margin)
        end do

        i_bln(n,l)=i_e
        j_bln(n,l)=j_e
        k_bln(n,l)=k_e
     end do l_forward

     l_backward: do l=lcen_bln-1,1,-1
        if (lmin_bln(n)>1) cycle

        i_s=i_bln(n,l+1)
        j_s=j_bln(n,l+1)
        k_s=k_bln(n,l+1)

        call interp_k(dx_s,k_s,dx,kx)
        call interp_k(dy_s,k_s,dy,kx)
        call interp_k(dz_s,k_s,dz,kx)

        call interp_ijk(bnx_s,i_s,j_s,k_s,bnx,ix,jx,kx,margin)
        call interp_ijk(bny_s,i_s,j_s,k_s,bny,ix,jx,kx,margin)
        call interp_ijk(bnz_s,i_s,j_s,k_s,bnz,ix,jx,kx,margin)

        ds_s=dz_s/real(nsubstepx,8)

        i_e=i_s-ds_s*bnx_s/dx_s
        j_e=j_s-ds_s*bny_s/dy_s
        k_e=k_s-ds_s*bnz_s/dz_s

        do nsubstep=1,nsubstepx
           if (k_e<fkb0.or.k_e>fkb1) then
              lmin_bln(n)=l+1
              cycle l_backward
           end if

           do niter=1,niterx
              call interp_k(dx_e,k_e,dx,kx)
              call interp_k(dy_e,k_e,dy,kx)
              call interp_k(dz_e,k_e,dz,kx)

              call interp_ijk(bnx_e,i_e,j_e,k_e,bnx,ix,jx,kx,margin)
              call interp_ijk(bny_e,i_e,j_e,k_e,bny,ix,jx,kx,margin)
              call interp_ijk(bnz_e,i_e,j_e,k_e,bnz,ix,jx,kx,margin)

              ds_e=dz_e/real(nsubstepx,8)

              i_e=i_s-0.5d0*(ds_s*bnx_s/dx_s+ds_e*bnx_e/dx_e)
              j_e=j_s-0.5d0*(ds_s*bny_s/dy_s+ds_e*bny_e/dy_e)
              k_e=k_s-0.5d0*(ds_s*bnz_s/dz_s+ds_e*bnz_e/dz_e)
           end do

           i_s=i_e
           j_s=j_e
           k_s=k_e

           call interp_k(dx_s,k_s,dx,kx)
           call interp_k(dy_s,k_s,dy,kx)
           call interp_k(dz_s,k_s,dz,kx)

           call interp_ijk(bnx_s,i_s,j_s,k_s,bnx,ix,jx,kx,margin)
           call interp_ijk(bny_s,i_s,j_s,k_s,bny,ix,jx,kx,margin)
           call interp_ijk(bnz_s,i_s,j_s,k_s,bnz,ix,jx,kx,margin)
        end do

        i_bln(n,l)=i_e
        j_bln(n,l)=j_e
        k_bln(n,l)=k_e
     end do l_backward
  end do

  return
end subroutine bbtobln

subroutine interp_ijk(ut,it,jt,kt,u,ix,jx,kx,margin)
  implicit none
  real(8),intent(OUT) :: ut
  real(8),intent(IN) :: it,jt,kt
  real(8),dimension(ix,jx,kx),intent(IN) :: u
  integer,intent(IN) :: ix,jx,kx,margin

  integer :: i,j,k
  integer :: i0,i1,j0,j1,k0,k1
  real(8) :: fi,fj,fk
  real(8) :: t000,t001,t010,t011,t100,t101,t110,t111
  integer :: is,js,ks
  integer :: ib0,jb0,kb0
  integer :: ib1,jb1,kb1
  integer :: inx,jnx,knx

  is=min(1,ix-1)
  js=min(1,jx-1)
  ks=min(1,kx-1)
  ib0=1+margin*is
  ib1=ix-margin*is
  jb0=1+margin*js
  jb1=jx-margin*js
  kb0=1+margin*ks
  kb1=kx-margin*ks
  inx=ib1-ib0+1
  jnx=jb1-jb0+1
  knx=kb1-kb0+1

  i0=floor(it)
  j0=floor(jt)
  k0=floor(kt)
  i1=i0+1
  j1=j0+1
  k1=k0+1

  i0=ib0+mod(mod(i0-ib0,inx)+inx,inx)
  i1=ib0+mod(mod(i1-ib0,inx)+inx,inx)
  j0=jb0+mod(mod(j0-jb0,jnx)+jnx,jnx)
  j1=jb0+mod(mod(j1-jb0,jnx)+jnx,jnx)
  k0=max(kb0,min(kb1,k0))
  k1=max(kb0,min(kb1,k1))

  fi=it-floor(it)
  fj=jt-floor(jt)
  fk=kt-floor(kt)

  t000=u(i0,j0,k0)
  t001=u(i0,j0,k1)
  t010=u(i0,j1,k0)
  t011=u(i0,j1,k1)
  t100=u(i1,j0,k0)
  t101=u(i1,j0,k1)
  t110=u(i1,j1,k0)
  t111=u(i1,j1,k1)

  ut=(1.d0-fi)*(1.d0-fj)*(1.d0-fk)*t000 + (1.d0-fi)*(1.d0-fj)*fk*t001 + &
      (1.d0-fi)*fj*(1.d0-fk)*t010 + (1.d0-fi)*fj*fk*t011 + &
      fi*(1.d0-fj)*(1.d0-fk)*t100 + fi*(1.d0-fj)*fk*t101 + &
      fi*fj*(1.d0-fk)*t110 + fi*fj*fk*t111

  return
end subroutine interp_ijk

subroutine interp_k(ut,kt,u,kx)
  implicit none
  real(8),intent(OUT) :: ut
  real(8),intent(IN) :: kt
  real(8),dimension(kx),intent(IN) :: u
  integer,intent(IN) :: kx

  integer :: k0,k1
  real(8) :: t0,t1

  k0=floor(kt)
  k1=k0+1

  k0=max(1,min(kx,k0))
  k1=max(1,min(kx,k1))

  t0=(real(k1,8)-kt)
  t1=(kt-real(k0,8))
  ut=t0*u(k0)+t1*u(k1)

  return
end subroutine interp_k

subroutine blntofopn(f_opn,i_bln,j_bln,k_bln,lmin_bln,lmax_bln, &
     & bzt,k_min,k_max,lx_bln,nx_bln,ix,jx,kx,margin)
  implicit none
  real(8),dimension(ix,jx,kx),intent(INOUT) :: f_opn
  real(8),dimension(nx_bln,lx_bln),intent(IN) :: i_bln,j_bln,k_bln
  integer,dimension(nx_bln),intent(IN) :: lmin_bln,lmax_bln
  real(8),dimension(ix,jx,kx),intent(IN) :: bzt
  integer,intent(IN) :: k_min,k_max
  integer,intent(IN) :: lx_bln,nx_bln
  integer,intent(IN) :: ix,jx,kx,margin

  integer,parameter :: di_nb=1
  integer,parameter :: dj_nb=di_nb
  real(8),dimension(3,3) :: kernel

  integer :: is,js
  integer :: ib0,jb0
  integer :: ib1,jb1
  integer :: inx,jnx
  integer :: l,n,i,j
  logical :: is_closed
  integer :: k,k_bln0,k_bln1,ks_bln
  real(8) :: t0,t1
  integer :: i_bln0,j_bln0
  real(8) :: bzt_sum,fbzt_sum

  is=min(1,ix-1)
  js=min(1,jx-1)
  kernel=(1.d0/16.d0)*reshape( &
       & (/1.d0,2.d0,1.d0,2.d0,4.d0,2.d0,1.d0,2.d0,1.d0/),(/3,3/))
  ib0=1+margin*is
  ib1=ix-margin*is
  jb0=1+margin*js
  jb1=jx-margin*js
  inx=ib1-ib0+1
  jnx=jb1-jb0+1

  f_opn(1:ix,1:jx,1:kx)=0.d0

  n_loop: do n=1,nx_bln
     is_closed=(maxval(k_bln(n,lmin_bln(n):lmax_bln(n)))<k_max)
     if (is_closed) cycle n_loop

     l_loop: do l=lmin_bln(n)+1,lmax_bln(n)
        k_bln0=floor(k_bln(n,l-1))
        k_bln1=floor(k_bln(n,l))
        ks_bln=sign(1,k_bln1-k_bln0)
        if (k_min>max(k_bln0,k_bln1)) cycle l_loop
        if (k_max<min(k_bln0,k_bln1)) cycle l_loop
        do k=k_bln0,k_bln1,ks_bln
           t0=(k_bln(n,l)-real(k,8)) &
                & /max(1.d-40,k_bln(n,l)-k_bln(n,l-1))
           t0=max(0.d0,min(1.d0,t0))
           t1=1.d0-t0
           i_bln0=nint(t0*i_bln(n,l)+t1*i_bln(n,l-1))
           j_bln0=nint(t0*j_bln(n,l)+t1*j_bln(n,l-1))
           i_bln0=ib0+mod(mod(i_bln0-ib0,inx)+inx,inx)
           j_bln0=jb0+mod(mod(j_bln0-jb0,jnx)+jnx,jnx)
           do j=-dj_nb,dj_nb
              do i=-di_nb,di_nb
                 f_opn(i_bln0+i,j_bln0+j,k)= &
                      & f_opn(i_bln0+i,j_bln0+j,k) &
                      & +kernel(1+di_nb+i,1+dj_nb+j)
              end do
           end do
        end do
     end do l_loop
  end do n_loop

  do k=k_min,k_max
     do i=1,margin
        f_opn(ib1-i+1,jb0:jb1,k)= &
             & f_opn(ib1-i+1,jb0:jb1,k)+f_opn(ib0-i,jb0:jb1,k)
     end do
     do i=1,margin
        f_opn(ib1+i,jb0:jb1,k)= &
             & f_opn(ib1+i,jb0:jb1,k)+f_opn(ib0+i-1,jb0:jb1,k)
     end do
     do j=1,margin
        f_opn(ib0:ib1,jb0-j,k)= &
             & f_opn(ib0:ib1,jb0-j,k)+f_opn(ib0:ib1,jb1-j+1,k)
     end do
     do j=1,margin
        f_opn(ib0:ib1,jb1+j,k)= &
             & f_opn(ib0:ib1,jb1+j,k)+f_opn(ib0:ib1,jb0+j-1,k)
     end do
  end do

  f_opn(1:ix,1:jx,1:kx)=max(0.d0,min(1.d0,f_opn(1:ix,1:jx,1:kx)))
  do k=1,k_min-1
     f_opn(1:ix,1:jx,k)=1.d0
  end do
  do k=k_max-1,kx
     f_opn(1:ix,1:jx,k)=1.d0
  end do
  do k=k_max,k_min,-1
     if (maxval(f_opn(1:ix,1:jx,k))<1.d-5) then
        f_opn(ib0:ib1,jb0:jb1,k)=f_opn(ib0:ib1,jb0:jb1,k+1)
     end if
  end do
  do k=k_min,k_max
     bzt_sum=sum(bzt(ib0:ib1,jb0:jb1,k))
     fbzt_sum=sum(f_opn(ib0:ib1,jb0:jb1,k) &
          & *bzt(ib0:ib1,jb0:jb1,k))
     f_opn(1:ix,1:jx,k)=f_opn(1:ix,1:jx,k) &
          & *max(0.d0,min(2.d0,bzt_sum/fbzt_sum))
  end do

  return
end subroutine blntofopn

subroutine blntobmap(i_obs,j_obs,dk_obs,i_bln,j_bln,k_bln, &
     & lmin_bln,lmax_bln,lcen_bln,k_obs,lx_bln,nx_bln)
  implicit none
  real(8),dimension(nx_bln),intent(INOUT) :: i_obs,j_obs,dk_obs
  real(8),dimension(nx_bln,lx_bln),intent(IN) :: i_bln,j_bln,k_bln
  integer,dimension(nx_bln),intent(IN) :: lmin_bln,lmax_bln
  integer,intent(IN) :: lcen_bln,k_obs,lx_bln,nx_bln

  integer,parameter :: dl_min=1
  real(8),parameter :: dk_errtol=1.2d0

  integer,dimension(nx_bln) :: l_obs
  integer :: l,n,kcen
  logical :: search_forward

  n_loop: do n = 1,nx_bln
     kcen=nint(k_bln(n,lcen_bln))

     if (k_obs>kcen) then
        search_forward=(k_bln(n,lcen_bln+dl_min)>kcen)
     else if (k_obs<kcen) then
        search_forward=(k_bln(n,lcen_bln+dl_min)<kcen)
     else
        l_obs(n)=-1
        cycle n_loop
     end if

     if (search_forward) then
        l_obs(n)=lcen_bln+dl_min
        do l=lcen_bln+dl_min+1,lmax_bln(n)
           if (abs(k_bln(n,l)-k_obs)<dk_errtol) then
              l_obs(n)=l
              cycle n_loop
           end if
           if (abs(k_bln(n,l_obs(n))-k_obs)>abs(k_bln(n,l)-k_obs)) then
              l_obs(n)=l
           end if
           if (abs(k_bln(n,l)-kcen)<dk_errtol) then
              l_obs(n)=l
              cycle n_loop
           end if
        end do
        l_obs(n)=-1
     else
        l_obs(n)=lcen_bln-dl_min
        do l=lcen_bln-dl_min-1,lmin_bln(n),-1
           if (abs(k_bln(n,l)-k_obs)<dk_errtol) then
              l_obs(n)=l
              cycle n_loop
           end if
           if (abs(k_bln(n,l_obs(n))-k_obs)>abs(k_bln(n,l)-k_obs)) then
              l_obs(n)=l
           end if
           if (abs(k_bln(n,l)-kcen)<dk_errtol) then
              l_obs(n)=l
              cycle n_loop
           end if
        end do
        l_obs(n)=-1
     end if
  end do n_loop

  do n=1,nx_bln
     if (l_obs(n)>0) then
        i_obs(n)=i_bln(n,l_obs(n))
        j_obs(n)=j_bln(n,l_obs(n))
        dk_obs(n)=k_bln(n,l_obs(n))-k_obs
     else
        i_obs(n)=1.d100
        j_obs(n)=1.d100
        dk_obs(n)=1.d100
     end if
  end do

  return
end subroutine blntobmap
