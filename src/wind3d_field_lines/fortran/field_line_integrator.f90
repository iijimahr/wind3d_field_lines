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
