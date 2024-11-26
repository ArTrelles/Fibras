      program read_test
      implicit none
      real(8) :: ao(6), ae(6), po(2), pe(2), apo(2), ape(2)
      real(8) :: xo_p,yo_p,zo_p,vx_pn,vy_pn, vz_pn, vxo_p, vyo_p,vzo_p
      real(8) :: mgen, z, PAge, priot, r7, gr7,z7,lo_p,bo_p,r_n,gr_n
      integer :: iden

      open(unit=10, file="GeneratedStarsle03.bin", form="unformatted",&
       status="old")

      do
      read(10, end=100) ao, ae, po, pe, apo, ape, xo_p, yo_p, zo_p,&
      vx_pn, vy_pn, vz_pn, &
      vxo_p, vyo_p, vzo_p, mgen, z, PAge, priot, iden, &
      r7, gr7, z7, lo_p,bo_p, r_n, gr_n

        print *, 'Valores leídos:'
        print *, 'ao:', ao
        print *, 'ae:', ae
        ! Agregar más impresiones si es necesario
      end do

100   continue
      close(10)
      end program read_test
