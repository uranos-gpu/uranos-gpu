module post_output_module
use storage_module
use post_storage_module
use post_output_gnu_module
use post_output_vtk_module
use m_npy
use FileModule
  
implicit none
private
public write_post_treatments

contains
subroutine write_post_treatments
        implicit none
        character(7)   :: iteration
        character(100) :: path

        write(iteration,'(i7.7)') it

        call write_screen
        if(gnuplot) then
          call write_gnuplot
          call write_vs_time
        endif

        if(VTK)   call write_vtk(dims,'VTK')
        if(VTK2D) call write_vtk(2,'VTK2D')

        if(NPY) then
            path = 'DATA/'//trim(data_dir)//'/NPY/' 
            call save_npy(trim(path)//'phi'//trim(iteration)//'.npy',phi(sx:ex,sy:ey,sz:ez,:))
            !if(it.eq.itmax) then
            !call save_npy(trim(path)//'x.npy',x(sx:ex))
            !call save_npy(trim(path)//'y.npy',y(sy:ey))
            !call save_npy(trim(path)//'z.npy',z(sz:ez))
            !endif
        endif

        call write_xz_plane_vtk(3)

        return
end subroutine write_post_treatments


subroutine write_screen
        implicit none
        real(rp) :: M_max, r_max, T_max, p_max, c_max, t_min, rho_min
        real(rp) :: rho_, i_rho, u_, v_, w_, T_, p_

        M_max = 0._rp
        r_max = 0._rp
        T_max = 0._rp
        p_max = 0._rp
        c_max = 0._rp

        t_min   = 10E+14
        rho_min = 10E+14
        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 ! compute speeds
                 rho_  = phi(i,j,k,1)
                 i_rho = 1._rp/rho_

                 u_ = phi(i,j,k,2)*i_rho
                 v_ = phi(i,j,k,3)*i_rho
                 w_ = phi(i,j,k,4)*i_rho

                 p_=(gamma0-1._rp)*(phi(i,j,k,5)-0.5_rp*rho_*(u_**2 + v_**2 + w_**2))

                 T_ = p_ * i_rho

                 rho_min = min(rho_, rho_min)
                   t_min = min(  t_,   t_min)

                 T_max = max(T_  ,T_max) 
                 r_max = max(rho_,r_max)
                 p_max = max(p_  ,p_max)
                        
                 M_max = max( (u_**2 + v_**2 + w_**2)/(gamma0*T_), M_max)

              enddo
           enddo
        enddo
        M_max = sqrt(M_max)
        c_max = sqrt(gamma0*T_max)

        ! ---- write screen output ---- !
        if(mod(it,10*itout)==0) then
               write(*,'(4x,13(A,4x))') ' it', ' file %','time %', 'dt ', '        M max',&
                       '      rho max',  '    rho_min', '    T max', '      T_min', '      p max', '      c max'
               write(*,*) ' -----------------------------------------------------------------------------------', &
                '------------------------------------------------------ post Uranos! '
        endif

        write(*,10) it ,real(f,rp)/real(size(main_dir%content),rp)*100,time/tmax*100 , dt, M_max, r_max, rho_min, T_max, &
                t_min, p_max, c_max
        10 format(I10,f8.3,f8.3,8(e15.3))

        return
end subroutine write_screen


end module post_output_module
