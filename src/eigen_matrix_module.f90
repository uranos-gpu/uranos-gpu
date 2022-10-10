module eigen_matrix_module
use parameters_module, only: rp, gamma0
use storage_module   , only: eqs, toll_equality
use mpi_module

implicit none
private
public eigenmatrix_x, eigenmatrix_y, eigenmatrix_z
contains

subroutine eigenmatrix_x(u_,v_,w_,c_,ek_,h_,right,left)
        implicit none
        real(rp)                    , intent(in ) :: u_, v_, w_, c_, ek_, h_
        real(rp), dimension(eqs,eqs), intent(out) :: right, left
        real(rp), parameter    :: hatGamma = gamma0-1.0_rp
        real(rp) :: cc, un_c2, vnc, hgU, hgV, hgW, hgE
        
#ifdef TIME
        call mpi_stime(s_eix_time)
#endif

        vnc = u_*c_
        hgU = hatGamma*u_
        hgV = hatGamma*v_
        hgW = hatGamma*w_
        hgE = hatGamma*ek_

        right(1,:)=[1.0_rp       , 1.0_rp , 1.0_rp       , 0.0_rp , 0.0_rp]
        right(2,:)=[u_-c_     , u_  , u_+c_     , 0.0_rp , 0.0_rp]
        right(3,:)=[v_        , v_  , v_        ,-1.0_rp , 0.0_rp]
        right(4,:)=[w_        , w_  , w_        , 0.0_rp , 1.0_rp]
        right(5,:)=[h_-vnc    , ek_ , h_+vnc    , -v_ , w_ ]

        cc    = c_*c_
        un_c2 = 1.0_rp/cc

        left(1,:)=0.5_rp*un_c2*[hgE + vnc  , -hgU-c_ , -hgV , -hgW,  hatGamma]
        left(2,:)=    un_c2*[cc  - hgE  ,  hgU    ,  hgV ,  hgW, -hatGamma]
        left(3,:)=0.5_rp*un_c2*[hgE - vnc  , -hgU+c_ , -hgV , -hgW,  hatGamma]
        left(4,:)=          [v_         , 0.0_rp     , -1.0_rp ,  0.0_rp,  0.0_rp     ]
        left(5,:)=          [-w_        , 0.0_rp     ,  0.0_rp ,  1.0_rp,  0.0_rp     ]

#ifdef DEBUG
        call check_matrix(left,right)
#endif

#ifdef TIME
        call mpi_etime(s_eix_time,t_eix_calls,t_eix_time)
#endif

return
end subroutine eigenmatrix_x

subroutine eigenmatrix_y(u_,v_,w_,c_,ek_,h_,right,left)
! Right eigen matrix      
! REF: "Eigenvalues and eigenvectors of the euler equation". Axel Rohde
! attenzione, in fondo al paper dice quali matrici usare per non incorrere in esplosioni lungo le dir.
        implicit none
        real(rp)                    , intent(in ) :: u_, v_,w_, c_, ek_, h_
        real(rp), dimension(eqs,eqs), intent(out) :: right, left
        real(rp), parameter    :: hatGamma = gamma0-1.0_rp
        real(rp) :: cc, un_c2, vnc, hgU, hgV, hgW, hgE

#ifdef TIME
        call mpi_stime(s_eiy_time)
#endif
        vnc = v_*c_
        hgU = hatGamma*u_
        hgV = hatGamma*v_
        hgW = hatGamma*w_
        hgE = hatGamma*ek_

        right(1,:)=[1.0_rp    , 1.0_rp, 1.0_rp       ,  0.0_rp , 0.0_rp  ]
        right(2,:)=[u_     , u_ , u_        ,  1.0_rp , 0.0_rp  ]
        right(3,:)=[v_-c_  , v_ , v_+c_     ,  0.0_rp , 0.0_rp  ]
        right(4,:)=[w_     , w_ , w_        ,  0.0_rp , -1.0_rp ]
        right(5,:)=[h_-vnc , ek_, h_+vnc    , u_   , -w_  ]

        ! Left eigen matrix 
        cc    = c_*c_
        un_c2 = 1.0_rp/cc

        left(1,:)=0.5_rp*un_c2*[hgE+vnc , -hgU , -hgV-c_ , -hgW , hatGamma]
        left(2,:)=    un_c2*[cc -hgE ,  hgU ,  hgV    ,  hgW ,-hatGamma]
        left(3,:)=0.5_rp*un_c2*[hgE-vnc , -hgU , -hgV+c_ , -hgW , hatGamma]
        left(4,:)=          [-u_     ,  1.0_rp ,  0.0_rp    ,  0.0_rp , 0.0_rp     ]
        left(5,:)=          [ w_     ,  0.0_rp ,  0.0_rp    , -1.0_rp , 0.0_rp     ]

#ifdef DEBUG
        call check_matrix(left,right)
#endif
#ifdef TIME
        call mpi_etime(s_eiy_time,t_eiy_calls,t_eiy_time)
#endif
        return
end subroutine eigenmatrix_y

subroutine eigenmatrix_z(u_,v_,w_,c_,ek_,h_,right,left)
        implicit none
        real(rp)                    , intent(in ) :: u_, v_,w_, c_, ek_, h_
        real(rp), dimension(eqs,eqs), intent(out) :: right, left
        real(rp), parameter    :: hatGamma = gamma0-1.0_rp
        real(rp) :: cc, un_c2, vnc, hgU, hgV, hgW, hgE

#ifdef TIME
        call mpi_stime(s_eiz_time)
#endif
        vnc = w_*c_
        hgU = hatGamma*u_
        hgV = hatGamma*v_
        hgW = hatGamma*w_
        hgE = hatGamma*ek_

        right(1,:)=[1.0_rp          , 1.0_rp, 1.0_rp       ,  0.0_rp          , 0.0_rp   ]
        right(2,:)=[u_           , u_ , u_        , -1.0_rp          , 0.0_rp   ]
        right(3,:)=[v_           , v_ , v_        ,  0.0_rp          , 1.0_rp   ]
        right(4,:)=[w_-c_        , w_ , w_+c_     ,  0.0_rp          , 0.0_rp   ]
        right(5,:)=[h_-vnc       , ek_, h_+vnc    , -u_           , v_    ]

        ! Left eigen matrix 
        cc    = c_*c_
        un_c2 = 1.0_rp/cc

        left(1,:)=0.5_rp*un_c2*[hgE + vnc    , -hgU, -hgV, -hgW-c_, hatGamma]
        left(2,:)=    un_c2*[cc  - hgE    ,  hgU,  hgV,  hgW   ,-hatGamma]
        left(3,:)=0.5_rp*un_c2*[hgE - vnc    , -hgU, -hgV, -hgW+c_, hatGamma]
        left(4,:)=          [u_           , -1.0_rp,  0.0_rp,  0.0_rp   , 0.0_rp     ]
        left(5,:)=          [-v_          ,  0.0_rp,  1.0_rp,  0.0_rp   , 0.0_rp     ]

#ifdef DEBUG
        call check_matrix(left,right)
#endif

#ifdef TIME
        call mpi_etime(s_eiz_time,t_eiz_calls,t_eiz_time)
#endif
        return
end subroutine eigenmatrix_z



subroutine check_matrix(left,right)
        real(rp), dimension(eqs,eqs), intent(in) :: left,right
        real(rp), dimension(eqs,eqs)             :: id
        real(rp), parameter                      :: toll = 1.0E-12_rp
        integer                                  :: h,k

        id = matmul(right, left)
        do k = 1,eqs
           do h = 1, eqs
        
              if(h.eq.k.and.abs(id(h,k)-1.0_rp).ge.toll) then
                              print*, 'in diagonal position', h,k, 'the WENO identity matrix does not value 1'
                              print*, 'identity value = ', id(h,k)
                              call secure_stop
              endif
              if(h.ne.k.and.abs(id(h,k)).ge.toll) then
                              print*, 'in extra diagonal position', h,k, 'the WENO identity matrix does not value 0'
                              print*, 'identity value = ', id(h,k)
                              call secure_stop
              endif
           
          enddo
        enddo
return
end subroutine check_matrix




end module eigen_matrix_module
