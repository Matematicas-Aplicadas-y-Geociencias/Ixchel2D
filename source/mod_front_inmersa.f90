!
! M\'odulo de frontera  inmersa
!
! Este m\'odulo contiene las subrutinas que permiten incluir s\'olidos en el dominio
! del fluido.
! Se utiliza la media geom\'etrica de los coeficientes de difusi\'on en las caras de los
! vol\'umenes de control como se describe en Patankar 1987.
!
! Autor: J.C. Cajas, K.Figueroa
!
module frontera_inmersa
  !
  use malla, only : mi, nj, DBL       ! dimensiones de la malla y tipo de doble
  use malla, only : xp, yp            ! coordenadas en la malla de la presi\'on
  use malla, only : xu, yv            ! coordenadas en la malla de la velocidad
  !
  ! Variables para fijar las fronteras inmersas
  !
  use ec_momento, only : u, u_ant, v, v_ant
  use ec_momento, only : au, av  
  use ec_momento, only : fuente_con_u, fuente_lin_u
  use ec_momento, only : fuente_con_v, fuente_lin_v
  use ec_momento, only : gamma_momen
  !
  ! Variables para fijar temperatura de las fronteras inmersas
  !
  use ec_energia, only : fuente_con_t, fuente_lin_t
  use ec_energia, only : gamma_energ
  !
  implicit none
  !
  ! Interfaces para poder usar funciones como argumentos en subrutinas
  !
  abstract interface
     !
     function funcion1d(x) result(y)
       !
       import DBL
       implicit none
       !
       real(kind=DBL), intent(in) :: x
       real(kind=DBL)             :: y
     end function funcion1d
     !
  end interface
  !  
contains
  !
  !*******************************************************************
  !
  ! lectura_cuerpo
  !
  ! Subrutina que abre el archivo de configuraci'on de fronteras inmersas.
  ! Se lee el n\'umero de archivos que definen las fronteras inmersas y
  ! sus nombres. Para cada uno se definen las propiedades indicadas
  !
  !*******************************************************************  
  !
  subroutine lectura_cuerpo()
    !
    implicit none
    !
    integer :: num_archivos   ! n'umero de archivos para definir la front. inmersa
    !
    character(len=64) :: archivo ! nombre del archivo a leer
    !r
    integer :: kk
    !
    open(82, file='front_inmersa.dat')
    !
    read(82,*) num_archivos
    !
    do kk=1, num_archivos
       !
       read(82,*) archivo
       !
       call define_cuerpo(archivo)
       !
    end do
    !
    close(82)
    !
  end subroutine lectura_cuerpo
  !
  !*******************************************************************
  !
  ! definir_cuerpo
  !
  ! Subrutina que abre el archivo de propiedades que definen los cuerpos
  ! para las fronteras inmersas. Se leen los 'indices ii, jj y las
  ! propiedades gamma_momen, gamma_energ, fuente_lin_t, fuente_con_t
  !
  !*******************************************************************
  !
  subroutine define_cuerpo(archivoo)
    !
    implicit none
    !
    character(len=64), intent(in) :: archivoo ! nombre del archivo de configuraci\'on
    !
    integer   :: num_puntos    ! n\'umero de puntos de la frontera inmersa
    integer   :: nn            ! enteros para iteraciones
    !
    ! Variables auxiliares para escribir las propiedades
    !
    real(kind=DBL) :: gam_mom, gam_ene, fue_con, fue_lin
    !
    integer   :: ii,jj
    !
    !------------------------
    !
    ! Se abre el archivo de configuraci'on
    !
    open(81, file=trim(archivoo))
    !
    ! Se lee el n'umero de puntos en el archivo
    !
    read(81,*) num_puntos
    !
    ! Se leen los puntos de la frontera inmersa
    !
    do nn = 1, num_puntos
       !
       read(81,*) ii, jj, gam_mom, gam_ene, fue_con, fue_lin
       gamma_momen(ii,jj) = real(gam_mom,DBL)
       gamma_energ(ii,jj) = real(gam_ene,DBL)
       !
       fuente_con_t(ii,jj) = real(fue_con,DBL)
       fuente_lin_t(ii,jj) = real(fue_lin,DBL)
       !
    end do
    !
    close(81)
    !
  end subroutine define_cuerpo
  !
  !*******************************************************************
  !
  ! definir_cuerpo
  !
  ! Subrutina que define la regi\'on de la malla que est\'a dentro de
  ! un cuerpo s\'olido o no. Asigna los valores correspondientes a las
  ! constantes de difusi\'on gamma.
  !
  !*******************************************************************
  !
  subroutine definir_cuerpo(gamma_momeno, gamma_enero, opcion)
    !
    ! import xp, yp, DBL
    ! import en_region_tipo1
    ! import recta_m1xb1, recta_m2xb2
    implicit none
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(inout) :: gamma_momeno, gamma_enero
    character(5), intent(in)                            :: opcion
    !
    integer         :: ii, jj
    real(kind=DBL)  :: xv, yu, hh, hx, hy
    real(kind=DBL)  :: m1, m2 ! Pend. tri\'angulo
    ! procedure(fun    :: f1, f2
    integer         :: xmax, xmin, ymax, ymin !indices de los extremos de los puntos interiores
    !
    ! Mensaje
    !
    print*, " Se define una frontera inmersa con la opcion: ",opcion
    !
    !xmin = mi+1
    !ymin = nj+1
    !xmax = 0
    !ymax = 0
    !
    if( opcion == 'horno' )then
       !
       ! Rectangulo s'olido  con centro en xv, yv y lado hh
       !
       yu = 3.0_DBL
       xv = 1.5_DBL
       hx = 2.8_DBL
       hy = 5.8_DBL
       !
       do jj = 1, nj+1
          !
          if( yu-(hy/2.0_DBL) <= yp(jj) .and. yp(jj) <= yu+(hy/2.0_DBL) )then
                !
             do ii = 1, mi+1
                if( xv-(hx/2.0_DBL) <= xp(ii) .and. xp(ii) <= xv+(hx/2.0_DBL) )then
                   !
                   gamma_momeno(ii,jj) = 10.0e20_DBL
                   !
                   gamma_energ(ii,jj) = 10.0e-10_DBL
                   !
                   fuente_lin_t(ii,jj) =-10.0e40_DBL
                   fuente_con_t(ii,jj) = 10.0e40_DBL*10.0e-5_DBL
                   !
                end if
                !
             end do
             !
          end if
          !
       end do
       !
       ! Fuente de calor
       !
       yu = 1.50_DBL
       xv = 0.05_DBL
       hx = 0.1_DBL
       hy = 1.0_DBL
       !
       do jj = 1, nj+1
          if( yu-(hy/2.0_DBL) <= yp(jj) .and. yp(jj) < yu+(hy/2.0_DBL) )then
                !
             do ii = 1, mi+1
                if( xv-(hx/2.0_DBL) <= xp(ii) .and. xp(ii) <= xv+(hx/2.0_DBL) )then
                   gamma_momeno(ii,jj) = sqrt(7.0_DBL/10.0e7_DBL)*5.0_DBL
                   !gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   fuente_lin_t(ii,jj) = 0.0_DBL
                   fuente_con_t(ii,jj) = 1.0e-2_DBL
                   !
                end if
                !
             end do
             !
          end if
          !
       end do
       !
       ! Cavidad (hornito)
       !
       yu = 4.0_DBL
       xv = 2.4_DBL
       hx = 0.9_DBL
       hy = 1.0_DBL
       !
       do jj = 1, nj+1
          if( yu-(hy/2.0_DBL) <= yp(jj) .and. yp(jj) <= yu+(hy/2.0_DBL) )then
                !
             do ii = 1, mi+1
                if( xv-(hx/2.0_DBL) <= xp(ii) .and. xp(ii) <= xv+(hx/2.0_DBL) )then
                   !
                   gamma_momeno(ii,jj) = sqrt(7.0_DBL/2.94e8_DBL)
                   gamma_energ(ii,jj)  = sqrt(1.0_DBL/(7.0_DBL*2.94e8_DBL))
                   ! gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   fuente_lin_t(ii,jj) = 0.0_DBL
                   fuente_con_t(ii,jj) = 0.0_DBL
                   !
                end if
                !
             end do
             !
          end if
          !
       end do
       !
       ! Pared de la cavidad (hornito)
       !
       ! yu = 8.5_DBL
       xv = 2.875_DBL
       hx = 0.05_DBL
       hy = 1.0_DBL
       !
       do jj = 1, nj+1
          !
          if( yu-(hy/2.0_DBL) <= yp(jj) .and. yp(jj) <= yu+(hy/2.0_DBL) )then
             !
             do ii = 1, mi+1
                !
                if( xv-(hx/2.0_DBL) <= xp(ii) .and. xp(ii) <= xv+(hx/2.0_DBL) )then
                   !
                   !
                   gamma_momeno(ii,jj) = 10.0e20_DBL
                   gamma_energ(ii,jj)  = 1.0e-3_DBL
                   !
                   fuente_lin_t(ii,jj) = 0.0_DBL
                   fuente_con_t(ii,jj) = 0.0_DBL
                   !
                end if
                !
             end do
             !
          end if
          !
       end do
    !
    else if( opcion == 'rectn' )then
       !
       ! Rectangulo con centro en xv, yv y lado hh
       !
       yu = 5.0_DBL
       xv = 2.5_DBL
       hx = 4.8_DBL
       hy = 9.8_DBL
       !
       do jj = 1, nj+1
          if( yu-(hy/2.0_DBL) <= yp(jj) .and. yp(jj) <= yu+(hy/2.0_DBL) )then
                !
             do ii = 1, mi+1
                if( xv-(hx/2.0_DBL) <= xp(ii) .and. xp(ii) <= xv+(hx/2.0_DBL) )then
                   !gamma_momeno(ii,jj) = 2.0_DBL
                   gamma_momeno(ii,jj) = 10.0e20_DBL
                   !
		   gamma_energ(ii,jj) = 10.0e-10_DBL
		   !
                   fuente_lin_t(ii,jj) =-10.0e40_DBL
                   fuente_con_t(ii,jj) = 10.0e40_DBL*10.0e-3_DBL
                   !
                end if
                !
             end do
             !
          end if
          !
       end do
       !
       yu = 5.0_DBL
       xv = 0.05_DBL
       hx = 0.1_DBL
       hy = 1.0_DBL
       !
       do jj = 1, nj+1
          if( yu-(hy/2.0_DBL) <= yp(jj) .and. yp(jj) <= yu+(hy/2.0_DBL) )then
                !
             do ii = 1, mi+1
                if( xv-(hx/2.0_DBL) <= xp(ii) .and. xp(ii) <= xv+(hx/2.0_DBL) )then
                   gamma_momeno(ii,jj) = sqrt(7.0_DBL/10.0e7_DBL)*5.0_DBL
                   !gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   fuente_lin_t(ii,jj) =-10.0e40_DBL
                   fuente_con_t(ii,jj) = 10.0e40_DBL*1.0_DBL
                   !
                end if
                !
             end do
             !
          end if
          !
       end do
       !
    else if( opcion == 'trian' )then
       !
       ! Tri\'angulo con v\'ertice en p(xv,yv) y altura hh
       ! definido usando las funciones recta_mxb1 y recta_mxb2
       !
       yu = 0.3_DBL
       hh = 0.3_DBL
       do jj = 1, nj+1
          if( yu-hh <= yp(jj) .and. yp(jj) <= yu )then
             do ii = 1, mi+1
                ! xx = yp(jj) Se utiliza como una regi\'on de tipo II
                ! yy = xp(ii)
                if( en_region_tipo1(yp(jj),xp(ii),yu-hh,yu,recta_mxb1,recta_mxb2) )then
                   ! print*, "DEBUG: dentro", xp(ii), yp(jj)
                   !gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   fuente_lin_u(ii,jj) =-10.0e50_DBL
                   fuente_con_u(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   !
                   fuente_lin_v(ii,jj) =-10.0e50_DBL
                   fuente_con_v(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   !
                end if
             end do
          end if
       end do
       !
    else if( opcion == 'cuadr' )then
       !
       ! Cuadrado con centro en xv, yv y lado hh
       !
       yu = 5.0_DBL
       xv = 5.0_DBL
       hh = 9.8_DBL
       !
       do jj = 1, nj+1
          if( yu-hh/2.0_DBL <= yp(jj) .and. yp(jj) <= yu+hh/2.0_DBL )then
                !
                !if(jj < ymin) .or. (jj > ymax) then
                !        ymin = jj-1
                !        ymax = jj
                !end if
                !
             do ii = 1, mi+1
                ! xx = yp(jj) Se utiliza como una regi\'on de tipo II
                ! yy = xp(ii)
                if( xv-hh/2.0_DBL <= xp(ii) .and. xp(ii) <= xv+hh/2.0_DBL )then
                   !gamma_momeno(ii,jj) = 2.0_DBL
                   gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   !if(ii < xmin .or. ii > xmax) then
                   !     xmin = ii-1
                   !     xmax = ii
                   !end if
                   !
                   ! u(ii,jj)            = 10.0_DBL
                   ! v(ii,jj)            = 5.0_DBL
                   !                   
                   ! u_ant(ii,jj)        = 10.0_DBL
                   ! v_ant(ii,jj)        = 5.0_DBL
                   !
                   !fuente_lin_u(ii,jj) =-10.0e50_DBL
                   !fuente_con_u(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   !
                   !fuente_lin_v(ii,jj) =-10.0e50_DBL
                   !fuente_con_v(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   !
                   !
                   !fuente_lin_t(ii,jj) =-10.0e40_DBL
                   !fuente_con_t(ii,jj) = 10.0e40_DBL*1.0_DBL
                   !
                   !print*, "DEBUG: dentro", ii,jj, u(ii,jj)
                end if
                !
             end do
             !
          end if
          !
       end do
       yu = 9.0_DBL
       xv = 0.05_DBL
       hx = 0.1_DBL
       hy = 2.0_DBL
       !
       do jj = 1, nj+1
          if( yu-(hy/2.0_DBL) <= yp(jj) .and. yp(jj) <= yu+(hy/2.0_DBL) )then
                !
             do ii = 1, mi+1
                if( xv-(hx/2.0_DBL) <= xp(ii) .and. xp(ii) <= xv+(hx/2.0_DBL) )then
                   gamma_momeno(ii,jj) = sqrt(7.0_DBL/10.0e7_DBL)*5.0_DBL
                   !gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   fuente_lin_t(ii,jj) =-10.0e40_DBL
                   fuente_con_t(ii,jj) = 10.0e40_DBL*1.0_DBL
                   !
                end if
                !
             end do
             !
          end if
          !
       end do
       !
       !
    else if( opcion == 'chime' )then
       !
       ! varillas inmersas para chimenea solar
       !
       yu = 0.5_DBL
       xv = 8.5_DBL
       hh = 0.3_DBL
       !
       do jj = 1, nj+1
          if( yu-hh/2.0_DBL <= yp(jj) .and. yp(jj) <= yu+hh/2.0_DBL )then
             do ii = 1, mi+1
                ! xx = yp(jj) Se utiliza como una regi\'on de tipo II
                ! yy = xp(ii)
                if( xv-hh/2.0_DBL <= xp(ii) .and. xp(ii) <= xv+hh/2.0_DBL )then
                   ! gamma_momeno(ii,jj) = 10.0e6_DBL
                   !
                   ! u(ii,jj)            = 10.0_DBL
                   ! v(ii,jj)            = 5.0_DBL
                   !                   
                   ! u_ant(ii,jj)        = 10.0_DBL
                   ! v_ant(ii,jj)        = 5.0_DBL
                   !
                   fuente_lin_u(ii,jj) =-10.0e50_DBL
                   fuente_con_u(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   au(ii,jj)           = 10.0e40_DBL
                   !
                   fuente_lin_v(ii,jj) =-10.0e50_DBL
                   fuente_con_v(ii,jj) = 10.0e50_DBL*10.0e-12_DBL
                   av(ii,jj)           = 10.0e40_DBL
                   !
                   fuente_lin_t(ii,jj) =-10e50_DBL
                   fuente_con_t(ii,jj) = 10e50_DBL*1.0_DBL
                   !print*, "DEBUG: dentro", ii,jj, u(ii,jj)
                end if
                !
             end do
             !
          end if
          !
       end do
       !
    else if ( opcion == 'tubov' )then
       !
       ! rectangulo
       !
       yu = 3.0_DBL
       hh = 2.0_DBL
       !
       do jj = 1, nj+1
          if( yu-(hh/2.0_DBL) <= yp(jj) .and. yp(jj) <= yu+(hh/2.0_DBL))then
             !
             do ii = 1, mi+1
                !
                gamma_momeno(ii,jj) = 4.0_DBL*10.0e-4_DBL
                !gamma_momeno(ii,jj) = 10.0e6_DBL
                !
             end do
             !
          end if
          !
       end do
       !
    end if ! Selecci\'on del caso del cuerpo inmerso
    !
  end subroutine definir_cuerpo
  !
  !
  !*******************************************************************
  !
  ! cond_front_inmersa
  !
  !*******************************************************************
  !
  subroutine cond_front_inmersa(AI_o,AC_o,AD_o,Rx_o, opcion)
    !
    implicit none
    !
    real(kind=DBL), dimension(mi+1,nj+1), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    character(5), intent(in)                          :: opcion
    !
    integer         :: ii, jj
    real(kind=DBL)  :: xv, yu, hh
    real(kind=DBL)  :: m1, m2 ! Pend. tri\'angulo
    !
    !
    if( opcion == 'cuadr' )then
       !
       ! Cuadrado con centro en xv, yu y lado hh
       !
       yu = 0.5_DBL
       xv = 6.0_DBL
       hh = 0.2_DBL
       !
       do jj = 1, nj+1
          if( yu-hh/2.0_DBL <= yp(jj) .and. yp(jj) <= yu+hh/2.0_DBL )then
             do ii = 1, mi+1
                if( xv-hh/2.0_DBL <= xp(ii) .and. xp(ii) <= xv+hh/2.0_DBL )then
                   AI_o(ii,jj) = 0.0_DBL
                   AD_o(ii,jj) = 0.0_DBL
                   AC_o(ii,jj) = 1.0_DBL
                   Rx_o(ii,jj) = 0.0_DBL
                   !print*, "DEBUG: dentro", ii,jj, u(ii,jj)
                end if
                !
             end do
             !
          end if
          !
       end do
       !
!     else if( opcion == 'trian' )then
!        !
!        ! Tri\'angulo con v\'ertice en p(xv,yu) y altura hh
!        ! definido usando las funciones recta_mxb1 y recta_mxb2
!        !
!        yv = 0.3_DBL
!        hh = 0.3_DBL
!        do jj = 1, nj
!           if( yu-hh <= yv(jj) .and. yv(jj) <= yu )then
!              do ii = 1, mi
!                 if( en_region_tipo1(yv(jj),xu(ii),yu-hh,yu,recta_mxb1,recta_mxb2) )then
!                    ! print*, "DEBUG: dentro", xp(ii), yp(jj)
!                    au_o(ii,jj) = 1.e40_DBL
!                    !
!                 end if
!              end do
!           end if
!        end do
       !
    end if ! Selecci\'on del caso del cuerpo inmerso
    !
  end subroutine cond_front_inmersa
  !
  !*******************************************************************
  !
  ! en_region_tipo1
  !
  ! Funci\'on para decidir si un punto est\'a en una regi\'on en el
  ! plano entre dos valores fijos a <= x <= b y dos funciones
  ! f1(x) <= y <= f2(x)
  !
  !*******************************************************************
  !
  function en_region_tipo1(xx,yy,aa,bb,f1,f2)
    !
    ! import DBL
    ! implicit none
    !
    logical                    :: en_region_tipo1 ! true para dentro de la region
    real(kind=DBL), intent(in) :: xx, yy          ! coordenadas del punto de prueba
    real(kind=DBL), intent(in) :: aa, bb          ! extremos para x
    !
    procedure(funcion1d)       :: f1, f2          ! extremos para y
    !
    en_region_tipo1=.false.
    !
    if( aa <= xx .and. xx <= bb )then
       !
       if( f1(xx) <= yy .and. yy<= f2(xx) )then
          !
          en_region_tipo1 = .true.
          !
       end if
       !
    end if
    !    
  end function en_region_tipo1
  !*******************************************************************
  !
  ! funciones para definir cuerpos inmersos
  !
  !*******************************************************************
  !
  ! Recta de pendiente mm y ordenada al origen bb
  !
  real(kind=DBL) function recta_mxb1(xx)
    !
    ! import DBL
    ! implicit none
    !
    real(kind=DBL), intent(in) :: xx
    real(kind=DBL)             :: mm, bb
    !
    mm = 1.0_DBL/2.0_DBL
    bb =-( 0.3_DBL - mm * 6.0_DBL ) / mm
    recta_mxb1 = xx/mm + bb
    !
  end function recta_mxb1
  !
  ! Recta de pendiente mm y ordenada al origen bb
  !
  real(kind=DBL) function recta_mxb2(xx)
    !
    ! import DBL
    implicit none
    !
    real(kind=DBL), intent(in) :: xx
    real(kind=DBL)             :: mm, bb
    !
    mm =-1.0_DBL/2.0_DBL
    bb =-( 0.3_DBL - mm * 6.0_DBL ) / mm
    recta_mxb2 = xx/mm + bb
    !
  end function recta_mxb2
  !  
end module frontera_inmersa
