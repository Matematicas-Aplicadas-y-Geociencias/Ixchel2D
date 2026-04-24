!
!
! M\'odulo de herramientas para la lectura de condiciones de frontera del programa SIMPLE 2D
!
! Este m\'odulo contiene herramientas que se utilizan para leer las condiciones de frontera
! en un dominio rectangular
!
MODULE cond_frontera
  !
  use malla, only : mi, nj, DBL
  !
  implicit none
  !
  ! Formato de escritura
  !
  character(len=16), parameter :: form36="(3A,1I2,2A)"
  character(len=32), parameter :: form37="(3A,F8.3,1A,F8.3,2A)"
  !
  ! Nombres de los archivos para condiciones de frontera
  !
  character(len=18) :: entrada_front_u
  character(len=18) :: entrada_front_v
  character(len=18) :: entrada_front_t
  !
  type tipo_cond_front
     !
     ! Estructura de datos para determinar condiciones de frontera mediante archivos
     ! Se hace una estructura r\'igida capaz de manejar 15 divisiones por lado. Si se
     ! requieren m\'as es posible aumentar las dimensiones de los arreglos o crear
     ! un puntero asignable (allocatable)
     ! 
     character(len=1)                :: lado_front  ! Se usan 4 lados en 2D: a,b,c,d
     integer                         :: ndivis      ! n'umero de divisiones
     character(len=4), dimension(15) :: tipo_condi  ! tipo de condici\'on de frontera
     real(kind=DBL), dimension(15)   :: valor_cond  ! valor de la condici\'on de frontera
     character(len=4), dimension(15) :: func_condi  ! Funci'on espacio-tiempo de la c.f. 
     integer, dimension(14)          :: indice_div  ! indices iniciales de las divisiones
     !
  end type tipo_cond_front
  !
contains
  !
  !-----------------------------------------------------------------------------
  !
  ! inicializa_cond_front
  !
  ! Esta subrutina inicializa la estructura de condiciones de frontera
  !
  subroutine inicializa_cond_front(cond_front_uu)
    !
    implicit none
    !
    class( tipo_cond_front), intent(inout) :: cond_front_uu
    !
    cond_front_uu % lado_front = 'z'
    cond_front_uu % ndivis     =  0
    cond_front_uu % tipo_condi = 'desc'
    cond_front_uu % valor_cond =  0.0_DBL
    cond_front_uu % func_condi = 'desc'
    cond_front_uu % indice_div = -14
    !
  end subroutine inicializa_cond_front
  !
  !-----------------------------------------------------------------------------
  !
  ! lectura_cond_frontera
  !
  ! Esta subrutina abre los archivos para leer las condiciones de frontera
  ! y crea los arreglos necesarios para imponerlas durante la construcci\'on
  ! de las matrices
  !
  subroutine lectura_cond_frontera(entrada_front_uuo,&
       &cond_front_uua,&
       &cond_front_uub,&
       &cond_front_uuc,&
       &cond_front_uud,&
       &xx,yy,&
       &mm,nn &
       &)
    !
    implicit none
    !
    character(len=18), intent(in) :: entrada_front_uuo ! archivo de entrada
    !
    class(tipo_cond_front), intent(inout) :: cond_front_uua    ! condici'on lado a
    class(tipo_cond_front), intent(inout) :: cond_front_uub    ! condici'on lado b
    class(tipo_cond_front), intent(inout) :: cond_front_uuc    ! condici'on lado c
    class(tipo_cond_front), intent(inout) :: cond_front_uud    ! condici'on lado d
    !
    real(kind=DBL), dimension(mm), intent(in) :: xx
    real(kind=DBL), dimension(nn), intent(in) :: yy
    integer, intent(in)                       :: mm,nn
    !
    character(len=1) :: lado
    character(len=4) :: variable
    character(len=4) :: tipo_condicion
    character(len=4) :: funcion
    !
    real(kind=DBL)   :: x0, x1, valor
    !
    integer          :: divisiones
    integer          :: indice
    integer          :: ii, jj, kk
    !
    ! *********************************************************
    !
    ! Apertura de los archivos y lectura de los arreglos
    !
    write(*,*) " "
    write(*,*) "Ixchel2D: Inicio de lectura de condiciones de frontera"
    write(*,*) " "    
    open(unit=111,file=entrada_front_uuo)
    !
    ! Se esperan 4 lados en el rect\'angulo en el archivo
    ! de condiciones de frontera
    !
    lectura: do ii = 1, 4
       !
       read(111,*) variable, lado, divisiones
       !
       write(*,form36) ' Cond. de frontera para ', variable, &
            &' con ', divisiones,' division(es) en lado ', lado
       !
       sel_lado: select case( lado )
          !
       case( 'a' )
          !
          cond_front_uua % indice_div(1) = 2 ! indice inicial para bucles de frontera
          !
          divisiones_xa: do jj = 1, divisiones
             !
             read(111,*) x0, x1, tipo_condicion, valor, funcion
             !
             write(*,form37) "    Tipo ", tipo_condicion, &
                  &" entre ", x0, " y ", x1," en lado ", lado
             !
             cond_front_uua % lado_front     = 'a'
             cond_front_uua % ndivis         = divisiones
             cond_front_uua % tipo_condi(jj) = tipo_condicion
             cond_front_uua % valor_cond(jj) = valor
             cond_front_uua % func_condi(jj) = funcion
             !
             indices_conda: do kk = 1, nn
                !
                if( yy(kk) > x1 ) exit indices_conda
                !
             end do indices_conda
             !
             cond_front_uua % indice_div(jj+1) = kk
             !
          end do divisiones_xa
          !
          cond_front_uua % indice_div(divisiones+1) = nn-1 ! indice final bucles frontera
          !
       case( 'b' )
          !
          cond_front_uub % indice_div(1) = 2 ! indice inicial para bucles de frontera
          !
          divisiones_yb: do jj = 1, divisiones
             !
             read(111,*) x0, x1, tipo_condicion, valor, funcion
             !
             write(*,form37) "    Tipo ", tipo_condicion, &
                  &" entre ", x0, " y ", x1," en lado ", lado
             !
             cond_front_uub % lado_front     = 'b'
             cond_front_uub % ndivis         = divisiones
             cond_front_uub % tipo_condi(jj) = tipo_condicion
             cond_front_uub % valor_cond(jj) = valor
             cond_front_uub % func_condi(jj) = funcion
             !
             indices_condb: do kk = 1, mm
                !
                if( xx(kk) > x1 ) exit indices_condb
                !
             end do indices_condb
             !
             cond_front_uub % indice_div(jj+1) = kk
             !
          end do divisiones_yb
          !
          cond_front_uub % indice_div(divisiones+1) = mm-1 ! indice final bucles frontera
          !
       case( 'c' )
          !
          cond_front_uuc % indice_div(1) = 2 ! indice inicial para bucles de frontera
          !
          divisiones_xc: do jj = 1, divisiones
             !
             read(111,*) x0, x1, tipo_condicion, valor, funcion
             !
             write(*,form37) "    Tipo ", tipo_condicion, &
                  &" entre ", x0, " y ", x1," en lado ", lado
             !
             cond_front_uuc % lado_front     = 'c'
             cond_front_uuc % ndivis         = divisiones
             cond_front_uuc % tipo_condi(jj) = tipo_condicion
             cond_front_uuc % valor_cond(jj) = valor
             cond_front_uuc % func_condi(jj) = funcion
             !
             indices_condc: do kk = 1, nn
                !
                if( yy(kk) > x1 ) exit indices_condc
                !
             end do indices_condc
             !
             cond_front_uuc % indice_div(jj+1) = kk
             !
          end do divisiones_xc
          !
          cond_front_uuc % indice_div(divisiones+1) = nn-1 ! indice final bucles frontera
          !
       case( 'd' )
          !
          cond_front_uud % indice_div(1) = 2 ! indice inicial para bucles de frontera
          !
          !
          divisiones_yd: do jj = 1, divisiones
             !
             read(111,*) x0, x1, tipo_condicion, valor, funcion
             !
             write(*,form37) "    Tipo ", tipo_condicion, &
                  &" entre ", x0, " y ", x1," en lado ", lado
             !
             cond_front_uud % lado_front     = 'd'
             cond_front_uud % ndivis         = divisiones
             cond_front_uud % tipo_condi(jj) = tipo_condicion
             cond_front_uud % valor_cond(jj) = valor
             cond_front_uud % func_condi(jj) = funcion
             !
             indices_condd: do kk = 1, mm
                !
                if( xx(kk) > x1 ) exit indices_condd
                !
             end do indices_condd
             !
             cond_front_uud % indice_div(jj+1) = kk
             !
          end do divisiones_yd
          !
          cond_front_uud % indice_div(divisiones+1) = mm-1 ! indice final bucles frontera
          !
       end select sel_lado
          !
    end do lectura
    !
    close(unit=111)
    !
    write(*,*) " "
    write(*,*) "Ixchel2D: Fin de lectura de condiciones de frontera"
    write(*,*) "-------------------------------------------------------"  
  end subroutine lectura_cond_frontera
  !
  !***************************************************************************
  !
  ! impone_cond_frontera
  !
  ! Esta subrutina recibe estructuras con condiciones de frontera y las impone
  ! en los coeficientes de las matrices para la direcci\'on x
  !
  subroutine impone_cond_frontera(cond_front_uu,&
       & AI_o,AC_o,AD_o,Rx_o, &
       & mm,nn,               &
       & kk,ll,               &
       & tiem,                &
       & au_o )
    !
    !$acc routine seq
    !
    implicit none
    !
    class( tipo_cond_front ), intent(in)          :: cond_front_uu
    !
    real(kind=DBL), dimension(mm,nn), intent(out) :: AI_o, AC_o, AD_o, Rx_o
    real(kind=DBL), intent(in)                    :: tiem
    !
    real(kind=DBL), dimension(kk,ll), intent(out), optional :: au_o
    !
    integer, intent(in)                           :: mm, nn, kk, ll
    !
    integer :: ldiv, ii, jj
    !
    !-------------------------------
    lado: select case( cond_front_uu % lado_front )
       !
       !-------------------------------
       !
       ! lado a
       !
    case( 'a' )
       !
       bucle_segmento_au: do ldiv = 1, cond_front_uu % ndivis
          !
          if( cond_front_uu % tipo_condi(ldiv) == 'diri' )then
             !
             do jj = cond_front_uu % indice_div(ldiv), cond_front_uu % indice_div(ldiv+1)
                !
                AI_o(1,jj) = 0.0_DBL
                AC_o(1,jj) = 1.0_DBL
                AD_o(1,jj) = 0.0_DBL
                Rx_o(1,jj) = cond_front_uu % valor_cond(ldiv) * &
                     dirichlet_ini( tiem, cond_front_uu % func_condi(ldiv) )
                if ( present(au_o) ) au_o(1,jj) = 1.e40_DBL
                ! print*, "DEBUG: Dirichlet en a"
             end do
             !
          else if( cond_front_uu % tipo_condi(ldiv) == 'neum' )then
             !
             ! print*, "DEBUG: ", cond_front_uu % valor_cond(ldiv)
             do jj = cond_front_uu % indice_div(ldiv), cond_front_uu % indice_div(ldiv+1)
                !
                AI_o(1,jj) = 0.0_DBL
                AC_o(1,jj) =-1.0_DBL
                AD_o(1,jj) = 1.0_DBL
                Rx_o(1,jj) = cond_front_uu % valor_cond(ldiv)
                if ( present(au_o) ) au_o(1,jj) = 1.e40_DBL
                ! print*, "DEBUG: neumann en a"
             end do
             !
          end if
          !
       end do bucle_segmento_au
       !
       !-------------------------------
       !
       ! lado b
       !
    case( 'b' )
       !
       bucle_segmento_bu: do ldiv = 1, cond_front_uu % ndivis
          !
          if( cond_front_uu % tipo_condi(ldiv) == 'diri' )then
             !
             do jj = cond_front_uu % indice_div(ldiv), cond_front_uu % indice_div(ldiv+1)
                !
                AI_o(1,jj) = 0.0_DBL
                AC_o(1,jj) = 1.0_DBL
                AD_o(1,jj) = 0.0_DBL
                Rx_o(1,jj) = cond_front_uu % valor_cond(ldiv) * &
                     dirichlet_ini( tiem, cond_front_uu % func_condi(ldiv) )
                if ( present(au_o) ) au_o(jj,1) = 1.e40_DBL
                ! print*, "DEBUG: Dirichlet en a"
             end do
             !
          else if( cond_front_uu % tipo_condi(ldiv) == 'neum' )then
             !
             ! print*, "DEBUG: ", cond_front_uu % valor_cond(ldiv)
             do jj = cond_front_uu % indice_div(ldiv), cond_front_uu % indice_div(ldiv+1)
                !
                AI_o(1,jj) = 0.0_DBL
                AC_o(1,jj) =-1.0_DBL
                AD_o(1,jj) = 1.0_DBL
                Rx_o(1,jj) = cond_front_uu % valor_cond(ldiv) * &
                     dirichlet_ini( tiem, cond_front_uu % func_condi(ldiv) )
                if ( present(au_o) ) au_o(jj,1) = 1.e40_DBL
                ! print*, "DEBUG: neumann en a"
             end do
             !
          end if
          !
       end do bucle_segmento_bu
       !
       !-------------------------------
       !
       ! lado c
       !
    case( 'c' )
       ! 
       ! $acc parallel loop vector
       bucle_segmento_cu: do ldiv = 1, cond_front_uu % ndivis
          if( cond_front_uu % tipo_condi(ldiv) == 'diri' )then
             !
             do jj = cond_front_uu % indice_div(ldiv), cond_front_uu % indice_div(ldiv+1)
                !
                AI_o(kk,jj) = 0.0_DBL
                AC_o(kk,jj) = 1.0_DBL
                AD_o(kk,jj) = 0.0_DBL
                Rx_o(kk,jj) = cond_front_uu % valor_cond(ldiv) * &
                     dirichlet_ini( tiem, cond_front_uu % func_condi(ldiv) )
                if ( present(au_o) ) au_o(kk,jj) = 1.e40_DBL
                ! print*, "DEBUG: Dirichlet en a"
             end do
             !
          else if( cond_front_uu % tipo_condi(ldiv) == 'neum' )then
             !
             do jj = cond_front_uu % indice_div(ldiv), cond_front_uu % indice_div(ldiv+1)
                !
                AI_o(kk,jj) =-1.0_DBL
                AC_o(kk,jj) = 1.0_DBL
                AD_o(kk,jj) = 0.0_DBL
                Rx_o(kk,jj) = cond_front_uu % valor_cond(ldiv)
                if ( present(au_o) ) au_o(kk,jj) = 1.e40_DBL
                ! print*, "DEBUG: neumann en a"
             end do
             !
          end if
          !
       end do bucle_segmento_cu
       !
       !-------------------------------
       !
       ! lado d
       !
    case( 'd' )
       ! 
       ! $acc parallel loop vector
       bucle_segmento_du: do ldiv = 1, cond_front_uu % ndivis
          if( cond_front_uu % tipo_condi(ldiv) == 'diri' )then
             !
             do jj = cond_front_uu % indice_div(ldiv), cond_front_uu % indice_div(ldiv+1)
                !
                AI_o(ll,jj) = 0.0_DBL
                AC_o(ll,jj) = 1.0_DBL
                AD_o(ll,jj) = 0.0_DBL
                Rx_o(ll,jj) = cond_front_uu % valor_cond(ldiv) * &
                     dirichlet_ini( tiem, cond_front_uu % func_condi(ldiv) )
                if ( present(au_o) ) au_o(jj,ll) = 1.e40_DBL
                ! print*, "DEBUG: Dirichlet en a"
             end do
             !
          else if( cond_front_uu % tipo_condi(ldiv) == 'neum' )then
             !
             do jj = cond_front_uu % indice_div(ldiv), cond_front_uu % indice_div(ldiv+1)
                !
                AI_o(ll,jj) =-1.0_DBL
                AC_o(ll,jj) = 1.0_DBL
                AD_o(ll,jj) = 0.0_DBL
                Rx_o(ll,jj) = cond_front_uu % valor_cond(ldiv)
                if ( present(au_o) ) au_o(jj,ll) = 1.e40_DBL
                ! print*, "DEBUG: neumann en a"
             end do
             !
          end if
          !
       end do bucle_segmento_du
       !
    end select lado
    !
  end subroutine impone_cond_frontera
  !
  !*************************************************************
  !
  ! funciones para definir condiciones dependientes del tiempo
  !
  !*************************************************************
  !
  ! funci'on del tiempo para arrancar simulaciones 
  !
  real(kind=DBL) function dirichlet_ini(xx,nombre)
    !
    real(kind=DBL), intent(in)   :: xx
    !
    character(len=4), intent(in) :: nombre
    !
    sel_funcion: select case( nombre )
       !
    case( 'cons' )
       !
       dirichlet_ini = 1.0_DBL
       !
    case( 'tanh' )
       !
       dirichlet_ini = dtanh(2.0_DBL*xx/3.0_DBL)
       !
    case( 'ramp' )
       !
       !       print*,"DEBUG: Rampa ", xx
       if( xx .le. 1.0_DBL )then
          !
          dirichlet_ini = xx
          !
       else
          !
          dirichlet_ini = 1.0_DBL
          !
       end if
       !
    end select sel_funcion
    !
  end function dirichlet_ini
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
end MODULE cond_frontera
