!
! M\'odulo de postproceso
!
! Este m\'odulo contiene las subrutinas que calculan cantidades 
! de salida y archivos de postproceso
!
! Autor: J.C. Cajas, K. Figueroa
!
module postproceso
  !
  use malla, only : mi, nj, DBL
  use malla, only : deltaxp, xp
  use malla, only : deltayp, yp
  use malla, only : deltaxu, xu
  use malla, only : deltayv, yv
  !
  implicit none
  !
  ! Formatos de salida
  !
  character(len=12), parameter :: form44="(33E25.15E3)"
  !
  type tipo_promedio_perfil
     !
     ! Estructura de datos para calcular cantidades promedio en cortes horizontales
     ! o verticales. Se crea una estructura r\'igida con capacidad de almacenar 20
     ! cortes. Si se necesitan m'as hay que modificar la dimensi'on del arreglo
     ! posicion del tipo
     ! 
     character(len=5)                :: orienta     ! Orientacion (horiz/verti)
     integer                         :: nposi       ! n'umero de posiciones
     real(kind=DBL), dimension(16,2) :: valor_prom  ! posicion,cantidad promediada
     integer, dimension(16)          :: indi_posi   ! indices entero de la posicion
     !
  end type tipo_promedio_perfil
  !
  ! Se declaran variables para calcular promedios en perfiles por defecto
  !
  type( tipo_promedio_perfil )      :: temp_promedio_perfilh, temp_promedio_perfilv
  type( tipo_promedio_perfil )      :: velu_promedio_perfilh, velu_promedio_perfilv
  type( tipo_promedio_perfil )      :: velv_promedio_perfilh, velv_promedio_perfilv
  !
contains
  !-----------------------------------------------------------------------------
  !
  ! inicializa_promedio_perfil
  !
  ! Esta subrutina inicializa la estructura de datos promedio_perfil
  !
  subroutine inicializa_promedio_perfil(prom_perf)
    !
    implicit none
    !
    class( tipo_promedio_perfil ), intent(inout) :: prom_perf
    !
    prom_perf % orienta    = 'aaaaz'
    prom_perf % nposi      = -2
    prom_perf % indi_posi  = -2
    prom_perf % valor_prom(:,:) = -444.0_DBL
    !
  end subroutine inicializa_promedio_perfil
  !
  !************************************************************
  ! lectura_archivo_prom
  !
  ! subrutina que lee el archivo que contiene el numero de integrales
  ! y las alturas a las que se desea hacerlas.
  !
  !************************************************************
  !
  subroutine lectura_archivo_prom()
    !
    implicit none
    !
    character(len=64) :: alturas ! nombre del archivo a leer
    integer           :: kk, ii
    !
    ! Inicializaci'on de los valores para las estructuras posibles
    !
    call inicializa_promedio_perfil(temp_promedio_perfilh)
    call inicializa_promedio_perfil(temp_promedio_perfilv)
    call inicializa_promedio_perfil(velu_promedio_perfilh)
    call inicializa_promedio_perfil(velu_promedio_perfilv)
    call inicializa_promedio_perfil(velv_promedio_perfilh)
    call inicializa_promedio_perfil(velv_promedio_perfilv)
    !
    temp_promedio_perfilh % orienta = 'horiz'
    temp_promedio_perfilv % orienta = 'verti'
    velu_promedio_perfilh % orienta = 'horiz'
    velu_promedio_perfilv % orienta = 'verti'
    velv_promedio_perfilh % orienta = 'horiz'
    velv_promedio_perfilv % orienta = 'verti'
    !
    ! Se abre el archivo de entrada de los perfiles a promediar, y se leen las
    ! posiciones de los perfiles
    !
    open(77, file='perfil_promedio.dat')
    !
    read(77,*) temp_promedio_perfilh % nposi
    !
    do kk=1, temp_promedio_perfilh % nposi
       !
       read(77,*) temp_promedio_perfilh % valor_prom(kk,1)
       !
    end do
    !
    read(77,*) temp_promedio_perfilv % nposi
    !
    do kk=1, temp_promedio_perfilv % nposi
       !
       read(77,*) temp_promedio_perfilv % valor_prom(kk,1)
       !
    end do
    !
    read(77,*) velu_promedio_perfilh % nposi
    !
    do kk=1, velu_promedio_perfilh % nposi
       !
       read(77,*) velu_promedio_perfilh % valor_prom(kk,1)
       !
    end do
    !
    read(77,*) velu_promedio_perfilv % nposi
    !
    do kk=1, velu_promedio_perfilv % nposi
       !
       read(77,*) velu_promedio_perfilv % valor_prom(kk,1)
       !
    end do
    !
    read(77,*) velv_promedio_perfilh % nposi
    !
    do kk=1, velv_promedio_perfilh % nposi
       !
       read(77,*) velv_promedio_perfilh % valor_prom(kk,1)
       !
    end do
    !
    read(77,*) velv_promedio_perfilv % nposi
    !
    do kk=1, velv_promedio_perfilv % nposi
       !
       read(77,*) velu_promedio_perfilv % valor_prom(kk,1)
       !
    end do
    !
    close(77)
    !
    ! Se determinan los \'indices para las posiciones deseadas
    !
    call determina_indices_horizontal(temp_promedio_perfilh,xp,mi+1)
    call determina_indices_vertical(temp_promedio_perfilv,yp,nj+1)
    call determina_indices_horizontal(velu_promedio_perfilh,xu,mi)
    call determina_indices_vertical(velu_promedio_perfilv,yp,nj+1)
    call determina_indices_horizontal(velv_promedio_perfilh,xp,mi+1)
    call determina_indices_vertical(velv_promedio_perfilv,yv,nj)
    !
  end subroutine lectura_archivo_prom
  !
  !***************************************************************
  ! determina_indices
  !
  ! Subrutina que determina los \'indices m\'as cercanos para las
  ! posiciones deseadas de los perfiles en direcci'on horizontal.
  ! Se usa un arreglo xx con dimensi\'on mi+1 para usarse con xp y xu 
  ! en el caso de xu, est'a sobredimensionado
  !
  !***************************************************************
  subroutine determina_indices_horizontal(prom_perf,xx,nn)
    !
    implicit none
    !
    class( tipo_promedio_perfil ), intent(inout) :: prom_perf
    !
    real(kind=DBL), dimension(mi+1), intent(in)  :: xx
    integer, intent(in) :: nn
    !
    integer :: ii, kk
    !
    ! Se usa la variable kk para recorrer las posiciones de los perfiles
    !
    kk = 1
    !
    do ii = 1, nn
       !
       ! Se comparan las alturas de la malla con las alturas deseadas
       ! que est'an en el primer 'indice del arreglo prom de la estructura.
       ! Se guardan los 'indices de los nodos m'as cercanos a la altura
       ! indicada (cercano por arriba)
       !
       if( prom_perf % valor_prom(kk,1) < xx(ii) .and. kk <= &
            & prom_perf % nposi) then
          prom_perf % indi_posi(kk) = ii
          kk = kk+1
       end if
       !
    end do
    !
  end subroutine determina_indices_horizontal
  !
  !
  !***************************************************************
  ! determina_indices_vertical
  !
  ! Subrutina que determina los \'indices m\'as cercanos para las
  ! posiciones deseadas de los perfiles en direcci'on vertical.
  ! Se usa un arreglo xx con dimensi\'on nj+1 para usarse con yp y yv 
  ! en el caso de yv, est'a sobredimensionado
  !
  !***************************************************************
  subroutine determina_indices_vertical(prom_perf,xx,nn)
    !
    implicit none
    !
    class( tipo_promedio_perfil ), intent(inout) :: prom_perf
    !
    real(kind=DBL), dimension(nj+1), intent(in)  :: xx
    integer, intent(in) :: nn
    !
    integer :: ii, kk
    !
    ! Se usa la variable kk para recorrer las posiciones de los perfiles
    !
    kk = 1
    !
    do ii = 1, nn
       !
       ! Se comparan las alturas de la malla con las alturas deseadas
       ! que est'an en el primer 'indice del arreglo prom de la estructura.
       ! Se guardan los 'indices de los nodos m'as cercanos a la altura
       ! indicada (cercano por arriba)
       !
       if( prom_perf % valor_prom(kk,1) < xx(ii) .and. kk <= &
            & prom_perf % nposi) then
          prom_perf % indi_posi(kk) = ii
          kk = kk+1
       end if
       !
    end do
    !
  end subroutine determina_indices_vertical
  !
  !************************************************************
  ! postpro_promedio
  !
  ! subrutina que calcula el promedio de la temperatura a lo
  ! largo del tiempo. Crea un archivo de slida con los datos
  !
  !************************************************************
  !
  subroutine postpro_promedio(opcion, tiempo, temp_o, file_name)
    !
    use malla, only :  mic, njc
    implicit none
    !
    real(kind=DBL), DIMENSION(mi+1,nj+1), intent(in) :: temp_o
    character(32), intent(in)     :: file_name
    character(5), intent(in)      :: opcion
    real(kind=DBL), intent(in)    :: tiempo
    integer                       :: kk, ii
    !
    if (opcion == 'horiz') then
       !
       do kk=1, temp_promedio_perfilh%nposi
          !
          call promedio_horizontal(temp_promedio_perfilh%valor_prom(kk,2),&
               & temp_o(1:mi+1,temp_promedio_perfilh % indi_posi(kk)))
          !
       end do
       !
    else if (opcion == 'verti' ) then
       do kk=1, temp_promedio_perfilh%nposi
          !
          call promedio_vertical(temp_promedio_perfilv%valor_prom(kk,2),&
               & temp_o(temp_promedio_perfilv % indi_posi(kk),1:nj+1))
          !
       end do
       !
    end if
    !
    open(unit = 76,file=file_name,access='append')
    !
    write(76,form44) tiempo, temp_promedio_perfilh%valor_prom(:,:)
    !
    close(76)
    !
  end subroutine postpro_promedio
  !
  !************************************************************
  ! promedio_horizontal
  !
  ! subrutina que calcula la integral de lineas horizontales
  !
  !************************************************************
  !
  subroutine promedio_horizontal(integral, variable)
    !
    implicit none
    !
    real(kind=DBL), intent(out)                 :: integral
    real(kind=DBL), DIMENSION(mi+1), intent(in) :: variable
    !
    integer :: ii
    !
    integral = 0.0_DBL
    integral = variable(1)*(deltaxp(1)/2.0_DBL)
    !
    do ii = 2, mi
       integral= integral + ( variable(ii) + variable(ii+1) ) * deltaxp(ii) * 0.5_DBL
    end do
    !
    integral= integral + variable(mi+1)*(deltaxp(mi)/2.0_DBL)
    !
  end subroutine promedio_horizontal
  !
  !************************************************************
  ! promedio_horizontal
  !
  ! subrutina que calcula la integral de lineas horizontales
  !
  !************************************************************
  !
  subroutine promedio_vertical(integral, variable)
    !
    implicit none
    !
    real(kind=DBL), intent(out)                 :: integral
    real(kind=DBL), DIMENSION(nj+1), intent(in) :: variable
    !
    integer :: jj
    !
    integral = 0.0_DBL
    integral = variable(1)*(deltayp(1)/2.0_DBL)
    !
    do jj = 2, nj
       !
       integral= integral + ( variable(jj) + variable(jj+1) ) * deltayp(jj) * 0.5_DBL
       !
    end do
    !
    integral= integral + variable(nj+1)*(deltayp(nj)/2.0_DBL)
    !
  end subroutine promedio_vertical
  !
  !************************************************************
  !
  ! postprocesa_parametros
  !
  ! Subrutina de postproceso que muestra los parámetros empleados
  ! y verifica que el directorio esté creado
  !
  !************************************************************
  !
  subroutine postprocesa_parametros(&
       &Ra,&
       &Pr,&
       &dt,&
       &itermax,&
       &paq_itera,&
       &Ri_1,&
       &rel_pres,&
       &rel_vel,&
       &rel_ener,&
       &conv_u,&
       &conv_t,&
       &conv_p,&
       &conv_resi,&
       &conv_paso,&
       &iter_simple_max,&
       &iter_ecuaci_max,&
       &entrada_u,&
       &entrada_v,&
       &entrada_tp,&
       &flujo_ini,&
       &tempe_ini,&
       &postpro,&
       &fron_inm,&
       &directorio&
       &)
    !
    use malla, only : mi, nj, DBL, mic, njc
    implicit none
    INTEGER          :: itermax, paq_itera, iter_simple_max, iter_ecuaci_max
    REAL(kind=DBL)   :: Ra,Pr,dt,Ri_1,rel_pres,rel_vel,rel_ener
    REAL(kind=DBL)   :: conv_u,conv_p,conv_t,conv_resi,conv_paso
    CHARACTER(len=28):: entrada_u,entrada_v,entrada_tp
    CHARACTER(len=36):: directorio
    character(len=8) :: flujo_ini, tempe_ini
    logical          :: postpro, fron_inm
    !
    ! Se escribe la informaci\'on con la que se realiza la ejecuci\'on que
    ! produce los archivos de salida en el directorio nxxxmxxxRxxx/
    ! Esto ayuda a detectar errores de ejecuci\'on por la ausencia de este
    ! directorio
    !
    open(unit=10, file=directorio)
    write (10,*) 'numero de Rayleigh                    ', Ra
    write (10,*) 'numero de Prandtl                     ', Pr
    write (10,*) 'incremento de tiempo                  ', dt
    write (10,*) 'iteraciones maximas                   ', itermax
    write (10,*) 'paquete de iteraciones                ', paq_itera
    write (10,*) 'numero de Richardson                  ', Ri_1
    write (10,*) 'relajacion de la presion              ', rel_pres
    write (10,*) 'relajacion de la velocidad            ', rel_vel
    write (10,*) 'relajacion de la temperatura          ', rel_ener
    write (10,*) 'convergencia de la velocidad          ', conv_u
    write (10,*) 'convergencia de la temperatura        ', conv_t
    write (10,*) 'convergencia de la presion            ', conv_p
    write (10,*) 'convergencia del residuo              ', conv_resi
    write (10,*) 'convergencia del paso de tiempo       ', conv_paso
    write (10,*) 'iteraciones maximas de SIMPLE         ', iter_simple_max
    write (10,*) 'iteraciones maximas de las ecuaciones ', iter_ecuaci_max
    write (10,*) 'archivo de entrada para u             ', entrada_u
    write (10,*) 'archivo de entrada para v             ', entrada_v
    write (10,*) 'archivo de entrada para t y p         ', entrada_tp
    write (10,*) 'opcion de flujo inicial               ', flujo_ini
    write (10,*) 'opcion de temperatura inicial         ', tempe_ini
    write (10,*) 'opcion de postproceso                 ', postpro
    write (10,*) 'opcion de frontera inmersa            ', fron_inm
    close(unit=10)
    !
    !
  end subroutine postprocesa_parametros
  !
  !*******************************************************************
  !
  ! nusselt_promedio
  !
  ! Subrutina que calcula el n\'umero de nusselt promedio
  ! en paredes verticales
  !
  !*******************************************************************
  subroutine nusselt_promedio_y(&
       &xpo,ypo,deltaxpo,deltaypo,&
       &temp_o,nusselt0_o,nusselt1_o,&
       &i_oo,i_1o&
       &)
    implicit none
    ! $acc routine
    !
    !-------------------------------------
    !
    ! Variables de malla, nusselt y temperatura
    ! 
    real(kind=DBL), dimension(mi+1,nj+1), intent(in) :: temp_o
    real(kind=DBL), dimension(mi+1), intent(in)      :: xpo
    real(kind=DBL), dimension(nj+1), intent(in)      :: ypo
    real(kind=DBL), dimension(mi),   intent(in)      :: deltaxpo
    real(kind=DBL), dimension(nj),   intent(in)      :: deltaypo
    real(kind=DBL), intent(out)                      :: nusselt0_o, nusselt1_o
    integer, intent(in)                              :: i_oo, i_1o
    !
    ! Variables de interpolaci\'on para derivadas
    !
    real(kind=DBL)                  :: a,b,c,dx1,dx2,dy1,dy2
    real(kind=DBL), dimension(mi+1) :: derivada
    !
    ! Variables de interpolaci\'on para integrales
    !
    real(kind=DBL) :: alpha,beta,gamma
    real(kind=DBL) :: x_o,x_1,x_2,y_o,y_1,y_2
    !
    ! Variables auxiliares
    !
    integer :: ii,jj
    !----------------------------------------
    !
    ! C\'alculo para jj = 1 usando interpolaci\'on cuadr\'atica
    !
    derivada = 0.0_DBL
    jj = 1
    do ii = i_oo, i_1o
       ! jj = 1
       dy1 = temp_o(ii,jj+1)-temp_o(ii,jj)
       dy2 = temp_o(ii,jj+2)-temp_o(ii,jj+1)
       dx1 = ypo(jj+1)-ypo(jj)
       dx2 = ypo(jj+2)-ypo(jj+1)
       a = (dy2/dx2-dy1/dx1)/(dx1+dx2)
       b = dy2/dx2-a*(ypo(jj+1)+ypo(jj+2))
       c = temp_o(ii,jj)-a*ypo(jj)*ypo(jj)-b*ypo(jj)
       derivada(ii) = 2._DBL*a*ypo(jj)+b
    end do
    !
    ! Integra con interpolaci\'on cuadr\'atica
    !
    nusselt0_o = 0.0_DBL
    do ii = i_oo, i_1o, 2
       x_o   = xpo(ii)
       x_1   = xpo(ii+1)
       x_2   = xpo(ii+2)
       y_o   = derivada(ii)
       y_1   = derivada(ii+1)
       y_2   = derivada(ii+2)
       alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
       beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
       gamma = y_o-alpha*x_o*x_o-beta*x_o
       nusselt0_o = nusselt0_o+&
            alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+&
            &beta/2._DBL*(x_2*x_2-x_o*x_o)+&
            &gamma*(x_2-x_o)
    end do
    !
    ! Signo por la ley de Fourier
    !
    nusselt0_o=-nusselt0_o
    !
    !---------------------------------------------------------
    !
    ! C\'alculo para jj = nj usando interpolaci\'on cuadr\'atica
    !
    derivada = 0.0_DBL
    jj = nj-1
    do ii = i_oo, i_1o
       ! jj = nj-1
       dy1 = temp_o(ii,jj+1)-temp_o(ii,jj)
       dy2 = temp_o(ii,jj+2)-temp_o(ii,jj+1)
       dx1 = ypo(jj+1)-ypo(jj)
       dx2 = ypo(jj+2)-ypo(jj+1)
       a = (dy2/dx2-dy1/dx1)/(dx1+dx2)
       b = dy2/dx2-a*(ypo(jj+1)+ypo(jj+2))
       c = temp_o(ii,jj)-a*ypo(jj)*ypo(jj)-b*ypo(jj)
       derivada(ii) = 2._DBL*a*ypo(jj+2)+b
    end do
    !
    ! Integra con interpolaci\'on cuadr\'atica
    !
    nusselt1_o = 0.0_DBL
    do ii = i_oo, i_1o, 2
       x_o   = xpo(ii)
       x_1   = xpo(ii+1)
       x_2   = xpo(ii+2)
       y_o   = derivada(ii)
       y_1   = derivada(ii+1)
       y_2   = derivada(ii+2)
       alpha = ((y_2-y_1)/(x_2-x_1)-(y_1-y_o)/(x_1-x_o))/(x_2-x_o)
       beta  = (y_2-y_1)/(x_2-x_1)-alpha*(x_1+x_2)
       gamma = y_o-alpha*x_o*x_o-beta*x_o
       nusselt1_o = nusselt1_o+&
            alpha/3._DBL*(x_2*x_2*x_2-x_o*x_o*x_o)+&
            &beta/2._DBL*(x_2*x_2-x_o*x_o)+&
            &gamma*(x_2-x_o)
    end do
    !
    ! Signo por la ley de Fourier
    !
    nusselt1_o=-nusselt1_o   
    !
  end subroutine nusselt_promedio_y
  !
  !************************************************************
  !
  ! postprocess_vtk
  !
  ! Subrutina de postproceso en formato vtk (archivos binarios)
  !
  !************************************************************
  !
  subroutine postproceso_vtk(&
       &xo,yo,uo,vo,presso,tempo,bo,archivoo&
       &)
    use malla, only : mi, nj, DBL, mic, njc, zkc
    implicit none
    INTEGER :: i,j,k
    REAL(kind=DBL), DIMENSION(mi+1), INTENT(in)           :: xo
    REAL(kind=DBL), DIMENSION(nj+1), INTENT(in)           :: yo
    REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: uo, vo
    REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: tempo,presso
    REAL(kind=DBL), DIMENSION(mi+1,nj+1),   INTENT(in)    :: bo
    CHARACTER(46), INTENT(in)                             :: archivoo
    character(64)                                         :: mico,njco,zkco
    character(128)                                        :: npuntosc
    !
    ! Creaci\'on de cadenas de caracteres para el contenido de los archivos
    !
    write(mico,*) mi+1
    write(njco,*) nj+1
    write(zkco,*) 1
    write(npuntosc,*) (mi+1)*(nj+1)*1
    !
    !************************************
    ! VTK
    open(78, file = trim(archivoo), access='stream', convert="big_endian")

    write(78) '# vtk DataFile Version 2.3'//new_line(' ')
    write(78) '3D Mesh'//new_line(' ')
    write(78) 'BINARY'//new_line(' ')
    write(78) 'DATASET STRUCTURED_GRID'//new_line(' ')
    write(78) 'DIMENSIONS '//trim(mico)//trim(njco)//trim(zkco)//new_line('a')
    write(78) 'POINTS '//trim(npuntosc)//' float',new_line('a')
    do k = 1, 1
       do j = 1, nj+1
          do i = 1, mi+1
             write(78) real(xo(i)),real(yo(j)),0.0
          enddo
       enddo
    end do
    write(78) new_line('a')//'POINT_DATA '//trim(npuntosc)
    write(78) 'SCALARS PRESS float',new_line('a')
    write(78) 'LOOKUP_TABLE default',new_line('a')
    do k = 1, 1
       do j =1, nj+1
          do i =1, mi+1
             write(78) real(presso(i,j))
             ! write(78) real(bo(i,j))
          end do
       end do
    end do
    write(78) new_line('a')//'SCALARS TEMPER float',new_line('a')
    write(78) 'LOOKUP_TABLE default',new_line('a')
    do k = 1, 1
       do j =1, nj+1
          do i =1, mi+1
             write(78) real(tempo(i,j))
          end do
       end do
    end do
    write(78) new_line('a')//'VECTORS VELOCITY float',new_line('a')
    do k = 1, 1
       do j =1, nj+1
          do i =1, mi+1
             write(78) real(uo(i,j)),real(vo(i,j)),0.0 
          end do
       end do
    end do
    close(78)
    !
    ! 100 FORMAT(3(f12.6));
    ! 110 FORMAT(A);
    ! 111 FORMAT(A,/);
    ! 120 FORMAT(A,I4,I4,I4);
    ! 130 FORMAT(A,I10,A);
    ! 140 FORMAT(A,I10);
    !
  end subroutine postproceso_vtk
  !
  !************************************************************
  !
  ! postprocess_bin
  !
  ! Subrutina de postproceso en formato binario
  !
  !************************************************************
  !
  subroutine postproceso_bin(xuo,yvo,xpo,ypo,&
       &uo,vo,presso,tempo,bo, &
       &Rxc                    &
       )
    use malla, only : mi, nj, DBL, mic, njc, zkc 
    implicit none
    integer :: i,j,k
    real(kind=DBL), DIMENSION(mi), INTENT(in)           :: xuo
    real(kind=DBL), DIMENSION(nj), INTENT(in)           :: yvo
    real(kind=DBL), DIMENSION(mi+1), INTENT(in)         :: xpo
    real(kind=DBL), DIMENSION(nj+1), INTENT(in)         :: ypo   
    real(kind=DBL), DIMENSION(mi,nj+1),   INTENT(in)    :: uo
    real(kind=DBL), DIMENSION(mi+1,nj),   INTENT(in)    :: vo
    real(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)    :: tempo,presso
    real(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)    :: bo
    character(6),   intent(in)                          :: Rxc
    ! character(64)                                       :: Rxc=repeat(' ',64)
    ! character(46),  INTENT(in)                          :: archivoo
    ! character(64)                                       :: mico,njco,zkco
    !
    ! Creaci\'on de cadenas de caracteres para el contenido de los archivos
    !
    ! write(mico,*) mi+1
    ! write(njco,*) nj+1
    ! write(zkco,*) 1
    ! Rxc = entero_caracter(ceiling(Rx))
    !********************************
    !*** Formato de escritura dat ***
    open(unit=2,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rxc)//'u.bin',access='stream')
    !write(2) placa_min,placa_max,itera_total,ao
    do j = 1, nj+1
       do i = 1, mi
          write(2) xuo(i),ypo(j),uo(i,j)
       end do
    end do
    close(unit=2)
    ! -----------
    open(unit=3,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rxc)//'v.bin',access='stream')
    ! WRITE(3) placa_min,placa_max,itera_total,ao
    DO j = 1, nj
       DO i = 1, mi+1
          WRITE(3) xpo(i),yvo(j),vo(i,j)
       END DO
    END DO
    CLOSE(unit=3)
    ! -----------
    OPEN(unit=4,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rxc)//'p.bin',access='stream')
    ! WRITE(4,*) placa_min,placa_max,itera_total,ao
    DO j = 1, nj+1
       DO i = 1, mi+1
          WRITE(4) tempo(i,j),presso(i,j)
       END DO
    END DO
    CLOSE(unit=4)

  end subroutine postproceso_bin
  !
  !************************************************************
  !
  ! entero_caracter
  !
  ! Subrutina que devuelve una cadena a partir de un entero
  !
  !************************************************************
  !
  function entero_caracter(entero)

    implicit none

    character(6)        :: entero_caracter 
    integer, intent(in) :: entero

    integer             :: uni, dec, cen, mil, dmi
    character(1)        :: un,  de,  ce,  mi,  dm

    dmi = entero/10000
    mil = ( entero-dmi*10000 ) / 1000
    cen = ( entero-dmi*10000-mil*1000 ) / 100
    dec = ( entero-dmi*10000-mil*1000-cen*100 ) / 10
    uni = ( entero-dmi*10000-mil*1000-cen*100-dec*10 )

    write(un,16) uni; 16 format(I1)
    write(de,16) dec
    write(ce,16) cen
    write(mi,16) mil
    write(dm,16) dmi

    entero_caracter = dm//mi//ce//de//un

  end function entero_caracter

end module postproceso
