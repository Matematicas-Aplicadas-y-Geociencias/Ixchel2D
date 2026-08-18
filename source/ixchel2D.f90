!
! Programa que implementa el algoritmo SIMPLE para resolver las ecuaciones
! de Navier-Stokes y la energ\'ia
!
! autor: J.C. Cajas
!
!
PROGRAM IXCHEL2D
  !
  use omp_lib
  !
  ! Variables de la malla, volumen de control y factores de interpolaci\'on
  !
  use malla, only : mi, nj, DBL
  use malla, only : mic, njc, zkc
  use malla, only : xu, yv, xp, yp
  use malla, only : deltaxp, deltayp, deltaxu
  use malla, only : deltayu, deltaxv, deltayv
  use malla, only : fexp, feyp, fexu, feyv
  use malla, only : form24, form25, form26, form27
  use malla, only : lectura_mallas_escalonadas
  !
  ! Subrutina para imponer condiciones de frontera
  !
  use cond_frontera, only : impone_cond_frontera
  !
  ! Componentes de velocidad, presi\'on
  ! residuos de las ecuaciones de momento y correcci\'on de la presi\'on
  ! coeficientes de difusi\'on y criterios e convergencia
  !
  use ec_momento, only : u, u_ant, du, au, Resu, Ri, Riy, fu
  use ec_momento, only : v, v_ant, dv, av, Resv, fv
  use ec_momento, only : fuente_con_u, fuente_lin_u
  use ec_momento, only : fuente_con_v, fuente_lin_v
  use ec_momento, only : pres, corr_pres, dcorr_pres, fcorr_pres
  use ec_momento, only : uf, vf, b_o
  use ec_momento, only : maxbo, conv_u, conv_p, conv_paso, rel_pres, rel_vel
  use ec_momento, only : gamma_momen
  !
  ! Rutinas de ensamblaje de la ec. de momento, correcci\'on de presi\'on y residuo
  !
  use ec_momento, only : ensambla_velu_x, ensambla_velu_y
  use ec_momento, only : ensambla_velv_x, ensambla_velv_y
  use ec_momento, only : ensambla_corr_pres_x, ensambla_corr_pres_y
  use ec_momento, only : residuo_u
  use ec_momento, only : ini_frontera_uv
  use ec_momento, only : cond_front_ua, cond_front_ub  
  use ec_momento, only : cond_front_uc, cond_front_ud
  use ec_momento, only : cond_front_va, cond_front_vb  
  use ec_momento, only : cond_front_vc, cond_front_vd
  use ec_momento, only : condicion_inicial_flujo
  !
  ! Variables para la ecuaci\'on de la energ\'ia
  ! temperatura, coeficiente de difusi\'on y criterios de convergencia
  !
  use ec_energia, only : temp, temp_ant, dtemp, Restemp
  use ec_energia, only : fuente_con_t, fuente_lin_t
  use ec_energia, only : ftemp, gamma_energ, conv_t,rel_ener
  use ec_energia, only : ini_frontera_t
  use ec_energia, only : cond_front_ta, cond_front_tb  
  use ec_energia, only : cond_front_tc, cond_front_td
  use ec_energia, only : condicion_inicial_tempe
  ! 
  ! Rutina de ensamblaje de la ec. de energ\'ia
  !
  use ec_energia, only : ensambla_energia_x, ensambla_energia_y
  !
  ! Rutina que determina viscosidades para fronteras inmersas
  !
  use frontera_inmersa, only : lectura_cuerpo
  use frontera_inmersa, only : definir_cuerpo
  use frontera_inmersa, only : cond_front_inmersa
  !
  ! Rutinas de soluci\'on de ecuaciones
  !
  use solucionador, only : tridiagonal
  !
  ! Subrutinas de entrada y postproceso
  !
  use postproceso, only  : inicializa_promedio_perfil
  use postproceso, only  : finaliza_promedio_perfil
  use postproceso, only  : lectura_archivo_parametros
  use postproceso, only  : nusselt_promedio_y, postproceso_vtk
  use postproceso, only  : postproceso_bin, entero_caracter
  use postproceso, only  : postprocesa_parametros, postpro_promedio
  use postproceso, only  : lectura_archivo_prom
  !
  ! Variables auxiliares para bucles y n\'umero de iteraciones
  !
  implicit none
  !include  'omp_lib.h'
  !
  INTEGER :: itera_total,itera,itera_inicial,i_1,paq_itera,itermax
  integer :: iter_ecuaci, iter_ecuaci_max
  integer :: iter_simple, iter_simple_max
  INTEGER :: millar,centena,decena,unidad,decima,id,nthreads
  integer :: ii, jj, kk, ll, iter, auxiliar, ldiv
  integer :: stream1 = 1, stream2 = 2, stream3 = 3
  ! ------------------------------------------------------------
  !
  ! Variables para los archivos de la entrada de datos
  !
  character(len=28) :: entrada_u,entrada_v,entrada_tp
  character(len=36) :: directorio
  character(len=36) :: file_name
  character(len=7)  :: adimen
  !
  ! Variables para los archivos de postproceso
  !
  character(len=46) :: archivo=repeat(' ',46)
  logical           :: postprocesar = .false.
  !
  ! Variables para opciones de inicializaci\'on
  !
  character(len=8)  :: flujo_ini, tempe_ini
  !
  ! Variables para opci'on de frontera inmersa
  !
  logical           :: front_inmersa = .false.
  !
  ! *******************************************
  !
  REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: entropia_calor,entropia_viscosa,entropia,gamma_t
  REAL(kind=DBL) :: temp_med,nusselt0,nusselt1,entropia_int,temp_int,gamma_s,residuo,error
  REAL(kind=DBL) :: conv_resi
  !
  ! Tamanio del dominio y ubicaci\'on de las fuentes de calor
  !
  REAL(kind=DBL)   :: ao
  integer          :: placa_min, placa_max
  !
  ! Coeficientes para las matrices 
  !
  real(kind=DBL), dimension(mi+1,nj+1) :: AI, AC, AD, Rx
  real(kind=DBL), dimension(nj+1,mi+1) :: BS, BC, BN, Ry
  !
  REAL(kind=DBL)   :: tiempo,tiempo_inicial,dt,Ra,Pr,Ri_1
  REAL(kind=DBL)   :: a_ent,lambda_ent
  CHARACTER(len=1) :: dec,un,de,ce,m
  ! CHARACTER(len=3) :: njc,mic,Rec
  CHARACTER(len=6) :: Rec=repeat(' ',6)
  CHARACTER(len=5) :: sample

  !****************************************
  !Variables de caracterizaci'on del fluido
  REAL(kind=DBL) :: temp_ref,visc_cin,dif_term,cond_ter,cons_gra
  real(kind=DBL) :: coef_exp,long_ref,dens_ref
  !****************************
  !declaraci´on de variable DBL
  REAL(kind=DBL) :: var2=0.0_DBL
  !*******************************************
  !
  ! Auxiliares de interpolaci\'on
  !
  real(kind=DBL) :: ui, ud, vs, vn
  real(kind=DBL) :: di, dd, ds, dn
  real(kind=DBL) :: gammai, gammad
  real(kind=DBL) :: gammas, gamman
  real(kind=DBL) :: deltax, deltay
  ! real(kind=DBL) :: temp_int
  !
  !***********************************************************
  !
  ! Mensaje de bienvenida
  !
  write(*,*) "-------------------------------------------------------"
  write(*,*) " "
  write(*,*) "                      IXCHEL2D"
  write(*,*) " "
  write(*,*) "           Simulaciones en termofluidos 2D      "
  write(*,*) "-------------------------------------------------------"
  write(*,*) "              "   
  !
  ! Valor por defecto de la variable de control de postproceso
  !
  postprocesar = .false.
  !
  ! Par'ametros para la simulaci'on
  !
  call lectura_archivo_parametros(&
       &adimen,          &
       &Ra,              &
       &Pr,              &
       &Ri_1,            &
       &dt,              &
       &paq_itera,       &
       &itermax,         &
       &rel_pres,        &
       &rel_vel,         &
       &rel_ener,        &
       &conv_u,          &
       &conv_t,          &
       &conv_p,          &
       &conv_resi,       &
       &conv_paso,       &
       &iter_ecuaci_max, &
       &iter_simple_max, &
       &entrada_u,       &
       &entrada_v,       &
       &entrada_tp,      &
       &flujo_ini,       &
       &tempe_ini,       &
       &postprocesar,    &
       &front_inmersa )
  !
  !--------------------------------------------------------------
  !
  ! Definici'on de caracteres para nombres de archivos y mensajes
  !
  mic = entero_caracter(mi)
  njc = entero_caracter(nj)
  !
  !------------------------------------------
  !
  ! Inicializaci\'on de los t\'erminos fuente
  !
  fuente_con_u = 0.0_DBL
  fuente_lin_u = 0.0_DBL
  !
  fuente_con_v = 0.0_DBL
  fuente_lin_v = 0.0_DBL
  !
  fuente_con_t = 0.0_DBL
  fuente_lin_t = 0.0_DBL
  !
  Ri           = Ri_1
  Riy          = 0.0_DBL
  !
  !-------------------------
  !
  !Selecci\'on de la adimensionalizaci\'on
  !
  if( trim(adimen) == 'natural' )then
     ! --------------------
     !
     ! Convecci\'on natural
     ! 
     gamma_momen = sqrt(Pr/Ra) 
     gamma_energ = sqrt(1._DBL/(Pr*Ra))
     Ri          = 0.0_DBL
     Riy         =-1.0_DBL
     Rec         = entero_caracter(ceiling(sqrt(Ra)))
     !
  else if( trim(adimen) == 'mixta' )then
     ! ------------------
     !
     ! Convecci\'on mixta
     ! 
     ! gamma_s     = 10._DBL*(1._DBL/(Ra*Pr))
     gamma_momen = 1.0_DBL/(Ra)    ! n\'umero de Reynolds
     gamma_energ = 1.0_DBL/(Ra*Pr) ! N'umero de P'eclet
     Ri          = 0.0_DBL
     Riy         = Ri_1
     Rec         = entero_caracter(ceiling(Ra))
     !
  end if
  !
  !------------------------------------------------
  !
  ! Verificación de existencia de directorio
  !
  directorio = 'n'//trim(njc)//'m'//trim(mic)//'R'//trim(Rec)//'/info_entrada.out'
  call postprocesa_parametros(&
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
       &postprocesar,&
       &front_inmersa,&
       &directorio)
  ! ----------------------------------------------------------------
  !
  ! Lectura de las mallas escalonadas e inicializaci\'on de arreglos
  !
  write(*,*) "Ixchel2D: Inicia lectura de mallas"
  call lectura_mallas_escalonadas(entrada_u,entrada_v,entrada_tp,&
       &u_ant,v_ant,pres,temp_ant,&
       &xp,yp,xu,yv,&
       &deltaxp,deltayp,&
       &deltaxu,deltayu,&
       &deltaxv,deltayv,&
       &fexp,feyp,fexu,feyv,&
       &ao,placa_min,placa_max,itera_inicial)
  write(*,*) "Ixchel2D: Finaliza lectura de mallas"
  write(*,*) "------------------------------------"
  ! !*****************
  !valores iniciales
  tiempo_inicial = itera_inicial*dt
  itera_total    = itera_inicial
  call condicion_inicial_flujo(flujo_ini)
  ! $acc parallel
  call condicion_inicial_tempe(tempe_ini)
  ! $acc end parallel
  !
  ! u_ant = 1.0_DBL
  ! u_ant = 0.0_DBL
  ! v_ant = 0.0_DBL
  ! temp_ant = 0.0_DBL
  u         = u_ant
  v         = v_ant
  uf        = 0.0_DBL
  vf        = 0.0_DBL
  temp      = temp_ant
  Resu      = 0.0_DBL
  corr_pres = 0.0_DBL !* dfloat(mi+1-ii)/dfloat(mi)
  au    = 1.e40_DBL
  av    = 1.e40_DBL
  b_o   = 0.0_DBL
  itera = 0
  iter_ecuaci = 0
  iter_simple = 0
  maxbo   = 0.0_DBL
  error   = 0.0_DBL
  residuo = 0.0_DBL
  ! --------------------------------------
  !
  ! Lectura de las condiciones de frontera
  !
  call ini_frontera_uv()
  call ini_frontera_t()
  !
  ! ------------------------------------------------
  !
  ! Construcci\'on de s\'olidos con frontera inmersa 
  !
  if( front_inmersa ) call lectura_cuerpo()
  !
  ! call definir_cuerpo(gamma_momen, gamma_energ, 'horno')
  !
  !--------------------------------------------------
  !
  ! Lectura de archivo para los promedios de perfiles
  !
  call lectura_archivo_prom()
  !
  file_name = 'n'//trim(njc)//'m'//trim(mic)//'R'//trim(Rec)//&
       &'/perfil_prome.out'
  !
  ! Se abre el archivo para almacenar los promedios en los
  ! perfiles
  !
  call inicializa_promedio_perfil(file_name)
  !
  !----------------------------------------------
  !
  !************************************************
  !escribe las caracter´isticas de las variable DBL
  write(*,100) 'Ixchel2D: Doble ',kind(var2),precision(var2),range(var2)
  write(*,*)' '
  !escribe informaci'on de los parametros usados
  write(*,101) Rec,Pr,Ri_1,rel_pres,rel_vel
  write(*,102) itera_inicial,mi,nj
  write(*,*)' '

  !--------------------------------------------
  !
  ! Inicio del repetidor principal
  !
  do ll = 1, itermax/paq_itera
     !
     ! Inicio del paquete de iteraciones
     !
     !------------------------------------------
     !
     ! Apertura de la regi\'on de datos paralela
     !
     !$omp     target data map(tofrom:                                               &
     !$omp     u(1:mi,1:nj+1),v(1:mi+1,1:nj),                            &
     !$omp     pres(1:mi+1,1:nj+1),temp(1:mi+1,1:nj+1),                  &
     !$omp     corr_pres(1:mi+1,1:nj+1),                                 &
     !$omp     u_ant(1:mi,1:nj+1),v_ant(1:mi+1,1:nj),                    &
     !$omp     temp_ant(1:mi+1,1:nj+1),b_o(1:mi+1,1:nj+1))               &
     !$omp     map(to:                                                   &
     !$omp     tiempo_inicial,                                           &
     !$omp     Resu(1:mi,1:nj+1),                                        &
     !$omp     au(1:mi,1:nj+1),av(1:mi+1,1:nj),                          &
     !$omp     gamma_momen(1:mi+1,1:nj+1),gamma_energ(1:mi+1,1:nj+1),    &
     !$omp     deltaxp(1:mi),deltayp(1:nj),deltaxu(1:mi),deltayu(1:nj),  &
     !$omp     deltaxv(1:mi),deltayv(1:nj),                              &
     !$omp     fexp(1:mi),feyp(1:nj),fexu(1:mi),feyv(1:nj),              &
     !$omp     Ri(1:mi,1:nj+1),Riy(1:mi+1,1:nj+1),dt,                    &
     !$omp     fuente_con_u(1:mi,1:nj+1), fuente_lin_u(1:mi,1:nj+1),     &
     !$omp     fuente_con_v(1:mi+1,1:nj), fuente_lin_v(1:mi+1,1:nj),     &
     !$omp     fuente_con_t(1:mi+1,1:nj+1), fuente_lin_t(1:mi+1,1:nj+1), &
     !$omp     rel_vel,conv_u,conv_p,                                    &
     !$omp     rel_ener,conv_t,placa_min,placa_max,                      &
     !$omp     cond_front_ua,cond_front_ub,cond_front_uc, cond_front_ud, &
     !$omp     cond_front_va,cond_front_vb,cond_front_vc, cond_front_vd, &
     !$omp     cond_front_ta,cond_front_tb,cond_front_tc, cond_front_td) &
     !$omp     map(alloc:                                                &
     !$omp     AI(1:mi+1,1:nj+1),AC(1:mi+1,1:nj+1),                      &
     !$omp     AD(1:mi+1,1:nj+1),Rx(1:mi+1,1:nj+1),                      &
     !$omp     BS(1:nj+1,1:mi+1),BC(1:nj+1,1:mi+1),                      &
     !$omp     BN(1:nj+1,1:mi+1),Ry(1:nj+1,1:mi+1),                      &
     !$omp     fu(1:mi,1:nj+1),du(1:mi,1:nj+1),                          &
     !$omp     fv(1:mi+1,1:nj),dv(1:mi+1,1:nj),                          &
     !$omp     fcorr_pres(1:mi+1,1:nj+1),dcorr_pres(1:mi+1,1:nj+1),      &
     !$omp     ftemp(1:mi+1,1:nj+1),dtemp(1:mi+1,1:nj+1)                 &
     !$omp     )
     !
     do kk=1,paq_itera
        !
        ! Inicio del algoritmo SIMPLE
        !
        ALGORITMO_SIMPLE: do
           !
           !-------------------------------------------
           !-------------------------------------------
           !
           ! Soluci\'on de la ecuaci\'on de momento
           !
           !-------------------------------------------
           !-------------------------------------------           
           ecuacion_momento: do
              !
              !$omp target teams distribute parallel do collapse(2)
              inicializacion_fu: do jj=1, nj+1
                 do ii = 1, mi
                    fu(ii,jj) = u(ii,jj)
                 end do
              end do inicializacion_fu
              !$omp end target teams distribute parallel do
              !------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de momento
              ! para u en direcci\'on y´
              !
              !---------------------------------------
              !
              !$omp target teams distribute parallel do collapse(2)
              bucle_uy_direccion_y: do jj = 2, nj
                 !
                 ! Llenado de la matriz
                 !
                 !
                 bucle_uy_direccion_x: do ii = 2, mi-1
                    call ensambla_velu_y(deltaxu,deltayu,deltaxp,&
                         &deltayv,fexp,feyp,fexu,gamma_momen,&
                         &fuente_con_u,fuente_lin_u,&
                         &u,u_ant,v,&
                         &temp,pres,Ri,dt,rel_vel,&
                         &BS,BC,BN,Ry,&
                         &jj,ii&
                         &)
                 end do bucle_uy_direccion_x
              end do bucle_uy_direccion_y
              !$omp end target teams distribute parallel do
              !
              ! Condiciones de frontera para u 
              !
              !-------------------------------------------------
              !
              ! Region paralela para imponer las cond. de front.
              ! en la GPU
              !
              !$omp target
              !-----------------------------------------------
              !
              ! lado b
              !
              call impone_cond_frontera(cond_front_ub,&
                   & BS,BC,BN,Ry, &
                   & nj+1,mi+1,   &
                   & mi,nj+1,     &
                   & au )
              !-----------------------------------------------
              !
              ! lado d
              !
              call impone_cond_frontera(cond_front_ud,&
                   & BS,BC,BN,Ry, &
                   & nj+1,mi+1,   &
                   & mi,nj+1,     &
                   & au )
              !$omp end target
              !
              !$omp target teams distribute parallel do
              solucion_momento_uy: do ii = 2, mi-1
                 
                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),&
                      &Ry(1:nj+1,ii),nj+1)

              end do solucion_momento_uy
              !$omp end target teams distribute parallel do
              !
              !----------------------------------
              !
              ! Actualizaci\'on de la velocidad u
              !
              !$omp target teams distribute parallel do collapse(2)
              do ii = 2, mi-1
                 do jj = 1, nj+1
                    u(ii,jj) = Ry(jj,ii)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !              
              !------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de momento
              ! para u en direcci\'on x
              !
              !$omp target teams distribute parallel do collapse(2)
              bucle_ux_direccion_y: do jj = 2, nj
                 !
                 ! Llenado de la matriz
                 !
                 !
                 bucle_ux_direccion_x: do ii = 2, mi-1
                    call ensambla_velu_x(deltaxu,deltayu,deltaxp,&
                         &deltayv,fexp,feyp,fexu,gamma_momen,&
                         &fuente_con_u,fuente_lin_u,&
                         &u,u_ant,v,&
                         &temp,pres,Ri,dt,rel_vel,&
                         &AI,AC,AD,Rx,au,&
                         &ii,jj&
                         &)
                 end do bucle_ux_direccion_x
              end do bucle_ux_direccion_y
              !$omp end target teams distribute parallel do
              !
              ! Condiciones de frontera para u 
              !
              !-------------------------------------------------
              !
              ! Region paralela para imponer las cond. de front.
              ! en la GPU
              !
              ! lado a
              !
              !$omp target
              call impone_cond_frontera(cond_front_ua,&
                   & AI,AC,AD,Rx, &
                   & mi+1,nj+1,   &
                   & mi,nj+1,     &
                   & au )
              !
              !-------------------------------------
              !
              ! lado c
              !
              call impone_cond_frontera(cond_front_uc,&
                   & AI,AC,AD,Rx, &
                   & mi+1,nj+1,   &
                   & mi,nj+1,     &
                   & au )
              !$omp end target
              !
              !-------------------------------------
              !
              ! Soluci\'on del sistema de ecuaciones
              !
              !$omp target teams distribute parallel do
              solucion_momento_ux: do jj = 2, nj
                 
                 call tridiagonal(AI(1:mi,jj),AC(1:mi,jj),AD(1:mi,jj),Rx(1:mi,jj),mi)
                 
              end do solucion_momento_ux
              !$omp end target teams distribute parallel do
              !----------------------------------
              !
              ! Actualizaci\'on de la velocidad u
              !
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj
                 do ii = 1, mi
                    u(ii,jj) = Rx(ii,jj)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              !$omp target teams distribute parallel do collapse(2)
              inicializacion_fv: do jj=1, nj
                 do ii = 1, mi+1
                    fv(ii,jj) = v(ii,jj)
                 end do
              end do inicializacion_fv
              !$omp end target teams distribute parallel do
              !---------------------------------------------
              !
              ! Se ensambla la velocidad v en direcci\'on y
              !
              !---------------------------------------------
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj-1
                 !
                 do ii = 2, mi
                    call ensambla_velv_y(deltaxv,deltayv,deltaxu,&
                         &deltayp,fexp,feyp,feyv,gamma_momen,&
                         &fuente_con_v,fuente_lin_v,&
                         &v,v_ant,u,&
                         &temp,pres,Riy,dt,rel_vel,&
                         &BS,BC,BN,Ry,&
                         &jj,ii&
                         &)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              ! Condiciones de frontera para v
              !
              !-------------------------------------------------
              !
              ! Region paralela para imponer las cond. de front.
              ! en la GPU
              !              
              !
              ! lado b
              !
              !$omp target
              call impone_cond_frontera(cond_front_vb,&
                   & BS,BC,BN,Ry, &
                   & nj+1,mi+1,   &
                   & mi+1,nj,     &
                   & av )
              
              ! !$acc parallel loop vector !async(stream1)
              ! bucle_direccionyv: do ii = 2, mi
              !    !***********************
              !    !Condiciones de frontera
              !    BC(1,ii)     = 1._DBL
              !    BN(1,ii)     = 0.0_DBL
              !    Ry(1,ii)     = 0.0_DBL
              !    av(1,ii)     = 1.0e40_DBL
              !    !
              ! end do bucle_direccionyv
              ! $acc parallel
              !-----------------------------------------------
              !
              ! lado d
              !
              call impone_cond_frontera(cond_front_vd,&
                   & BS,BC,BN,Ry, &
                   & nj+1,mi+1,   &
                   & mi+1,nj,     &
                   & av )         
              !$omp end target
              !
              !$omp target teams distribute parallel do
              solucion_momento_vy: do ii = 2, mi

                 call tridiagonal(BS(1:nj,ii),BC(1:nj,ii),BN(1:nj,ii),Ry(1:nj,ii),nj)
                
              end do solucion_momento_vy
              !$omp end target teams distribute parallel do
              !----------------------------------
              !
              ! Actualizaci\'on de la velocidad v
              !
              !$omp target teams distribute parallel do collapse(2)
              do ii = 2, mi
                 do jj = 1, nj
                    v(ii,jj) = Ry(jj,ii)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              !---------------------------
              !
              ! Se ensambla la velocidad v en direcci\'on x
              !
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj-1
                 ! $acc loop vector
                 do ii = 2, mi
                    call ensambla_velv_x(deltaxv,deltayv,deltaxu,&
                         &deltayp,fexp,feyp,feyv,gamma_momen,&
                         &fuente_con_v,fuente_lin_v,&
                         &v,v_ant,u,&
                         &temp,pres,Riy,dt,rel_vel,&
                         &AI,AC,AD,Rx,av,&
                         &ii,jj&
                         &)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              ! Condiciones de frontera para v
              !
              !-----------------------------------------------
              !
              ! Region paralela para imponer las cond. de front.
              ! en la GPU
              ! 
              ! lado a
              !
              !$omp target
              call impone_cond_frontera(cond_front_va,&
                   & AI,AC,AD,Rx, &
                   & mi+1,nj+1,   &
                   & mi+1,nj,     &
                   & av )
       
              !-----------------------------------------------
              !
              ! lado c
              !
              call impone_cond_frontera(cond_front_vc,&
                   & AI,AC,AD,Rx, &
                   & mi+1,nj+1,   &
                   & mi+1,nj,     &
                   & av )            
              !$omp end target
              !
              !------------------------------------
              !
              ! Soluci\'on de las ecs. de momento v
              !
              !$omp target teams distribute parallel do
              solucion_momento_vx: do jj = 2, nj-1

                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),&
                      &Rx(1:mi+1,jj),mi+1)
                 
              end do solucion_momento_vx
              !$omp end target teams distribute parallel do
              !----------------------------------
              !
              ! Actualizaci\'on de la velocidad v
              !
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj-1
                 do ii = 1, mi+1
                    v(ii,jj) = Rx(ii,jj)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              ! error de la ecuacion de momento
              !
              error = 0.0_DBL
              !$omp target teams distribute parallel do collapse(2) reduction(+:error)
              calculo_diferencias_dv: do jj=2, nj-1
                 do ii = 2, mi

                    error = error + dabs(v(ii,jj)-fv(ii,jj))*deltaxp(ii)*deltayv(jj)

                 end do
              end do calculo_diferencias_dv
              !$omp end target teams distribute parallel do
              ! error = sqrt(error)
              !
              ! Criterio de convergencia de la velocidad
              !
              if ( error < conv_u ) then
                 iter_ecuaci = 0
                 ! write(101,*) 'velocidad ', error
                 exit
              else if (iter_ecuaci > iter_ecuaci_max) then
                 iter_ecuaci = 0
                 ! write(*,*) "ADVER. MOMEN: convergencia no alcanzada ", error
                 exit
              else
                 iter_ecuaci = iter_ecuaci+1
                 ! write(101,*) 'velocidad ', error
              end if
              !            
           end do ecuacion_momento
           !$acc wait
           !-----------------------------------------
           !-----------------------------------------
           !
           ! Se calcula la correcci'on de la presi'on
           !
           !-----------------------------------------
           !-----------------------------------------
           !
           !$omp target teams distribute parallel do collapse(2)
           inicializa_corrector_presion: do jj = 1, nj+1
              do ii = 1, mi+1
                 corr_pres(ii,jj) = 0.0_DBL
                 fcorr_pres(ii,jj)= 0.0_DBL
              end do
           end do inicializa_corrector_presion
           !$omp end target teams distribute parallel do
           !
           correccion_presion: do
              !
              !$omp target teams distribute parallel do collapse(2) &
              !$omp map(from: fcorr_pres) map(to:corr_pres)
              inicializa_fcorr_press: do jj=2, nj
                 do ii = 2, mi
                    fcorr_pres(ii,jj) = corr_pres(ii,jj)
                 end do
              end do inicializa_fcorr_press
              !$omp end target teams distribute parallel do
              !---------------------------------------------------
              !frontera inmersa
              !
              !call cond_front_inmersa(au, av, 'cuadr')
              !-----------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de la presi\'on en y
              !
              !$omp target teams distribute parallel do collapse(2)
              do ii = 2, mi
                 !
                 do jj = 2, nj
                    call ensambla_corr_pres_y(deltaxp,deltayp,&
                         &deltaxu,deltayv,&
                         &u,v,&
                         &corr_pres,rel_pres,&
                         &BS,BC,BN,Ry,au,av,&
                         &ii,jj)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !-------------------------
              !
              ! Condiciones de frontera
              !AC_o
              !$omp target teams distribute parallel do
              bucle_direccionxe: do ii = 2, mi
                 !***********************
                 !Condiciones de frontera
                 BC(1,ii)     = 1._DBL
                 BN(1,ii)     = 0.0_DBL
                 Ry(1,ii)     = 0.0_DBL
                 !
                 BC(nj+1,ii)  = 1._DBL
                 BS(nj+1,ii)  = 0.0_DBL
                 Ry(nj+1,ii)  = 0.0_DBL
              end do bucle_direccionxe
              !$omp end target teams distribute parallel do
              !---------------------------------------------------
              !
              ! imponer correccion de la presion en frontera inmersa
              !
              !call cond_front_inmersa(BS, BC, BN, Ry, 'cuadr')
              !---------------------------------------------------
              !
              ! Soluci\'on de la correcci\'on de la presi\'on en y
              !
              !
              !$omp target teams distribute parallel do
              solucion_presion_y: do ii = 2, mi

                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),&
                      &Ry(1:nj+1,ii),nj+1)
  
              end do solucion_presion_y
              !$omp end target teams distribute parallel do
              !----------------------------------------------------
              !
              ! Actualizaci\'on del corrector de la presi\'on en y
              !
              !$omp target teams distribute parallel do collapse(2)
              do ii = 2, mi
                 do jj = 1, nj+1
                    corr_pres(ii,jj) = Ry(jj,ii)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              !-----------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de la presi\'on en x
              !
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj
                 !
                 do ii = 2, mi
                    call ensambla_corr_pres_x(deltaxp,deltayp,&
                         &deltaxu,deltayv,&
                         &u,v,b_o,&
                         &corr_pres,rel_pres,&
                         &AI,AC,AD,Rx,au,av,&
                         &jj,ii)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              !$omp target teams distribute parallel do
              do jj = 2, nj
                 !------------------------
                 ! Condiciones de frontera
                 AC(1,jj) = 1._DBL
                 AD(1,jj) = 0._DBL
                 Rx(1,jj) = 0._DBL
                 !
                 AI(mi+1,jj) = 0.0_DBL
                 AC(mi+1,jj) = 1.0_DBL
                 Rx(mi+1,jj) = 0.0_DBL
              end do
              !$omp end target teams distribute parallel do
              !---------------------------------------------------
              !
              ! imponer correccion de la presion en frontera inmersa
              !
              !call cond_front_inmersa(AI, AC, AD, Rx, 'cuadr')
              !---------------------------------------------------
              !
              ! Soluci\'on de la correcci\'on de la presi\'on en x
              !       
              !$omp target teams distribute parallel do
              solucion_presion_x: do jj = 2, nj

                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),&
                      &Rx(1:mi+1,jj),mi+1)
                 
              end do solucion_presion_x
              !$omp end target teams distribute parallel do
              !
              !----------------------------------------------------
              !
              ! Actualizaci\'on del corrector de la presi\'on en x
              !
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj
                 do ii = 1, mi+1
                    corr_pres(ii,jj) = Rx(ii,jj)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              ! C\'alculo de diferencias y criterio de convergencia
              !
              error = 0.0_DBL
              maxbo = 0.0_DBL
              !
              !$omp target teams distribute parallel do collapse(2) reduction(+:error)
              calculo_dif_corr_pres: do jj=2, nj
                 do ii=2, mi

                    error = error + dabs(corr_pres(ii,jj)-fcorr_pres(ii,jj))*&
                         & deltaxp(ii)*deltayp(jj)

                 end do
              end do calculo_dif_corr_pres
              !$omp end target teams distribute parallel do
              ! error=sqrt(error)
              !
              !$omp target teams distribute parallel do collapse(2) reduction(max:maxbo)
              calculo_dif_maxbo: do jj=2, nj
                 do ii=2, mi

                    maxbo = max(maxbo,dabs(b_o(ii,jj)))
                    
                 end do
              end do calculo_dif_maxbo
              !$omp end target teams distribute parallel do
              ! maxbo = sqrt(maxbo)
              !
              !-----------------------------------------------------
              !
              ! Critero de convergencia del corrector de la presi'on
              !
              !
              if( error<conv_p )then
                 !write(*,*)"PRES: conver ", error," con ",iter_ecuaci," iteraciones"
                 iter_ecuaci = 0
                 exit
              else if (iter_ecuaci > iter_ecuaci_max) then
                 iter_ecuaci = 0
                 ! write(*,*) "ADVER. PRES: convergencia no alcanzada ", error
                 exit
              else
                 iter_ecuaci = iter_ecuaci+1
                 ! write(*,*) 'corrector presion ',error,maxbo
              end if
              !
           end do correccion_presion
           !
           !--------------------------------------------
           !
           ! Se corrige la presion
           !
           !$omp target teams distribute parallel do collapse(2)
           do jj = 2, nj
              do ii = 2, mi
                 pres(ii,jj) = pres(ii,jj) + 0.6_DBL*corr_pres(ii,jj)
              end do
           end do
           !$omp end target teams distribute parallel do
           !---------------------------------
           !
           ! Se corrigen las velocidades
           !
           !$omp target teams distribute parallel do collapse(2)
           do jj = 2, nj-1
              do ii = 2, mi-1
                 u(ii,jj) = u(ii,jj)+deltayu(jj)*&
                      &(corr_pres(ii,jj)-corr_pres(ii+1,jj))/au(ii,jj)
                 v(ii,jj) = v(ii,jj)+deltaxv(ii)*&
                      &(corr_pres(ii,jj)-corr_pres(ii,jj+1))/av(ii,jj)
              end do
           end do
           !$omp end target teams distribute parallel do
           !
           !$omp target teams distribute parallel do
           do ii = 2, mi-1
              u(ii,nj) = u(ii,nj)+deltayu(nj)*&
                   &(corr_pres(ii,nj)-corr_pres(ii+1,nj))/au(ii,nj)
           end do
           !$omp end target teams distribute parallel do
           !
           !$omp target teams distribute parallel do
           do jj = 2, nj-1
              v(mi,jj) = v(mi,jj)+deltaxv(mi)*&
                   &(corr_pres(mi,jj)-corr_pres(mi,jj+1))/av(mi,jj)
           end do
           !$omp end target teams distribute parallel do
           !
           !--------------------------------------------------
           !--------------------------------------------------
           !
           ! Se resuelve la ecuaci\'on de balance de energ\'ia
           !
           !--------------------------------------------------
           !--------------------------------------------------
           solucion_energia: do
              !
              !$omp target teams distribute parallel do collapse(2)
              inicializacion_ftemp: do jj=2, nj
                 do ii = 2, mi
                    ftemp(ii,jj) = temp(ii,jj)
                 end do
              end do inicializacion_ftemp
              !$omp end target teams distribute parallel do
              !------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de la energ\'ia en y
              !
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj
                 !
                 do ii = 2, mi
                    call ensambla_energia_y(deltaxp,deltayp,&
                         &deltaxu,deltayv,fexu,feyv,gamma_energ,&
                         &fuente_con_t,fuente_lin_t,&
                         &u,v,&
                         &temp,temp_ant,dt,&
                         &rel_ener,placa_min,placa_max,&
                         &BS,BC,BN,Ry,&
                         &jj,ii&
                         &)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !-----------------------------------------------
              !
              ! Condiciones de frontera
              !
              !-----------------------------------------------
              !
              ! Region paralela para imponer las cond. de front.
              ! en la GPU
              ! 
              !------------
              ! lado b
              !
              !$omp target
              call impone_cond_frontera(cond_front_tb,&
                   & BS,BC,BN,Ry, &
                   & nj+1,mi+1,   &
                   & mi+1,nj+1) 
              !-----------------------------------------------
              !
              ! lado d
              !
              call impone_cond_frontera(cond_front_td,&
                   & BS,BC,BN,Ry, &
                   & nj+1,mi+1,   &
                   & mi+1,nj+1)
              !$omp end target
              !---------------------------------------------
              !
              ! Soluci\'on de la ecuaci\'on de la energ\'ia en y
              !
              !$omp target teams distribute parallel do
              solucion_energia_y: do ii = 2, mi

                 call tridiagonal(BS(1:nj+1,ii),BC(1:nj+1,ii),BN(1:nj+1,ii),&
                      &Ry(1:nj+1,ii),nj+1)

              end do solucion_energia_y
              !$omp end target teams distribute parallel do
              !
              !----------------------------------------
              !
              ! Actualizacion de la energia
              !
              !$omp target teams distribute parallel do collapse(2)
              do ii = 2, mi
                 do jj = 1, nj+1
                    temp(ii,jj) = Ry(jj,ii)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              !-----------------------------------------------
              !
              ! Se ensambla la ecuaci\'on de la energ\'ia en x
              !
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj
                 !
                 do ii = 2, mi
                    call ensambla_energia_x(deltaxp,deltayp,&
                         &deltaxu,deltayv,fexu,feyv,gamma_energ,&
                         &fuente_con_t,fuente_lin_t,&
                         &u,v,&
                         &temp,temp_ant,dt,&
                         &rel_ener,placa_min,placa_max,&
                         &AI,AC,AD,Rx,&
                         &ii,jj&
                         &)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !------------------------------------------
              !
              ! Condiciones de frontera
              !
              !------------------------------------------
              !
              ! Regi\'on paralela para imponer las
              ! condiciones de frontera en la GPU
              !
              !------------
              !
              ! lado a
              !
              !$omp target
              call impone_cond_frontera(cond_front_ta,&
                   & AI,AC,AD,Rx, &
                   & mi+1,nj+1,   &
                   & mi+1,nj+1)
              !-----------------------------------------------
              !
              ! lado c
              !
              call impone_cond_frontera(cond_front_tc,&
                   & AI,AC,AD,Rx, &
                   & mi+1,nj+1,   &
                   & mi+1,nj+1)
              !$omp end target
              !
              !---------------------------------------------
              !
              ! Soluci\'on de la ecuaci\'on de la energ\'ia en x
              !
              !$omp target teams distribute parallel do
              solucion_energia_x: do jj = 2, nj

                 call tridiagonal(AI(1:mi+1,jj),AC(1:mi+1,jj),AD(1:mi+1,jj),&
                      &Rx(1:mi+1,jj),mi+1)

              end do solucion_energia_x
              !$omp end target teams distribute parallel do
              !----------------------------------------
              !
              ! Actualizacion de la energia
              !
              !$omp target teams distribute parallel do collapse(2)
              do jj = 2, nj
                 do ii = 1, mi+1
                    temp(ii,jj) = Rx(ii,jj)
                 end do
              end do
              !$omp end target teams distribute parallel do
              !
              !----------------------------------------
              !
              ! error de la ecuaci\'on de la energ\'ia
              !
              error = 0.0_DBL
              !
              !$omp target teams distribute parallel do collapse(2) reduction(+:error)
              calculo_diferencias_dtemp: do jj = 2, nj
                 do ii = 2, mi
                    error = error + dabs(temp(ii,jj)-ftemp(ii,jj))*&
                         & deltaxp(ii)*deltayp(jj)
                 end do
              end do calculo_diferencias_dtemp
              ! error = sqrt(error)
              !$omp end target teams distribute parallel do
              !
              !------------------------------------------
              !
              ! Criterio de convergencia de la energ\'ia
              !
              !
              if( error < conv_t )then
                 iter_ecuaci = 0
                 ! write(*,*) 'temp', error
                 exit
              else if( iter_ecuaci > iter_ecuaci_max )then
                 iter_ecuaci = 0
                 ! write(*,*) "Adver. ENERG: convergencia no alcanzada ", &
                 ! error
                 exit
              else
                 iter_ecuaci = iter_ecuaci + 1
                 ! write(*,*) 'temp', error
              end if
              !
           end do solucion_energia
           !$acc wait
           !
           !--------------------------------------------
           !--------------------------------------------
           !
           ! Criterio de convergencia del paso de tiempo
           !
           !--------------------------------------------
           !--------------------------------------------
           !
           !$omp target teams distribute parallel do collapse(2)
           bucle_residuo_direccion_y: do jj = 2, nj
              !
              ! Llenado de la matriz
              !
              bucle_residuo_direccion_x: do ii = 2, mi-1
                 !
                 call residuo_u(deltaxu,deltayu,deltaxp,&
                      &deltayv,fexp,feyp,fexu,gamma_momen,&
                      &fuente_con_u,fuente_lin_u,&
                      &u,u_ant,v,&
                      &temp,pres,Ri,dt,rel_vel,&
                      &Resu,&
                      &ii,jj)
                 !
              end do bucle_residuo_direccion_x
              !
           end do bucle_residuo_direccion_y
           !$omp end target teams distribute parallel do
           !
           !--------------------------------
           !
           ! residuo del algoritmo
           !
           residuo = 0.0_DBL
           !$omp target teams distribute parallel do collapse(2) reduction(max:residuo)
           calculo_maximo_residuou: do jj=2, nj
              do ii = 2, mi-1
                 residuo = max(residuo, dabs(Resu(ii,jj)) )
                 ! residuo = residuo + Resu(ii,jj)*Resu(ii,jj)
              end do
           end do calculo_maximo_residuou
           ! residuo = sqrt(residuo)
           !$omp end target teams distribute parallel do
           !
           if ( maxbo<conv_paso .and. residuo<conv_resi)then
              !
              iter_simple = iter_simple + 1
              write(102,*) 'SIMPLE', iter_simple, maxbo, residuo
              !
              iter_simple = 0
              exit
              !
           else if ( iter_simple > iter_simple_max ) then
              !
              iter_simple = iter_simple + 1
              write(102,*) 'SIMPLE', iter_simple, maxbo, residuo
              !
              iter_simple = 0
              exit
              !
           else
              !
              iter_simple = iter_simple + 1
              write(102,*) 'SIMPLE', iter_simple, maxbo, residuo
              !
           end if
           !
        end do ALGORITMO_SIMPLE  !final del algoritmo SIMPLE
        !
        ! -------------------------------------
        !
        ! Escritura de mensajes y postprocesos
        !
        itera  = itera + 1
        tiempo = tiempo_inicial+itera*dt
        !
        ! --------------------------------------------
        !
        ! Se calculan promedios en perfiles definidos.
        ! Se hacen tres por cada unidad de tiempo 
        !
        if( mod(itera, ceiling( 1._DBL/(3._DBL*dt) ) ) == 0 ) then
           !
           ! Mensaje de avance de simulaci'on
           !
           write(*,106) tiempo,maxbo,residuo
           !
           call postpro_promedio( tiempo, temp, u, v )
           !
        end if
        !
        ! ********************************************************
        !
        ! Se actualizan los  arreglos para paso de tiempo anterior
        !
        !$omp target teams distribute parallel do collapse(2)
        do jj=1, nj+1
           do ii=1, mi+1
              temp_ant(ii,jj) = temp(ii,jj)
           end do
        end do
        !$omp end target teams distribute parallel do
        !
        !$omp target teams distribute parallel do collapse(2)
        do jj=1, nj+1
           do ii=1, mi
              u_ant(ii,jj) = u(ii,jj)
           end do
        end do
        !$omp end target teams distribute parallel do
        !
        !$omp target teams distribute parallel do collapse(2)
        do jj=1, nj
           do ii=1, mi+1
              v_ant(ii,jj) = v(ii,jj)
           end do
        end do
        !$omp end target teams distribute parallel do
     end do
     !--------------------------------------------
     !
     ! Se cierra la regi\'on paralela de datos
     !
     !$omp end target data
     !
     !*************       termina el paquete de iteraciones
     !*****************************************************
     !*****************************************************
     !
     itera_total = itera_inicial+itera
     millar      = itera_total/(1000*paq_itera)
     centena     = (itera_total-millar*1000*paq_itera)/(100*paq_itera)
     decena      = (itera_total-millar*1000*paq_itera-centena*100*paq_itera)&
          &/(10*paq_itera)
     unidad      = (itera_total-millar*1000*paq_itera-centena*100*paq_itera-&
          &decena*10*paq_itera)/(paq_itera)
     decima      = (itera_total-millar*1000*paq_itera-centena*100*paq_itera-&
          &decena*10*paq_itera-unidad*paq_itera)/(paq_itera/10)
     !
     WRITE(dec,16)decima;16 format(I1)
     WRITE(un,16) unidad
     WRITE(de,16) decena
     WRITE(ce,16) centena
     WRITE(m,16)  millar
     DO jj = 2, nj
        DO ii = 2, mi
           uf(ii,jj) = (u(ii,jj)+u(ii-1,jj))/2._DBL
        END DO
        vf(1,jj)    = v(1,jj)
        vf(mi+1,jj) = v(mi+1,jj)
     END DO
     DO ii = 2, mi
        DO jj = 2, nj
           vf(ii,jj) = (v(ii,jj)+v(ii,jj-1))/2._DBL
        END DO
     END DO
     DO jj = 1, nj+1
        uf(1,jj)    = u(1,jj)
        uf(mi+1,jj) = u(mi,jj)
     END DO
     DO ii = 1, mi+1
        vf(ii,1)    = v(ii,1)
        vf(ii,nj+1) = v(ii,nj)
     END DO
     !************************************
     !
     ! Mensaje de paquete de iteraciones
     !
     write(*,*)' '
     WRITE(*,104) itera_total
     write(*,*)' '
     ! WRITE(*,105) MAXVAL(ABS(Restemp)),MAXVAL(ABS(Resv))
     ! WRITE(*,*)' '
     !********************************
     !*** Formato de escritura dat ***
     !--------------------------------
     if( postprocesar )then
        !
        OPEN(unit=2,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rec)//'u.dat')
        WRITE(2,*) placa_min,placa_max,itera_total,ao
        DO jj = 1, nj+1
           DO ii = 1, mi
              WRITE(2,form24) xu(ii),yp(jj),u(ii,jj)
           END DO
        END DO
        CLOSE(unit=2)
        !
        OPEN(unit=3,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rec)//'v.dat')
        WRITE(3,*) placa_min,placa_max,itera_total,ao
        DO jj = 1, nj
           DO ii = 1, mi+1
              WRITE(3,form24) xp(ii),yv(jj),v(ii,jj)
           END DO
        END DO
        CLOSE(unit=3)
        !
        OPEN(unit=4,file='out_n'//trim(njc)//'m'//trim(mic)//'_R'//trim(Rec)//'p.dat')
        WRITE(4,*) placa_min,placa_max,itera_total,ao
        DO jj = 1, nj+1
           DO ii = 1, mi+1
              WRITE(4,form25) xp(ii),yp(jj),temp(ii,jj),pres(ii,jj)
           END DO
        END DO
        CLOSE(unit=4)
        !
        ! *************************************
        ! *** Formato escritura VTK ***********
        !
        archivo = 'n'//trim(njc)//'m'//trim(mic)//'R'//trim(Rec)//'/t_'//m//ce//&
             &de//un//dec//'.vtk'
        call postproceso_bin(xu,yv,xp,yp,u,v,pres,temp,b_o,Rec)
        call postproceso_vtk(xp,yp,uf,vf,pres,temp,b_o,archivo)
        !
     end if ! Postprocesar
     !
  end do !*********** final del repetidor principal
  !
  ! Cierre de archivos de postproceso
  !
  call finaliza_promedio_perfil()
  !
  ! Formatos
  !
100 format(1X,A,'kind= ',I2,', Precision= ',I2,' Rango= ',I3)
101 format(1X,'R =',A,', Pr=',F8.3', Ri=',F8.3', rel_pres=',F5.2', rel_vel=',F5.2)
102 format(1X,'Iteracion inicial=',I7,', mi=',I5,', nj=',I5)
103 format(1X,'N_Izq=',D23.15,', N_Der=',D23.15)
104 format(1X,"Ixchel2D: iter_total=", I7)
105 format(1X,'Res_T=',D23.15,', Res_v=',D23.15)
106 format(1X,"Ixchel2D: tiempo= ",E10.4,", maxbo= ",E9.3,", res_u= ",&
         &E9.3 )
  !
end program IXCHEL2D

