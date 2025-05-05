!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  USER FUNCTIONS FOR TETE ROUSSE CAVITY PROBLEM
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
! User Function Water Pressure Restricted
! Fix some area where the pressure is set to 0 to ovoid unphysical lost of
! contact far from the cavity. 
! At the beginning, the water pressure only apply in a restricted area and then
! can extend with the Mask. The zone x < 947980 i still unreachable by the water
! (frozen bed)
!!------------------------------------------------------------------------------!!
FUNCTION WaterPressureRestricted ( Model, nodenumber, Mask) RESULT(pw)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver

   TYPE(Variable_t), POINTER :: WaterPressureVariable
   REAL(KIND=dp), POINTER :: WaterPressureValues(:)
   INTEGER, POINTER :: WaterPressurePerm(:)

   INTEGER :: nodenumber  
   LOGICAL :: WaterZero = .FALSE.
   REAL(KIND=dp) :: Mask, pw   
   REAL(KIND=dp) :: x, y, z, SeaPressure    
   
   ! Get Water Pressure Variable -> Must be exported somewhere on Body Bottom Surface
   WaterPressureVariable => VariableGet( Model % Variables, 'WaterPressure' )
   IF ( ASSOCIATED( WaterPressureVariable ) ) THEN
      WaterPressurePerm    => WaterPressureVariable % Perm
      WaterPressureValues  => WaterPressureVariable % Values
   ELSE
      CALL FATAL('USF_TR WaterPressureRestricted', 'Variable WaterPressure not found. Must be exported somewhere !')
   END IF

   x = Model % Nodes % x( nodenumber )
   y = Model % Nodes % y( nodenumber )
   z = Model % Nodes % z( nodenumber )

   WaterZero = .FALSE. 
   IF (x<947980.0) THEN
     WaterZero = .TRUE.
   ELSE 
     IF (SQRT((y - 2105060.0)**2.0)>20.0) WaterZero = .TRUE.  
     ! Go back to False is mask is <= 0 (in the cavity)
     IF (Mask < 0.5) WaterZero = .FALSE. 
   END IF

   pw = 0.0
   IF (.Not.WaterZero) pw = SeaPressure(Model, nodenumber, z) 
   
   WaterPressureValues(WaterPressurePerm(nodenumber)) = pw

END FUNCTION WaterPressureRestricted 
   

! Applying water pressure in cavity only. The zone x < 947980 i still unreachable by the water
! (frozen bed)
!!------------------------------------------------------------------------------!!
FUNCTION WaterPressureCavityOnly ( Model, nodenumber, Mask) RESULT(pw)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver

   TYPE(Variable_t), POINTER :: WaterPressureVariable
   REAL(KIND=dp), POINTER :: WaterPressureValues(:)
   INTEGER, POINTER :: WaterPressurePerm(:)

   INTEGER :: nodenumber  
   LOGICAL :: WaterZero = .TRUE.
   REAL(KIND=dp) :: Mask, pw   
   REAL(KIND=dp) :: x, y, z, SeaPressure    
   
   ! Get Water Pressure Variable -> Must be exported somewhere on Body Bottom Surface
   WaterPressureVariable => VariableGet( Model % Variables, 'WaterPressure' )
   IF ( ASSOCIATED( WaterPressureVariable ) ) THEN
      WaterPressurePerm    => WaterPressureVariable % Perm
      WaterPressureValues  => WaterPressureVariable % Values
   ELSE
      CALL FATAL('USF_TR WaterPressureCavityOnly', 'Variable WaterPressure not found. Must be exported somewhere !')
   END IF

   x = Model % Nodes % x( nodenumber )
   y = Model % Nodes % y( nodenumber )
   z = Model % Nodes % z( nodenumber )

   WaterZero = .TRUE. 
   IF (x>947980.0 .AND. Mask < 0.5) THEN
     WaterZero = .FALSE.
   END IF

   pw = 0.0
   IF (.Not.WaterZero) THEN
      pw = SeaPressure(Model, nodenumber, z) 
!      PRINT*, "Node Number =", nodenumber
!      PRINT*, "Water Pressure =", pw
   END IF   

   WaterPressureValues(WaterPressurePerm(nodenumber)) = pw

END FUNCTION WaterPressureCavityOnly    

! Same as above for elastic solver pressure need a minus signs !!!
! 
!!------------------------------------------------------------------------------!!
FUNCTION WaterPressureCavityOnly_Elastic ( Model, nodenumber, Mask) RESULT(pw)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver

   TYPE(Variable_t), POINTER :: WaterPressureVariable
   REAL(KIND=dp), POINTER :: WaterPressureValues(:)
   INTEGER, POINTER :: WaterPressurePerm(:)

   INTEGER :: nodenumber  
   LOGICAL :: WaterZero = .TRUE.
   REAL(KIND=dp) :: Mask, pw   
   REAL(KIND=dp) :: x, y, z, SeaPressure    
   
   ! Get Water Pressure Variable -> Must be exported somewhere on Body Bottom Surface
   WaterPressureVariable => VariableGet( Model % Variables, 'WaterPressure' )
   IF ( ASSOCIATED( WaterPressureVariable ) ) THEN
      WaterPressurePerm    => WaterPressureVariable % Perm
      WaterPressureValues  => WaterPressureVariable % Values
   ELSE
      CALL FATAL('USF_TR WaterPressureCavityOnly', 'Variable WaterPressure not found. Must be exported somewhere !')
   END IF

   x = Model % Nodes % x( nodenumber )
   y = Model % Nodes % y( nodenumber )
   z = Model % Nodes % z( nodenumber )

   WaterZero = .TRUE. 
   IF (x>947980.0 .AND. Mask < 0.5) THEN
     WaterZero = .FALSE.
   END IF

   pw = 0.0
   IF (.Not.WaterZero) THEN
!      PRINT*,'Node of cavity'
      pw = SeaPressure(Model, nodenumber, z) 
!      PRINT*, "Node Number =", nodenumber
!      PRINT*, "Water Pressure =", pw
   END IF

   WaterPressureValues(WaterPressurePerm(nodenumber)) = pw

END FUNCTION WaterPressureCavityOnly_Elastic    
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
! User Function Water Pressure
! Fix some area where the pressure is set to 0 to ovoid unphysical lost of
! contact far from the cavity
!!------------------------------------------------------------------------------!!
FUNCTION WaterPressure ( Model, nodenumber, z) RESULT(pw)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver

   TYPE(Variable_t), POINTER :: WaterPressureVariable
   REAL(KIND=dp), POINTER :: WaterPressureValues(:)
   INTEGER, POINTER :: WaterPressurePerm(:)

   INTEGER :: nodenumber  
   REAL(KIND=dp) :: z, pw   
   REAL(KIND=dp) :: x, y, SeaPressure    
   
   ! Get Water Pressure Variable -> Must be exported somewhere on Body Bottom Surface
   WaterPressureVariable => VariableGet( Model % Variables, 'WaterPressure' )
   IF ( ASSOCIATED( WaterPressureVariable ) ) THEN
      WaterPressurePerm    => WaterPressureVariable % Perm
      WaterPressureValues  => WaterPressureVariable % Values
   ELSE
      CALL FATAL('USF_TR WaterPressure', 'Variable WaterPressure not found. Must be exported somewhere !')
   END IF

   x = Model % Nodes % x( nodenumber )
   IF (x<947980.0) THEN
     pw = 0.0
   ELSE
     pw = SeaPressure(Model, nodenumber, z) 
   END IF
   
   WaterPressureValues(WaterPressurePerm(nodenumber)) = pw

END FUNCTION WaterPressure    

!!!!!!!!!!!!!!!!!!!
! User Function MinZsBottom
! Return the bedrock altitude for a given nodes
!!------------------------------------------------------------------------------!!
FUNCTION MinZsBottom ( Model, nodenumber, znode) RESULT(Zbed)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber  
   REAL(KIND=dp) :: znode, Zbed 
   INTEGER :: imin, Npt, t
   INTEGER :: NMAX, i, j,Nb, Nbx, Nby, ib, ix, iy
   REAL(KIND=dp) :: x, y, z, xb0, yb0, x1, x2, y1, y2, zi(2,2) 
   REAL(KIND=dp) :: R, Rmin, dbx, dby, lbx, lby
   REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE xb, yb, zb

   Nbx = 301
   Nby = 176
   xb0 = 947700.0d0
   yb0 = 2104850.0d0      
   lbx = 600.0
   lby = 350.0

   Nb = Nbx * Nby
   dbx = lbx / (Nbx-1.0)
   dby = lby / (Nby-1.0)

   IF (FirstTime) THEN
        FirstTime = .False.
        
! Open the DEM
! Grid specificity can be read in mesh_input.dat
! 
        OPEN(10,file="../../Data/Topo/MNTBedInf2014.dat")
        
        ALLOCATE(xb(Nb), yb(Nb), zb(Nb))
        READ(10,*)(xb(i), yb(i), zb(i), i=1,Nb)
        CLOSE(10)
   END IF
        
! Compute zbed for that point (x,y)

        x = Model % Nodes % x (nodenumber)
        y = Model % Nodes % y (nodenumber)
         
! Find zbed for that point from the Bedrock MNT 

        ix = INT((x-xb0)/dbx)+1
        iy = INT((y-yb0)/dbx)+1
        ib = Nbx * (iy - 1) + ix
        
        x1 = xb(ib)
        x2 = xb(ib+1)
        y1 = yb(ib)
        y2 = yb(ib + Nbx)
        
        zi(1,1) = zb(ib)
        zi(2,1) = zb(ib+1)
        zi(2,2) = zb(ib + Nbx + 1)
        zi(1,2) = zb(ib + Nbx)
        
        
        IF ((zi(1,1)<-9990.0).OR.(zi(1,2)<-9990.0).OR.(zi(2,1)<-9990.0).OR.(zi(2,2)<-9990.0)) THEN
           IF ((zi(1,1)<-9990.0).AND.(zi(1,2)<-9990.0).AND.(zi(2,1)<-9990.0).AND.(zi(2,2)<-9990.0)) THEN
           ! Find the nearest point avalable
             Rmin = 9999.0
             DO i=1, Nb
               IF (zb(i)>0.0) THEN
                 R = SQRT((x-xb(i))**2.0+(y-yb(i))**2.0)
                 IF (R<Rmin) THEN
                   Rmin = R
                   imin = i
                 END IF
               END IF
             END DO
            zbed = zb(imin)
                        
           ELSE
            ! Mean value over the avalable data
             zbed = 0.0
             Npt = 0
             DO i=1, 2
               DO J=1, 2
                  IF (zi(i,j) > 0.0) THEN 
                     zbed = zbed + zi(i,j)
                     Npt = Npt + 1
                  END IF   
               END DO
             END DO
             zbed = zbed / Npt
             
           END IF
        ELSE
          zbed = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dbx*dby)      
        END IF

END FUNCTION MinZsBottom


!!!!!!!!!!!!!!!!!!!
! User Function MaskCavity
! Mask = -1 for cavity nodes (zb > b), Mask = 1 for bedrock (zb=b)
!!------------------------------------------------------------------------------!!
FUNCTION MaskCavity ( Model, nodenumber, znode) RESULT(Mask)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber  
   REAL(KIND=dp) :: znode, Mask
   INTEGER :: imin, Npt, t
   INTEGER :: NMAX, i, j,Nb, Nbx, Nby, ib, ix, iy
   REAL(KIND=dp) :: x, y, z, xb0, yb0, x1, x2, y1, y2, zi(2,2) 
   REAL(KIND=dp) :: R, Rmin, dbx, dby, lbx, lby, zbed
   REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE xb, yb, zb


   Nbx = 301
   Nby = 176
   xb0 = 947700.0d0
   yb0 = 2104850.0d0      
   lbx = 600.0
   lby = 350.0

   Nb = Nbx * Nby
   dbx = lbx / (Nbx-1.0)
   dby = lby / (Nby-1.0)

   IF (FirstTime) THEN
        FirstTime = .False.
        
! Open the DEM
! Grid specificity can be read in mesh_input.dat
! 
        OPEN(10,file="../../Data/Topo/MNTBedInf2014.dat")
        
        ALLOCATE(xb(Nb), yb(Nb), zb(Nb))
        READ(10,*)(xb(i), yb(i), zb(i), i=1,Nb)
        CLOSE(10)
   END IF
        
! Compute zbed for that point (x,y)

        x = Model % Nodes % x (nodenumber)
        y = Model % Nodes % y (nodenumber)
         
! Find zbed for that point from the Bedrock MNT 

        ix = INT((x-xb0)/dbx)+1
        iy = INT((y-yb0)/dbx)+1
        ib = Nbx * (iy - 1) + ix
        
        x1 = xb(ib)
        x2 = xb(ib+1)
        y1 = yb(ib)
        y2 = yb(ib + Nbx)
        
        zi(1,1) = zb(ib)
        zi(2,1) = zb(ib+1)
        zi(2,2) = zb(ib + Nbx + 1)
        zi(1,2) = zb(ib + Nbx)
        
        
        IF ((zi(1,1)<-9990.0).OR.(zi(1,2)<-9990.0).OR.(zi(2,1)<-9990.0).OR.(zi(2,2)<-9990.0)) THEN
           IF ((zi(1,1)<-9990.0).AND.(zi(1,2)<-9990.0).AND.(zi(2,1)<-9990.0).AND.(zi(2,2)<-9990.0)) THEN
           ! Find the nearest point avalable
             Rmin = 9999.0
             DO i=1, Nb
               IF (zb(i)>0.0) THEN
                 R = SQRT((x-xb(i))**2.0+(y-yb(i))**2.0)
                 IF (R<Rmin) THEN
                   Rmin = R
                   imin = i
                 END IF
               END IF
             END DO
            zbed = zb(imin)
                        
           ELSE
            ! Mean value over the avalable data
             zbed = 0.0
             Npt = 0
             DO i=1, 2
               DO J=1, 2
                  IF (zi(i,j) > 0.0) THEN 
                     zbed = zbed + zi(i,j)
                     Npt = Npt + 1
                  END IF   
               END DO
             END DO
             zbed = zbed / Npt
             
           END IF
        ELSE
          zbed = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dbx*dby)      
        END IF
                  
   
   znode = Model % Nodes % z(nodenumber) 
   IF (znode > Zbed+0.1) THEN
      Mask = -1.0
   ELSE
      Mask = 1.0
   END IF

END FUNCTION MaskCavity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Read the water level in a file
!  and interpolate linearly for the given time 
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION WaterLevel ( Model, nodenumber, time ) RESULT(value)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: value, time

   REAL(KIND=dp) :: t, t0
   TYPE(ValueList_t), POINTER :: Material
   TYPE(Variable_t), POINTER :: Timevar
   REAL(KIND=dp), ALLOCATABLE :: data_value(:)
   REAL(KIND=dp), ALLOCATABLE :: Coordinate(:)
   INTEGER :: Ns, i  
   LOGICAL :: FirstTime = .TRUE., GotIt 
   CHARACTER(LEN=MAX_NAME_LEN) :: FileDataName

   SAVE FirstTime, Ns
   SAVE data_value, Coordinate, t0

   IF (FirstTime) THEN
      FirstTime = .FALSE.
!------------------------------------
! Get the name of the input data file         
!------------------------------------
      FileDataName = '../../Data/WaterLevel/Evol_Niveau2010-2013-moyen.dat'

      OPEN(1,file=TRIM(FileDataName))
      READ(1,*)Ns

      ALLOCATE (Coordinate(Ns), data_value(Ns))

      READ(1,*)(Coordinate(i), data_value(i), i=1,Ns)
      CLOSE(1) 
  
      Material => GetMaterial( Model % CurrentElement)
      t0 = GetConstReal(Material, 'StartDay', GotIt)
      IF (.Not.GotIt) CALL FATAL( 'WaterLevel', 'StarDay not given' )
   ENDIF ! End First Time


!-----------------------
! Interpolate for that time
!-----------------------
! time in year in Elmer
! time in days in the file, t=1 01/07/2011
    Timevar => VariableGet( Model % Variables, 'Time')
    time = TimeVar % Values(1)
    t = t0 + time*365.25
    value = InterpolateCurve(coordinate,data_value,t)
    PRINT*, 'time in days =', t, 'water level =', value
END FUNCTION WaterLevel  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Read the cumulated SMB from the input file
!  and interpolate linearly for the given time and given altitude
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION TRsmbTime ( Model, nodenumber, time ) RESULT(value)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: value, time

   TYPE(ValueList_t), POINTER :: Material
   TYPE(Variable_t), POINTER :: Timevar, TimeStepVar 
   REAL(KIND=dp), ALLOCATABLE :: data_value(:,:)
   REAL(KIND=dp), ALLOCATABLE :: Coordinate(:)
   REAL(KIND=dp) :: Altitude(8), t0, t1, dt, z
   REAL(KIND=dp) :: val0, val1, valt0, valt1, day0
   INTEGER :: Ns, i, j0, j1 
   LOGICAL :: FirstTime = .TRUE., GotIt  
   CHARACTER(LEN=MAX_NAME_LEN) :: FileDataName

   SAVE FirstTime, Ns, Altitude
   SAVE data_value, Coordinate, day0

   IF (FirstTime) THEN
      FirstTime = .FALSE.
!------------------------------------
! Get the name of the input data file         
!------------------------------------
      FileDataName = '../../Data/SMB/SMB_TeteRousse.dat'

      OPEN(1,file=TRIM(FileDataName))
      READ(1,*)Ns

      ALLOCATE (Coordinate(Ns), data_value(Ns,8))

      READ(1,*)(Coordinate(i), data_value(i,1), data_value(i,2), &
       data_value(i,3), data_value(i,4), data_value(i,5), data_value(i,6), &
       data_value(i,7), data_value(i,8), i=1,Ns)
      CLOSE(1) 

      Altitude=(/3125.0,3138.0,3154.0,3169.0,3185.0,3200.0,3226.0,3251.0/)

      Material => GetMaterial( Model % CurrentElement)
      day0 = GetConstReal(Material, 'StartDay', GotIt)
      IF (.Not.GotIt) CALL FATAL( 'WaterLevel', 'StarDay not given' )
   ENDIF ! End First Time

!------------------------------------
! Find the altitude to interpolate to
!------------------------------------
   j1 = 0
   z = Model % Nodes % z ( nodenumber)
   IF (z < Altitude(1)) THEN
     j0 = 1
   ELSE IF (z < Altitude(2)) THEN
     j0 = 1
     j1 = 2
   ELSE IF (z < Altitude(3)) THEN
     j0 = 2
     j1 = 3
   ELSE IF (z < Altitude(4)) THEN
     j0 = 3
     j1 = 4
   ELSE IF (z < Altitude(5)) THEN
     j0 = 4
     j1 = 5
   ELSE IF (z < Altitude(6)) THEN
     j0 = 5
     j1 = 6
   ELSE IF (z < Altitude(7)) THEN
     j0 = 6
     j1 = 7
   ELSE IF (z < Altitude(8)) THEN
     j0 = 7
     j1 = 8
   ELSE 
     j0 = 8
   END IF
!-----------------------
! Interpolate for that time
!-----------------------
    Timevar => VariableGet( Model % Variables, 'Time')
    time = TimeVar % Values(1)
    TimeStepvar => VariableGet( Model % Variables, 'timestep size')
    dt = TimeStepVar % Values(1)
! Time in year in Elmer
! time in days in the file
    t0 = day0 + (time - dt)*365.25_dp
    t1 = day0 + time*365.25_dp

!    PRINT*, 'time (d) lower bound SMB USF=', t0
!    PRINT*, 'time (d) upper bound SMB USF=', t1
! Cumulated SMB at t-dt
    val0 = InterpolateCurve(coordinate,data_value(:,j0),t0)
!    PRINT*, 'SMB on day', t0, 'at altitude', Altitude(j0), 'is ', val0
    IF (j1>0) val1 = InterpolateCurve(coordinate,data_value(:,j1),t0)
!    IF (j1>0) PRINT*, 'SMB on day', t0, 'at altitude', Altitude(j1), 'is ', val1
    IF (j1>0) THEN
      valt0 = val0 + (val1-val0)*(z-Altitude(j0))/(Altitude(j1)-Altitude(j0)) 
!      PRINT*, 'SMB on day', t0, 'at altitude', z, 'is ', valt0
    ELSE
      Valt0 = val0
!      PRINT*, 'SMB on day', t0, 'at altitude', z, 'is ', valt0
    ENDIF 
! Cumulated SMB at t 
    val0 = InterpolateCurve(coordinate,data_value(:,j0),t1)
!    PRINT*, 'SMB on day', t1, 'at altitude', Altitude(j0), 'is ', val0
    IF (j1>0) val1 = InterpolateCurve(coordinate,data_value(:,j1),t1)
!    IF (j1>0) PRINT*, 'SMB on day', t1, 'at altitude', Altitude(j1), 'is ', val1
    IF (j1>0) THEN
      valt1 = val0 + (val1-val0)*(z-Altitude(j0))/(Altitude(j1)-Altitude(j0)) 
!      PRINT*, 'SMB on day', t1, 'at altitude', z, 'is ', valt1
    ELSE
      Valt1 = val0
!      PRINT*, 'SMB on day', t1, 'at altitude', z, 'is ', valt1
    ENDIF 

! Valt1 and Valt0 are Cumulated SMB in cm water eq. 
! SMB between t-dt and t in m/a.ice.eq
    value = (valt1 - valt0) / (0.917 * 100.0 * dt)
   
!    PRINT*, 'SMB in m/a on day', t1, 'at altitude', z, 'is ', value
  ! write(*,*)t0,t1,dt,z,j0,j1,value
END FUNCTION TRsmbTime              

