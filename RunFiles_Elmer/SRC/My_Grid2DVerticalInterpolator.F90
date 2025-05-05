!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *  
! *  
! *  This is an adaptation of the Elmer/Ice solver Grid2DInterpolator
! *  to produce a vertical 2d interpolation of any field.
! *  It war originally coded by J. Brondex (Jan 2025) to interpolate
! *  2010 temperature measurements from boreholes 2, 4, 5, 10, 13, 17, 18,
! *  Assuming no lateral variaiton (along y) of the field
! * 
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: J. Brondex 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Interpolate data given on a regular 2D regular grid in an ASCII file (x z Value)
!    in the mesh nodes using bilinear interpolation
!    The data are ordered such that   
!    x1 z1 val11
!    x2 z1 val21
!    ...
!    xn z1 valn1
!    x1 z2 val12
!    ...
!    xn zn valnn 
!    
!    The grid is described by giving:
!    (x0, z0) the left-bottom corner coordinate
!    (lx, lz) the x and z lengths of the covered domain
!    (Nx, Nz) the number of cells in x and z directions 
!    No data are given by -9999 with a tolerance of 0.001
!    These can be over-ridden in the sif by 'no data' and 'no data tol'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE My_Grid2DVerticalInterpolator( Model,Solver,dt,TransientSimulation )

   USE DefUtils

   IMPLICIT NONE
   TYPE(Solver_t), TARGET :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation

   TYPE(ValueList_t), POINTER :: Params
   TYPE(Variable_t), POINTER :: Var
   REAL(KIND=dp), POINTER :: Values(:)
   INTEGER, POINTER :: Perm(:)

   REAL(KIND=DP) :: Rmin, Rmax
   REAL(KIND=DP) :: x, z, Tinterp, x0, z0, lx, lz, dx, dz
   REAL(KIND=DP), ALLOCATABLE :: xb(:), zb(:), Tb(:), xbaux(:), zbaux(:), Tbaux(:)
   REAL(KIND=dp) :: noDataVal, noDataTol, posTol
   REAL(KIND=dp), PARAMETER :: noDataValDefault = -9999.0, noDataTolDefault = 0.001, posTolDefault= 1.0D-06

   INTEGER,parameter :: io=20
   INTEGER :: ok, Nx, Nz, Nb, Nbaux, OutNode
   INTEGER :: i, j, k, l, kmin, NoVar

   CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, DataF
   CHARACTER(LEN=MAX_NAME_LEN) :: Name, FName, ParaName
   CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='My_Grid2DVerticalInterpolator'

   LOGICAL :: GotVar, Found, InvertOrder, FillIn, UnFoundFatal=.TRUE.

   NULLIFY(Params,Var,Values,Perm)

   Params => GetSolverParams()

   ! Read variable to initialize and Data
   NoVar=0
   GotVar=.True.

   DO WHILE(GotVar)
      NoVar = NoVar + 1
      WRITE (Name,'(A,I0)') 'Variable ',NoVar

      VariableName = ListGetString( Params, TRIM(Name), GotVar )
      IF (.NOT.GotVar) EXIT

      Var => VariableGet(Model %  Mesh % Variables, VariableName,UnFoundFatal=UnFoundFatal )
      Values => Var % Values
      Perm => Var % Perm

      WRITE (FName,'(A,I0,A)') 'Variable ',NoVar,' Data File'
      DataF = ListGetString( Params, TRIM(FName), Found, UnFoundFatal )

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' x0'
      x0 = ListGetConstReal( Params, TRIM(ParaName), Found,UnFoundFatal=UnFoundFatal)

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' z0'
      z0 = ListGetConstReal( Params, TRIM(ParaName), Found,UnFoundFatal=UnFoundFatal )
            
      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' lx'
      lx = ListGetConstReal( Params, TRIM(ParaName), Found,UnFoundFatal=UnFoundFatal )
            
      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' lz'
      lz = ListGetConstReal( Params, TRIM(ParaName), Found,UnFoundFatal=UnFoundFatal )

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' Nx'
      Nx = ListGetInteger( Params, TRIM(ParaName), Found ,UnFoundFatal=UnFoundFatal)

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' Nz'
      Nz = ListGetInteger( Params, TRIM(ParaName), Found ,UnFoundFatal=UnFoundFatal)

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' Invert'
      InvertOrder = GetLogical( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) THEN
         InvertOrder = .FALSE.
      END IF
      IF (InvertOrder) THEN
         WRITE(message,'(A,A,I0)')'Inverting order (row major) for variable ', 'Variable ',NoVar
         CALL INFO(Trim(SolverName),Trim(message),Level=1)
      END IF

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' Fill'
      FillIn = GetLogical( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) THEN
         FillIn = .FALSE.
      END IF
      IF (FillIn) THEN
         WRITE(message,'(A,A,I0)')'Filling empty entries for ', 'Variable ',NoVar
         CALL INFO(Trim(SolverName),Trim(message),Level=1)
      END IF

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' no data'
      noDataVal = ListGetConstReal( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) then
         noDataVal = noDataValDefault
         WRITE(message,'(A,A,A,e15.8)')'Keyword <',Trim(ParaName), & 
              '> not found, using default ',noDataValDefault
         CALL INFO(SolverName, Message, Level=3)
      END IF

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' no data tol'
      noDataTol = ListGetConstReal( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) then
         noDataTol = noDataTolDefault
         WRITE(message,'(A,A,A,e15.8)')'Keyword <',Trim(ParaName), & 
              '> not found, using default ',noDataTolDefault
         CALL INFO(SolverName, Message, Level=3)
      END IF

      WRITE (ParaName,'(A,I0,A)') 'Variable ',NoVar,' position tol'
      posTol = ListGetConstReal( Params, TRIM(ParaName), Found )
      IF (.NOT.Found) then
         posTol = posTolDefault
         WRITE(message,'(A,A,A,e15.8)')'Keyword <',Trim(ParaName), & 
              '> not found, using default ',posTolDefault
         CALL INFO(SolverName, Message, Level=3)
      END IF

      OPEN(unit = io, file = TRIM(DataF), status = 'old',iostat = ok)

      IF (ok /= 0) THEN
         WRITE(message,'(A,A)') 'Unable to open file ',TRIM(DataF)
         CALL FATAL(Trim(SolverName),Trim(message))
      END IF
            
      Nb = Nx*Nz 
          
      ALLOCATE(xb(Nb), zb(Nb), Tb(Nb), xbaux(Nb), zbaux(Nb), Tbaux(Nb))

      ! read data
      DO i = 1, Nb 
         READ(io,*,iostat = ok, end=100) xbaux(i), zbaux(i), Tbaux(i)
      END DO
100   Nbaux = Nb - i
      IF (Nbaux > 0) THEN
         WRITE(message,'(I0,A,I0,A,A)') Nbaux,' out of ',Nb,' datasets in file ', TRIM(DataF)
         CALL INFO(Trim(SolverName),Trim(message))         
      END IF
      CLOSE(io)

      ! Make some verifications and - in case - manipulation 
      !on the DEM structure
      dx = lx / (Nx-1.0)
      dz = lz / (Nz-1.0)
      k = 0 
      l = 0
      IF (.NOT.InvertOrder) THEN
         DO j = 1, Nz
            z = z0 + dz*(j-1)
            DO i = 1, Nx 
               k = k + 1
               x = x0 + dx*(i-1)
               IF (.NOT.FillIn) THEN
                  xb(k) = xbaux(k)
                  zb(k) = zbaux(k)
                  Tb(k) = Tbaux(k)
                  IF ((ABS(x-xbaux(k))>posTol*dx).OR.(ABS(z-zbaux(k))>posTol*dz)) THEN
                     
                     WRITE(Message,'(A,A)')'Structure of the DEM is not conforming to what is given in the sif for ',TRIM(FName) 
                     CALL INFO(SolverName, Message, Level=1)
                     WRITE(Message,'(A,i4,A,i4,A,e15.8,2x,e15.8,A,e15.8,2x,e15.8,A)') &
                          'Variable', NoVar, ': Found that point ',k,&
                          ' coordinate is (',xbaux(k),zbaux(k),&
                          '), whereas it should be (',x,z,')' 
                     CALL FATAL(SolverName, Message)                      
                  END IF
               ELSE
                  IF ((ABS(x-xbaux(l+1))>posTol*dx).OR.(ABS(z-zbaux(l+1))>posTol*dz)) THEN
                     xb(k) = x
                     zb(k) = z
                     Tb(k) = noDataVal ! setting to NaN                   
                  ELSE
                     l=l+1
                     xb(k) = xbaux(l)
                     zb(k) = zbaux(l)
                     Tb(k) = Tbaux(l)
                  END IF
               END IF
            END DO
         END DO
      ELSE ! inverse order
         DO i = 1, Nx 
            x = x0 + dx*(i-1) 
            DO j = 1, Nz
               k = k + 1
               z = z0 + dz*(j-1)

               IF (.NOT.FillIn) THEN
                  xb((j-1)*Nx + i) = xbaux(k)
                  zb((j-1)*Nx + i) = zbaux(k)
                  Tb((j-1)*Nx + i) = Tbaux(k)
                  IF ((ABS(x-xb((j-1)*Nx + i))>posTol*dx).OR.(ABS(z-zb((j-1)*Nx + i))>posTol*dz)) THEN
                     
                     WRITE(Message,'(A,A)')'Structure of the DEM is not conforming to what is given in the sif for ',TRIM(FName) 
                     CALL INFO(SolverName, Message, Level=1)
                     WRITE(Message,'(A,i4,A,i4,A,e15.8,2x,e15.8,A,e15.8,2x,e15.8,A)') &
                          'Variable', NoVar, ':Found that point ',k,&
                          ' coordinate is (',xb((j-1)*Nx),zb((j-1)*Nx + i),'),&
                          whereas 3 it should be (',x,z,')' 
                     CALL FATAL(SolverName, Message)                      
                  END IF
               ELSE
                  IF ((ABS(x-xbaux(l+1))>posTol*dx).OR.(ABS(z-zbaux(l+1))>posTol*dz)) THEN
                     xb((j-1)*Nx + i) = x
                     zb((j-1)*Nx + i) = z
                     Tb((j-1)*Nx + i) = noDataVal ! setting to NaN
                  ELSE
                     l=l+1
                     xb((j-1)*Nx + i) = xbaux(l)
                     zb((j-1)*Nx + i) = zbaux(l)
                     Tb((j-1)*Nx + i) = Tbaux(l)
                  END IF
               END IF

            END DO
         END DO
      END IF

      OutNode = 0
      Rmax = 0.0
      DO i=1,Model % Mesh % NumberOfNodes
         x = Model % Mesh % Nodes % x(i)
         z = Model % Mesh % Nodes % z(i)
         Rmin = 0.0
         CALL InterpolateDEM(x,z,xb,zb,Tb,Nx,Nz,x0,z0,lx,lz,Rmin,Tinterp,noDataVal,noDataTol)
         if ( perm(i) .eq. 0 ) CYCLE
         Values(Perm(i)) = Tinterp
         IF (Rmin > 0.0) THEN
            OutNode = OutNode + 1
            IF (Rmin > Rmax) Rmax = Rmin
         END IF
      END DO
          
      ! Give information on the number of Nodes which are outside of the
      ! DEM domain
      IF (OutNode > 0) THEN
         WRITE( Message, '(I0,A,A)' )OutNode,' nodes where found outside of &
                 the DEM domain in ',TRIM(DataF)
         CALL Info( TRIM(SolverName), Message, Level=3 )
         WRITE( Message, '(A,e15.8)' )'The farthest DEM point used to evaluate & 
                 the nodal value was: ', Rmax
         CALL Info( TRIM(SolverName), Message, Level=3 )
      END IF
            
      DEALLOCATE(xb, zb, Tb, xbaux, zbaux, Tbaux)
   END DO

   CALL INFO(Trim(SolverName), '----------ALL DONE----------',Level=5)

END SUBROUTINE My_Grid2DVerticalInterpolator


!!!!!!!!!!!!!!!!!!!
! Subroutine InterpolateDEM
!!------------------------------------------------------------------------------!!
SUBROUTINE InterpolateDEM (x, z, xb, zb, Tb, Nbx, Nbz, xb0, zb0, lbx, lbz, Rmin, Tinterp, noDataVal, noDataTol)
  USE DefUtils
  IMPLICIT NONE
  REAL(KIND=dp),INTENT(IN) :: noDataVal, noDataTol
  INTEGER :: imin, Npt, t
  INTEGER :: NMAX, i, j, Nb, Nbx, Nbz, ib, ix, iz
  REAL(KIND=dp) :: x, z, Tinterp, xb0, zb0, x1, x2, z1, z2, Ti(2,2) 
  REAL(KIND=dp) :: R, Rmin, lbx, lbz, dbx, dbz
  REAL(KIND=dp) :: xb(Nbx*Nbz), zb(Nbx*Nbz), Tb(Nbx*Nbz)       

  ! Find interpolated T for that point from the Temperature map 
  dbx = lbx / (Nbx-1.0)
  dbz = lbz / (Nbz-1.0)
  Nb = Nbx*Nbz

  ix = INT((x-xb0)/dbx)+1
  iz = INT((z-zb0)/dbz)+1
  ib = Nbx * (iz - 1) + ix

  ! if we are already at the end of the domain then collapse the 2 by 2 interpolation 
  ! square to just 2 points at the end of the domain (else we get interpolation involving 
  ! points at the beginning of the domain).  This comment refers to the x direction.
  IF (MOD(ib,Nbx) .eq. 0.0) THEN
     Ti(2,1) = noDataVal
     Ti(2,2) = noDataVal
  ELSE
     IF ( (ib+1).gt.size(Tb) ) THEN
        Ti(2,1) = noDataVal
     ELSE
        Ti(2,1) = Tb(ib+1)
     END IF
     IF ( (ib+Nbx+1).gt.size(Tb) ) THEN 
        Ti(2,2) = noDataVal
     ELSE
        Ti(2,2) = Tb(ib + Nbx + 1)
     END IF
  END IF

  x1 = xb(ib)
  IF ( (ib+1).gt.size(xb) ) THEN 
     x2 = noDataVal
  ELSE
     x2 = xb(ib+1)
  END IF

  z1 = zb(ib)
  IF ( (ib+Nbx).gt.size(zb) ) THEN 
     z2 = noDataVal
  ELSE
     z2 = zb(ib + Nbx)
  END IF
  
  IF ( (ib).gt.size(Tb) ) THEN
     Ti(1,1) = noDataVal
  ELSE  
     Ti(1,1) = Tb(ib)
  END IF
  IF ( (ib+Nbx).gt.size(Tb) ) THEN
     Ti(1,2) = noDataVal
  ELSE
     Ti(1,2) = Tb(ib + Nbx)
  END IF

  IF ( (isNoData(Ti(1,1))).OR. &
       (isNoData(Ti(1,2))).OR. &
       (isNoData(Ti(2,1))).OR. &
       (isNoData(Ti(2,2))) ) THEN

     IF ( (isNoData(Ti(1,1))).AND. &
          (isNoData(Ti(1,2))).AND. &
          (isNoData(Ti(2,1))).AND. &
          (isNoData(Ti(2,2))) ) THEN

        ! Find the nearest point available if all neighbouring points have noData
        Rmin = 9999999.0
        DO i=1, Nb
           IF (.NOT.isNoData(Tb(i))) THEN
              R = SQRT((x-xb(i))**2.0+(z-zb(i))**2.0)
              IF (R<Rmin) THEN
                 Rmin = R
                 imin = i
              END IF
           END IF
        END DO
        Tinterp = Tb(imin)

     ELSE
        ! Mean value over the available data if only some points have noData
        Tinterp = 0.0
        Npt = 0
        DO i=1, 2
           DO J=1, 2
              IF (.NOT. isNoData(Ti(i,j))) THEN 
                 Tinterp = Tinterp + Ti(i,j)
                 Npt = Npt + 1
              END IF
           END DO
        END DO
        Tinterp = Tinterp / Npt
     END IF
  ELSE
     ! linear interpolation is only carried out if all 4 neighbouring points have data.
     Tinterp = (Ti(1,1)*(x2-x)*(z2-z)+Ti(2,1)*(x-x1)*(z2-z)+Ti(1,2)*(x2-x)*(z-z1)+Ti(2,2)*(x-x1)*(z-z1))/(dbx*dbz)      
  END IF


CONTAINS

  LOGICAL FUNCTION isNoData(val)

    IMPLICIT NONE
    REAL(KIND=dp),INTENT(IN) :: val

    IF ((val .GT. noDataVal-noDataTol) .AND. (val .LT. noDataVal+noDataTol)) THEN
       isNoData = .TRUE.
    ELSE
       isNoData = .FALSE.
    END IF

    RETURN 

  END FUNCTION isNoData

END SUBROUTINE InterpolateDEM

