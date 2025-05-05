! *****************************************************************************
!> Export in a file for each time step the Min, Mean and Max surface vertical 
!> velocity above the cavity
SUBROUTINE ExportVelo( Model, Solver, dt, TransientSimulation )
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material, Params
  TYPE(variable_t), POINTER :: TimeVar
  TYPE(Nodes_t), SAVE :: Nodes

  TYPE(Variable_t), POINTER :: FlowVariable
  REAL(KIND=dp), POINTER :: FlowValues(:)
  INTEGER, POINTER :: FlowPerm(:), NodeIndexes(:)

  LOGICAL :: stat, FirstTime = .TRUE., GotIt

  INTEGER :: i, j, n, M, t, istat, Nac, nodenumber, DIM

  REAL(KIND=dp) :: xc, yc, Rx, Ry, x, y, MeanW, MinW, MaxW, time, t0 
  REAL(KIND=dp), ALLOCATABLE :: Wac(:) 
  INTEGER, ALLOCATABLE :: NodesAboveCavity(:), NodesAC(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'ExportVelo', FileName

  SAVE NodesAboveCavity, NodesAC, Wac, t0, FirstTime, SolverName, FileName  
  SAVE xc, yc, Rx, Ry, NodeIndexes, DIM, Nac
  !------------------------------------------------------------------------------

! At the first time, store the node number above the cavity
  IF (FirstTime) THEN
    FirstTime = .FALSE.
    xc = 948000.0_dp
    yc = 2105060.0_dp
    Rx = 19.4235_dp
    Ry = 51.43_dp
    DIM = CoordinateSystemDimension()
    n = Solver % Mesh % MaxElementNodes
    M = Model % Mesh % NumberOfNodes

    ALLOCATE(NodesAC(M), Nodes % x(n), Nodes % y(n), Nodes % z(n), STAT=istat)

    IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
    END IF
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )

!! Get name of file in which output will be written
    Params => GetSolverParams()
    FileName = ListGetString(Params, 'Filename', GotIt )
    IF (.NOT. GotIt) THEN
      CALL FATAL( SolverName, 'Name of output file must be provided as "Filename" in solver section.' )
    END IF

! loop over active elements
! locate the nodes above the cavity
   NodesAC = 0
    DO t = 1, Solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      !IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      n = GetElementNOFNodes( Element )
      CALL GetElementNodes( Nodes )
      NodeIndexes => Element % NodeIndexes

      DO i=1,n
        nodenumber = NodeIndexes(i)
        x = Nodes % x(i) 
        y = Nodes % y(i) 
        x = (xc - x) / Rx
        y = (yc - y) / Ry
        
        IF ((x*x+y*y).LE.1.0_dp) THEN
          IF (NodesAC(nodenumber)==0) NodesAC(nodenumber)=1
        END IF
      END DO
    END DO
    Nac = COUNT(NodesAC>0)

    ALLOCATE(NodesAboveCavity(Nac), Wac(Nac), STAT=istat)
    IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
    END IF
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )

! fill a vector with the Nac node numbers of the nodes above the cavity
    j = 0
    DO i=1,M
      IF (NodesAC(i)>0) THEN
        j = j + 1
        NodesAboveCavity(j)=i
      END IF
    END DO

    Material => GetMaterial( Element )
    t0 = GetConstReal(Material, 'StartDay', GotIt)
    IF (.Not.GotIt) CALL FATAL( SolverName, 'StarDay not given' )

  END IF ! first time

! Get the velocity solution
   FlowVariable => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
   IF ( ASSOCIATED( FlowVariable ) ) THEN
      FlowPerm    => FlowVariable % Perm
      FlowValues  => FlowVariable % Values
   ELSE
      CALL Info(SolverName, &
       & 'No variable for velocity associated.',Level=4)
   END IF

   Wac(1:Nac) = FlowValues((DIM+1)*(FlowPerm(NodesAboveCavity(1:Nac))-1) + 3)
! In mm/day   
   MinW = Minval(Wac)*1000.0/365.25
   MaxW = MaxVal(Wac)*1000.0/365.25
   MeanW = SUM(Wac(1:Nac))/Nac*1000.0/365.25

   write(*,*)'Mean, min, max, partition',MeanW,MinW,MaxW,ParEnv % myPe

   TimeVar => VariableGet( Model % Variables,'Time')
   time = TimeVar % Values(1) 

! In days t=1 01/07/2011
   time = t0 + time*365.25_dp 
!    time = time*365.25_dp - 309.0_dp 

   OPEN(10, file = Filename, Access = 'APPEND')
   WRITE(10,1000)time, MeanW, MinW, MaxW, Nac, ParEnv % myPe 
   CLOSE(10)

1000 Format(4(ES20.11E3),2(i6))

    CALL INFO( SolverName , 'Done')
    

!------------------------------------------------------------------------------
END SUBROUTINE ExportVelo 
!------------------------------------------------------------------------------


