!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! USF_Damage.f90
! Computes the evolution of damage and computes the Enhancement factor of the Glen's law
! as a function of damage evolution 
!
! Last modif : 31 October 2024
! by J. Brondex: can now deal with various damage criterion: Max. princip. stress, von Mises, Hayurst, and Coulomb
!
! (1) Enhancement Factor 
! Need some inputs in the sif file.
! Parameters: 
! Glen Exponent
!             
!
! (2) SourceDamage
! Need some inputs in the sif file.
! Parameters: 
! Damage Enhancement Factor
! Damage Parameter sigmath
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! EnhancementFactor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION EnhancementFactor ( Model, nodenumber, D) RESULT(E)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(ValueList_t), POINTER :: Material
   TYPE(Solver_t), TARGET :: Solver
   REAL(KIND=dp) :: D, E, n  
   INTEGER :: nodenumber
   LOGICAL :: FirstTime=.TRUE., GotIt

   SAVE FirstTime, n 

   IF (FirstTime) THEN
   FirstTime = .False.
    
      Material => GetMaterial()
      n = GetConstReal( Material, 'Glen Exponent', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Glen Exponent not found. &
              &Setting to 3.0'
         CALL INFO('Damage Enhancement Factor', Message, level=2)
         n = 3.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'n = ', n 
         CALL INFO('Damage EnhancementFactor', Message, level=2)
      END IF
   END IF

   !Damage must remain <0 to avoid infinit E
   D = MIN(0.7_dp, D)
   D = MAX(0.0_dp, D)
   E = (1.0 - D)**(-n) 
END FUNCTION EnhancementFactor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SourceDamage 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION SourceDamage (Model, nodenumber, D) RESULT(Source)


   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   USE GeneralUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: D, Source          
   INTEGER :: nodenumber
 
   TYPE(Solver_t):: Solver 
   TYPE(ValueList_t), POINTER :: Material, Constants

   TYPE(Variable_t), POINTER :: SigmaVariable
   TYPE(Variable_t), POINTER :: StressVariable, FlowVariable, ChiVariable

   REAL(KIND=dp), POINTER :: SigmaValues(:)
   REAL(KIND=dp), POINTER :: StressValues(:), FlowValues(:), ChiValues(:)

   INTEGER, POINTER :: StressPerm(:), FlowPerm(:), ChiPerm(:), SigmaPerm(:)

   INTEGER :: Ind(3,3), DIM, i, j, indice(3), infor
   REAL (KIND=dp) :: Sig(3,3), SigDev(3,3), EigVect(3,3), EigValues(3),tmp, Sigma(3)
   REAL (KIND=dp) :: SigmaI, SigmaII, SigmaIII, SigmaEq, Chi, B, sigmath, lambdah
   REAL (KIND=dp) :: alphaH, betaH, mu
   REAL (KIND=DP) :: EI(3),Dumy(1),Work(24), sigmath_var, stress_threshold
  ! REAL (KIND=DP) :: u, v, nbrPi, s
   LOGICAL :: GotIt, FirstTime = .TRUE., Cauchy, UnFoundFatal=.TRUE.
   CHARACTER*20 :: USF_Name='SourceDamage'
   CHARACTER(LEN=MAX_NAME_LEN) :: Criterion

   SAVE :: Criterion
   SAVE :: Ind, DIM, sigmath, B, lambdah, alphaH, betaH, mu
   SAVE :: FirstTime, Cauchy   

   !!!!############################## START INITIALIZATION  #######################################!!!!
   IF (FirstTime) THEN
      FirstTime = .FALSE.  
      DIM = CoordinateSystemDimension()

      DO i=1, 3
         Ind(i,i) = i
      END DO
      Ind(1,2) = 4
      Ind(2,1) = 4
      Ind(2,3) = 5
      Ind(3,2) = 5
      Ind(3,1) = 6
      Ind(1,3) = 6

      Constants => GetConstants()
      Material => GetMaterial()

     !Read in material which damage criterion needs to be applied: Max Principal Stress, von Mises, Hayurst or Coulomb
      Criterion = GetString( Material, 'Damage Criterion', GotIt ) !Be Ware: If tests ons tring are casse sensitive and Criterion is returned in lower casse (even if written differently in sif)
      IF (.NOT. GotIt) THEN
         CALL FATAL('Damage Source', 'String "Damage Criterion" not found in Material. Must be set to "Max Principal Stress", "von Mises", "Hayurst", or "Coulomb"')
      ELSE
         WRITE(Message,'(A,A)') 'Damage Criterion = ', Criterion
         CALL INFO('Damage Source', Message, level=2)
      END IF 
 
     !For all criterion, read in constants the value of enhancement factor B 
      B = GetConstReal( Constants, 'Damage Enhancement Factor', GotIt )
      IF (.NOT.GotIt) THEN
         CALL FATAL('Damage Source', 'Damage Enhancement Factor B not Found')
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Enhancement Factor = ', B
         CALL INFO('Damage Source', Message, level=2)
      END IF

     !For all criterion, read in constants the value of stress threshold sigmath  
      sigmath = GetConstReal( Constants, 'Damage Parameter sigmath', GotIt )
      IF (.NOT.GotIt) THEN
         CALL FATAL('Damage Source', 'Damage Parameter Sigmath not Found')
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Parameter sigmath = ', sigmath 
         CALL INFO('Damage Source', Message, level=2)
      END IF

     !For all criterion, read in constants the value of healing parameter lambdah  
      lambdah = GetConstReal( Constants, 'Damage Healing Parameter', GotIt)
      IF (.NOT.GotIt) THEN
         lambdah = 0.0_dp
         WRITE(Message,'(A,F10.4)') 'No "Damage Healing Parameter" Given, set to ', lambdah
         CALL INFO('Damage Source', Message, level=2)
      ELSE
         WRITE(Message,'(A,F10.4)') '"Damage Healing Parameter" Given, set to ', lambdah
         CALL INFO('Damage Source', Message, level=2)
      END IF
     
     !Now read parameters specific to some criterion: 
     !alphaH and betaH for Hayurst:
     IF ( Criterion == "hayurst" ) THEN !!!ATTENTION, test on string is case sensitive
         alphaH = GetConstReal( Constants, 'alphaH', GotIt )
         IF (.NOT.GotIt) THEN
            CALL FATAL('Damage Source', 'Hayurst criterion prescribed but "alphaH" not Found in Constants.')
         ELSE
            WRITE(Message,'(A,F10.4)') 'alphaH = ', alphaH
            CALL INFO('Damage Source', Message, level=2)
         END IF
         betaH = GetConstReal( Constants, 'betaH', GotIt )
         IF (.NOT.GotIt) THEN
            CALL FATAL('Damage Source', 'Hayurst criterion prescribed but "betaH" not Found in Constants.')
         ELSE
            WRITE(Message,'(A,F10.4)') 'betaH = ', betaH
            CALL INFO('Damage Source', Message, level=2)
         END IF
     !mu for Coulomb
      ELSE IF ( Criterion == "coulomb") THEN
         mu = GetConstReal( Constants, 'mu', GotIt )
         IF (.NOT.GotIt) THEN
            CALL FATAL('Damage Source', 'Coulomb criterion prescribed but "mu" not Found in Constants.')
         ELSE
            WRITE(Message,'(A,F10.4)') 'mu = ', mu
            CALL INFO('Damage Source', Message, level=2)
         END IF
      END IF

      
     !Cauchy or deviatoric stresses ? -> Always deviatoric when exported in Porous
      Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )
      WRITE(Message,'(A,L1)') 'Cauchy stress tensor computed ? ', Cauchy 
      CALL INFO('Damage Source', Message, level=2)
   END IF ! FirstTime
   !!!!############################## END INITIALIZATION  #######################################!!!!


   !!!! ########## GET ALL NEEDED VARIABLES ######################!!
   ! Get the Stress -> Must be calculated somewhere, i.e. Porous (DeviatoricStress) or ComputeDevStress (Stress)                    
   StressVariable => VariableGet( Model % Variables, 'Stress' )
   IF ( ASSOCIATED( StressVariable ) ) THEN
      StressPerm    => StressVariable % Perm
      StressValues  => StressVariable % Values
   ELSE
      StressVariable => VariableGet( Model % Variables, 'DeviatoricStressDedicSolv' )
      IF ( ASSOCIATED( StressVariable ) ) THEN
         IF ( Cauchy ) THEN
            CALL FATAL('Damage Source', 'Associated Stresses are DeviatoricStress but Cauchy is set to true : Not consistent !')
         ELSE
            StressPerm    => StressVariable % Perm
            StressValues  => StressVariable % Values
         END IF
      ELSE
         CALL FATAL('Damage Source', 'Stresses not associated, Need ComputeDevStress or Porous Solver !')
      END IF
   END IF 
   
   ! Get Chi variable (positive where damage increases) for visualization -> Must be exported somewhere
   ChiVariable => VariableGet( Model % Variables, 'Chi' )
   IF ( ASSOCIATED( ChiVariable ) ) THEN
      ChiPerm    => ChiVariable % Perm
      ChiValues  => ChiVariable % Values
   ELSE
      CALL FATAL('Damage Source', 'Variable Chi not found. Must be exported somewhere !')
   END IF

   ! Get Sigma variable (Principal Stress) for visualization -> Must be exported somewhere
   SigmaVariable => VariableGet( Model % Variables, 'PrincipleStress',UnFoundFatal=UnFoundFatal)
   IF ( ASSOCIATED( SigmaVariable ) ) THEN
      SigmaPerm    => SigmaVariable % Perm
      SigmaValues  => SigmaVariable % Values
   ELSE
      CALL FATAL('Damage Source', 'Variable PrincipleStress not found. Must be exported somewhere !')
   END IF


   ! Get the variables to compute the hydrostatic pressure  
   FlowVariable => VariableGet( Model % Variables, 'Flow Solution' )
   IF ( ASSOCIATED( FlowVariable ) ) THEN
      FlowPerm    => FlowVariable % Perm
      FlowValues  => FlowVariable % Values
   ELSE
      FlowVariable => VariableGet( Model % Variables, 'Porous' )
      IF ( ASSOCIATED( FlowVariable ) ) THEN
         FlowPerm    => FlowVariable % Perm
         FlowValues  => FlowVariable % Values
      ELSE 
         CALL FATAL('Damage Source', 'Flow Solution not associated, Need NS or Porous Solver !')
      END IF
   END IF

   !!!! ########## MAKE SOME MANIPULATIONS ON VARIABLES ######################!!
   !!! Create the stress tensor from the stress vector
   Sig = 0.0
   DO i=1, DIM
      DO j= 1, DIM
         Sig(i,j) =  &
              StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(i,j) )
      END DO
   END DO
   IF (DIM==2) Sig(3,3) = StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(3,3))

  !!! Convert Cauchy stress in deviatoric stress or the other way round 
  ! S = sigma + p
   IF (.NOT.Cauchy) THEN ! If Deviatoric Stress is computed, then, get the
                         ! Cauchy Stress
       SigDev = Sig
       DO i=1,3  
       !!To access pressure: FlowValues((DIM+1)*(FlowPerm(Nodenumber)-1) + (DIM+1))
           Sig(i,i) = SigDev(i,i) - FlowValues((DIM+1)*FlowPerm(Nodenumber))
       END DO
   ELSE ! If the Cauchy Stress is computed, then get the Deviatoric Stress 
       DO i=1,3  
           SigDev(i,i) = Sig(i,i) + FlowValues((DIM+1)*FlowPerm(Nodenumber))
            !write(*,*)Sig(i,1), Sig(i,2), Sig(i,3)
       END DO
   END IF

   ! Compute the principal stresses:

   ! Get the principal stresses by computing the Eigen value of the Cauchy stress tensor
   CALL DGEEV('N','N',3,Sig,3,EigValues,EI,Dumy,1,Dumy,1,Work,24,infor )
   IF (infor.ne.0) &
   CALL FATAL('Compute EigenValues', 'Failed to compute EigenValues') 

   ! Get the eigenvectors (if necessary)
   CALL DGEEV('N','V',3,Sig,3,EigValues,EI,Dumy,1,EigVect,3,Work,24,infor ) 
   IF (infor.ne.0) &
   CALL FATAL('Compute EigenVectors', 'Failed to compute EigenVectors') 
   
   indice = (/(i,i=1,3)/)
   CALL sortd(3,EigValues,indice)
 
   ! 3 Eigen Values for a 3D Stress State (Sigma I, II and III)
   ! This is achieved using the solver ComputeEigenValues.F90
   ! Values are ordered following Ev1 < Ev2 < Ev3 (Sigma I is therefore the larger)
   SigmaI = EigValues(3)
   SigmaII = EigValues(2)
   SigmaIII = EigValues(1)

   ! We save Sigma I, II and III for possible visualisation
   Sigma(1) = SigmaI
   Sigma(2) = SigmaII
   Sigma(3) = SigmaIII
   
   !!!! ########## COMPUTE EQUIVALENT STRESS DEPENDING ON DAMAGE CRITERION ######################!!
   IF ( Criterion == 'max principal stress') THEN
      SigmaEq = SigmaI
   ELSE IF (Criterion == 'von mises') THEN
      SigmaEq = SQRT(0.5_dp * ((SigmaI - SigmaII)**2.0_dp + (SigmaI - SigmaIII)**2.0_dp + (SigmaII - SigmaIII)**2.0_dp))
   ELSE IF (Criterion == 'hayurst') THEN
      SigmaEq = alphaH * SigmaI + betaH * SQRT(0.5_dp * ((SigmaI - SigmaII)**2.0_dp + (SigmaI - SigmaIII)**2.0_dp + (SigmaII - SigmaIII)**2.0_dp)) &
                                - (1.0_dp - alphaH - betaH) * 3.0_dp * FlowValues((DIM+1)*(FlowPerm(Nodenumber)-1) + (DIM+1)) !!!First invariant of Cauchy stress tensor is -3*pressure
   ELSE IF (Criterion == 'coulomb') THEN
      SigmaEq = 0.5_dp * (SigmaI - SigmaIII) - mu * 0.5_dp * (SigmaI + SigmaIII)
   END IF 


   !!!! ########## PUT SOME RANDOM NOISE ON STRESS THRESHOLD ####################!!
   !!Define Pi for the damage treshold
   !nbrPi = 3.141592_dp
   !
   !u = EvenRandom()
   !v = EvenRandom()
   ! 
   !! Determination of the stress threshold   
   !Constants => GetConstants()
   !s = GetConstReal( Constants, 'Dev Tensile Strength Modifier', GotIt) 
   ! IF (.NOT.GotIt) THEN
   !   CALL FATAL('USF_Damage','No "Dev tensile strength modifier" given, set &
   !   to 0.05')
   ! END IF
   !
   !! Get a normal distribution of mean 0 and std dev "s" using two random numbers
   !! "u" and "v" 
   !sigmath_var = ABS(0+s*SQRT(-2.0_dp*LOG((1.0_dp-u)))*COS(2.0_dp*nbrPi*v))
   !stress_threshold = sigmath*(1+sigmath_var)
   
   stress_threshold = sigmath

   !!!############# CALCULATE CHI AND DEDUCE SOURCE TERM ###################!!!
   Chi = (1.0_dp / (1.0_dp - D)) * SigmaEq  - stress_threshold

   ! Save Chi in ChiVariable
   ChiValues(ChiPerm(nodenumber)) = Chi

   !Here we fill up the newly created PrincipalStress variable with convention SigmaI>SigmaII>SigmaIII
   ! Note that principal stress are also reachable through the ComputeEigenValues.F90 solver which manipulate the EigenStress variable (but with convention SigmaIII>SigmaII>SigmaI)
   SigmaValues(DIM*(SigmaPerm(nodenumber)-1)+1) = Sigma(1)
   SigmaValues(DIM*(SigmaPerm(nodenumber)-1)+2) = Sigma(2)
   SigmaValues(DIM*(SigmaPerm(nodenumber)-1)+3) = Sigma(3)
   
   !!!! Below we calculate the damage source from Chi for all damage criterion. This is ok to have same formulation for the 4 damage criterion as:
   !! * SigmaEq is always positive with von Mises so, as formulated here, healing never occur
   !! * SigmaEq with Coulomb is positive if SigmaI > (1+mu)/(1-mu)*SigmaIII which is unlikely for usual values of mu (mu=0.1) -> healing unlikely to happen (very small risk at the very surface)
   !! * SigmaEq can be negative with max principal stress criterion and healing is wanted in that case (means there is healing for compression)
   !! * SigmaEq can be negative with Hayurst criterion -> Means pressure is winning the game despite always positive von Mises contribution and healing is wanted in that case 
   IF ( Chi .GE. 0.0_dp ) THEN !!SigmaI/(1-D) => Sigmath -> Damage increases at rate B*Chi
     Source = B * Chi
   ELSE IF (SigmaEq > 0.0_dp) THEN !! Chi<0 but SigmaEq > 0 -> No damage and no healing
     Source = 0.0_dp
   ELSE !! SigmaI < 0 (compressive stress) -> healing at rate B*Lambdah*Chi (normally only occurs for Hayurst and Max. principal stress criteria)
     Source = B * lambdah * Chi
   END IF
        
END FUNCTION SourceDamage   
