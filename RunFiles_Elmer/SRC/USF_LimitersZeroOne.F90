FUNCTION ApplyLimiters(Model, n, f) RESULT(g)
USE DefUtils
TYPE(Model_t) :: Model
INTEGER :: n
REAL(KIND=dp) :: f, g

g = MIN( MAX(f,0.0_dp),1.0_dp )

END FUNCTION ApplyLimiters
