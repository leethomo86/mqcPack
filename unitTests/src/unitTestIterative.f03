Program unitTestIterative
  Use iso_fortran_env,Only:int64,real64
  Use MQC_General
  Use MQC_Algebra
  Implicit None

  Integer(Kind=int64),Parameter::Dimension=24_int64,NumberOfRoots=4_int64
  Type(MQC_Matrix)::RealMatrix,ComplexMatrix,DenseVectors,DavidsonVectors
  Type(MQC_Vector)::DenseValues,DavidsonValues,ResidualNorms
  Type(MQC_Explicit_Matrix_Operator)::Operator
  Real(Kind=real64)::Value,ErrorValue,MaxEigenvalueError,MaxResidual
  Complex(Kind=real64)::ComplexValue
  Integer(Kind=int64)::i,j,Iterations
  Logical::Converged

  Call RealMatrix%Init(Dimension,Dimension)
  Do j=1,Dimension
    Do i=j,Dimension
      Value=Sin(Real(3*i+5*j,Kind=real64))+&
        0.1_real64*Cos(Real(7*i-2*j,Kind=real64))
      If(i.eq.j) Value=Value+Real(i,Kind=real64)
      Call RealMatrix%Put(Value,i,j)
      Call RealMatrix%Put(Value,j,i)
    EndDo
  EndDo
  Write(*,'(1x,A)') 'Iterative test: dense real diagonalization'
  Call RealMatrix%Diag(DenseValues,DenseVectors)
  Write(*,'(1x,A)') 'Iterative test: real Davidson'
  Call Operator%Init(RealMatrix)
  Call MQC_Davidson(Operator,NumberOfRoots,DavidsonValues,DavidsonVectors,&
    Tolerance=1.0e-10_real64,MaximumSubspace=12_int64,&
    MaximumIterations=100_int64,Converged=Converged,Iterations=Iterations,&
    ResidualNorms=ResidualNorms)
  Call AssertTrue(Converged,'Real symmetric Davidson calculation did not converge')
  MaxEigenvalueError=0.0_real64
  Do i=1,NumberOfRoots
    ErrorValue=DavidsonValues%At(i)-DenseValues%At(i)
    MaxEigenvalueError=Max(MaxEigenvalueError,Abs(ErrorValue))
  EndDo
  MaxResidual=ResidualNorms%MaxVal()
  Call AssertTrue(MaxEigenvalueError.lt.1.0e-9_real64,&
    'Real symmetric Davidson eigenvalues disagree with dense diagonalization')
  Call AssertTrue(MaxResidual.lt.1.0e-9_real64,&
    'Real symmetric Davidson residual is too large')

  Call ComplexMatrix%Init(Dimension,Dimension,&
    Scalar=(0.0_real64,0.0_real64))
  Do j=1,Dimension
    Do i=j,Dimension
      If(i.eq.j) Then
        ComplexValue=Cmplx(Real(i,Kind=real64)+&
          0.2_real64*Sin(Real(i,Kind=real64)),0.0_real64,Kind=real64)
      Else
        ComplexValue=Cmplx(0.2_real64*Sin(Real(2*i+j,Kind=real64)),&
          0.15_real64*Cos(Real(i-3*j,Kind=real64)),Kind=real64)
      EndIf
      Call ComplexMatrix%Put(ComplexValue,i,j)
      Call ComplexMatrix%Put(Conjg(ComplexValue),j,i)
    EndDo
  EndDo
  Write(*,'(1x,A)') 'Iterative test: dense complex diagonalization'
  Call ComplexMatrix%Diag(DenseValues,DenseVectors)
  Write(*,'(1x,A)') 'Iterative test: complex Davidson'
  Call Operator%Init(ComplexMatrix)
  Call MQC_Davidson(Operator,NumberOfRoots,DavidsonValues,DavidsonVectors,&
    Tolerance=1.0e-10_real64,MaximumSubspace=12_int64,&
    MaximumIterations=100_int64,Converged=Converged,Iterations=Iterations,&
    ResidualNorms=ResidualNorms)
  Call AssertTrue(Converged,'Complex Hermitian Davidson calculation did not converge')
  MaxEigenvalueError=0.0_real64
  Do i=1,NumberOfRoots
    ErrorValue=DavidsonValues%At(i)-DenseValues%At(i)
    MaxEigenvalueError=Max(MaxEigenvalueError,Abs(ErrorValue))
  EndDo
  MaxResidual=ResidualNorms%MaxVal()
  Call AssertTrue(MaxEigenvalueError.lt.1.0e-9_real64,&
    'Complex Hermitian Davidson eigenvalues disagree with dense diagonalization')
  Call AssertTrue(MaxResidual.lt.1.0e-9_real64,&
    'Complex Hermitian Davidson residual is too large')

  Write(*,'(1x,A,I0,A,ES12.5,A,ES12.5)') 'unitTestIterative: PASS (',&
    Iterations,' iterations, eigenvalue error ',MaxEigenvalueError,&
    ', residual ',MaxResidual

Contains

  Subroutine AssertTrue(Condition,Message)
    Logical,Intent(In)::Condition
    Character(Len=*),Intent(In)::Message

    If(.not.Condition) Then
      Write(*,'(1x,A)') 'FAIL: '//Trim(Message)
      Error Stop 1
    EndIf
  End Subroutine AssertTrue

End Program unitTestIterative
