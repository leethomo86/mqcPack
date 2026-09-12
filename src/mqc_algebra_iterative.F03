!> \brief Implementations of the MQC_Algebra matrix-free eigensolver API.
Submodule(MQC_Algebra) MQC_Algebra_Iterative
  Implicit None
Contains

!> Store a dense matrix in the explicit linear-operator adapter.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Explicit_Operator_Init
    If(Size(Matrix,1).ne.Size(Matrix,2)) &
      Call MQC_Error('Explicit linear operator requires a square matrix')
    Operator%Matrix=Matrix
  End Procedure MQC_Explicit_Operator_Init

!> Return the dimension of an explicit linear operator.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Explicit_Operator_Dimension
    Dimension=Size(Operator%Matrix,1)
  End Procedure MQC_Explicit_Operator_Dimension

!> Apply the explicit matrix to a block of trial vectors.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Explicit_Operator_Apply
    Complex(Kind=real64)::Element
    Integer(Kind=int64)::i,j

    If(Size(Vectors,1).ne.Operator%Dimension()) &
      Call MQC_Error('Vector dimension disagrees with explicit operator')
    If(Any(Shape(Products).ne.Shape(Vectors))) &
      Call MQC_Error('Product shape disagrees with input vector block')
    Products=(0.0_real64,0.0_real64)
    Do j=1,Operator%Dimension()
      Do i=1,Operator%Dimension()
        Element=Operator%Matrix%At(i,j)
        Products(i,:)=Products(i,:)+Element*Vectors(j,:)
      EndDo
    EndDo
  End Procedure MQC_Explicit_Operator_Apply

!> Return the diagonal of an explicit matrix operator.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Explicit_Operator_Diagonal
    Type(MQC_Vector)::MatrixDiagonal
    Integer(Kind=int64)::i

    If(Size(Diagonal).ne.Operator%Dimension()) &
      Call MQC_Error('Diagonal dimension disagrees with explicit operator')
    MatrixDiagonal=Operator%Matrix%Diagonal()
    Do i=1,Operator%Dimension()
      Diagonal(i)=MatrixDiagonal%At(i)
    EndDo
  End Procedure MQC_Explicit_Operator_Diagonal

!> Compute the lowest algebraic eigenpairs of a Hermitian linear operator
!> using thick-restarted block Davidson iteration.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Davidson
    Complex(Kind=real64),Allocatable::Basis(:,:),OperatorBasis(:,:),&
      ProjectedVectors(:,:),RitzVectors(:,:),OperatorRitz(:,:),Residuals(:,:),&
      Corrections(:,:),Diagonal(:),WorkVector(:)
    Real(Kind=real64),Allocatable::ProjectedValues(:),RootResiduals(:)
    Type(MQC_Matrix)::ProjectedMatrix,ProjectedEigenvectors
    Type(MQC_Vector)::ProjectedEigenvalues,DiagonalValues,DiagonalOrder
    Complex(Kind=real64)::Denominator
    Real(Kind=real64)::ConvergenceTolerance,NormValue,DenominatorThreshold
    Integer(Kind=int64)::Dimension,MaxSubspace,MaxIterations,SubspaceSize,&
      CorrectionCount,IterationCount,i,j,k,InitialCount,NewCount
    Logical::AllConverged

    Dimension=Operator%Dimension()
    If(NumberOfRoots.lt.1.or.NumberOfRoots.gt.Dimension) &
      Call MQC_Error('Invalid number of roots in MQC_Davidson')
    ConvergenceTolerance=1.0e-8_real64
    If(Present(Tolerance)) ConvergenceTolerance=Tolerance
    MaxSubspace=Min(Dimension,Max(4_int64*NumberOfRoots,20_int64))
    If(Present(MaximumSubspace)) MaxSubspace=Min(Dimension,MaximumSubspace)
    MaxSubspace=Min(Dimension,Max(MaxSubspace,NumberOfRoots+1_int64))
    MaxIterations=100_int64
    If(Present(MaximumIterations)) MaxIterations=MaximumIterations
    If(ConvergenceTolerance.le.0.0_real64.or.MaxIterations.lt.1) &
      Call MQC_Error('Invalid convergence controls in MQC_Davidson')

    Allocate(Basis(Dimension,MaxSubspace),OperatorBasis(Dimension,MaxSubspace),&
      ProjectedVectors(MaxSubspace,MaxSubspace),RitzVectors(Dimension,NumberOfRoots),&
      OperatorRitz(Dimension,NumberOfRoots),Residuals(Dimension,NumberOfRoots),&
      Corrections(Dimension,NumberOfRoots),Diagonal(Dimension),&
      WorkVector(Dimension),ProjectedValues(MaxSubspace),&
      RootResiduals(NumberOfRoots))
    Basis=(0.0_real64,0.0_real64)
    OperatorBasis=(0.0_real64,0.0_real64)
    Call Operator%Diagonal(Diagonal)
    Call DiagonalValues%Init(Dimension)
    Do i=1,Dimension
      Call DiagonalValues%Put(Real(Diagonal(i),Kind=real64),i)
    EndDo
    DiagonalOrder=DiagonalValues%ArgSort()

    InitialCount=NumberOfRoots
    If(Present(InitialVectors)) Then
      If(Size(InitialVectors,1).ne.Dimension.or.&
        Size(InitialVectors,2).lt.NumberOfRoots) &
        Call MQC_Error('Invalid initial-vector dimensions in MQC_Davidson')
      Do j=1,NumberOfRoots
        Do i=1,Dimension
          Basis(i,j)=InitialVectors%At(i,j)
        EndDo
      EndDo
    Else
      Do j=1,NumberOfRoots
        i=DiagonalOrder%At(j)
        Basis(i,j)=(1.0_real64,0.0_real64)
      EndDo
    EndIf
    Call MQC_Iterative_Orthonormalize_Block(Basis,InitialCount)
    If(InitialCount.lt.NumberOfRoots) &
      Call MQC_Error('Initial Davidson vectors are linearly dependent')
    SubspaceSize=InitialCount
    Call Operator%Apply(Basis(:,1:SubspaceSize),OperatorBasis(:,1:SubspaceSize))
    AllConverged=.False.

    Do IterationCount=1,MaxIterations
      Call MQC_Iterative_Projected_Matrix(Basis(:,1:SubspaceSize),&
        OperatorBasis(:,1:SubspaceSize),ProjectedMatrix)
      Call ProjectedMatrix%Diag(ProjectedEigenvalues,ProjectedEigenvectors)
      Do i=1,SubspaceSize
        ProjectedValues(i)=ProjectedEigenvalues%At(i)
        Do j=1,SubspaceSize
          ProjectedVectors(j,i)=ProjectedEigenvectors%At(j,i)
        EndDo
      EndDo
      RitzVectors=MatMul(Basis(:,1:SubspaceSize),&
        ProjectedVectors(1:SubspaceSize,1:NumberOfRoots))
      OperatorRitz=MatMul(OperatorBasis(:,1:SubspaceSize),&
        ProjectedVectors(1:SubspaceSize,1:NumberOfRoots))
      CorrectionCount=0
      Do k=1,NumberOfRoots
        Residuals(:,k)=OperatorRitz(:,k)-ProjectedValues(k)*RitzVectors(:,k)
        RootResiduals(k)=Sqrt(Real(Dot_Product(Residuals(:,k),Residuals(:,k)),&
          Kind=real64))
        If(RootResiduals(k).le.ConvergenceTolerance) Cycle
        WorkVector=-Residuals(:,k)
        DenominatorThreshold=Max(1.0e-12_real64,&
          1.0e-10_real64*Max(1.0_real64,Abs(ProjectedValues(k))))
        Do i=1,Dimension
          Denominator=Diagonal(i)-ProjectedValues(k)
          If(Abs(Denominator).lt.DenominatorThreshold) Then
            If(Real(Denominator,Kind=real64).lt.0.0_real64) Then
              Denominator=Cmplx(-DenominatorThreshold,0.0_real64,Kind=real64)
            Else
              Denominator=Cmplx(DenominatorThreshold,0.0_real64,Kind=real64)
            EndIf
          EndIf
          WorkVector(i)=WorkVector(i)/Denominator
        EndDo
        Call MQC_Iterative_Orthogonalize(WorkVector,Basis(:,1:SubspaceSize))
        If(CorrectionCount.gt.0) Call MQC_Iterative_Orthogonalize(WorkVector,&
          Corrections(:,1:CorrectionCount))
        NormValue=Sqrt(Real(Dot_Product(WorkVector,WorkVector),Kind=real64))
        If(NormValue.gt.1.0e-12_real64) Then
          CorrectionCount=CorrectionCount+1
          Corrections(:,CorrectionCount)=WorkVector/NormValue
        EndIf
      EndDo
      AllConverged=All(RootResiduals.le.ConvergenceTolerance)
      If(AllConverged.or.CorrectionCount.eq.0) Exit

      If(SubspaceSize+CorrectionCount.gt.MaxSubspace) Then
        Basis=(0.0_real64,0.0_real64)
        Basis(:,1:NumberOfRoots)=RitzVectors
        SubspaceSize=NumberOfRoots
        Call MQC_Iterative_Orthonormalize_Block(Basis,SubspaceSize)
        NewCount=0
        Do k=1,Min(CorrectionCount,MaxSubspace-SubspaceSize)
          WorkVector=Corrections(:,k)
          Call MQC_Iterative_Orthogonalize(WorkVector,Basis(:,1:SubspaceSize))
          If(NewCount.gt.0) Call MQC_Iterative_Orthogonalize(WorkVector,&
            Basis(:,SubspaceSize+1:SubspaceSize+NewCount))
          NormValue=Sqrt(Real(Dot_Product(WorkVector,WorkVector),Kind=real64))
          If(NormValue.gt.1.0e-12_real64) Then
            NewCount=NewCount+1
            Basis(:,SubspaceSize+NewCount)=WorkVector/NormValue
          EndIf
        EndDo
        SubspaceSize=SubspaceSize+NewCount
        Call Operator%Apply(Basis(:,1:SubspaceSize),OperatorBasis(:,1:SubspaceSize))
      Else
        Basis(:,SubspaceSize+1:SubspaceSize+CorrectionCount)=&
          Corrections(:,1:CorrectionCount)
        Call Operator%Apply(Basis(:,SubspaceSize+1:SubspaceSize+CorrectionCount),&
          OperatorBasis(:,SubspaceSize+1:SubspaceSize+CorrectionCount))
        SubspaceSize=SubspaceSize+CorrectionCount
      EndIf
    EndDo

    Call Eigenvalues%Init(NumberOfRoots)
    Call Eigenvectors%Init(Dimension,NumberOfRoots,&
      Scalar=(0.0_real64,0.0_real64))
    Do k=1,NumberOfRoots
      Call Eigenvalues%Put(ProjectedValues(k),k)
      Do i=1,Dimension
        Call Eigenvectors%Put(RitzVectors(i,k),i,k)
      EndDo
    EndDo
    If(Present(ResidualNorms)) Then
      Call ResidualNorms%Init(NumberOfRoots)
      Do k=1,NumberOfRoots
        Call ResidualNorms%Put(RootResiduals(k),k)
      EndDo
    EndIf
    If(Present(Converged)) Converged=AllConverged
    If(Present(Iterations)) Iterations=Min(IterationCount,MaxIterations)
  End Procedure MQC_Davidson

!> Construct the Hermitian Ritz matrix in the current Davidson subspace.
!> \author L. M. Thompson
!> \date 2026
  Subroutine MQC_Iterative_Projected_Matrix(Basis,OperatorBasis,ProjectedMatrix)
    Complex(Kind=real64),Intent(In)::Basis(:,:),OperatorBasis(:,:)
    Type(MQC_Matrix),Intent(Out)::ProjectedMatrix
    Complex(Kind=real64)::Element
    Integer(Kind=int64)::i,j,Dimension

    Dimension=Size(Basis,2)
    Call ProjectedMatrix%Init(Dimension,Dimension,&
      Scalar=(0.0_real64,0.0_real64))
    Do j=1,Dimension
      Do i=j,Dimension
        Element=Dot_Product(Basis(:,i),OperatorBasis(:,j))
        If(i.eq.j) Element=Cmplx(Real(Element,Kind=real64),0.0_real64,Kind=real64)
        Call ProjectedMatrix%Put(Element,i,j)
        Call ProjectedMatrix%Put(Conjg(Element),j,i)
      EndDo
    EndDo
  End Subroutine MQC_Iterative_Projected_Matrix

!> Orthonormalize a block and discard numerically dependent columns.
!> \author L. M. Thompson
!> \date 2026
  Subroutine MQC_Iterative_Orthonormalize_Block(Vectors,NumberOfVectors)
    Complex(Kind=real64),Intent(InOut)::Vectors(:,:)
    Integer(Kind=int64),Intent(InOut)::NumberOfVectors
    Complex(Kind=real64),Allocatable::WorkVector(:)
    Real(Kind=real64)::NormValue
    Integer(Kind=int64)::Source,Target

    Allocate(WorkVector(Size(Vectors,1)))
    Target=0
    Do Source=1,NumberOfVectors
      WorkVector=Vectors(:,Source)
      If(Target.gt.0) Call MQC_Iterative_Orthogonalize(WorkVector,Vectors(:,1:Target))
      NormValue=Sqrt(Real(Dot_Product(WorkVector,WorkVector),Kind=real64))
      If(NormValue.gt.1.0e-12_real64) Then
        Target=Target+1
        Vectors(:,Target)=WorkVector/NormValue
      EndIf
    EndDo
    NumberOfVectors=Target
  End Subroutine MQC_Iterative_Orthonormalize_Block

!> Reorthogonalize one vector against a basis using two MGS passes.
!> \author L. M. Thompson
!> \date 2026
  Subroutine MQC_Iterative_Orthogonalize(Vector,Basis)
    Complex(Kind=real64),Intent(InOut)::Vector(:)
    Complex(Kind=real64),Intent(In)::Basis(:,:)
    Complex(Kind=real64)::Overlap
    Integer(Kind=int64)::Pass,i

    Do Pass=1,2
      Do i=1,Size(Basis,2)
        Overlap=Dot_Product(Basis(:,i),Vector)
        Vector=Vector-Basis(:,i)*Overlap
      EndDo
    EndDo
  End Subroutine MQC_Iterative_Orthogonalize

End Submodule MQC_Algebra_Iterative
