#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

namespace Eigen {

namespace internal {

/** \internal Low-level conjugate gradient algorithm
  * \param mat The matrix A
  * \param rhs The right hand side vector b
  * \param x On input and initial solution, on output the computed solution.
  * \param precond A preconditioner being able to efficiently solve for an
  *                approximation of Ax=b (regardless of b)
  * \param iters On input the max number of iteration, on output the number of performed iterations.
  * \param tol_error On input the tolerance error, on output an estimation of the relative error.
  */
template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
EIGEN_DONT_INLINE
void conjugate_projected_gradient(const MatrixType& mat, const Rhs& rhs, Dest& x, const MatrixType& projectionMatrix,
                        const Preconditioner& precond, int& iters,
                        typename Dest::RealScalar& tol_error)
{

  using std::sqrt;
  using std::abs;
  typedef typename Dest::RealScalar RealScalar;
  typedef typename Dest::Scalar Scalar;
  typedef Matrix<Scalar,Dynamic,1> VectorType;

  RealScalar tol = tol_error;
  int maxIters = iters;

  int n = mat.cols();

  VectorType residual = rhs - mat * x; //initial residual

  RealScalar rhsNorm2 = rhs.squaredNorm();
  if(rhsNorm2 == 0)
  {
    x.setZero();
    iters = 0;
    tol_error = 0;
    return;
  }
  RealScalar threshold = tol*tol*rhsNorm2;
  RealScalar residualNorm2 = residual.squaredNorm();
  if (residualNorm2 < threshold)
  {
    iters = 0;
    tol_error = sqrt(residualNorm2 / rhsNorm2);
    return;
  }

  VectorType w = projectionMatrix * precond.solve(residual);      //initial search direction
  VectorType p = w;

//  MatrixXd pInitial(p.rows(), maxIters+1);  // test for reorthogonlization
//  pInitial.col(0) = w;



  VectorType z(n), tmp(n);
  RealScalar absNew = numext::real(w.dot(w));  // the square of the absolute value of r scaled by invM
  RealScalar wNorm2;
  int i = 0;
  while(i < maxIters)
  {
    tmp.noalias() = mat * p;              // the bottleneck of the algorithm

    Scalar alpha = absNew / p.dot(tmp);   // the amount we travel on dir
    x += alpha * p;                       // update solution
    residual -= alpha * tmp;              // update residue

    residualNorm2 = residual.squaredNorm();
    if(residualNorm2 < threshold)
      break;

    w = projectionMatrix * residual;

//   wNorm2 = w.squaredNorm();
//    if (wNorm2 < threshold)
//        break;

    z = precond.solve(w);          // approximately solve for "A z = residual"

    RealScalar absOld = absNew;
    absNew = numext::real(w.dot(z));     // update the absolute value of r
    RealScalar beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
    p = z + beta * p;                             // update search direction

//    pInitial.col(i+1) = z;
//    std::cout << "bla" << std::endl;
//    for(int k = 0; k < i; ++k)
//    {
//        pInitial.col(i+1) -= z.dot(mat*pInitial.col(k))/pInitial.col(k).dot(mat*pInitial.col(k))*pInitial.col(k);
//        std::cout << "bla" << std::endl;
//    }
//    p = pInitial.col(i+1);



    i++;
  }
  tol_error = sqrt(residualNorm2 / rhsNorm2);
  iters = i;
}

}


template< typename _MatrixType, int _UpLo=Lower,
          typename _Preconditioner = DiagonalPreconditioner<typename _MatrixType::Scalar> >
class ConjugateProjectedGradient;

namespace internal {

template< typename _MatrixType, int _UpLo, typename _Preconditioner>
struct traits<ConjugateProjectedGradient<_MatrixType,_UpLo,_Preconditioner> >
{
  typedef _MatrixType MatrixType;
  typedef _Preconditioner Preconditioner;
};

}


template< typename _MatrixType, int _UpLo, typename _Preconditioner>
class ConjugateProjectedGradient : public IterativeSolverBase<ConjugateProjectedGradient<_MatrixType,_UpLo,_Preconditioner> >
{
  typedef IterativeSolverBase<ConjugateProjectedGradient> Base;
  using Base::mp_matrix;
  using Base::m_error;
  using Base::m_iterations;
  using Base::m_info;
  using Base::m_isInitialized;
public:
  typedef _MatrixType MatrixType;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::Index Index;
  typedef typename MatrixType::RealScalar RealScalar;
  typedef _Preconditioner Preconditioner;

  enum {
    UpLo = _UpLo
  };

public:

  /** Default constructor. */
  ConjugateProjectedGradient() : Base() {}

  /** Initialize the solver with matrix \a A for further \c Ax=b solving.
    *
    * This constructor is a shortcut for the default constructor followed
    * by a call to compute().
    *
    * \warning this class stores a reference to the matrix A as well as some
    * precomputed values that depend on it. Therefore, if \a A is changed
    * this class becomes invalid. Call compute() to update it with the new
    * matrix A, or modify a copy of A.
    */
  ConjugateProjectedGradient(const MatrixType& A) : Base(A) {}

  ~ConjugateProjectedGradient() {}



  /** \internal */
  template<typename Rhs,typename Dest>
  void _solveWithGuess(const Rhs& b, Dest& x, const MatrixType& projectionMatrix) const
  {
    m_iterations = Base::maxIterations();
    m_error = Base::m_tolerance;

    for(int j=0; j<b.cols(); ++j)
    {
      m_iterations = Base::maxIterations();
      m_error = Base::m_tolerance;

      typename Dest::ColXpr xj(x,j);
      internal::conjugate_projected_gradient(mp_matrix->template selfadjointView<UpLo>(), b.col(j), xj, projectionMatrix.template selfadjointView<UpLo>(), Base::m_preconditioner, m_iterations, m_error);
    }

    m_isInitialized = true;
    m_info = m_error <= Base::m_tolerance ? Success : NoConvergence;
  }

protected:

};




} // end namespace Eigen
