#ifndef KRON_LAPLACE_HPP
#define KRON_LAPLACE_HPP

#include "mfem.hpp"
#include "kron_mult.hpp"

namespace mfem
{

class KroneckerLaplacian : public Operator
{
   ParBilinearForm Ls, Ms;
   BilinearForm Lt, Mt;
   const int nVs, nVt;
   const double mass_coeff;

   OperatorHandle Ls_op, Ms_op, Lt_op, Mt_op;

   Array<int> empty;

   KronMult kron_mult;
public:
   KroneckerLaplacian(ParFiniteElementSpace &Vs, FiniteElementSpace &Vt,
                      const double mass_coeff_ = 0.0);

   void Mult(const Vector &x, Vector &y) const override;

   void AssembleDiagonal(Vector &diag) const override;
};

}

#endif
