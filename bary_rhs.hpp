#ifndef BARY_RHS_HPP
#define BARY_RHS_HPP

#include "mfem.hpp"
#include "kron_quad.hpp"

namespace mfem
{

class KroneckerLinearForm : public Vector
{
   ParFiniteElementSpace &Vs, &Vt;
   const int dim_s;

   QuadratureSpace &Qs, &Qt;

   ParLinearForm Ls_grad, Ls_interp, Lt_grad, Lt_interp;
   KroneckerQuadratureFunction qf_s; // Dot product with the spatial gradient
   KroneckerQuadratureFunction qf_t; // Coefficient of the time derivative
   KroneckerQuadratureFunction qf_scalar; // Coefficient of the source term

   QuadratureFunction qs, qs_vec, qt;
   QuadratureFunctionCoefficient qs_coeff, qt_coeff;
   VectorQuadratureFunctionCoefficient qs_vec_coeff, qt_vec_coeff;

   Vector z1, z2, z3, z4;

   void Assemble();

public:
   KroneckerLinearForm(ParFiniteElementSpace &Vs_, ParFiniteElementSpace &Vt_,
                       QuadratureSpace &Qs_, QuadratureSpace &Qt_);

   void Update(KroneckerQuadratureFunction &u_qf,
               KroneckerQuadratureFunction &q_qf, 
               KroneckerQuadratureFunction &phi_qf,
               int k, int N, double sigma_phi);

   using Vector::operator=;
};

} // namespace mfem

#endif
