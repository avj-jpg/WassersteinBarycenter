#include "kron_laplace.hpp"

namespace mfem
{

KroneckerLaplacian::KroneckerLaplacian(
   ParFiniteElementSpace &Vs, FiniteElementSpace &Vt, const double mass_coeff_)
   : Ls(&Vs),
     Ms(&Vs),
     Lt(&Vt),
     Mt(&Vt),
     nVs(Vs.GetTrueVSize()),
     nVt(Vt.GetTrueVSize()),
     mass_coeff(mass_coeff_)
{
   height = width = nVs*nVt;
   const int dim = Vs.GetMesh()->Dimension();

   const bool pa = dim >= 2 && UsesTensorBasis(Vs);

   Ls.AddDomainIntegrator(new DiffusionIntegrator);
   if (pa) { Ls.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   Ls.Assemble();
   Ls.FormSystemMatrix(empty, Ls_op);

   Ms.AddDomainIntegrator(new MassIntegrator);
   if (pa) { Ms.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   Ms.Assemble();
   Ms.FormSystemMatrix(empty, Ms_op);

   Lt.AddDomainIntegrator(new DiffusionIntegrator);
   Lt.Assemble();
   Lt.FormSystemMatrix(empty, Lt_op);

   Mt.AddDomainIntegrator(new MassIntegrator);
   Mt.Assemble();
   Mt.FormSystemMatrix(empty, Mt_op);
}

void KroneckerLaplacian::Mult(const Vector &x, Vector &y) const
{
   MFEM_ASSERT(x.Size() == nVs*nVt, "");
   MFEM_ASSERT(y.Size() == nVs*nVt, "");

   kron_mult.Mult(*Mt_op, *Ms_op, x, y);
   y *= mass_coeff;
   kron_mult.AddMult(*Mt_op, *Ls_op, x, y);
   kron_mult.AddMult(*Lt_op, *Ms_op, x, y);
}

void KroneckerLaplacian::AssembleDiagonal(Vector &diag) const
{
   Vector diag_Ms(nVs), diag_Ls(nVs), diag_Mt(nVt), diag_Lt(nVt);

   Ms_op->AssembleDiagonal(diag_Ms);
   Ls_op->AssembleDiagonal(diag_Ls);
   Mt_op.As<SparseMatrix>()->GetDiag(diag_Mt);
   Lt_op.As<SparseMatrix>()->GetDiag(diag_Lt);

   diag.SetSize(nVs*nVt);

   for (int it = 0; it < nVt; ++it)
   {
      const double dMt = diag_Mt[it];
      const double dLt = diag_Lt[it];
      for (int is = 0; is < nVs; ++is)
      {
         const double dMs = diag_Ms[is];
         const double dLs = diag_Ls[is];
         diag[is + it*nVs] = mass_coeff*dMt*dMs + dMt*dLs + dLt*dMs;
      }
   }
}

} // namespace mfem
