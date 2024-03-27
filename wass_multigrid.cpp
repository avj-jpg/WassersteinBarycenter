#include "wass_multigrid.hpp"
#include "wass_laplace.hpp"

namespace mfem
{

KroneckerMultigrid::KroneckerMultigrid(ParFiniteElementSpaceHierarchy &Hs,
                                       ParFiniteElementSpaceHierarchy &Ht,
                                       const double mass_coeff,
                                       const int bdr_offset)
{
   for (int i = 0; i < Hs.GetNumLevels() - 1; ++i)
   {
      const Operator &Ps = *Hs.GetProlongationAtLevel(i);
      const Operator &Pt = *Ht.GetProlongationAtLevel(i);
      prolongations.emplace_back(new KroneckerProlongation(Pt, Ps));
   }

   for (int i = 0; i < Hs.GetNumLevels(); ++i)
   {
      const double coarse = i == 0;
      FormLevel(Hs.GetFESpaceAtLevel(i), Ht.GetFESpaceAtLevel(i), coarse, mass_coeff, bdr_offset);
   }
}

void KroneckerMultigrid::FormLevel(
   ParFiniteElementSpace &Vs,
   ParFiniteElementSpace &Vt,
   const bool coarse,
   const double mass_coeff,
   const int bdr_offset)
{
   auto *A = new KroneckerLaplacian(Vs, Vt, mass_coeff, bdr_offset);
   Vector diag;
   A->AssembleDiagonal(diag);

   Solver *S;
   if (coarse)
   {
      coarse_jacobi.reset(new OperatorJacobiSmoother(diag, empty));
      CGSolver *cg = new CGSolver(Vs.GetComm());
      cg->SetMaxIter(10);
      cg->SetRelTol(1e-4);
      cg->SetOperator(*A);
      cg->SetPreconditioner(*coarse_jacobi);
      S = cg;
   }
   else
   {
      S = new OperatorChebyshevSmoother(*A, diag, empty, 2, Vs.GetComm());
   }
   AddLevel(A, S, true, true);
}

const Operator *KroneckerMultigrid::GetProlongationAtLevel(int level) const
{
   return prolongations[level].get();
}

} // namespace mfem
