#include "kron_multigrid.hpp"
#include "kron_laplace.hpp"

namespace mfem
{

class KroneckerCoarseSolver : public Solver
{
   const int nVs, nVt;
   ParBilinearForm a_s;
   ConstantCoefficient coeff_L, coeff_M;
   Array<int> empty;
   HypreParMatrix A_s;
   HypreBoomerAMG amg;

public:
   KroneckerCoarseSolver(ParFiniteElementSpace &Vs,
                         ParFiniteElementSpace &Vt,
                         const double mass_coeff)
      : nVs(Vs.GetTrueVSize()),
        nVt(Vt.GetTrueVSize()),
        a_s(&Vs)
   {
      height = width = nVs*nVt;

      auto avg_diag = [&](BilinearFormIntegrator *integ)
      {
         Vector diag(nVt);
         BilinearForm a_t(&Vt);
         a_t.AddDomainIntegrator(integ);
         a_t.Assemble();
         a_t.Finalize();
         a_t.AssembleDiagonal(diag);
         return diag.Sum()/nVt;
      };

      const double avg_L = avg_diag(new DiffusionIntegrator);
      const double avg_M = avg_diag(new MassIntegrator);

      coeff_L.constant = avg_M;
      coeff_M.constant = avg_L + mass_coeff*avg_M;
      a_s.AddDomainIntegrator(new DiffusionIntegrator(coeff_L));
      a_s.AddDomainIntegrator(new MassIntegrator(coeff_M));
      a_s.Assemble();
      a_s.FormSystemMatrix(empty, A_s);
      amg.SetOperator(A_s);
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      x.HostRead();
      y.HostWrite();

      Vector slice_x, slice_y;
      for (int i = 0; i < nVt; ++i)
      {
         slice_x.MakeRef(const_cast<Vector&>(x), i*nVs, nVs);
         slice_y.MakeRef(y, i*nVs, nVs);
         amg.Mult(slice_x, slice_y);
      }
   }

   void SetOperator(const Operator &op) override { }
};

KroneckerMultigrid::KroneckerMultigrid(ParFiniteElementSpaceHierarchy &Hs,
                                       ParFiniteElementSpaceHierarchy &Ht,
                                       const double mass_coeff)
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
      FormLevel(Hs.GetFESpaceAtLevel(i), Ht.GetFESpaceAtLevel(i), coarse, mass_coeff);
   }
}

void KroneckerMultigrid::FormLevel(ParFiniteElementSpace &Vs,
                                   ParFiniteElementSpace &Vt,
                                   const bool coarse,
                                   const double mass_coeff)
{
   auto *A = new KroneckerLaplacian(Vs, Vt, mass_coeff);
   Vector diag;
   A->AssembleDiagonal(diag);

   Solver *S;
   if (coarse)
   {
      S = new KroneckerCoarseSolver(Vs, Vt, mass_coeff);
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
