#ifndef KRON_MULTIGRID_HPP
#define KRON_MULTIGRID_HPP

#include "mfem.hpp"
#include "kron_mult.hpp"
#include <memory>

namespace mfem
{

class KroneckerProlongation : public Operator
{
   const Operator &A, &B;
   KronMult kron_mult;
public:
   KroneckerProlongation(const Operator &A_, const Operator &B_) : A(A_), B(B_) { }
   void Mult(const Vector &x, Vector &y) const { kron_mult.Mult(A, B, x, y); }
   void MultTranspose(const Vector &x, Vector &y) const { kron_mult.MultTranspose(A, B, x, y); }
};

class KroneckerMultigrid : public MultigridBase
{
   Array<int> empty; // no essential DOFs for now...
   std::vector<std::unique_ptr<KroneckerProlongation>> prolongations;
   std::unique_ptr<OperatorJacobiSmoother> coarse_jacobi;
   void FormLevel(ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt, const bool coarse, const double mass_coeff);
public:
   KroneckerMultigrid(ParFiniteElementSpaceHierarchy &Hs,
                      ParFiniteElementSpaceHierarchy &Ht,
                      const double mass_coeff = 0.0);
   const Operator *GetProlongationAtLevel(int level) const override;
};

class MultigridHierarchy
{
   H1_FECollection coarse_fec;
   ParFiniteElementSpace coarse_space;
   std::vector<std::unique_ptr<H1_FECollection>> fe_collections;
   ParFiniteElementSpaceHierarchy space_hierarchy;
   const int order;
public:
   MultigridHierarchy(ParMesh &coarse_mesh, int h_refs, int p_refs)
   : coarse_fec(1, coarse_mesh.Dimension()),
     coarse_space(&coarse_mesh, &coarse_fec),
     space_hierarchy(&coarse_mesh, &coarse_space, false, false),
     order(pow(2, p_refs))
   {
      const int dim = coarse_mesh.Dimension();
      for (int level = 0; level < h_refs; ++level)
      {
         space_hierarchy.AddUniformlyRefinedLevel();
      }
      for (int level = 0; level < p_refs; ++level)
      {
         fe_collections.push_back(std::unique_ptr<H1_FECollection>(
            new H1_FECollection(pow(2, level+1), dim)
         ));
         space_hierarchy.AddOrderRefinedLevel(fe_collections.back().get());
      }
   }

   ParFiniteElementSpace &GetFinestSpace()
   {
      return space_hierarchy.GetFinestFESpace();
   }

   ParMesh &GetFinestMesh()
   {
      return *space_hierarchy.GetFinestFESpace().GetParMesh();
   }

   ParFiniteElementSpaceHierarchy &GetSpaceHierarchy() { return space_hierarchy; }

   int GetOrder() const { return order; }
};

} // namespace mfem

#endif
