#include "bary_rhs.hpp"
//#include "bary_params.hpp"

namespace mfem
{

KroneckerLinearForm::KroneckerLinearForm(
   ParFiniteElementSpace &Vs_, ParFiniteElementSpace &Vt_,
   QuadratureSpace &Qs_, QuadratureSpace &Qt_)
   : Vector(Vs_.GetTrueVSize() * Vt_.GetTrueVSize()),
     Vs(Vs_),
     Vt(Vt_),
     dim_s(Vs.GetMesh()->SpaceDimension()), // background space dimension
     Qs(Qs_),
     Qt(Qt_),
     Ls_grad(&Vs),
     Ls_interp(&Vs),
     Lt_grad(&Vt),
     Lt_interp(&Vt),
     qf_s(Qs, Qt, dim_s),
     qf_t(Qs, Qt, 1),
     qf_scalar(Qs, Qt, 1),
     qs(Qs),
     qs_vec(Qs, dim_s),
     qt(Qt),
     qs_coeff(qs),
     qt_coeff(qt),
     qs_vec_coeff(qs_vec),
     qt_vec_coeff(qt),
     z1(Vs.GetTrueVSize()*Qt.GetSize()),
     z2(Vs.GetTrueVSize()*Qt.GetSize()),
     z3(Vt.GetTrueVSize()),
     z4(Vs.GetTrueVSize())
{
   auto set_int_rule = [&](ParLinearForm &lf, const IntegrationRule &ir)
   {
      for (auto *i : *lf.GetDLFI()) { i->SetIntRule(&ir); }
   };

   Ls_grad.AddDomainIntegrator(new DomainLFGradIntegrator(qs_vec_coeff));
   Ls_interp.AddDomainIntegrator(new DomainLFIntegrator(qs_coeff));

   Lt_grad.AddDomainIntegrator(new DomainLFGradIntegrator(qt_vec_coeff));
   Lt_interp.AddDomainIntegrator(new DomainLFIntegrator(qt_coeff));

   set_int_rule(Ls_grad, Qs.GetIntRule(0));
   set_int_rule(Ls_interp, Qs.GetIntRule(0));
   set_int_rule(Lt_grad, Qt.GetIntRule(0));
   set_int_rule(Lt_interp, Qt.GetIntRule(0));
}

void KroneckerLinearForm::Assemble()
{
   const int nQs = Qs.GetSize();
   const int nQt = Qt.GetSize();

   const int nVs = Vs.GetTrueVSize();
   const int nVt = Vt.GetTrueVSize();

   // Compute space integrals for each time quadrature point
   for (int it = 0; it < nQt; ++it)
   {
      qs_vec.MakeRef(qf_s, it*nQs*dim_s, nQs*dim_s);
      qs.MakeRef(qf_t, it*nQs, nQs);

      Vector z1_slice(z1, it*nVs, nVs);
      Vector z2_slice(z2, it*nVs, nVs);

      Ls_grad.Assemble();
      Ls_grad.ParallelAssemble(z1_slice);
      Ls_interp.Assemble();
      Ls_interp.ParallelAssemble(z2_slice);

      // Add source term
      qs.MakeRef(qf_scalar, it*nQs, nQs);
      Ls_interp.Assemble();
      Ls_interp.ParallelAssemble(z4);
      z1_slice += z4;
   }

   for (int is = 0; is < nVs; ++is)
   {
      for (int it = 0; it < nQt; ++it)
      {
         qt[it] = z1[is + it*nVs];
      }
      Lt_interp.Assemble();
      Lt_interp.ParallelAssemble(z3);

      for (int it = 0; it < nVt; ++it)
      {
         (*this)[is + it*nVs] = z3[it];
      }

      for (int it = 0; it < nQt; ++it)
      {
         qt[it] = z2[is + it*nVs];
      }
      Lt_grad.Assemble();
      Lt_grad.ParallelAssemble(z3);

      for (int it = 0; it < nVt; ++it)
      {
         (*this)[is + it*nVs] += z3[it];
      }
   }
}

void KroneckerLinearForm::Update(
   KroneckerQuadratureFunction &u_qf, KroneckerQuadratureFunction &s_qf,
   KroneckerQuadratureFunction &phi_qf, int k, int N, double sigma_phi)
{
   const int nq = qf_t.Size();
   const int dim = dim_s + 1;
         
   for (int i = 0; i < nq; ++i)
   {
      int ks0 = (k+N-1)%N, ks1 = k, kr0 = (k+N-1)%N, kr1 = (k+1)%N;
      qf_scalar[i] = sigma_phi * (s_qf[i*N+ks0] - s_qf[i*N+ks1])
          + phi_qf(i*N+kr0) + phi_qf(i*N+kr1);
      for (int d = 0; d < dim_s; ++d)
      {
         // m components
         qf_s[d + i*dim_s] = - sigma_phi * u_qf[i*dim*N + dim*k +d];
      }
      // rho component
      qf_t[i] = -sigma_phi * u_qf[i*dim*N+dim*k+dim_s];
   }

   // std::cout << "Scalar: \n";
   // qf_scalar.Print(std::cout, 1);
   // std::cout << "\n\n\n";
   // std::cout << "Space: \n";
   // qf_s.Print(std::cout, 1);
   // std::cout << "\n\n\n";
   // std::cout << "Time: \n";
   // qf_t.Print(std::cout, 1);
   // std::cout << "\n\n\n";

   Assemble();
}

} // namespace mfem
