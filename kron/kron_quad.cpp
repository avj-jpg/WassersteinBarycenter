#include "kron_quad.hpp"

namespace mfem

{

KroneckerQuadratureFunction::KroneckerQuadratureFunction(
   QuadratureSpace &Qs_, QuadratureSpace &Qt_, const int vdim_)
   : Vector(Qs_.GetSize()*Qt_.GetSize()*vdim_),
     Qs(Qs_),
     Qt(Qt_),
     vdim(vdim_)
{ }

void KroneckerQuadratureFunction::ProjectComponent(Coefficient &coeff, int c)
{
   Vector tvec, svec, xvec, values;

   for (int iel_t = 0; iel_t < Qt.GetNE(); ++iel_t)
   {
      const IntegrationRule &ir_t = Qt.GetElementIntRule(iel_t);
      ElementTransformation &T_t = *Qt.GetTransformation(iel_t);
      const int nq_t = ir_t.Size();
      for (int iq_t = 0; iq_t < nq_t; ++iq_t)
      {
         const IntegrationPoint &ip_t = ir_t[iq_t];
         T_t.Transform(ip_t, tvec);

         const double t = tvec[0];

         // Fill in the space QuadratureFunction qs
         for (int iel_s = 0; iel_s < Qs.GetNE(); ++iel_s)
         {
            const IntegrationRule &ir_s = Qs.GetElementIntRule(iel_s);
            ElementTransformation &T_s = *Qs.GetTransformation(iel_s);
            const int nq_s = ir_s.Size();
            for (int iq_s = 0; iq_s < nq_s; ++iq_s)
            {
               const IntegrationPoint &ip_s = ir_s[iq_s];

               const double val = coeff.Eval(T_s, ip_s, t);
               const int s_idx = iq_s + iel_s*nq_s;
               const int t_idx = iq_t + iel_t*nq_t;
               const int idx = s_idx + Qs.GetNE()*nq_s*t_idx;
               (*this)[c + idx*vdim] = val;
            }
         }
      }
   }
}

void KroneckerQuadratureFunction::Project(Coefficient &coeff)
{
   MFEM_VERIFY(vdim == 1, "Wrong vdim");
   ProjectComponent(coeff, 0);
}

void KroneckerQuadratureFunction::ProjectGradient(
   Vector &u, ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt)
{
   const ElementDofOrdering ordering = Vs.GetElementDofOrdering();
   const Operator *Rs = Vs.GetElementRestriction(ordering);
   const Operator *Rt = Vt.GetElementRestriction(ordering);

   const int nQs = Qs.GetSize();
   const int nQt = Qt.GetSize();
   const int nVs = Vs.GetTrueVSize();
   const int nVt = Vt.GetTrueVSize();
   // FIXME: for surface mesh dim_s is the space dimension
   const int dim_s = Vs.GetMesh()->SpaceDimension();
   const int dim = dim_s + 1;
   const int s_vdim = Vs.GetVDim();

   MFEM_VERIFY(vdim == s_vdim*dim, "Bad vdim");

   Vector tv_s, tv_t;
   ParGridFunction gf_s(&Vs);
   ParGridFunction gf_t(&Vt);
   Vector ev_s(Rs->Height());
   Vector ev_t(Rt->Height());
   Vector grad_s(nQs * s_vdim * dim_s);
   Vector val_s(nQs * s_vdim);
   Vector q_t(nQt);

   Vector z(nVt * nQs * vdim); // intermediate vector

   const QuadratureInterpolator *qi_s = Vs.GetQuadratureInterpolator(Qs);
   qi_s->SetOutputLayout(QVectorLayout::byVDIM);
   qi_s->DisableTensorProducts(false);

   for (int it = 0; it < nVt; ++it)
   {
      tv_s.MakeRef(u, it*nVs, nVs);
      gf_s.SetFromTrueDofs(tv_s);
      Rs->Mult(gf_s, ev_s);

      qi_s->PhysDerivatives(ev_s, grad_s);
      qi_s->Values(ev_s, val_s);

      grad_s.HostReadWrite();
      val_s.HostReadWrite();

      for (int is = 0; is < nQs; ++is)
      {
         for (int vd = 0; vd < s_vdim*dim_s; ++vd)
         {
            z[it + is*nVt + vd*nQs*nVt] = grad_s[vd + is*s_vdim*dim_s];
         }
         for (int vd = 0; vd < s_vdim; ++vd)
         {
            z[it + is*nVt + (vd + s_vdim*dim_s)*nQs*nVt] = val_s[vd + is*s_vdim];
         }
      }
   }

   const QuadratureInterpolator *qi_t = Vt.GetQuadratureInterpolator(Qt);
   qi_t->SetOutputLayout(QVectorLayout::byVDIM);
   qi_t->DisableTensorProducts(false);

   for (int is = 0; is < nQs; ++is)
   {
      // Interpolate spatial derivatives
      for (int vd = 0; vd < s_vdim*dim_s; ++vd)
      {
         tv_t.MakeRef(z, is*nVt + vd*nQs*nVt, nVt);
         gf_t.SetFromTrueDofs(tv_t);
         Rt->Mult(gf_t, ev_t);
         qi_t->Values(ev_t, q_t);

         q_t.HostReadWrite();

         for (int it = 0; it < nQt; ++it)
         {
            (*this)[vd + is*vdim + it*nQs*vdim] = q_t[it];
         }
      }

      for (int vd = 0; vd < s_vdim; ++vd)
      {
         const int offset = vd + s_vdim*dim_s;
         tv_t.MakeRef(z, is*nVt + offset*nQs*nVt, nVt);
         gf_t.SetFromTrueDofs(tv_t);
         Rt->Mult(gf_t, ev_t);
         qi_t->PhysDerivatives(ev_t, q_t);

         q_t.HostReadWrite();

         for (int it = 0; it < nQt; ++it)
         {
            (*this)[offset + is*vdim + it*nQs*vdim] = q_t[it];
         }
      }
   }
}

void KroneckerQuadratureFunction::ProjectDivergence(
   Vector &u, ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt)
{
   const ElementDofOrdering ordering = Vs.GetElementDofOrdering();
   const Operator *Rs = Vs.GetElementRestriction(ordering);
   const Operator *Rt = Vt.GetElementRestriction(ordering);

   const int nQs = Qs.GetSize();
   const int nQt = Qt.GetSize();
   const int nVs = Vs.GetTrueVSize();
   const int nVt = Vt.GetTrueVSize();

   MFEM_VERIFY(vdim == 1, "Bad vdim");

   Vector tv_s, tv_t;
   ParGridFunction gf_s(&Vs);
   ParGridFunction gf_t(&Vt);
   Vector ev_s(Rs->Height());
   Vector ev_t(Rt->Height());
   Vector div_s(nQs);
   Vector q_t(nQt);

   Vector z(nVt * nQs * vdim); // intermediate vector

   const QuadratureInterpolator *qi_s = Vs.GetQuadratureInterpolator(Qs);
   qi_s->SetOutputLayout(QVectorLayout::byVDIM);
   qi_s->DisableTensorProducts(false);

   for (int it = 0; it < nVt; ++it)
   {
      tv_s.MakeRef(u, it*nVs, nVs);
      gf_s.SetFromTrueDofs(tv_s);
      Rs->Mult(gf_s, ev_s);
      qi_s->PhysDivergence(ev_s, div_s);

      div_s.HostReadWrite();

      for (int is = 0; is < nQs; ++is)
      {
         z[it + is*nVt] = div_s[is];
      }
   }

   const QuadratureInterpolator *qi_t = Vt.GetQuadratureInterpolator(Qt);
   qi_t->SetOutputLayout(QVectorLayout::byVDIM);
   qi_t->DisableTensorProducts(false);

   for (int is = 0; is < nQs; ++is)
   {
      tv_t.MakeRef(z, is*nVt, nVt);
      gf_t.SetFromTrueDofs(tv_t);
      Rt->Mult(gf_t, ev_t);
      qi_t->Values(ev_t, q_t);

      q_t.HostReadWrite();

      for (int it = 0; it < nQt; ++it)
      {
         (*this)[is + it*nQs] = q_t[it];
      }
   }
}

void KroneckerQuadratureFunction::ProjectValue(
   Vector &u, ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt)
{
   const ElementDofOrdering ordering = Vs.GetElementDofOrdering();
   const Operator *Rs = Vs.GetElementRestriction(ordering);
   const Operator *Rt = Vt.GetElementRestriction(ordering);

   const int nQs = Qs.GetSize();
   const int nQt = Qt.GetSize();
   const int nVs = Vs.GetTrueVSize();
   const int nVt = Vt.GetTrueVSize();

   //MFEM_VERIFY(vdim == Vs.GetVDim() || vdim == Vs.GetFE(0)->GetRangeDim(), "Bad vdim");
   //std::cout << vdim << ": "<< Vs.GetVDim() << ": " << Vs.GetFE(0)->GetRangeDim() << std::endl;

   Vector tv_s, tv_t;
   ParGridFunction gf_s(&Vs);
   ParGridFunction gf_t(&Vt);
   Vector ev_s(Rs->Height());
   Vector ev_t(Rt->Height());
   Vector val_s(nQs*vdim);
   Vector q_t(nQt);

   Vector z(nVt * nQs * vdim); // intermediate vector

   const QuadratureInterpolator *qi_s = Vs.GetQuadratureInterpolator(Qs);
   qi_s->SetOutputLayout(QVectorLayout::byVDIM);
   qi_s->DisableTensorProducts(false);

   for (int it = 0; it < nVt; ++it)
   {
      tv_s.MakeRef(u, it*nVs, nVs);
      gf_s.SetFromTrueDofs(tv_s);
      Rs->Mult(gf_s, ev_s);

      qi_s->Values(ev_s, val_s);

      val_s.HostReadWrite();

      for (int is = 0; is < nQs; ++is)
      {
         for (int vd = 0; vd < vdim; ++vd)
         {
            z[it + is*nVt + vd*nQs*nVt] = val_s[vd + is*vdim];
         }
      }
   }

   const QuadratureInterpolator *qi_t = Vt.GetQuadratureInterpolator(Qt);
   qi_t->SetOutputLayout(QVectorLayout::byVDIM);
   qi_t->DisableTensorProducts(false);

   for (int is = 0; is < nQs; ++is)
   {
      for (int vd = 0; vd < vdim; ++vd)
      {
         tv_t.MakeRef(z, is*nVt + vd*nQs*nVt, nVt);
         gf_t.SetFromTrueDofs(tv_t);
         Rt->Mult(gf_t, ev_t);
         qi_t->Values(ev_t, q_t);

         q_t.HostReadWrite();

         for (int it = 0; it < nQt; ++it)
         {
            (*this)[vd + is*vdim + it*nQs*vdim] = q_t[it];
         }
      }
   }
}

// TODO
void KroneckerQuadratureFunction::ProjectSpaceGradient(
   Vector &u, ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt)
{
   const ElementDofOrdering ordering = Vs.GetElementDofOrdering();
   const Operator *Rs = Vs.GetElementRestriction(ordering);
   const Operator *Rt = Vt.GetElementRestriction(ordering);

   const int nQs = Qs.GetSize();
   const int nQt = Qt.GetSize();
   const int nVs = Vs.GetTrueVSize();
   const int nVt = Vt.GetTrueVSize();
   const int dim_s = Vs.GetMesh()->Dimension();
   const int s_vdim = Vs.GetVDim();

   MFEM_VERIFY(vdim == s_vdim*dim_s, "Bad vdim");

   Vector tv_s, tv_t;
   ParGridFunction gf_s(&Vs);
   ParGridFunction gf_t(&Vt);
   Vector ev_s(Rs->Height());
   Vector ev_t(Rt->Height());
   Vector grad_s(nQs * s_vdim * dim_s);
   Vector q_t(nQt);

   Vector z(nVt * nQs * vdim); // intermediate vector

   const QuadratureInterpolator *qi_s = Vs.GetQuadratureInterpolator(Qs);
   qi_s->SetOutputLayout(QVectorLayout::byVDIM);
   qi_s->DisableTensorProducts(false);

   for (int it = 0; it < nVt; ++it)
   {
      tv_s.MakeRef(u, it*nVs, nVs);
      gf_s.SetFromTrueDofs(tv_s);
      Rs->Mult(gf_s, ev_s);

      qi_s->PhysDerivatives(ev_s, grad_s);

      grad_s.HostReadWrite();

      for (int is = 0; is < nQs; ++is)
      {
         for (int vd = 0; vd < s_vdim*dim_s; ++vd)
         {
            z[it + is*nVt + vd*nQs*nVt] = grad_s[vd + is*s_vdim*dim_s];
         }
      }
   }

   const QuadratureInterpolator *qi_t = Vt.GetQuadratureInterpolator(Qt);
   qi_t->SetOutputLayout(QVectorLayout::byVDIM);
   qi_t->DisableTensorProducts(false);

   for (int is = 0; is < nQs; ++is)
   {
      // Interpolate spatial derivatives
      for (int vd = 0; vd < s_vdim*dim_s; ++vd)
      {
         tv_t.MakeRef(z, is*nVt + vd*nQs*nVt, nVt);
         gf_t.SetFromTrueDofs(tv_t);
         Rt->Mult(gf_t, ev_t);
         qi_t->Values(ev_t, q_t);

         q_t.HostReadWrite();

         for (int it = 0; it < nQt; ++it)
         {
            (*this)[vd + is*vdim + it*nQs*vdim] = q_t[it];
         }
      }
   }
}

double KroneckerQuadratureFunction::L1Norm() const
{
   MFEM_VERIFY(vdim == 1, "Only implemented for vdim == 1");

   const Vector &ws = Qs.GetWeights();
   const Vector &wt = Qt.GetWeights();

   const int nQs = Qs.GetSize();
   const int nQt = Qt.GetSize();

   double integ = 0.0;
   for (int it = 0; it < nQt; ++it)
   {
      const double wt_i = wt[it];
      for (int is = 0; is < nQs; ++is)
      {
         const double val = (*this)[is + it*nQs];
         integ += wt_i*ws[is]*std::abs(val);
      }
   }

   if (auto *pmesh = dynamic_cast<const ParMesh*>(Qs.GetMesh()))
   {
      MPI_Comm comm = pmesh->GetComm();
      MPI_Allreduce(MPI_IN_PLACE, &integ, 1, MPI_DOUBLE, MPI_SUM, comm);
   }

   return integ;
}

void GetTimeSlice(const Vector &X, Vector &X_slice, const double t,
                  const int vd, const int vdim,
                  const ParFiniteElementSpace &Vs,
                  const ParFiniteElementSpace &Vt,
                  bool nearest_nbr)
{
   const Mesh &mesh_t = *Vt.GetMesh();
   int found_elem = -1;
   IntegrationPoint ip;
   // Find time element
   for (int iel_t = 0; iel_t < mesh_t.GetNE(); ++iel_t)
   {
      Array<int> v;
      mesh_t.GetElementVertices(iel_t, v);
      const double a = *mesh_t.GetVertex(v[0]);
      const double b = *mesh_t.GetVertex(v[1]);
      const double t0 = std::min(a, b);
      const double t1 = std::max(a, b);

      if (t >= t0 && t <= t1)
      {
         found_elem = iel_t;
         ip.x = (t - t0)/(t1 - t0);
         break;
      }
   }
   MFEM_VERIFY(found_elem >= 0, "Could not find time element");

   Array<int> time_dofs;
   Vt.GetElementVDofs(found_elem, time_dofs);
   const int nt = time_dofs.Size();

   Vector shape(nt);
   Vt.GetFE(found_elem)->CalcShape(ip, shape);
   const IntegrationRule &time_ir = Vt.GetFE(found_elem)->GetNodes();

   int idx = -1;
   double min_diff = 2.0;
   for (int it = 0; it < nt; ++it)
   {
      const double diff = std::abs(time_ir[it].x - ip.x);
      if (diff < min_diff)
      {
         idx = it;
         min_diff = diff;
      }
   }

   const int nVs = Vs.GetTrueVSize();
   X_slice.SetSize(nVs);
   for (int is = 0; is < nVs; ++is)
   {
      if (nearest_nbr)
      {
         X_slice[is] = X[vd + vdim*(is + time_dofs[idx]*nVs)];
      }
      else
      {
         X_slice[is] = 0.0;
         for (int it = 0; it < nt; ++it)
         {
            X_slice[is] += X[vd + vdim*(is + time_dofs[it]*nVs)]*shape[it];
         }
      }
   }
}

KroneckerFaceInterpolator::KroneckerFaceInterpolator(ParFiniteElementSpace &Vs_,
                                                     QuadratureSpace &Qs_)
   : Vs(Vs_), Qs(Qs_), gf_s(&Vs)
{ }

void KroneckerFaceInterpolator::Project(const Vector &u_tdof,
                                        QuadratureFunction &qf, int offset) const
{
   const int nV = Vs.GetTrueVSize();
   const Vector u_tdof_slice(const_cast<Vector&>(u_tdof), offset*nV, nV);
   const ElementDofOrdering ordering = Vs.GetElementDofOrdering();
   const Operator *Rs = Vs.GetElementRestriction(ordering);
   const QuadratureInterpolator *qi_s = Vs.GetQuadratureInterpolator(Qs);

   gf_s.SetFromTrueDofs(u_tdof_slice);
   ev_s.SetSize(Rs->Height());
   Rs->Mult(gf_s, ev_s);

   qi_s->Values(ev_s, qf);
}

} // namespace mfem
