#ifndef KRON_QUAD_HPP
#define KRON_QUAD_HPP

#include "mfem.hpp"

namespace mfem
{

class KroneckerQuadratureFunction : public Vector
{
   QuadratureSpace &Qs, &Qt;
   const int vdim;
public:
   KroneckerQuadratureFunction(QuadratureSpace &Qs_, QuadratureSpace &Qt_, const int vdim_);
   void ProjectComponent(Coefficient &coeff, int c);
   void Project(Coefficient &coeff);
   void ProjectGradient(Vector &u, ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt);
   void ProjectSpaceGradient(Vector &u, ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt);
   void ProjectDivergence(Vector &u, ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt);
   void ProjectValue(Vector &u, ParFiniteElementSpace &Vs, ParFiniteElementSpace &Vt);
   using Vector::operator=;
   int GetVDim() const { return vdim; }
   double L1Norm() const;
};

void GetTimeSlice(const Vector &X, Vector &X_slice, const double t,
                  const int vd, const int vdim,
                  const ParFiniteElementSpace &Vs,
                  const ParFiniteElementSpace &Vt,
                  bool nearest_nbr = false);

class KroneckerFaceInterpolator
{
   ParFiniteElementSpace &Vs;
   QuadratureSpace &Qs;
   mutable ParGridFunction gf_s;
   mutable Vector ev_s;
public:
   KroneckerFaceInterpolator(ParFiniteElementSpace &Vs_, QuadratureSpace &Qs_);
   void Project(const Vector &u_tdof, QuadratureFunction &qf, int offset) const;
};

} // namespace mfem

#endif
