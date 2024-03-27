#ifndef KRON_MULT_HPP
#define KRON_MULT_HPP

#include "mfem.hpp"
#include <memory>

namespace mfem
{

class KronMult
{
   mutable Vector z, slice_x, slice_y, slice_z;
public:
   template <bool ADD, bool TRANSPOSE>
   void Mult_(const Operator &A, const Operator &B, const Vector &x, Vector &y) const;

   void Mult(const Operator &A, const Operator &B, const Vector &x, Vector &y) const;

   void AddMult(const Operator &A, const Operator &B, const Vector &x, Vector &y) const;

   void MultTranspose(const Operator &A, const Operator &B, const Vector &x, Vector &y) const;

   void AddMultTranspose(const Operator &A, const Operator &B, const Vector &x, Vector &y) const;
};

} // namespace mfem

#endif
