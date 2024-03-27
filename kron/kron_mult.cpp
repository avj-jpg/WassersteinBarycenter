#include "kron_mult.hpp"

namespace mfem
{

template <bool ADD, bool TRANSPOSE>
void KronMult::Mult_(const Operator &A, const Operator &B, const Vector &x, Vector &y) const
{
   const int ms = TRANSPOSE ? B.Width() : B.Height();
   const int mt = TRANSPOSE ? A.Width() : A.Height();

   const int ns = TRANSPOSE ? B.Height() : B.Width();
   const int nt = TRANSPOSE ? A.Height() : A.Width();

   // In the below, A really means op(A), where op is either transpose or not,
   // depending on TRANSPOSE, and likewsie for B.

   // A is of size (mt) x (nt)
   // B is of size (ms) x (ns)

   // (A (x) B) x = B X A^t.

   x.HostRead();
   x.Read();
   y.HostWrite();

   // Set Z = B X.
   z.UseDevice(true);
   z.SetSize(ms*nt);
   z.Write();
   z = 0.0;

   for (int i = 0; i < nt; ++i)
   {
      slice_x.MakeRef(const_cast<Vector&>(x), i*ns, ns);
      slice_z.MakeRef(z, i*ms, ms);
      if (TRANSPOSE) { B.MultTranspose(slice_x, slice_z); }
      else { B.Mult(slice_x, slice_z); }
   }
   z.GetMemory().Sync(slice_z.GetMemory());
   z.HostReadWrite();
   // Set Y = Z A^t (or add if requested)
   Vector Zi(nt), Yi(mt);
   for (int i = 0; i < ms; ++i)
   {
      Zi.HostWrite();
      Yi.Write();
      for (int j = 0; j < nt; ++j) { Zi(j) = z(i + j*ms); }
      if (TRANSPOSE) { A.MultTranspose(Zi, Yi); }
      else { A.Mult(Zi, Yi); }
      Yi.HostReadWrite();
      if (ADD)
      {
         for (int j = 0; j < mt; ++j) { y(i + j*ms) += Yi(j); }
      }
      else
      {
         for (int j = 0; j < mt; ++j) { y(i + j*ms) = Yi(j); }
      }
   }
}

void KronMult::Mult(const Operator &A, const Operator &B, const Vector &x, Vector &y) const
{ Mult_<false, false>(A, B, x, y); }

void KronMult::AddMult(const Operator &A, const Operator &B, const Vector &x, Vector &y) const
{ Mult_<true, false>(A, B, x, y); }

void KronMult::MultTranspose(const Operator &A, const Operator &B, const Vector &x, Vector &y) const
{ Mult_<false, true>(A, B, x, y); }

void KronMult::AddMultTranspose(const Operator &A, const Operator &B, const Vector &x, Vector &y) const
{ Mult_<true, true>(A, B, x, y); }

} // namespace mfem
