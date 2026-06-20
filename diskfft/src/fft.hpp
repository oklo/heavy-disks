// Minimal self-contained radix-2 Cooley-Tukey FFT (complex double).
// Sizes must be powers of two.  Used for the Kalnajs convolution; swap in FFTW
// later for speed if needed.
#ifndef DISKFFT_FFT_HPP
#define DISKFFT_FFT_HPP

#include <complex>
#include <vector>
#include <cmath>

namespace diskfft {

using cdouble = std::complex<double>;

// In-place 1D FFT. sign=-1 forward, sign=+1 inverse (no 1/N normalisation).
inline void fft1d(std::vector<cdouble> &a, int sign) {
  const int n = static_cast<int>(a.size());
  // bit-reversal permutation
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (int len = 2; len <= n; len <<= 1) {
    double ang = sign * 2.0 * M_PI / len;
    cdouble wlen(std::cos(ang), std::sin(ang));
    for (int i = 0; i < n; i += len) {
      cdouble w(1.0, 0.0);
      for (int k = 0; k < len / 2; ++k) {
        cdouble u = a[i + k];
        cdouble v = a[i + k + len / 2] * w;
        a[i + k] = u + v;
        a[i + k + len / 2] = u - v;
        w *= wlen;
      }
    }
  }
}

// 2D FFT of a row-major (n1 x n2) array, in place. sign as in fft1d.
inline void fft2d(std::vector<cdouble> &a, int n1, int n2, int sign) {
  std::vector<cdouble> row(n2);
  for (int i = 0; i < n1; ++i) {                 // transform each row (n2)
    for (int j = 0; j < n2; ++j) row[j] = a[i * n2 + j];
    fft1d(row, sign);
    for (int j = 0; j < n2; ++j) a[i * n2 + j] = row[j];
  }
  std::vector<cdouble> col(n1);
  for (int j = 0; j < n2; ++j) {                 // transform each column (n1)
    for (int i = 0; i < n1; ++i) col[i] = a[i * n2 + j];
    fft1d(col, sign);
    for (int i = 0; i < n1; ++i) a[i * n2 + j] = col[i];
  }
}

// inverse 2D FFT with 1/(n1*n2) normalisation
inline void ifft2d(std::vector<cdouble> &a, int n1, int n2) {
  fft2d(a, n1, n2, +1);
  double inv = 1.0 / (static_cast<double>(n1) * n2);
  for (auto &z : a) z *= inv;
}

}  // namespace diskfft
#endif
