// 183702374
#include <complex>
#include <iostream>
#include <vector>

namespace math {
const double kPI = 3.14159265358979323846;
}

template <class Field>
class FFTConverter;

template <typename Type>
class IField {
 public:
  using ObjectT = Type;

  virtual ObjectT PrimitiveRootOfUnity(long long) = 0;

  virtual ObjectT InvertedPrimitiveRootOfUnity(long long) = 0;

  friend class FFTConverter<IField<Type>>;
};

class ComplexField : public IField<std::complex<double>> {
 public:
  ObjectT PrimitiveRootOfUnity(long long len) override {
    double div = len;
    double angle = 2 * math::kPI / div;
    return ObjectT(std::cos(angle), std::sin(angle));
  }

  ObjectT InvertedPrimitiveRootOfUnity(long long len) override {
    double div = (-1) * len;
    double angle = 2 * math::kPI / div;
    return ObjectT(std::cos(angle), std::sin(angle));
  }
};

template <typename Type>
class Polynomial {
 public:
  std::vector<Type> coeffs;

  Polynomial() = default;

  Polynomial(const std::vector<Type>& input) : coeffs(input) {}
};

template <class Field = ComplexField>
class FFTConverter {
 public:
  using ObjectT = typename Field::ObjectT;

  template <typename Type>
  void FFTConvert(Polynomial<Type>& polynomial, bool is_inverted) {
    long long polynom_degree = (long long)polynomial.coeffs.size();
    for (long long i = 1, j = 0; i < polynom_degree; ++i) {
      long long bit = polynom_degree >> 1;
      for (; j >= bit; bit >>= 1) {
        j -= bit;
      }
      j += bit;
      if (i < j) {
        std::swap(polynomial.coeffs[i], polynomial.coeffs[j]);
      }
    }
    Field field;
    for (long long len = 2; len <= polynom_degree; len <<= 1) {
      int coeff = (is_inverted ? -1 : 1);
      ObjectT primitive_root = field.PrimitiveRootOfUnity(coeff * len);
      for (long long i = 0; i < polynom_degree; i += len) {
        ObjectT root_of_unity(1);
        for (long long j = 0; j < len / 2; ++j) {
          ObjectT even = polynomial.coeffs[i + j];
          ObjectT odd = polynomial.coeffs[i + j + len / 2] * root_of_unity;
          polynomial.coeffs[i + j] = even + odd;
          polynomial.coeffs[i + j + len / 2] = even - odd;
          root_of_unity *= primitive_root;
        }
      }
    }
  }

  template <typename Type>
  void Convert(Polynomial<Type>& polynomial) {
    FFTConvert(polynomial, false);
  }

  template <typename Type>
  void Invert(Polynomial<Type>& polynomial) {
    long long polynom_degree = (long long)polynomial.coeffs.size();
    FFTConvert(polynomial, true);
    for (size_t i = 0; i < polynomial.coeffs.size(); ++i) {
      polynomial.coeffs[i] /= polynom_degree;
    }
  }
};

template <typename Type>
std::istream& operator>>(std::istream& input, Polynomial<Type>& polynom) {
  size_t degree;
  input >> degree;
  polynom.coeffs.resize(degree + 1);
  for (size_t i = 0; i <= degree; ++i) {
    input >> polynom.coeffs[degree - i];
  }
  return input;
}

template <typename Type>
std::ostream& operator<<(std::ostream& out, const Polynomial<Type>& polynom) {
  // output: degree & coeffs
  size_t degree = polynom.coeffs.size() - 1;
  out << degree << ' ';
  for (size_t i = 0; i <= degree; ++i) {
    out << polynom.coeffs[degree - i] << ' ';
  }
  return out;
}

template <typename Type>
Polynomial<Type> operator*(const Polynomial<Type>& first,
                           const Polynomial<Type>& second) {
  const double kRounding = 0.5;
  std::vector<std::complex<double>> fvec(first.coeffs.begin(),
                                         first.coeffs.end());
  std::vector<std::complex<double>> svec(second.coeffs.begin(),
                                         second.coeffs.end());
  Polynomial<std::complex<double>> fft_first(fvec);
  Polynomial<std::complex<double>> fft_second(svec);
  Polynomial<Type> res;
  FFTConverter<ComplexField> fft;
  size_t degree = 1;
  while (degree < std::max(fvec.size(), svec.size())) {
    degree <<= 1;
  }
  degree <<= 1;
  fft_first.coeffs.resize(degree);
  fft_second.coeffs.resize(degree);
  fft.Convert(fft_first);
  fft.Convert(fft_second);
  for (size_t i = 0; i < degree; ++i) {
    fft_first.coeffs[i] *= fft_second.coeffs[i];
  }
  fft.Invert(fft_first);
  res.coeffs.resize(degree);
  for (size_t i = 0; i < degree; ++i) {
    res.coeffs[i] = Type(std::abs(fft_first.coeffs[i].real()) + kRounding);
    if (fft_first.coeffs[i].real() < 0) {
      res.coeffs[i] *= -1;
    }
  }
  size_t last = res.coeffs.size();
  while (last > 0 && std::abs(res.coeffs[--last]) == 0) {
    res.coeffs.pop_back();
  }
  return res;
}

int main() {
  Polynomial<int> first;
  Polynomial<int> second;
  std::cin >> first >> second;
  std::cout << first * second;
}
 
