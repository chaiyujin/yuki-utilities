#pragma once

#include "core.h"
#include <iostream>
#include <algorithm>
#include <memory>
#include "Eigen/Eigen"
#define _USE_MATH_DEFINES
#include <math.h>

namespace yuki
{
	namespace math
	{
		template <typename T>
		struct range_
		{
		public:
			T start, end;
			bool is_valid() const { return start <= end; }
		};

		typedef range_<float>  range_float;
		typedef range_<double> range_double;
		typedef range_<int>	   range_int;

		/* -----------------------------------------------vec234-------------------------------------------------- */
		template <typename T> struct vec2_;
		template <typename T> struct vec3_;
		template <typename T> struct vec4_;

		template <typename T>
		struct vec2_
		{
		public:
			T x, y;
			vec2_(T x = 0, T y = 0) : x(x), y(y) {}
			vec2_(const vec2_<T> &b) : x(b.x), y(b.y) {}
			vec2_(const vec3_<T> &b);
			vec2_(const vec4_<T> &b);
			bool isnan() const { return std::isnan(x) || std::isnan(y); }
			T * operator&() { return &x; }
			const T * operator&() const { return &x; }
			bool operator<(const vec2_ &b) const { return (x < b.x || (x == b.x && y < b.y)); }
			bool operator==(const vec2_ &b) const { return x == b.x && y == b.y; }
		};

		template <typename T>
		struct vec3_
		{
		public:
			T x, y, z;
			vec3_(T x = 0, T y = 0, T z = 0) : x(x), y(y), z(z) {}
			vec3_(const vec2_<T> &b) : x(b.x), y(b.y), z(0) {}
			vec3_(const vec3_<T> &b) : x(b.x), y(b.y), z(b.z) {}
			vec3_(const vec4_<T> &b);
			bool isnan() const { return std::isnan(x) || std::isnan(y) || std::isnan(z); }
			T &r() { return x; } const T &r() const { return x; }
			T &g() { return y; } const T &g() const { return y; }
			T &b() { return z; } const T &b() const { return z; }
			T * operator&() { return &x; }
			const T * operator&() const { return &x; }
			bool operator<(const vec3_ &b) const { return (x < b.x || (x == b.x && y < b.y) || (x == b.x && y == b.y && z < b.z)); }
			bool operator==(const vec3_ &b) const { return x == b.x && y == b.y && z == b.z; }
			T &operator[](int idx) { if (idx == 0) return x; else if (idx == 1) return y; else return z;}
			const T &operator[](int idx) const { if (idx == 0) return x; else if (idx == 1) return y; else return z;}
			vec3_ operator+(const vec3_ &b) const
			{
				return vec3_(x + b.x, y + b.y, z + b.z);
			}
			template <typename U>
			vec3_ operator/(const U v) const
			{
				return vec3_(x / (T)v, y / (T)v, z / (T)v);
			}
			template <typename U>
			vec3_ operator*(const U v) const
			{
				return vec3_(x * (T)v, y * (T)v, z * (T)v);
			}
			vec3_ normalize() const
			{
				T length = sqrt(x * x + y * y + z * z);
				return vec3_(x / length, y / length, z / length);
			}
		};

		template <typename T>
		struct vec4_
		{
		public:
			T x, y, z, w;
			vec4_(T x = 0, T y = 0, T z = 0, T w = 0) : x(x), y(y), z(z), w(w) {}
			vec4_(const vec2_<T> &b) : x(b.x), y(b.y), z(0), w(0) {}
			vec4_(const vec3_<T> &b) : x(b.x), y(b.y), z(b.z), w(0) {}
			vec4_(const vec4_<T> &b) : x(b.x), y(b.y), z(b.z), w(b.w) {}
			bool isnan() const { return std::isnan(x) || std::isnan(y) || std::isnan(z) || std::isnan(w); }
			T &r() { return x; } const T &r() const { return x; }
			T &g() { return y; } const T &g() const { return y; }
			T &b() { return z; } const T &b() const { return z; }
			T &a() { return w; } const T &a() const { return w; }
			T * operator&() { return &x; }
			const T * operator&() const { return &x; }
			bool operator<(const vec4_ &b) const { return (x < b.x || (x == b.x && y < b.y) || (x == b.x && y == b.y && z < b.z) || (x == b.x && y == b.y && z == b.z && w < b.w)); }
			bool operator==(const vec4_ &b) const { return x == b.x && y == b.y && z == b.z && w == b.w; }
		};

		template <class T> inline vec2_<T>::vec2_(const vec3_<T> &b) : x(b.x), y(b.y) {}
		template <class T> inline vec2_<T>::vec2_(const vec4_<T> &b) : x(b.x), y(b.y) {}
		template <class T> inline vec3_<T>::vec3_(const vec4_<T> &b) : x(b.x), y(b.y), z(b.z) {}

		template <class T> inline std::ostream & operator<<(std::ostream & out, const vec2_<T> &v) { out << v.x << " " << v.y; return out; }
		template <class T> inline std::ostream & operator<<(std::ostream & out, const vec3_<T> &v) { out << v.x << " " << v.y << " " << v.z; return out; }
		template <class T> inline std::ostream & operator<<(std::ostream & out, const vec4_<T> &v) { out << v.x << " " << v.y << " " << v.z << " " << v.w; return out; }
		template <class T> inline std::istream & operator >> (std::istream & in, const vec2_<T> &v) { in >> v.x >> v.y; return in; }
		template <class T> inline std::istream & operator >> (std::istream & in, const vec3_<T> &v) { in >> v.x >> v.y >> v.z; return in; }
		template <class T> inline std::istream & operator >> (std::istream & in, const vec4_<T> &v) { in >> v.x >> v.y >> v.z >> v.w; return in; }

		template <typename T>
		T cross(const vec2_<T> &a, const vec2_<T> &b)
		{
			return a.x * b.y - a.y * b.x;
		}
		template <typename T>
		vec3_<T> cross(const vec3_<T> &a, const vec3_<T> &b)
		{
			return vec3_<T>(
				a.y * b.z - a.z * b.y,
				a.z * b.x - a.x * b.z,
				a.x * b.y - a.y * b.x
				);
		}
		template <typename T>
		T dot(const vec2_<T> &a, const vec2_<T> &b)
		{
			return a.x * b.x + a.y * b.y;
		}
		template <typename T>
		T dot(const vec3_<T> &a, const vec3_<T> &b)
		{
			return a.x * b.x + a.y * b.y + a.z * b.z;
		}

		typedef vec2_<int> int2;
		typedef vec3_<int> int3;
		typedef vec2_<float> float2;
		typedef vec3_<float> float3;
		typedef vec2_<double> double2;
		typedef vec3_<double> double3;
		typedef vec3_<unsigned char> rgb24;
		typedef vec4_<unsigned char> rgba32;
		typedef vec4_<float> rgba_flt;

		template <typename T> Eigen::Matrix<T, 2, 1> to_eigen(const vec2_<T> &v) { return Eigen::Matrix<T, 2, 1>(v.x, v.y); }
		template <typename T> Eigen::Matrix<T, 3, 1> to_eigen(const vec3_<T> &v) { return Eigen::Matrix<T, 3, 1>(v.x, v.y, v.z); }
		template <typename T> Eigen::Matrix<T, 4, 1> to_eigen(const vec4_<T> &v) { return Eigen::Matrix<T, 4, 1>(v.x, v.y, v.z, v.w); }

		/* ------------------------------------------------------------------------------------------------- */

		/* --------------------------------------- easy function ------------------------------------------- */
		template <typename T>
		T clamp_value(const T &x, const T &a, const T &b)
		{
			if (x < a) return a;
			if (b < x) return b;
			return x;
		}
		
		template <typename T>
		T accumlate_mul(const std::vector<T> &vec, int s = 0, int e = INT32_MAX)
		{
			T ret(1);
			if (e < 0) e += (int)vec.size();
			e = clamp_value(e, 0, (int)vec.size());
			for (int i = s; i < e; ++i) ret *= vec[i];
			return ret;
		}
		template <typename T>
		T accumlate_add(const std::vector<T> &vec, int s = 0, int e = INT32_MAX)
		{
			T ret(0);
			if (e < 0) e += (int)vec.size();
			e = clamp_value(e, 0, (int)vec.size());
			for (int i = s; i < e; ++i) ret += vec[i];
			return ret;
		}
		/* ------------------------------------------------------------------------------------------------- */

		/* -------------------------------------------Tensor3------------------------------------------------ */

		template <int Dim>
		class Shape
		{
			std::vector<int> v_;
		public:
			Shape() : v_(Dim, 1) {}
			Shape(const std::vector<int> &v) : Shape()
			{
				assert(v.size() == Dim);
				for (size_t i = 0; i < v.size() && i < Dim; ++i)
					v_[i] = v[i];
			}
			Shape(const Shape &b) : v_(b.v_) {}

			size_t dims() const { return v_.size(); }
			size_t size() const { return accumlate_mul(v_); }
			int &operator[](size_t dim) { return v_[dim]; }
			const int &operator[](size_t dim) const { return v_[dim]; }
			std::vector<int> &data() { return v_; }
			const std::vector<int> &data() const { return v_; }
			bool operator==(const Shape<Dim> &b) const
			{
				for (size_t i = 0; i < v_.size(); ++i)
					if (v_[i] != b.v_[i]) return false;
				return true;
			}
			bool operator!=(const Shape<Dim> &b) const { return !((*this) == b); }
		};

		template <int Dim> inline std::ostream &operator<<(std::ostream &out, const Shape<Dim> &s)
		{
			out << "[";
			for (size_t i = 0; i < s.dims(); ++i)
				out << " " << s[i];
			out << " ]";
			return out;
		}

		class Tensor3
		{
			std::shared_ptr<double>			p_real_;
			double *						p_head_;
			Shape<3>						shape_;
			bool							is_sub_;
		public:
			Tensor3() : p_real_(NULL), p_head_(NULL), shape_(), is_sub_(false) {}
			Tensor3(const Tensor3 &b) : p_real_(b.p_real_), p_head_(b.p_head_), shape_(b.shape_), is_sub_(b.is_sub_) {}
			Tensor3(const std::vector<int> &s) : Tensor3() { resize(s); }
			Tensor3(const Shape<3> &s) : Tensor3() { resize(s.data()); }
			Tensor3(const Tensor3 &b, int start, int length=-1) : Tensor3()
			{
				assert(start < b.shape_[0]);
				is_sub_ = true;  // it's a part of a father Tensor3
				shape_ = b.shape_;
				shape_[0] = (length < 0) ? shape_[0] - start : length;
				p_head_ = p_head_ + start * accumlate_mul(b.shape_.data(), 1);
			}

			void resize(const std::vector<int> &s)
			{
				if (is_sub_) { std::cerr << "[Tensor3]: do not resize sub Tensor3.\n"; return; }
				Shape<3> newshape_(s);
				if (shape_ != newshape_)
				{
					shape_ = newshape_;
					p_real_.reset(new double[shape_.size()]);
					p_head_ = p_real_.get();
				}
			}

			const Shape<3> &shape() const { return shape_; }
			bool is_null() const { return p_head_ == NULL; }
			void set_zero() { if (p_head_) memset(p_head_, 0, sizeof(double) * shape_.size()); }
			double *data() { return p_head_; }
			const double*data() const { return p_head_; }

			template <typename T>
			void unfold_data(T *data, int unfold_mode, bool is_data_colmajor = true)
			{
				assert(unfold_mode >= 0 && unfold_mode < shape_.size());

				std::vector<int> indices(shape_.dims());
				std::vector<int> div(shape_.dims());  // idx -> i, j, k...
				std::vector<int> mul(shape_.dims());  // idx <- i, j, k...
				int D = indices.size();
	
				if (is_data_colmajor)
				{
					div.back() = 1;
					for (int i = D - 2; i >= 0; --i)
					{
						div[i] = div[i + 1] * shape_[i + 1];
					}
					for (int i = 0, d = unfold_mode, v = 1; i < D; ++i)
					{
						mul[d] = v;
						v *= shape_[d];
						d = (d + 1) % D;
					}
				}
				else
				{
					std::cerr << "not support so far.\n";
					exit(1);
					return;
				}

				int size_ = shape_.size();
				for (int i = 0; i < size_; ++i)
				{
					int ii = i;
					int j = 0;
					for (int d = 0; d < D; ++d)
					{
						indices[d] = ii / div[d];
						ii %= div[d];
						j += indices[d] * mul[d];
					}
					p_head_[i] = data[j];
				}
			}

			template <int Mode>
			void mul_vec(const double*vec, Tensor3 &result) const;
			void mul_vec(const double*vec1, const double*vec2, Tensor3 &result) const;
			template <int Mode>
			void mul_mat(const double*mat, int rows, int cols, Tensor3 &result) const;
			void mul(double value, Tensor3 &result) const;

			template <int Mode>
			void unfold(Tensor3 &result);

			void mesh_eigen(Eigen::MatrixXd &mat)
			{
				assert(shape_[1] == shape_[2] && shape_[1] == 1);
				mat.resize(3, shape_[0] / 3);
				memcpy(mat.data(), p_head_, sizeof(double) * mat.size());
			}
		};


		template <>
		inline void Tensor3::mul_vec<1>(const double*vec, Tensor3 &result) const
		{
			if (!(result.shape_[0] == shape_[0] && result.shape_[1] == 1 && result.shape_[2] == shape_[2]))
				result.resize({ shape_[0], 1, shape_[2] });
			result.set_zero();
#pragma omp parallel for
			for (int i = 0; i < shape_[0]; ++i)
			{
				int c = i * shape_[2];
				int c0 = i * shape_[1] * shape_[2];
				for (int j = 0; j < shape_[1]; ++j)
				{
					int c1 = j * shape_[2];
					for (int k = 0; k < shape_[2]; ++k)
					{
						result.p_head_[c + k] += p_head_[c0 + c1 + k] * vec[j];
					}
				}
			}
		}

		template <>
		inline void Tensor3::mul_vec<2>(const double*vec, Tensor3 &result) const
		{
			if (!(result.shape_[0] == shape_[0] && result.shape_[1] == shape_[1] && result.shape_[2] == 1))
				result.resize({ shape_[0], shape_[1], 1 });
			result.set_zero();
#pragma omp parallel for
			for (int i = 0; i < shape_[0]; ++i)
			{
				int c = i * shape_[1];
				int c0 = i * shape_[1] * shape_[2];
				for (int j = 0; j < shape_[1]; ++j)
				{
					int c1 = j * shape_[2];
					for (int k = 0; k < shape_[2]; ++k)
					{
						result.p_head_[c + j] += p_head_[c0 + c1 + k] * vec[k];
					}
				}
			}
		}

		inline void Tensor3::mul_vec(const double*vec1, const double*vec2, Tensor3 &result) const
		{
			if (!(result.shape_[0] == shape_[0] && result.shape_[1] == 1 && result.shape_[2] == 1))
				result.resize({ shape_[0], 1, 1 });
			result.set_zero();
#pragma omp parallel for
			for (int i = 0; i < shape_[0]; ++i)
			{
				int c0 = i * shape_[1] * shape_[2];
				for (int j = 0; j < shape_[1]; ++j)
				{
					int c1 = j * shape_[2];
					for (int k = 0; k < shape_[2]; ++k)
					{
						result.p_head_[i] += p_head_[c0 + c1 + k] * vec1[j] * vec2[k];
					}
				}
			}
		}

		template <>
		inline void Tensor3::mul_mat<1>(const double*mat, int rows, int cols, Tensor3 &result) const
		{
			assert(rows == shape_[1]);
			if (!(result.shape_[0] == shape_[0] && result.shape_[1] == cols && result.shape_[2] == shape_[2]))
				result.resize({ shape_[0], cols, shape_[2] });
			result.set_zero();
#pragma omp parallel for
			for (int i = 0; i < shape_[0]; ++i)
			{
				int c0 = i * shape_[1] * shape_[2];
				int r0 = i * cols * shape_[2];
				for (int j = 0; j < shape_[1]; ++j)
				{
					int c1 = j * shape_[2];
					int idx = j * cols;
					for (int m = 0; m < cols; ++m)
					{
						int r1 = m * shape_[2];
						for (int k = 0; k < shape_[2]; ++k)
						{
							result.p_head_[r0 + r1 + k] += p_head_[c0 + c1 + k] * mat[idx + m];
						}
					}
				}
			}
		}

		template <>
		inline void Tensor3::mul_mat<2>(const double*mat, int rows, int cols, Tensor3 &result) const
		{
			assert(rows == shape_[2]);
			if (!(result.shape_[0] == shape_[0] && result.shape_[1] == shape_[1] && result.shape_[2] == cols))
				result.resize({ shape_[0], shape_[1], cols });
			result.set_zero();
#pragma omp parallel for
			for (int i = 0; i < shape_[0]; ++i)
			{
				int c0 = i * shape_[1] * shape_[2];
				int r0 = i * shape_[1] * cols;
				for (int j = 0; j < shape_[1]; ++j)
				{
					int c1 = j * shape_[2];
					int r1 = j * cols;
					for (int k = 0; k < shape_[2]; ++k)
					{
						int idx = k * cols;
						for (int m = 0; m < cols; ++m)
						{
							result.p_head_[r0 + r1 + m] += p_head_[c0 + c1 + k] * mat[idx + m];
						}
					}
				}
			}
		}

		template <>
		inline void Tensor3::unfold<1>(Tensor3 &result)
		{
			if (!(result.shape_[0] == shape_[1] && result.shape_[1] == shape_[2] * shape_[0] && result.shape_[2] == 1))
				result.resize({ shape_[1], shape_[2] * shape_[0], 1 });
			result.set_zero();
#pragma omp parallel for
			for (int i = 0; i < shape_[0]; ++i)
			{
				int c0 = i * shape_[1] * shape_[2];
				for (int j = 0; j < shape_[1]; ++j)
				{
					int c1 = j * shape_[2];
					int r0 = j * shape_[2] * shape_[0];
					for (int k = 0; k < shape_[2]; ++k)
					{
						result.p_head_[r0 + k * shape_[0] + i] = p_head_[c0 + c1 + k];
					}
				}
			}
		}

		template <>
		inline void Tensor3::unfold<2>(Tensor3 &result)
		{
			if (!(result.shape_[0] == shape_[2] && result.shape_[1] == shape_[0] * shape_[1] && result.shape_[2] == 1))
				result.resize({ shape_[2], shape_[0] * shape_[1], 1 });
			result.set_zero();
			int shape01 = shape_[0] * shape_[1];
#pragma omp parallel for
			for (int i = 0; i < shape_[0]; ++i)
			{
				int c0 = i * shape_[1] * shape_[2];
				int r1 = i * shape_[1];
				for (int j = 0; j < shape_[1]; ++j)
				{
					int c1 = j * shape_[2];
					for (int k = 0; k < shape_[2]; ++k)
					{
						result.p_head_[k * shape01 + r1 + j] = p_head_[c0 + c1 + k];
					}
				}
			}
		}

		inline void Tensor3::mul(double scale, Tensor3 &result) const
		{
			if (!(result.shape_ == shape_))
				result.resize(shape_.data());
#pragma omp parallel for
			for (int i = 0; i < shape_[0]; ++i)
			{
				int c0 = i * shape_[1] * shape_[2];
				for (int j = 0; j < shape_[1]; ++j)
				{
					int c1 = j * shape_[2];
					for (int k = 0; k < shape_[2]; ++k)
					{
						result.p_head_[c0 + c1 + k] = p_head_[c0 + c1 + k] * scale;
					}
				}
			}
		}

		/* ------------------------------------------------------------------------------------------------- */

	}

	namespace math
	{
		/* used frequently in dsp */

		template <typename T>
		T nextpow2(T n)
		{
			if (n <= 1) { return 1; }
			double p = std::log2(n);
			if (p > (int)p) return std::pow(2, (int)p + 1);
			return n;
		}
		/*
			* The actual computation :
			*      - in    : the input vector which defines the toeplitz matrix
			*      - size  : size of in (ie number of elements)
			*      - order : size of the system to solve. order must be < size -1
			*      - acoeff: solution (ie ar coefficients). Size must be at last order+1
			*      - err   : *prediction* error (scalar)
			*      - kcoeff: reflexion coefficients. Size must be at last equal to equal to order.
			*      - tmp   : cache, must have at least order elements, if NULL, will be allocated and free in this function
			*
			* this function assume all arrays are allocated with the right size, and that
			* the parameters make sense. No checking is done, must be done before calling
			* this function: in particular, in[0] must be non zero.
			*
			* Returns 0 on success, -1 if a compuation error happened (overflow, underflow
			* for error calculation)
		*/
		int levinson(const double *in, int order, double *acoeff, double *err, double *kcoeff, double *tmp=NULL);

		/*
			* Return the roots of a polynomial with coefficients given in p.
			* The values in the rank-1 array `p` are coefficients of a polynomial.
			* If the length of `p` is n then the polynomial is described by::
			*     p[0] * x**(n - 1) + p[1] * x**(n-2) + ... + p[n-2]*x + p[n - 1]
		*/
		std::vector<std::complex<double>> roots(const double *p, int n);

		/*
			Evaluate a polynomial at specific values.
    		If `p` is of length N, this function returns the value:
        	``p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]`
		*/
		template <typename T, typename U>
		T polyval(const U *p, int n, T x, bool reverse=false)
		{
			T ret(0);
			for (int i = 0; i < n; ++i)
				ret = ret * x + ((reverse)? p[n - 1 - i] : p[i]);
			return ret;
		}
	}

}
