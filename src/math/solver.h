#pragma once

#include "mathutils.h"
#include <vector>
#include <complex>

namespace yuki
{
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