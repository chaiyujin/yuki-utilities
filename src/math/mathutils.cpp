#include "mathutils.h"

namespace yuki {
namespace math {


int levinson(const double *in, int order, double *acoeff, double *err, double *kcoeff, double *tmp)
{
    bool delete_ = false;
    if (tmp == NULL)
    {
        tmp = new double[order];
        delete_ = true;
    }

    int i, j;
    double acc;
    int ret = 0;

    /* order 0 */
    acoeff[0] = 1.0;
    *err = in[0];

    for (i = 1; i <= order; ++i)
    {
        acc = in[i];
        for (j = 1; j <= i - 1; ++j)
            acc += acoeff[j] * in[i - j];
        kcoeff[i - 1] = -acc / (*err);
        acoeff[i] = kcoeff[i - 1];

        for (j = 0; j < order; ++j)
            tmp[j] = acoeff[j];
        
        for (j = 1; j < i; ++j)
            acoeff[j] += kcoeff[i - 1] * tmp[i - j];
        *err *= (1 - kcoeff[i - 1] * kcoeff[i - 1]);
    }

    if (delete_)
        delete[] tmp;

    return ret;
}


}}