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


std::vector<std::complex<double>> roots(double *p, int n)
{
    Eigen::MatrixXd mat(n - 1, n - 1);
    mat.setZero();
    mat.bottomLeftCorner(n - 2, n - 2) = Eigen::MatrixXd::Identity(n - 2, n - 2);
    for (int i = 1; i < n; ++i)
        mat(n - i - 1, n - 2) = -p[i] / p[0];
    Eigen::ComplexEigenSolver<Eigen::MatrixXd> solver;
    solver.compute(mat);
    std::cout << solver.eigenvalues() << std::endl;
    return std::vector<std::complex<double>>();
}

}}