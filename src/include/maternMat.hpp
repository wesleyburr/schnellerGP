#ifndef __RHODLR_MATERN_MAT__
#define __RHODLR_MATERN_MAT__

#include "HODLR_Matrix.hpp"
#include "HODLR_Tree.hpp"

// Taking Matern Kernel for p = 2:
// K(r) = σ^2 * [1 + sqrt(5) * r / ρ + 5/3 * (r / ρ)^2] * exp[-sqrt(5) * r / ρ]
class Matern_Kernel : public HODLR_Matrix
{
  private:
    Mat x;
    double sigma_squared, rho;
    int D;

  public:

    // Constructor:
    Matern_Kernel(Mat tX, int N, double sigma, double rho) : HODLR_Matrix(N), x(tX)
    {
        this -> sigma_squared = sigma * sigma;
        this -> rho           = rho;
        this -> D             = tX.cols();
    };

    dtype getMatrixEntry(int i, int j)
    {
        double temp = 0;
        for(int d = 0; d < D; d++){
            double temp2 = x(i, d) - x(j, d);
            temp = temp + temp2 * temp2;
        }
        double R_by_rho = sqrt(temp) / rho;
        double nugget = (i == j)?1e-8:0.0;  // if i==j, return 1e-8, else, return 0.0
        return sigma_squared * (1.0 + sqrt(5.0) * R_by_rho + (5.0/3.0) * (R_by_rho * R_by_rho)) * exp(-sqrt(5.0) * R_by_rho) + nugget;
    }

    // Destructor:
    ~Matern_Kernel() {};
};

#endif
