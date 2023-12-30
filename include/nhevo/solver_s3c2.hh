#pragma once

#include <nhevo/solver.hh>

namespace nhevo {

class SolverS3C2 : public Solver
{
public:
  SolverS3C2();

  ~SolverS3C2();

  std::vector<double> solve(const std::vector<double>& xs,
                            const std::vector<double>& delta_t,
                            const double param) override;

private:
  Eigen::Matrix<double, 9, 11> getCoeffs(const double x, const double delta_t);

  std::vector<double> getDetCoeffs(const Eigen::Matrix<double, 9, 11>& C);

  Eigen::Matrix<double, -1, -1> coeffA(const double x,
                                       const double delta_t) override;
};

}
