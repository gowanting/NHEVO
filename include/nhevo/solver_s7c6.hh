#pragma once

#include <nhevo/solver.hh>

namespace nhevo {

class SolverS7C6 : public Solver
{
public:
  SolverS7C6();

  ~SolverS7C6();

  std::vector<double> solve(const std::vector<double>& xs,
                            const std::vector<double>& delta_t,
                            const double param) override;

private:
  Eigen::Matrix<double, 9, 27> getCoeffs(const double x, const double delta_t);

  std::vector<double> getDetCoeffs(const Eigen::Matrix<double, 9, 27>& C);

  Eigen::Matrix<double, -1, -1> coeffA(const double x,
                                       const double delta_t) override;
};

}
