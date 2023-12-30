#pragma once

#include <nhevo/solver.hh>

namespace nhevo {

class SolverS5C4 : public Solver
{
public:
  SolverS5C4();

  ~SolverS5C4();

  std::vector<double> solve(const std::vector<double>& xs,
                            const std::vector<double>& delta_t,
                            const double param) override;

private:
  Eigen::Matrix<double, 9, 19> getCoeffs(const double x, const double delta_t);

  std::vector<double> getDetCoeffs(const Eigen::Matrix<double, 9, 19>& C);

  Eigen::Matrix<double, -1, -1> coeffA(const double x,
                                       const double delta_t) override;
};

}
