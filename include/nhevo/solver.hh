#pragma once

#include <Eigen/Core>
#include <numeric>
#include <vector>

namespace nhevo {

class Solver
{
public:
  /**
   * @brief Computes 4x4 transformation matrix (including rotation, translation).
   *
   * @param[in] xss A set of horizontal bearing measurements.
   * @param[in] delta_ts A set of time intervals represents the time relative to the start time t0.
   * @param[in] param A numerical parameter that influences the calculation.
   * @param[in] elapse Elapsed time which used in the calculation of motion.
   *
   * @return A 4x4 transformation matrix (including rotation, translation).
   */
  Eigen::Matrix<double, 4, 4> getPose(
    const std::vector<std::vector<double>>& xss,
    const std::vector<std::vector<double>>& delta_ts,
    const double param,
    const double elapse);

  virtual ~Solver() = default;

  /**
   * @brief Sets the value of the 'tau' parameter for the Solver object.
   *
   * @param[in] tau Fixed time interval used to again fix the scale.
   *                For details, kindly refer to the paper :)
   */
  void setTau(const double tau);

  /**
   * @brief Sets the parameters for our bound-constrained minimization process in the Solver object.
   *
   * @param[in] fmin_eps The epsilon value for the minimization procedure,
   *                     controls the precision and stopping criteria.
   * @param[in] fmin_max_iter Maximum number of iterations for the minimization procedure,
   *                          control how long the algorithm runs before termination.
   * @param[in] fmin_init_radius The initial radius used in the minimization process,
   *                             controls the initial search space.
   */
  void setFMinBnd(const double fmin_eps,
                  const long   fmin_max_iter,
                  const double fmin_init_radius);

  /**
   * @brief Sets the parameters for the sturm root-finding algorithm in the Solver object.
   *
   * @param[in] find_root_lb The lower bound of the interval within which the root is to be searched.
   * @param[in] find_root_ub The upper bound of the interval within which the root is to be searched.
   * @param[in] find_root_eps_x Tolerance for determining how close to the actual root the solution needs to be.
   * @param[in] find_root_eps_val Tolerance for how close the polynomial's value must be to zero to consider a solution a root.
   * @param[in] find_root_max_iter The maximum number of iterations the algorithm will perform to find a root.
   */
  void setFindRoot(const double find_root_lb,
                   const double find_root_ub,
                   const double find_root_eps_x,
                   const double find_root_eps_val,
                   const double find_root_max_iter);

protected:
  virtual std::vector<double> solve(const std::vector<double>& xs,
                                    const std::vector<double>& delta_t,
                                    const double param) = 0;

  virtual Eigen::Matrix<double, -1, -1> coeffA(const double x,
                                               const double delta_t) = 0;

  std::vector<double> detSolver(const std::vector<double>& C_det);

  double multiA(const std::vector<std::vector<double>>& xss,
                const std::vector<std::vector<double>>& delta_ts,
                const double param);

  double tau_;

private:
  Eigen::Matrix<double, -1, -1> formulateMultiAc(
    const std::vector<std::vector<double>>& xss,
    const std::vector<std::vector<double>>& delta_t,
    const double zN);

  void histogram(const std::vector<double>& v,
                 std::vector<size_t>& sort_v,
                 size_t& lower_bound,
                 size_t& upper_bound);

  template<typename T>
  std::vector<size_t> sortIndexes(const std::vector<T>& v)
  {
    auto idx = std::vector<size_t>(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), [&v](const size_t& i1, const size_t& i2) {
      return v[i1] < v[i2];
    });

    return idx;
  }

  double fmin_eps_;
  long   fmin_max_iter_;
  double fmin_init_radius_;

  double find_root_lb_;
  double find_root_ub_;
  double find_root_eps_x_;
  double find_root_eps_val_;
  size_t find_root_max_iter_;
};

}
