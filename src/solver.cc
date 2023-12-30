#include <opengv/math/Sturm.hpp>
#include <dlib/optimization.h>

#include <nhevo/solver.hh>

namespace nhevo {

Eigen::Matrix<double, 4, 4>
Solver::getPose(const std::vector<std::vector<double>>& xss,
                const std::vector<std::vector<double>>& delta_ts,
                const double param,
                const double elapse)
{
  auto pose = Eigen::Matrix<double, 4, 4>();
  pose.setIdentity();

  const auto omega = multiA(xss, delta_ts, param);
  const auto theta = -omega * elapse;
  const auto C = cos(theta);
  const auto S = sin(theta);

  pose(0, 0) = pose(1, 1) = C;
  pose(0, 1) = -S;
  pose(1, 0) = S;

  if (theta > 0) {
    pose(0, 3) = S;
    pose(1, 3) = 1 - C;
  } else {
    pose(0, 3) = -S;
    pose(1, 3) = C - 1;
  }

  return pose;
}

double
Solver::multiA(const std::vector<std::vector<double>>& xss,
               const std::vector<std::vector<double>>& delta_ts,
               const double param)
{
  if (xss.empty()) {
    return 0;
  }

  auto sols = std::vector<double>();
  auto map = std::vector<double>();

  for (auto k = 0; k < xss.size(); k++) {
    const auto& xs = xss[k];
    const auto& delta_t = delta_ts[k];
    auto z = solve(xs, delta_t, param);

    if (!z.empty()) {
      for (auto i = 0; i < z.size(); i++) {
        map.push_back(k);
      }
      sols.insert(sols.end(), z.begin(), z.end());
    }
  }
  if (sols.empty()) {
    std::cerr << "sols is empty \n";
    return 0;
  }

  auto sort_v = std::vector<size_t>();
  size_t lower_bound;
  size_t upper_bound;
  histogram(sols, sort_v, lower_bound, upper_bound);

  double start_point =
    sols[sort_v[static_cast<size_t>((lower_bound + upper_bound) / 2)]];
  const double begin = sols[sort_v[lower_bound]];
  const double end = sols[sort_v[upper_bound]];

  auto minSingular = [&xss, &delta_ts, this](const double zN) {
    const auto A = formulateMultiAc(xss, delta_ts, zN);
    Eigen::BDCSVD<Eigen::MatrixXd> svd(A);
    auto s = svd.singularValues();
    return s.minCoeff();
  };

  dlib::find_min_single_variable(minSingular,
                                 start_point,
                                 begin,
                                 end,
                                 fmin_eps_,
                                 fmin_max_iter_,
                                 fmin_init_radius_);
  return start_point;
}

std::vector<double>
Solver::detSolver(const std::vector<double>& C_det)
{
  std::vector<double> dCoeffs, d2Coeffs, roots, sol;

  dCoeffs.resize(C_det.size() - 1);
  for (size_t i = 1; i < C_det.size(); i++) {
    dCoeffs[i - 1] = C_det[i] * i;
  }

  d2Coeffs.resize(C_det.size() - 2);
  for (size_t i = 1; i < dCoeffs.size(); i++) {
    d2Coeffs[i - 1] = dCoeffs[i] * i;
  }

  std::reverse(dCoeffs.begin(), dCoeffs.end());
  opengv::math::Sturm sturm(dCoeffs);
  sturm.findRoots3(roots,
                    find_root_lb_,
                    find_root_ub_,
                    find_root_eps_x_,
                    find_root_eps_val_,
                    find_root_max_iter_);
  // sturm.findRoots2(roots, 0.0045, 1e-8);

  for (auto z : roots) {
    double value = d2Coeffs[0];
    for (size_t i = 1; i < d2Coeffs.size(); i++) {
      value += d2Coeffs[i] * pow(z, i);
    }

    if (value > 0)
      sol.push_back(z);
  }
  return sol;
}

Eigen::Matrix<double, -1, -1>
Solver::formulateMultiAc(const std::vector<std::vector<double>>& xss,
                         const std::vector<std::vector<double>>& delta_t,
                         const double zN)
{
  const auto inlier_num = xss.size();
  auto rows = 0;
  for (size_t i = 0; i < inlier_num; i++) {
    rows = rows + xss[i].size();
  }
  auto cols = 2 * inlier_num + 1;

  auto multi_ac = Eigen::Matrix<double, -1, -1>();
  multi_ac.resize(rows, cols);
  multi_ac.setZero();

  auto count = 0;
  for (size_t i = 0; i < inlier_num; i++) {
    for (size_t j = 0; j < xss[i].size(); j++) {
      auto cur_C = coeffA(xss[i][j], delta_t[i][j]);

      auto c = Eigen::Matrix<double, 3, 1>();
      c.setZero();
      for (size_t kk = 0; kk < 3; kk++) {
        for (size_t k = 0; k < cur_C.cols(); k++) {
          c(kk) = c(kk) + cur_C(kk, k) * std::pow(zN, k - 1);
        }
      }

      multi_ac(count, 2 * i) = c(0);
      multi_ac(count, 2 * i + 1) = c(1);
      multi_ac(count, multi_ac.cols() - 1) = c(2);

      count++;
    }
  }

  multi_ac = multi_ac / rows;
  return multi_ac;
}

void
Solver::histogram(const std::vector<double>& v,
                  std::vector<size_t>& sort_v,
                  size_t& lower_bound,
                  size_t& upper_bound)
{
  lower_bound = 0;
  upper_bound = 0;
  sort_v = sortIndexes(v);

  const auto size = v.size();
  if (size <= 1) {
    return;
  }

  const auto IQR = v[sort_v[size * 3 / 4]] - v[sort_v[size / 4]];

  const auto min = v[sort_v[0]];
  const auto max = v[sort_v[size - 1]];
  const auto bin_width = 2 * IQR / (std::cbrt(static_cast<double>(size)));
  const int bin_num = std::ceil((max - min) / bin_width);

  if (!IQR || bin_num <= 0) {
    return;
  }

  auto N = std::vector<size_t>();
  auto values = std::vector<std::array<size_t, 3>>(
    bin_num, std::array<size_t, 3>{ 0, size, size });
  for (size_t i = 0; i < size; i++) {
    const size_t index = std::floor((v[sort_v[i]] - min) / bin_width);
    values[index][0]++;

    auto& min_idx = values[index][1];
    auto& max_idx = values[index][2];
    if (min_idx == size) {
      min_idx = i;
    }
    max_idx = i;
  }
  for (const auto& it : values) {
    N.push_back(it[0]);
  }

  const auto sort_N = sortIndexes(N);
  const auto& _1 = sort_N[bin_num - 1];
  const auto& _2 = sort_N[bin_num - 2];

  if (std::abs(static_cast<int>(_1) - static_cast<int>(_2)) == 1 &&
      1.1 * N[_2] > N[_1]) {
    if (_1 < _2) {
      lower_bound = values[_1][1];
      upper_bound = values[_2][2];
    } else {
      lower_bound = values[_2][1];
      upper_bound = values[_1][2];
    }
    return;
  }

  lower_bound = values[_1][1];
  upper_bound = values[_1][2];
}

void
Solver::setTau(const double tau)
{
  tau_ = tau;
}

void
Solver::setFMinBnd(const double fmin_eps,
                   const long fmin_max_iter,
                   const double fmin_init_radius)
{
  fmin_eps_ = fmin_eps;
  fmin_max_iter_ = fmin_max_iter;
  fmin_init_radius_ = fmin_init_radius;
}

void
Solver::setFindRoot(const double find_root_lb,
                    const double find_root_ub,
                    const double find_root_eps_x,
                    const double find_root_eps_val,
                    const double find_root_max_iter)
{
  find_root_lb_ = find_root_lb;
  find_root_ub_ = find_root_ub;
  find_root_eps_x_ = find_root_eps_x;
  find_root_eps_val_ = find_root_eps_val;
  find_root_max_iter_ = find_root_max_iter;
}

}
