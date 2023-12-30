#pragma once

#include "nhevo/trackers/tracker.hh"
#include "nhevo/detectors/detector.hh"

#include <nhevo/solver_s3c2.hh>
#include <nhevo/solver_s5c4.hh>
#include <nhevo/solver_s7c6.hh>

#include <Eigen/Eigen>
#include <memory>
#include <atomic>
#include <vector>

namespace nhevo {

class NHEVO
{
public:
  NHEVO();

  void calculate();

  void setEvent(const int ex,
                const int ey,
                const double et,
                const bool ep);

  void setGroundTruth(const double stamp, Eigen::Matrix4d& ground_truth);

private:
  void declParams();

private:
  int width_;
  int height_;
  Eigen::Matrix<double, 3, 3> K_;

  /* variables */
  std::atomic<double> scale_;
  std::atomic<double> gt_del_t_;
  std::atomic<double> curr_time_;

  Eigen::Matrix4d pose_;
  Eigen::Matrix4d ground_truth_;
  Eigen::Matrix4d rel_ground_truth_;

  /* detector */
  std::unique_ptr<nhevo::Detector> detector_;

  /* tracker */
  std::unique_ptr<nhevo::Tracker> tracker_;

  /* solver */
  std::string solver_type_;
  std::unique_ptr<nhevo::Solver> solver_;

  double solver_tau_;

  double solver_fmin_eps_;
  long   solver_fmin_max_iter_;
  double solver_fmin_init_radius_;

  double solver_find_root_lb_;
  double solver_find_root_ub_;
  double solver_find_root_eps_x_;
  double solver_find_root_eps_val_;
  double solver_find_root_max_iter_;
};

} // namespace nhevo
