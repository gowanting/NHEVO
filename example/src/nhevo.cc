#include "nhevo/nhevo.hh"
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <thread>


namespace nhevo {

NHEVO::NHEVO()
{
  declParams();

  scale_ = 0;
  pose_.setIdentity();
  ground_truth_.setIdentity();

  tracker_ = nullptr;
  detector_ = nullptr;

  if (solver_type_ == "s3c2") {
    solver_ = std::make_unique<nhevo::SolverS3C2>();
  } else if (solver_type_ == "s5c4") {
    solver_ = std::make_unique<nhevo::SolverS5C4>();
  } else if (solver_type_ == "s7c6") {
    solver_ = std::make_unique<nhevo::SolverS7C6>();
  }

  solver_->setTau(solver_tau_);

  solver_->setFMinBnd(solver_fmin_eps_,
                      solver_fmin_max_iter_,
                      solver_fmin_init_radius_);

  solver_->setFindRoot(solver_find_root_lb_,
                       solver_find_root_ub_,
                       solver_find_root_eps_x_,
                       solver_find_root_eps_val_,
                       solver_find_root_max_iter_);
}

void
NHEVO::setEvent(const int ex,
                const int ey,
                const double et,
                const bool ep)
{
  curr_time_ = et;

  if (detector_->isCorner(ex, ey, et, ep)) {
    tracker_->addEvent(ex, ey, et, ep);
  }
}

void
NHEVO::setGroundTruth(const double stamp, Eigen::Matrix4d& ground_truth)
{
  if (ground_truth.isIdentity()) {
    return;
  }

  curr_time_ = stamp;

  static double last = stamp;
  gt_del_t_ = stamp - last;

  rel_ground_truth_ = ground_truth_.inverse() * ground_truth;
  scale_ = rel_ground_truth_.block<3, 1>(0, 3).norm();

  last = stamp;
  ground_truth_ = std::move(ground_truth);
}

void
NHEVO::calculate()
{
  if (gt_del_t_ < 1e-2) {
    return;
  }

  auto xss = std::vector<std::vector<double>>();
  auto delta_ts = std::vector<std::vector<double>>();

  tracker_->getXss(xss, delta_ts);

  static double last   = curr_time_;
  const  double elapse = curr_time_ - last;
  last = curr_time_;

  auto pose = solver_->getPose(xss, delta_ts, 1e11, elapse);

  auto&& t = pose.block<3, 1>(0, 3);
  const auto norm = t.norm();

  if (abs(norm) > 1e-5) {
    t = t * scale_ / norm / (gt_del_t_ / elapse);
  } else {
    t(0) = scale_ / (gt_del_t_ / elapse);
    t(1) = 0;
  }

  pose_ = pose_ * pose;

  auto error = (pose.inverse() * rel_ground_truth_ - Eigen::Matrix4d::Identity()).norm();
  std::cout << "error: " << error << "\n";

  std::this_thread::sleep_for(std::chrono::milliseconds(20));
}

void
NHEVO::declParams()
{
  auto node = YAML::LoadFile("./config/config.yaml");

  width_ = node["general"]["width"].as<int>();
  height_ = node["general"]["height"].as<int>();

  std::vector<std::vector<double>> vK;
  vK = node["general"]["camera_matrix"].as<decltype(vK)>();
  K_ << vK[0][0], vK[0][1], vK[0][2], vK[1][0], vK[1][1], vK[1][2],
        vK[2][0], vK[2][1], vK[2][2];

  solver_type_ = node["solver"]["type"].as<std::string>();

  solver_tau_ = node["solver"]["tau"].as<double>();

  solver_fmin_eps_         = node["solver"]["fmin_eps"].as<double>();
  solver_fmin_max_iter_    = node["solver"]["fmin_max_iter"].as<long>();
  solver_fmin_init_radius_ = node["solver"]["fmin_init_radius"].as<double>();

  solver_find_root_lb_       = node["solver"]["find_root_lb"].as<double>();
  solver_find_root_ub_       = node["solver"]["find_root_ub"].as<double>();
  solver_find_root_eps_x_    = node["solver"]["find_root_eps_x"].as<double>();
  solver_find_root_eps_val_  = node["solver"]["find_root_eps_val"].as<double>();
  solver_find_root_max_iter_ = node["solver"]["find_root_max_iter"].as<double>();
}

} // namespace nhevo
