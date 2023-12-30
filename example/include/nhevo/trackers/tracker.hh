#pragma once

#include <vector>
#include <mutex>

namespace nhevo {

class Tracker
{
public:
  /**
   * @brief Add event into the tracker.
   *
   * @param[in] et The timestamp of the event, in second.
   * @param[in] ex The x-axis coordinate.
   * @param[in] ey The y-axis coordinate.
   * @param[in] ep The polarity.
   */
  virtual void addEvent(const int ex,
                        const int ey,
                        const double et,
                        const bool ep) noexcept = 0;

  /**
   * @brief Get feature tracks.
   *
   * @param[out] xss The features' x-axis coordinates, one `std::vector<double>` for each track.
   * @param[out] delta_ts The features' timestamps, start from `0` for each vector.
   */
  virtual void getXss(std::vector<std::vector<double>>& xss,
                      std::vector<std::vector<double>>& delta_ts) noexcept = 0;

  virtual ~Tracker() = default;

protected:
  std::mutex mutex_;
};

} // namespace nhevo
