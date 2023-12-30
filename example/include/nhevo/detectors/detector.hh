#pragma once

namespace nhevo {

class Detector
{
public:
  /**
   * @brief Determine if the event is a feature point.
   *
   * @param[in] et The timestamp of the event, in second.
   * @param[in] ex The x-axis coordinate.
   * @param[in] ey The y-axis coordinate.
   * @param[in] ep The polarity.
   *
   * @return bool.
   */
  virtual bool isCorner(const int ex,
                        const int ey,
                        const double et,
                        const bool ep) = 0;
  virtual ~Detector() = default;
};

} // namespace nhevo
