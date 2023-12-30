#include "nhevo/nhevo.hh"
#include <thread>

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <prophesee_event_msgs/EventArray.h>

using PoseStamped = geometry_msgs::PoseStamped;
using EventArray = prophesee_event_msgs::EventArray;

int
main(int argc, char* argv[])
{
  ros::init(argc, argv, "nhevo");
  ros::NodeHandle node;

  auto nhevo = nhevo::NHEVO();

  ros::Subscriber events_subscriber_ =
    node.subscribe<EventArray>("/prophesee/left/events", 100, [&](const EventArray::ConstPtr& msg){
          for (const auto& ev : msg->events) {
            nhevo.setEvent(ev.x,
                           ev.y,
                           ev.ts.toSec(),
                           ev.polarity);
          }
        });

  ros::Subscriber ground_truth_subscriber_ =
    node.subscribe<PoseStamped>("/gt/pose", 10, [&](const PoseStamped::ConstPtr& msg){
          Eigen::Matrix4d m;
          m.setIdentity();

          m.block<3, 3>(0, 0) = Eigen::Quaternion<double>(msg->pose.orientation.w,
                                                          msg->pose.orientation.x,
                                                          msg->pose.orientation.y,
                                                          msg->pose.orientation.z)
                                  .toRotationMatrix();
          m(0, 3) = msg->pose.position.x;
          m(1, 3) = msg->pose.position.y;
          m(2, 3) = msg->pose.position.z;

          nhevo.setGroundTruth(msg->header.stamp.toSec(), m);
        });

  auto thread = std::thread(
    [&]() {
      while (ros::ok()) {
        nhevo.calculate();
      }
    });

  ros::MultiThreadedSpinner spinner(2);
  spinner.spin();

  thread.join();
}
