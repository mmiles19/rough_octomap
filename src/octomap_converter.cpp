#include <ros/ros.h>
#include <rough_octomap/tree_conversions.h>

ros::Publisher out_pub;

void callback(const octomap_msgs::Octomap& tree_msg_in)
{
    octomap_msgs::Octomap tree_msg_out;
    if(octomap::convertOccupancy(tree_msg_in, tree_msg_out))
        out_pub.publish(tree_msg_out);
}

int main (int argc, char *argv[])
{
    ros::init(argc, argv, "octomap_converter");
    ros::NodeHandle n("~");
    ros::Subscriber in_sub = n.subscribe("input",1,&callback);
    out_pub = n.advertise<octomap_msgs::Octomap>("output", 1);
    ros::spin();
}