cmake_minimum_required(VERSION 2.8)
project(rough_octomap)

set(PACKAGE_DEPENDENCIES
  roscpp
  visualization_msgs
  std_msgs
  std_srvs
  octomap_ros
  octomap_msgs
  rviz
)

find_package(catkin REQUIRED COMPONENTS ${PACKAGE_DEPENDENCIES})

find_package(octomap REQUIRED)
add_definitions(-DOCTOMAP_NODEBUGOUT)

find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
set(QT_LIBRARIES Qt5::Widgets)
add_definitions(-DQT_NO_KEYWORDS)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${OCTOMAP_INCLUDE_DIRS}
)

QT5_WRAP_CPP(MOC_FILES
  include/rough_octomap/occupancy_grid_display.h
  OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED -DBOOST_LEXICAL_CAST_INCLUDED
)

catkin_package(
  INCLUDE_DIRS include ${catkin_INCLUDE_DIRS} ${OCTOMAP_INCLUDE_DIRS}
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS ${PACKAGE_DEPENDENCIES}
  DEPENDS OCTOMAP
)

set(LINK_LIBS
  ${OCTOMAP_LIBRARIES}
  ${catkin_LIBRARIES}
)

add_library(${PROJECT_NAME} src/RoughOcTree.cpp)
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})

add_library(rough_octomap_rviz_plugin src/occupancy_grid_display.cpp ${MOC_FILES})
target_link_libraries(rough_octomap_rviz_plugin ${PROJECT_NAME} ${LINK_LIBS} ${QT_LIBRARIES})

install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
        PATTERN ".svn" EXCLUDE)

install(FILES plugin_description.xml
        DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)
