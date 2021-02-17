/**
 * OctoMap ROS message conversions / [de-] serialization
 *
 * @author A. Hornung, University of Freiburg, Copyright (C) 2011-2013.
 * @see http://www.ros.org/wiki/octomap_ros
 * License: BSD
 */

/*
 * Copyright (c) 2011-2013, A. Hornung, University of Freiburg
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ROUGH_OCTOMAP_MSGS_CONVERT_MSGS_H
#define ROUGH_OCTOMAP_MSGS_CONVERT_MSGS_H

#include <octomap/octomap.h>
#include <octomap_msgs/Octomap.h>
#include <octomap/ColorOcTree.h>
#include <octomap/RoughOcTree.h>

// new conversion functions  
namespace octomap_msgs{  
  /**
   * @brief Creates a new octree by deserializing from msg,
   * e.g. from a message or service (binary: only free and occupied .bt file format).
   * This creates a new OcTree object and returns a pointer to it.
   * You will need to free the memory when you're done.
   */
  static inline octomap::AbstractOcTree* binaryMsgToMap(const Octomap& msg){
    if (!msg.binary)
      return NULL;

    octomap::AbstractOcTree* tree;
    if (msg.id == "ColorOcTree"){
      octomap::ColorOcTree* octree = new octomap::ColorOcTree(msg.resolution);    
      readTree(octree, msg);
      tree = octree;
    }
    else if (msg.id == "RoughOcTree"){
      octomap::RoughOcTree* octree = new octomap::RoughOcTree(msg.resolution);    
      readTree(octree, msg);
      tree = octree;
    }
    else if (msg.id == "RoughOcTreeStamped"){
      octomap::RoughOcTreeStamped* octree = new octomap::RoughOcTreeStamped(msg.resolution);    
      readTree(octree, msg);
      tree = octree;
    } else {
      octomap::OcTree* octree = new octomap::OcTree(msg.resolution);    
      readTree(octree, msg);
      tree = octree;
    }
    return tree;      
  }
}


#endif

