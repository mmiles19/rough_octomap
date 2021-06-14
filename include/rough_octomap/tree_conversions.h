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

#ifndef ROUGH_OCTOMAP_CONVERT_TREES_H
#define ROUGH_OCTOMAP_CONVERT_TREES_H

#include <octomap/octomap.h>
#include <octomap_msgs/Octomap.h>
#include <octomap/ColorOcTree.h>
#include <rough_octomap/RoughOcTree.h>
#include <rough_octomap/conversions.h>

// new conversion functions
namespace octomap{
  
  // /**
  //  * @brief Serialization of an octree into binary data e.g. for messages and services.
  //  * Full probability version (stores complete state of tree, .ot file format).
  //  * The data will be much smaller if you call octomap.toMaxLikelihood() and octomap.prune()
  //  * before.
  //  * @return success of serialization
  //  */
  // template <class InputTreeType, class OutputTreeType>
  // static inline bool convert(const InputTreeType& tree_in, OutputTreeType& tree_out){
  //   std::stringstream datastream;
  //   if (!tree_in.write(datastream))
  //     return false;

  //   switch (InputTreeType) {
  //     case RoughOcTree:
  //       convertFromRoughOcTree(tree_in, tree_out)
  //       break;
  //     default:
  //       OCTOMAP_ERROR("Conversion of %s to %s is undefined.", tree_in.getTreeType().c_str(), tree_out.getTreeType().c_str());
  //       break;
  //   }

  //   std::string datastring = datastream.str();
  //   mapData = std::vector<int8_t>(datastring.begin(), datastring.end());

  //   return true;
  // }

  // /**
  //  * @brief Serialization of an octree into binary data e.g. for messages and services.
  //  * Full probability version (stores complete state of tree, .ot file format).
  //  * The data will be much smaller if you call octomap.toMaxLikelihood() and octomap.prune()
  //  * before.
  //  * @return success of serialization
  //  */
  // template <class InputNodeType, class OutputNodeType>
  // static inline bool convert(const InputNodeType& tree_in, OutputNodeType& tree_out){
  //   std::stringstream datastream;
  //   if (!tree_in.write(datastream))
  //     return false;

  //   switch (InputTreeType) {
  //     case RoughOcTree:
  //       convertFromRoughOcTree(tree_in, tree_out)
  //       break;
  //     default:
  //       OCTOMAP_ERROR("Conversion of %s to %s is undefined.", tree_in.getTreeType().c_str(), tree_out.getTreeType().c_str());
  //       break;
  //   }

  //   std::string datastring = datastream.str();
  //   mapData = std::vector<int8_t>(datastring.begin(), datastring.end());

  //   return true;
  // }

  // typedef octomap::convert<octomap::RoughOcTreeNode, octomap::OcTreeNode> RoughOcTreeNode2OcTreeNode;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  // trees

  template <class TREE>
  static inline std::ostream& writeOccupancy(const TREE* tree, std::ostream &s) {
    auto root = tree->getRoot();
    if (root)
      writeOccupancyRecurs(tree, root, s);
    return s;
  }

  template <class TREE, class NODE>
  static inline std::ostream& writeOccupancyRecurs(const TREE* tree, const NODE* node, std::ostream &s) {
    writeOccupancyData(node, s);

    // 1 bit for each children; 0: empty, 1: allocated
    std::bitset<8> children;
    for (unsigned int i=0; i<8; i++) {
      if (tree->nodeChildExists(node, i))
        children[i] = 1;
      else
        children[i] = 0;
    }

    char children_char = (char) children.to_ulong();

    // recursively write children
    for (unsigned int i=0; i<8; i++) {
      if (children[i] == 1) {
        writeNOccupancyRecurs(tree, tree->getNodeChild(node, i), s);
      }
    }

    return s;
  }

  template <class NODE>
  static inline std::ostream& writeOccupancyData(NODE* node, std::ostream &s)  {
    s.write((const char*) &(node->value), sizeof(node->value)); // occupancy
    return s;
  }

  template <class TREE>
  static inline std::istream& readOccupancy(TREE* tree, std::istream &s) {

    if (!s.good()){
      OCTOMAP_WARNING_STR(__FILE__ << ":" << __LINE__ << "Warning: Input filestream not \"good\"");
    }

    tree->tree_size = 0;
    tree->size_changed = true;

    // tree needs to be newly created or cleared externally
    auto root = tree->getRoot();
    if (root) {
      OCTOMAP_ERROR_STR("Trying to read into an existing tree.");
      return s;
    }

    root = new TREE::NodeType();
    readOccupancyRecurs(tree, root, s);

    tree->tree_size = tree->calcNumNodes();  // compute number of nodes
    return s;
  }

  template <class TREE, class NODE>
  static inline std::istream& readOccupancyRecurs(TREE* tree, NODE* node, std::istream &s) {

    readOccupancyData(node, s);

    char children_char;
    s.read((char*)&children_char, sizeof(char));
    std::bitset<8> children ((unsigned long long) children_char);

    for (unsigned int i=0; i<8; i++) {
      if (children[i] == 1){
        NODE* newNode = tree->createNodeChild(node, i);
        readOccupancyRecurs(tree, newNode, s);
      }
    }

    return s;
  }

  template <class NODE>
  static inline std::istream& readOccupancyData(NODE* node, std::istream &s) {
    s.read((char*) &(node->value), sizeof(node->value)); // occupancy
    return s;
  }

  // /**
  //  * @brief Serialization of an octree into binary data e.g. for messages and services.
  //  * Full probability version (stores complete state of tree, .ot file format).
  //  * @return success of serialization
  //  */
  // static inline bool convertOccupancy(const octomap::RoughOcTree& tree_in, octomap::OcTree& tree_out){
  //   std::stringstream datastream;
  //   if (!writeOccupancy(&tree_in, datastream))
  //     return false;
  //   if (!readOccupancy(&tree_out, datastream))
  //     return false;
  //   return true;
  // }

  /**
   * @brief Templated tree conversion. Only encountered when conversion is not defined elsewhere.
   * @return false
   */
  template <class InputTreeType, class OutputTreeType>
  static inline bool convertOccupancy(const InputTreeType& tree_in, OutputTreeType& tree_out){
    std::stringstream datastream;
    if (!writeOccupancy(&tree_in, datastream))
      return false;
    if (!readOccupancy(&tree_out, datastream))
      return false;
    return true;
  }

  static inline bool convertOccupancy(const octomap_msgs::Octomap& tree_msg_in, octomap_msgs::Octomap& tree_msg_out){
    if (tree_msg_in.id=="RoughOcTree" && tree_msg_out.id=="OcTree") {
      RoughOcTree* tree_in = (RoughOcTree*)octomap_msgs::msgToMap(tree_msg_in);
      std::vector<int8_t> mapData;
      if(!fullMapToMsgData(tree_in, mapData))
        return false;
      
    } else {
      OCTOMAP_ERROR("Conversion of %s to %s is undefined.", tree_msg_in.id.c_str(), tree_msg_in.id.c_str());
      return false;
    }
    return true;
  }

  // nodes

  // /**
  //  * @brief Serialization of an octree node into binary data e.g. for messages and services.
  //  * Full probability version (stores complete state of node, .ot file format).
  //  * @return success of serialization
  //  */
  // static inline bool convert(const octomap::RoughOcTreeNode& node_in, octomap::OcTreeNode& node_out){
  //   std::stringstream datastream;
  //   if (!node_in.writeData(datastream))
  //     return false;

  //   switch (InputTreeType) {
  //     case RoughOcTree:
  //       convertFromRoughOcTree(tree_in, tree_out)
  //       break;
  //     default:
  //       OCTOMAP_ERROR("Conversion of %s to %s is undefined.", tree_in.getTreeType().c_str(), tree_out.getTreeType().c_str());
  //       break;
  //   }

  //   std::string datastring = datastream.str();
  //   mapData = std::vector<int8_t>(datastring.begin(), datastring.end());

  //   return true;
  // }

}

#endif
