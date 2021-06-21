#ifndef OCTOMAP_ROUGH_OCTREE_H
#define OCTOMAP_ROUGH_OCTREE_H

#include <iostream>
#include <boost/dynamic_bitset.hpp>

#include <octomap/OcTreeNode.h>
#include <octomap/OcTreeStamped.h>
#include <octomap/OccupancyOcTreeBase.h>

struct RGBColor{ float r, g, b; };
inline RGBColor getRGBColor(float ratio)
{
    RGBColor out;
    if(isnan(ratio))
    {
      out.r = 0.0;
      out.g = 0.0;
      out.b = 0.0;
    }
    else
    {
      //we want to normalize ratio so that it fits in to 6 regions
      //where each region is 256 units long
      int normalized = int(ratio * 255 * 5);

      //find the distance to the start of the closest region
      int x = normalized % 256;

      int red = 0, grn = 0, blu = 0;
      switch(normalized / 256)
      {
          case 0: red = 255;      grn = x;        blu = 0;       break;//red
          case 1: red = 255 - x;  grn = 255;      blu = 0;       break;//yellow
          case 2: red = 0;        grn = 255;      blu = x;       break;//green
          case 3: red = 0;        grn = 255 - x;  blu = 255;     break;//cyan
          case 4: red = x;        grn = 0;        blu = 255;     break;//blue
      }

      out.r = (float)red/255.0;
      out.g = (float)grn/255.0;
      out.b = (float)blu/255.0;
    }
    return out;
}
inline RGBColor getBWColor(float ratio)
{
    RGBColor out;
    if(isnan(ratio))
    {
      out.r = 1.0;
      out.g = 0.0;
      out.b = 0.0;
    }
    else
    {
      out.r = ratio;
      out.g = ratio;
      out.b = ratio;
    }
    return out;
}

namespace octomap {
  enum RoughBinaryEncodingMode {
    THRESHOLDING,
    BINNING
  };
}

namespace octomap {

  // forward declaraton for "friend"
  class RoughOcTree;

  // node definition
  class RoughOcTreeNode : public OcTreeNode {
  public:
    friend class RoughOcTree; // needs access to node children (inherited)

  public:
    RoughOcTreeNode() : OcTreeNode(), rough(NAN), agent(0) {}

    RoughOcTreeNode(const RoughOcTreeNode& rhs) : OcTreeNode(rhs), rough(rhs.rough), agent(rhs.agent) {}

    bool operator==(const RoughOcTreeNode& rhs) const{
      return (rhs.value == value && rhs.rough == rough && rhs.agent == agent);
    }

    void copyData(const RoughOcTreeNode& from){
      OcTreeNode::copyData(from);
      this->rough =  from.getRough();
      this->agent =  from.getAgent();
    }

    inline float getRough() const { return rough; }
    inline void  setRough(float c) {this->rough = c; }

    inline char getAgent() const { return agent; }
    inline void setAgent(char a) { this->agent = a; }

    // Potential future use but would have to publish the agent data in the messages, which we don't currently
    RGBColor getAgentColor(/*float rough_=NAN*/) {
      return getBWColor(getRough());
    }

    RGBColor getRoughColor(/*float rough_=NAN*/) {
      return getBWColor(getRough());
    }

    // has any color been integrated? (pure white is very unlikely...)
    inline bool isRoughSet() const {
      return !isnan(rough);
    }

    void updateRoughChildren();

    float getAverageChildRough() const;

    // file I/O
    std::istream& readData(std::istream &s);
    std::ostream& writeData(std::ostream &s) const;

  protected:
    float rough;
    char agent;
  };


  // tree definition
  class RoughOcTree : public OccupancyOcTreeBase <RoughOcTreeNode> {

  public:
    /// Default constructor, sets resolution of leafs
    RoughOcTree(double resolution);

    /// virtual constructor: creates a new object of same type
    /// (Covariant return type requires an up-to-date compiler)
    RoughOcTree* create() const {return new RoughOcTree(resolution); }

    std::string getTreeType() const {return "RoughOcTree";}

    inline bool getRoughEnabled() const { return roughEnabled; }
    inline void setRoughEnabled(bool e) {
      this->roughEnabled = e;
      if (!e) this->num_binary_bins = 0;
    }

    inline uint getNumBins() const { return num_binary_bins; }
    inline void setNumBins(uint n) { this->num_binary_bins = n; }

     /**
     * Prunes a node when it is collapsible. This overloaded
     * version only considers the node occupancy for pruning,
     * different colors of child nodes are ignored.
     * @return true if pruning was successful
     */
    virtual bool pruneNode(RoughOcTreeNode* node);

    virtual bool isNodeCollapsible(const RoughOcTreeNode* node) const;

    RoughOcTreeNode* updateNodeRough(RoughOcTreeNode* node, const OcTreeKey& key, bool occupied, char agent);

    // Set agent for the given node by key or coordinate
    RoughOcTreeNode* setNodeAgent(const OcTreeKey& key, char agent);

    RoughOcTreeNode* setNodeAgent(float x, float y,
                                 float z, char agent) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeAgent(key,agent);
    }

    RoughOcTreeNode* setNodeAgent(point3d pt, char agent) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(pt, key)) return NULL;
      return setNodeAgent(key,agent);
    }

    // set node roughness at given key or coordinate. Replaces previous roughness.
    RoughOcTreeNode* setNodeRough(const OcTreeKey& key, float rough);

    RoughOcTreeNode* setNodeRough(float x, float y,
                                 float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeRough(key,rough);
    }

    RoughOcTreeNode* setNodeRough(point3d pt, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(pt, key)) return NULL;
      return setNodeRough(key,rough);
    }

    float getNodeRough(const OcTreeKey& key);

    float getNodeRough(float x, float y, float z) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NAN;
      return getNodeRough(key);
    }

    float getNodeRough(point3d pt) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(pt, key)) return NAN;
      return getNodeRough(key);
    }

    // integrate roughness measurement at given key or coordinate. Average with previous roughness
    RoughOcTreeNode* averageNodeRough(const OcTreeKey& key, float rough);

    RoughOcTreeNode* averageNodeRough(float x, float y,
                                      float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return averageNodeRough(key,rough);
    }

    RoughOcTreeNode* averageNodeRough(point3d pt, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(pt, key)) return NULL;
      return averageNodeRough(key,rough);
    }

    // integrate roughness measurement at given key or coordinate. Average with previous roughness
    RoughOcTreeNode* integrateNodeRough(const OcTreeKey& key, float rough);

    RoughOcTreeNode* integrateNodeRough(float x, float y,
                                      float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return integrateNodeRough(key,rough);
    }

    RoughOcTreeNode* integrateNodeRough(point3d pt, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(pt, key)) return NULL;
      return integrateNodeRough(key,rough);
    }

    // update inner nodes, sets roughness to average child roughness
    void updateInnerOccupancy();

    // uses gnuplot to plot a RGB histogram in EPS format
    void writeRoughHistogram(std::string filename);

    // binary io overloaded from OcTreeBase
    std::istream& readBinaryData(std::istream &s);
    std::ostream& writeBinaryData(std::ostream &s) const;
    std::istream& readBinaryNode(std::istream &s, RoughOcTreeNode* node);
    std::ostream& writeBinaryNode(std::ostream &s, const RoughOcTreeNode* node) const;
    std::istream& readBinaryNodeViaThresholding(std::istream &s, RoughOcTreeNode* node);
    std::ostream& writeBinaryNodeViaThresholding(std::ostream &s, const RoughOcTreeNode* node) const;
    std::istream& readBinaryNodeViaBinning(std::istream &s, RoughOcTreeNode* node);
    std::ostream& writeBinaryNodeViaBinning(std::ostream &s, const RoughOcTreeNode* node) const;

    RoughBinaryEncodingMode binary_encoding_mode;
    float rough_binary_thres; // must be between 0 and 1
    uint num_binary_bins; // must be power of 2

  protected:
    bool roughEnabled = false;

    void updateInnerOccupancyRecurs(RoughOcTreeNode* node, unsigned int depth);

    /**
     * Static member object which ensures that this OcTree's prototype
     * ends up in the classIDMapping only once. You need this as a
     * static member in any derived octree class in order to read .ot
     * files through the AbstractOcTree factory. You should also call
     * ensureLinking() once from the constructor.
     */
    class StaticMemberInitializer{
       public:
         StaticMemberInitializer() {
           RoughOcTree* tree = new RoughOcTree(0.1);
           tree->clearKeyRays();
           AbstractOcTree::registerTreeType(tree);
         }

         /**
         * Dummy function to ensure that MSVC does not drop the
         * StaticMemberInitializer, causing this tree failing to register.
         * Needs to be called from the constructor of this octree.
         */
         void ensureLinking() {};
    };
    /// static member to ensure static initialization (only once)
    static StaticMemberInitializer roughOcTreeMemberInit;

  };
}

#endif
