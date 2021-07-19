#ifndef OCTOMAP_ROUGH_OCTREE_H
#define OCTOMAP_ROUGH_OCTREE_H

#include <iostream>
#include <boost/dynamic_bitset.hpp>

#include <octomap/OcTreeNode.h>
#include <octomap/OcTreeStamped.h>
#include <octomap/OccupancyOcTreeBase.h>

struct RGBColor{ float r, g, b; };
inline RGBColor HSVtoRGB(double h, double s, double v) {
  RGBColor color;
  int i;
  double m, n, p, f;

  h -= floor(h);
  h *= 6;
  i = floor(h);
  f = h - i;
  m = v * (1 - s);
  n = v * (1 - s * f);
  p = v * (1 - s * (1 - f));

  switch (i) {
    case 6:
    case 0:
      color.r = v; color.g = p; color.b = m;
      break;
    case 1:
      color.r = n; color.g = v; color.b = m;
      break;
    case 2:
      color.r = m; color.g = v; color.b = p;
      break;
    case 3:
      color.r = m; color.g = n; color.b = v;
      break;
    case 4:
      color.r = p; color.g = m; color.b = v;
      break;
    case 5:
      color.r = v; color.g = m; color.b = n;
      break;
    default:
      color.r = 1; color.g = 0.5; color.b = 0.5;
      break;
  }

  return color;
}

inline RGBColor ratioToBW(float ratio)
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

inline RGBColor ratioToRGB(float ratio)
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
    RoughOcTreeNode() : OcTreeNode(), rough(NAN), stair_logodds(0), agent(0) {}

    RoughOcTreeNode(const RoughOcTreeNode& rhs) : OcTreeNode(rhs), rough(rhs.rough), stair_logodds(rhs.stair_logodds), agent(rhs.agent) {}

    bool operator==(const RoughOcTreeNode& rhs) const{
      return (rhs.value == value && rhs.rough == rough && rhs.stair_logodds == stair_logodds && rhs.agent == agent);
    }

    void copyData(const RoughOcTreeNode& from){
      OcTreeNode::copyData(from);
      this->rough =  from.getRough();
      this->stair_logodds =  from.getStairLogOdds();
      this->agent =  from.getAgent();
    }

    inline float getRough() const { return rough; }
    inline void  setRough(float c) {this->rough = c; }

    inline float getStairLogOdds() const { return stair_logodds; }
    inline void  setStairLogOdds(float s) {this->stair_logodds = s; }
    inline float getStairProbability() const { return probability(stair_logodds); }
    inline void updateStairChildren() { this->setStairLogOdds(this->getMaxChildStairLogOdds()); }
    float getMeanChildStairLogOdds() const;
    float getMaxChildStairLogOdds() const;
    void addStairValue(const float& p);

    inline char getAgent() const { return agent; }
    inline void setAgent(char a) { this->agent = a; }

    // Only used for markers right now but can be used for binary/full messages if the agent field is included
    RGBColor getAgentColor(double atZ, double minZ, double maxZ, bool adjustAgent) {
      double h, s, v, sb, vb, sm, vm, split;
      char agent = getAgent();
      if (adjustAgent && (agent > 0)) agent--;
      // Find the standardized height of the voxel
      double z = std::min(std::max((atZ - minZ) / (maxZ - minZ), 0.0), 1.0);

      // Restrict the agents to our preselected colors
      agent = agent % 6;
      sb = 0.2;
      switch (agent) {
        case 0:
          // Black / Green
          h = 0.47;
          sb = 0.1;
          vb = 0.0;
          break;
        case 1:
          // Dark Blue
          h = 0.666;
          vb = 0.55;
          break;
        case 2:
          // Purple
          h = 0.833;
          vb = 0.44;
          break;
        case 3:
          // Green
          h = 0.422;
          vb = 0.53;
          break;
        case 4:
          // Yellow
          h = 0.133;
          vb = 0.48;
          break;
        case 5:
          // Red
          h = 0.0;
          vb = 0.55;
          break;
        case 6:
          // Light Blue
          h = 0.544;
          vb = 0.42;
          break;
      }

      // Multipliers
      sm = 1.0 - sb;
      vm = 1.0 - vb;

      // Build HSV color
      // Center the hue around 0
      h = h + (z - 0.5) * (1.0 / 6.0);

      if (agent == 0) {
        // For merged maps, slowly increase value and decrease saturation to bottom
        s = sb + (1 - z) * sm;
        v = z * z;
      } else {
        // For regular agents, raise saturation, then raise value, then lower saturation
        split = 1.0 / 3.0;
        if (z < split) {
          s = sb + (z / split) * sm;
          v = vb;
        } else if (z < split * 2) {
          s = 1.0;
          v = vb + ((z - split) / split) * vm;
        } else {
          s = sb + (1 - (z - 2 * split) / split) * sm;
          v = 1.0;
        }
      }

      // Convert the HSV to RGB
      return HSVtoRGB(h, s, v);
    }

    RGBColor getRoughColor(/*float rough_=NAN*/) {
      return ratioToBW(getRough());
      // return ratioToRGB(getRough());
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
    float stair_logodds;
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

    // Set stair values for the given node by key or coordinate
    RoughOcTreeNode* setNodeStairValue(const OcTreeKey& key, float logodds);

    RoughOcTreeNode* setNodeStairValue(float x, float y, float z, float logodds) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeStairValue(key,logodds);
    }

    RoughOcTreeNode* setNodeStairValue(point3d pt, float logodds) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(pt, key)) return NULL;
      return setNodeStairValue(key,logodds);
    }

    /// queries whether a node is stairs according to the tree's parameter for "occupancy"
    inline bool isNodeStairs(const RoughOcTreeNode* node) const{
      return (node->getStairLogOdds() >= this->occ_prob_thres_log);
    }

    inline bool isNodeStairs(const RoughOcTreeNode& node) const{
      return (node.getStairLogOdds() >= this->occ_prob_thres_log);
    }

    RoughOcTreeNode* updateNodeStairs(const OcTreeKey& key, float log_odds_update);
    RoughOcTreeNode* updateNodeStairs(const point3d& value, float log_odds_update);
    RoughOcTreeNode* updateNodeStairs(const OcTreeKey& key, bool is_stairs);
    RoughOcTreeNode* updateNodeStairs(const point3d& value, bool is_stairs);
    RoughOcTreeNode* updateNodeStairsRecurs(RoughOcTreeNode* node, const OcTreeKey& key,
                                                    unsigned int depth, const float& log_odds_update);
    void updateNodeStairLogOdds(RoughOcTreeNode* occupancyNode, const float& update) const;

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
