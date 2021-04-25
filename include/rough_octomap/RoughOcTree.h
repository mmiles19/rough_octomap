#ifndef OCTOMAP_ROUGH_OCTREE_H
#define OCTOMAP_ROUGH_OCTREE_H


#include <iostream>
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
  
  // forward declaraton for "friend"
  class RoughOcTree;
  
  // node definition
  class RoughOcTreeNode : public OcTreeNode {    
  public:
    friend class RoughOcTree; // needs access to node children (inherited)

  public:
    RoughOcTreeNode() : OcTreeNode(), rough(NAN) {}

    RoughOcTreeNode(const RoughOcTreeNode& rhs) : OcTreeNode(rhs), rough(rhs.rough) {}

    bool operator==(const RoughOcTreeNode& rhs) const{
      return (rhs.value == value && rhs.rough == rough);
    }
    
    void copyData(const RoughOcTreeNode& from){
      OcTreeNode::copyData(from);
      this->rough =  from.getRough();
    }
        
    inline float getRough() const { return rough; }
    inline void  setRough(float c) {this->rough = c; }

    float getRough() { return rough; }

    // RGBColor getRoughColor() { return getRGBColor(getRough()); }
    RGBColor getRoughColor(/*float rough_=NAN*/) { 
      // if (isnan(rough_)) rough_=getRough();
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
    
     /**
     * Prunes a node when it is collapsible. This overloaded
     * version only considers the node occupancy for pruning,
     * different colors of child nodes are ignored.
     * @return true if pruning was successful
     */
    virtual bool pruneNode(RoughOcTreeNode* node);
    
    virtual bool isNodeCollapsible(const RoughOcTreeNode* node) const;
       
    // set node color at given key or coordinate. Replaces previous color.
    RoughOcTreeNode* setNodeRough(const OcTreeKey& key, float rough);

    RoughOcTreeNode* setNodeRough(float x, float y, 
                                 float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeRough(key,rough);
    }

    float getNodeRough(const OcTreeKey& key);

    float getNodeRough(float x, float y, float z) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NAN;
      return getNodeRough(key);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    RoughOcTreeNode* averageNodeRough(const OcTreeKey& key, float rough);
    
    RoughOcTreeNode* averageNodeRough(float x, float y, 
                                      float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return averageNodeRough(key,rough);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    RoughOcTreeNode* integrateNodeRough(const OcTreeKey& key, float rough);
    
    RoughOcTreeNode* integrateNodeRough(float x, float y, 
                                      float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return integrateNodeRough(key,rough);
    }

    // update inner nodes, sets color to average child color
    void updateInnerOccupancy();

    // uses gnuplot to plot a RGB histogram in EPS format
    void writeRoughHistogram(std::string filename);

    // binary io overloaded from OcTreeBase
    // std::istream& readBinaryData(std::istream &s);
    // std::ostream& writeBinaryData(std::ostream &s) const;
    std::istream& readBinaryNode(std::istream &s, RoughOcTreeNode* node);
    std::ostream& writeBinaryNode(std::ostream &s, const RoughOcTreeNode* node) const;
    float rough_binary_thres;
    
  protected:
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
  //////////////////// Stamped
namespace octomap {
  
  // forward declaraton for "friend"
  class RoughOcTreeStamped;
  
  // node definition
  class RoughOcTreeNodeStamped : public OcTreeNodeStamped {    
  public:
    friend class RoughOcTreeStamped; // needs access to node children (inherited)

  public:
    RoughOcTreeNodeStamped() : OcTreeNodeStamped(), rough(NAN) {}

    RoughOcTreeNodeStamped(const RoughOcTreeNodeStamped& rhs) : OcTreeNodeStamped(rhs), rough(rhs.rough) {}

    bool operator==(const RoughOcTreeNodeStamped& rhs) const{
      return (rhs.value == value && rhs.timestamp == timestamp && rhs.rough == rough);
    }
    
    void copyData(const RoughOcTreeNodeStamped& from){
      OcTreeNode::copyData(from);
      this->timestamp = from.getTimestamp();
      this->rough =  from.getRough();
    }
        
    inline float getRough() const { return rough; }
    inline void  setRough(float c) {this->rough = c; }

    float getRough() { return rough; }

    // RGBColor getRoughColor() { return getRGBColor(getRough()); }
    RGBColor getRoughColor() { 
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
  };


  // tree definition
  class RoughOcTreeStamped : public OccupancyOcTreeBase <RoughOcTreeNodeStamped> {

  public:
    /// Default constructor, sets resolution of leafs
    RoughOcTreeStamped(double resolution);
      
    /// virtual constructor: creates a new object of same type
    /// (Covariant return type requires an up-to-date compiler)
    RoughOcTreeStamped* create() const {return new RoughOcTreeStamped(resolution); }

    std::string getTreeType() const {return "RoughOcTreeStamped";}
    
     /**
     * Prunes a node when it is collapsible. This overloaded
     * version only considers the node occupancy for pruning,
     * different colors of child nodes are ignored.
     * @return true if pruning was successful
     */
    virtual bool pruneNode(RoughOcTreeNodeStamped* node);
    
    virtual bool isNodeCollapsible(const RoughOcTreeNodeStamped* node) const;
       
    // set node color at given key or coordinate. Replaces previous color.
    RoughOcTreeNodeStamped* setNodeRough(const OcTreeKey& key, float rough);

    RoughOcTreeNodeStamped* setNodeRough(float x, float y, 
                                 float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeRough(key,rough);
    }

    float getNodeRough(const OcTreeKey& key);

    float getNodeRough(float x, float y, float z) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NAN;
      return getNodeRough(key);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    RoughOcTreeNodeStamped* averageNodeRough(const OcTreeKey& key, float rough);
    
    RoughOcTreeNodeStamped* averageNodeRough(float x, float y, 
                                      float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return averageNodeRough(key,rough);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    RoughOcTreeNodeStamped* integrateNodeRough(const OcTreeKey& key, float rough);
    
    RoughOcTreeNodeStamped* integrateNodeRough(float x, float y, 
                                      float z, float rough) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return integrateNodeRough(key,rough);
    }

    // update inner nodes, sets color to average child color
    void updateInnerOccupancy();

    // uses gnuplot to plot a RGB histogram in EPS format
    void writeRoughHistogram(std::string filename);
    
  protected:
    void updateInnerOccupancyRecurs(RoughOcTreeNodeStamped* node, unsigned int depth);

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
           RoughOcTreeStamped* tree = new RoughOcTreeStamped(0.1);
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
    static StaticMemberInitializer roughOcTreeStampedMemberInit;

  };

  //! user friendly output in format (r g b)
  // std::ostream& operator<<(std::ostream& out, RoughOcTree::Color const& c);

} // end namespace

#endif
