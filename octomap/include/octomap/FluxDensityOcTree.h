/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
 * http://octomap.github.com/
 *
 * Copyright (c) 2009-2013, K.M. Wurm and A. Hornung, University of Freiburg
 * All rights reserved.
 * License: New BSD
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


#ifndef OCTOMAP_FLUX_DENSITY_OCTREE_H
#define OCTOMAP_FLUX_DENSITY_OCTREE_H


#include <iostream>
#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>

/* OcTree implementation with annotation of 3d flux density measurements (for 
 * example magnetic field measurements). Heavily boroughed / copied from the 
 * ColorOcTree example.
 *
 * Copyright (c) 2015, Leif Christensen <leif.christensen@gmx.de>
 */

namespace octomap {
  
  // node definition
  class FluxDensityOcTreeNode : public OcTreeNode {    
  public:
    
    class FluxDensity {
    public:
    FluxDensity() : x(0.0), y(0.0), z(0.0){}
    FluxDensity(float _x, float _y, float _z) 
      : x(_x), y(_y), z(_z){}
// TODO equality/unequality operators for float
//      inline bool operator== (const Color &other) const {
//        return (r==other.r && g==other.g && b==other.b);
//      }
//      inline bool operator!= (const Color &other) const {
//        return (r!=other.r || g!=other.g || b!=other.b);
//      }
    float x,y,z;
    };

  public:
    FluxDensityOcTreeNode() : OcTreeNode() {}

    FluxDensityOcTreeNode(const FluxDensityOcTreeNode& rhs) : OcTreeNode(rhs), flux(rhs.flux) {}

//    bool operator==(const ColorOcTreeNode& rhs) const{
//      return (rhs.value == value && rhs.color == color);
//    }
    
    // children
    inline FluxDensityOcTreeNode* getChild(unsigned int i) {
      return static_cast<FluxDensityOcTreeNode*> (OcTreeNode::getChild(i));
    }
    inline const FluxDensityOcTreeNode* getChild(unsigned int i) const {
      return static_cast<const FluxDensityOcTreeNode*> (OcTreeNode::getChild(i));
    }

    bool createChild(unsigned int i) {
      if (children == NULL) allocChildren();
      children[i] = new FluxDensityOcTreeNode();
      return true;
    }

    bool pruneNode();
    void expandNode();
    
    inline FluxDensity getFluxDensity() const { return flux; }
    inline void  setFluxDensity(FluxDensity fd) {this->flux = fd; }
    inline void  setFluxDensity(float x, float y, float z) {
      this->flux = FluxDensity(x,y,z); 
    }

    FluxDensity& getFluxDensity() { return flux; }

    // has any flux density been integrated? TODO all zero not a good measure of being unlikely 
    inline bool isFluxDensitySet() const { 
      return ((flux.x != 0.0) || (flux.y != 0.0) || (flux.z != 0.0)); 
    }

    void updateFluxDensityChildren();


    FluxDensityOcTreeNode::FluxDensity getAverageChildFluxDensity() const;
  
    // file I/O
    std::istream& readValue (std::istream &s);
    std::ostream& writeValue(std::ostream &s) const;
    
  protected:
    FluxDensity flux;
  };


  // tree definition
  class FluxDensityOcTree : public OccupancyOcTreeBase <FluxDensityOcTreeNode> {

  public:
    /// Default constructor, sets resolution of leafs
    FluxDensityOcTree(double resolution) : OccupancyOcTreeBase<FluxDensityOcTreeNode>(resolution) {};  
      
    /// virtual constructor: creates a new object of same type
    /// (Covariant return type requires an up-to-date compiler)
    FluxDensityOcTree* create() const {return new FluxDensityOcTree(resolution); }

    std::string getTreeType() const {return "ColorOcTree";}
   
    // set node color at given key or coordinate. Replaces previous color.
    FluxDensityOcTreeNode* setNodeFluxDensity(const OcTreeKey& key, const float& fdx, const float& fdy, const float& fdz); 

    FluxDensityOcTreeNode* setNodeFluxDensity(const float& x, const float& y, 
                                 const float& z, const float& fdx, const float& fdy, const float& fdz) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeFluxDensity(key,fdx, fdy, fdz);
    }

    // integrate flux density measurement at given key or coordinate. Average with previous flux
    FluxDensityOcTreeNode* averageNodeFluxDensity(const OcTreeKey& key, const float& fdx, 
                                  const float& fdy, const float& fdz);
    
    FluxDensityOcTreeNode* averageNodeFluxDensity(const float& x, const float& y, 
                                      const float& z, const float& fdx, 
                                      const float& fdy, const float& fdz) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return averageNodeFluxDensity(key,fdx,fdy,fdz);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    FluxDensityOcTreeNode* integrateNodeFluxDensity(const OcTreeKey& key, const float& fdx, 
                                      const float& fdy, const float& fdz);
    
    FluxDensityOcTreeNode* integrateNodeFluxDensity(const float& x, const float& y, 
                                      const float& z, const float& fdx, 
                                      const float& fdy, const float& fdz) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return integrateNodeFluxDensity(key, fdx, fdy, fdz);
    }

    // update inner nodes, sets color to average child color
    void updateInnerOccupancy();

    // uses gnuplot to plot a RGB histogram in EPS format
    void writeColorHistogram(std::string filename);
    
  protected:
    void updateInnerOccupancyRecurs(FluxDensityOcTreeNode* node, unsigned int depth);

    /**
     * Static member object which ensures that this OcTree's prototype
     * ends up in the classIDMapping only once
     */
    class StaticMemberInitializer{
       public:
         StaticMemberInitializer() {
           FluxDensityOcTree* tree = new FluxDensityOcTree(0.1);
           AbstractOcTree::registerTreeType(tree);
         }
    };
    /// static member to ensure static initialization (only once)
    static StaticMemberInitializer fluxDensityOcTreeMemberInit;

  };

  //! user friendly output in format (r g b)
  std::ostream& operator<<(std::ostream& out, FluxDensityOcTreeNode::FluxDensity const& c);

} // end namespace

#endif
