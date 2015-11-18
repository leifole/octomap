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

#include <octomap/FluxDensityOcTree.h>

namespace octomap {


  // node implementation  --------------------------------------
  std::ostream& FluxDensityOcTreeNode::writeValue (std::ostream &s) const {
    // 1 bit for each children; 0: empty, 1: allocated
    std::bitset<8> children;
    for (unsigned int i=0; i<8; i++) {
      if (childExists(i)) children[i] = 1;
      else                children[i] = 0;
    }
    char children_char = (char) children.to_ulong();
    
    // write node data
    s.write((const char*) &value, sizeof(value)); // occupancy
    s.write((const char*) &flux, sizeof(FluxDensity)); // flux
    s.write((char*)&children_char, sizeof(char)); // child existence

    // write existing children
    for (unsigned int i=0; i<8; ++i) 
      if (children[i] == 1) this->getChild(i)->writeValue(s);    
    return s;
  }

  std::istream& FluxDensityOcTreeNode::readValue (std::istream &s) {
    // read node data
    char children_char;
    s.read((char*) &value, sizeof(value)); // occupancy
    s.read((char*) &flux, sizeof(FluxDensity)); // flux
    s.read((char*)&children_char, sizeof(char)); // child existence

    // read existing children
    std::bitset<8> children ((unsigned long long) children_char);
    for (unsigned int i=0; i<8; i++) {
      if (children[i] == 1){
        createChild(i);
        getChild(i)->readValue(s);
      }
    }
    return s;
  }

  FluxDensityOcTreeNode::FluxDensity FluxDensityOcTreeNode::getAverageChildFluxDensity() const {
    float mx(0.0), my(0.0), mz(0.0);
    int c(0);
    for (int i=0; i<8; i++) {
      if (childExists(i) && getChild(i)->isFluxDensitySet()) {
        mx += getChild(i)->getFluxDensity().x;
        my += getChild(i)->getFluxDensity().y;
        mz += getChild(i)->getFluxDensity().z;
        ++c;
      }
    }
    if (c) {
      mx /= c;
      my /= c;
      mz /= c;
      return FluxDensity(mx, my,mz);
    }
    else { // no child has flux density set
      return FluxDensity(0.0,0.0,0.0);//TODO return FluxDensity with values of type UNSET or similar
    }
  }


  void FluxDensityOcTreeNode::updateFluxDensityChildren() {      
    flux = getAverageChildFluxDensity();
  }

  // pruning =============

  bool FluxDensityOcTreeNode::pruneNode() {
    // checks for equal occupancy only, flux density ignored
    if (!this->collapsible()) return false;
    // set occupancy value 
    setLogOdds(getChild(0)->getLogOdds());
    // set flux density to average flux density
    if (isFluxDensitySet()) flux = getAverageChildFluxDensity();
    // delete children
    for (unsigned int i=0;i<8;i++) {
      delete children[i];
    }
    delete[] children;
    children = NULL;
    return true;
  }

  void FluxDensityOcTreeNode::expandNode() {
    assert(!hasChildren());
    for (unsigned int k=0; k<8; k++) {
      createChild(k);
      children[k]->setValue(value);
      getChild(k)->setFluxDensity(flux);
    }
  }

  // tree implementation  --------------------------------------

  FluxDensityOcTreeNode* FluxDensityOcTree::setNodeFluxDensity(const OcTreeKey& key, 
                                             const float& x, 
                                             const float& y, 
                                             const float& z) {
    FluxDensityOcTreeNode* n = search (key);
    if (n != 0) {
      n->setFluxDensity(x, y, z); 
    }
    return n;
  }

  FluxDensityOcTreeNode* FluxDensityOcTree::averageNodeFluxDensity(const OcTreeKey& key, 
                                                 const float& x, 
                                                 const float& y, 
                                                 const float& z) {
    FluxDensityOcTreeNode* n = search (key);
    if (n != 0) {
      if (n->isFluxDensitySet()) {
        FluxDensityOcTreeNode::FluxDensity prev_flux = n->getFluxDensity();
        n->setFluxDensity((prev_flux.x + x)/2, (prev_flux.y + y)/2, (prev_flux.z + z)/2); 
      }
      else {
        n->setFluxDensity(x, y, z);
      }
    }
    return n;
  }

  FluxDensityOcTreeNode* FluxDensityOcTree::integrateNodeFluxDensity(const OcTreeKey& key, 
                                                   const float& x, 
                                                   const float& y, 
                                                   const float& z) {
    FluxDensityOcTreeNode* n = search (key);
    if (n != 0) {
      if (n->isFluxDensitySet()) {
        FluxDensityOcTreeNode::FluxDensity prev_flux = n->getFluxDensity();
        double node_prob = n->getOccupancy();
        float new_x = (float) ((double) prev_flux.x * node_prob 
                                               +  (double) x * (0.99-node_prob));
        float new_y = (float) ((double) prev_flux.y * node_prob 
                                               +  (double) y * (0.99-node_prob));
        float new_z = (float) ((double) prev_flux.z * node_prob 
                                               +  (double) z * (0.99-node_prob));
        n->setFluxDensity(new_x, new_y, new_z); 
      }
      else {
        n->setFluxDensity(x, y, z);
      }
    }
    return n;
  }
  
  
  void FluxDensityOcTree::updateInnerOccupancy() {
    this->updateInnerOccupancyRecurs(this->root, 0);
  }

  void FluxDensityOcTree::updateInnerOccupancyRecurs(FluxDensityOcTreeNode* node, unsigned int depth) {
    // only recurse and update for inner nodes:
    if (node->hasChildren()){
      // return early for last level:
      if (depth < this->tree_depth){
        for (unsigned int i=0; i<8; i++) {
          if (node->childExists(i)) {
            updateInnerOccupancyRecurs(node->getChild(i), depth+1);
          }
        }
      }
      node->updateOccupancyChildren();
      node->updateFluxDensityChildren();
    }
  }

//  void ColorOcTree::writeColorHistogram(std::string filename) {
//
//#ifdef _MSC_VER
//    fprintf(stderr, "The color histogram uses gnuplot, this is not supported under windows.\n");
//#else
//    // build RGB histogram
//    std::vector<int> histogram_r (256,0);
//    std::vector<int> histogram_g (256,0);
//    std::vector<int> histogram_b (256,0);
//    for(ColorOcTree::tree_iterator it = this->begin_tree(),
//          end=this->end_tree(); it!= end; ++it) {
//      if (!it.isLeaf() || !this->isNodeOccupied(*it)) continue;
//      ColorOcTreeNode::Color& c = it->getColor();
//      ++histogram_r[c.r];
//      ++histogram_g[c.g];
//      ++histogram_b[c.b];
//    }
//    // plot data
//    FILE *gui = popen("gnuplot ", "w");
//    fprintf(gui, "set term postscript eps enhanced color\n");
//    fprintf(gui, "set output \"%s\"\n", filename.c_str());
//    fprintf(gui, "plot [-1:256] ");
//    fprintf(gui,"'-' w filledcurve lt 1 lc 1 tit \"r\",");
//    fprintf(gui, "'-' w filledcurve lt 1 lc 2 tit \"g\",");
//    fprintf(gui, "'-' w filledcurve lt 1 lc 3 tit \"b\",");
//    fprintf(gui, "'-' w l lt 1 lc 1 tit \"\",");
//    fprintf(gui, "'-' w l lt 1 lc 2 tit \"\",");
//    fprintf(gui, "'-' w l lt 1 lc 3 tit \"\"\n");
//
//    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_r[i]);    
//    fprintf(gui,"0 0\n"); fprintf(gui, "e\n");
//    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_g[i]);    
//    fprintf(gui,"0 0\n"); fprintf(gui, "e\n");
//    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_b[i]);    
//    fprintf(gui,"0 0\n"); fprintf(gui, "e\n");
//    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_r[i]);    
//    fprintf(gui, "e\n");
//    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_g[i]);    
//    fprintf(gui, "e\n");
//    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_b[i]);    
//    fprintf(gui, "e\n");
//    fflush(gui);
//#endif
//  }

  std::ostream& operator<<(std::ostream& out, FluxDensityOcTreeNode::FluxDensity const& f) {
    return out << '(' << f.x << ' ' << f.y << ' ' << f.z << ')';
  }


  FluxDensityOcTree::StaticMemberInitializer FluxDensityOcTree::fluxDensityOcTreeMemberInit;

} // end namespace

