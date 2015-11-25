/*
 * This file is part of OctoMap - An Efficient Probabilistic 3D Mapping
 * Framework Based on Octrees
 * http://octomap.github.io
 *
 * Copyright (c) 2009-2014, K.M. Wurm and A. Hornung, University of Freiburg
 * All rights reserved. License for the viewer octovis: GNU GPL v2
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 */

#include <octovis/FluxDensityOcTreeDrawer.h>
#include <octomap/FluxDensityOcTree.h>

namespace octomap {

  FluxDensityOcTreeDrawer::FluxDensityOcTreeDrawer() 
    : OcTreeDrawer() {
  }

  FluxDensityOcTreeDrawer::~FluxDensityOcTreeDrawer() {
  }

  void FluxDensityOcTreeDrawer::setOcTree(const AbstractOcTree& tree_pnt,
                                    const octomap::pose6d& origin_,
                                    int map_id_) {

    const FluxDensityOcTree& tree = ((const FluxDensityOcTree&) tree_pnt);

    this->map_id = map_id_;

    // save origin used during cube generation
    this->initial_origin = octomap::pose6d(octomap::point3d(0,0,0), origin_.rot());
    // origin is in global coords
    this->origin = origin_;
    
    // maximum size to prevent crashes on large maps: (should be checked in a better way than a constant)
    bool showAll = (tree.size() < 5 * 1e6);
    bool uses_origin = ( (origin_.rot().x() != 0.) && (origin_.rot().y() != 0.)
                         && (origin_.rot().z() != 0.) && (origin_.rot().u() != 1.) );

    // walk the tree one to find the number of nodes in each category
    // (this is used to set up the OpenGL arrays)
    // TODO: this step may be left out, if we maintained the GLArrays in std::vectors instead...
    unsigned int cnt_occupied(0), cnt_occupied_thres(0), cnt_free(0), cnt_free_thres(0);
    for(FluxDensityOcTree::tree_iterator it = tree.begin_tree(this->m_max_tree_depth),
          end=tree.end_tree(); it!= end; ++it) {
      if (it.isLeaf()) { 
        if (tree.isNodeOccupied(*it)){ // occupied voxels
          if (tree.isNodeAtThreshold(*it)) ++cnt_occupied_thres;
          else                               ++cnt_occupied;
        }
        else if (showAll) { // freespace voxels
          if (tree.isNodeAtThreshold(*it)) ++cnt_free_thres;
          else                               ++cnt_free;
        }
      }        
    }    
    // setup GL arrays for cube quads and cube colors
    initGLArrays(cnt_occupied      , m_occupiedSize     , &m_occupiedArray     , &m_occupiedColorArray);
    initGLArrays(cnt_occupied_thres, m_occupiedThresSize, &m_occupiedThresArray, &m_occupiedThresColorArray);
    initGLArrays(cnt_free          , m_freeSize         , &m_freeArray, NULL);
    initGLArrays(cnt_free_thres    , m_freeThresSize    , &m_freeThresArray, NULL);

    std::vector<octomath::Vector3> cube_template;
    initCubeTemplate(origin, cube_template);

    unsigned int idx_occupied(0), idx_occupied_thres(0), idx_free(0), idx_free_thres(0);
    unsigned int color_idx_occupied(0), color_idx_occupied_thres(0);

    m_grid_voxels.clear();
    OcTreeVolume voxel; // current voxel, possibly transformed 
    for(FluxDensityOcTree::tree_iterator it = tree.begin_tree(this->m_max_tree_depth),
          end=tree.end_tree(); it!= end; ++it) {

      if (it.isLeaf()) { // voxels for leaf nodes
          // Color depending on flux density
          FluxColor c = fluxDensityToRGB(it->getFluxDensity(),0.0,1.0);
        if (uses_origin) 
          voxel = OcTreeVolume(origin.rot().rotate(it.getCoordinate()), it.getSize());
        else 
          voxel = OcTreeVolume(it.getCoordinate(), it.getSize());
        
        if (tree.isNodeOccupied(*it)){ // occupied voxels
          if (tree.isNodeAtThreshold(*it)) {
            idx_occupied_thres = generateCube(voxel, cube_template, idx_occupied_thres, &m_occupiedThresArray);
            color_idx_occupied_thres = setCubeColorHeightmap(voxel, color_idx_occupied_thres, &m_occupiedThresColorArray);
            color_idx_occupied_thres =  setCubeColorRGBA(c.r, c.g, c.b, 
                                                         (unsigned char) (it->getOccupancy() * 255.),
                                                         color_idx_occupied_thres, &m_occupiedThresColorArray);
          }
          else {
            idx_occupied = generateCube(voxel, cube_template, idx_occupied, &m_occupiedArray);
            color_idx_occupied = setCubeColorHeightmap(voxel, color_idx_occupied, &m_occupiedColorArray);
            color_idx_occupied = setCubeColorRGBA(c.r, c.g, c.b, 
                                                  (unsigned char)(it->getOccupancy() * 255.),
                                                  color_idx_occupied, &m_occupiedColorArray);
          }
        }
        //TODO also color freespace voxels using FluxColor from FluxDensity
        else if (showAll) { // freespace voxels
          if (tree.isNodeAtThreshold(*it)) {
            idx_free_thres = generateCube(voxel, cube_template, idx_free_thres, &m_freeThresArray);
          }
          else {
            idx_free = generateCube(voxel, cube_template, idx_free, &m_freeArray);
          }
        }

        // grid structure voxel
        if (showAll) m_grid_voxels.push_back(voxel);        
      }
      
      else { // inner node voxels (for grid structure only)
        if (showAll) {
          if (uses_origin)
            voxel = OcTreeVolume(origin.rot().rotate(it.getCoordinate()), it.getSize());
          else
            voxel = OcTreeVolume(it.getCoordinate(), it.getSize());
          m_grid_voxels.push_back(voxel);
        }
      }      
    } // end for all voxels

    m_octree_grid_vis_initialized = false;

    if(m_drawOcTreeGrid)
      initOctreeGridVis();    
  }

  // Hot-to-Cold colormap after Paul Bourke http://paulbourke.net/texture_color/colourspace
  octomap::FluxColor FluxDensityOcTreeDrawer::fluxDensityToRGB(const octomap::FluxDensityOcTreeNode::FluxDensity& fd,float vmin, float vmax){
     float v = sqrtf(fd.x * fd.x + fd.y * fd.y + fd.z * fd.z); 

     FluxColor c = {1.0,1.0,1.0}; // white
     double dv;

     if (v < vmin)
         v = vmin;
     if (v > vmax)
         v = vmax;
     dv = vmax - vmin;

     if (v < (vmin + 0.25 * dv)) {
         c.r = 0;
         c.g = 4 * (v - vmin) / dv;
     } else if (v < (vmin + 0.5 * dv)) {
         c.r = 0;
         c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
     } else if (v < (vmin + 0.75 * dv)) {
         c.r = 4 * (v - vmin - 0.5 * dv) / dv;
         c.b = 0;
     } else {
         c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
         c.b = 0;
     }
     return c;
  }


} // end namespace
