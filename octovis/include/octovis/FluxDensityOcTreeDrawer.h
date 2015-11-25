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

#ifndef OCTOVIS_FLUXDENSITY_OCTREEDRAWER_H_
#define OCTOVIS_FLUXDENSITY_OCTREEDRAWER_H_

#include <octovis/OcTreeDrawer.h>
#include <octomap/FluxDensityOcTree.h>

namespace octomap {

    struct FluxColor{
        double  r;
        double  g;
        double  b;
    };

    class FluxDensityOcTreeDrawer : public OcTreeDrawer {
        public:
            FluxDensityOcTreeDrawer();
            virtual ~FluxDensityOcTreeDrawer();

            virtual void setOcTree(const AbstractOcTree& tree_pnt, const pose6d& origin, int map_id_);
            FluxColor fluxDensityToRGB(const FluxDensityOcTreeNode::FluxDensity& fd,float min, float max);

        protected:


    };


} // end namespace

#endif
