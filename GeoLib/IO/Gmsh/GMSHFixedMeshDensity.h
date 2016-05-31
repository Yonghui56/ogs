/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GMSHFIXEDMESHDENSITY_H_
#define GMSHFIXEDMESHDENSITY_H_

#include "GMSHMeshDensityStrategy.h"

namespace GeoLib
{
namespace IO
{
namespace GMSH
{

class GMSHFixedMeshDensity : public GMSHMeshDensityStrategy
{
public:
    GMSHFixedMeshDensity(double mesh_density);
    void init(std::vector<GeoLib::Point const*> const& vec);
    double getMeshDensityAtPoint(GeoLib::Point const*const) const;
    double getMeshDensityAtStation(GeoLib::Point const*const) const;
    virtual ~GMSHFixedMeshDensity() {}

private:
    double _mesh_density;
};

}
} // end namespace IO
} // end namespace GeoLib

#endif /* GMSHFIXEDMESHDENSITY_H_ */
