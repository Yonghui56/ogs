/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateJacobianAssembler.h"
#include "BaseLib/Error.h"

#include "AnalyticalJacobianAssembler.h"
#include "CentralDifferencesJacobianAssembler.h"

namespace ProcessLib
{
std::unique_ptr<AbstractJacobianAssembler> createJacobianAssembler(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
        return std::unique_ptr<AbstractJacobianAssembler>(
            new AnalyticalJacobianAssembler);

    auto const type = config->peekConfigParameter<std::string>("type");

    if (type == "Analytical") {
        config->ignoreConfigParameter("type");
        return std::unique_ptr<AbstractJacobianAssembler>(
            new AnalyticalJacobianAssembler);
    } else if (type == "CentralDifferences") {
        return createCentralDifferencesJacobianAssembler(*config);
    }

    OGS_FATAL("Unknown Jacobian assembler type: `%s'.", type.c_str());
}
}  // ProcessLib
