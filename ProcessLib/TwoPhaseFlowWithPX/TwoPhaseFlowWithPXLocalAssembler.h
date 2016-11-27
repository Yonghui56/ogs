/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEFLOWWITHPXLOCALASSEMBLER_H
#define OGS_TWOPHASEFLOWWITHPXLOCALASSEMBLER_H

#include <vector>
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "TwoPhaseFlowWithPXProcessData.h"
#include "EoSBase.h"
struct IntegrationPointData final
{
	explicit IntegrationPointData(
		ProcessLib::TwoPhaseFlowWithPX::EoSBase& EoS_material)
		: _EoS_material(EoS_material)
	{
	}

	double _sw, _sw_pre;
	double _x_m, _x_m_prev;

	ProcessLib::TwoPhaseFlowWithPX::EoSBase& _EoS_material;

	double _detJ;
	double _integralMeasure;
	//copy the previous to current
	void pushBackState()
	{
		_sw_pre = _sw;
		_x_m_prev = _x_m;
	}
};
namespace ProcessLib
{
namespace TwoPhaseFlowWithPX
{
const unsigned NUM_NODAL_DOF = 2;

class TwoPhaseFlowWithPXLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtWettingPressure(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class TwoPhaseFlowWithPXLocalAssembler
    : public TwoPhaseFlowWithPXLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using LocalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using LocalVectorType = typename LocalAssemblerTraits::LocalVector;

public:
    TwoPhaseFlowWithPXLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        TwoPhaseFlowWithPXProcessData const& process_data)
        : _element(element),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _process_data(process_data),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_wetting(
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const override
    {
        assert(_saturation.size() > 0);
        return _saturation;
    }

    std::vector<double> const& getIntPtWettingPressure(
        std::vector<double>& /*cache*/) const override
    {
        assert(_pressure_wetting.size() > 0);
        return _pressure_wetting;
    }

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices> _shape_matrices;

    TwoPhaseFlowWithPXProcessData const& _process_data;
	std::vector<IntegrationPointData> _ip_data;
    // Note: currently only isothermal case is considered, so the temperature is
    // assumed to be const
    // the variation of temperature will be taken into account in future
	const double molar_mass_h2o =
		MaterialLib::PhysicalConstant::MolarMass::Water;
	const double molar_mass_h2 = MaterialLib::PhysicalConstant::MolarMass::H2;
	double _temperature = 293.15;
    std::vector<double> _saturation;
    std::vector<double> _pressure_wetting;
    static const int nonwet_pressure_coeff_index = 0;
    static const int cap_pressure_coeff_index = 1;

    static const int nonwet_pressure_matrix_index = 0;
    static const int cap_pressure_matrix_index = ShapeFunction::NPOINTS;

    static const int nonwet_pressure_size = ShapeFunction::NPOINTS;
    static const int cap_pressure_size = ShapeFunction::NPOINTS;
};

}  // end of namespace
}  // end of namespace

#include "TwoPhaseFlowWithPXLocalAssembler-impl.h"

#endif /* TWOPHASEFLOWWITHPXLOCALASSEMBLER_H */
