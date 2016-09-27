/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_RICHARDSFLOW_FEM_H_
#define PROCESS_LIB_RICHARDSFLOW_FEM_H_

#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "RichardsFlowProcessData.h"
#include "MeshLib/CoordinateSystem.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace ProcessLib
{

namespace RichardsFlow
{

const unsigned NUM_NODAL_DOF = 1;

class RichardsFlowLocalAssemblerInterface
        : public ProcessLib::LocalAssemblerInterface
        , public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        std::vector<double>& /*cache*/) const = 0;

};

template <typename ShapeFunction,
         typename IntegrationMethod,
         unsigned GlobalDim>
class LocalAssemblerData
        : public RichardsFlowLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       RichardsFlowProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, _integration_method)),
		  _saturation(std::vector<double>(_integration_method.getNumberOfPoints()))
    {
		assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
		const MeshLib::CoordinateSystem coordsystem(element);
		//const MeshLib::ElementCoordinatesMappingLocal ele_local_coord(element, coordsystem);// wp.getCoords());
		dim = coordsystem.getDimension();
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

		double Sw(0.0);//water saturation
		double Pc(0.0);//capillary pressure
		double k_rel = 0.0;  // relative permeability
		double drhow_dp(0.0);
		double dSwdPc(0.0);
		
		Eigen::MatrixXd mass_mat_coeff = Eigen::MatrixXd::Zero(1, 1);
		Eigen::MatrixXd K_mat_coeff = Eigen::MatrixXd::Zero(1, 1);

		MathLib::PiecewiseLinearInterpolation const&  interP_Pc = *_process_data.curves.at("curveA");
		MathLib::PiecewiseLinearInterpolation const&  interP_Kr = *_process_data.curves.at("curveB");
		
		auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
			local_M_data, local_matrix_size, local_matrix_size);
		auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
		auto local_b = MathLib::createZeroedVector<NodalVectorType>(
			local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());
		double P_int_pt = 0.0;
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
			pos.setIntegrationPoint(ip);
			auto const& sm = _shape_matrices[ip];
			auto const& wp = _integration_method.getWeightedPoint(ip);

			//NumLib::shapeFunctionInterpolate(local_x, sm.N, int_pt_array);
			NumLib::shapeFunctionInterpolate(local_x, sm.N, P_int_pt);
			
			auto const K = _process_data.intrinsic_permeability(t, pos)[0];//permeability
			auto const poro = _process_data.porosity(t, pos)[0];//porosity
			auto const mu = _process_data.viscosity(t, pos)[0];

			Pc = -P_int_pt;
			//Sw = getSwbyPc_van(Pc);

			Sw = interP_Pc.getValue(Pc);//read from Pc-S curve
										//dSwdPc = getdSwdPc_van(Pc);
			_saturation[ip] = Sw;
			//dSwdPc = interP_Pc.getSlope(Pc);//read from slope of Pc-S curve
			//k_rel = getKrelbySw_van(Sw,0);
			dSwdPc = interP_Pc.PressureSaturationDependency(Pc, true);
			k_rel = interP_Kr.getValue(Sw);//read from S-Kr curve

			mass_mat_coeff(0, 0) = storage * Sw + poro * Sw * drhow_dp - poro * dSwdPc;
			K_mat_coeff(0, 0) = K*k_rel / mu;

			local_K.noalias() += sm.dNdx.transpose() *
				K_mat_coeff(0, 0) * sm.dNdx *
				sm.detJ * wp.getWeight();

			//std::cout << t << "  " << _localA << "\n";
			local_M.noalias() += sm.N.transpose() *
				mass_mat_coeff(0, 0) * sm.N *
				sm.detJ * wp.getWeight();//Eigen::Map<Eigen::VectorXd>
										 //std::cout << t << "  " << _localM << "\n";
			if (_process_data.has_gravity) {
				typename ShapeMatricesType::GlobalDimVectorType vec_g;
				vec_g.resize(dim);
				vec_g(dim - 1) = -9.81;
				//const double vec_g(-9.81);//1D
				// since no primary vairable involved
				// directly assemble to the Right-Hand-Side
				// F += dNp^T * K * gz
				local_b.noalias() += sm.dNdx.transpose()*K_mat_coeff(0, 0) * rho_w * vec_g * sm.detJ * wp.getWeight();
			} // end of if hasGravityEffect
			  //std::cout << t << "  " << _localRhs << "\n";
        }//end of GP
		for (int idx_ml = 0; idx_ml < local_M.cols(); idx_ml++)
		{
			double mass_lump_val;
			mass_lump_val = local_M.col(idx_ml).sum();
			local_M.col(idx_ml).setZero();
			local_M(idx_ml, idx_ml) = mass_lump_val;
		}
    }

    Eigen::Map<const Eigen::RowVectorXd>
    getShapeMatrix(const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const&
    getIntPtSaturation(std::vector<double>& /*cache*/) const override
    {
        assert(_saturation.size() > 0);
        return _saturation;
    }

private:
    MeshLib::Element const& _element;
    RichardsFlowProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices> _shape_matrices;

	unsigned dim;

	std::vector<double> _saturation;
        //= std::vector<std::vector<double>>(
            //GlobalDim, std::vector<double>(_integration_method.getNumberOfPoints()));

	const double rho_w = 1000.;//water density
	const double storage = 0.0;
};


}   // namespace RichardsFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_RichardsFlow_FEM_H_
