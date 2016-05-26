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

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "RichardsFlowProcessData.h"

namespace ProcessLib
{

namespace RichardsFlow
{
    enum class IntegrationPointValue {
        Saturation
    };

    const unsigned NUM_NODAL_DOF = 1;
	
    template <typename GlobalMatrix, typename GlobalVector>
    class RichardsFlowLocalAssemblerInterface
        : public ProcessLib::LocalAssemblerInterface<GlobalMatrix, GlobalVector>
        , public NumLib::Extrapolatable<GlobalVector, IntegrationPointValue>
	{
	public:
		virtual void postSetInitialConditions(
			std::vector<double> const& local_x) = 0;
	};
	
template <typename ShapeFunction,
         typename IntegrationMethod,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData : public RichardsFlowLocalAssemblerInterface<GlobalMatrix, GlobalVector>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;

public:
    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       RichardsFlowProcessData const& process_data)
        : _element(element)
        , _shape_matrices(
              initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod, GlobalDim>(
                  element, integration_order))
        , _process_data(process_data)
        , _localM(local_matrix_size, local_matrix_size)
        , _localA(local_matrix_size, local_matrix_size) // TODO narrowing conversion
        , _localRhs(local_matrix_size)
		, _localJac(local_matrix_size, local_matrix_size)
        , _integration_order(integration_order)
    {
		
        const MeshLib::CoordinateSystem coordsystem(element);
        const MeshLib::ElementCoordinatesMappingLocal ele_local_coord(element, coordsystem);// wp.getCoords());
		const int node = element.getNNodes();
		dim = coordsystem.getDimension();
		//local_x0=
    }
	void postSetInitialConditions(std::vector<double> const& local_x) override
	{
		local_x0 = local_x;
	}
    void assemble(double const t, std::vector<double> const& local_x) override
    {
        _localA.setZero();
        _localM.setZero();
        _localRhs.setZero();
        const double rho_w = 1000.;//water density
        const double storage = 0.0;
        double Sw(0.0);//water saturation
        double Pc(0.0);//capillary pressure
        double k_rel = 0.0;  // relative permeability
        double drhow_dp(0.0);
        double dSwdPc(0.0);

        /*if (hasGravityEffect) {
            vec_g = LocalVectorType::Zero(_problem_coordinates.getDimension());
            vec_g[_problem_coordinates.getIndexOfY()] = -9.81;
        }
        */
        Eigen::MatrixXd mass_mat_coeff = Eigen::MatrixXd::Zero(1, 1);
        Eigen::MatrixXd K_mat_coeff = Eigen::MatrixXd::Zero(1, 1);
        
        MathLib::PiecewiseLinearInterpolation const&  interP_Pc = *_process_data.curves.at("curveA");
        MathLib::PiecewiseLinearInterpolation const&  interP_Kr = *_process_data.curves.at("curveB");

        //const bool hasGravityEffect = false;

        IntegrationMethod integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();//return gauss point number

        double P_int_pt = 0.0;
        std::array<double*, 1> const int_pt_array = { &P_int_pt };

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);

            //NumLib::shapeFunctionInterpolate(local_x, sm.N, int_pt_array);
            NumLib::shapeFunctionInterpolate(local_x, sm.N, P_int_pt);
            Pc = -P_int_pt;
            //Sw = getSwbyPc_van(Pc);
            Sw = interP_Pc.getValue(Pc);//read from Pc-S curve
                                        //dSwdPc = getdSwdPc_van(Pc);
            _saturation[0][ip] = Sw;
            //dSwdPc = interP_Pc.getSlope(Pc);//read from slope of Pc-S curve
                                            //k_rel = getKrelbySw_van(Sw,0);
            dSwdPc = interP_Pc.PressureSaturationDependency(Pc,true);
            k_rel = interP_Kr.getValue(Sw);//read from S-Kr curve

            mass_mat_coeff(0, 0) = storage * Sw + _process_data.porosity(_element) * Sw * drhow_dp - _process_data.porosity(_element) * dSwdPc;
            K_mat_coeff(0, 0) = _process_data.intrinsic_permeability(_element)*k_rel / _process_data.viscosity(_element);

            _localA.noalias() += sm.dNdx.transpose() *
                K_mat_coeff(0, 0) * sm.dNdx *
                sm.detJ * wp.getWeight();

            //std::cout << t << "  " << _localA << "\n";
            _localM.noalias() += sm.N *
                mass_mat_coeff(0, 0) * sm.N.transpose() *
                sm.detJ * wp.getWeight();//Eigen::Map<Eigen::VectorXd>
            //std::cout << t << "  " << _localM << "\n";
            if (_process_data.has_gravity) {

                //Eigen::Vector3d const vec_g(0, 0, -9.81); //2D X-Z 
                //Eigen::Vector2d const vec_g(0, -9.81);//2D X-Y
                typename ShapeMatricesType::GlobalDimVectorType vec_g;
                vec_g.resize(dim);
                vec_g(dim-1) = -9.81;
                //const double vec_g(-9.81);//1D
                // since no primary vairable involved
                // directly assemble to the Right-Hand-Side
                // F += dNp^T * K * gz
                _localRhs.noalias() += sm.dNdx.transpose()*K_mat_coeff(0, 0) * rho_w * vec_g * sm.detJ * wp.getWeight();
            } // end of if hasGravityEffect
            //std::cout << t << "  " << _localRhs << "\n";
        }
    }
    
	void assembleJacobian(double const t,
        std::vector<double> const& local_x) override
    {
        size_t el(0);
        const size_t n_nodes = _element.getNNodes();

        const size_t n_dof = local_x.size();
        
        Eigen::MatrixXd TMP_M_1 = Eigen::MatrixXd::Zero(n_dof, n_dof);
        Eigen::MatrixXd TMP_M_2 = Eigen::MatrixXd::Zero(n_dof, n_dof);
        
		std::vector<double> _local_x(n_dof);
		_local_x = local_x;
        double eps = 1e-7;
        double _Theta=0.5;
        for (size_t u_idx = 0; u_idx < n_dof; u_idx++)
        {
            //clear M1,K1,F1,M2,K2,F2
            Eigen::MatrixXd M1 = Eigen::MatrixXd::Zero(n_dof, n_dof);
            Eigen::MatrixXd K1 = Eigen::MatrixXd::Zero(n_dof, n_dof);
            Eigen::VectorXd F1 = Eigen::VectorXd::Zero(n_dof);
            Eigen::MatrixXd M2 = Eigen::MatrixXd::Zero(n_dof, n_dof);
            Eigen::MatrixXd K2 = Eigen::MatrixXd::Zero(n_dof, n_dof);
            Eigen::VectorXd F2 = Eigen::VectorXd::Zero(n_dof);
            ////clear local_r
            Eigen::VectorXd local_r_1 = Eigen::VectorXd::Zero(n_dof);
            Eigen::VectorXd local_r_2 = Eigen::VectorXd::Zero(n_dof);
            ////clear the epsilon vector 
            std::vector<double> e_vec(n_dof);
			Eigen::VectorXd e_vec_ = Eigen::VectorXd::Zero(n_dof);
			//Eigen::Vector4d local_x_(local_x.data());
			Eigen::Map<const Eigen::VectorXd> local_x_(local_x.data(), local_x.size());
			Eigen::Map<const Eigen::VectorXd> local_x0_(local_x0.data(), local_x0.size());
			//Eigen::Map<Eigen::VectorXd> local_x_map(local_x.data(), local_x.size());
            e_vec[u_idx] = eps*(1 + std::abs(local_x[u_idx]));
			e_vec_[u_idx] = eps*(1 + std::abs(local_x[u_idx]));
			for (size_t n_idx = 0; n_idx < n_dof; n_idx++)
				_local_x[n_idx] = local_x[n_idx] - e_vec[n_idx];
            assemble(t, _local_x);
			//assemble(t, local_x_-e_vec_);//second argument should be std::vec
            TMP_M_1 = (1.0 /  _process_data.dt)*_localM + _Theta*_localA;
			
            local_r_1 = TMP_M_1*(local_x_ - e_vec_);
			//local_r_1 = TMP_M_1*_local_x;
            TMP_M_1 = (1.0 /  _process_data.dt) * _localM - (1. - _Theta) * _localA;
			local_r_1.noalias() -= TMP_M_1 * local_x0_;//u0
            local_r_1.noalias() -= _localRhs;
			for (size_t n_idx = 0; n_idx < n_dof; n_idx++)
				_local_x[n_idx] = local_x[n_idx] + e_vec[n_idx];
            assemble(t, _local_x);
			TMP_M_2 = (1.0 / _process_data.dt)*_localM + _Theta*_localA;
            local_r_2 = TMP_M_2*(local_x_ + e_vec_);
			TMP_M_2 = (1.0 / _process_data.dt) * _localM - (1. - _Theta) * _localA;
			local_r_2.noalias() -= TMP_M_2 *local_x0_;//u0
            local_r_2.noalias() -= _localRhs;
            _localJac.block(0, u_idx, local_x.size(), 1) = -(local_r_1 - local_r_2) / 2 / eps / (1 + std::abs(local_x_[u_idx]));
  //          debugging--------------------------
  //          std::cout << "_localJac: \n";
  //          std::cout << _localJac << std::endl;
  //          end of debugging-------------------
  //          _localJac=
        }
    }

    void addJacobianToGlobal(
        AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& Jac) const override
    {
        Jac.add(indices, _localJac);
    }
	void preTimestep(std::vector<double> const& /*local_x*/) override
	{
		DBUG("Richards: preTimestep");
	}
	void postTimestep(std::vector<double> const& local_x,
		double const t) override
	{
		local_x0 = local_x;
		DBUG("Richards: postTimestep");
	}
    void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
        const override
    {
        M.add(indices, _localM);
        K.add(indices, _localA);
        b.add(indices.rows, _localRhs);
    }
    
    Eigen::Map<const Eigen::VectorXd>
        getShapeMatrix(const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::VectorXd>(N.data(), N.size());
    }
    

    
    std::vector<double> const&
        getIntegrationPointValues(IntegrationPointValue const property,
        std::vector<double>& /*cache*/) const override
    {
        switch (property)
        {
        case IntegrationPointValue::Saturation:
            return _saturation[0];
    
        }

        std::abort();
    }

private:
    MeshLib::Element const& _element;
    unsigned dim;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>> _shape_matrices;
    RichardsFlowProcessData const& _process_data;

    NodalMatrixType _localA;
    NodalMatrixType _localM;
    NodalMatrixType _localJac;
    NodalVectorType _localRhs;

	std::vector<double> local_x0;

    unsigned const _integration_order;
    //std::vector<double> _saturation=std::vector<double>()
    std::vector<std::vector<double>> _saturation
        = std::vector<std::vector<double>>(1,
        std::vector<double>(ShapeFunction::NPOINTS));
    public:
#ifdef OGS_USE_EIGEN
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#endif
};


}   // namespace  RichardsFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_RICHARDSFLOW_FEM_H_
