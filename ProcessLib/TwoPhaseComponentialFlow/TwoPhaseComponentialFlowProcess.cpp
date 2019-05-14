/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "TwoPhaseComponentialFlowProcess.h"

#include <cassert>
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "TwoPhaseComponentialFlowLocalAssembler.h"

namespace ProcessLib
{
    namespace TwoPhaseComponentialFlow
    {
        TwoPhaseComponentialFlowProcess::TwoPhaseComponentialFlowProcess(
            MeshLib::Mesh& mesh,
            std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
            std::vector<std::unique_ptr<ParameterBase>> const& parameters,
            unsigned const integration_order,
            std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
            TwoPhaseComponentialFlowProcessData&& process_data,
            SecondaryVariableCollection&& secondary_variables,
            NumLib::NamedFunctionCaller&& named_function_caller,
            BaseLib::ConfigTree const& /*config*/,
            std::map<std::string,
            std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        /*curves*/)
            : Process(mesh, std::move(jacobian_assembler), parameters,
                integration_order, std::move(process_variables),
                std::move(secondary_variables), std::move(named_function_caller)),
            _process_data(std::move(process_data))
        {
            DBUG("Create Two Phase Componential Flow Process.");
        }

        void TwoPhaseComponentialFlowProcess::initializeConcreteProcess(
            NumLib::LocalToGlobalIndexMap const& dof_table,
            MeshLib::Mesh const& mesh,
            unsigned const integration_order)
        {
            const int process_id = 0;
            ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
            ProcessLib::createLocalAssemblers<TwoPhaseComponentialFlowLocalAssembler>(
                mesh.getDimension(), mesh.getElements(), dof_table,
                pv.getShapeFunctionOrder(), _local_assemblers,
                mesh.isAxiallySymmetric(), integration_order, _process_data);

            _secondary_variables.addSecondaryVariable(
                "saturation",
                makeExtrapolator(1,
                    getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtSaturation));

            _secondary_variables.addSecondaryVariable(
                "pressure_wetting",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::
                    getIntPtWettingPressure));

            _secondary_variables.addSecondaryVariable(
                "mol_frac_nonwet_vapor",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::
                    getIntPtMolFracNonwetVapor));

            _secondary_variables.addSecondaryVariable(
                "mol_frac_nonwet_air",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtMolFracNonwetAir));

            _secondary_variables.addSecondaryVariable(
                "co2_concentration",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::
                    getIntPtCO2Concentration));

            _secondary_variables.addSecondaryVariable(
                "porosity",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtPorosityValue));

            _secondary_variables.addSecondaryVariable(
                "pH_value",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtpHValue));

            _secondary_variables.addSecondaryVariable(
                "co2_cumulated_prev",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtRhoMolCo2CumulTotalPrev));

            _secondary_variables.addSecondaryVariable(
                "mole_amount_organic_waste_cellulose",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtMolAmountWasteCell));

            _secondary_variables.addSecondaryVariable(
                "mole_amount_organic_waste_polystyrene",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtMolAmountWastePoly));


            _secondary_variables.addSecondaryVariable(
                "sio2_cumulated_prev",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtMolRhoSiO2CumulPrev));

            _secondary_variables.addSecondaryVariable(
                "moldensity_gasphase",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::
                    getIntPtMolDensityGasPhase));

            _secondary_variables.addSecondaryVariable(
                "moldensity_liquidphase",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::
                    getIntPtMolDensityLiquidPhase));

            _secondary_variables.addSecondaryVariable(
                "total_velocity_gas_x",
                makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtTotalVelocityGasX));

            if (mesh.getDimension() > 1)
            {
                _secondary_variables.addSecondaryVariable(
                    "total_velocity_gas_y",
                    makeExtrapolator(
                        1, getExtrapolator(), _local_assemblers,
                        &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtTotalVelocityGasY));
            }
            if (mesh.getDimension() > 2)
            {
                _secondary_variables.addSecondaryVariable(
                    "total_velocity_gas_z",
                    makeExtrapolator(
                        1, getExtrapolator(), _local_assemblers,
                        &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtTotalVelocityGasZ));
            }

            _secondary_variables.addSecondaryVariable(
                "total_velocity_liquid_x",
                makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtTotalVelocityLiquidX));

            if (mesh.getDimension() > 1)
            {
                _secondary_variables.addSecondaryVariable(
                    "total_velocity_y",
                    makeExtrapolator(
                        1, getExtrapolator(), _local_assemblers,
                        &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtTotalVelocityLiquidY));
            }
            if (mesh.getDimension() > 2)
            {
                _secondary_variables.addSecondaryVariable(
                    "total_velocity_z",
                    makeExtrapolator(
                        1, getExtrapolator(), _local_assemblers,
                        &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtTotalVelocityLiquidZ));
            }


            _secondary_variables.addSecondaryVariable(
                "darcy_volumetric_flux_total_gas_phase",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtOverallDarcyVolumetricFluxGas));

            _secondary_variables.addSecondaryVariable(
                "darcy_volumetric_flux_total_liquid_phase",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtOverallDarcyVolumetricFluxLiquid));

            _secondary_variables.addSecondaryVariable(
                "overall_gas_generation_rate",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasGenerationRate));

            _secondary_variables.addSecondaryVariable(
                "gas_hydrogen_generation_rate",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasHydrogenGenerateRate));

            _secondary_variables.addSecondaryVariable(
                "gas_hydrogen_generation_source_rate",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasHydrogenGenerateSourceRate));

            _secondary_variables.addSecondaryVariable(
                "gas_methane_generation_rate",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasMethaneGenerateRate));

            _secondary_variables.addSecondaryVariable(
                "gas_carbon_degradation_rate",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasCarbonDegradationRate));

            _secondary_variables.addSecondaryVariable(
                "gas_carbon_generation_rate",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasCarbonGenerateRate));

            _secondary_variables.addSecondaryVariable(
                "gas_co2_darcy_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasCO2DarcyVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_co2_diffusive_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasCO2DiffusiveVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_hydrogen_darcy_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasHydrogenDarcyVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_hydrogen_diffusive_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasHydrogenDiffusiveVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_methane_darcy_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasMethaneDarcyVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_methane_diffusive_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasMethaneDiffusiveVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "CO2_Consumed_For_Current_Step",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtCO2ConsumedcurrentStep));

            _secondary_variables.addSecondaryVariable(
                "gas_water_vapor_darcy_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtVaporDarcyVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_water_vapor_diffusive_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtVaporDiffusiveVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_nitrogen_darcy_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasNitrogenDarcyVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_nitrogen_diffusive_volumetric_flux",
                makeExtrapolator(mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasNitrogenDiffusiveVolumetricFlux));

            _secondary_variables.addSecondaryVariable(
                "gas_nitrogen_diffusive_volumetric_flux_chainrule",
                makeExtrapolator(
                    mesh.getDimension(), getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::
                        getIntPtGasNitrogenDiffusiveVolumetricFluxChainRule));

            _secondary_variables.addSecondaryVariable(
                "relative_humidity",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtRelHumidity));

            _secondary_variables.addSecondaryVariable(
                "bazant_power",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtReactivity));

            _secondary_variables.addSecondaryVariable(
                "water_consumption_rate",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtWaterConsumpRate));

            _secondary_variables.addSecondaryVariable(
                "gas_n2_variation_rate_porosity",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtGasN2VariationPorosity));

            _secondary_variables.addSecondaryVariable(
                "molar_fraction_water_in_liquid_phase",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::getIntPtmolfracwaterinliquid));
            //secondary variable on each cell
            auto mesh_prop_saturation = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "saturation_cell",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_saturation->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_saturation = mesh_prop_saturation;

            auto mesh_prop_porosity=MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "porosity_cell",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_porosity->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_porosity = mesh_prop_porosity;

            auto mesh_prop_pHvalue = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "pH_value_cell",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_pHvalue->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_pHvalue = mesh_prop_pHvalue;

            auto mesh_prop_bazant_power = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "bazant_power_cell",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_bazant_power->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_bazant_power = mesh_prop_bazant_power;

            auto mesh_prop_mol_density_gas = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "molar_density_gas_cell",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_mol_density_gas->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_mol_density_gas = mesh_prop_mol_density_gas;

            auto mesh_prop_mol_density_liquid = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "molar_density_liquid_cell",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_mol_density_liquid->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_mol_density_liquid = mesh_prop_mol_density_liquid;

            auto mesh_prop_co2_cumulate_consume
                = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "consume_co2_cumulate_cell",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_co2_cumulate_consume->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_co2_cumulate_consume = mesh_prop_co2_cumulate_consume;

            auto mesh_prop_sio2_cumulate_consume
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "consume_sio2_cumulate_cell",
                    MeshLib::MeshItemType::Cell, 1);
            mesh_prop_sio2_cumulate_consume->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_sio2_cumulate_consume = mesh_prop_sio2_cumulate_consume;

            auto mesh_prop_total_liquid_darcy_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "total_liquid_phase_darcy_volumetric_flux_cell",
                MeshLib::MeshItemType::Cell, 3);
            mesh_prop_total_liquid_darcy_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_overall_liquid_darcy_volumetric_flux
                = mesh_prop_total_liquid_darcy_volumetric_flux;

            auto mesh_prop_total_gas_darcy_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "total_gas_phase_darcy_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_total_gas_darcy_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_overall_gas_darcy_volumetric_flux
                = mesh_prop_total_gas_darcy_volumetric_flux;

            auto mesh_prop_gas_co2_darcy_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "co2_gas_darcy_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_gas_co2_darcy_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_gas_co2_darcy_volumetric_flux
                = mesh_prop_gas_co2_darcy_volumetric_flux;

            auto mesh_prop_gas_co2_diffusive_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "co2_gas_diffusive_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_gas_co2_diffusive_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_gas_co2_diffusive_volumetric_flux
                = mesh_prop_gas_co2_diffusive_volumetric_flux;


            auto mesh_prop_gas_h2_darcy_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "hydrogen_gas_darcy_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_gas_h2_darcy_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_gas_hydrogen_darcy_volumetric_flux
                = mesh_prop_gas_h2_darcy_volumetric_flux;

            auto mesh_prop_gas_h2_diffusive_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "hydrogen_gas_diffusive_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_gas_h2_diffusive_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_gas_hydrogen_diffusive_volumetric_flux
                = mesh_prop_gas_h2_diffusive_volumetric_flux;

            auto mesh_prop_gas_ch4_darcy_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "methane_gas_darcy_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_gas_ch4_darcy_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_gas_methane_darcy_volumetric_flux
                = mesh_prop_gas_ch4_darcy_volumetric_flux;

            auto mesh_prop_gas_ch4_diffusive_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "methane_gas_diffusive_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_gas_ch4_diffusive_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_gas_methane_diffusive_volumetric_flux
                = mesh_prop_gas_ch4_diffusive_volumetric_flux;

            auto mesh_prop_gas_h2o_vapor_darcy_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "water_vapor_gas_darcy_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_gas_h2o_vapor_darcy_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_gas_water_vapor_darcy_volumetric_flux
                = mesh_prop_gas_h2o_vapor_darcy_volumetric_flux;

            auto mesh_prop_gas_h2o_vapor_diffusive_volumetric_flux
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "water_vapor_gas_diffusive_volumetric_flux_cell",
                    MeshLib::MeshItemType::Cell, 3);
            mesh_prop_gas_h2o_vapor_diffusive_volumetric_flux
                ->resize(mesh.getNumberOfElements() * 3);
            _process_data.mesh_prop_gas_water_vapor_diffusive_volumetric_flux
                = mesh_prop_gas_h2o_vapor_diffusive_volumetric_flux;

            auto mesh_prop_mol_fraction_h2o_vaor
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "molar_fraction_h2o_vapor_cell",
                    MeshLib::MeshItemType::Cell, 1);
            mesh_prop_mol_fraction_h2o_vaor->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_mol_frac_h2o_vapor = mesh_prop_mol_fraction_h2o_vaor;

            auto mesh_prop_mol_fraction_n2
                = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "molar_fraction_n2_cell",
                    MeshLib::MeshItemType::Cell, 1);
            mesh_prop_mol_fraction_n2->resize(mesh.getNumberOfElements() * 1);
            _process_data.mesh_prop_mol_frac_n2 = mesh_prop_mol_fraction_n2;


        }

        void TwoPhaseComponentialFlowProcess::assembleConcreteProcess(const double t,
            GlobalVector const& x,
            GlobalMatrix& M,
            GlobalMatrix& K,
            GlobalVector& b)
        {
            DBUG("Assemble TwoPhaseFlowWithPrhoProcess.");
            std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
                dof_table = { std::ref(*_local_to_global_index_map) };
            // Call global assembler for each local assembly item.
            GlobalExecutor::executeMemberDereferenced(
                _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
                dof_table, t, x, M, K, b, _coupled_solutions);
        }

        void TwoPhaseComponentialFlowProcess::assembleWithJacobianConcreteProcess(
            const double t, GlobalVector const& x, GlobalVector const& xdot,
            const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
            GlobalVector& b, GlobalMatrix& Jac)
        {
            DBUG("AssembleWithJacobian TwoPhaseFlowWithPrhoProcess.");


            std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
                dof_table = { std::ref(*_local_to_global_index_map) };
            // Call global assembler for each local assembly item.
            GlobalExecutor::executeMemberDereferenced(
                _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
                _local_assemblers, dof_table, t, x, xdot, dxdot_dx,
                dx_dx, M, K, b, Jac, _coupled_solutions);
        }

    }  // end of namespace
}  // end of namespace
