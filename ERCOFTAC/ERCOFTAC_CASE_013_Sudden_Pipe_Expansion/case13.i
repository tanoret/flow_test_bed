##########################################################
# ERCOFTAC test case
# Case Number: 013
# Author: Emmanuel Gault
# Last Update: May, 2024
# Turbulent model using:
# k-epsilon model
# Equilibrium + Newton wall treatement
# SIMPLE solve
##########################################################


Re = 400000

rho = 998
D = 0.1524
d = 0.0783
mu = 0.00100
bulk_u =  '${fparse Re * mu / rho / (3 * d)}'


pressure_tag = "pressure_grad"

### k-epslilon Closure Parameters ###
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09

### Initial and Boundary Conditions ###
intensity = 0.0319
k_init = '${fparse 1.5*(intensity * bulk_u/2)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / d}'

### Modeling parameters ###
non_equilibrium_treatment = true
bulk_wall_treatment = false
walls = 'wall-side_to_3 top wall-side_to_2'
walls_eps = 'wall-side_to_3'
noeps_walls = 'wall-side_to_2 top'
max_mixing_length = 1e10
linearized_yplus_mu_t = false
wall_treatment = 'eq_incremental' # Options: eq_newton, eq_incremental, eq_linearized, neq


[Functions]
    [u_fonc]
      type = PiecewiseLinear
      data_file = data_u.csv
      scale_factor = '${fparse bulk_u}'
      format = columns
    []
    [v_fonc]
        type = PiecewiseLinear
        data_file = data_v.csv
        scale_factor = '${fparse bulk_u}'
        format = columns
    []
[]

[Mesh]
    [pipe]
        type = CartesianMeshGenerator
        dim = 2
        dx = '${fparse 6.0*D} ${fparse 3*D}'
        dy = '${fparse 0.5*d} ${fparse 0.5*(D-d)}'
        ix = '50 25'
        iy = '4  4'
        subdomain_id = '
                        1 1
                        2 1
                    '
    []
    [corner_walls]
        type = SideSetsBetweenSubdomainsGenerator
        input = pipe
        primary_block = '1'
        paired_block = '2'
        new_boundary = 'wall-side'
    []
    [delete_bottom]
        type = BlockDeletionGenerator
        input = corner_walls
        block = '2'
    []
    [subdomain1]
        input = delete_bottom
        type = SubdomainBoundingBoxGenerator
        bottom_left = '0 ${fparse 0.5*d+0.0001} 0'
        top_right = '1 1 0'
        block_id = 2
    []
    [break_boundary]
        input = subdomain1
        type = BreakBoundaryOnSubdomainGenerator
        boundaries = 'wall-side'
    []
    [subdomain2]
        input = break_boundary
        type = SubdomainBoundingBoxGenerator
        bottom_left = '0 0 0'
        top_right = '1 ${fparse 0.5*d+0.0001} 0'
        block_id = 3
    []
    [bb2]
        input = subdomain2
        type = BreakBoundaryOnSubdomainGenerator
        boundaries = 'wall-side'
    []
    coord_type = RZ
    rz_coord_axis = X
    uniform_refine = 2
[]

[Problem]
    nl_sys_names = 'u_system v_system pressure_system TKE_system TKED_system'
    previous_nl_solution_required = true
[]

[GlobalParams]
    rhie_chow_user_object = 'rc'
    velocity_interp_method = 'rc'
[]

[UserObjects]
    [rc]
        type = INSFVRhieChowInterpolatorSegregated
        u = vel_x
        v = vel_r
        pressure = pressure
    []
[]

[Variables]
    [vel_x]
        type = INSFVVelocityVariable
        solver_sys = u_system
        two_term_boundary_expansion = false
        initial_condition = 1e-6
    []
    [vel_r]
        type = INSFVVelocityVariable
        solver_sys = v_system
        two_term_boundary_expansion = false
        initial_condition = 1e-6
    []
    [pressure]
        type = INSFVPressureVariable
        solver_sys = pressure_system
        two_term_boundary_expansion = false
        initial_condition = 1e-6
    []
    [TKE]
        type = INSFVEnergyVariable
        solver_sys = TKE_system
        initial_condition = ${k_init}
    []
    [TKED]
        type = INSFVEnergyVariable
        solver_sys = TKED_system
        initial_condition = ${eps_init}
    []
[]

[FVKernels]
    [u_advection]
        type = INSFVMomentumAdvection
        variable = vel_x
        advected_interp_method = 'average'
        rho = ${rho}
        momentum_component = 'x'
    []
    [u_viscosity]
        type = INSFVMomentumDiffusion
        variable = vel_x
        mu = ${mu}
        momentum_component = 'x'
    []
    [u_viscosity_turbulent]
        type = INSFVMomentumDiffusion
        variable = vel_x
        mu = 'mu_t'
        momentum_component = 'x'
        complete_expansion = true
        u = vel_x
        v = vel_r
    []
    [u_pressure]
        type = INSFVMomentumPressure
        variable = vel_x
        momentum_component = 'x'
        pressure = pressure
        extra_vector_tags = ${pressure_tag}
    []

    [v_advection]
        type = INSFVMomentumAdvection
        variable = vel_r
        advected_interp_method = 'average'
        rho = ${rho}
        momentum_component = 'y'
    []
    [v_viscosity]
        type = INSFVMomentumDiffusion
        variable = vel_r
        mu = ${mu}
        momentum_component = 'y'
    []
    [v_viscosity_turbulent]
        type = INSFVMomentumDiffusion
        variable = vel_r
        mu = 'mu_t'
        momentum_component = 'y'
        complete_expansion = true
        u = vel_x
        v = vel_r
    []

    [v_pressure]
        type = INSFVMomentumPressure
        variable = vel_r
        momentum_component = 'y'
        pressure = pressure
        extra_vector_tags = ${pressure_tag}
    []

    [p_diffusion]
        type = FVAnisotropicDiffusion
        variable = pressure
        coeff = "Ainv"
        coeff_interp_method = 'average'
    []
    [p_source]
        type = FVDivergence
        variable = pressure
        vector_field = "HbyA"
        force_boundary_execution = true
    []

    [TKE_advection]
        type = INSFVTurbulentAdvection
        variable = TKE
        advected_interp_method = 'upwind'
        rho = ${rho}
    []
    [TKE_diffusion]
        type = INSFVTurbulentDiffusion
        variable = TKE
        coeff = ${mu}
    []
    [TKE_diffusion_turbulent]
        type = INSFVTurbulentDiffusion
        variable = TKE
        coeff = 'mu_t'
        scaling_coef = ${sigma_k}
    []
    [TKE_source_sink]
        type = INSFVTKESourceSink
        variable = TKE
        u = vel_x
        v = vel_r
        epsilon = TKED
        rho = ${rho}
        mu = ${mu}
        mu_t = 'mu_t'
        walls = ${walls}
        non_equilibrium_treatment = ${non_equilibrium_treatment}
        max_mixing_length = ${max_mixing_length}
    []

    [TKED_advection]
        type = INSFVTurbulentAdvection
        variable = TKED
        advected_interp_method = 'upwind'
        rho = ${rho}
        walls = ${walls_eps}
    []
    [TKED_diffusion]
        type = INSFVTurbulentDiffusion
        variable = TKED
        coeff = ${mu}
        walls = ${walls_eps}
    []
    [TKED_diffusion_turbulent]
        type = INSFVTurbulentDiffusion
        variable = TKED
        coeff = 'mu_t'
        scaling_coef = ${sigma_eps}
        walls = ${walls_eps}
    []
    [TKED_source_sink]
        type = INSFVTKEDSourceSink
        variable = TKED
        u = vel_x
        v = vel_r
        k = TKE
        rho = ${rho}
        mu = ${mu}
        mu_t = 'mu_t'
        C1_eps = ${C1_eps}
        C2_eps = ${C2_eps}
        walls = ${noeps_walls}
        non_equilibrium_treatment = ${non_equilibrium_treatment}
        max_mixing_length = ${max_mixing_length}
    []
[]

[FVBCs]
    [inlet-u]
        type = INSFVInletVelocityBC
        boundary = 'left'
        variable = vel_x
        functor = u_fonc
    []
    [inlet-v]
        type = INSFVInletVelocityBC
        boundary = 'left'
        variable = vel_r
        functor = v_fonc
    []

    [inlet_TKE]
        type = INSFVInletIntensityTKEBC
        boundary = 'left'
        variable = TKE
        u = vel_x
        v = vel_r
        intensity = ${intensity}
    []
    [inlet_TKED]
        type = INSFVMixingLengthTKEDBC
        boundary = 'left'
        variable = TKED
        k = TKE
        characteristic_length = '${fparse d}'
    []
    [outlet_p]
        type = INSFVOutletPressureBC
        boundary = 'right'
        variable = pressure
        functor = 0
    []
    [walls-u]
        type = FVDirichletBC
        boundary = ${walls}
        variable = vel_x
        value = 0
    []
    [walls-v]
        type = FVDirichletBC
        boundary = ${walls}
        variable = vel_r
        value = 0
    []

    [walls_mu_t]
        type = INSFVTurbulentViscosityWallFunction
        boundary = ${walls}
        variable = mu_t
        u = vel_x
        v = vel_r
        rho = ${rho}
        mu = ${mu}
        mu_t = 'mu_t'
        k = TKE
        wall_treatment = ${wall_treatment}
    []
    [axis-u]
        type = INSFVSymmetryVelocityBC
        boundary = 'bottom'
        variable = vel_x
        u = vel_x
        v = vel_r
        mu = ${mu}
        momentum_component = x
    []
    [axis-v]
        type = INSFVSymmetryVelocityBC
        boundary = 'bottom'
        variable = vel_r
        u = vel_x
        v = vel_r
        mu = ${mu}
        momentum_component = y
    []
    [axis-p]
        type = INSFVSymmetryPressureBC
        boundary = 'bottom'
        variable = pressure
    []
    [axis-k]
        type = INSFVSymmetryScalarBC
        boundary = 'bottom'
        variable = TKE
    []
    [axis-eps]
        type = INSFVSymmetryScalarBC
        boundary = 'bottom'
        variable = TKED
    []
[]

[AuxVariables]
    [mu_t]
        type = MooseVariableFVReal
        initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
        two_term_boundary_expansion = false
    []
[]

[AuxKernels]
    [compute_mu_t]
        type = kEpsilonViscosityAux
        variable = mu_t
        C_mu = ${C_mu}
        k = TKE
        epsilon = TKED
        mu = ${mu}
        rho = ${rho}
        u = vel_x
        v = vel_r
        bulk_wall_treatment = ${bulk_wall_treatment}
        walls = ${walls}
        linearized_yplus = ${linearized_yplus_mu_t}
        non_equilibrium_treatment = ${non_equilibrium_treatment}
        execute_on = 'NONLINEAR'
    []
[]

[Postprocessors]
    [inlet_pressure]
        type = SideAverageValue
        boundary = 'left'
        variable = 'pressure'
        execute_on = 'FINAL'
        outputs = 'csv'
    []
    [average_TKE]
        type = ElementAverageValue
        variable = 'TKE'
        execute_on = 'FINAL'
        outputs = 'csv'
    []
    [average_TKED]
        type = ElementAverageValue
        variable = 'TKED'
        execute_on = 'FINAL'
        outputs = 'csv'
    []
    [average_outlet_velocity]
        type = SideAverageValue
        boundary = 'right'
        variable = 'vel_x'
        execute_on = 'FINAL'
        outputs = 'csv'
    []
[]



[Executioner]
    type = SIMPLE
    rhie_chow_user_object = 'rc'
    momentum_systems = 'u_system v_system'
    pressure_system = 'pressure_system'
    turbulence_systems = 'TKED_system TKE_system'

    pressure_gradient_tag = ${pressure_tag}
    momentum_equation_relaxation = 0.7
    pressure_variable_relaxation = 0.3
    turbulence_equation_relaxation = '0.3 0.3'
    num_iterations = 2000
    pressure_absolute_tolerance = 1e-12
    momentum_absolute_tolerance = 1e-12
    turbulence_absolute_tolerance = '1e-12 1e-12'
    momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
    momentum_petsc_options_value = 'hypre boomeramg'
    pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
    pressure_petsc_options_value = 'hypre boomeramg'

    momentum_l_abs_tol = 1e-14
    pressure_l_abs_tol = 1e-14
    turbulence_l_abs_tol = 1e-14
    momentum_l_max_its = 30
    pressure_l_max_its = 30
    momentum_l_tol = 0.0
    pressure_l_tol = 0.0
    turbulence_l_tol = 0.0
    print_fields = false
[]

[Outputs]
    exodus = true
    [csv]
      type = CSV
      execute_on = FINAL
    []
[]
