##########################################################
# ERCOFTAC test case
# Case Number: 060
# Author: Emmanuel Gault
# Last Update: June, 2024
# Turbulent model using:
# k-epsilon model
# Equilibrium + Newton wall treatement
# SIMPLE solve
##########################################################


#Re = 202000

rho = 1.225
mu = 1.79e-5
bulk_u = 11.6
D = 0.26


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
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / D}'

### Modeling parameters ###
bulk_wall_treatment = false
wall_treatment = 'eq_newton' # Options: eq_newton, eq_incremental, eq_linearized, neq




[Mesh]
    [fmg]
        type = FileMeshGenerator
        file = 'mesh_ERCOFTAC_case_60_normal_with_bcs1.e'
    []
    [subdomain1]
        input = fmg
        type = SubdomainBoundingBoxGenerator
        bottom_left = '0 -1 -1'
        top_right = '0.1 1 1'
        block_id = 2
    []
    [break_boundary]
        input = subdomain1
        type = BreakBoundaryOnSubdomainGenerator
        boundaries = 'wall'
    []
    [subdomain2]
        input = break_boundary
        type = SubdomainBoundingBoxGenerator
        bottom_left = '0.1 -1 -1'
        top_right = '1 1 1'
        block_id = 3
    []
    [bb2]
        input = subdomain2
        type = BreakBoundaryOnSubdomainGenerator
        boundaries = 'wall'
    []
    uniform_refine = 0
[]

[Functions]
    [v_fonc]
      type = PiecewiseLinear
      data_file = testswirlv.csv
      scale_factor = 53
      format = columns
      axis = z
    []
    [w_fonc]
      type = PiecewiseLinear
      data_file = testswirlw.csv
      scale_factor = 53
      format = columns
      axis = y
    []
    [R_fonc]
        type = PiecewiseLinear
        data_file = curvR.csv
        scale_factor = 1
        format = columns
        axis = x
    []
[]

[Problem]
    nl_sys_names = 'u_system v_system w_system pressure_system TKE_system TKED_system'
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
        v = vel_y
        w = vel_z
        pressure = pressure
    []
[]

[Variables]
    [vel_x]
        type = INSFVVelocityVariable
        solver_sys = u_system
        two_term_boundary_expansion = false
        initial_condition = 2
    []
    [vel_y]
        type = INSFVVelocityVariable
        solver_sys = v_system
        two_term_boundary_expansion = false
        initial_condition = 0
    []
    [vel_z]
        type = INSFVVelocityVariable
        solver_sys = w_system
        two_term_boundary_expansion = false
        initial_condition = 0
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
        advected_interp_method = 'upwind'
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
        v = vel_y
        w = vel_z
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
        variable = vel_y
        advected_interp_method = 'upwind'
        rho = ${rho}
        momentum_component = 'y'
    []
    [v_viscosity]
        type = INSFVMomentumDiffusion
        variable = vel_y
        mu = ${mu}
        momentum_component = 'y'
    []
    [v_viscosity_turbulent]
        type = INSFVMomentumDiffusion
        variable = vel_y
        mu = 'mu_t'
        momentum_component = 'y'
        complete_expansion = true
        u = vel_x
        v = vel_y
        w = vel_z
    []
    [v_pressure]
        type = INSFVMomentumPressure
        variable = vel_y
        momentum_component = 'y'
        pressure = pressure
        extra_vector_tags = ${pressure_tag}
    []

    [w_advection]
        type = INSFVMomentumAdvection
        variable = vel_z
        advected_interp_method = 'upwind'
        rho = ${rho}
        momentum_component = 'z'
    []
    [w_viscosity]
        type = INSFVMomentumDiffusion
        variable = vel_z
        mu = ${mu}
        momentum_component = 'z'
    []
    [w_viscosity_turbulent]
        type = INSFVMomentumDiffusion
        variable = vel_z
        mu = 'mu_t'
        momentum_component = 'z'
        complete_expansion = true
        u = vel_x
        v = vel_y
        w = vel_z
    []
    [w_pressure]
        type = INSFVMomentumPressure
        variable = vel_z
        momentum_component = 'z'
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
        v = vel_y
        w = vel_z
        epsilon = TKED
        rho = ${rho}
        mu = ${mu}
        mu_t = 'mu_t'
        walls = 'wall'
        wall_treatment = ${wall_treatment}
    []

    [TKED_advection]
        type = INSFVTurbulentAdvection
        variable = TKED
        advected_interp_method = 'upwind'
        rho = ${rho}
        walls = 'wall_to_2'
    []
    [TKED_diffusion]
        type = INSFVTurbulentDiffusion
        variable = TKED
        coeff = ${mu}
        walls = 'wall_to_2'
    []
    [TKED_diffusion_turbulent]
        type = INSFVTurbulentDiffusion
        variable = TKED
        coeff = 'mu_t'
        scaling_coef = ${sigma_eps}
        walls = 'wall_to_2'
    []
    [TKED_source_sink]
        type = INSFVTKEDSourceSink
        variable = TKED
        u = vel_x
        v = vel_y
        w = vel_z
        k = TKE
        rho = ${rho}
        mu = ${mu}
        mu_t = 'mu_t'
        C1_eps = ${C1_eps}
        C2_eps = ${C2_eps}
        walls = 'wall' 
        wall_treatment = ${wall_treatment}
    []
[]

[FVBCs]
    [inlet-u]
        type = INSFVInletVelocityBC
        boundary = 'inlet'
        variable = vel_x
        functor = '${fparse bulk_u}'
    []
    [inlet-v]
        type = INSFVInletVelocityBC
        boundary = 'inlet'
        variable = vel_y
        functor = v_fonc
    []
    [inlet-w]
        type = INSFVInletVelocityBC
        boundary = 'inlet'
        variable = vel_z
        functor = w_fonc
    []

    [inlet_TKE]
        type = INSFVInletIntensityTKEBC
        boundary = 'inlet'
        variable = TKE
        u = vel_x
        v = vel_y
        w = vel_z
        intensity = ${intensity}
    []
    [inlet_TKED]
        type = INSFVMixingLengthTKEDBC
        boundary = 'inlet'
        variable = TKED
        k = TKE
        characteristic_length = '${fparse D}'
    []
    [outlet_p]
        type = INSFVOutletPressureBC
        boundary = 'outlet'
        variable = pressure
        functor = 0
    []
    [walls-u]
        type = FVDirichletBC
        boundary = 'wall'
        variable = vel_x
        value = 0
    []
    [walls-v]
        type = FVDirichletBC
        boundary = 'wall'
        variable = vel_y
        value = 0
    []
    [walls-w]
        type = FVDirichletBC
        boundary = 'wall'
        variable = vel_z
        value = 0
    []

    [walls_mu_t]
        type = INSFVTurbulentViscosityWallFunction
        boundary = 'wall'
        variable = mu_t
        u = vel_x
        v = vel_y
        w = vel_z
        rho = ${rho}
        mu = ${mu}
        mu_t = 'mu_t'
        k = TKE
        wall_treatment = ${wall_treatment}
        curv_R = R_fonc
        x_curvature_axis = 1
        y_curvature_axis = 0
        z_curvature_axis = 0
        alpha_curv = 2
    []
[]

[AuxVariables]
    [mu_t]
        type = MooseVariableFVReal
        initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
        two_term_boundary_expansion = false
    []
    [mu_eff]
        type = MooseVariableFVReal
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
        v = vel_y
        bulk_wall_treatment = ${bulk_wall_treatment}
        walls = 'wall'
        wall_treatment = ${wall_treatment}
        execute_on = 'NONLINEAR'
    []
    [compute_mu_eff]
        type = ParsedAux
        variable = mu_eff
        coupled_variables = 'mu_t'
        expression = '${mu} + mu_t'
    []
[]



[Executioner]
    type = SIMPLENonlinearAssembly
    rhie_chow_user_object = 'rc'
    momentum_systems = 'u_system v_system w_system'
    pressure_system = 'pressure_system'
    turbulence_systems = 'TKED_system TKE_system'

    pressure_gradient_tag = ${pressure_tag}
    momentum_equation_relaxation = 0.6
    pressure_variable_relaxation = 0.4
    turbulence_equation_relaxation = '0.2 0.2'
    num_iterations = 800
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
[]
