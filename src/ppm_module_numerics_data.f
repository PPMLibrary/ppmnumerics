      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_numerics_data
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_numerics_data
      !!! Declares global data types and variables

      !-------------------------------------------------------------------------
      !  Define quadrature rules for the unit triangle (BEM)
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_center = 1
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_nodes  = 2
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_edges  = 3
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_cne    = 4
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_stroud = 5
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_hammer3= 6
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_hammer4= 7
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_hammer7= 8
      INTEGER, PARAMETER :: ppm_param_bem_quadrule_hammer12= 9

      !-------------------------------------------------------------------------
      !  Define basis functions for the unit triangle (BEM)
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_bem_basis_const     = 1
      INTEGER, PARAMETER :: ppm_param_bem_basis_linear    = 2
      INTEGER, PARAMETER :: ppm_param_bem_basis_quad      = 3

      !-------------------------------------------------------------------------
      !  Define field solvers
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_solver_mg           = 0
      INTEGER, PARAMETER :: ppm_param_solver_sor          = 1
      INTEGER, PARAMETER :: ppm_param_solver_fft          = 2

      !------------------------------------------------------------------------
      !  Define equation ppm solvers
      !------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_eq_poisson          = 1

      !------------------------------------------------------------------------
      !  Define order for finite difference
      !------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_order_1             = 1
      INTEGER, PARAMETER :: ppm_param_order_2             = 2
      INTEGER, PARAMETER :: ppm_param_order_3             = 3
      INTEGER, PARAMETER :: ppm_param_order_4             = 4

      !------------------------------------------------------------------------
      !  Define smoother
      !------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_smooth_rbsor        = 1

      !-------------------------------------------------------------------------
      !  ODE parameters
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER :: ppm_param_ode_scheme_eulerf   = 1
      INTEGER, PARAMETER :: ppm_param_ode_scheme_tvdrk2   = 2
      INTEGER, PARAMETER :: ppm_param_ode_scheme_midrk2   = 3
      INTEGER, PARAMETER :: ppm_param_ode_scheme_rk4      = 4
      !INTEGER, PARAMETER :: ppm_param_ode_scheme_trapez   = 5
      !INTEGER, PARAMETER :: ppm_param_ode_scheme_tvdrk3   = 6
      INTEGER, PARAMETER :: ppm_param_ode_scheme_sts      = 7

      ! nothing here for now

      END MODULE ppm_module_numerics_data
