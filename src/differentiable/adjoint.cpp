#include <polyfem/solver/AdjointTools.hpp>
#include <polyfem/utils/MatrixUtils.hpp>
#include <polyfem/State.hpp>
#include "binding.hpp"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace polyfem;
using namespace polyfem::solver;

void define_adjoint(py::module_ &m)
{
  m.def(
      "shape_derivative",
      [](State &state) {
        Eigen::VectorXd term;
        if (state.problem->is_time_dependent())
          AdjointTools::dJ_shape_transient_adjoint_term(
              state, state.get_adjoint_mat(1), state.get_adjoint_mat(0), term);
        else
          AdjointTools::dJ_shape_static_adjoint_term(
              state, state.diff_cached.u(0), state.get_adjoint_mat(0), term);
        return utils::unflatten(term, state.mesh->dimension());
      },
      py::arg("solver"));

  m.def(
      "elastic_material_derivative",
      [](State &state) {
        Eigen::VectorXd term;
        if (state.problem->is_time_dependent())
          AdjointTools::dJ_material_transient_adjoint_term(
              state, state.get_adjoint_mat(1), state.get_adjoint_mat(0), term);
        else
          AdjointTools::dJ_material_static_adjoint_term(
              state, state.diff_cached.u(0), state.get_adjoint_mat(0), term);

        return utils::unflatten(term, state.bases.size());
      },
      py::arg("solver"));

  m.def(
      "friction_coefficient_derivative",
      [](State &state) {
        Eigen::VectorXd term;
        if (state.problem->is_time_dependent())
          AdjointTools::dJ_friction_transient_adjoint_term(
              state, state.get_adjoint_mat(1), state.get_adjoint_mat(0), term);
        else
          log_and_throw_adjoint_error(
              "Friction coefficient derivative is only supported for transient problems!");

        return term(0);
      },
      py::arg("solver"));

  m.def(
      "initial_velocity_derivative",
      [](State &state) {
        const int dim = state.mesh->dimension();

        Eigen::VectorXd term;
        if (!state.problem->is_time_dependent())
          log_and_throw_adjoint_error(
              "Initial condition derivative is only supported for transient problems!");

        AdjointTools::dJ_initial_condition_adjoint_term(
            state, state.get_adjoint_mat(1), state.get_adjoint_mat(0), term);

        std::unordered_map<int, Eigen::VectorXd> map;
        for (int e = 0; e < state.mesh->n_elements(); e++)
        {
          const int id = state.mesh->get_body_id(e);
          if (map.find(id) == map.end())
            map[id] = Eigen::VectorXd::Zero(dim);
        }

        Eigen::Matrix<bool, 1, -1> visited(state.n_bases);
        visited.setConstant(false);
        for (int e = 0; e < state.mesh->n_elements(); e++)
        {
          const int bid = state.mesh->get_body_id(e);
          auto &vec = map[bid];
          for (const auto &b : state.bases[e].bases)
            for (const auto &g : b.global())
            {
              if (visited(g.index))
                continue;
              visited(g.index) = true;
              vec += term.segment(state.ndof() + g.index * dim, dim);
            }
        }

        return map;
      },
      py::arg("solver"));

  m.def(
      "initial_displacement_derivative",
      [](State &state) {
        const int dim = state.mesh->dimension();

        Eigen::VectorXd term;
        if (!state.problem->is_time_dependent())
          log_and_throw_adjoint_error(
              "Initial condition derivative is only supported for transient problems!");

        AdjointTools::dJ_initial_condition_adjoint_term(
            state, state.get_adjoint_mat(1), state.get_adjoint_mat(0), term);

        std::unordered_map<int, Eigen::VectorXd> map;
        for (int e = 0; e < state.mesh->n_elements(); e++)
        {
          const int id = state.mesh->get_body_id(e);
          if (map.find(id) == map.end())
            map[id] = Eigen::VectorXd::Zero(dim);
        }

        Eigen::Matrix<bool, 1, -1> visited(state.n_bases);
        visited.setConstant(false);
        for (int e = 0; e < state.mesh->n_elements(); e++)
        {
          const int bid = state.mesh->get_body_id(e);
          auto &vec = map[bid];
          for (const auto &b : state.bases[e].bases)
            for (const auto &g : b.global())
            {
              if (visited(g.index))
                continue;
              visited(g.index) = true;
              vec += term.segment(g.index * dim, dim);
            }
        }

        return map;
      },
      py::arg("solver"));

  m.def(
      "dirichlet_derivative",
      [](State &state) {
        const int dim = state.mesh->dimension();

        Eigen::VectorXd term;
        if (state.problem->is_time_dependent())
          log_and_throw_adjoint_error(
              "Dirichlet derivative is only supported for static problems!");

        AdjointTools::dJ_dirichlet_static_adjoint_term(
            state, state.get_adjoint_mat(0), term);
        return utils::unflatten(term, state.mesh->dimension());
      },
      py::arg("solver"));
}
