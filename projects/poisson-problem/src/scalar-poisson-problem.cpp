// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"

#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ScalarPoissonProblem::ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh,
                                           VertexPositionGeometry* inputGeo) {
  mesh = inputMesh;
  geometry = inputGeo;
  this->A = geometry->laplaceMatrix();
  this->M = geometry->massMatrix();
  this->totalArea = geometry->totalArea();
}

/*
 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), where A
 * is the POSITIVE DEFINITE Laplace matrix and M is the mass matrix.
 *
 * Input: <rho>, the density of vertices in the mesh.
 * Returns: The solution vector.
 */
Vector<double> ScalarPoissonProblem::solve(const Vector<double>& rho) const {
  // Note: Geometry Central has linear solvers:
  // https://geometry-central.net/numerical/linear_solvers/
  // Page 101 of the Notes
  // rhoBar = (1/totalArea) * \int rho dV = (1/totalArea) * M * rho
  auto rhoBar = Vector<double>::Constant(
      rho.size(), (this->M * rho).sum() / this->totalArea);
  SparseMatrix<double> laplace{this->A};
  Vector<double> rhs = -this->M * (rho - rhoBar);
  // Solve Ax = -M(rho - rhoBar)
  return solvePositiveDefinite(laplace, rhs);
}