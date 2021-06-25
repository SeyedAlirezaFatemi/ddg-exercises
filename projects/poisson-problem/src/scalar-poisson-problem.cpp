// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"

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
  // TODO
  // Note: Geometry Central has linear solvers:
  // https://geometry-central.net/numerical/linear_solvers/
  return Vector<double>::Zero(rho.rows());  // placeholder
}