// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"

#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh,
                                     VertexPositionGeometry* inputGeo) {
  // Build member variables: mesh, geometry
  mesh = inputMesh;
  geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(
    const SparseMatrix<double>& M, double h) const {
  auto const vertexCount = mesh->nVertices();
  auto identity = identityMatrix<double>(vertexCount);
  // (M-h\Delta)f_h=Mf_0
  // In the code: (M+h\Delta)f_h=Mf_0 (our discrete Laplace matrix is the
  // negative of the actual Laplacian.)
  return M + h * geometry->laplaceMatrix();
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {
  auto const vertexCount = mesh->nVertices();
  auto mass = geometry->massMatrix();
  auto flowOperator = buildFlowOperator(mass, h);
  DenseMatrix<double> f{vertexCount, 3};
  for (const auto& vertex : mesh->vertices()) {
    f(vertex.getIndex(), 0) = geometry->inputVertexPositions[vertex].x;
    f(vertex.getIndex(), 1) = geometry->inputVertexPositions[vertex].y;
    f(vertex.getIndex(), 2) = geometry->inputVertexPositions[vertex].z;
  }
  // We need to solve Af_h=mf_0
  f = mass * f;
  PositiveDefiniteSolver<double> solver(flowOperator);
  f.col(0) = solver.solve(f.col(0));
  f.col(1) = solver.solve(f.col(1));
  f.col(2) = solver.solve(f.col(2));
  // Note: Update positions via geometry->inputVertexPositions
  for (const auto& vertex : mesh->vertices()) {
    geometry->inputVertexPositions[vertex] =
        Vector3{f(vertex.getIndex(), 0), f(vertex.getIndex(), 1),
                f(vertex.getIndex(), 2)};
  }
}