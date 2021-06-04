// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class
// in Geometry Central. Because we are "inside" the class, we no longer have to
// call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using
// geometry->buildHodgeStar0Form(), etc. where "geometry" is a pointer to a
// VertexPositionGeometry. This avoids having to declare a GeometryRoutines
// object in every project, and also mimics the way that geometry routines are
// normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be
 * applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
  // Diagonal Hodge Star. |V|*|V|
  auto vertexCount = mesh.nVertices();
  std::vector<Eigen::Triplet<double>> tripletList;
  for (const auto &vertex : mesh.vertices()) {
    tripletList.emplace_back(vertex.getIndex(), vertex.getIndex(),
                             barycentricDualArea(vertex));
  }
  auto star0 = geometrycentral::SparseMatrix<double>(vertexCount, vertexCount);
  star0.setFromTriplets(tripletList.begin(), tripletList.end());
  return star0;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be
 * applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
  // |E|*|E|
  auto edgeCount = mesh.nEdges();
  std::vector<Eigen::Triplet<double>> tripletList;
  for (const auto &edge : mesh.edges()) {
    double cotanWeight =
        0.5 * (cotan(edge.halfedge()) + cotan(edge.halfedge().twin()));
    tripletList.emplace_back(edge.getIndex(), edge.getIndex(), cotanWeight);
  }
  auto star1 = geometrycentral::SparseMatrix<double>(edgeCount, edgeCount);
  star1.setFromTriplets(tripletList.begin(), tripletList.end());
  return star1;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be
 * applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
  // |F|*|F|
  auto faceCount = mesh.nFaces();
  std::vector<Eigen::Triplet<double>> tripletList;
  for (const auto &face : mesh.faces()) {
    tripletList.emplace_back(face.getIndex(), face.getIndex(),
                             1.0 / faceArea(face));
  }
  auto star2 = geometrycentral::SparseMatrix<double>(faceCount, faceCount);
  star2.setFromTriplets(tripletList.begin(), tripletList.end());
  return star2;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be
 * applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form()
    const {
  // |E|*|V| * |V|*1 = |E|*1
  auto edgeCount = mesh.nEdges();
  auto vertexCount = mesh.nVertices();
  auto d0 = geometrycentral::SparseMatrix<double>(edgeCount, vertexCount);
  std::vector<Eigen::Triplet<double>> tripletList;
  for (const auto &edge : mesh.edges()) {
    auto firstVertexIndex = edge.firstVertex().getIndex();
    auto secondVertexIndex = edge.secondVertex().getIndex();
    tripletList.emplace_back(edge.getIndex(), firstVertexIndex, 1);
    tripletList.emplace_back(edge.getIndex(), secondVertexIndex, -1);
  }
  d0.setFromTriplets(tripletList.begin(), tripletList.end());
  return d0;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be
 * applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form()
    const {
  // |F|*|E| * |E|*1 = |F|*1
  auto edgeCount = mesh.nEdges();
  auto faceCount = mesh.nFaces();
  auto d1 = geometrycentral::SparseMatrix<double>(faceCount, edgeCount);
  std::vector<Eigen::Triplet<double>> tripletList;
  for (const auto &face : mesh.faces()) {
    for (const auto &halfedge : face.adjacentHalfedges()) {
      auto edge = halfedge.edge();
      if (edge.halfedge() == halfedge) {
        tripletList.emplace_back(face.getIndex(), edge.getIndex(), 1);
      } else {
        tripletList.emplace_back(face.getIndex(), edge.getIndex(), -1);
      }
    }
  }
  d1.setFromTriplets(tripletList.begin(), tripletList.end());
  return d1;
}

}  // namespace surface
}  // namespace geometrycentral