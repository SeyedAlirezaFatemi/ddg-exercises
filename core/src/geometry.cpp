// PLEASE READ:
//
// This file implements additional geometry routines for the
// VertexPositionGeometry class in Geometry Central. Because we are "inside" the
// class, we no longer have to call
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
// geometry->cotan(he), geometry->barycentricDualArea(v), etc. where "geometry"
// is a pointer to a VertexPositionGeometry. This avoids having to declare a
// GeometryRoutines object in every project, and also mimics the way that
// geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include <complex>

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
  return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {
  double total = 0.0;
  for (Edge e : mesh.edges()) {
    total += edgeLength(e);
  }
  return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {
  double total = 0.0;
  for (Face f : mesh.faces()) {
    total += faceArea(f);
  }
  return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use
 * built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {
  if (!he.isInterior()) {
    return 0.0;
  }
  // cot(a, b) = (a * b) / |a x b|
  // Get the outgoing halfedges
  auto firstHalfedgeVector = halfedgeVector(he.next().twin());
  auto secondHalfedgeVector = halfedgeVector(he.next().next());
  return dot(firstHalfedgeVector, secondHalfedgeVector) /
         cross(firstHalfedgeVector, secondHalfedgeVector).norm();
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {
  // The barycentric dual area associated with a vertex i is equal to one-third
  // the area of all triangles ijk touching i.
  double totalArea = 0.0;
  for (const auto& face : v.adjacentFaces()) {
    totalArea += faceArea(face);
  }
  return totalArea / 3;
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in
 * function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {
  // Get the two outgoing vectors from corner c
  // c.halfedge() returns the halfedge whose tail touches this corner
  auto firstVector = halfedgeVector(c.halfedge());
  auto secondVector = halfedgeVector(c.halfedge().next().next().twin());
  return acos(dot(firstVector, secondVector) /
              (firstVector.norm() * secondVector.norm()));
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT
 * use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral
 * angle is computed. Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {
  auto firstNormal = faceNormal(he.face());
  auto secondNormal = faceNormal(he.twin().face());
  auto edgeVector = halfedgeVector(he);
  return atan2(
      dot(edgeVector, cross(firstNormal, secondNormal)) / edgeVector.norm(),
      dot(firstNormal, secondNormal));
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {
  Vector3 normal{0., 0., 0.};
  for (const auto& face : v.adjacentFaces()) {
    normal += faceNormal(face);
  }
  return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {
  Vector3 normal{0., 0., 0.};
  for (const auto& corner : v.adjacentCorners()) {
    normal += angle(corner) * faceNormal(corner.face());
  }
  return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {
  Vector3 normal{0., 0., 0.};
  for (const auto& corner : v.adjacentCorners()) {
    auto firstEdge = halfedgeVector(corner.halfedge());
    auto secondEdge = halfedgeVector(corner.halfedge().next().next().twin());
    normal +=
        cross(firstEdge, secondEdge) / (firstEdge.norm2() * secondEdge.norm2());
  }
  return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {
  Vector3 normal{0., 0., 0.};
  for (const auto& face : v.adjacentFaces()) {
    normal += faceArea(face) * faceNormal(face);
  }
  return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {
  Vector3 normal{0., 0., 0.};
  for (const auto& halfedge : v.outgoingHalfedges()) {
    auto edgeVector = halfedgeVector(halfedge);
    normal += dihedralAngle(halfedge) * edgeVector / edgeVector.norm();
  }
  return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent
 * to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {
  Vector3 normal{0., 0., 0.};
  for (const auto& halfedge : v.outgoingHalfedges()) {
    normal += edgeCotanWeight(halfedge.edge()) * halfedgeVector(halfedge);
  }
  return normal.normalize();
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {
  auto angleSum{0.0};
  for (const auto& corner : v.adjacentCorners()) {
    angleSum += cornerAngle(corner);
  }
  return 2 * PI - angleSum;
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {
  auto totalDefect{0.0};
  for (const auto& vertex : mesh.vertices()) {
    totalDefect += angleDefect(vertex);
  }
  return totalDefect;
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {
  // Scalar mean curvature = 0.5 * \sum_{ij \in E} l_{ij}\theta_{ij}
  auto meanCurvature{0.0};
  for (const auto& edge : v.adjacentEdges()) {
    meanCurvature += edgeDihedralAngle(edge) * edgeLength(edge);
  }
  return 0.5 * meanCurvature;
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {
  auto dualArea{0.0};
  for (const auto& corner : v.adjacentCorners()) {
    auto firstHalfedge = corner.halfedge();
    auto secondHalfedge = corner.halfedge().next().next();
    dualArea += halfedgeVector(firstHalfedge).norm2() * cotan(firstHalfedge) +
                halfedgeVector(secondHalfedge).norm2() * cotan(secondHalfedge);
  }
  return 0.125 * dualArea;
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a
 * vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature
 * values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(
    Vertex v) const {
  auto area = circumcentricDualArea(v);
  auto meanCurvature = scalarMeanCurvature(v);
  auto gaussianCurvature = vertexGaussianCurvature(v);
  auto temp1 = meanCurvature / area;
  auto temp2 = sqrt(pow(temp1, 2) - gaussianCurvature / area);
  return std::make_pair(temp1 - temp2, temp1 + temp2);
}

/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the
 * negative semidefinite Laplace matrix, multiplying by -1, and shifting the
 * diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {
  auto vertexCount = mesh.nVertices();
  auto laplace =
      geometrycentral::SparseMatrix<double>(vertexCount, vertexCount);
  std::vector<Eigen::Triplet<double>> tripletList;
  for (const auto& vertex : mesh.vertices()) {
    auto vertexWeight{0.0};
    auto vertexIndex{vertex.getIndex()};
    for (const auto& halfedge : vertex.outgoingHalfedges()) {
      // We use -edgeCotanWeight instead of edgeCotanWeight to produce a
      // positive definite matrix
      auto adjacentVertexWeight = -edgeCotanWeight(halfedge.edge());
      tripletList.emplace_back(vertexIndex, halfedge.tipVertex().getIndex(),
                               adjacentVertexWeight);
      vertexWeight -= adjacentVertexWeight;
    }
    vertexWeight += 1e-8;  // To make this matrix strictly positive definite
    tripletList.emplace_back(vertexIndex, vertexIndex, vertexWeight);
  }
  laplace.setFromTriplets(tripletList.begin(), tripletList.end());
  return laplace;
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area
 * of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {
  auto vertexCount = mesh.nVertices();
  auto mass = geometrycentral::SparseMatrix<double>(vertexCount, vertexCount);
  std::vector<Eigen::Triplet<double>> tripletList;
  for (const auto& vertex : mesh.vertices()) {
    tripletList.emplace_back(vertex.getIndex(), vertex.getIndex(),
                             barycentricDualArea(vertex));
  }
  mass.setFromTriplets(tripletList.begin(), tripletList.end());
  return mass;
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by
 * building the negative semidefinite Laplace matrix, multiplying by -1, and
 * shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>>
VertexPositionGeometry::complexLaplaceMatrix() const {
  // TODO
  return identityMatrix<std::complex<double>>(1);  // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {
  // Compute center of mass.
  Vector3 center = {0.0, 0.0, 0.0};
  for (Vertex v : mesh.vertices()) {
    center += inputVertexPositions[v];
  }
  center /= mesh.nVertices();

  return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {
  // Compute center of mass.
  Vector3 center = centerOfMass();

  // Translate to origin [of original mesh].
  double radius = 0;
  for (Vertex v : mesh.vertices()) {
    inputVertexPositions[v] -= center;
    radius = std::max(radius, inputVertexPositions[v].norm());
  }

  // Rescale.
  if (rescale) {
    for (Vertex v : mesh.vertices()) {
      inputVertexPositions[v] /= radius;
    }
  }

  // Translate to origin [of original mesh].
  for (Vertex v : mesh.vertices()) {
    inputVertexPositions[v] += origin;
  }
}

}  // namespace surface
}  // namespace geometrycentral