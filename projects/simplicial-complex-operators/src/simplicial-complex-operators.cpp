// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer
 * to the mesh via <mesh>. Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {
  // Needed to access geometry->vertexIndices, etc. as cached quantities.
  // Not needed if you're just using v->getIndex(), etc.
  geometry->requireVertexIndices();
  geometry->requireEdgeIndices();
  geometry->requireFaceIndices();

  // You can set the index field of a vertex via geometry->vertexIndices[v],
  // where v is a Vertex object (or an integer). Similarly you can do edges and
  // faces via geometry->edgeIndices, geometry->faceIndices, like so:

  for (size_t i = 0; i < mesh->nVertices(); i++) {
    geometry->vertexIndices[i] = i;
  }

  // This is also a valid way to access indices (directly by mesh element).
  size_t idx = 0;
  for (Vertex v : mesh->vertices()) {
    geometry->vertexIndices[v] = idx;
    idx++;
  }

  for (size_t i = 0; i < mesh->nEdges(); i++) {
    geometry->edgeIndices[i] = i;
  }

  for (size_t i = 0; i < mesh->nFaces(); i++) {
    geometry->faceIndices[i] = i;
  }

  // You can more easily get the indices of mesh elements using the function
  // getIndex(), like so:
  //
  //      v.getIndex()
  //
  // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

  for (const Vertex &v : mesh->vertices()) {
    idx = v.getIndex();  // == geometry->vertexIndices[v])
  }

  // Geometry Central already sets the indices for us, though, so this function
  // is just here for demonstration. You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the
 * global variable A0.
 */
SparseMatrix<size_t>
SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
  auto edgeCount = mesh->nEdges();
  auto vertexCount = mesh->nVertices();
  std::vector<Eigen::Triplet<size_t>> tripletList;
  tripletList.reserve(2 * edgeCount);
  auto mat = geometrycentral::SparseMatrix<size_t>(edgeCount, vertexCount);
  for (auto const &edge : mesh->edges()) {
    auto edgeIndex = edge.getIndex();
    auto firstVertexIndex = edge.firstVertex().getIndex();
    auto secondVertexIndex = edge.secondVertex().getIndex();
    tripletList.emplace_back(edgeIndex, firstVertexIndex, 1);
    tripletList.emplace_back(edgeIndex, secondVertexIndex, 1);
  }
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return mat;
  // Note: You can build an Eigen sparse matrix from triplets, then return it as
  // a Geometry Central SparseMatrix. See
  // <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for
  // documentation.
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the
 * global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix()
    const {
  auto edgeCount = mesh->nEdges();
  auto faceCount = mesh->nFaces();
  std::vector<Eigen::Triplet<size_t>> tripletList;
  tripletList.reserve(3 * faceCount);
  auto mat = geometrycentral::SparseMatrix<size_t>(faceCount, edgeCount);
  for (auto const &face : mesh->faces()) {
    auto faceIndex = face.getIndex();
    for (auto const &edge : face.adjacentEdges()) {
      tripletList.emplace_back(faceIndex, edge.getIndex(), 1);
    }
  }
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(
    const MeshSubset &subset) const {
  auto vertexCount = mesh->nVertices();
  Vector<size_t> vec = Vector<size_t>::Zero(vertexCount);
  for (auto const &vertexIndex : subset.vertices) {
    vec(vertexIndex) = 1;
  }
  return vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(
    const MeshSubset &subset) const {
  auto edgeCount = mesh->nEdges();
  Vector<size_t> vec = Vector<size_t>::Zero(edgeCount);
  for (auto const &edgeIndex : subset.edges) {
    vec(edgeIndex) = 1;
  }
  return vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(
    const MeshSubset &subset) const {
  auto faceCount = mesh->nFaces();
  Vector<size_t> vec = Vector<size_t>::Zero(faceCount);
  for (auto const &faceIndex : subset.faces) {
    vec(faceIndex) = 1;
  }
  return vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active
 * vertices, edges, and faces, respectively. Returns: The star of the given
 * subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset &subset) const {
  auto vertices = subset.vertices;
  auto edges = subset.edges;
  auto faces = subset.faces;
  //  SparseMatrix<size_t> faceVectorAdjacencyMatrix = A1 * A0;
  for (auto const &vertexIndex : subset.vertices) {
    for (int edgeIndex = 0; edgeIndex < A0.rows(); ++edgeIndex) {
      // Check if the edge contains the vertex
      if (A0.coeff(edgeIndex, vertexIndex)) {
        edges.emplace(edgeIndex);
        // First solution for finding faces that contain the vertex
        // Find the Faces that contain the edge
        for (int faceIndex = 0; faceIndex < A1.rows(); ++faceIndex) {
          // Check if face contains the edge
          if (A1.coeff(faceIndex, edgeIndex)) {
            faces.emplace(faceIndex);
          }
        }
      }
      /**
      // Second solution for finding faces that contain the vertex. SLOW!
      for (int faceIndex = 0; faceIndex < faceVectorAdjacencyMatrix.rows();
      ++faceIndex) {
        // Check if face contains the vertex
        if (faceVectorAdjacencyMatrix.coeff(faceIndex, vertexIndex)) {
          faces.emplace(faceIndex);
        }
      }
      */
    }
  }
  for (auto const &edgeIndex : subset.edges) {
    for (int faceIndex = 0; faceIndex < A1.rows(); ++faceIndex) {
      // Check if face contains the edge
      if (A1.coeff(faceIndex, edgeIndex)) {
        faces.emplace(faceIndex);
      }
    }
  }
  return MeshSubset(vertices, edges, faces);
}

/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active
 * vertices, edges, and faces, respectively. Returns: The closure of the given
 * subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset &subset) const {
  auto vertices = subset.vertices;
  auto edges = subset.edges;
  auto faces = subset.faces;
  // First, add all the edges of faces
  for (auto const &faceIndex : subset.faces) {
    auto const &face = mesh->face(faceIndex);
    for (const auto &edge : face.adjacentEdges()) {
      edges.emplace(edge.getIndex());
    }
  }
  // Then, add all the vertices of edges
  for (auto const &edgeIndex : edges) {
    auto const &edge = mesh->edge(edgeIndex);
    vertices.emplace(edge.firstVertex().getIndex());
    vertices.emplace(edge.secondVertex().getIndex());
  }
  return MeshSubset(vertices, edges, faces);
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active
 * vertices, edges, and faces, respectively. Returns: The link of the given
 * subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset &subset) const {
  // TODO
  return subset;  // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active
 * vertices, edges, and faces, respectively. Returns: True if given subset is a
 * simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset &subset) const {
  // TODO
  return false;  // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the
 * degree of the complex. Otherwise, return -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active
 * vertices, edges, and faces, respectively. Returns: int representing the
 * degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset &subset) const {
  // TODO
  return -1;  // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected
 * subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active
 * vertices, edges, and faces, respectively. Returns: The boundary of the given
 * subset.
 */
MeshSubset SimplicialComplexOperators::boundary(
    const MeshSubset &subset) const {
  // TODO
  return subset;  // placeholder
}