#pragma once

#include "SavingSystem.h"

#include <Eigen/Core>

#include <Eigen/SparseCore>
#include<Eigen/SparseCholesky>	

#include <vector>

#include <iostream>
#include <numeric>

enum class Boundary
{
	CANTILEVER,
	EULER,
	FLOWING,
	DEFAULT
};

class Rod3D
{
	//typedef Eigen::Matrix<long double, 2, 1> Vector2D; 
	//typedef Eigen::Matrix<long double, 1, 1> Scalar;

	typedef long double Scalar;
	Scalar PI{ EIGEN_PI };

public:
	//types
	typedef Eigen::Matrix<Scalar, 3, 1> Vector3D;
	typedef Eigen::Matrix<Scalar, 2, 1> Vector2D;
	typedef Eigen::Matrix<Scalar, 3, 3> Matrix3D;
	typedef Eigen::Matrix<Scalar, 2, 2> Matrix2D;
	typedef Eigen::Matrix<Scalar, 3, 2> Matrix32;

	typedef Eigen::SparseMatrix<Scalar> SparseMat; // declares a column-major sparse matrix type of double
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXd;

	//methods

	Rod3D(int numOfVerts = 1, Boundary bType = Boundary::DEFAULT); //makes a straight rod at rest with vertices equally spaces from zero to 1 in the x-direction

	Rod3D(const std::vector<Vector3D>& positions, Boundary bType = Boundary::DEFAULT);
	Rod3D(const std::vector<Vector3D>& positions, const std::vector<Vector3D>& velocities, Boundary bType = Boundary::DEFAULT);

	Rod3D(const std::vector<Vector3D>& positions, const std::vector<Vector3D>& velocities, const std::vector<Scalar>& masses, const std::vector<Scalar>& tEdgeLengths, Boundary bType = Boundary::DEFAULT);


	void print(int index) { std::cout << m_vertexPositions[index] << "\n"; } // some print functionality for testing



	int getNumVertices() { return m_numOfVertices; }
	int getNumEdges() { return m_numOfEdges; }
	std::vector<Vector3D> getVertexPositions() { return m_vertexPositions; }
	std::vector<Vector3D> getVertexVelocities() { return m_vertexVelocities; }
	std::vector<Matrix2D> getFrameCurvatures() { return m_frameCurvatures; }
	std::vector<Vector3D> getCurvatureBinormals() { return m_curvatureBinormals; }


	std::vector<Scalar> getEdgeLengths() { calculateEdgeLengths(); return m_EdgeLengths; }

	void setPositions(std::vector<Vector3D> vertexPositions) { m_vertexPositions = vertexPositions; }

	void setCoefficients(Scalar bendingStiffness, Scalar stretchingStiffness, Scalar dragCoefficient);

	void SetFluidVelocity(Vector3D fluidVelocity) { m_fluidVelocity = fluidVelocity; }

	void setMass(Scalar density);

	void addGravity(Scalar gGravity, Vector3D direction = { 1, 0, 0 });

	void addPrestrain(Scalar prestrain);

	void saveState(Scalar time);

	void addReserve(int numOfReserveVerts);

	void addVertex(Vector3D newPosition, Vector3D newVelocity, Scalar density, Scalar targetLength);

	void removeVertex(Scalar density);

	void setEnterSpeed(Scalar enterSpeed) { m_enterSpeed = enterSpeed; m_vertexVelocities.assign(m_numOfVertices, { enterSpeed, 0 , 0 }); }

	void setRadius(Scalar radius) { m_rodRadius = radius; }

	void setGravityAcc(Scalar gravityAcc) { m_gravityAcc = gravityAcc;}

	void setWallPosition(Scalar wallPosition) { m_wallPosition = wallPosition; }

	void setRandomForceMag(Scalar randomForceMag) {m_randomForceMag = randomForceMag;}

	void setTargetCurvatures(std::vector<Matrix2D> barFrameCurvatures) { m_barFrameCurvatures = barFrameCurvatures; }
	// std::fill(m_barFrameCurvatures.begin(), m_barFrameCurvatures.end(), targetCurvature);

	void updatePositions(Scalar dtime); //updates positions and velocities after calculating the force

	void crossLink(Scalar timeScale);

	void startSaveFile(std::string saveName);

	void testing();


private:


	//methods

	void initializeBoundaryConds(Boundary boundaryType); // saves the initial values used in the boundary conditions

	void initializePhysicsProperties(Boundary bType);

	void calculateEdgeLengths(); // given the vertex positions, calculates the distances between consecutive vertices

	void updateKinemticFrames(); // updates the bishop frame and other kinematic quantities for given the vertex position

	void kinematicUpdate(int index);

	void boundaryConditions(Scalar dtime);

	void addExternalForce(std::vector<Vector3D>);


	void calculateForce(int index);

	Vector3D calculateBendingForce(int index);

	Vector3D calculateDragForce(int index);

	Rod3D::Vector3D calculateGravityForce(int index);

	Vector3D calculateJayVector(int index1, int index2);

	Vector3D calculateDragForce(Vector3D velocity, int index);

	//Vector3D getStretchingForce(int index);

	void constraintForce(int index);

	void inextensibility(Scalar deltaT);

	Vector3D projectOut(Vector3D vector1, Vector3D vector2);

	Matrix3D rotateTowards(Vector3D vector1, Vector3D vector2);

	Matrix3D updateFrame(Matrix3D oldFrame, Vector3D newTangent);
	Matrix3D updateBishopFrame(int index);

	Vector3D getCurvatureBinormal(int index);

	Matrix2D getFrameCurvatures(int index);

	Matrix3D getHolonomyGrads(int index);

	Matrix32 BishopAngleGrad(int index1, int index2, int index3);

	Matrix32 CurvatureGrad(int index1, int index2, int index3); // term in Eq.11 of discrete elastic rods

	Matrix3D makeTwistMatrix(Vector3D vector); // returns the matrix that correspons to M.V = vector.cross(V)

	Matrix3D randOrthoMat();


	//SavingSystem<Scalar> m_saver0{ "Run-Results/times.dat" };
	SavingSystem<Vector3D> m_saver1{ "Run-Results/positions.dat" };
	SavingSystem<Vector3D> m_saver2{ "Run-Results/velocities.dat" };
	


	//variables
	int m_numOfVertices{ 0 };
	int m_numOfEdges{ 0 };

	int m_numDimensions{ 3 };

	Scalar m_bendingStiffness{ 1 }, m_stretchingStiffness{ 1 }, m_dragCoefficient{ 1 };

	//Scalar m_dragTaylorParam{ 1 };

	std::vector<Vector3D> m_vertexPositions{}; // the 3d positions of the vertices (n by 3).

	std::vector<Vector3D> m_vertexVelocities{}; // the 3d velocities of the vertices (n by 3).

	std::vector<Scalar> m_masses{};  // masses (per unit length) at the vertices. (n by 1)

	std::vector<Vector3D> m_edgeDirections{}; // the 3d unit vectors pointing between the vertices (n-1 by 3).

	std::vector<Matrix3D> m_bishopFrames{}; // the 3by3 matrices (for each edge) represting the bishop frames.

	std::vector<Vector3D> m_curvatureBinormals{}; // the discrete curvature binormal vectors (n by 3).

	std::vector<Matrix2D> m_frameCurvatures{}; // the components of the curvature vector in the material frame (n by 2 by 2).
	std::vector<Matrix2D> m_barFrameCurvatures{}; // the components of the curvature vector in the material frame (n by 2 by 2).

	std::vector<Matrix3D> m_holonomyGradients{}; // the gradients of the angles psi that correpond to the discrerte holonomy (n by 2 by 2).

	std::vector<Scalar> m_EdgeLengths{};  // edges lengths from the edge vectors. There are n - 1 of them (n-1 by 1)

	std::vector<Scalar> m_targetEdgeLengths{};  // the preferred edge lengths. There are n - 1 of them (n-1 by 1)

	std::vector<Vector3D> m_forces{}; // same dimensionality as vertices
	std::vector<Vector3D> m_externalForces{}; // same dimensionality as vertices

	std::vector<Vector3D> m_boundaryData{}; // records the relevant quantities that are fixed with the boundary conditions

	Boundary m_boundaryType{ Boundary::DEFAULT };

	Vector3D m_fluidVelocity{ Vector3D::Zero() }; //background fluid velocity

	Scalar m_meanMasses{ 0.0 };

	Scalar m_enterSpeed{ 0 }; // speed at the boundary point when the rod is flowing in

	Scalar m_wallPosition{ 0.0 };

	bool m_hasReachedWall{ false };
	int m_floorVertexIndex{ -1 }; // index of vertex that reached the floor.


	Scalar m_terminalSpeed{ 0 };
	Scalar m_randomForceMag{ 0 };

	Scalar m_rodRadius{0};
	Scalar m_rodLength{ 0 };
	Scalar m_gravityAcc{ 9.8 };


	//std::string m_rodName{ "" };

};

