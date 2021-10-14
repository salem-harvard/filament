
//#include "SavingSystem.h"
#include "Rod3D.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <vector>
#include "SimulationMaker.h"



Rod3D::Rod3D(int numOfVerts, Boundary bType)
{
	m_numOfVertices = numOfVerts;
	m_numOfEdges = m_numOfVertices - 1;

	if (m_numOfVertices > 1) {
		for (int index{ 0 }; index < m_numOfVertices; ++index)
		{
			m_vertexPositions.push_back({ static_cast<Scalar>(index) / (m_numOfEdges), 0 , 0});
			m_vertexVelocities.push_back({ 0,0,0 });
		}
	}
	else {
		m_vertexPositions.push_back({ 0, 0, 0 });
		m_vertexVelocities.push_back({ 0,0,0 });
	}

	initializePhysicsProperties(bType);

}


Rod3D::Rod3D(const std::vector<Vector3D>& positions, Boundary bType)
{
	//TODO test whether the positions array can be modified from here. 
	// TODO check whether the variables I wanted to change like m_vertexPositions are behaving as expected
	m_vertexPositions = positions;

	m_numOfVertices = m_vertexPositions.size();
	m_numOfEdges = m_numOfVertices - 1; // is true when the rod is not closed

	m_vertexVelocities.assign(m_numOfVertices, { 0,0,0 });

	initializePhysicsProperties(bType);

	std::cout << "number of initial vertices:" << m_numOfVertices << std::endl;
}

Rod3D::Rod3D(const std::vector<Vector3D>& positions, const std::vector<Vector3D>& velocities, Boundary bType)
{

	m_vertexPositions = positions;
	m_vertexVelocities = velocities;

	m_numOfVertices = m_vertexPositions.size();
	m_numOfEdges = m_numOfVertices - 1; // is true when the rod is not closed

	initializePhysicsProperties(bType);

}


Rod3D::Rod3D(const std::vector<Vector3D>& positions, const std::vector<Vector3D>& velocities, const std::vector<Scalar>& masses, const std::vector<Scalar>& tEdgeLengths, Boundary bType)
{

	m_vertexPositions = positions;
	m_vertexVelocities = velocities;

	m_numOfVertices = m_vertexPositions.size();
	m_numOfEdges = m_numOfVertices - 1; // is true when the rod is not closed

	m_masses = masses;

	m_targetEdgeLengths = tEdgeLengths;

	calculateEdgeLengths();

	m_externalForces.assign(m_numOfVertices, { 0,0,0 });
	m_forces.assign(m_numOfVertices, { 0,0,0 }); //doesn't initialize to the correct elastic force.

	initializeBoundaryConds(bType);
}


void Rod3D::initializePhysicsProperties(Boundary bType)
{
	m_masses.assign(m_numOfVertices, 0);

	m_EdgeLengths.assign(m_numOfEdges, 0);
	m_edgeDirections.assign(m_numOfEdges, { 0,0,0 });
	calculateEdgeLengths();
	m_targetEdgeLengths = m_EdgeLengths;

	m_rodLength = std::accumulate(m_EdgeLengths.begin(), m_EdgeLengths.end(), decltype(m_EdgeLengths)::value_type(0));

	//Matrix3D initFrame = rotateTowards(Matrix3D::Identity().col(0), m_edgeDirections[0]) * Matrix3D::Identity();
	m_bishopFrames.assign(m_numOfEdges, randOrthoMat());
	m_curvatureBinormals.assign(m_numOfVertices, Vector3D::Zero());
	m_frameCurvatures.assign(m_numOfVertices, Matrix2D::Zero());
	m_barFrameCurvatures.assign(m_numOfVertices, Matrix2D::Zero());
	m_holonomyGradients.assign(m_numOfVertices, Matrix3D::Zero());


	updateKinemticFrames();

	//m_barFrameCurvatures = m_frameCurvatures;
	//for (int i = 0; i < m_barFrameCurvatures.size(); ++i) { m_barFrameCurvatures[i] = m_frameCurvatures[i] / 2.0; }
	

	m_externalForces.assign(m_numOfVertices, Vector3D::Zero());
	m_forces.assign(m_numOfVertices, Vector3D::Zero()); //doesn't initialize to the correct elastic force, still needs to be calculated.

	initializeBoundaryConds(bType);
}

void Rod3D::calculateEdgeLengths()
{
	for (int index{ 0 }; index < m_numOfEdges; ++index)
	{
		//not yet unit vectors first extract the edge lengths from this
		//m_edgeDirections.push_back({ m_vertexPositions[index + 1][1] - m_vertexPositions[index][1], m_vertexPositions[index + 1][2] - m_vertexPositions[index][2]});

		m_edgeDirections[index] = (m_vertexPositions[index + 1] - m_vertexPositions[index]);


		m_EdgeLengths[index] = (m_edgeDirections[index].norm());

		//make sure the we have unit vectors
		//m_edgeDirections[index] = { m_edgeDirections[index][1]/ m_EdgeLengths[index], m_edgeDirections[index][2] / m_EdgeLengths[index] };
		m_edgeDirections[index].normalize();
	}
}


void Rod3D::updateKinemticFrames()
{
	

	for (int index{ 0 }; index < m_numOfEdges; ++index)
	{

		kinematicUpdate(index);

	}
}

void Rod3D::kinematicUpdate(int index) {
	m_bishopFrames[index] = updateBishopFrame(index);

	m_curvatureBinormals[index] = getCurvatureBinormal(index);

	m_frameCurvatures[index] = getFrameCurvatures(index);

	m_holonomyGradients[index] = getHolonomyGrads(index);
}

void Rod3D::setCoefficients(Scalar bendingStiffness, Scalar stretchingStiffness, Scalar dragCoefficient)
{
	m_stretchingStiffness = stretchingStiffness;
	m_bendingStiffness = bendingStiffness;
	m_dragCoefficient = dragCoefficient;
}

void Rod3D::setMass(Scalar density)
{
	//calculateEdgeLengths();

	m_masses[0] = (density * m_targetEdgeLengths[0] * 0.5);

	for (int index{ 1 }; index < m_numOfVertices - 1; ++index)
	{
		long double edgeLength{ (m_targetEdgeLengths[index] + m_targetEdgeLengths[index - 1]) * 0.5 };
		m_masses[index] = (density * edgeLength);
	}


	m_masses[m_numOfVertices - 1] = (density * m_EdgeLengths[m_numOfEdges - 1] * 0.5);

	m_meanMasses = std::accumulate(m_masses.begin(), m_masses.end(), 0.0) / m_masses.size();
}


// saves the initial values, like end point positions, used in the boundary conditions
void Rod3D::initializeBoundaryConds(Boundary boundaryType)
{
	m_boundaryType = boundaryType;

	switch (m_boundaryType)
	{
	case Boundary::EULER:
		m_boundaryData.push_back(m_vertexPositions[0]);
		m_boundaryData.push_back(m_vertexPositions[m_numOfVertices - 1]);
		//m_boundaryData.push_back(m_vertexVelocities[0]);
		//m_boundaryData.push_back(m_vertexVelocities[m_numOfVertices - 1]);
		break;

	case Boundary::CANTILEVER: // cantilever
		m_boundaryData.push_back(m_vertexPositions[0]);
		m_boundaryData.push_back(m_vertexPositions[1]);
		break;
	}
}

void Rod3D::addExternalForce(std::vector<Vector3D> externalForce)
{
	for (int index{ 0 }; index < m_numOfVertices; ++index)
	{
		m_externalForces[index] += externalForce[index];
	}
}

void Rod3D::addGravity(Scalar gGravity, Vector3D direction)
{
	direction = direction.normalized();

	for (int index{ 0 }; index < m_numOfVertices; ++index)
	{
		//m_externalForces[index][0] += m_masses[index] * gGravity * direction[0];
		//m_externalForces[index][1] += m_masses[index] * gGravity * direction[1];
		m_externalForces[index] += m_masses[index] * gGravity * direction;
	}
}


// change the target edge lengths by the given factor
void Rod3D::addPrestrain(Scalar prestrain)
{
	calculateEdgeLengths();

	for (int index{ 0 }; index < m_numOfEdges; ++index)
	{
		m_targetEdgeLengths[index] = m_EdgeLengths[index] * prestrain;
	}
}

// update velocity based on force and then update positions
void Rod3D::updatePositions(Scalar dtime)
{
	//std::cout << m_forces[10] << "\n";
	for (int index{ m_floorVertexIndex + 1 }; index < m_numOfVertices; ++index)
	{

		calculateForce(index);
		
		Vector3D totalForce = m_forces[index];// +m_randomForceMag * Vector3D::Random();// + m_extenalForces[index];


		//this is a step for approximating the implicit symplectic euler
		//Vector3D velocityDifference = dtime * (totalForce / m_masses[index]);
		//totalForce += calculateDragForce(velocityDifference, index);

		m_vertexVelocities[index] += (dtime / m_masses[index]) * (totalForce);

		constraintForce(index);

		m_vertexPositions[index] += m_vertexVelocities[index] * dtime;




	}

	boundaryConditions(dtime); // applies the boundary condition

	calculateEdgeLengths();
	inextensibility(dtime);
	calculateEdgeLengths();

	updateKinemticFrames(); // recalculate the edges



	// update the bishop (material) frame
}

/*
void Rod3D::calculateForce(int index)
{
	Vector3D bendingForce, stretchingForce, tangent, dragForce; // tangent is average of the two tangents at each vertex

	if (index >= 2 && index < m_numOfVertices - 2)
	{
		bendingForce = (projectOut(m_edgeDirections[index - 2], m_edgeDirections[index - 1]) / (m_EdgeLengths[index - 1] * (m_targetEdgeLengths[index - 2] + m_targetEdgeLengths[index - 1]) / 2) +
			projectOut(m_edgeDirections[index], m_edgeDirections[index - 1]) / (m_EdgeLengths[index - 1] * (m_targetEdgeLengths[index - 1] + m_targetEdgeLengths[index]) / 2) -
			projectOut(m_edgeDirections[index - 1], m_edgeDirections[index]) / (m_EdgeLengths[index] * (m_targetEdgeLengths[index - 1] + m_targetEdgeLengths[index]) / 2) -
			projectOut(m_edgeDirections[index + 1], m_edgeDirections[index]) / (m_EdgeLengths[index] * (m_targetEdgeLengths[index] + m_targetEdgeLengths[index + 1]) / 2)) * m_bendingStiffness;

		stretchingForce = -(-(m_EdgeLengths[index] - m_targetEdgeLengths[index]) * m_edgeDirections[index] / m_targetEdgeLengths[index] +
			(m_EdgeLengths[index - 1] - m_targetEdgeLengths[index - 1]) * m_edgeDirections[index - 1] / m_targetEdgeLengths[index - 1]) * m_stretchingStiffness;

		tangent = (m_edgeDirections[index - 1] + m_edgeDirections[index]).normalized();

	}
	else if (index == 0)
	{
		bendingForce = (-
			projectOut(m_edgeDirections[index + 1], m_edgeDirections[index]) / (m_EdgeLengths[index] * (m_targetEdgeLengths[index] + m_targetEdgeLengths[index + 1]) / 2)) * m_bendingStiffness;

		stretchingForce = -(-(m_EdgeLengths[index] - m_targetEdgeLengths[index]) * m_edgeDirections[index] / m_targetEdgeLengths[index]) * m_stretchingStiffness;

		tangent = m_edgeDirections[0];

	}
	else if (index == 1)
	{
		bendingForce = (projectOut(m_edgeDirections[index], m_edgeDirections[index - 1]) / (m_EdgeLengths[index - 1] * (m_targetEdgeLengths[index - 1] + m_targetEdgeLengths[index]) / 2) -
			projectOut(m_edgeDirections[index - 1], m_edgeDirections[index]) / (m_EdgeLengths[index] * (m_targetEdgeLengths[index - 1] + m_targetEdgeLengths[index]) / 2) -
			projectOut(m_edgeDirections[index + 1], m_edgeDirections[index]) / (m_EdgeLengths[index] * (m_targetEdgeLengths[index] + m_targetEdgeLengths[index + 1]) / 2)) * m_bendingStiffness;


		stretchingForce = -(-(m_EdgeLengths[index] - m_targetEdgeLengths[index]) * m_edgeDirections[index] / m_targetEdgeLengths[index] +
			(m_EdgeLengths[index - 1] - m_targetEdgeLengths[index - 1]) * m_edgeDirections[index - 1] / m_targetEdgeLengths[index - 1]) * m_stretchingStiffness;

		tangent = (m_edgeDirections[index - 1] + m_edgeDirections[index]).normalized();

	}
	else if (index == m_numOfVertices - 2)
	{
		//edge indices for the edge directions cannot exceed m_numOfVertices - 2 = m_numOfEdges - 1
		bendingForce = (projectOut(m_edgeDirections[index - 2], m_edgeDirections[index - 1]) / (m_EdgeLengths[index - 1] * (m_targetEdgeLengths[index - 2] + m_targetEdgeLengths[index - 1]) / 2) +
			projectOut(m_edgeDirections[index], m_edgeDirections[index - 1]) / (m_EdgeLengths[index - 1] * (m_targetEdgeLengths[index - 1] + m_targetEdgeLengths[index]) / 2) -
			projectOut(m_edgeDirections[index - 1], m_edgeDirections[index]) / (m_EdgeLengths[index] * (m_targetEdgeLengths[index - 1] + m_targetEdgeLengths[index]) / 2)) * m_bendingStiffness;


		stretchingForce = -(-(m_EdgeLengths[index] - m_targetEdgeLengths[index]) * m_edgeDirections[index] / m_targetEdgeLengths[index] +
			(m_EdgeLengths[index - 1] - m_targetEdgeLengths[index - 1]) * m_edgeDirections[index - 1] / m_targetEdgeLengths[index - 1]) * m_stretchingStiffness;


		tangent = (m_edgeDirections[index - 1] + m_edgeDirections[index]).normalized();

	}
	else if (index == m_numOfVertices - 1)
	{
		bendingForce = (projectOut(m_edgeDirections[index - 2], m_edgeDirections[index - 1]) / (m_EdgeLengths[index - 1] * (m_targetEdgeLengths[index - 2] +
			m_targetEdgeLengths[index - 1]) / 2)) * m_bendingStiffness;

		stretchingForce = -((m_EdgeLengths[index - 1] - m_targetEdgeLengths[index - 1]) * m_edgeDirections[index - 1] / m_targetEdgeLengths[index - 1]) * m_stretchingStiffness;

		tangent = (m_edgeDirections[index - 1]);

	}

	Vector3D velocity = m_vertexVelocities[index] - m_fluidVelocity;

	dragForce = -m_dragCoefficient * m_masses[index] * (2 * velocity - (velocity.dot(tangent)) * tangent);


	Vector3D gravityForce = { (0.0 / 1054)*(m_masses[index] * 9.8), 0, 0 }; // 30.0 / 1054.


	m_forces[index] = bendingForce;// +stretchingForce + dragForce + gravityForce;

}
*/

void Rod3D::calculateForce(int index)
{

	Vector3D bendingForce{Vector3D::Zero()}, dragForce{ Vector3D::Zero() }, gravityForce{ Vector3D::Zero() }; 



	bendingForce = calculateBendingForce(index);

	dragForce = calculateDragForce(index);

	gravityForce = calculateGravityForce(index); 



	m_forces[index] = bendingForce + dragForce + gravityForce;//+ stretchingForce;// +getStretchingForce(index);// +dragForce + gravityForce;

	//std::cout << "\n" << m_forces[index] << "\n";

}

Rod3D::Vector3D Rod3D::calculateBendingForce(int index) {

	Vector3D bendingForce{ Vector3D::Zero() };


	for (int vIndx{ 1 }; vIndx < m_numOfEdges; vIndx++) {
		bendingForce = bendingForce - m_bendingStiffness * (1 / (m_EdgeLengths[vIndx] + m_EdgeLengths[vIndx - 1])) * (
			(CurvatureGrad(index, vIndx - 1, vIndx)) * (m_frameCurvatures[vIndx].col(0) - m_barFrameCurvatures[vIndx].col(0)) +
			(CurvatureGrad(index, vIndx, vIndx)) * (m_frameCurvatures[vIndx].col(1) - m_barFrameCurvatures[vIndx].col(1)));
			//(CurvatureGrad(index, vIndx - 1, vIndx) - BishopAngleGrad(index, vIndx - 1, vIndx)) * (m_frameCurvatures[vIndx].col(0) - m_barFrameCurvatures[vIndx].col(0)) +
			//(CurvatureGrad(index, vIndx, vIndx) - BishopAngleGrad(index, vIndx, vIndx)) * (m_frameCurvatures[vIndx].col(1) - m_barFrameCurvatures[vIndx].col(1)));
	}

	return bendingForce; 
}

Rod3D::Vector3D Rod3D::calculateDragForce(int index)
{
	Vector3D jayVec{ Vector3D::Zero() }, tangent{ Vector3D::Zero() };

	Scalar vorLength{0};

	if (index == 0) {
		tangent = m_edgeDirections[0];
		vorLength = m_EdgeLengths[0] * 0.5;
	}
	else if (index == m_numOfVertices - 1) {
		tangent = m_edgeDirections[index - 1];
		vorLength = m_EdgeLengths[index - 1] * 0.5;
	}
	else {
		tangent = (m_edgeDirections[index - 1] + m_edgeDirections[index]).normalized();
		vorLength = (m_EdgeLengths[index - 1] + m_EdgeLengths[index]) * 0.5;
	}

	Vector3D velocity{  m_vertexVelocities[index] - m_fluidVelocity };

	Matrix3D M1{ 2*Matrix3D::Identity() -  tangent * tangent.transpose() };

	//Scalar totalLength = std::accumulate(m_EdgeLengths.begin(), m_EdgeLengths.end(),
	//	decltype(m_EdgeLengths)::value_type(0));

	Scalar dragTaylorParam = 1/log(m_rodLength / m_rodRadius);


	Vector3D dragForce{ Vector3D::Zero() };
	//add leading order
	dragForce += -2 * PI * m_dragCoefficient * dragTaylorParam* M1 * velocity;

	//return dragForce * vorLength;

	//add next order
	//****************
	for (int vIndx{ 0 }; vIndx < m_numOfVertices; vIndx++) {
		jayVec += calculateJayVector(index, vIndx);
	}
	Matrix3D M2{ Matrix3D::Identity() - 1.5 * tangent * tangent.transpose() };

	dragForce += 2 * PI * dragTaylorParam * dragTaylorParam * m_dragCoefficient * M1 * (jayVec + velocity*log(vorLength));
	dragForce += 2 * PI * dragTaylorParam * dragTaylorParam * m_dragCoefficient * M2 * velocity;
	//****************

	return dragForce*vorLength;

}

Rod3D::Vector3D Rod3D::calculateGravityForce(int index) {
	return { m_gravityAcc * m_masses[index], 0, 0 };
}

Rod3D::Vector3D Rod3D::calculateJayVector(int index1, int index2) {

	if (index1 == index2) { return Vector3D::Zero(); }

	Vector3D tangent2{ Vector3D::Zero() };

	Scalar vorLength{ 0 };
	if (index2 == 0) { tangent2 = m_edgeDirections[0];
	vorLength = m_EdgeLengths[0] * 0.5;
	}
	else if (index2 == m_numOfVertices - 1) { tangent2 = (m_edgeDirections[index2 - 1]);
	vorLength = m_EdgeLengths[index2 - 1] * 0.5;
	}
	else { tangent2 = (m_edgeDirections[index2 - 1] + m_edgeDirections[index2]).normalized();
	vorLength = (m_EdgeLengths[index2 - 1] + m_EdgeLengths[index2]) * 0.5;
	}

	Vector3D velocity2{ m_vertexVelocities[index2] - m_fluidVelocity };

	Vector3D displacement{ m_vertexPositions[index1] - m_vertexPositions[index2] };

	Matrix3D M1{ Matrix3D::Identity() + displacement.normalized() * displacement.normalized().transpose() };
	Matrix3D M2{ Matrix3D::Identity() - 0.5 * tangent2 * tangent2.transpose() };

	return 0.5 * (M1 * M2 * velocity2)*(vorLength / displacement.norm());

}

Rod3D::Vector3D Rod3D::calculateDragForce(Vector3D velocity, int index)
{

	Vector3D  dragForce{ Vector3D::Zero() }, tangent{ Vector3D::Zero() }; 

	if (index == 0) {tangent = m_edgeDirections[0];}
	else if (index == m_numOfVertices - 1) {tangent = (m_edgeDirections[index - 1]);}
	else {tangent = (m_edgeDirections[index - 1] + m_edgeDirections[index]).normalized();}



	dragForce = -m_dragCoefficient * m_masses[index] * (2 * velocity - (velocity.dot(tangent)) * tangent);

	return dragForce;

}

/*
Rod3D::Vector3D Rod3D::getStretchingForce(int index)
{
	Vector3D  stretchingForce{ Vector3D::Zero() }; // tangent is average of the two tangents at each vertex


	if (index == 0)
	{

		stretchingForce = -(-(m_EdgeLengths[index] - m_targetEdgeLengths[index]) * m_edgeDirections[index] / m_targetEdgeLengths[index]) * m_stretchingStiffness;

	}
	else if (index == m_numOfVertices - 1)
	{
		
		stretchingForce = -((m_EdgeLengths[index - 1] - m_targetEdgeLengths[index - 1]) * m_edgeDirections[index - 1] / m_targetEdgeLengths[index - 1]) * m_stretchingStiffness;


	}
	else {


			stretchingForce = -(-(m_EdgeLengths[index] - m_targetEdgeLengths[index]) * m_edgeDirections[index] / m_targetEdgeLengths[index] +
				(m_EdgeLengths[index - 1] - m_targetEdgeLengths[index - 1]) * m_edgeDirections[index - 1] / m_targetEdgeLengths[index - 1]) * m_stretchingStiffness;

	}

	return stretchingForce;

}
*/

void Rod3D::constraintForce(int index)
{
	Scalar epsilon{ 0.001 };

	if (m_vertexPositions[index][0] > m_wallPosition - epsilon) {
		
		//std::cout << "crossing wall" << std::endl;

		// for hitting the floor. 
		m_vertexPositions[index][0] = m_wallPosition;
		m_vertexVelocities[index] = {0,0,0};

		m_floorVertexIndex = index; 

		/*
		if (!m_hasReachedWall) 
		{ 
			m_hasReachedWall = true; 
			m_terminalSpeed = m_vertexVelocities[index][0];
		}

		// for reaching steady state 
		m_vertexVelocities[index] = { m_terminalSpeed, 0,0}; */
	}

	/*
	//enforcing inextensibility constraint
	if (index > 0) {
		Vector3D edgeDirection{ (m_vertexPositions[index] - m_vertexPositions[index - 1]).normalized() };
		m_vertexVelocities[index] = projectOut(m_vertexVelocities[index], edgeDirection) + (m_vertexVelocities[index - 1].dot(edgeDirection)) * edgeDirection;
		m_vertexPositions[index] = m_vertexPositions[index - 1] + edgeDirection*m_EdgeLengths[index - 1];
	}
	*/
	

}

void Rod3D::inextensibility(Scalar deltaT) {
	// construct the relevant matrices, solve the sparse system
	// add the solution to the vertices. 

	SparseMat massMatrix(m_numDimensions * m_numOfVertices, m_numDimensions * m_numOfVertices);
	SparseMat constraintGrad(m_numOfEdges, m_numDimensions * m_numOfVertices); 
	VectorXd constraints(m_numOfEdges, 1), deltaLambda(m_numOfEdges, 1), deltaX(m_numDimensions * m_numOfVertices, 1);

	for (int indx{ 0 }; indx < m_numOfEdges; ++indx) {
		massMatrix.insert(indx * m_numDimensions, indx * m_numDimensions) = m_meanMasses/m_masses[indx];
		massMatrix.insert(indx * m_numDimensions + 1, indx * m_numDimensions + 1) = m_meanMasses / m_masses[indx];
		massMatrix.insert(indx * m_numDimensions + 2, indx * m_numDimensions + 2) = m_meanMasses / m_masses[indx];


		constraints[indx] = m_EdgeLengths[indx] * m_EdgeLengths[indx]/m_targetEdgeLengths[indx] - m_targetEdgeLengths[indx];

		constraintGrad.insert(indx, indx * m_numDimensions) = -2*m_bishopFrames[indx].col(0)[0];
		constraintGrad.insert(indx, indx * m_numDimensions + 1) = -2 * m_bishopFrames[indx].col(0)[1];
		constraintGrad.insert(indx, indx * m_numDimensions + 2) = -2 * m_bishopFrames[indx].col(0)[2];
		constraintGrad.insert(indx, (indx + 1) * m_numDimensions) = 2 * m_bishopFrames[indx].col(0)[0];
		constraintGrad.insert(indx, (indx + 1) * m_numDimensions + 1) = 2 * m_bishopFrames[indx].col(0)[1];
		constraintGrad.insert(indx, (indx + 1) * m_numDimensions + 2) = 2 * m_bishopFrames[indx].col(0)[2];

	}
	massMatrix.insert(m_numOfEdges * m_numDimensions, m_numOfEdges * m_numDimensions) = m_meanMasses / m_masses[m_numOfEdges];
	massMatrix.insert(m_numOfEdges * m_numDimensions + 1, m_numOfEdges * m_numDimensions + 1) = m_meanMasses / m_masses[m_numOfEdges];
	massMatrix.insert(m_numOfEdges * m_numDimensions + 2, m_numOfEdges * m_numDimensions + 2) = m_meanMasses / m_masses[m_numOfEdges];

	//add the boundary condition constraint
	/*
	constraints[m_numOfEdges] = m_edgeDirections[0].dot(m_edgeDirections[1]) - 1;

	constraintGrad.insert(m_numOfEdges, 0) = - m_bishopFrames[1].col(0)[0]/m_targetEdgeLengths[0];
	constraintGrad.insert(m_numOfEdges, 1) = -m_bishopFrames[1].col(0)[1] / m_targetEdgeLengths[0];
	constraintGrad.insert(m_numOfEdges, 2) = -m_bishopFrames[1].col(0)[2] / m_targetEdgeLengths[0];

	constraintGrad.insert(m_numOfEdges, 3) = m_bishopFrames[1].col(0)[0] / m_targetEdgeLengths[0] - m_bishopFrames[0].col(0)[0] / m_targetEdgeLengths[1];
	constraintGrad.insert(m_numOfEdges, 4) = m_bishopFrames[1].col(0)[1] / m_targetEdgeLengths[0] -m_bishopFrames[0].col(0)[1] / m_targetEdgeLengths[1];
	constraintGrad.insert(m_numOfEdges, 5) = m_bishopFrames[1].col(0)[2] / m_targetEdgeLengths[0] -m_bishopFrames[0].col(0)[2] / m_targetEdgeLengths[1];

	constraintGrad.insert(m_numOfEdges, 6) = -m_bishopFrames[0].col(0)[0] / m_targetEdgeLengths[1];
	constraintGrad.insert(m_numOfEdges, 7) = -m_bishopFrames[0].col(0)[1] / m_targetEdgeLengths[1];
	constraintGrad.insert(m_numOfEdges, 8) = -m_bishopFrames[0].col(0)[2] / m_targetEdgeLengths[01];
	*/


	SparseMat A = deltaT * deltaT * constraintGrad * massMatrix * (constraintGrad.transpose());

	Eigen::SimplicialLLT<SparseMat > solver;
	solver.compute(A);
	
	deltaLambda = solver.solve(constraints);

	deltaX = -deltaT * deltaT * massMatrix * (constraintGrad.transpose())* deltaLambda;


	for (int indx{ m_floorVertexIndex + 1 }; indx < m_numOfVertices; ++indx) {
		m_vertexVelocities[indx] += Vector3D{deltaX[indx * m_numDimensions], deltaX[indx * m_numDimensions + 1], deltaX[indx * m_numDimensions +2]}/ deltaT;
		m_vertexPositions[indx] += Vector3D{ deltaX[indx * m_numDimensions], deltaX[indx * m_numDimensions + 1], deltaX[indx * m_numDimensions + 2] } ;
	}
	

}

void Rod3D::boundaryConditions(Scalar deltaT)
{
	
	switch (m_boundaryType)
	{
	case Boundary::EULER:
	{
		m_vertexPositions[0] = m_boundaryData[0];
		m_vertexPositions[m_numOfVertices - 1] = m_boundaryData[1];

		m_vertexVelocities[0] = { 0, 0 ,0};
		m_vertexVelocities[m_numOfVertices - 1] = { 0, 0, 0};
	}
		break;

	case Boundary::FLOWING: {

		Vector3D velocityDifference = Vector3D{ m_enterSpeed, 0, 0 } - m_vertexVelocities[m_numOfVertices - 1];

		//m_vertexPositions[m_numOfVertices - 1] += velocityDifference * deltaT;
		m_vertexVelocities[m_numOfVertices - 1] += velocityDifference;
		

	}
	break;

	case Boundary::CANTILEVER: // 
	{
		m_vertexPositions[0] = m_boundaryData[0];
		m_vertexPositions[1] = m_boundaryData[1];

		m_vertexVelocities[0] = { 0,0,0 };
		m_vertexVelocities[1] = { 0,0,0 };
	}
	break;

	case Boundary::DEFAULT: // 
	{
		//free boundary
	}
	break;
	}

}

//takes 2 unit vectors and finds the part of unitVector1 that is normal to unitVector2
Rod3D::Vector3D Rod3D::projectOut(Vector3D vector1, Vector3D unitVector2)
{
	return vector1 - vector1.dot(unitVector2) * unitVector2;
}


// a matrix that rotates vector1 to vector2 about an axis perpendicular to both of them (discrete parallel transport).
Rod3D::Matrix3D Rod3D::rotateTowards(Vector3D unitVec1, Vector3D unitVec2)
{
	Vector3D crossProduct{ unitVec1.cross(unitVec2) };

	Scalar sinTheta{crossProduct.norm()}, cosTheta{ unitVec1.dot(unitVec2) };

	Scalar epsilon{ 0.00000000001 }; // if vectors are nearly parallel return the identity
	if (abs(sinTheta) < epsilon && cosTheta > epsilon) { return Matrix3D::Identity(); }
	if (abs(sinTheta) < epsilon && cosTheta < epsilon) { std::cerr << "Sharp turning angle" << std::endl; }

	//normal axes to unitVec1, normal1 is the rotation axis
	Vector3D normal1{crossProduct.normalized()};
	//perpendicular to both the axes of rotation and unitVec1
	Vector3D normal2{ normal1.cross(unitVec1) };

	return normal1*normal1.transpose() + cosTheta*(unitVec1 * unitVec1.transpose() + normal2 * normal2.transpose()) + 
		sinTheta*(normal2 * unitVec1.transpose() - unitVec1 * normal2.transpose());
}

Rod3D::Matrix3D Rod3D::updateFrame(Matrix3D oldFrame, Vector3D newTangent) {
	return rotateTowards(oldFrame.col(0), newTangent)*oldFrame;
}

Rod3D::Matrix3D Rod3D::updateBishopFrame(int index) {

	
	Vector3D newTangent = (m_vertexPositions[index + 1] - m_vertexPositions[index]).normalized();
	Matrix3D oldFrame{};

	if (index == 0) {
		oldFrame = m_bishopFrames[0];	
	}
	else {
		oldFrame = m_bishopFrames[index - 1];
	}

	Matrix3D newFrame{ rotateTowards(oldFrame.col(0), newTangent) * oldFrame };

	newFrame.col(0) = newTangent;

	return newFrame;
}

Rod3D::Vector3D Rod3D::getCurvatureBinormal(int index)
{
	if (index == 0) { return Vector3D::Zero(); }
	
	Vector3D tcross{Vector3D::Zero()};
	Scalar tdot{ 0 };
	
	tcross = (m_bishopFrames[index - 1].col(0)).cross(m_bishopFrames[index].col(0));
	tdot = (m_bishopFrames[index - 1].col(0)).dot(m_bishopFrames[index].col(0));

	return 2*tcross/(tdot + 1);
}

Rod3D::Matrix2D Rod3D::getFrameCurvatures(int index)
{
	if(index == 0) { return Matrix2D::Zero(); }
	

	Vector3D curvBinormal{ m_curvatureBinormals[index] };
	Matrix2D frameCurv{};

	

	frameCurv(0, 0) =  curvBinormal.dot(m_bishopFrames[index - 1].col(2));
	frameCurv(1, 0) = -curvBinormal.dot(m_bishopFrames[index - 1].col(1));
	frameCurv(0, 1) =  curvBinormal.dot(m_bishopFrames[index].col(2));
	frameCurv(1, 1) = -curvBinormal.dot(m_bishopFrames[index].col(1));
	

	return frameCurv;

}

Rod3D::Matrix3D Rod3D::getHolonomyGrads(int index)
{
	if (index == 0) { return Matrix3D::Zero(); }

	Matrix3D hGrad{};

	hGrad.col(0) = 0.5*m_curvatureBinormals[index] / m_EdgeLengths[index - 1];
	hGrad.col(2) = -0.5*m_curvatureBinormals[index] / m_EdgeLengths[index];
	hGrad.col(1) = -(hGrad.col(0) + hGrad.col(2));

	return hGrad;
}


//grad_{index1} Psi_{index_2}
Rod3D::Matrix32 Rod3D::BishopAngleGrad(int index1, int index2, int index3) 
{
	Vector3D angleGrad{ Vector3D::Zero()};
	
	if ((index2 < index1 - 1) || index2 == 0) { angleGrad = Vector3D(); }
	else if (index2 == index1 - 1) { angleGrad =  m_holonomyGradients[index1 - 1].col(2); }
	else if (index2 == index1 - 1) { angleGrad =  m_holonomyGradients[index1 - 1].col(2) + m_holonomyGradients[index1].col(1); }
	else { angleGrad = m_holonomyGradients[index1 - 1].col(2) + m_holonomyGradients[index1].col(1) + m_holonomyGradients[index1 + 1].col(0); }

	Matrix2D quarterRotation{}; 
	quarterRotation << 0, -1,
					  1, 0;

	Vector2D frameCurvature{};

	if (index3 == index2 + 1) { frameCurvature = quarterRotation*m_frameCurvatures[index3].col(0); }
	if (index3 == index2) { frameCurvature = quarterRotation * m_frameCurvatures[index3].col(1); }
	
	return angleGrad * (frameCurvature.transpose());

}

Rod3D::Matrix32 Rod3D::CurvatureGrad(int index1, int index2, int index3)
{
	Matrix3D binormalGrad{ Matrix3D::Zero()};

	if (index3 > index1 + 1 || index3 < index1 - 1 || index3 == 0 || index3 == m_numOfVertices -  1) {
		return Matrix32::Zero();
	}
	else if (index3 == index1 - 1) {
		binormalGrad = (2 * makeTwistMatrix(m_bishopFrames[index3 - 1].col(0)) - m_curvatureBinormals[index3] * (m_bishopFrames[index3 - 1].col(0).transpose())) / (
			m_EdgeLengths[index3] + m_EdgeLengths[index3] * m_bishopFrames[index3 - 1].col(0).dot(m_bishopFrames[index3].col(0)));
	}
	else if (index3 == index1 + 1) {
		binormalGrad = (2 * makeTwistMatrix(m_bishopFrames[index3].col(0)) + m_curvatureBinormals[index3] * (m_bishopFrames[index3].col(0).transpose())) / (
			m_EdgeLengths[index1] + m_EdgeLengths[index1] * m_bishopFrames[index1].col(0).dot(m_bishopFrames[index3].col(0)));
	}
	else if (index3 == index1) {
		binormalGrad = -(2 * makeTwistMatrix(m_bishopFrames[index3].col(0)) + m_curvatureBinormals[index3] * (m_bishopFrames[index3].col(0).transpose())) / (
			m_EdgeLengths[index3 - 1] + m_EdgeLengths[index3 - 1] * m_bishopFrames[index3-1].col(0).dot(m_bishopFrames[index3].col(0)));

		binormalGrad += -(2 * makeTwistMatrix(m_bishopFrames[index3-1].col(0)) - m_curvatureBinormals[index3] * (m_bishopFrames[index3-1].col(0).transpose())) / (
			m_EdgeLengths[index3] + m_EdgeLengths[index3] * m_bishopFrames[index3 - 1].col(0).dot(m_bishopFrames[index3].col(0)));
	}

	Matrix32 frameMat{};

	frameMat.col(0) = m_bishopFrames[index2].col(2); frameMat.col(1) = -m_bishopFrames[index2].col(1);
	
	return binormalGrad.transpose() * frameMat;
}

//update natural curvature and rigidity
void Rod3D::crossLink(Scalar timeScale) {

	for (int index{ 2 }; index < m_numOfVertices - 2; ++index) {
		m_barFrameCurvatures[index] = m_barFrameCurvatures[index] - timeScale * (m_barFrameCurvatures[index] - m_frameCurvatures[index]);
	}
	return;
}


//matrix that corresponds M.x is vector.cross(x)
Rod3D::Matrix3D Rod3D::makeTwistMatrix(Vector3D vector)
{
	Matrix3D twistMat{};

	twistMat << 0, -vector(2), vector(1),
		vector(2), 0, -vector(0),
		-vector(1), vector(0), 0;
	
	return  twistMat;
}

Rod3D::Matrix3D Rod3D::randOrthoMat() {
	Eigen::Quaterniond q = Eigen::Quaterniond::UnitRandom();

	Matrix3D mat;

	mat = q.toRotationMatrix().cast<Scalar>();

	//return Matrix3D::Random();

	return mat;
}

void Rod3D::addReserve(int numOfReserveVerts)
{
	m_vertexPositions.reserve(numOfReserveVerts);
	m_vertexVelocities.reserve(numOfReserveVerts);
	m_masses.reserve(numOfReserveVerts);

	m_edgeDirections.reserve(numOfReserveVerts - 1);
	m_EdgeLengths.reserve(numOfReserveVerts - 1);
	m_targetEdgeLengths.reserve(numOfReserveVerts - 1);


	m_bishopFrames.reserve(numOfReserveVerts - 1);
	m_curvatureBinormals.reserve(numOfReserveVerts);
	m_frameCurvatures.reserve(numOfReserveVerts);
	m_barFrameCurvatures.reserve(numOfReserveVerts);
	m_holonomyGradients.reserve(numOfReserveVerts);


	m_forces.reserve(numOfReserveVerts);
	m_externalForces.reserve(numOfReserveVerts);
}

void Rod3D::addVertex(Vector3D newPosition, Vector3D newVelocity, Scalar density, Scalar targetLength)
{

	// add the edge information
	Vector3D edgeVector{ newPosition - m_vertexPositions[m_numOfVertices - 1] };

	m_vertexPositions.push_back(newPosition);
	m_vertexVelocities.push_back(newVelocity);

	Scalar newEdgeLength{ edgeVector.norm() };

	m_EdgeLengths.push_back(newEdgeLength);
	m_edgeDirections.push_back(edgeVector.normalized());


	// whether the target length can have initial incompatibility
	//targetLength =  newEdgeLength;//Scalar targetLength{newEdgeLength};
	m_targetEdgeLengths.push_back(newEdgeLength);
	m_rodLength += newEdgeLength;

	//change the mass. The previous vertex needs to be changed as well due to voronoi cell change
	m_masses[m_numOfVertices - 1] += 0.5 * targetLength * density;
	m_masses.push_back(0.5 * targetLength * density);

	//not adding external force so far
	m_externalForces.push_back({ 0,0,0 });

	//extending force size, will be calculated in the update method
	m_forces.push_back({ 0,0,0 });

	
	m_bishopFrames.push_back(randOrthoMat());
	m_curvatureBinormals.push_back(Vector3D::Zero());
	m_frameCurvatures.push_back(Matrix2D::Zero());
	m_holonomyGradients.push_back(Matrix3D::Zero());


	m_numOfVertices += 1;
	m_numOfEdges += 1;

	kinematicUpdate(m_numOfEdges - 1);

	m_barFrameCurvatures.push_back(Matrix2D::Zero());
	//m_barFrameCurvatures.push_back(m_frameCurvatures.back());
}

//removes a vertex from the end of the array
void Rod3D::removeVertex(Scalar density) {


	m_vertexPositions.erase(m_vertexPositions.begin());
	m_vertexVelocities.erase(m_vertexVelocities.begin());


	//The first remaining mass needs to be changed as well before removing the target length
	m_masses.erase(m_masses.begin());
	m_masses[0] -= 0.5 * m_targetEdgeLengths[0] * density;

	m_rodLength -= m_targetEdgeLengths[0];


	m_EdgeLengths.erase(m_EdgeLengths.begin());
	m_targetEdgeLengths.erase(m_targetEdgeLengths.begin());
	m_edgeDirections.erase(m_edgeDirections.begin());



	m_externalForces.erase(m_externalForces.begin());
	m_forces.erase(m_forces.begin());


	m_bishopFrames.erase(m_bishopFrames.begin());
	m_curvatureBinormals.erase(m_curvatureBinormals.begin());
	m_frameCurvatures.erase(m_frameCurvatures.begin());
	m_barFrameCurvatures.erase(m_barFrameCurvatures.begin());
	m_holonomyGradients.erase(m_holonomyGradients.begin());


	m_numOfVertices -= 1;
	m_numOfEdges -= 1;

}


void Rod3D::saveState(Scalar time)
{
	//m_saver0.saveElement(time);
	m_saver1.saveVectors(m_vertexPositions, m_numOfVertices);

	m_saver2.saveVectors(m_vertexVelocities, m_numOfVertices);
	
}


void Rod3D::startSaveFile(std::string saveName) {
	//m_saver0.startFile("Run-Results/times" + saveName + ".dat");
	m_saver1.startFile("Run-Results/positions" + saveName + ".dat");
	m_saver2.startFile("Run-Results/velocities" + saveName + ".dat");
}

void Rod3D::testing()
{
	//std::cout << m_curvatureBinormals.size() <<std::endl;



	//Vector3D vec1{ m_bishopFrames[indx - 1].col(0) }, vec2{ m_bishopFrames[indx].col(0) };

	
	


	int index{ 10 }, vIndx{ 40 };

	//Vector3D knormal{ 2 * vec1.cross(vec2) / (1 + vec1.dot(vec2)) };

	Scalar vorLength = 0.5 * (m_EdgeLengths[index-1] + m_EdgeLengths[index]);

	std::cout << (m_frameCurvatures[index] - m_barFrameCurvatures[index])/vorLength << std::endl;

	std::cout << "force: \n" << m_forces[index]/vorLength << std::endl << std::endl;

	//std::cout <<  m_bendingStiffness << std::endl << std::endl;
	//std::cout << (m_EdgeLengths[vIndx] + m_EdgeLengths[vIndx - 1]) << std::endl << std::endl;

	//std::cout << m_vertexPositions[index] << std::endl << std::endl;


	//std::cout <<( CurvatureGrad(index, index, index) + CurvatureGrad(index, index - 1, index) +
	//	CurvatureGrad(index, index - 1, index-1) + CurvatureGrad(index, index - 2, index - 1) + 
	//	CurvatureGrad(index, index + 1, index + 1) + CurvatureGrad(index, index, index + 1)
	//	)* vorLength << std::endl << std::endl;
	

}


