#include "SimulationMaker.h"
#include "Rod3D.h"
#include "SavingSystem.h"


#include <cmath>
#include <vector>
#include <string>

#include <cstdlib>






void SimulationMaker::flowingIntoFluid(Scalar finalTime, int numOfTimeSteps, int vertsPerLength, int maxNumOfVerts, Scalar enterSpeed, Scalar crossLinkRate, Scalar wallPosition, Scalar density, Scalar YModulus,
											Scalar muDragCoeff, Scalar radius, Scalar gravityAcc, int numSaves, Scalar bias, std::string saveName)
{
	
	Scalar length {enterSpeed * finalTime };


	int numOfVerts{ static_cast<int>(vertsPerLength * length) };

	assert(numOfTimeSteps > numOfVerts);

	int reducedNumVerts{ static_cast<int>(numOfVerts - 2) }; // since we start with having two vertices

	//change the number of time steps to accomodate all the inserted vertices.
	numOfTimeSteps += reducedNumVerts - numOfTimeSteps % reducedNumVerts;

	Scalar targetLength{ length / (numOfVerts - 1.0) }; // target length of every incoming edge, assumed constant

	//{8*targetLength, 0, 0},{7 * targetLength, 0 ,0},{6 * targetLength, 0, 0},{5 * targetLength, 0 ,0},{4 * targetLength, 0, 0}, { 3 * targetLength, 0 ,0 },
	Rod3D rod{ {{8 * targetLength, 0, 0},{7 * targetLength, 0 ,0},{6 * targetLength, 0, 0},{5 * targetLength, 0 ,0},
		{4 * targetLength, 0, 0}, { 3 * targetLength, 0 ,0 },{2 * targetLength, 0, 0},{targetLength, 0 ,0}},
		{{enterSpeed, 0, 0},{enterSpeed, 0, 0},{enterSpeed, 0, 0},{enterSpeed, 0, 0},{enterSpeed, 0, 0},
		{enterSpeed, 0, 0},{enterSpeed, 0, 0},{enterSpeed, 0, 0}}, Boundary::FLOWING };

	//when testing gravity and drag on a 4 vertex rod
	//Rod3D rod{ {{4 * targetLength, 0, 0} ,{3 * targetLength, 0, 0}, {2 * targetLength, 0, 0}, {targetLength, 0 ,0}}, Boundary::FLOWING };

	rod.addReserve(numOfVerts); // the final number of vertices

	rod.setMass(density);
	rod.setRadius(radius);
	rod.setGravityAcc(gravityAcc);

	rod.setEnterSpeed(enterSpeed); // speed of incoming material
	rod.setWallPosition(wallPosition);

	Scalar skinThickness{ 0 };

	Scalar bendingStiffness{ YModulus * PI * pow(radius, 4.0)  / 4.0 };
	Scalar stretchingStiffness{ YModulus * PI * pow(radius, 2.0) };
	Scalar gDragCoefficient{ muDragCoeff};
	rod.setCoefficients(bendingStiffness, stretchingStiffness, gDragCoefficient);
	rod.setRandomForceMag(bias);


	Scalar dt{ finalTime / numOfTimeSteps }; //time step size

	int saveRate{ numOfTimeSteps / numSaves}; // number of time steps per save

	int enterRate{ numOfTimeSteps / reducedNumVerts }; // number of time steps per addition of vertex. total will be equal to numOfVerts

	targetLength = enterRate * dt * enterSpeed;


	//std::cout << "Total simulation time: " << finalTime << std::endl;

	rod.startSaveFile(saveName);


	std::cout << "\nStarting Simulation ...\n" << std::endl;
	for (int index{ 0 }; index < numOfTimeSteps; ++index)
	{

		if (index % saveRate == 0)
		{
			std::cout << index / saveRate + 1 << " - " << "saving data at time = " << m_time << std::endl;

			rod.saveState(m_time);
		}

		if ((index % enterRate == 0))
		{
			if (rod.getNumVertices() >= maxNumOfVerts) {
				rod.removeVertex(density);
			}

			//std::cout  << "adding new vertex at time: " << m_time << std::endl;
			//Scalar relativeTime = (finalTime - m_time) / finalTime;

			//Scalar biasMag{ bias * pow(relativeTime, 3) };
			Scalar rotRate{ 1 };
			Scalar rx = bias * cos(2 * PI * rotRate * m_time) / (1 + m_time * m_time) + (2 * (static_cast <Scalar> (rand()) / static_cast <Scalar> (RAND_MAX)) - 0.5) * bias/ (1 + m_time * m_time);//
			Scalar ry = bias * sin(2 * PI * rotRate * m_time) / (1 + m_time * m_time) + (2 * (static_cast <Scalar> (rand()) / static_cast <Scalar> (RAND_MAX)) - 0.5) * bias/ (1 + m_time * m_time);//
			//rod.addVertex({ 0, bias*pow(relativeTime, 3.96), 0.96 * bias * pow(relativeTime, 4.1) } , { enterSpeed, 0, 0 }, density, targetLength);
			//rod.addVertex({ 0, biasMag * cos(10 * PI * relativeTime), biasMag * sin(10 * PI * relativeTime) }, {enterSpeed, 0, 0}, density, targetLength);
			rod.addVertex({0, rx, ry}, {enterSpeed, 0, 0}, density, targetLength);
			
			
		}

		rod.updatePositions(dt);
		rod.crossLink(dt * crossLinkRate);
		rod.setCoefficients(bendingStiffness, stretchingStiffness, gDragCoefficient);
		m_time += dt;

	}

}



void SimulationMaker::terminalVelocity(Scalar finalTime, int numOfTimeSteps, int numOfVerts, Scalar wallPosition, Scalar density, Scalar YModulus, Scalar muDragCoeff, 
									Scalar radius, Scalar length, Scalar gravityAcc, int numSaves, std::string saveName)
{
	

	//{8*targetLength, 0, 0},{7 * targetLength, 0 ,0},{6 * targetLength, 0, 0},{5 * targetLength, 0 ,0},{4 * targetLength, 0, 0}, { 3 * targetLength, 0 ,0 },
	Rod3D rod{ initializehLine(numOfVerts, length), Boundary::DEFAULT };



	rod.setMass(density);
	rod.setRadius(radius);
	rod.setWallPosition(wallPosition);
	rod.setGravityAcc(gravityAcc);


	Scalar bendingStiffness{ YModulus * PI * pow(radius, 4.0) / 4.0 };
	Scalar stretchingStiffness{ YModulus * PI * pow(radius, 2.0) };
	Scalar gDragCoefficient{ muDragCoeff };
	rod.setCoefficients(bendingStiffness, stretchingStiffness, gDragCoefficient);



	Scalar dt{ finalTime / numOfTimeSteps }; //time step size

	int saveRate{ numOfTimeSteps / numSaves }; // number of time steps per save

	
	//std::cout << "Total simulation time: " << finalTime << std::endl;

	rod.startSaveFile(saveName);


	std::cout << "\nStarting Simulation ...\n" << std::endl;
	for (int index{ 0 }; index < numOfTimeSteps; ++index)
	{

		if (index % saveRate == 0)
		{
			std::cout << index / saveRate + 1 << " - " << "saving data at time = " << m_time << std::endl;

			rod.saveState(m_time);
		}


		rod.updatePositions(dt);
		m_time += dt;

	}

}





void SimulationMaker::rodRelaxation(int numOfTimeSteps, Scalar totalTime, int numOfVerts, Scalar density, Scalar YModulus, Scalar muDragCoeff, Scalar radius, 
	Scalar wallPosition, int numSaves, Scalar bias, std::string saveName)
{
	std::vector<Rod3D::Vector3D> vertexPositions{};

	vertexPositions.assign(numOfVerts, { 0,0,0 });
	vertexPositions = initializeCircle(numOfVerts, bias);
	//vertexPositions = initializeLine(numOfVerts, bias);


	Rod3D rod{ vertexPositions, Boundary::CANTILEVER };


	rod.setMass(density);

	Scalar bendingStiffness{ YModulus * PI * pow(radius, 4.0) / 4.0 };
	Scalar stretchingStiffness{ YModulus * PI * pow(radius,2.0) };
	Scalar gDragCoefficient{ muDragCoeff};
	//rod.setCoefficients(bendingStiffness, stretchingStiffness, gDragCoefficient);
	rod.setCoefficients(1, 100, 0.1);
	rod.setWallPosition(wallPosition);

	
	rod.setTargetCurvatures(helixCurvatures(numOfVerts));

	Scalar dt{ totalTime / numOfTimeSteps };
	int saveRate{ numOfTimeSteps / numSaves }; // number of time steps per save\

	rod.startSaveFile(saveName);


	std::cout << "\nStarting Simulation ...\n" << std::endl;
	for (int index{ 0 }; index < numOfTimeSteps; ++index)
	{

		if (index % saveRate == 0)
		{
			std::cout << index / saveRate + 1 << " - " << "saving data at time = " << m_time << std::endl;

			rod.saveState(m_time);
			rod.testing();
		}


		rod.updatePositions(dt);
		m_time += dt;

	}
}




void SimulationMaker::EulerBuckling(Scalar preStrain,  int numOfTimeSteps, Scalar dt, int numOfVerts, Scalar density, Scalar YModulus, Scalar muDragCoeff,
										Scalar radius, int numSaves, Scalar bias)
{
	//Scalar bias{ 0.01 };

	//Rod3D rod{ EulerInitPositions(numOfVerts, bias), Boundary::EULER };

	Rod3D rod{numOfVerts, Boundary::EULER };

	rod.setMass(density);

	rod.addPrestrain(preStrain);

	rod.setPositions(EulerInitPositions(numOfVerts, bias));

	//define external force which is at the tip
	//std::vector<Rod3D::Vector3D> externalForces{ static_cast<unsigned int>(numOfVerts), { 0,0 } };

	//a uniform load pointing downwards to bias the buckling
	//rod.addGravity(bias);

	Scalar bendingStiffness{ YModulus * PI * pow(radius, 4.0) / 4.0 };
	Scalar stretchingStiffness{ YModulus * PI * pow(radius,2.0) };
	Scalar gDragCoefficient{ muDragCoeff};


	rod.setCoefficients(bendingStiffness, stretchingStiffness, gDragCoefficient);

	int saveRate{ numOfTimeSteps / numSaves };

	//std::cout << "saving data at time = " << m_time << std::endl;
	//rod.saveState();

	std::cout << "Starting Simulation \n";
	for (int index{ 0 }; index < numOfTimeSteps; ++index)
	{

		if (index % saveRate == 0)
		{
			std::cout << "saving data at time = " << m_time << std::endl;
			rod.saveState(m_time);
		}

		rod.updatePositions(dt);
		m_time += dt;

		
	}
}


std::vector<Rod3D::Vector3D> SimulationMaker::EulerInitPositions(int numOfVerts, Scalar bias)
{
	std::vector<Rod3D::Vector3D> positions{};

	for (int index{ 0 }; index < numOfVerts; ++index)
	{
		Scalar u{ static_cast<long double>(index) / (numOfVerts - 1.0) };

		positions.push_back({ u , bias * u * (u - 1) });
	}

	return positions;
}




void SimulationMaker::fluidBuckling(int numOfTimeSteps, Scalar dt, int numOfVerts, Scalar length, Scalar density, Scalar YModulus, Scalar muDragCoeff,
	Scalar radius, Rod3D::Vector3D fluidVelocity, Scalar bias, int numSaves,  std::string saveName)
{
	Rod3D rod{ fluidInitPositions(numOfVerts, length, bias), Boundary::CANTILEVER };

	//rod.setPositions(fluidInitPositions(numOfVerts, length, bias));

	rod.setMass(density);
	rod.setRadius(radius);
	rod.SetFluidVelocity(fluidVelocity);




	Scalar bendingStiffness{ YModulus * PI * pow(radius, 4.0) / 4.0 };
	Scalar stretchingStiffness{ YModulus * PI * pow(radius,2.0) };
	Scalar gDragCoefficient{ muDragCoeff };
	rod.setCoefficients(bendingStiffness, stretchingStiffness, gDragCoefficient);
	rod.setGravityAcc(0.0);
	rod.setWallPosition(10 * length);

	int saveRate{ numOfTimeSteps / numSaves };


	rod.startSaveFile(saveName);

	std::cout << "Total simulation time: " << numOfTimeSteps * dt << std::endl;

	std::cout << "\nStarting Simulation ...\n" << std::endl;
	for (int index{ 0 }; index < numOfTimeSteps; ++index)
	{

		if (index % saveRate == 0)
		{
			std::cout << index / saveRate + 1 << " - " << "saving data at time = " << m_time << std::endl;
			rod.saveState(m_time);
		}

		rod.updatePositions(dt);
		m_time += dt;


	}
}


std::vector<Rod3D::Vector3D> SimulationMaker::fluidInitPositions(int numOfVerts, Scalar length, Scalar bias)
{
	std::vector<Rod3D::Vector3D> positions{};

	for (int index{ 0 }; index < numOfVerts; ++index)
	{
		Scalar u{ (length * index) / (numOfVerts - 1.0) };

		positions.push_back({ u , bias * u * u, 0 });
	}

	return positions;
}

std::vector<Rod3D::Vector3D> SimulationMaker::initializeCircle(int numOfVerts, Scalar bias=0)
{
	
	std::vector<Rod3D::Vector3D> vertexPositions{};

	vertexPositions.assign(numOfVerts, { 0,0,0 });

	for (int index{ 0 }; index < numOfVerts; index++) {
		Scalar theta{ 2 * PI * index / (numOfVerts - 1.0) };
		vertexPositions[index] = { cos(theta), sin(theta), generateRandomScalar(bias) };
	}
	return vertexPositions;
}

std::vector<Rod3D::Matrix2D> SimulationMaker::helixCurvatures(int numOfVerts)
{

	std::vector<Rod3D::Vector3D> vertices{};
	vertices.assign(numOfVerts, Rod3D::Vector3D::Zero());


	Scalar torsion{ 1 }, radius{ 1 };

	for (int index{ 0 }; index < numOfVerts; index++) {
		Scalar theta{ 2 * PI * index / (numOfVerts - 1.0) };
		vertices[index] = { radius*cos(theta), radius*sin(theta), theta*torsion };
	}


	Rod3D rod{ vertices,  Boundary::CANTILEVER };


	std::vector<Rod3D::Vector3D> curvatures{ rod.getCurvatureBinormals() };
	std::vector<Rod3D::Matrix2D> Fcurvatures{ rod.getFrameCurvatures() };
	std::vector<Scalar> edgeLs{ rod.getEdgeLengths() };

	int indx{ 50 };
	Scalar vorLength = 0.5 * (edgeLs[indx - 1] + edgeLs[indx]);

	std::cout << curvatures[indx].norm()/vorLength << std::endl << std::endl;
	std::cout << Fcurvatures[indx] / vorLength << std::endl;
	return rod.getFrameCurvatures();
}

std::vector<Rod3D::Vector3D> SimulationMaker::initializeLine(int numOfVerts, Scalar bias)
{
	std::vector<Rod3D::Vector3D> vertexPositions{};

	vertexPositions.assign(numOfVerts, { 0,0,0 });

	for (int index{ 0 }; index < numOfVerts; index++) {
		Scalar x{ 2 * PI * index / (numOfVerts - 1.0) };
		vertexPositions[index] = { x, generateRandomScalar(bias), 0};
	}
	return vertexPositions;
}


std::vector<Rod3D::Vector3D> SimulationMaker::initializehLine(int numOfVerts, Scalar length)
{
	std::vector<Rod3D::Vector3D> vertexPositions{};

	vertexPositions.assign(numOfVerts, { 0,0,0 });

	for (int index{ 0 }; index < numOfVerts; index++) {
		Scalar y{ length * index / (numOfVerts - 1.0) };
		vertexPositions[index] = { 0, y, 0 };
	}
	return vertexPositions;
}

//generates a random scalar between -bias and bias
SimulationMaker::Scalar SimulationMaker::generateRandomScalar(Scalar bias) {

	return 2 * bias * ((Scalar(rand()) / (Scalar(RAND_MAX) + 1.0)) - 0 * 0.5);
}





