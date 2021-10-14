#include "SavingSystem.h"
#include "Rod3D.h"
#include "SimulationMaker.h"


#include <Eigen/Core>

#include <string>

//#include <cmath>
//#include <fstream>
#include <iostream>

#include <vector>
//using namespace Eigen;

void testing();
void rodRelaxation();
void flowingIntoFluid();
void terminalVelocity();
void fluidBuckling();

//main fro testing, the commented out one below is for the simulation
int main()
{
	flowingIntoFluid();
}



void flowingIntoFluid()
{
	typedef long double Scalar;
	Scalar PI{ EIGEN_PI };

	std::string saveName{ "-0.4mm-1mLperMin-50Pasc-20s-5g-75cm-1secInv-1dr" }; //"-0.4mm-6mLperMin-2700Pasc-4s-3g-2s^-1"

	Scalar radius{ 0.4 * pow(10.,-3) }; //
	Scalar volumeFlowRate{ 1 * pow(10.,-6.) / 60. }; //convert to M^3/second
	Scalar finalTime{ 20 };

	
	int numOfTimeSteps{ 500000 };  int vertsPerLength{ 500 }; int maxNumOfVerts{ 250 };  Scalar bias{ 0.001 }; int numSaves{ 300 };   Scalar wallPosition{ 0.75 };

	Scalar YModulus{ 2.7 * pow(10.,3) }; Scalar muDragCoeff{ 0.001 }; Scalar volumeDensity{ 1054. }; Scalar area{ PI * pow(radius, 2.0) };
	Scalar dragRatio{ 1 }; Scalar muEffective{ dragRatio * muDragCoeff };
	Scalar density{ volumeDensity * area };
	Scalar enterSpeed{ volumeFlowRate / area };
	Scalar crossLinkRate{ 1.0 };

	//effective gravitational acceleration (ratio of boyant to inertial mass times g = 9.8)
	Scalar g{9.8}; Scalar densityDifference{ 5. };  Scalar gravityAcc{ g * densityDifference / volumeDensity };

	SimulationMaker sim;
	

	sim.flowingIntoFluid(finalTime, numOfTimeSteps, vertsPerLength, maxNumOfVerts, enterSpeed, crossLinkRate, wallPosition, density, 
							YModulus, muEffective, radius, gravityAcc, numSaves, bias, saveName);
}

void terminalVelocity()
{
	typedef long double Scalar;
	Scalar PI{ EIGEN_PI };

	std::string saveName{ "-falling-rod" }; //"-0.5mm-10mLperMin-4800Pasc-10s-0g"

	Scalar radius{ 1.8 * pow(10.,-3) }; //
	
	Scalar finalTime{ 30 };

	//500000000
	int numOfTimeSteps{ 50000 };  int numOfVerts{ 10 };  int numSaves{ 100 };   Scalar wallPosition{ 0.5 };

	Scalar YModulus{ 2700 }; Scalar muDragCoeff{ 0.001 }; Scalar volumeDensity{ 1054. }; Scalar area{ PI * pow(radius, 2.0) };
	Scalar density{ volumeDensity * area }; Scalar length{ 0.02 }; 

	Scalar densityDifference{ 0 };  Scalar gravityAcc{ densityDifference / volumeDensity }; //effective gravitational acceleration (ratio of boyant to inertial mass times g = 9.8)


	SimulationMaker sim;

	sim.terminalVelocity(finalTime, numOfTimeSteps, numOfVerts, wallPosition, density, YModulus, muDragCoeff, radius, length, gravityAcc, numSaves, saveName);

}

void fluidBuckling()
{
	typedef long double Scalar;
	Scalar PI{ EIGEN_PI };
	typedef Eigen::Matrix<long double, 3, 1> Vector3D;

	std::string saveName{ "-fluidBuckling" };

	// parameters for fluid buckling
	int numOfTimeSteps{ 1000000 }; Scalar dt{ 0.02}; int numOfVerts{ 50 }; int numSaves{ 100 };  Scalar bias{ 0.001 }; Scalar density{ 10.0 };
	Scalar muDragCoeff{ 0.001 }; Scalar radius{ 0.5 * pow(10.,-3) }; Scalar length = 0.0245;  Scalar YModulus{ 2700 }; 	Vector3D fluidVelocity{ -0.02, 0, 0 };


	



	//std::cout << enterSpeed << std::endl;

	SimulationMaker sim;
	

	sim.fluidBuckling(numOfTimeSteps, dt, numOfVerts, length, density, YModulus, muDragCoeff, radius, fluidVelocity, bias, numSaves, saveName);

}

void rodRelaxation() {
	typedef long double Scalar;
	Scalar PI{ EIGEN_PI };
	typedef Eigen::Matrix<long double, 3, 1> Vector3D;

	Scalar radius{ 0.1 }; Scalar dt{ 0.005 };

	int numOfTimeSteps{ 850000 }, numOfVerts{ 100 }, numSaves{ 400 };

	Scalar totalTime{ dt * numOfTimeSteps };

	Scalar YModulus{ 480000 }; Scalar muDragCoeff{ 0.001 }; Scalar volumeDensity{ 1054. }; Scalar area{ PI * pow(radius, 2.0) }; Scalar density{ 10 };//{ volumeDensity * area };

	Scalar wallPosition{ 100 };

	SimulationMaker sim;

	std::string saveName{ "-relaxation-ToHelix3" };

	Scalar bias{ 0.0001 };
	sim.rodRelaxation(numOfTimeSteps, totalTime, numOfVerts, density, YModulus, muDragCoeff, radius, wallPosition, numSaves, bias, saveName);
}


/*int main()
{

	typedef long double Scalar;
	Scalar PI{ EIGEN_PI };

	std::string saveName{ "-0.5mm-10mLperMin-4800Pasc-10s-0g" };

	Scalar radius{ 0.5 * pow(10.,-3) }; //
	Scalar volumeFlowRate{10 * pow(10.,-6.) / 60. }; //convert to M^3/second
	//Scalar enterSpeed{ 30 * pow(10.,-3) };
	Scalar finalTime{ 10 };


	int numOfTimeSteps{ 500000000 };  int vertsPerLength{ 400 }; int maxNumOfVerts{ 1000 };  Scalar bias{ 0.00001 }; int numSaves{ 400 };   Scalar wallPosition{ 100 };
	
	Scalar YModulus{ 4800 }; Scalar muDragCoeff{ 0.001 }; Scalar volumeDensity{ 1054. }; Scalar area{ PI * pow(radius, 2.0) }; Scalar density{ volumeDensity * area };  
	 Scalar enterSpeed{ volumeFlowRate / area };



	// when dropping a rope without drag
	//Scalar YModulus{ 8.5 * pow(10.,6) }; Scalar muDragCoeff{ 0.000 }; Scalar density{ 2.9 * pow(10.,-3) }; 
	//Scalar enterSpeed{ 0.03 };  radius = 0.375 * pow(10., -3);



	//std::cout << enterSpeed << std::endl;

	SimulationMaker sim;
	//sim.EulerBuckling(prestrain, numOfTimeSteps,dt, numOfVerts, density, YModulus, muDragCoeff, radius, numSaves, bias);

	//sim.fluidBuckling(numOfTimeSteps, dt, numOfVerts, length, density, YModulus, muDragCoeff, radius, numSaves, bias, fluidVelocity);

	sim.flowingIntoFluid(finalTime, numOfTimeSteps, vertsPerLength, maxNumOfVerts, enterSpeed, wallPosition, density, YModulus, muDragCoeff, radius, numSaves, bias, saveName);





	// parameters for fluid buckling
	//int numOfTimeSteps{ 100000000 }; Scalar dt{ 0.0001 }; int numOfVerts{ 300 }; Scalar length = 6; Scalar density{ 50.0 }; Scalar YModulus{ 50 }; Scalar muDragCoeff{ 0.1 };
	//Scalar radius{ 0.03 }; int numSaves{ 100 };  Scalar bias{0.00000001 };
	//Rod2D::Vector2D fluidVelocity{ -0.002, 0 };

	//for euler buckling
	//Scalar prestrain{ 1.007 };
}*/


void testing() {
	typedef long double Scalar;
	Scalar PI{ EIGEN_PI };
	typedef Eigen::Matrix<long double, 3, 1> Vector3D;


	std::vector<Vector3D> vertexPositions{};
	int numOfVertices{ 1000 };

	vertexPositions.assign(numOfVertices, { 0,0,0 });

	for (int index{ 0 }; index < numOfVertices; index++) {
		Scalar theta{ 2 * PI * index / (numOfVertices - 1.0) };
		vertexPositions[index] = { cos(theta), sin(theta), 0 };
	}

	Rod3D rod{ vertexPositions, Boundary::CANTILEVER };
	rod.setWallPosition(100);
	rod.setMass(1);
	rod.setCoefficients(1, 1, 1);

	//Rod3D::Matrix2D barCurvature{}; barCurvature << 1, 1,
	//	0, 0;
	//rod.setTargetCurvatures(barCurvature);

	rod.testing();

	/*
	rod.startSaveFile("-rod-relaxation");

	int numOfTimeSteps{ 1000 }, numSaves{100};
	Scalar dt{ 0.0001 };

	int saveRate{ numOfTimeSteps / numSaves }; // number of time steps per save

	for (int tIndex{ 0 }; tIndex < numOfTimeSteps; ++tIndex) {

		if (tIndex % saveRate == 0)
		{
			//std::cout << tIndex / saveRate + 1 << " - " << "saving data at time = " << m_time << std::endl;

			rod.saveState();
		}

		rod.updatePositions(dt);

	}
	*/

}