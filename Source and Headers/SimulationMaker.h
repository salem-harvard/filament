#pragma once

#include "SavingSystem.h"
#include "Rod3D.h"


#include <Eigen/Core>



#include <vector>


class SimulationMaker
{
	typedef long double Scalar;
	Scalar PI{ EIGEN_PI };

public:

	/*enum class SimulationTypes
	{
		StraightRodTipForce = 1
	};*/

	//SimulationMaker(SimulationTypes whichCase, Scalar params);

	void EulerBuckling(Scalar preStrain, int numOfTimeSteps, Scalar dt, int numOfVerts, Scalar density, Scalar YModulus, Scalar muDragCoeff,
		Scalar radius, int numSaves, Scalar bias);


	void fluidBuckling(int numOfTimeSteps, Scalar dt, int numOfVerts, Scalar length, Scalar density, Scalar YModulus, Scalar muDragCoeff,
		Scalar radius, Rod3D::Vector3D fluidVelocity, Scalar bias, int numSaves, std::string saveName);

	void flowingIntoFluid(Scalar finalTime, int numOfTimeSteps, int vertsPerLength, int maxNumOfVerts, Scalar enterSpeed, Scalar crossLinkRate, Scalar wallPosition, Scalar density, Scalar YModulus,
		Scalar muDragCoeff, Scalar radius, Scalar gravityAcc, int numSaves, Scalar bias, std::string saveName);

	void terminalVelocity(Scalar finalTime, int numOfTimeSteps, int numOfVerts, Scalar wallPosition, Scalar density, Scalar YModulus, 
							Scalar muDragCoeff, Scalar radius, Scalar length, Scalar gravityAcc, int numSaves,std::string saveName);

	void rodRelaxation(int numOfTimeSteps, Scalar totalTime, int numOfVerts, Scalar density, Scalar YModulus, Scalar muDragCoeff,
		Scalar radius, Scalar wallPosition, int numSaves, Scalar bias, std::string saveName);


private:

	std::vector<Rod3D::Vector3D> EulerInitPositions(int numOfVerts, Scalar bias);

	std::vector<Rod3D::Vector3D> fluidInitPositions(int numOfVerts, Scalar length, Scalar bias);

	std::vector<Rod3D::Vector3D> initializeCircle(int numOfVerts, Scalar bias);

	std::vector<Rod3D::Matrix2D> helixCurvatures(int numOfVerts);

	std::vector<Rod3D::Vector3D> initializeLine(int numOfVerts, Scalar bias);

	std::vector<Rod3D::Vector3D> initializehLine(int numOfVerts, Scalar length=1.0);

	Scalar generateRandomScalar(Scalar bias);

	Scalar m_time{ 0 };
	
};
