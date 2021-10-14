#pragma once

//#include "Rod3D.h"

#include <Eigen/Core>

#include <string>

//#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>


template <class T>
class SavingSystem
{
	std::string m_fileName;

	std::ofstream m_outputFile{ m_fileName };


public:

	SavingSystem(std::string fileName)
	{
		startFile(fileName);
	}

	void saveArray(const std::vector<T> array, int arraySize)
	{
		m_outputFile.open(m_fileName, std::ofstream::out | std::ofstream::app);

		if (!m_outputFile)
		{
			std::cerr << "dang it! could not open file" << std::endl;
		}


		for (int index{ 0 }; index < arraySize; ++index)
		{

			assert(!std::isnan(array[index]));
			m_outputFile << array[index] << ", ";
		}

		m_outputFile << "\n";

		m_outputFile.close();
	}


	typedef long double Scalar;
	void saveElement(Scalar element)
	{
		m_outputFile.open(m_fileName, std::ofstream::out | std::ofstream::app);

		if (!m_outputFile)
		{
			std::cerr << "dang it! could not open file" << std::endl;
		}


		m_outputFile << element << ", ";

		m_outputFile << "\n";

		m_outputFile.close();
	}


	typedef Eigen::Matrix<long double, 3, 1> Vector3D;
	void saveVectors(const std::vector<Vector3D> array, int arraySize)
	{
		m_outputFile.open(m_fileName, std::ofstream::out | std::ofstream::app);

		if (!m_outputFile)
		{
			std::cerr << "dang it! could not open file" << std::endl;
		}


		for (int index{ 0 }; index < arraySize; ++index)
		{
			assert(!(std::isnan(array[index][0]) || std::isnan(array[index][1]) || std::isnan(array[index][2])));
			
			
			m_outputFile << array[index][0] << ", " << array[index][1] << ", " << array[index][2] << ", ";
		}

		m_outputFile << "\n";

		m_outputFile.close();
	}

	void startFile(std::string fileName) {
		m_fileName = fileName;
		std::ofstream m_outputFile{ m_fileName };

		m_outputFile.open(m_fileName, std::ofstream::out | std::ofstream::trunc);

		m_outputFile.close();
	}

};
