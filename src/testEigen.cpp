#include <iostream>
#include "../Eigen/Dense"
#include "../Eigen/LU" 

using Eigen::MatrixXd;

using namespace Eigen;
using namespace std;


int main(){

	cout << "states * Wout = desiredOutputs\n" << endl;
	
	MatrixXd states(2,3);
	states(0,0) = 1;
	states(0,1) = 2;
	states(0,2) = 1;
	states(1,0) = 5;
	states(1,1) = 1;
	states(1,2) = 3;

	MatrixXd desiredOutputs(2,1);
	desiredOutputs(0,0) = 5;
	desiredOutputs(1,0) = -4;
	
	cout << "\nHere is the matrix states:\n" << states << endl;
	cout << "\nHere is the matrix desiredOutputs:\n" << desiredOutputs << endl;

	MatrixXd Wout = states.colPivHouseholderQr().solve(desiredOutputs);

	cout << "\nThe solution (Wout) is:\n" << Wout << endl;
	
	cout << "\nTest states * Wout:\n" <<  states * Wout << endl;

/*
	MatrixXd R = states.transpose() * states;
	MatrixXd P = states.transpose() * desiredOutputs;
	MatrixXd Wout = (R.inverse() * P).transpose();


	std::cout <<  R << std::endl;
	std::cout << "\n\n" << std::endl;
	std::cout <<  R.inverse() << std::endl;
*/
	

}
