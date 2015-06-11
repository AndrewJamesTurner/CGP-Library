// http://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

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

	// error in the inversion 
	double relative_error = (states*Wout - desiredOutputs).norm() / desiredOutputs.norm();
	cout << "The relative error is:\n" << relative_error << endl;


	// solve using least squares
	Wout = states.jacobiSvd(ComputeThinU | ComputeThinV).solve(desiredOutputs);
	cout << "\nThe solution (Wout) using least squares is:\n" << Wout << endl;

	// error in the inversion 
	relative_error = (states*Wout - desiredOutputs).norm() / desiredOutputs.norm();
	cout << "The relative error is:\n" << relative_error << endl;

}
