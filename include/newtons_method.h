#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value


// f(x) is the E(x) in the slide
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd& x0, Objective& f, Jacobian& g, Hessian& H, unsigned int maxSteps, Eigen::VectorXd& tmp_g, Eigen::SparseMatrixd& tmp_H) {
	for (int i = 0; i < maxSteps; i++) {
		// gradient
		g(tmp_g, x0);
		// hessian
		H(tmp_H, x0);

		// check for convergence
		// we choose tol to be 1e-12
		// here means our first guess is good, exit loop
		double tol = 1e-10;
		if (tmp_g.norm() < tol) {
			return 0.0;
		}
		

		Eigen::SimplicialLDLT < Eigen::SparseMatrix < double > > solver;

		//Eigen::SparseMatrixd A;
		//Eigen::VectorXd b;


		// lecture 5 slide 25 and 27
		// A = tmp_H
		// b = - tmp_g
		
		solver.compute(tmp_H);
		Eigen::VectorXd search_direction;
		search_direction = solver.solve(-1.0 * tmp_g);

		// now add the line search 
		double a = 1.0;
		double c = 1e-8;
		double p = 0.5;
		// lecture 5 slide 39
		while (1) {

			if (a < tol || f(x0 + (a * search_direction)) <= (f(x0) + c * (search_direction.transpose() * tmp_g)(0)) ) {
				break;
			}
			a = a * p;
		}
		
		// Use search direction to update current guess
		x0 = x0 + a * search_direction;
	}


	return 0.0;
}
