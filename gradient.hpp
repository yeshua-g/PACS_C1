#include<vector>
#include<string>
#include "GetPot"
#include "muParser.h"
#include <functional>
#include <cmath>


enum class StepSizeStrategy {
    ExponentialDecay,
    InverseDecay,
    ApproximateLineSearch
};

struct OptimizationParameters {
    std::string expression_f; // Mathematical expression to minimize
    std::vector<std::string> expression_grad_f; // Mathematical expressions for gradient components
    std::vector<double> initial_condition; // Initial guess
    double epsilon_r; // Residual tolerance
    double epsilon_s; // Step length tolerance
    double alpha_0; // Initial step size
    int max_iterations; // Max number of iterations
    double mu; // Decay parameter for step size strategies
    double sigma; // Armijo rule parameter
    StepSizeStrategy step_size_strategy; // Step size strategy
    int dimension;
};

// Function to compute the norm of a vector
double vector_norm(const std::vector<double>& vec) {
    double norm = 0.0;
    for (size_t i=0; i<vec.size();++i) {
        norm += vec[i] * vec[i];
    }
    return std::sqrt(norm);
}

//scalar times vector produt
std::vector<double> product(const double s,const  std::vector<double>& v){
    std::vector<double> res(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        res[i] = v[i] * s;
    }
    return res;
}

//vector difference
std::vector<double> diff(const std::vector<double>& v1,const std::vector<double>& v2){
    std::vector<double> res(v1.size());
    for(auto i=0; i<v1.size(); ++i){
        res[i]=(v1[i]-v2[i]);
    }
    return res;
}

//print vector
void print(const std::vector<double>& v){
    std::cout<<"The computed minimum is:"<<std::endl;
    std::cout<<" [ " ;
    for(auto elem:v){
        std::cout<< elem << "  ";
    }
    std::cout<<"] "<<std::endl ;
    return;
}

// Function to evaluate mathematical expressions using muparser
double evaluate_expression(const std::string& expression, const std::vector<double>& variables) {
    mu::Parser parser;
    parser.SetExpr(expression);
    
    // Define variables
    std::vector<double> x_values(variables); // Make a copy to ensure they remain valid
    for (size_t i = 0; i < variables.size(); ++i) {
        parser.DefineVar("x" + std::to_string(i + 1), &x_values[i]);
    }

    return parser.Eval();
}



// Armijo rule for step size selection
double armijo_rule(double alpha_0, const std::vector<double>& xk, const std::vector<double>& grad_fk, const std::function<double(const std::vector<double>&)>& f, double sigma) {
    double alpha = alpha_0;
    
    while (f(xk) - f(diff(xk,product(alpha,grad_fk))) < sigma * alpha * vector_norm(grad_fk) * vector_norm(grad_fk)) {
        alpha /= 2.0; // Reduce alpha if Armijo condition is not satisfied
    }
    return alpha;
}

template<StepSizeStrategy Strategy>
std::vector<double> gradient_descent(const OptimizationParameters& params){
   
    std::vector<double> xk = params.initial_condition;
    std::vector<double> xk_n = params.initial_condition;
    double alpha = params.alpha_0;
    std::vector<double> grad_fk(params.dimension);
    bool converged=false;

    std::cout << "Applying gradient descent method for " << params.expression_f << std::endl;

    for (int k = 0; k < params.max_iterations; ++k) {

        // Compute gradient of f at xk
        for(int i=0; i< params.dimension; ++i){
            grad_fk[i]=evaluate_expression(params.expression_grad_f[i],xk);
        }

        // Update xk using gradient descent
        xk_n=diff(xk,product(alpha,grad_fk));

        // Check convergence criteria
        if (vector_norm(diff(xk_n,xk)) < params.epsilon_s || std::abs(evaluate_expression(params.expression_f, xk_n) - evaluate_expression(params.expression_f, xk)) < params.epsilon_r) {
            std::cout<< "Convergence achieved in " << k << " iterations" << std::endl;
            converged=true;
            break; // Convergence achieved
        }

        // Update step size alpha using the chosen strategy
        if constexpr (Strategy == StepSizeStrategy::ExponentialDecay) {
            alpha = params.alpha_0 * std::exp(-params.mu * (k+1));
        } else if constexpr (Strategy == StepSizeStrategy::InverseDecay) {
            alpha = params.alpha_0 / (1 + params.mu * (k+1));
        } else {
            // Approximate line search using Armijo rule
            alpha = armijo_rule(params.alpha_0, xk, grad_fk, [&](const std::vector<double>& variables) {
                return evaluate_expression(params.expression_f, variables);
            }, params.sigma);
        }
        xk=xk_n;
    }
    if(converged==false)
    std::cout << "Method not converged, max_iterations reached"<<std::endl;
    return xk; // Return the final result
}

OptimizationParameters read_optimization_parameters(const std::string& filename) {
    GetPot config(filename.c_str()); // Create GetPot object with the filename
    OptimizationParameters params;

    // Read parameters from the configuration file
    params.expression_f = config("Function/Expression", "");
    params.dimension = config("Function/Dimension", 2);
    std::string gradstring;
    std::string initstring;
    for(int i=0; i<params.dimension; ++i){
    gradstring="Function/ExactGrad"+std::to_string(i);
    initstring="Function/x_0"+std::to_string(i);
    params.expression_grad_f.push_back(config(gradstring.c_str(),""));
    params.initial_condition.push_back(config(initstring.c_str(),0.));
    }
    params.epsilon_r = config("GradientParameters/EpsilonR", 1e-6);
    params.epsilon_s = config("GradientParameters/EpsilonS", 1e-6);
    params.alpha_0 = config("GradientParameters/Alpha0",0.15);
    params.max_iterations = config("GradientParameters/MaxIterations", 1000);
    params.mu = config("GradientParameters/Mu", 0.2);
    params.sigma = config("GradientParameters/Sigma", 0.1);

    std::string strategy_str = config("GradientParameters/StepSizeMethod", "ExponentialDecay");
    /*std::cout<<strategy_str<<std::endl;*/
    if (strategy_str == "ExponentialDecay")
        params.step_size_strategy = StepSizeStrategy::ExponentialDecay;
    else if (strategy_str == "InverseDecay")
        params.step_size_strategy = StepSizeStrategy::InverseDecay;
    else if (strategy_str == "ApproximateLineSearch")
        params.step_size_strategy = StepSizeStrategy::ApproximateLineSearch;
    else
        throw std::invalid_argument("Invalid step size strategy");

    return params;
}


