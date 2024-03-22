#ifndef Minimization_HPP
#define Minimization_HPP


#include<vector>
#include<string>
#include "GetPot"
#include "muParser.h"
#include <functional>
#include <cmath>

//enumerator class to choose step strategy
enum class StepSizeStrategy {
    ExponentialDecay,
    InverseDecay,
    ApproximateLineSearch,
    FixedAlpha
};


//Struct that contains all data
struct OptimizationParameters {
    std::string Method;// Chosen method to apply for the problem
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
    int dimension; //dimension of the problem
    bool usenumGrad; //finite difference gradient option
    double eta; //parameter for heavy ball,nesterov or Adam
    double epsilon; //parameter for Adam
    double beta_1; //parameter for Adam
    double beta_2; //parameter for Adam;
};

// Function to compute the norm of a vector
double vector_norm(const std::vector<double>& vec) {
    double norm = 0.0;
    for (size_t i=0; i<vec.size();++i) {
        norm += vec[i] * vec[i];
    }
    return std::sqrt(norm);
}

//Function to compute the elementwise product between vectors, used in Adam function 
std::vector<double> elementwise_prod(const std::vector<double>& vec1,const std::vector<double>& vec2) {
    std::vector<double> res(vec1.size());
    for (size_t i=0; i<vec1.size();++i) {
         res[i] = vec1[i] * vec2[i];
    }
    return res;
}

//scalar times vector produt
std::vector<double> product(const double s,const  std::vector<double>& v){
    std::vector<double> res(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        res[i] = v[i] * s;
    }
    return res;
}

//Vector division by a scalar
std::vector<double> div(const double s,const  std::vector<double>& v){
    std::vector<double> res(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        res[i] = v[i] / s;
    }
    return res;
}

//vector sum
std::vector<double> sum(const std::vector<double>& v1,const std::vector<double>& v2){
    std::vector<double> res(v1.size());
    for(size_t i=0; i<v1.size(); ++i){
        res[i]=(v1[i]+v2[i]);
    }
    return res;
}

//vector difference
std::vector<double> diff(const std::vector<double>& v1,const std::vector<double>& v2){
    std::vector<double> res(v1.size());
    for(size_t i=0; i<v1.size(); ++i){
        res[i]=(v1[i]-v2[i]);
    }
    return res;
}

//print vectorial result
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
    std::vector<double> x_values(variables); 
    for (size_t i = 0; i < variables.size(); ++i) {
        parser.DefineVar("x" + std::to_string(i + 1), &x_values[i]);
    }

    return parser.Eval();
}



// Armijo rule for step size selection
double armijo_rule(double alpha_0, const std::vector<double>& xk, const std::vector<double>& grad_fk, const std::function<double(const std::vector<double>&)>& f, double sigma) {
    double alpha = alpha_0;
    
    while (f(xk) - f(diff(xk,product(alpha,grad_fk))) < sigma * alpha * vector_norm(grad_fk) * vector_norm(grad_fk)) {
        alpha /= 2.0; 
    }
    return alpha;
}

//function to numerically evaluate the gradient by finite difference
std::vector<double> finiteDifferenceGradient(const std::string& f, const std::vector<double>& x, double epsilon = 1e-5) {
    int n = x.size();
    std::vector<double> gradient(n);

    for (int i = 0; i < n; ++i) {
        std::vector<double> x_plus_epsilon = x;
        x_plus_epsilon[i] += epsilon;
        gradient[i] = (evaluate_expression(f,x_plus_epsilon) - evaluate_expression(f,x)) / epsilon;
    }

    return gradient;
}

//main function thath brings together the gradient descent method
template<StepSizeStrategy Strategy>
std::vector<double> gradient_descent(const OptimizationParameters& params){
   
    //setting initial data
    std::vector<double> xk = params.initial_condition;
    std::vector<double> xk_n = params.initial_condition;
    double alpha = params.alpha_0;

    //initializing an empty vector for the gradient and a bool to check convergence
    std::vector<double> grad_fk(params.dimension);
    bool converged=false;

    std::cout << "Applying gradient descent method for " << params.expression_f;

    //using the Exact Graadient
    if(!params.usenumGrad){

        std::cout<<" using the exact gradient"<<std::endl;
        //looping till kmax if convergence is not reached
        for (int k = 0; k < params.max_iterations; ++k) {
            
            // Compute gradient of f at xk by evaluating it using the exact formula
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
                alpha = params.alpha_0 * std::exp(-params.mu * (k));
                } else if constexpr (Strategy == StepSizeStrategy::InverseDecay) {
                    alpha = params.alpha_0 / (1 + params.mu * (k));
                    } else {// Approximate line search using Armijo rule
                        alpha = armijo_rule(params.alpha_0, xk, grad_fk, [&](const std::vector<double>& variables) {
                        return evaluate_expression(params.expression_f, variables);}, params.sigma);
                        }
            
            xk=xk_n;
         }
     } //using the Finite difference Graadient
    else{

           std::cout<<" using the fd gradient"<<std::endl;
           
           //looping till k_max if convergence is not reached
           for (int k = 0; k < params.max_iterations; ++k) {
            // Compute gradient of f at xk by finite differences
            grad_fk=finiteDifferenceGradient(params.expression_f,xk);
            
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
                alpha = params.alpha_0 * std::exp(-params.mu * (k));
                } else if constexpr (Strategy == StepSizeStrategy::InverseDecay) {
                    alpha = params.alpha_0 / (1 + params.mu * (k));
                    } else {// Approximate line search using Armijo rule
                        alpha = armijo_rule(params.alpha_0, xk, grad_fk, [&](const std::vector<double>& variables) {
                        return evaluate_expression(params.expression_f, variables);}, params.sigma);
                        }
            
            xk=xk_n;
             }
        }
    
    if(converged==false)
    std::cout << "Method not converged, max_iterations reached"<<std::endl;
    return xk_n; // Return the final result
}

OptimizationParameters read_optimization_parameters(const std::string& filename) {
    GetPot config(filename.c_str()); // Create GetPot object with the filename
    OptimizationParameters params;

    // Read parameters from the configuration file
    params.Method = config("Method/choice","");
    params.expression_f = config("Function/Expression", "");
    params.dimension = config("Function/Dimension", 2);
    std::string gradstring;
    std::string initstring;
    for(int i=0; i<params.dimension; ++i){
    gradstring="Function/ExactGrad_"+std::to_string(i);
    initstring="Function/x_0_"+std::to_string(i);
    params.expression_grad_f.push_back(config(gradstring.c_str(),""));
    params.initial_condition.push_back(config(initstring.c_str(),0.));
    }

    //Reading specific parameters for Gradient descent
    if(params.Method=="Gradient"){
    params.epsilon_r = config("GradientParameters/EpsilonR", 1e-6);
    params.epsilon_s = config("GradientParameters/EpsilonS", 1e-6);
    params.alpha_0 = config("GradientParameters/Alpha0",0.05);
    params.max_iterations = config("GradientParameters/MaxIterations", 1000);
    params.mu = config("GradientParameters/Mu", 0.2);
    params.sigma = config("GradientParameters/Sigma", 0.1);
    params.usenumGrad= config("GradientParameters/UseNumericalGradient",false);

    std::string strategy_str = config("GradientParameters/StepSizeMethod", "ApproximateLineSearch");
    
    if (strategy_str == "ExponentialDecay")
        params.step_size_strategy = StepSizeStrategy::ExponentialDecay;
    else if (strategy_str == "InverseDecay")
        params.step_size_strategy = StepSizeStrategy::InverseDecay;
    else if (strategy_str == "ApproximateLineSearch")
        params.step_size_strategy = StepSizeStrategy::ApproximateLineSearch;
    else
        throw std::invalid_argument("Invalid step size strategy");    
    }
    //reading specific parameters for Heavyball
    else if(params.Method=="HeavyBall"){
    params.epsilon_r = config("HeavyballParameters/EpsilonR", 1e-6);
    params.epsilon_s = config("HeavyballParameters/EpsilonS", 1e-6);
    params.alpha_0 = config("HeavyballParameters/Alpha0",0.05);
    params.mu = config("HeavyballParameters/Mu", 0.2);
    params.max_iterations = config("HeavyballParameters/MaxIterations", 1000);
    params.eta = config("HeavyballParameters/Eta", 0.9);
    params.usenumGrad= config("HeavyballParameters/UseNumericalGradient",false);

    std::string strategy_str = config("HeavyballParameters/StepSizeMethod", "ExponentialDecay");
    
    if (strategy_str == "ExponentialDecay")
        params.step_size_strategy = StepSizeStrategy::ExponentialDecay;
    else if (strategy_str == "InverseDecay")
        params.step_size_strategy = StepSizeStrategy::InverseDecay;
    else if (strategy_str == "FixedAlpha")
        params.step_size_strategy = StepSizeStrategy::FixedAlpha;
    else
        throw std::invalid_argument("Invalid step size strategy");
    }
     //reading specific parameters for Nesterov
    else if(params.Method=="Nesterov"){
    params.epsilon_r = config("Nesterov/EpsilonR", 1e-6);
    params.epsilon_s = config("Nesterov/EpsilonS", 1e-6);
    params.alpha_0 = config("Nesterov/Alpha0",0.05);
    params.mu = config("Nesterov/Mu", 0.2);
    params.max_iterations = config("Nesterov/MaxIterations", 1000);
    params.eta = config("Nesterov/Eta", 0.9);
    params.usenumGrad= config("Nesterov/UseNumericalGradient",false);

    std::string strategy_str = config("Nesterov/StepSizeMethod", "ExponentialDecay");
    
    if (strategy_str == "ExponentialDecay")
        params.step_size_strategy = StepSizeStrategy::ExponentialDecay;
    else if (strategy_str == "InverseDecay")
        params.step_size_strategy = StepSizeStrategy::InverseDecay;
    else if (strategy_str == "FixedAlpha")
        params.step_size_strategy = StepSizeStrategy::FixedAlpha;
    else
        throw std::invalid_argument("Invalid step size strategy");
    }
     //reading specific parameters for Adam
    else if(params.Method=="Adam"){
    params.epsilon = config("Adam/Epsilon", 1e-8);
    params.epsilon_r = config("Adam/EpsilonR", 1e-6);
    params.epsilon_s = config("Adam/EpsilonS", 1e-6);
    params.eta = config("Adam/Eta", 0.01);
    params.beta_1 = config("Adam/Beta_1", 0.2);
    params.beta_2 = config("Adam/Beta_2", 0.2);
    params.max_iterations = config("Adam/MaxIterations", 1000);
    params.usenumGrad= config("Adam/UseNumericalGradient",false);
    }

    else throw std::invalid_argument("Invalid Method");

    return params;
}

//heavy ball method
template<StepSizeStrategy Strategy>
std::vector<double> heavy_ball(const OptimizationParameters& params){
   
    //setting initial data
    std::vector<double> xk = params.initial_condition;
    std::vector<double> xk_n(xk.size());
    double alpha = params.alpha_0;
    std::vector<double> grad_fk(xk.size());
    bool converged=false;
     std::vector<double> d(xk.size());

    std::cout << "Applying  heavy ball method for " << params.expression_f;

    //using the Exact Graadient
    if(!params.usenumGrad){

        std::cout<<" using the exact gradient"<<std::endl;
        
        //initializing the exact gradient and d
        for(int i=0; i< params.dimension; ++i){
                grad_fk[i]=evaluate_expression(params.expression_grad_f[i],xk);
                }
        d=product(-alpha,grad_fk);

        //looping till kmax if convergence is not reached
        for (int k = 0; k < params.max_iterations; ++k) {
            
            //update xk
            xk_n=sum(xk,d);
            
            //compute new gradient
             for(int i=0; i< params.dimension; ++i){
                grad_fk[i]=evaluate_expression(params.expression_grad_f[i],xk_n);
                }
            
            // Update step size alpha using the chosen strategy
             if constexpr (Strategy == StepSizeStrategy::ExponentialDecay) {
                alpha = params.alpha_0 * std::exp(-params.mu * (k));
                } else if constexpr (Strategy == StepSizeStrategy::InverseDecay) {
                    alpha = params.alpha_0 / (1 + params.mu * (k));
                    }

            //update d
            d=diff(product(params.eta,d),product(alpha,grad_fk));
            
            // Check convergence criteria
            if (vector_norm(diff(xk_n,xk)) < params.epsilon_s || std::abs(evaluate_expression(params.expression_f, xk_n) - evaluate_expression(params.expression_f, xk)) < params.epsilon_r) {
                 std::cout<< "Convergence achieved in " << k << " iterations" << std::endl;
                 converged=true;
                 break; // Convergence achieved
                 }
                 
            xk=xk_n;
         }
     } //using the Finite difference Graadient
    else{
         std::cout<<" using the fd gradient"<<std::endl;

         //initializing the fd gradient and d
         grad_fk=finiteDifferenceGradient(params.expression_f,xk);
         d=product(-alpha,grad_fk);

        //  //looping till kmax if convergence is not reached
        for (int k = 0; k < params.max_iterations; ++k) {


            // updating x and the gradient
            xk_n=sum(xk,d);
            grad_fk=finiteDifferenceGradient(params.expression_f,xk_n);
            
            
            // Update step size alpha using the chosen strategy
            if constexpr (Strategy == StepSizeStrategy::ExponentialDecay) {
                alpha = params.alpha_0 * std::exp(-params.mu * (k));
                } else if constexpr (Strategy == StepSizeStrategy::InverseDecay) {
                    alpha = params.alpha_0 / (1 + params.mu * (k));
                    } 

            //update d
            d=diff(product(params.eta,d),product(alpha,grad_fk));
            
            // Check convergence criteria
            if (vector_norm(diff(xk_n,xk)) < params.epsilon_s || std::abs(evaluate_expression(params.expression_f, xk_n) - evaluate_expression(params.expression_f, xk)) < params.epsilon_r) {
                 std::cout<< "Convergence achieved in " << k << " iterations" << std::endl;
                 converged=true;
                 break; // Convergence achieved
                 }
                 
            xk=xk_n;
         }
           
        }
    
    if(converged==false)
    std::cout << "Method not converged, max_iterations reached"<<std::endl;
    return xk_n; // Return the final result
}

template<StepSizeStrategy Strategy>
std::vector<double> Nesterov(const OptimizationParameters& params){
   
    //setting initial data
    std::vector<double> xk = params.initial_condition;
    std::vector<double> xk_n(xk.size());
    std::vector<double> y(xk.size());
    double alpha = params.alpha_0;
    std::vector<double> grad_fk(xk.size());
    bool converged=false;

    std::cout << "Applying  Nesterov method for " << params.expression_f;

    //using the Exact Graadient
    if(!params.usenumGrad){

        std::cout<<" using the exact gradient"<<std::endl;
        
        //initializing the gradient and xk
        for(int i=0; i< params.dimension; ++i){
                grad_fk[i]=evaluate_expression(params.expression_grad_f[i],xk);
                }
        xk_n=diff(xk,product(alpha,grad_fk));
        
        //looping till kmax if convergence is not reached
        for (int k = 0; k < params.max_iterations; ++k) {
            
            y=sum(xk_n,product(params.eta,diff(xk_n,xk)));

             for(int i=0; i< params.dimension; ++i){
                grad_fk[i]=evaluate_expression(params.expression_grad_f[i],y);
                }
            
            // Update step size alpha using the chosen strategy
            if constexpr (Strategy == StepSizeStrategy::ExponentialDecay) {
                alpha = params.alpha_0 * std::exp(-params.mu * (k));
                } else if constexpr (Strategy == StepSizeStrategy::InverseDecay) {
                    alpha = params.alpha_0 / (1 + params.mu * (k));
                    }

            //updating xk
            xk_n=diff(y,product(alpha,grad_fk));
            
            // Check convergence criteria
            if (vector_norm(diff(xk_n,xk)) < params.epsilon_s || std::abs(evaluate_expression(params.expression_f, xk_n) - evaluate_expression(params.expression_f, xk)) < params.epsilon_r) {
                 std::cout<< "Convergence achieved in " << k << " iterations" << std::endl;
                 converged=true;
                 break; // Convergence achieved
                 }
                 
            xk=xk_n;
         }
     } //using the Finite difference Graadient
    else{
         std::cout<<" using the fd gradient"<<std::endl;
        
        //Initializing gradient by fd and xk
        grad_fk=finiteDifferenceGradient(params.expression_f,xk);
        xk_n=diff(xk,product(alpha,grad_fk));

        //looping till kmax if convergence is not reached
        for (int k = 0; k < params.max_iterations; ++k) {
            
            y=sum(xk_n,product(params.eta,diff(xk_n,xk)));

            grad_fk=finiteDifferenceGradient(params.expression_f,y);
            
            // Update step size alpha using the chosen strategy
             if constexpr (Strategy == StepSizeStrategy::ExponentialDecay) {
                alpha = params.alpha_0 * std::exp(-params.mu * (k));
                } else if constexpr (Strategy == StepSizeStrategy::InverseDecay) {
                    alpha = params.alpha_0 / (1 + params.mu * (k));
                    }

            //updating xk
            xk_n=diff(y,product(alpha,grad_fk));
            
            // Check convergence criteria
            if (vector_norm(diff(xk_n,xk)) < params.epsilon_s || std::abs(evaluate_expression(params.expression_f, xk_n) - evaluate_expression(params.expression_f, xk)) < params.epsilon_r) {
                 std::cout<< "Convergence achieved in " << k << " iterations" << std::endl;
                 converged=true;
                 break; // Convergence achieved
                 }
                 
            xk=xk_n;
         }
           
        }
    
    if(converged==false)
    std::cout << "Method not converged, max_iterations reached"<<std::endl;
    return xk_n; // Return the final result
}


std::vector<double> Adam(const OptimizationParameters& params){
   
    //setting initial data
    std::vector<double> xk = params.initial_condition;
    std::vector<double> xk_n(xk.size());
    std::vector<double> m(xk.size(),0);
    std::vector<double> m_hat(xk.size(),0);
    std::vector<double> v(xk.size(),0);
    std::vector<double> v_hat(xk.size(),0);
    double beta_1 = params.beta_1;
    double beta_2 = params.beta_2;
    double beta_1_t = params.beta_1;
    double beta_2_t = params.beta_2;
    std::vector<double> grad_fk(xk.size());
    bool converged=false;

    std::cout << "Applying Adam method for " << params.expression_f;

    //using the Exact Graadient
    if(!params.usenumGrad){

        std::cout<<" using the exact gradient"<<std::endl;
        
    //looping till kmax if convergence is not reached
    for (int k = 0; k < params.max_iterations; ++k) {

        //computing the gradient
        for(int i=0; i< params.dimension; ++i){
                grad_fk[i]=evaluate_expression(params.expression_grad_f[i],xk);
                }
        
        //Computing m and v
        m=sum(product(beta_1,m),product((1-beta_1),grad_fk));
        v=sum(product(beta_2,v),product((1-beta_2),elementwise_prod(grad_fk,grad_fk)));

        //Computing m_hat and v_hat
        m_hat=div((1-beta_1_t),m);
        v_hat=div((1-beta_2_t),v);
        
        for (size_t i = 0; i < xk_n.size(); ++i) {
            
            // Update parameters
            xk_n[i] = xk[i] -( params.eta/ (std::sqrt(v_hat[i]) + params.epsilon)) * m_hat[i];
        }
            
            // Check convergence criteria
            if (vector_norm(diff(xk_n,xk)) < params.epsilon_s || std::abs(evaluate_expression(params.expression_f, xk_n) - evaluate_expression(params.expression_f, xk)) < params.epsilon_r) {
                 std::cout<< "Convergence achieved in " << k << " iterations" << std::endl;
                 converged=true;
                 break; // Convergence achieved
                 }

            //"elevating" betas
            beta_1_t*=beta_1;
            beta_2_t*=beta_2;     
            xk=xk_n;
         }
     } //using the Finite difference Graadient
    else{
          std::cout<<" using the fd gradient"<<std::endl;
        
    //looping till kmax if convergence is not reached
    for (int k = 0; k < params.max_iterations; ++k) {
     
        //computing gradient by df, m and v
        grad_fk=finiteDifferenceGradient(params.expression_f,xk);
        m=sum(product(beta_1,m),product((1-beta_1),grad_fk));
        v=sum(product(beta_2,v),product((1-beta_2),elementwise_prod(grad_fk,grad_fk)));

        //computing m_hat and v_hat
        m_hat=div((1-beta_1_t),m);
        v_hat=div((1-beta_2_t),v);
        
        for (size_t i = 0; i < xk_n.size(); ++i) {
            
            // Update parameters
            xk_n[i] = xk[i] -( params.eta/ (std::sqrt(v_hat[i]) + params.epsilon)) * m_hat[i];
        }
            
            // Check convergence criteria
            if (vector_norm(diff(xk_n,xk)) < params.epsilon_s || std::abs(evaluate_expression(params.expression_f, xk_n) - evaluate_expression(params.expression_f, xk)) < params.epsilon_r) {
                 std::cout<< "Convergence achieved in " << k << " iterations" << std::endl;
                 converged=true;
                 break; // Convergence achieved
                 }
                 
            //"elevating" betas
            beta_1_t*=beta_1;
            beta_2_t*=beta_2;
            xk=xk_n;
        }
    }
    if(converged==false)
    std::cout << "Method not converged, max_iterations reached"<<std::endl;
    return xk_n; // Return the final result
}

#endif