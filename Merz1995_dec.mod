/*
 * This file replicates the model studied in:
 * Merz (1995): "Search in labor market and the real business cycle", 
 * Journal of Monetary Economics, 36 (1995), pp. 269-300.
 * 
 * It provides a replication of the decentralized version of the model from the paper  
 *of the model from the paper
 *
 * This implementation was written by Artem Shramko. 
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2018 Artem Shramko
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
 */



 //endogenous variables
var varrho r q G W C Y I S V M U Z K N p theta y_dev u_dev v_dev;

 //predetermined variables
predetermined_variables K N;

 //exogenous variables
varexo e_Z;

 //parameters
parameters lambda alpha rho cs etas cv etav delta_prime psi beta sig_e nu_est;

 //initialize parameters
lambda=0.4;   % weight of search effort and unemployed in matching function
alpha=0.36;   % weight of capital in production function
rho=0.95;     % persistence of technology shock
cs=0.005;    % search costs (c0)
etas=10;      % eta in search cost function
cv=0.05;      % vacancy posting costs
etav=1;       % curvature of vacancy posting costs (=1 linear costs)
delta=0.022;  % depreciation rate of capital
mu=0.004;     % growth rate
delta_prime=1-(1-delta)*(exp(-mu));
psi=0.07;     % constant separation rate
beta= 1/(1.04^(1/4));  % discount factor
nu_est=-1.25;  % Frisch elasticity
sig_e=0.007;  % standard devaition of shock



 //model equations
 
model;
 //----------- Household -----------------------------------

// 1. Marginal Untility from Consumption
[varrho]
exp(varrho) = (1/exp(C));

 // 2. Marginal Disutility from Work
[G]
exp(G)=(exp(N)^(1/nu_est));

 // 3. Consumption Euler equation
[C]
    exp(varrho)=beta*(exp(varrho(+1))*(exp(r(+1))+1-delta_prime));

 // 4. Search intensity Euler equation
[S]
    exp(varrho)*cs=exp(p)*beta*(exp(varrho(+1))*(exp(W(+1))+cs)-exp(G(+1))+(exp(varrho(+1))*cs/exp(p(+1)))*(1+psi-exp(p(+1))*exp(S(+1))));
	 
 //------------- Firm --------------------------------------

// 5. Production function
[Y]
    exp(Y)=exp(Z)*exp(K)^alpha*(exp(N))^(1-alpha);
    
 // 6. Interest rate
[r]
    exp(r)=alpha*exp(Y)/exp(K);
    
 // 7. Vacancy creation Euler equation [cv - vacancy posting costs]
[V]
    exp(varrho)*cv=exp(q)*beta*exp(varrho(+1))*((1-alpha)*exp(Y(+1))/exp(N(+1))-exp(W(+1))+(cv/exp(q(+1)))*(1-psi));
    
 // 8. Wage
[W]
    exp(W) = lambda*((1-alpha)*exp(Y)/exp(N)+cv*exp(V)/exp(U))+(1-lambda)*(0.8*exp(G-varrho)-cs);
    
 // --------Equlibrium market clearing ----------------------
 
// 9. Aggregate Reource constraint
[I]
    exp(Y)=exp(C)+exp(I)+cs*exp(S)^etas*exp(U)+cv*exp(V)^etav;
    
 //----------Simplifictions------------------------------------
 
 // 10. Unemployment
[U]
    exp(U)=1-exp(N);
    
 // 11. Matching function
[M]
    exp(M)=exp(V)^(1-lambda)*(exp(S+U))^lambda;
    
 // 12. Kapital dynamics
[K]
    exp(K(+1))=(1-delta_prime)*exp(K)+exp(I);
    
 // 13. Labor dynamics
[N]
    exp(N(+1))=(1-psi)*exp(N)+exp(M);

 // 14. Probability to find job
[P]
    //exp(p)=exp(theta)^(1-lambda);
    
    exp(p)=exp(M-S-U);
    
 // 15. Job match probability
[q]
    exp(q)=(1-lambda)*exp(M-V);
    
 // 16. Labor augmenting factor
[Z]
    Z=rho*Z(-1)+e_Z;
        

 // 17. Output deviation from St.St.
[y_dev]
    y_dev = Y-STEADY_STATE(Y) ;
 // 18. Unmployment deviation from St.St.
[u_dev]
    u_dev = U-STEADY_STATE(U) ;
 //19. Vacancy deviation from St.St.
[v_dev]
    v_dev = V-STEADY_STATE(V) ;
end;

 //initial values for parameters
initval;

    W=0.8513;
    G=-0.0688;
    varrho=-0.9097;
    r=-3.3298;
    C=0.909695877517769;
    Y=1.21297419220172; 
    I=-0.131035097064121; 
    S=7.88090800314440e-10; 
    V=-2.91164322016667; 
    M=-2.74529037907900;
    U=-2.49576111823558;
    Z=0;
    K=3.52231558659806;
    N=-0.0860303421462215;
    //P=1.29900453434794;
    p=-0.2495;
    q=-0.3445;
    theta=-0.415882101931087;
    y_dev=0;
    u_dev=0;
    v_dev=0;
end;

check;

 //calculate steady state
steady;

 //variance of shocks
shocks;
var e_Z=sig_e^2;
end;

stoch_simul(irf=120,order=1,hp_filter=1600, periods=2100);
