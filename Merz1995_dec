% The model of Merz (1995)
% In log deviations


//endogenous variables
var varrho r q G W C Y I S V M U Z K N P y_dev u_dev v_dev;

//predetermined variables
predetermined_variables K N;

//exogenous variables
varexo e_Z;

//parameters
parameters lambda alpha rho cs etas cv etav delta_prime psi beta sig_e nu_est;

//initialize parameters
lambda=0.4;
alpha=0.36; 
rho=0.95;
cs=0.00000000001; 
etas=1000000000;
cv=0.05; 
etav=1; 
delta=0.022;
mu=0.004;
delta_prime=1-(1-delta)*(exp(-mu));
psi=0.07; 
beta= 1/(1.04^(1/4));
nu_est=1.25;
sig_e=0.007;

//model equations
model;

//----------- Household -----------------------------------
// 1. Marginal Untility from Consumption
[varrho]
exp(varrho) = (1/exp(C));

// 2. Marginal Disutility from Work
[G]
exp(G)=(exp(N)^(1-1/nu_est))/(1-1/nu_est);

// 3. Consumption Euler equation
[C]
    exp(varrho)=beta*(exp(varrho(+1))*(exp(r(+1))+1-delta_prime));

// 4. Search intensity Euler equation
[S]
    exp(varrho)*cs=exp(P)*beta*(exp(varrho(+1))*(W(+1)+cs)-exp(G(+1))+(exp(varrho(+1))*cs/exp(P(+1)))*(1+psi-exp(P(+1))*S(+1)));
	 

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
    exp(W) = lambda*((1-alpha)*exp(Y)/exp(N)+cv*exp(V)/exp(U))+(1-lambda)*((exp(N)^(-1/nu_est))/exp(varrho)-cs);

// --------Equlibrium market clearing ----------------------
// 9. Aggregate Reource constraint
[I]
    exp(Y)=exp(C)+exp(I)+cs*exp(S)^etas*exp(U)+cv*exp(V)^etav;


//----------Declarations------------------------------------

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

// 14. Avg. labor productivity
[P]
    exp(P)=exp(Y-N);

// 15. Job match probability
[q]
    exp(q)=(1-lambda)*exp(M)/exp(V);

// 16. Labor augmenting factor
[Z]
    Z=rho*Z(-1)+e_Z;

       


//--------- Cental Planner - [commented out] --------------
//MC of searching = MB from working
//    (cv*etav*exp(V)^(etav-1)/(exp(C)*(1-lambda)*exp(M-V)))=beta*((1/exp(C(+1)))*((1-alpha)*exp(Y(+1)-N(+1))+cs*exp(S(+1))^etas)-exp(N(+1))^(1/nu_est)+(cv*etav*exp(V(+1))^(etav-1)/(exp(C(+1))*(1-lambda)*exp(M(+1)-V(+1))))*(1-psi-lambda*exp(M(+1))/exp(U(+1))));

//MC of searching = MC of posting vacancies
//    (cs*etas*exp(S)^(etas)*exp(U)/((lambda)))=(cv*etav*exp(V)^(etav)/((1-lambda)));

//Labor market tightness 
//    exp(theta)=exp(V-U);

//Budget constraint
//    exp(C)=exp(Y)-exp(I)-cs*exp(S)^etas*exp(U)-cv*exp(V)^etav;


// 16. Output deviation from St.St.
[y_dev]
    y_dev = Y-STEADY_STATE(Y) ;

// 17. Unmployment deviation from St.St.
[u_dev]
    u_dev = U-STEADY_STATE(U) ;

//18. Vacancy deviation from St.St.
[v_dev]
    v_dev = V-STEADY_STATE(V) ;
end;

//initial values for parameters
initval;
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
    P=1.29900453434794;
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

//observed model variables (names must coincide with data variables in datafile)
//varobs y_dev u_dev v_dev; 


//specify parameter to be estimated, including their priors
//estimated_params;
//lambda, 0.75,0.01,0.99, beta_pdf, 0.5, 0.1;
//psi, 0.0765, 0.001,1, beta_pdf, 0.07, 0.05;
//nu_est, 1.25, 0.001, 10, gamma_pdf, 1.25, 0.2;
//stderr u_dev, 0.25*0.115881787006766, gamma_pdf,0.25*0.115881787006766,0.01;
//stderr v_dev, 0.25*0.136135889718564, gamma_pdf,0.25*0.136135889718564,0.01;
//stderr e_Z, 0.007, 0, 0.5,inv_gamma_pdf,0.01, 0.5;
//rho,0.9534, 0.0001,.999, beta_pdf, 0.95, 0.005;
//end;


//estimate model
//estimation(mode_compute=1,lik_init=2,order=1,plot_priors=1,datafile=data_dev,mh_replic=10000,mh_nblocks=2,mh_drop=0.2,mh_jscale=1,mode_check); 

//stochastic simulation using the estimated parameters
stoch_simul(irf=120,order=1,hp_filter=1600, periods=2100);
