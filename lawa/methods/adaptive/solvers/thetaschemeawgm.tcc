namespace lawa {

template <typename Index, typename ThetaTimeStepSolver>
ThetaSchemeAWGM<Index,ThetaTimeStepSolver>::ThetaSchemeAWGM(ThetaTimeStepSolver &_timestep_solver)
: timestep_solver(_timestep_solver), theta(1.), timestep(0.1), numOfTimesteps(10), timestep_eps(0.1),
  maxiterations(1), init_cgtol(1e-2)
{

}

template <typename Index, typename ThetaTimeStepSolver>
void
ThetaSchemeAWGM<Index,ThetaTimeStepSolver>::setParameters(T _theta, T _timestep, int _numOfTimesteps,
                                                          T _timestep_eps, int _maxiterations,
                                                          T _init_cgtol)
{
    theta = _theta;
    timestep = _timestep;
    numOfTimesteps = _numOfTimesteps;
    timestep_eps = _timestep_eps;
    maxiterations = _maxiterations;
    init_cgtol = _init_cgtol;
    timestep_solver.Op.setThetaTimeStepParameters(theta, timestep);
    timestep_solver.Prec.setThetaTimeStepParameters(theta, timestep);
}

template <typename Index, typename ThetaTimeStepSolver>
void
ThetaSchemeAWGM<Index,ThetaTimeStepSolver>::applyPreconditioner
(Coefficients<Lexicographical,T,Index> &v)
{
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= timestep_solver.Prec[(*it).first];
    }
}

template <typename Index, typename ThetaTimeStepSolver>
void
ThetaSchemeAWGM<Index,ThetaTimeStepSolver>::applyInvPreconditioner
(Coefficients<Lexicographical,T,Index> &v)
{
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./timestep_solver.Prec[(*it).first];
    }
}

template <typename Index, typename ThetaTimeStepSolver>
void
ThetaSchemeAWGM<Index,ThetaTimeStepSolver>::solve(Coefficients<Lexicographical,T,Index> &u)
{
    size_t hms;
    hms = u.bucket_count();

    T discrete_timepoint = 0.;
    Coefficients<Lexicographical,T,Index> propagated_u_k(hms);

    int maxDoF = 0;
    T avDoF = 0.;

    for (int i=1; i<=numOfTimesteps; ++i) {
        discrete_timepoint += timestep;

        T current_theta = theta;
        if (i<=4 && theta!=1.) {
            current_theta = 1.;
        }
        timestep_solver.Op.setThetaTimeStepParameters(current_theta, timestep);
        timestep_solver.Prec.setThetaTimeStepParameters(current_theta, timestep);

        propagated_u_k.setToZero();
        if (current_theta!=1.) {
            timestep_solver.Op.evalA(u, propagated_u_k, "galerkin");
            propagated_u_k *= (current_theta-1.)*timestep;
        }
        timestep_solver.Op.evalM(u, propagated_u_k, "galerkin");
        std::cerr << "i = " << i << ", timestep = " << timestep << ", theta = " << current_theta << std::endl;
        timestep_solver.F.setThetaTimeStepParameters(current_theta,timestep,discrete_timepoint,&propagated_u_k);
        /*
        if (discrete_timepoint <= 0.1) {
            u.clear();
            getSparseGridVector(timestep_solver.basis, u, 7, 0.L);
        }
        else {
            maxiterations = 1;
            init_cgtol = 1e-6;
        }
        */
        this->applyInvPreconditioner(u);
        timestep_solver.cg_solve(u,timestep_eps,maxiterations,init_cgtol);
        maxDoF = std::max(maxDoF, (int)u.size());
        avDoF += u.size();
        this->applyPreconditioner(u);
    }
    avDoF /= numOfTimesteps;
    std::cerr << "Theta timestepping with theta = " << theta << " and timestep " << timestep
              << " finished with maxDof = " << maxDoF << " and avDof = " << avDoF << std::endl;
}

}   // namespace lawa

