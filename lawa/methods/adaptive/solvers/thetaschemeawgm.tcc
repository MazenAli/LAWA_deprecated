namespace lawa {

template <typename Index, typename ThetaTimeStepSolver>
ThetaSchemeAWGM<Index,ThetaTimeStepSolver>::ThetaSchemeAWGM(ThetaTimeStepSolver &_timestep_solver)
: timestep_solver(_timestep_solver), theta(1.), timestep(0.1), numOfTimesteps(10), timestep_eps(0.1),
  maxiterations(1), init_cgtol(1e-2), strategy(0)
{

}

template <typename Index, typename ThetaTimeStepSolver>
void
ThetaSchemeAWGM<Index,ThetaTimeStepSolver>::setParameters(T _theta, T _timestep, int _numOfTimesteps,
                                                          T _timestep_eps, int _maxiterations,
                                                          T _init_cgtol, int _strategy)
{
    theta = _theta;
    timestep = _timestep;
    numOfTimesteps = _numOfTimesteps;
    timestep_eps = _timestep_eps;
    maxiterations = _maxiterations;
    init_cgtol = _init_cgtol;
    strategy   = _strategy;

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
ThetaSchemeAWGM<Index,ThetaTimeStepSolver>::solve(Coefficients<Lexicographical,T,Index> &u,
                                                  int &avDof, int &maxDof, int &terminalDof)
{
    size_t hms;
    hms = u.bucket_count();

    T discrete_timepoint = 0.;
    Coefficients<Lexicographical,T,Index> propagated_u_k(hms), u0(hms), tmp(hms), tmp2(hms);
    u0 = u;
    T target_residual = 0.001*timestep_eps;
    bool strategy1_applied = false, strategy2_applied = false;


    if (strategy==2) {
        tmp.clear();
        tmp2.clear();
        tmp  = u;
        tmp2 = u;
        //T threshold = exp(3*discrete_timepoint)*0.001*timestep_eps*u.norm(2.);
        //T threshold = 0.000001*timestep_eps*u.norm(2.);
        T threshold = 1e-6;
        std::cerr << "   Size of u before threshold: " << u.size() << std::endl;
        std::cerr << "   Threshold is: " << threshold << std::endl;
        int maxLevelSum = 0;
        bool erasedNode = false;
        for (const_coeff_it it=tmp.begin(); it!=tmp.end(); ++it) {
           maxLevelSum = std::max(maxLevelSum, (*it).first.index1.j+(*it).first.index2.j);
           if ((fabs((*it).second)<threshold)) {
               if ((*it).first.index1.xtype==XBSpline || (*it).first.index2.xtype==XBSpline)
                       continue;
               else {
                   std::cerr << "   Erasing " << (*it).first << std::endl;
                   u.erase((*it).first);
                   erasedNode = true;
               }
           }
        }
        std::cerr << "   Max level sum = " << maxLevelSum << std::endl;
        tmp.clear();
        if (erasedNode) {
            tmp = u;
            for (coeff_it it=tmp.begin(); it!=tmp.end(); ++it) {
                completeMultiTree(timestep_solver.basis, (*it).first, u, 0, true, false);
            }
            tmp.clear();
        }
        for (coeff_it it=u.begin(); it!=u.end(); ++it) {
            (*it).second = tmp2[(*it).first];
        }
        std::cerr << "   Size of u after threshold: " << u.size() << std::endl;
        strategy2_applied = true;
    }

    int N = u.size();

    maxDof = 0;
    T _avDof = 0.;
    T current_theta = theta;
    T current_timestep = timestep;

    int refinement_timesteps_mod = numOfTimesteps / 8;

    std::stringstream filename;
    filename << "conv_awgm_in_thetascheme_" << u.size() << ".txt";
    std::ofstream file(filename.str().c_str());


    for (int i=1; i<=numOfTimesteps+2; ++i) {
//    for (int i=1; i<=numOfTimesteps; ++i) {


        current_theta    = theta;
        current_timestep = timestep;
        if (i<=4) {
            current_theta = 1.;
            current_timestep = 0.5*timestep;
            //current_timestep = timestep;
        }
        discrete_timepoint += current_timestep;


        std::cerr << "******* THETA scheme iteration i = " << i << " ****** " << std::endl;
        std::cerr << "   timestep = " << timestep << ", theta = " << current_theta << std::endl;

        timestep_solver.Op.setThetaTimeStepParameters(current_theta, current_timestep);
        timestep_solver.Prec.setThetaTimeStepParameters(current_theta, current_timestep);
        timestep_solver.F.setThetaTimeStepParameters(current_theta,current_timestep,discrete_timepoint,u);


        this->applyInvPreconditioner(u);
        T res = 0.;
        if (strategy==1 && discrete_timepoint>=0.1) {
            if (!strategy1_applied) {
                int sg_j = 0;
                for (sg_j=0; sg_j<20; ++sg_j) {
                    tmp.clear();
                    getSparseGridVector(timestep_solver.basis, tmp, sg_j, (T)0.);
                    if (tmp.size()>u.size()) break;
                }
                sg_j -= 1;
                u.clear();
                getSparseGridVector(timestep_solver.basis, u, sg_j, (T)0.);
                strategy1_applied = true;
            }
            res = timestep_solver.cg_solve(u,target_residual,1,0.001*timestep_eps,0.,
                                                         "conv.dat", "coeff.dat", N);
        }
        /*
        else if (strategy==2 && discrete_timepoint>=0.1) {
            res = timestep_solver.cg_solve(u,target_residual,1,0.0001*timestep_eps,0.,
                                             "conv.dat", "coeff.dat", N);
            target_residual = res;
            std::cerr << "   Residual = " << res << ", target residual was " << target_residual << std::endl;
        }
        */
        else if (strategy==3 && discrete_timepoint>=0.1 && i%refinement_timesteps_mod==0) {
            int sg_j = 0;
            int maxN = 0;
            for (sg_j=0; sg_j<20; ++sg_j) {
                tmp.clear();
                getSparseGridVector(timestep_solver.basis, tmp, sg_j, (T)0.);
                if      (tmp.size()>u.size()) break;
                else    maxN = tmp.size();
            }
            sg_j -= 1;
            u.clear();
            T target_res = 0.001*timestep_eps;
            getSparseGridVector(timestep_solver.basis, u, 2, (T)0.);
            res = timestep_solver.cg_solve(u,target_res,50,0.00001*target_res,0.,
                                                         "conv.dat", "coeff.dat", maxN);
            std::cerr << "   Adaptive solver: " << u.size() << " " << res << " " << target_res << " " << u.norm(2.) << std::endl;

        }
        else {
            res = timestep_solver.cg_solve(u,target_residual,1,1e-8,0.,
                                             "conv.dat", "coeff.dat", N);
        }
        file << u.size() << " " << res << " " << u.norm(2.) << std::endl;

        this->applyPreconditioner(u);

        IndexSet<Index1D> Lambda_x, Lambda_y;
        split(supp(u),Lambda_x,Lambda_y);
        int jmin_x=0, jmax_x=0, jmin_y=0, jmax_y=0;
        getMinAndMaxLevel(Lambda_x,jmin_x,jmax_x);
        getMinAndMaxLevel(Lambda_y,jmin_y,jmax_y);
        std::cout << "   Max level = (" << jmax_x << ", " << jmax_y << ")" << std::endl;

        if (strategy==2 && discrete_timepoint>=0.1) {
           tmp.clear();
           tmp2.clear();
           tmp   = u;
           tmp2  = u;
           //T threshold = exp(3*discrete_timepoint)*0.001*timestep_eps*u.norm(2.);
           //T threshold = 0.000001*timestep_eps*u.norm(2.);
           T threshold = 1e-6;
           std::cerr << "   Size of u before threshold: " << u.size() << std::endl;
           std::cerr << "   Threshold is: " << threshold << std::endl;
           int maxLevelSum = 0;
           bool erasedNode = false;
           for (const_coeff_it it=tmp.begin(); it!=tmp.end(); ++it) {
               maxLevelSum = std::max(maxLevelSum, (*it).first.index1.j+(*it).first.index2.j);
               if ((fabs((*it).second)<threshold)) {
                   if ((*it).first.index1.xtype==XBSpline || (*it).first.index2.xtype==XBSpline)
                           continue;
                   else  {
                       //std::cerr << "   Erasing " << (*it).first << std::endl;
                       u.erase((*it).first);
                       erasedNode = true;
                   }
               }
           }
           std::cerr << "   Size of u after threshold and no completion: " << u.size() << std::endl;
           std::cerr << "   Max level sum = " << maxLevelSum << std::endl;
           tmp.clear();
           if (erasedNode) {
               tmp = u;
               for (coeff_it it=tmp.begin(); it!=tmp.end(); ++it) {
                   completeMultiTree(timestep_solver.basis, (*it).first, u, 0, true, false);
               }
               tmp.clear();
           }
           for (coeff_it it=u.begin(); it!=u.end(); ++it) {
               (*it).second = tmp2[(*it).first];
           }
           std::cerr << "   Size of u after threshold: " << u.size() << std::endl;
           strategy2_applied = true;
       }

        maxDof = std::max(maxDof, (int)u.size());
        _avDof += u.size();

        std::cerr << std::endl;
        std::cerr << "   supp u = " << u.size() << std::endl << std::endl;
    }
    terminalDof = (int)u.size();
    _avDof /= (numOfTimesteps+2);
    avDof = (int)_avDof;
    std::cerr << "Theta timestepping with theta = " << theta << " and timestep " << timestep
              << " finished at t = " << discrete_timepoint << " with maxDof = " << maxDof
              << " and avDof = " << avDof << std::endl;
}

}   // namespace lawa

/*
if (current_theta!=1.) {
    timestep_solver.Op.evalA(u, propagated_u_k, "galerkin");
    propagated_u_k *= (current_theta-1.)*current_timestep;
}
timestep_solver.Op.evalM(u, propagated_u_k, "galerkin");
std::cerr << "DEBUG: propagated_u_k = " << propagated_u_k << std::endl;
*/
//timestep_solver.F.setThetaTimeStepParameters(current_theta,current_timestep,discrete_timepoint,&propagated_u_k);



/*
int sg_j=0;
if (strategy==1 && discrete_timepoint>=0.1 && !strategy1_applied) {
    std::cerr << "APPLYING strategy 1" << std::endl;
    std::cerr << "   Current size of index set before restriction: " << propagated_u_k.size() << std::endl;
    for (sg_j=0; sg_j<20; ++sg_j) {
        propagated_u_k.clear();
        getSparseGridVector(timestep_solver.basis, propagated_u_k, sg_j, (T)0.);
        std::cout << "    sg_j = " << sg_j << " : " <<  propagated_u_k.size() << " " << u.size() << std::endl;
        if (propagated_u_k.size()>u.size()) break;
    }
    sg_j -= 1;
    propagated_u_k.clear();
    getSparseGridVector(timestep_solver.basis, propagated_u_k, sg_j, (T)0.);
    std::cout << "   Final sg_j = " << sg_j << " " << propagated_u_k.size() << std::endl;
    std::cerr << "   Current size of index set after restriction: " << propagated_u_k.size() << std::endl;
}

if (strategy==1 && discrete_timepoint>=0.1 && !strategy1_applied) {
    std::cerr << "APPLYING strategy 1" << std::endl;
    u.clear();
    getSparseGridVector(timestep_solver.basis, u, sg_j, (T)0.);
    strategy1_applied = true;
}

*/

/*
 if (strategy==2 && i%refinement_timesteps_mod==0) {
    std::cerr << "Size of u before threshold: " << u.size() << ", ||u||_2" << u.norm(2.) << std::endl;
    u = MULTITREE_THRESH(u,(T)0.01*timestep_eps*u.norm(2.)); //0.001*timestep_eps*u.norm(2.)
    tmp = u;
    std::cerr << "Size of u after threshold: " << u.size() << std::endl;
    for (coeff_it it=tmp.begin(); it!=tmp.end(); ++it) {
        completeMultiTree(timestep_solver.basis, (*it).first, u, 0, true);
    }
    std::cerr << "Size of u after extension to multitree: " << u.size() << std::endl;
    propagated_u_k.clear();
    propagated_u_k = u;
    propagated_u_k.setToZero();
}
    //this->applyPreconditioner(u);
 */


/*
std::cerr << "Comparison of coefficients at t=0 and t=1:" << std::endl;
for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
    std::cerr << (*it).first << " : " << (*it).second << " " << u0[(*it).first] << std::endl;
}
*/

