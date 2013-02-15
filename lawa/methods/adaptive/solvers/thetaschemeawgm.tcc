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
                                                  int &avDof, int &maxDof, int &terminalDof, int j)
{
    size_t hms;
    hms = u.bucket_count();

    T discrete_timepoint = 0.;
    Coefficients<Lexicographical,T,Index> propagated_u_k(hms), u0(hms), tmp(hms), tmp2(hms), tmp3(hms);
    u0 = u;

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
        std::cerr << "   timestep = " << timestep << ", theta = " << current_theta
                  << ", time point = " << discrete_timepoint << std::endl;

        timestep_solver.Op.setThetaTimeStepParameters(current_theta, current_timestep);
        timestep_solver.Prec.setThetaTimeStepParameters(current_theta, current_timestep);
        timestep_solver.F.setThetaTimeStepParameters(current_theta,current_timestep,discrete_timepoint,u);

        this->applyInvPreconditioner(u);
        T res = 0.;
        if (strategy==3 && discrete_timepoint>=0.1) {
            res = timestep_solver.cg_solve(u,1e-12,1,1e-10,0.,
                                             "conv.dat", "coeff.dat", N);
        }
        else {
            res = timestep_solver.cg_solve(u,1e-12,1,1e-10,0.,
                                             "conv.dat", "coeff.dat", N);
        }

        this->applyPreconditioner(u);


        IndexSet<Index1D> Lambda_x, Lambda_y;
        split(supp(u),Lambda_x,Lambda_y);
        int jmin_x=0, jmax_x=0, jmin_y=0, jmax_y=0;
        getMinAndMaxLevel(Lambda_x,jmin_x,jmax_x);
        getMinAndMaxLevel(Lambda_y,jmin_y,jmax_y);
        std::cout << "   Max level = (" << jmax_x << ", " << jmax_y << ")" << std::endl;


        tmp.clear();
        tmp2.clear();
        tmp2  = u;
        bool erasedNode = false;


        if (strategy==1 && i>4) {
            //T threshold = j <= 5 ? 0.0001*std::pow((T)2.,(T)(-2.*j)) : 0.00001*std::pow((T)2.,(T)(-2.*j));//T threshold = 1e-12; //ex1
            T threshold = j <= 5 ? 0.0005*std::pow((T)2.,(T)(-2.*j)) : 0.00005*std::pow((T)2.,(T)(-2.*j));//T threshold = 1e-12; //ex2
            std::cerr << "   APPLYING strategy 1" << std::endl;
            int sizeBeforeThresh = u.size(), sizeAfterThresh = 0;
            u = MULTITREE_THRESH(u,threshold*u.norm(2.));
            erasedNode = true;
            sizeAfterThresh = u.size();
            (sizeBeforeThresh > sizeAfterThresh) ? erasedNode = true : erasedNode = false;
            std::cerr << "      Threshold: " << threshold << " -> (" << sizeBeforeThresh << ", "
                      << sizeAfterThresh << ")" << std::endl;
       }

        if (strategy==2 && discrete_timepoint>=0.125) {
           //T threshold = j <= 5 ? 0.0001*std::pow((T)2.,(T)(-2.*j)) : 0.00001*std::pow((T)2.,(T)(-2.*j));
           //T threshold = j <= 5 ? 0.001*std::pow((T)2.,(T)(-2.*j)) : 0.0001*std::pow((T)2.,(T)(-2.*j));

            T threshold =  0.00001*std::pow((T)2.,(T)(-2.*j))*u.norm(2.);
//            T threshold =  0.000001*std::pow((T)2.,(T)(-2.*j))*u.norm(2.);
            int sizeBeforeThresh = u.size(), sizeAfterThresh = 0;
            std::cerr << "   APPLYING strategy 2" << std::endl;
            int maxLevelSum = 0;
            tmp = u;
            for (const_coeff_it it=tmp.begin(); it!=tmp.end(); ++it) {
               maxLevelSum = std::max(maxLevelSum, (*it).first.index1.j+(*it).first.index2.j);

               if ((fabs((*it).second)<threshold)) {
                   if ( (((*it).first.index1.xtype==XBSpline) || ((*it).first.index2.xtype==XBSpline))
                        && (fabs((*it).second)>1e-10)   )
                           continue;
                   //else  {
                       //std::cerr << "   Erasing " << (*it).first << std::endl;
                   u.erase((*it).first);
                   //}
               }
            }
            sizeAfterThresh = u.size();
            (sizeBeforeThresh > sizeAfterThresh) ? erasedNode = true : erasedNode = false;
            std::cerr << "      Threshold: " << threshold << " -> (" << sizeBeforeThresh << ", "
                     << sizeAfterThresh << ")" << std::endl;
            std::cerr << "   Max level sum = " << maxLevelSum << std::endl;
            tmp.clear();
        }

        if ((strategy==1 || strategy==2) && erasedNode) {
            std::cerr << "   Completing multitree." << std::endl;
            tmp = u;
            for (coeff_it it=tmp.begin(); it!=tmp.end(); ++it) {
                completeMultiTree(timestep_solver.basis, (*it).first, u, 0, true, false);
            }
            tmp.clear();

            for (coeff_it it=u.begin(); it!=u.end(); ++it) {
                (*it).second = tmp2[(*it).first];
            }
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

/*
int  j_x=(*it).first.index1.j,    j_y=(*it).first.index2.j;
long k_x=(*it).first.index1.k,    k_y=(*it).first.index2.k;
XType xtype_x=(*it).first.index1.xtype, xtype_y=(*it).first.index2.xtype;

Support<T> supp_x = timestep_solver.basis.first.generator(xtype_x).support(j_x,k_x);
T center_x  = 0.5*(supp_x.l1 + supp_x.l2);
center_x *= 4.; center_x += -2.;

Support<T> supp_y = timestep_solver.basis.second.generator(xtype_y).support(j_y,k_y);
T center_y  = 0.5*(supp_y.l1 + supp_y.l2);
center_y *= 4.; center_y += -2.;
T factor = exp(5.*(center_x*center_x + center_y*center_y));

 *
 */

