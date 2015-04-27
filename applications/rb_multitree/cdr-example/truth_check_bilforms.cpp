#include "cdr_problem.h"

typedef AdaptiveSpaceTimePDEOperator1D<T, Basis2D_Trial,
NoPrec2D, NoPrec2D, NoInitialCondition>                                 SpaceTimeOp2D;
typedef AdaptiveSpaceTimeTimeDerivOperator1D<T, Basis2D_Trial,
NoPrec2D, NoPrec2D, NoInitialCondition>                                 SpaceTimeTimeDerivOp2D;
typedef AdaptiveSpaceTimeLaplaceOperator1D<T, Basis2D_Trial,
NoPrec2D, NoPrec2D, NoInitialCondition>                                 SpaceTimeLaplaceOp2D;
typedef AdaptiveSpaceTimeWeightedConvectionOperator1D<T, Basis2D_Trial,
NoPrec2D, NoPrec2D, NoInitialCondition>                                 SpaceTimeWeightedConvectionOp2D;
typedef AdaptiveSpaceTimeReactionOperator1D<T, Basis2D_Trial,
NoPrec2D, NoPrec2D, NoInitialCondition>                                 SpaceTimeReactionOp2D;


int main () {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//


    int d   = 2;
    int d_  = 2;
    int j0  = 2;

    //getchar();

    /// Basis initialization
    TrialBasis_Time      basis_per(d,d_,j0);
    TestBasis_Time       basis_int(d,d_,j0);
    Basis_Space 		 basis_intbc(d,d_,j0);
    basis_intbc.enforceBoundaryCondition<DirichletBC>();

    Basis2D_Trial basis2d_trial(basis_per,basis_intbc);
    Basis2D_Test  basis2d_test(basis_int,basis_intbc);

    /// Initialization of operators
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);
    Function<T> w_conv_Fct(weight_convection,no_singPts);

    // Bilinear Forms
    Convection1D_Time			ConvectionBil_t(basis_per, basis_int);
    Identity1D_Time 		    IdentityBil_t(basis_per, basis_int);
    Identity1D_Space 	        IdentityBil_x(basis_intbc, basis_intbc);
    Laplace1D_Space 	        LaplaceBil_x(basis_intbc, basis_intbc);
    Convection1D_Space 	        ConvectionBil_x(basis_intbc, basis_intbc,
    											zero_Fct, w_conv_Fct, zero_Fct, 10,true, false, true);
    Identity1D_Space 	        IdentityBil_Test_t(basis_int, basis_int);
    Identity1D_Time_Trial		IdentityBil_Trial_t(basis_per, basis_per);
    Convection1D_Time_Trial		ConvectionBil_Trial_t(basis_per, basis_per);


    RefConvection1D_Time 		RefConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Time 		    RefIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Space 	    RefIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefLaplace1D_Space 	        RefLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefConvection1D_Space 	    RefConvectionBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis,
													zero_Fct, w_conv_Fct, zero_Fct, 10,true, false, true);
    RefIdentity1D_Space 	    RefIdentityBil_Test_t(basis_int.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Time_Trial	RefIdentityBil_Trial_t(basis_per.refinementbasis, basis_per.refinementbasis);
    RefConvection1D_Time_Trial	RefConvectionBil_Trial_t(basis_per.refinementbasis, basis_per.refinementbasis);

    // Transposed Bilinear Forms
    TranspConvection1D_Time 	TranspConvectionBil_t(basis_per, basis_int);
    TranspIdentity1D_Time 		TranspIdentityBil_t(basis_per, basis_int);
    TranspIdentity1D_Space 	    TranspIdentityBil_x(basis_intbc, basis_intbc);
    TranspLaplace1D_Space 	    TranspLaplaceBil_x(basis_intbc, basis_intbc);
    TranspConvection1D_Space 	TranspConvectionBil_x(basis_intbc, basis_intbc,
    													zero_Fct, w_conv_Fct, zero_Fct, 10,true, false, true);
    
    RefTranspConvection1D_Time 	RefTranspConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefTranspIdentity1D_Time 	RefTranspIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefTranspIdentity1D_Space 	RefTranspIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefTranspLaplace1D_Space 	RefTranspLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefTranspConvection1D_Space RefTranspConvectionBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis,
    													zero_Fct, w_conv_Fct, zero_Fct, 10, true, false, true);

    /// Initialization of local operator
    LOp_Conv1D_Time			lOp_Conv1D_t	(basis_int, basis_per, RefConvectionBil_t, ConvectionBil_t);
    LOp_Id1D_Time			lOp_Id1D_t  	(basis_int, basis_per, RefIdentityBil_t, IdentityBil_t);
    LOp_Id1D_Space			lOp_Id1D_x  	(basis_intbc, basis_intbc, RefIdentityBil_x, IdentityBil_x);
    LOp_Lapl1D_Space		lOp_Lapl1D_x	(basis_intbc, basis_intbc, RefLaplaceBil_x, LaplaceBil_x);
    LOp_Conv1D_Space		lOp_Conv1D_x	(basis_intbc, basis_intbc, RefConvectionBil_x, ConvectionBil_x);
    LOp_Id1D_Space			lOp_Id1D_Test_t (basis_int, basis_int, RefIdentityBil_Test_t, IdentityBil_Test_t);
    LOp_Id1D_Time_Trial		lOp_Id1D_Trial_t(basis_per, basis_per, RefIdentityBil_Trial_t, IdentityBil_Trial_t);
    LOp_Conv1D_Time_Trial 	lOp_Conv1D_Trial_t(basis_per, basis_per, RefConvectionBil_Trial_t, ConvectionBil_Trial_t);
    
    LOpT_Conv1D_Time		lOpT_Conv1D_t(basis_per, basis_int, RefTranspConvectionBil_t, TranspConvectionBil_t);
    LOpT_Id1D_Time			lOpT_Id1D_t  (basis_per, basis_int, RefTranspIdentityBil_t, TranspIdentityBil_t);
    LOpT_Id1D_Space			lOpT_Id1D_x  (basis_intbc, basis_intbc, RefTranspIdentityBil_x, TranspIdentityBil_x);
    LOpT_Lapl1D_Space		lOpT_Lapl1D_x(basis_intbc, basis_intbc, RefTranspLaplaceBil_x, TranspLaplaceBil_x);
    LOpT_Conv1D_Space		lOpT_Conv1D_x(basis_intbc, basis_intbc, RefTranspConvectionBil_x, TranspConvectionBil_x);

    LOp_Conv_Id_2D			localConvectionIdentityOp2D		(lOp_Conv1D_t, 		lOp_Id1D_x);
    LOp_Id_Id_2D			localIdentityIdentityOp2D		(lOp_Id1D_t, 		lOp_Id1D_x);
    LOp_Id_Lapl_2D			localIdentityLaplaceOp2D		(lOp_Id1D_t, 		lOp_Lapl1D_x);
    LOp_Id_Conv_2D			localIdentityConvectionOp2D		(lOp_Id1D_t, 		lOp_Conv1D_x);
    LOp_Test_Id_Id_2D		localIdentityIdentityOp2D_Test	(lOp_Id1D_Test_t, 	lOp_Id1D_x);
    LOp_Test_Id_Lapl_2D		localIdentityLaplaceOp2D_Test	(lOp_Id1D_Test_t, 	lOp_Lapl1D_x);
    LOp_Trial_Id_Id_2D		localIdentityIdentityOp2D_Trial	(lOp_Id1D_Trial_t, 	lOp_Id1D_x);
    LOp_Trial_Id_Lapl_2D 	localIdentityLaplaceOp2D_Trial	(lOp_Id1D_Trial_t, 	lOp_Lapl1D_x);
    LOp_Trial_Id_Conv_2D 	localIdentityConvectionOp2D_Trial(lOp_Id1D_Trial_t, lOp_Conv1D_x);
    LOp_Trial_Conv_Id_2D 	localConvectionIdentityOp2D_Trial(lOp_Conv1D_Trial_t,lOp_Id1D_x);

    LOpT_Conv_Id_2D		transpLocalConvectionIdentityOp2D	(lOpT_Conv1D_t, lOpT_Id1D_x);
    LOpT_Id_Id_2D		transpLocalIdentityIdentityOp2D		(lOpT_Id1D_t, 	lOpT_Id1D_x);
    LOpT_Id_Lapl_2D		transpLocalIdentityLaplaceOp2D		(lOpT_Id1D_t, 	lOpT_Lapl1D_x);
    LOpT_Id_Conv_2D		transpLocalIdentityConvectionOp2D	(lOpT_Id1D_t, 	lOpT_Conv1D_x);

    localConvectionIdentityOp2D.setJ(9);
    localIdentityIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    localIdentityConvectionOp2D.setJ(9);
    transpLocalConvectionIdentityOp2D.setJ(9);
    transpLocalIdentityIdentityOp2D.setJ(9);
    transpLocalIdentityLaplaceOp2D.setJ(9);
    transpLocalIdentityConvectionOp2D.setJ(9);
    localIdentityIdentityOp2D_Test.setJ(9);
    localIdentityLaplaceOp2D_Test.setJ(9);
    localIdentityIdentityOp2D_Trial.setJ(9);
    localIdentityConvectionOp2D_Trial.setJ(9);
    localIdentityLaplaceOp2D_Trial.setJ(9);
    localConvectionIdentityOp2D_Trial.setJ(9);

    /// Initialization of preconditioner
    LeftPrec2D leftPrec(basis2d_test);
    RightPrec2D rightPrec(basis2d_trial);

    NoPrec2D noPrec;


    //======= SETUP SPACETIME-OPERATORS ===============//


    /* Affine Decomposition Left Hand Side */
    DenseVectorT nosingpts(2);
    nosingpts = 0., 1.;
    Function<T>  w_conv(weight_convection, nosingpts);
    
    SpaceTimeOp2D                   spacetimenorm(basis2d_trial,noPrec, noPrec, 1., 0., 1., 0.);
    SpaceTimeTimeDerivOp2D          A1(basis2d_trial, noPrec, noPrec);          // time derivative
    SpaceTimeLaplaceOp2D            A2(basis2d_trial, noPrec, noPrec);          // space laplacian
    SpaceTimeWeightedConvectionOp2D A3(basis2d_trial, noPrec, noPrec, w_conv);  // space convection (weighted)
    SpaceTimeReactionOp2D           A4(basis2d_trial, noPrec, noPrec);          // space reaction



	//===============================================================//
	//===============  RB SETUP =====================================//
	//===============================================================//

    // Affine Decompositions:
    // 	Left Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct>	lhs_theta_fcts;
    lhs_theta_fcts.push_back(no_theta);
    lhs_theta_fcts.push_back(no_theta);
    lhs_theta_fcts.push_back(theta_conv);
    lhs_theta_fcts.push_back(theta_reac);
    ThetaStructure<ParamType> lhs_theta(lhs_theta_fcts);

    vector<AbstractLocalOperator2D<T>* > lhs_ops, lhs_opsT;
    lhs_ops.push_back(&localConvectionIdentityOp2D);
    lhs_ops.push_back(&localIdentityLaplaceOp2D);
    lhs_ops.push_back(&localIdentityConvectionOp2D);
    lhs_ops.push_back(&localIdentityIdentityOp2D);
    lhs_opsT.push_back(&transpLocalConvectionIdentityOp2D);
    lhs_opsT.push_back(&transpLocalIdentityLaplaceOp2D);
    lhs_opsT.push_back(&transpLocalIdentityConvectionOp2D);
    lhs_opsT.push_back(&transpLocalIdentityIdentityOp2D);

    Affine_Op_2D affine_lhs(lhs_theta, lhs_ops);
    Affine_Op_2D affine_lhs_T(lhs_theta, lhs_opsT);

    COp_TestInnProd innprod_Y	 (localIdentityIdentityOp2D_Test, localIdentityLaplaceOp2D_Test);
    TrialInnProdY 	innprod_Y_u_u(localIdentityIdentityOp2D_Trial, localIdentityLaplaceOp2D_Trial);

    vector<AbstractLocalOperator2D<T>* > A_u_u_ops;
    A_u_u_ops.push_back(&localConvectionIdentityOp2D_Trial);
    A_u_u_ops.push_back(&localIdentityLaplaceOp2D_Trial);
    A_u_u_ops.push_back(&localIdentityConvectionOp2D_Trial);
    A_u_u_ops.push_back(&localIdentityIdentityOp2D_Trial);
    Flex_COp_2D A_u_u(A_u_u_ops);


	//===============================================================//
	//=============== COMPARISON =====================================//
	//===============================================================//


    IndexSet<Index2D> LambdaTrial;
    //getFullIndexSet(basis2d_trial, LambdaTrial, 2,2,0);
    getScalingFctIndexSet(basis2d_trial, LambdaTrial, 5,5);

/*
    for(auto& index1 : LambdaTrial){
    	for(auto& index2 : LambdaTrial){
    		cout << "A1:" << (((A_u_u.eval(0,index1, index2) - A1(index1, index2))>1e-14)?"ERROR":"CHECK") << " -- " <<  index1 << " " << index2 << endl;
    		cout << "A2:" << (((A_u_u.eval(1,index1, index2) - A2(index1, index2))>1e-14)?"ERROR":"CHECK") << " -- " <<  index1 << " " << index2 << endl;
    		cout << "A3:" << (((A_u_u.eval(2,index1, index2) - A3(index1, index2))>1e-14)?"ERROR":"CHECK") << " -- " <<  index1 << " " << index2 << " "<< A_u_u.eval(2,index1, index2) - A3(index1, index2) << endl;
    		cout << "A4:" << (((A_u_u.eval(3,index1, index2) - A4(index1, index2))>1e-14)?"ERROR":"CHECK") << " -- " <<  index1 << " " << index2 << endl;
    		cout << "I:" << (((innprod_Y_u_u.eval(index1, index2) - spacetimenorm(index1, index2))>1e-14)?"ERROR":"CHECK") << " -- " <<  index1 << " " << index2 << endl;
    	}
    }
    */
    
    Index1D i1(5,2,XBSpline);
    Index1D i1_b(5,3,XBSpline);
    Index2D i2(i1,i1);
    Index2D i2_b(i1_b, i1);
    
    cout.precision(8);
    
    cout << "Index : " << i2 << " x " << i2 << endl;
    cout << "A1: " << A_u_u.eval(0,i2, i2) << " " << A1(i2, i2) << endl;
    cout << "A2: " << A_u_u.eval(1,i2, i2) << " " << A2(i2, i2) << endl;
    cout << "A3: " << A_u_u.eval(2,i2, i2) << " " << A3(i2, i2) << endl;
    cout << "A4: " << A_u_u.eval(3,i2, i2) << " " << A4(i2, i2) << endl;
    cout << "I : " << innprod_Y_u_u.eval(i2, i2) << " " << spacetimenorm(i2, i2) << endl;
    
    cout << "Index : " << i2 << " x " << i2_b << endl;
    cout << "A1: " << A_u_u.eval(0,i2, i2_b) << " " << A1(i2, i2_b) << endl;
    cout << "A2: " << A_u_u.eval(1,i2, i2_b) << " " << A2(i2, i2_b) << endl;
    cout << "A3: " << A_u_u.eval(2,i2, i2_b) << " " << A3(i2, i2_b) << endl;
    cout << "A4: " << A_u_u.eval(3,i2, i2_b) << " " << A4(i2, i2_b) << endl;
    cout << "I : " << innprod_Y_u_u.eval(i2, i2_b) << " " << spacetimenorm(i2, i2_b) << endl;
    
    cout << "Index : " << i2_b << " x " << i2 << endl;
    cout << "A1: " << A_u_u.eval(0,i2_b, i2) << " " << A1(i2_b, i2) << endl;
    cout << "A2: " << A_u_u.eval(1,i2_b, i2) << " " << A2(i2_b, i2) << endl;
    cout << "A3: " << A_u_u.eval(2,i2_b, i2) << " " << A3(i2_b, i2) << endl;
    cout << "A4: " << A_u_u.eval(3,i2_b, i2) << " " << A4(i2_b, i2) << endl;
    cout << "I : " << innprod_Y_u_u.eval(i2_b, i2) << " " << spacetimenorm(i2_b, i2) << endl;

    return 0;
}

