int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DM             da;
  AppCtx         appctx;
  PetscBool      forwardonly = PETSC_FALSE,implicitform = PETSC_FALSE;
  PetscInt       perturbic = 1;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetBool(NULL,NULL,"-forwardonly",&forwardonly,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-implicitform",&implicitform,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-perturbic",&perturbic,NULL);CHKERRQ(ierr);

  appctx.D1    = 8.0e-5;
  appctx.D2    = 4.0e-5;
  appctx.gamma = .024;
  appctx.kappa = .06;

  /* Create distributed array (DMDA) to manage parallel grid and vectors */
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,DMDA_STENCIL_STAR,64,64,PETSC_DECIDE,PETSC_DECIDE,2,1,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,0,"u");CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,1,"v");CHKERRQ(ierr);

  /* Extract global vectors from DMDA; then duplicate for remaining
     vectors that are the same types */
  ierr = DMCreateGlobalVector(da,&appctx.U);CHKERRQ(ierr);

  /* Create timestepping solver context */
  ierr = TSCreate(PETSC_COMM_WORLD,&appctx.ts);CHKERRQ(ierr);
  ierr = TSSetType(appctx.ts,TSCN);CHKERRQ(ierr);
  ierr = TSSetDM(appctx.ts,da);CHKERRQ(ierr);
  ierr = TSSetProblemType(appctx.ts,TS_NONLINEAR);CHKERRQ(ierr);
  if (!implicitform) {
    ierr = TSSetRHSFunction(appctx.ts,NULL,RHSFunction,&appctx);CHKERRQ(ierr);
    ierr = TSSetRHSJacobian(appctx.ts,NULL,NULL,RHSJacobian,&appctx);CHKERRQ(ierr);
  } else {
    ierr = TSSetIFunction(appctx.ts,NULL,IFunction,&appctx);CHKERRQ(ierr);
    ierr = TSSetIJacobian(appctx.ts,NULL,NULL,IJacobian,&appctx);CHKERRQ(ierr);
  }

  /* Set initial conditions */
  ierr = InitialConditions(da,appctx.U);CHKERRQ(ierr);
  ierr = TSSetSolution(appctx.ts,appctx.U);CHKERRQ(ierr);

  /* Set solver options */
  ierr = TSSetMaxTime(appctx.ts,2000.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(appctx.ts,0.5);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(appctx.ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetFromOptions(appctx.ts);CHKERRQ(ierr);

  ierr = GenerateOBs(appctx.ts,appctx.U,&appctx);CHKERRQ(ierr);

  if (!forwardonly) {
    Tao tao;
    Vec P;
    Vec lambda[1];

    if (perturbic == 1) {
      ierr = PerturbedInitialConditions(da,appctx.U);CHKERRQ(ierr);
    } else if (perturbic == 2) {
      ierr = PerturbedInitialConditions2(da,appctx.U);CHKERRQ(ierr);
    } else if (perturbic == 3) {
      ierr = PerturbedInitialConditions3(da,appctx.U);CHKERRQ(ierr);
    }

    ierr = VecDuplicate(appctx.U,&lambda[0]);CHKERRQ(ierr);
    ierr = TSSetCostGradients(appctx.ts,1,lambda,NULL);CHKERRQ(ierr);

    /* Have the TS save its trajectory needed by TSAdjointSolve() */
    ierr = TSSetSaveTrajectory(appctx.ts);CHKERRQ(ierr);

    /* Create TAO solver and set desired solution method */
    ierr = TaoCreate(PETSC_COMM_WORLD,&tao);CHKERRQ(ierr);
    ierr = TaoSetType(tao,TAOBLMVM);CHKERRQ(ierr);

    /* Set initial guess for TAO */
    ierr = VecDuplicate(appctx.U,&P);CHKERRQ(ierr);
    ierr = VecCopy(appctx.U,P);CHKERRQ(ierr);
    ierr = TaoSetInitialVector(tao,P);CHKERRQ(ierr);

    /* Set routine for function and gradient evaluation */
    ierr = TaoSetObjectiveAndGradientRoutine(tao,FormFunctionAndGradient,&appctx);CHKERRQ(ierr);

    /* Check for any TAO command line options */
    ierr = TaoSetFromOptions(tao);CHKERRQ(ierr);

    /*
    ierr = TaoGetKSP(tao,&ksp);CHKERRQ(ierr);
    if (ksp) {
      ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
      ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
    }
    */

    ierr = TaoSolve(tao);CHKERRQ(ierr);
    ierr = TaoDestroy(&tao);CHKERRQ(ierr);
    ierr = VecDestroy(&lambda[0]);CHKERRQ(ierr);
    ierr = VecDestroy(&P);CHKERRQ(ierr);
  }

  /* Free work space.  All PETSc objects should be destroyed when they
     are no longer needed. */
  ierr = VecDestroy(&appctx.U);CHKERRQ(ierr);
  ierr = TSDestroy(&appctx.ts);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}