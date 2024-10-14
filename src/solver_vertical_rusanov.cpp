// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Initial version: A. Ridley, June 2023

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Project the values to the faces
// --------------------------------------------------------------------------

void calc_facevalues_alts_rusanov(Grid &grid,
                                  arma_cube &inVar,
                                  arma_cube &outLeft,
                                  arma_cube &outRight) {

  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nZs = grid.get_nZ();
  int64_t nGCs = grid.get_nGCs();
  int64_t iX, iY, iZ;

  precision_t factor1 = 0.625;
  precision_t factor2 = 0.0416667;

  precision_t beta = 1.5;

  // inverse delta-alt at cell edge
  arma_mat ida(nXs, nYs);
  // Projection of variables to the cell edge:
  arma_mat dVarUp(nXs, nYs), dVarDown(nXs, nYs);

  arma_cube dVarLimited(nXs, nYs, nZs);

  dVarLimited.zeros();

  // Only do calculation on physical cells
  for (iZ = nGCs; iZ < nZs - nGCs; iZ++) {
    ida = 2.0 / grid.dalt_lower_scgc.slice(iZ + 1);
    dVarUp = ida %
             (factor1 * (inVar.slice(iZ + 1) - inVar.slice(iZ)) -
              factor2 * (inVar.slice(iZ + 2) - inVar.slice(iZ - 1)));

    ida = 2.0 / grid.dalt_lower_scgc.slice(iZ);
    dVarDown = ida %
               (factor1 * (inVar.slice(iZ) - inVar.slice(iZ - 1)) -
                factor2 * (inVar.slice(iZ + 1) - inVar.slice(iZ - 2)));

    for (iX = nGCs; iX < nXs - nGCs; iX++)
      for (iY = nGCs; iY < nYs - nGCs; iY++)
        dVarLimited(iX, iY, iZ) =
          limiter_mc(dVarUp(iX, iY), dVarDown(iX, iY), beta);

  }

  // Ghostcell closest to the bottom physical cell:
  iZ = nGCs - 1;
  ida = 1.0 / grid.dalt_lower_scgc.slice(iZ + 1);
  dVarUp = ida % (inVar.slice(iZ + 1) - inVar.slice(iZ));
  ida = 1.0 / grid.dalt_lower_scgc.slice(iZ);
  dVarDown = ida % (inVar.slice(iZ) - inVar.slice(iZ - 1));

  for (iX = nGCs; iX < nXs - nGCs; iX++)
    for (iY = nGCs; iY < nYs - nGCs; iY++)
      dVarLimited(iX, iY, iZ) =
        limiter_mc(dVarUp(iX, iY), dVarDown(iX, iY), beta);

  // Ghostcell closest to the top physical cell:
  iZ = nZs - nGCs;
  ida = 1.0 / grid.dalt_lower_scgc.slice(iZ + 1);
  dVarUp = ida % (inVar.slice(iZ + 1) - inVar.slice(iZ));
  ida = 1.0 / grid.dalt_lower_scgc.slice(iZ);
  dVarDown = ida % (inVar.slice(iZ) - inVar.slice(iZ - 1));

  for (iX = nGCs; iX < nXs - nGCs; iX++)
    for (iY = nGCs; iY < nYs - nGCs; iY++)
      dVarLimited(iX, iY, iZ) =
        limiter_mc(dVarUp(iX, iY), dVarDown(iX, iY), beta);
    
  for (iZ = nGCs; iZ < nZs - nGCs + 1; iZ++) {
    outLeft.slice(iZ) =
      inVar.slice(iZ - 1) +
      0.5 * dVarLimited.slice(iZ - 1) % grid.dalt_lower_scgc.slice(iZ);
    outRight.slice(iZ) =
      inVar.slice(iZ) -
      0.5 * dVarLimited.slice(iZ) % grid.dalt_lower_scgc.slice(iZ);
  }
  /*
  if (iProc == 11)
    std::cout << "facevalues : "
	      << inVar(7,19,17) << " "
	      << inVar(7,19,18) << " "
	      << inVar(7,19,19) << " "
	      << inVar(7,19,20) << " "
	      << dVarLimited(7,19,18) << " "
	      << grid.dalt_lower_scgc(7,19,17) << " "
	      << outRight(7, 19, 17) << " "
	      << outRight(7, 19, 18) << " "
	      << outLeft(7, 19, 17) << " "
	      << outLeft(7, 19, 18) << " "
	      << dVarUp(7, 19) << " "
	      << dVarDown(7, 19) << "\n";
  */
  
  return;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

void calc_grad_and_diff_alts_rusanov(Grid &grid,
                                     arma_cube &inVar,
                                     arma_cube &cMax,
                                     arma_cube &outGrad,
                                     arma_cube &outDiff) {

  std::string function = "calc_grad_and_diff_alts_rusanov";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nZs = grid.get_nZ();
  int64_t nGCs = grid.get_nGCs();
  arma_cube varLeft(nXs, nYs, nZs), varRight(nXs, nYs, nZs);
  arma_cube diffFlux(nXs, nYs, nZs);
  int64_t iX, iY, iZ;
  precision_t cMaxLocal, diffFluxLocal;

  outDiff.zeros();
  outGrad.zeros();

  report.print(3, "before facevalues");

  calc_facevalues_alts_rusanov(grid, inVar, varLeft, varRight);

  report.print(3, "after facevalues");

  for (iZ = nGCs; iZ < nZs - nGCs; iZ++)
    outGrad.slice(iZ) = 0.5 *
                        (varLeft.slice(iZ + 1) + varRight.slice(iZ + 1) -
                         varLeft.slice(iZ) - varRight.slice(iZ)) /
                        grid.dalt_center_scgc.slice(iZ);
  /*
  if (iProc == 11)
    std::cout << "calc_grad : "
	      << varLeft(7, 19, 17) << " "
	      << varLeft(7, 19, 18) << " "
	      << varRight(7, 19, 17) << " "
	      << varRight(7, 19, 18) << " "
	      << grid.dalt_center_scgc(7, 19, 17) << " "
	      << outGrad(7, 19, 17) << "\n";
  */
  for (iZ = nGCs; iZ < nZs - nGCs + 1; iZ++) {
    for (iX = nGCs; iX < nXs - nGCs; iX++)
      for (iY = nGCs; iY < nYs - nGCs; iY++) {
        if (cMax(iX, iY, iZ - 1) > cMax(iX, iY, iZ))
          cMaxLocal = cMax(iX, iY, iZ - 1);
        else
          cMaxLocal = cMax(iX, iY, iZ);

        diffFlux(iX, iY, iZ) =
          0.5 * cMaxLocal * (varRight(iX, iY, iZ) - varLeft(iX, iY, iZ));
        //if (iZ <= 10 && iX == 4 && iY == 4)
        //  std::cout << "diff flux : " << diffFlux(iX, iY, iZ)
        //      << " " << cMaxLocal
        //      << " " << varRight(iX, iY, iZ)
        //      << " " << varLeft(iX, iY, iZ) << "\n";
      }
  }

  report.print(3, "after diff flux");

  for (iZ = nGCs; iZ < nZs - nGCs; iZ++)
    outDiff.slice(iZ) =
      (diffFlux.slice(iZ + 1) - diffFlux.slice(iZ)) /
      grid.dalt_center_scgc.slice(iZ);

  report.exit(function);
  return;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

void Neutrals::solver_vertical_rusanov(Grid grid,
                                       Times time) {

  std::string function = "Neutrals::solver_vertical_rusanov";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nXs = grid.get_nX(), iX;
  int64_t nYs = grid.get_nY(), iY;
  int64_t nZs = grid.get_nZ(), iZ;
  int64_t nGCs = grid.get_nGCs();
  int iDir, iSpecies;

  precision_t dt = time.get_dt();

  // -----------------------------------------------------------
  // Bulk Variables:
  arma_cube gradTemp(nXs, nYs, nZs), diffTemp(nXs, nYs, nZs);
  std::vector<arma_cube> gradVel, diffVel;
  gradVel = make_cube_vector(nXs, nYs, nZs, 3);
  diffVel = make_cube_vector(nXs, nYs, nZs, 3);

  calc_grad_and_diff_alts_rusanov(grid,
                                  temperature_scgc,
                                  cMax_vcgc[2],
                                  gradTemp,
                                  diffTemp);

  arma_cube gradDummy(nXs, nYs, nZs), diffDummy(nXs, nYs, nZs);

  for (iDir = 0; iDir < 3; iDir++) {
    calc_grad_and_diff_alts_rusanov(grid,
                                    velocity_vcgc[iDir],
                                    cMax_vcgc[2],
                                    gradDummy,
                                    diffDummy);
    gradVel[iDir] = gradDummy;
    diffVel[iDir] = diffDummy;
  }

  // Vertical direction needs extra term:
  arma_cube divBulkVertVel(nXs, nYs, nZs);
  divBulkVertVel = gradVel[2] + 2 * velocity_vcgc[2] / grid.radius_scgc;

  // -----------------------------------------------------------
  // species dependent variables:
  std::vector<arma_cube> gradLogN_s, diffLogN_s;
  std::vector<arma_cube> gradVertVel_s, diffVertVel_s, divVertVel_s;
  gradLogN_s = make_cube_vector(nXs, nYs, nZs, nSpecies);
  diffLogN_s = make_cube_vector(nXs, nYs, nZs, nSpecies);
  gradVertVel_s = make_cube_vector(nXs, nYs, nZs, nSpecies);
  diffVertVel_s = make_cube_vector(nXs, nYs, nZs, nSpecies);
  divVertVel_s = make_cube_vector(nXs, nYs, nZs, nSpecies);

  arma_cube log_s(nXs, nYs, nZs), vv_s(nXs, nYs, nZs);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect) {

      // Log(number density):
      log_s = log(species[iSpecies].density_scgc);

      calc_grad_and_diff_alts_rusanov(grid,
                                      log_s,
                                      cMax_vcgc[2],
                                      gradDummy,
                                      diffDummy);

      gradLogN_s[iSpecies] = gradDummy;
      diffLogN_s[iSpecies] = diffDummy;

      // Vertical Velocity for each species:
      vv_s = species[iSpecies].velocity_vcgc[2];
      calc_grad_and_diff_alts_rusanov(grid,
                                      vv_s,
                                      cMax_vcgc[2],
                                      gradDummy,
                                      diffDummy);
      gradVertVel_s[iSpecies] = gradDummy;
      diffVertVel_s[iSpecies] = diffDummy;
      divVertVel_s[iSpecies] = gradDummy + 2 * vv_s / grid.radius_scgc;
    } else {
      gradVertVel_s[iSpecies].zeros();
      diffVertVel_s[iSpecies].zeros();
      divVertVel_s[iSpecies].zeros();
    }
  }

  // v2or = (Ve^2 + Vn^2)/R term:
  arma_cube v2or(nXs, nYs, nZs);
  v2or.zeros();
  v2or = (velocity_vcgc[0] % velocity_vcgc[0] +
          velocity_vcgc[1] % velocity_vcgc[1]) / grid.radius_scgc;

  // calculate vertical momentum due to eddy diffusion:
  vertical_momentum_eddy(grid);

  bool useImplicitFriction = true;

  if (useImplicitFriction) {
    calc_neutral_friction_coefs();
  }

  // -----------------------------------------------------------
  // Now calculate new states:
  precision_t mass;
  arma_cube one(nXs, nYs, nZs);
  one.ones();
  arma_cube gmo(nXs, nYs, nZs);
  gmo = gamma_scgc - one;

  arma_cube accTotal(nXs, nYs, nZs);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect) {

      mass = species[iSpecies].mass;

      // densities:
      log_s =
        log(species[iSpecies].density_scgc)
        - dt * (divVertVel_s[iSpecies] +
                species[iSpecies].velocity_vcgc[2] % gradLogN_s[iSpecies])
        + dt * diffLogN_s[iSpecies];
      species[iSpecies].newDensity_scgc = exp(log_s);

      accTotal = 
            dt * grid.gravity_vcgc[2]
          - dt * temperature_scgc % gradLogN_s[iSpecies] * cKB / mass
          + dt * diffVertVel_s[iSpecies]
          - dt * species[iSpecies].velocity_vcgc[2] % gradVertVel_s[iSpecies]
          + dt * v2or
          + dt * species[iSpecies].acc_eddy
          + dt * acc_coriolis[2]
          + dt * grid.cent_acc_vcgc[2];

      // vertical velocities:
      if (useImplicitFriction) {
        species[iSpecies].newVelocity_vcgc[2] =
          (species[iSpecies].velocity_vcgc[2] + accTotal + 
           dt * species[iSpecies].acc_neutral_friction[2]) /
          (1.0 + dt * species[iSpecies].neutral_friction_coef);
      } else {
        species[iSpecies].newVelocity_vcgc[2] =
          (species[iSpecies].velocity_vcgc[2] + accTotal);
      }

    } else {
      species[iSpecies].newVelocity_vcgc[2].zeros();
      species[iSpecies].newDensity_scgc = species[iSpecies].density_scgc;
    }
  }
  if (!useImplicitFriction)
    calc_neutral_friction(dt);

  newTemperature_scgc =
    temperature_scgc
    - dt * (velocity_vcgc[2] % gradTemp
            + gmo % (temperature_scgc % divBulkVertVel))
    + dt * diffTemp;

  precision_t maxVerticalVelocity = 100.0;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    if (species[iSpecies].DoAdvect)
      species[iSpecies].newVelocity_vcgc[2].clamp(-maxVerticalVelocity, maxVerticalVelocity);
  newTemperature_scgc.clamp(150, 1e32);

  for (iX = nGCs; iX < nXs - nGCs; iX++)
    for (iY = nGCs; iY < nYs - nGCs; iY++)
      for (iZ = nGCs; iZ < nZs - nGCs; iZ++) {
        temperature_scgc(iX, iY, iZ) = newTemperature_scgc(iX, iY, iZ);

        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          if (species[iSpecies].DoAdvect) {
            species[iSpecies].density_scgc(iX, iY, iZ) =
              species[iSpecies].newDensity_scgc(iX, iY, iZ);
            species[iSpecies].velocity_vcgc[2](iX, iY, iZ) =
              species[iSpecies].newVelocity_vcgc[2](iX, iY, iZ);
          } else {
            // assign bulk vertical velocity to the non-advected species:
            species[iSpecies].velocity_vcgc[2](iX, iY, iZ) =
              velocity_vcgc[2](iX, iY, iZ);
          }
        }
      }

  bool doPrintThis = false;
  if (doPrintThis) {
    iX = 2;
    iY = 2;
    iSpecies = 4;
    mass = species[iSpecies].mass;

    for (int iAlt = 0; iAlt < 20; iAlt++) {
      std::cout << iAlt << " " 
        << log(species[iSpecies].density_scgc(iX, iY,iAlt)) << " "
        << temperature_scgc(iX,iY,iAlt) << " "
        << species[iSpecies].velocity_vcgc[2](iX,iY,iAlt) << " "
        << temperature_scgc(iX,iY,iAlt) * gradLogN_s[iSpecies](iX,iY,iAlt) * cKB / mass << " "
        << gradTemp(iX,iY,iAlt) * cKB / mass << " "
        << grid.gravity_vcgc[2](iX,iY,iAlt) << "\n";
    }
  }
  //calc_neutral_friction();
/*
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect) {
        species[iSpecies].velocity_vcgc[2] = 
          species[iSpecies].velocity_vcgc[2] + dt * 
          species[iSpecies].acc_neutral_friction[2];
    }
  }
*/
  calc_mass_density();
  // Calculate bulk vertical winds:
  velocity_vcgc[2].zeros();
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    if (species[iSpecies].DoAdvect) {
      velocity_vcgc[2] = velocity_vcgc[2] + 
        species[iSpecies].mass * species[iSpecies].density_scgc % 
        species[iSpecies].velocity_vcgc[2] / rho_scgc;
    }

  report.exit(function);
  return;
}


// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

void Ions::solver_vertical_rusanov(Grid grid,
                                  Times time) {

  std::string function = "Ions::solver_vertical_rusanov";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nXs = grid.get_nX(), iX;
  int64_t nYs = grid.get_nY(), iY;
  int64_t nZs = grid.get_nZ(), iZ;
  int64_t nGCs = grid.get_nGCs();
  int iDir, iSpecies;

  precision_t dt = time.get_dt();

  // -----------------------------------------------------------
  // Bulk Variables:
  std::vector<arma_cube> gradVel, diffVel;
  gradVel = make_cube_vector(nXs, nYs, nZs, 3);
  diffVel = make_cube_vector(nXs, nYs, nZs, 3);

  arma_cube gradDummy(nXs, nYs, nZs), diffDummy(nXs, nYs, nZs);

  // -----------------------------------------------------------
  // species dependent variables:
  std::vector<arma_cube> gradLogN_s, diffLogN_s;
  std::vector<arma_cube> gradVertVel_s, diffVertVel_s, divVertVel_s;
  gradLogN_s = make_cube_vector(nXs, nYs, nZs, nSpecies);
  diffLogN_s = make_cube_vector(nXs, nYs, nZs, nSpecies);
  gradVertVel_s = make_cube_vector(nXs, nYs, nZs, nSpecies);
  diffVertVel_s = make_cube_vector(nXs, nYs, nZs, nSpecies);
  divVertVel_s = make_cube_vector(nXs, nYs, nZs, nSpecies);

  arma_cube log_s(nXs, nYs, nZs), vv_s(nXs, nYs, nZs);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect) {

      // Log(number density):
      log_s = log(species[iSpecies].density_scgc);

      calc_grad_and_diff_alts_rusanov(grid,
                                      log_s,
                                      cMax_vcgc[2],
                                      gradDummy,
                                      diffDummy);

      gradLogN_s[iSpecies] = gradDummy;
      diffLogN_s[iSpecies] = diffDummy;

      // Vertical Velocity for each species:
      vv_s = species[iSpecies].velocity_vcgc[2];
      calc_grad_and_diff_alts_rusanov(grid,
                                      vv_s,
                                      cMax_vcgc[2],
                                      gradDummy,
                                      diffDummy);
      gradVertVel_s[iSpecies] = gradDummy;
      diffVertVel_s[iSpecies] = diffDummy;
      divVertVel_s[iSpecies] = gradDummy + 2 * vv_s / grid.radius_scgc;
    } else {
      gradVertVel_s[iSpecies].zeros();
      diffVertVel_s[iSpecies].zeros();
      divVertVel_s[iSpecies].zeros();
    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect) {

      // densities:
      log_s =
        log(species[iSpecies].density_scgc)
        - dt * (divVertVel_s[iSpecies] +
                species[iSpecies].velocity_vcgc[2] % gradLogN_s[iSpecies])
        + dt * diffLogN_s[iSpecies];
      species[iSpecies].newDensity_scgc = exp(log_s);

    } else {
      species[iSpecies].newDensity_scgc = species[iSpecies].density_scgc;
    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    if (species[iSpecies].DoAdvect) 
      species[iSpecies].density_scgc = species[iSpecies].newDensity_scgc;

  fill_electrons();

  report.exit(function);
  return;
}
