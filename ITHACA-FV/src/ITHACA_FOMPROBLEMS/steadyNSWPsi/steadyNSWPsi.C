#include "steadyNSWPsi.H"
#include "steadyNS.H"
#include "viscosityModel.H"


steadyNSWPsi::steadyNSWPsi() : steadyNS() {}
steadyNSWPsi::steadyNSWPsi(int argc, char* argv[])  :
    steadyNS()
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();
#include "createFields.H"
#include "createFvOptions.H"
    turbulence->validate();
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty" || bcMethod == "none",
             "The BC method must be set to lift or penalty or none in ITHACAdict");
    fluxMethod = ITHACAdict->lookupOrDefault<word>("fluxMethod", "inconsistent");
    M_Assert(fluxMethod == "inconsistent" || bcMethod == "consistent",
             "The flux method must be set to inconsistent or consistent in ITHACAdict");
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
};

// Method to perform a truthSolve
void steadyNSWPsi::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U(); 
    volScalarField& W = _W();
    volScalarField& Psi_z = _Psi_z();
    volVectorField& Psi = _Psi();
    volVectorField& temp = _temp();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    #include "NLsolvesteadyNS.H"
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/"); 
    ITHACAstream::exportSolution(W, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(Psi_z, name(counter), "./ITHACAoutput/Offline/");
    Wfield.append(W.clone());
    Psi_zfield.append(Psi_z.clone());
    counter++;
    writeMu(mu_now);
    // --- Fill in the mu_samples with parameters (mu) to be used for the PODI sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (label i = 0; i < mu_now.size(); i++)
    {
        mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}

///////////////////////////////
void steadyNSWPsi::projectSUP(fileName folder, label NWmodes, label NPsi_zmodes)
{
    /*NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP; */
    NWmodes = NWmodes;
    NPsi_zmodes = NPsi_zmodes;
    /* L_U_SUPmodes.resize(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(liftfield[k].clone());
        }
    }  */
////////MODIFIED WITH W
    /* if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(Umodes[k].clone());
        }
    } 

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            L_U_SUPmodes.append(supmodes[k].clone());
        }
    }  */

    if (NWmodes != 0)
    {
        for (label k = 0; k < NWmodes; k++)
        {
            W_proj_modes.append(W_proj_modes[k].clone());    ///il primo era L_U_SUPmodes
        }
    } 

    

     if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        /* word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        }
        else
        {
            B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        }

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        }
        else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        }

        word P_str = "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
        }
        else
        {
            P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        }

        word M_str = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
        {
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
        }
        else
        {
            M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        }
        else
        {
            C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        }  */

        word AWPsi_str = "AWPsi_" + name(NWmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + AWPsi_str))
        {
            ITHACAstream::ReadDenseMatrix(AWPsi_matrix, "./ITHACAoutput/Matrices/", AWPsi_str);
        }
        else
        {
            AWPsi_matrix = diffusiveW_term(NWmodes, NPsi_zmodes);;
        }

        word BWPsi_str = "BWPsi_" + name(NWmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + BWPsi_str))
        {
            ITHACAstream::ReadDenseMatrix(BWPsi_matrix, "./ITHACAoutput/Matrices/", BWPsi_str);
        }
        else
        {
            BWPsi_matrix = diffusivePsi_z_term(NWmodes, NPsi_zmodes);
        }

        word MW_str = "MW_" + name(NWmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + MW_str))
        {
            ITHACAstream::ReadDenseMatrix(MW_matrix, "./ITHACAoutput/Matrices/", MW_str);
        }
        else
        {
            MW_matrix = massW_term(NWmodes, NPsi_zmodes);
        }

        word MPsi_str = "MPsi_" + name(NPsi_zmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + MPsi_str))
        {
            ITHACAstream::ReadDenseMatrix(MPsi_matrix, "./ITHACAoutput/Matrices/", MPsi_str);
        }
        else
        {
            MPsi_matrix = massPsi_z_term(NWmodes, NPsi_zmodes);
        }

        word GWPsi_str = "GWPsi_" + name(NWmodes)+ "_" +name(NPsimodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + GWPsi_str))
        {
            ITHACAstream::ReadDenseTensor(GWPsi_tensor, "./ITHACAoutput/Matrices/", GWPsi_str);
        }
        else
        {
            GWPsi_tensor = convective_term_tens(NWmodes, NPsi_zmodes);
        }


        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        } 
    }
    else
    {
        /* B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes); */
 
        AWPsi_matrix = diffusiveW_term(NWmodes, NPsi_zmodes);
        BWPsi_matrix = diffusivePsi_z_term(NWmodes, NPsi_zmodes);
        GWPsi_tensor = convective_term_tens(NWmodes, NPsi_zmodes);
        //KWPsi_matrix = pressure_gradient_term(NWmodes, NPsi_zmodes);
        //PWPsi_matrix = divergence_term(NWmodes, NPsi_zmodes);
        MW_matrix = massW_term(NWmodes, NPsi_zmodes);
        MPsi_matrix = massPsi_z_term(NWmodes, NPsi_zmodes);

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

}


Eigen::MatrixXd steadyNSWPsi::diffusiveW_term(label NWmodes, label NPsi_zmodes)
{
    label AWPsi_size = NWmodes;
    Eigen::MatrixXd AWPsi_matrix;
    AWPsi_matrix.resize(AWPsi_size, AWPsi_size);

    /* //Project everything
    for (label i = 0; i < Bsize; i++)
    {
        for (label j = 0; j < Bsize; j++)
        {
            B_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), L_U_SUPmodes[j])).value();
        }
    } */

    for (label i = 0; i < AWPsi_size; i++)     //L_U_SUPmodes = modes W
    {
        for (label j = 0; j < AWPsi_size; j++)
        {
            AWPsi_matrix(i, j) = fvc::domainIntegrate(W_proj_modes[i] & fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), W_proj_modes[j])).value();
        }
    }


    /* if (Pstream::parRun())
    {
        reduce(B_matrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/",
                                      "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    }  */

    return AWPsi_matrix;
}

Eigen::MatrixXd steadyNSWPsi::diffusivePsi_z_term(label NWmodes, label NPsi_zmodes)
{
    label BWPsi_size = NPsi_zmodes;
    Eigen::MatrixXd BWPsi_matrix;
    BWPsi_matrix.resize(BWPsi_size, BWPsi_size);

    /* //Project everything
    for (label i = 0; i < Bsize; i++)
    {
        for (label j = 0; j < Bsize; j++)
        {
            B_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), L_U_SUPmodes[j])).value();
        }
    }  */

    for (label i = 0; i < BWPsi_size; i++)     //L_U_SUPmodes = modes W
    {
        for (label j = 0; j < BWPsi_size; j++)
        {
            BWPsi_matrix(i, j) = fvc::domainIntegrate(Psi_z_proj_modes[i] & fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), Psi_z_proj_modes[j])).value();
        }
    }


    /* if (Pstream::parRun())
    {
        reduce(B_matrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/",
                                      "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    }  */

    return BWPsi_matrix;
}


Eigen::Tensor<double, 3> steadyNSWPsi::convective_term_tens(label NWmodes,
        label NPsi_zmodes)

{
    label CW_size = NWmodes;
    label CPsi_size = NPsi_zmodes;
    Eigen::Tensor<double, 3> C_tensor;
    C_tensor.resize(CW_size, CPsi_size, CW_size);
    const fvMesh& mesh = L_U_SUPmodes[0].mesh();   
    /* volVectorField Uflux = _Uflux();
    volVectorField Vflux = _Uflux();
    volVectorField temp = _temp();
    surfaceScalarField& phi = _phi();  */

    for (label i = 0; i < CW_size; i++)
    {
        for (label j = 0; j < CPsi_size; j++)
        {   
            for (label k = 0; k < CW_size; k++)
            {
                if (fluxMethod == "consistent")
                {
                    C_tensor(i, j, k) = fvc::domainIntegrate(W_proj_modes[i] & fvc::div(
                                            fvc::flux(fvc::curl(Psi_zmodes[j])),
                                            W_proj_modes[k])).value();
                }
                else
                {
                    C_tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::div(
                                            linearInterpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                            L_U_SUPmodes[k])).value();
                }
            }
        }
    }

    /* if (Pstream::parRun())
    {
        reduce(C_tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    if (Pstream::master())
    {
        // Export the tensor
        ITHACAstream::SaveDenseTensor(C_tensor, "./ITHACAoutput/Matrices/",
                                      "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                          NSUPmodes) + "_t");
    }
  */
    return C_tensor;
}

Eigen::MatrixXd steadyNSWPsi::massW_term(label NWmodes, label NPsi_zmodes)
//modified
{
    label MW_size = NWmodes;
    Eigen::MatrixXd MW_matrix(MW_size, MW_size);

    // Project everything
    for (label i = 0; i < MW_size; i++)
    {
        for (label j = 0; j < MW_size; j++)
        {
            MW_matrix(i, j) = fvc::domainIntegrate(W_proj_modes[i] &
                                                  W_proj_modes[j]).value();
        }
    }

     /*if (Pstream::parRun())
    {
        reduce(MW_matrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(MW_matrix, "./ITHACAoutput/Matrices/",
                                      "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    } 
 */
    return MW_matrix;
}

Eigen::MatrixXd steadyNSWPsi::massPsi_z_term(label NWmodes, label NPsi_zmodes)
{
    label MW_size = NWmodes;
    label MPsi_size = NPsi_zmodes;
    Eigen::MatrixXd MPsi_matrix(MPsi_size, MW_size);

    // Project everything
    for (label i = 0; i < MPsi_size; i++)
    {
        for (label j = 0; j < MPsi_size; j++)
        {
            MPsi_matrix(i, j) = fvc::domainIntegrate(Psi_zmodes[i] &
                                                  W_proj_modes[j]).value();
        }
    }

    /* if (Pstream::parRun())
    {
        reduce(MPsi_matrix, sumOp<Eigen::MatrixXd>());
    }
 */
    /* if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(MPsi_matrix, "./ITHACAoutput/Matrices/",
                                      "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    }  */

    return MPsi_matrix;
}
