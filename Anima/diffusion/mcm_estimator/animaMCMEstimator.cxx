#include <animaMCMEstimatorImageFilter.h>
#include <animaGradientFileReader.h>
#include <animaVectorOperations.h>
#include <animaReadWriteFunctions.h>

#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

#include <fstream>

// Update progression of the process
void eventCallback(itk::Object *caller, const itk::EventObject &event, void *clientData)
{
    itk::ProcessObject *processObject = (itk::ProcessObject *)caller;
    std::cout << "\033[K\rProgression: " << (int)(processObject->GetProgress() * 100) << "%" << std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    // Required arguments
    TCLAP::ValueArg<std::string> dwiArg(
        "i", "dwi",
        "DWI 4D Volume",
        true, "", "DWI images", cmd);
    TCLAP::ValueArg<std::string> gradsArg(
        "g", "grads",
        "Gradient table",
        true, "", "gradients", cmd);
    TCLAP::ValueArg<std::string> bvalsArg(
        "b", "bvals",
        "B-value list",
        true, "", "b-values", cmd);
    TCLAP::SwitchArg bvalueScaleArg(
        "B", "b-no-scale",
        "Do not scale b-values according to gradient norm",
        cmd);
    TCLAP::ValueArg<double> smallDeltaArg(
        "", "small-delta",
        "Diffusion small delta (in seconds)",
        false, anima::DiffusionSmallDelta, "small delta", cmd);
    TCLAP::ValueArg<double> bigDeltaArg(
        "", "big-delta",
        "Diffusion big delta (in seconds)",
        false, anima::DiffusionBigDelta, "big delta", cmd);

    // Outputs
    TCLAP::ValueArg<std::string> outArg(
        "o", "out-mcm",
        "MCM output volume",
        true, "", "MCM output", cmd);
    TCLAP::ValueArg<std::string> aicArg(
        "a", "out-aic",
        "Output estimated AICu image",
        false, "", "AICu output image", cmd);
    TCLAP::ValueArg<std::string> outB0Arg(
        "", "out-b0",
        "Output estimated B0 image",
        false, "", "B0 output image", cmd);
    TCLAP::ValueArg<std::string> outSigmaArg(
        "", "out-sig",
        "Output estimated noise sigma square image",
        false, "", "noise sigma square output image", cmd);
    TCLAP::ValueArg<std::string> outMoseArg(
        "", "out-mose",
        "Output model selection map",
        false, "", "model selection output image", cmd);
    TCLAP::ValueArg<std::string> inMoseArg(
        "I", "in-mose",
        "Input model selection map (overrides model selection step)",
        false, "", "model selection input image", cmd);

    // Optional arguments
    TCLAP::ValueArg<std::string> computationMaskArg(
        "m", "mask",
        "Computation mask",
        false, "", "computation mask", cmd);
    TCLAP::ValueArg<double> b0thrArg(
        "", "b0-thr",
        "Background threshold on B0 value (default: 10)",
        false, 10.0, "B0 theshold", cmd);

    TCLAP::ValueArg<unsigned int> nbFasciclesArg(
        "n", "nb-fascicles",
        "Number of computed fascicles (default: 2)",
        false, 2, "number of fascicles", cmd);
    TCLAP::ValueArg<unsigned int> cylinderTypeArg(
        "c", "cylinder-type",
        "Compartment type for cylinders: 1: Stick, 2: Zeppelin, 3: Tensor, 4: NODDI, 5: DDI, 6: CHARMED (default: 3)",
        false, 3, "cylinder type", cmd);
    TCLAP::SwitchArg aicSelectNbCompartmentsArg(
        "M", "opt-nb-comp",
        "Activate AICC-based number of compartments selection",
        cmd, false);
    TCLAP::ValueArg<unsigned int> sphereTypeArg(
        "s", "sphere-type",
        "Compartment type for spheres: 0: none, 1: SphereGPDPulsedGradient, 2: PlaneSGPPulsedGradient (default: 1)",
        false, 3, "sphere type", cmd);

    TCLAP::SwitchArg freeWaterCompartmentArg(
        "F", "free-water",
        "Model with free water",
        cmd, false);
    
    TCLAP::SwitchArg smallGlialCellCompartmentArg(
        "", "with-small-glial-cells",
        "Model with small glial cells",
        cmd, false);
    
    TCLAP::SwitchArg mediumGlialCellCompartmentArg(
        "", "with-medium-glial-cells",
        "Model with medium glial cells",
        cmd, false);
    
    TCLAP::SwitchArg largeGlialCellCompartmentArg(
        "", "with-large-glial-cells",
        "Model with large glial cells",
        cmd, false);
    
    TCLAP::SwitchArg optSphereRadiusArg(
        "", "opt-sphere-radius",
        "Optimize sphere radius value",
        cmd, false);
    TCLAP::SwitchArg optSphereDiffArg(
        "", "opt-sphere-diff",
        "Optimize sphere diffusivity value",
        cmd, false);
    TCLAP::SwitchArg optCylinderRadiusArg(
        "", "opt-cylinder-radius",
        "Optimize cylinder radius value",
        cmd, false);

    TCLAP::SwitchArg fixDiffArg(
        "", "fix-diff",
        "Fix diffusivity value",
        cmd, false);
    TCLAP::SwitchArg fixKappaArg(
        "", "fix-kappa",
        "Fix orientation concentration values",
        cmd, false);
    TCLAP::SwitchArg fixEAFArg(
        "", "fix-eaf",
        "Fix extra axonal fraction values",
        cmd, false);

    TCLAP::SwitchArg commonDiffusivitiesArg(
        "", "common-diffusivities",
        "Share diffusivity values among compartments",
        cmd, false);
    TCLAP::SwitchArg commonKappaArg(
        "", "common-kappa",
        "Share orientation concentration values among compartments",
        cmd, false);
    TCLAP::SwitchArg commonEAFArg(
        "", "common-eaf",
        "Share extra axonal fraction values among compartments",
        cmd, false);

    // Initial values for sphere compartment
    TCLAP::ValueArg<double> initSphereDiffArg(
        "", "init-sphere-diff",
        "Initial sphere diffusivity in mm2/s (default: 1.71e-3 mm2/s)",
        false, 1.71e-3, "initial sphere diffusivity", cmd);
    TCLAP::ValueArg<double> initSphereRadiusArg(
        "", "init-sphere-radius",
        "Initial sphere radius in mm (default: 0.015 mm)",
        false, 0.015, "initial sphere radius", cmd);
    
    // Initial values for cylinder compartments
    TCLAP::ValueArg<double> initAxialDiffArg(
        "", "init-axial-diff",
        "Initial axial diffusivity (default: 1.71e-3)",
        false, 1.71e-3, "initial axial diffusivity", cmd);
    TCLAP::ValueArg<double> initRadialDiff1Arg(
        "", "init-radial-diff1",
        "Initial first radial diffusivity (default: 1.9e-4)",
        false, 1.9e-4, "initial first radial diffusivity", cmd);
    TCLAP::ValueArg<double> initRadialDiff2Arg(
        "", "init-radial-diff2",
        "Initial second radial diffusivity (default: 1.5e-4)",
        false, 1.5e-4, "initial second radial diffusivity", cmd);
    TCLAP::ValueArg<double> initCylinderRadiusArg(
        "", "init-cylinder-radius",
        "Initial cylinder radius (default: 0.0005)",
        false, 0.0005, "initial cylinder radius", cmd);

    // Optimization parameters
    TCLAP::ValueArg<std::string> optimizerArg(
        "", "optimizer",
        "Optimizer for estimation: bobyqa (default), ccsaq, bfgs or levenberg",
        false, "bobyqa", "optimizer", cmd);
    TCLAP::ValueArg<double> absCostChangeArg(
        "", "abs-cost-change",
        "Cost function change to stop estimation (default: 0.01)",
        false, 0.01, "cost change threshold", cmd);
    TCLAP::ValueArg<unsigned int> mlModeArg(
        "", "ml-mode",
        "ML estimation strategy: marginal likelihood (0), profile likelihood (1, default), Variable projection (2)",
        false, 1, "ML mode", cmd);
    TCLAP::ValueArg<double> xTolArg(
        "x", "x-tol",
        "Tolerance for relative position in optimization (default: 0 -> 1.0e-4 or 1.0e-7 for bobyqa)",
        false, 0, "position relative tolerance", cmd);
    TCLAP::ValueArg<double> fTolArg(
        "", "f-tol",
        "Tolerance for relative cost in optimization (default: 0 -> function of position tolerance)",
        false, 0, "cost relative tolerance", cmd);
    TCLAP::ValueArg<unsigned int> maxEvalArg(
        "e", "max-eval",
        "Maximum evaluations (default: 0 -> function of number of unknowns)",
        false, 0, "max evaluations", cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg(
        "T", "nb-threads",
        "Number of threads to run on (default: all cores)",
        false, itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using FilterType = anima::MCMEstimatorImageFilter<double, double>;
    using InputImageType = FilterType::InputImageType;
    using MaskImageType = FilterType::MaskImageType;
    using FilterPointer = FilterType::Pointer;

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    FilterPointer filter = FilterType::New();

    filter->SetUseConstrainedOrientationConcentration(fixKappaArg.isSet());
    if (!fixKappaArg.isSet())
        filter->SetUseCommonConcentrations(commonKappaArg.isSet());
    else
        filter->SetUseCommonConcentrations(false);

    filter->SetUseConstrainedExtraAxonalFraction(fixEAFArg.isSet());
    if (!fixEAFArg.isSet())
        filter->SetUseCommonExtraAxonalFractions(commonEAFArg.isSet());
    else
        filter->SetUseCommonExtraAxonalFractions(false);

    std::cout << "Loading input DWI image..." << std::endl;

    anima::setMultipleImageFilterInputsFromFileName<InputImageType, FilterType>(dwiArg.getValue(), filter);

    // Load gradient table and b-value list
    std::cout << "Importing gradient table and b-values..." << std::endl;

    using GFReaderType = anima::GradientFileReader<vnl_vector_fixed<double, 3>, double>;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradsArg.getValue());
    gfReader.SetBValueBaseString(bvalsArg.getValue());
    gfReader.SetGradientIndependentNormalization(bvalueScaleArg.isSet());
    gfReader.SetSmallDelta(smallDeltaArg.getValue());
    gfReader.SetBigDelta(bigDeltaArg.getValue());
    gfReader.Update();

    GFReaderType::GradientVectorType directions = gfReader.GetGradients();
    for (unsigned int i = 0; i < directions.size(); ++i)
        filter->AddGradientDirection(i, directions[i]);

    GFReaderType::BValueVectorType mb = gfReader.GetGradientStrengths();

    filter->SetGradientStrengths(mb);
    filter->SetSmallDelta(smallDeltaArg.getValue());
    filter->SetBigDelta(bigDeltaArg.getValue());

    if (computationMaskArg.getValue() != "")
        filter->SetComputationMask(anima::readImage<MaskImageType>(computationMaskArg.getValue()));

    if (inMoseArg.getValue() != "")
        filter->SetMoseVolume(anima::readImage<FilterType::MoseImageType>(inMoseArg.getValue()));

    filter->SetB0Threshold(b0thrArg.getValue());

    filter->SetModelWithFreeWaterComponent(freeWaterCompartmentArg.isSet());
    filter->SetModelWithSmallGlialCellComponent(smallGlialCellCompartmentArg.isSet());
    filter->SetModelWithMediumGlialCellComponent(mediumGlialCellCompartmentArg.isSet());
    filter->SetModelWithLargeGlialCellComponent(largeGlialCellCompartmentArg.isSet());

    bool modelWithSphereComponent = true;
    switch(sphereTypeArg.getValue())
    {
        case 0:
        default:
        modelWithSphereComponent = false;
        break;

        case 1:
        filter->SetSphereCompartmentType(anima::SphereGPDPulsedGradient);
        break;

        case 2:
        filter->SetSphereCompartmentType(anima::PlaneSGPPulsedGradient);
        break;
    }
    filter->SetModelWithSphereComponent(modelWithSphereComponent);

    switch (cylinderTypeArg.getValue())
    {
    case 1:
        filter->SetCylinderCompartmentType(anima::Stick);
        break;

    case 2:
        filter->SetCylinderCompartmentType(anima::Zeppelin);
        break;

    case 3:
        filter->SetCylinderCompartmentType(anima::Tensor);
        break;

    case 4:
        filter->SetCylinderCompartmentType(anima::NODDI);
        break;

    case 5:
        filter->SetCylinderCompartmentType(anima::DDI);
        break;

    case 6:
        filter->SetCylinderCompartmentType(anima::CHARMED);
        break;

    default:
        std::cerr << "Unsupported compartment type" << std::endl;
        return EXIT_FAILURE;
    }

    filter->SetAxialDiffusivityValue(initAxialDiffArg.getValue());
    filter->SetRadialDiffusivity1Value(initRadialDiff1Arg.getValue());
    filter->SetRadialDiffusivity2Value(initRadialDiff2Arg.getValue());
    filter->SetSphereDiffusivityValue(initSphereDiffArg.getValue());
    filter->SetSphereRadiusValue(initSphereRadiusArg.getValue());
    filter->SetCylinderRadiusValue(initCylinderRadiusArg.getValue());

    filter->SetNumberOfCompartments(nbFasciclesArg.getValue());
    filter->SetFindOptimalNumberOfCompartments(aicSelectNbCompartmentsArg.isSet());

    filter->SetOptimizer(optimizerArg.getValue());
    filter->SetAbsoluteCostChange(absCostChangeArg.getValue());

    filter->SetNoiseType(FilterType::Gaussian);
    filter->SetMLEstimationStrategy((FilterType::MaximumLikelihoodEstimationMode)mlModeArg.getValue());

    filter->SetXTolerance(xTolArg.getValue());
    filter->SetFTolerance(fTolArg.getValue());
    filter->SetMaxEval(maxEvalArg.getValue());

    filter->SetUseConstrainedDiffusivity(fixDiffArg.isSet());
    filter->SetUseConstrainedSphereDiffusivity(!optSphereDiffArg.isSet());
    filter->SetUseConstrainedSphereRadius(!optSphereRadiusArg.isSet());
    filter->SetUseConstrainedCylinderRadius(!optCylinderRadiusArg.isSet());

    if (!fixDiffArg.isSet())
        filter->SetUseCommonDiffusivities(commonDiffusivitiesArg.isSet());
    else
        filter->SetUseCommonDiffusivities(false);

    filter->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    filter->AddObserver(itk::ProgressEvent(), callback);

    std::cout << "Estimating MCM..." << std::endl;

    itk::TimeProbe tmpTimer;
    tmpTimer.Start();

    try
    {
        filter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    tmpTimer.Stop();

    std::cout << "\nEstimation done in " << tmpTimer.GetTotal() << " s" << std::endl;
    std::cout << "Writing MCM to: " << outArg.getValue() << std::endl;

    try
    {
        filter->WriteMCMOutput(outArg.getValue());
    }
    catch (std::exception &e)
    {
        std::cerr << "Error writing output " << e.what() << std::endl;
    }

    if (aicArg.getValue() != "")
    {
        std::cout << "Writing AICu image to: " << aicArg.getValue() << std::endl;
        anima::writeImage(aicArg.getValue(), filter->GetAICcVolume());
    }

    if (outB0Arg.getValue() != "")
    {
        std::cout << "Writing B0 image to: " << outB0Arg.getValue() << std::endl;
        anima::writeImage(outB0Arg.getValue(), filter->GetB0Volume());
    }

    if (outSigmaArg.getValue() != "")
    {
        std::cout << "Writing noise sigma square image to: " << outSigmaArg.getValue() << std::endl;
        anima::writeImage(outSigmaArg.getValue(), filter->GetSigmaSquareVolume());
    }

    if (outMoseArg.getValue() != "")
    {
        std::cout << "Writing model selection image to: " << outMoseArg.getValue() << std::endl;
        anima::writeImage(outMoseArg.getValue(), filter->GetMoseVolume());
    }

    return EXIT_SUCCESS;
}
