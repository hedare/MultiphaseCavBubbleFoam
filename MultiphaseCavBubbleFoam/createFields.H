Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p_temp
(
    IOobject
    (
        "p_temp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p_rgh
);

volScalarField p_g
(
    IOobject
    (
        "p_g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p_rgh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Constructing multiphaseMixtureThermo\n" << endl;
multiphaseMixtureThermo mixture(U, phi);

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoMin("rhoMin", dimDensity, transportProperties);

dimensionedScalar rhoMax("rhoMax", dimDensity, transportProperties);

dimensionedScalar pFact("pFact", dimless, transportProperties);

IOdictionary thermophysicalProperties_air
(
    IOobject
    (
        "thermophysicalProperties.air",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
dimensionedScalar B_air	("B", 		dimPressure, 	thermophysicalProperties_air.subDict("mixture").subDict("equationOfState"));
dimensionedScalar rho0_air	("rho0", 	dimDensity, 	thermophysicalProperties_air.subDict("mixture").subDict("equationOfState"));
dimensionedScalar p0_air	("p0", 	dimPressure, 	thermophysicalProperties_air.subDict("mixture").subDict("equationOfState"));
dimensionedScalar gamma_air	("gamma", 	dimless, 	thermophysicalProperties_air.subDict("mixture").subDict("equationOfState"));
bool correctBubbleMass = mesh.solutionDict().subDict("mysolver").lookupOrDefault<Switch>("correctBubbleMass", false);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mixture.rho()
);

volScalarField rho_water
(
    IOobject
    (
        "rho_water",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);

volScalarField rho_oil
(
    IOobject
    (
        "rho_oil",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);

volScalarField rho_air
(
    IOobject
    (
        "rho_air",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);

volVectorField inkMap
(
    IOobject
    (
        "inkMap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar pMin("pMin", dimPressure, mixture);

dimensionedScalar pMinGas("pMinGas", dimPressure, mixture);

mesh.setFluxRequired(p_rgh.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// Construct compressible turbulence model
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        mixture.rhoPhi(),
        mixture
    )
);

#include "createK.H"

const word kinematicCloudName
(
    args.getOrDefault<word>("cloud", "kinematicCloud")
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.nu()*rho
);

Info << "Constructing kinematicCloud " << kinematicCloudName << endl;
basicKinematicCloud kinematicCloud
(
    kinematicCloudName,
    rho,
    U,
    mu,
    g
);
