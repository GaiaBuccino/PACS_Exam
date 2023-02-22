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

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

Info<< "Reading field W\n" << endl;
volScalarField W
(
    IOobject
    (
        "W",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Psi_z\n" << endl;
volScalarField Psi_z
(
    IOobject
    (
        "Psi_z",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField temp
    (
    	IOobject
		(
			"temp",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
	    mesh,
		vector(0,0,1)
    );


Info<< "Reading vector field Psi\n" << endl;
volVectorField Psi
(
    IOobject
    (
        "Psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,   //si puo' mettere no_read??
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading vector field U\n" << endl;
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

//Remember to compute U = rot Psi
#include "createPhi.H"


