/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    initbubble

Description
    Initalization utility for example case

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include <iostream>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    volScalarField alpha_water
    (
        IOobject
        (
            "alpha.water",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField alpha_air
    (
        IOobject
        (
            "alpha.air",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField alpha_oil
    (
        IOobject
        (
            "alpha.oil",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField alphas
    (
        IOobject
        (
            "alphas",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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

    IOdictionary thermophysicalProperties_oil
    (
        IOobject
        (
            "thermophysicalProperties.oil",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

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

    dimensionedScalar B_oil	("B", 		dimPressure, 	thermophysicalProperties_oil.subDict("mixture").subDict("equationOfState"));
    dimensionedScalar rho0_oil	("rho0", 	dimDensity, 	thermophysicalProperties_oil.subDict("mixture").subDict("equationOfState"));
    dimensionedScalar p0_oil	("p0", 	dimPressure, 	thermophysicalProperties_oil.subDict("mixture").subDict("equationOfState"));
    dimensionedScalar gamma_oil	("gamma", 	dimless, 	thermophysicalProperties_oil.subDict("mixture").subDict("equationOfState"));

    dimensionedScalar B_air	("B", 		dimPressure, 	thermophysicalProperties_air.subDict("mixture").subDict("equationOfState"));
    dimensionedScalar rho0_air	("rho0", 	dimDensity, 	thermophysicalProperties_air.subDict("mixture").subDict("equationOfState"));
    dimensionedScalar p0_air	("p0", 	dimPressure, 	thermophysicalProperties_air.subDict("mixture").subDict("equationOfState"));
    dimensionedScalar gamma_air	("gamma", 	dimless, 	thermophysicalProperties_air.subDict("mixture").subDict("equationOfState"));

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

    Info<< "Reading \n" << endl;
    
    dimensionedScalar epsDtB_air("epsDtB_air", dimLength*dimLength, transportProperties);
    
    dimensionedScalar epsDtB_water("epsDtB_water", dimLength*dimLength, transportProperties);

    dimensionedVector BubbleCenter("BubbleCenter", dimLength, transportProperties);

    dimensionedScalar BubbleRadius("BubbleRadius", dimLength, transportProperties);

    dimensionedVector dropCenter("dropCenter", dimLength, transportProperties);

    dimensionedScalar dropRadius("dropRadius", dimLength, transportProperties);

    dimensionedScalar epsionBubble("epsionBubble", dimless, transportProperties);

    dimensionedScalar pliquid("pliquid", dimPressure, transportProperties);

    dimensionedScalar pBubble("pBubble", dimPressure, transportProperties);

    dimensionedScalar rho = rho0_oil*pow((pliquid+B_oil)/(p0_oil+B_oil),1/gamma_oil);
    pBubble = (p0_air + B_air)*pow(rho/rho0_air,gamma_air)-B_air;
	
    forAll(alpha_water, cellI)
    {
    	vector x = mesh.C()[cellI];

        scalar kb = magSqr((x[0]-BubbleCenter.value()[0])/BubbleRadius.value())+magSqr((x[1]-BubbleCenter.value()[1])/BubbleRadius.value())+magSqr((x[2]-BubbleCenter.value()[2])/BubbleRadius.value());

        alpha_air[cellI] = (1-Foam::tanh((kb-1)/epsionBubble.value()))/2;
    	
        scalar kb2 = magSqr((x[0]-dropCenter.value()[0])/dropRadius.value())+magSqr((x[1]-dropCenter.value()[1])/dropRadius.value())+magSqr((x[2]-dropCenter.value()[2])/dropRadius.value());

        alpha_water[cellI] = (1-Foam::tanh((kb2-1)/epsionBubble.value()))/2;

        inkMap[cellI] = x - BubbleCenter.value();
    }


//smear the interface to avoid the numerical instability grows in the beginning

    volScalarField alpha_air0(alpha_air);
    fvScalarMatrix alpha_airSEqn
    (
        fvm::Sp(scalar(1),alpha_air) - fvm::laplacian(epsDtB_air,alpha_air) == alpha_air0
    );
    alpha_airSEqn.solve();

    if (epsDtB_water.value() != 0)
    {
        volScalarField alpha_water0(alpha_water);
        fvScalarMatrix alpha_waterSEqn
        (
            fvm::Sp(scalar(1),alpha_water) - fvm::laplacian(epsDtB_water,alpha_water) == alpha_water0
        );
        alpha_waterSEqn.solve();
    }

    alpha_oil = max(min(1.0-alpha_air-alpha_water,1.0),0.0);
    
    alphas = alpha_oil + 2*alpha_air;

    p = pliquid*(1-alpha_air) + pBubble*alpha_air;

    p_rgh = p;
         
    alpha_water.write();
    
    alpha_air.write();
    
    alpha_oil.write();
    
    alphas.write();

    p_rgh.write();

    p.write();
    
    inkMap.write();
    
    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
