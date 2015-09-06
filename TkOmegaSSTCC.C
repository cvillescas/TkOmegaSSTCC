/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/
// TkOmegaSSTCC-incompressible
#include "TkOmegaSSTCC.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"
#include "cellZoneMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(TkOmegaSSTCC, 0);
addToRunTimeSelectionTable(RASModel, TkOmegaSSTCC, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TkOmegaSSTCC::TkOmegaSSTCC
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName

)
:
    kOmegaSST(U, phi, transport, turbulenceModelName, modelName), 

    cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr1",
            coeffDict_,
            1.0
        )
    ),
    cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr2",
            coeffDict_,
            2.0
        )
    ),
    cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr3",
            coeffDict_,
            1.0
        )
    )
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);
    nut_ =
    (
        a1_*k_
      / max
        (
            a1_*omega_,
            b1_*F23()*sqrt(2.0)*mag(symm(fvc::grad(U_)))
        )
    );
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void TkOmegaSSTCC::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    // Compute rStar
    tmp<volTensorField> tSkew = skew(fvc::grad(U_)); 
    tmp<volSymmTensorField> tSymm = symm(fvc::grad(U_));
    volScalarField symInnerProduct(2.0*tSymm() && tSymm()); 
    volScalarField asymInnerProduct
    (
        max(2.0*tSkew() && tSkew(),
        dimensionedScalar("0", dimensionSet(0, 0, -2, 0, 0), 0.0))
    );
    volScalarField w
    (
        atan(dimensionedScalar("4",dimensionSet(0,0,2,0,0),1.0e-02)*asymInnerProduct)*2.0/(constant::mathematical::pi)*(asymInnerProduct-symInnerProduct)+symInnerProduct
    );
    volScalarField rStar(sqrt(symInnerProduct/w));


    // Compute rTilda 
    volScalarField D(sqrt(max(symInnerProduct, 0.09*omega_*omega_)));
    tmp<volSymmTensorField> divS =
    (
        fvc::ddt(tSymm())
       +fvc::div
        (
            phi_, tSymm()
        )
    );
    volScalarField rT((tSkew().T() & tSymm) && divS);

    divS.clear();
    tSkew.clear();
    tSymm.clear();
    volScalarField w2
    (
        atan(dimensionedScalar("1",dimensionSet(0,0,2,0,0),1.0e-2)*asymInnerProduct)*2.0/(constant::mathematical::pi)*(sqrt(asymInnerProduct)-D)+D //T
    );
    volScalarField rTilda(2.0*rT/w2/D/D/D);

    // Compute Frot
    volScalarField Frot
    (
        IOobject
        (
            "Frot",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        max
        (
            min
            (
                (scalar(1.0)+cr1_)*2.0*rStar/(scalar(1)+rStar)*(scalar(1.0)
               -cr3_*atan(cr2_*rTilda))-cr1_,
               scalar(1.25)
            ),
        scalar(0.0))
    );
    rStar.clear();
    rTilda.clear();
    rT.clear();
    D.clear();

    // Compute kOmegaSST variables
    volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    const volScalarField G(type() + ".G", nut_*S2);

    omega_.boundaryField().updateCoeffs();
    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );
    volScalarField F1(this->F1(CDkOmega));


    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)*S2*Frot
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();
    omegaEqn().boundaryManipulate(omega_.boundaryField());
    solve(omegaEqn);
    bound(omega_, omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G*Frot, c1_*betaStar_*k_*omega_)
      - fvm::Sp(betaStar_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
