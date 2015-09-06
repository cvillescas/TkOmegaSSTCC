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

Class
    Foam::incompressible::RASModels::TkOmegaSSTCC

Description
    Implementation of the k-omega-SST turbulence model for compressible flows.

    Turbulence model described in:
    \verbatim
        Menter, F., Kuntz, M. and Langtry. R.
        "Ten years of Industrial Experience with the SST Turbulence Model"
        Turbulence, Heat and Mass Transfer 4.
        2003.
    \endverbatim

    With curvature correction term, described in:
    \verbatim
        Shur, M. L., Strelets, M. K., Travin, A. K., Spalart, P. R.
        "Turbulence Modeling in Rotating and Curved Channels:
        Assessing the Spalart-Shur Correction"
        AIAA Journal Vol. 38, No. 5, 2000, pp. 784-792.
        http://arc.aiaa.org/doi/abs/10.2514/2.1058
    \endverbatim

    and reformulated according to the SST model in:
    \verbatim
        Menter F.R. Smirnov P.E.
        "Sensitization of the sst turbulence model to rotation and curvature
        by applying the spalart–shur correction term"
        Journ. of Turbomachinery, 131(4), 2009.
    \endverbatim

    Note that this implementation is written in terms of alpha diffusion
    coefficients rather than the more traditional sigma (alpha = 1/sigma) so
    that the blending can be applied to all coefficuients in a consistent
    manner.  The paper suggests that sigma is blended but this would not be
    consistent with the blending of the k-epsilon and k-omega models.

    Also note that the error in the last term of equation (2) relating to
    sigma has been applied.

    Wall-functions are applied in this implementation by using equations (14)
    to specify the near-wall omega as appropriate.

    The blending functions (15) and (16) are not currently used because of the
    uncertainty in their origin, range of applicability and that is y+ becomes
    sufficiently small blending u_tau in this manner clearly becomes nonsense.

    The default model coefficients correspond to the following:
    \verbatim
        bTkOmegaSSTCCCoeffs
        {
            alphaK1     0.85034;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.85616;
            Prt         1.0;    // only for compressible
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      0.5532;
            gamma2      0.4403;
            a1          0.31;
            c1          10.0;
            cr1         1.0;
            cr2         2.0;
            cr3         1.0;
        }
    \endverbatim

SourceFiles
    TkOmegaSSTCC.C
\*---------------------------------------------------------------------------*/

#ifndef incompressibleTkOmegaSSTCC_H //?
#define incompressibleTkOmegaSSTCC_H //?

#include "kOmegaSST.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class bTkOmegaSSTCC Declaration
\*---------------------------------------------------------------------------*/

class TkOmegaSSTCC
:
    public kOmegaSST
{

protected:

    // Protected data

        // Model coefficients
            dimensionedScalar cr1_;
            dimensionedScalar cr2_;
            dimensionedScalar cr3_;

public:

    //- Runtime type information
    TypeName("TkOmegaSSTCC");


    // Constructors

		//- Construct from components
		TkOmegaSSTCC
		(
		    const volVectorField& U,
		    const surfaceScalarField& phi,
                    transportModel& transport,
		    const word& turbulenceModelName = turbulenceModel::typeName,
		    const word& modelName = typeName
		);

    //- Destructor
    virtual ~TkOmegaSSTCC()
    {}


	// Member Functions

	    //- Solve the turbulence equations and correct the turbulence viscosity
		virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
