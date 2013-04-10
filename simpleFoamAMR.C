/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    simpleFoamAMR

Description
    Steady-state solver for incompressible, turbulent flow with Automatic 
    Mesh Refinement for Tetrahedral Meshes
    
    NOTE: Highly Experimental!!!!

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "dynamicFvMesh.H"

#include "errorEstimate.H"
#include "resError.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"
#       include "initConvergenceCheck.H"

        // ----- Begin - Dynamic part of the solver ------- //
        fvc::makeAbsolute(phi, U);

        bool meshChanged = mesh.update();

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);
        // ----- End - Dynamic part of the solver --------- //

        p.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        {
#           include "UEqn.H"
#           include "pEqn.H"
        }

        turbulence->correct();

        errorEstimate<vector> ee
        (
            resError::div(phi,U)
          - resError::laplacian(turbulence->nuEff(),U)
          - fvc::div(turbulence->nuEff()*dev(fvc::grad(U)().T()))
         ==
          -fvc::grad(p)
        );

        errField = mag(ee.error());
        pGrad = mag(fvc::grad(p));

        runTime.write();

        Info<< endl
            << " errField Max = " << max(errField)
            << " errField Min = " << min(errField)
            << nl << endl;

        Info<< " pGrad Max = " << max(pGrad)
            << " pGrad Min = " << min(pGrad)
            << nl << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

#       include "convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
