/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2023 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "internal/XtbForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <cmath>

using namespace XtbPlugin;
using namespace OpenMM;
using namespace std;

XtbForceImpl::XtbForceImpl(const XtbForce& owner) : CustomCPPForceImpl(owner), owner(owner), env(nullptr), calc(nullptr), res(nullptr), mol(nullptr), hasInitializedMolecule(false) {
}

XtbForceImpl::~XtbForceImpl() {
    if (res != nullptr)
        xtb_delResults(&res);
    if (calc != nullptr)
        xtb_delCalculator(&calc);
    if (mol != nullptr)
        xtb_delMolecule(&mol);
    if (env != nullptr)
        xtb_delEnvironment(&env);
}

void XtbForceImpl::initialize(ContextImpl& context) {
    CustomCPPForceImpl::initialize(context);
    indices = owner.getParticleIndices();
    numbers = owner.getAtomicNumbers();
    if (indices.size() != numbers.size())
        throw OpenMMException("Different numbers of particle indices and atomic numbers are specified");
    charge = owner.getCharge();
    multiplicity = owner.getMultiplicity();
    electrostaticEmbedding = owner.hasElectrostaticEmbedding();
    if (electrostaticEmbedding) {
        pointCharges = owner.getPointCharges();
        qmParticleIndices = owner.getQMParticleIndices();
        pcCutoff2 = std::pow(owner.getPointChargeCutoff(), 2);

        if (owner.getMethod() == XtbForce::GFNFF) {
            throw OpenMMException("GFNFF method does not support external charges");
        }
    }
    env = xtb_newEnvironment();
    calc = xtb_newCalculator();
    res = xtb_newResults();
    xtb_setVerbosity(env, XTB_VERBOSITY_MUTED);
    checkErrors();
    int numParticles = indices.size();
    positionVec.resize(3*numParticles, 0.0);
    forceVec.resize(3*numParticles);
  }

double XtbForceImpl::computeForce(ContextImpl& context, const vector<Vec3>& positions, vector<Vec3>& forces) {
    const double distanceScale = 18.897261246257703; // Convert nm to bohr
    const double energyScale = 2625.4996394798254; // Convert Hartree to kJ/mol
    const double forceScale = 49614.75258920568; // Convert Hartree/bohr to kJ/mol/nm

    // Pass the current state to XTB.

    int numParticles = indices.size();
    for (int i = 0; i < numParticles; i++) {
        positionVec[3*i] = distanceScale*positions[indices[i]][0];
        positionVec[3*i+1] = distanceScale*positions[indices[i]][1];
        positionVec[3*i+2] = distanceScale*positions[indices[i]][2];
    }
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double boxVectors[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            boxVectors[3*i+j] = distanceScale*box[i][j];
    if (hasInitializedMolecule)
        xtb_updateMolecule(env, mol, positionVec.data(), boxVectors);
    else {
        bool periodic[3] = {owner.usesPeriodicBoundaryConditions(), owner.usesPeriodicBoundaryConditions(), owner.usesPeriodicBoundaryConditions()};
        mol = xtb_newMolecule(env, &numParticles, numbers.data(), positionVec.data(), &charge, &multiplicity, boxVectors, periodic);
        checkErrors();
        if (owner.getMethod() == XtbForce::GFN1xTB)
            xtb_loadGFN1xTB(env, mol, calc, NULL);
        else if (owner.getMethod() == XtbForce::GFN2xTB)
            xtb_loadGFN2xTB(env, mol, calc, NULL);
        else if (owner.getMethod() == XtbForce::GFNFF)
            xtb_loadGFNFF(env, mol, calc, NULL);
        checkErrors();
        hasInitializedMolecule = true;
    }
    checkErrors();

    if (electrostaticEmbedding) {
        boundaryPCIndices.clear();
        boundaryPCCharges.clear();
        boundaryPCNumbers.clear();
        boundaryPCPositions.clear();

        auto chargeGroupInBoundaryRegion = [&, periodic = owner.usesPeriodicBoundaryConditions()](const std::vector<XtbPointCharge>& chargeGroup) {
            for (auto [index, number, charge]: chargeGroup) {
                for (auto j: qmParticleIndices) {
                    // see https://github.com/openmm/openmm/blob/master/platforms/reference/src/SimTKReference/ReferenceForce.cpp#L80
                    // works only for cubic boxes
                    Vec3 diff = positions[index] - positions[j];
                    if (periodic) {
                        const Vec3 base(std::floor(diff[0] / box[0][0] + 0.5) * box[0][0],
                                        std::floor(diff[1] / box[1][1] + 0.5) * box[1][1],
                                        std::floor(diff[2] / box[2][2] + 0.5) * box[2][2]);
                        diff -= base;
                    }
                    if (diff.dot(diff) <= pcCutoff2) {
                        return true;
                    }

                }
            }
            return false;
        };


        for (const auto& chargeGroup : pointCharges) {
            if (chargeGroupInBoundaryRegion(chargeGroup)) {
                for (auto [index, number, charge]: chargeGroup) {
                    boundaryPCIndices.push_back(index);
                    boundaryPCNumbers.push_back(number);
                    boundaryPCCharges.push_back(charge);

                    boundaryPCPositions.push_back(distanceScale*positions[index][0]);
                    boundaryPCPositions.push_back(distanceScale*positions[index][1]);
                    boundaryPCPositions.push_back(distanceScale*positions[index][2]);
                }
            }
        }
        numBoundaryPC = boundaryPCIndices.size();
        pcForceVec.resize(3 * numBoundaryPC);

        xtb_setExternalCharges(env, calc, &numBoundaryPC, boundaryPCNumbers.data(), boundaryPCCharges.data(), boundaryPCPositions.data());
        checkErrors();
    }

    // Perform the computation.

    xtb_singlepoint(env, mol, calc, res);
    checkErrors();
    double energy;
    xtb_getEnergy(env, res, &energy);
    checkErrors();
    xtb_getGradient(env, res, forceVec.data());
    for (int i = 0; i < positions.size(); i++)
        forces[i] = Vec3();
    for (int i = 0; i < numParticles; i++)
        forces[indices[i]] = -forceScale*Vec3(forceVec[3*i], forceVec[3*i+1], forceVec[3*i+2]);
    if (electrostaticEmbedding) {
        xtb_getPCGradient(env, res, pcForceVec.data());
        for (int i = 0; i < numBoundaryPC; i++)
            forces[boundaryPCIndices[i]] = -forceScale*Vec3(pcForceVec[3*i], pcForceVec[3*i+1], pcForceVec[3*i+2]);
        xtb_releaseExternalCharges(env, calc);
    }
    return energyScale*energy;
}

void XtbForceImpl::checkErrors() {
    if (xtb_checkEnvironment(env)) {
        vector<char> buffer(1000);
        int maxLength = buffer.size();
        xtb_getError(env, buffer.data(), &maxLength);
        throw OpenMMException(string(buffer.data()));
    }
}