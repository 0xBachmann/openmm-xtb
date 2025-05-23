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

#include "openmm/OpenMMException.h"
#include "XtbForce.h"
#include "internal/XtbForceImpl.h"

using namespace XtbPlugin;
using namespace OpenMM;
using namespace std;

XtbPointCharge::XtbPointCharge(int index, int number, double charge) : index(index), number(number), charge(charge) {
}


XtbForce::XtbForce(XtbForce::Method method, double charge, int multiplicity, bool periodic,
                   const vector<int> &particleIndices, const vector<int> &atomicNumbers) : method(method),
    charge(charge), multiplicity(multiplicity), periodic(periodic), electrostaticEmbedding(false),
    particleIndices(particleIndices), atomicNumbers(atomicNumbers) {
}

XtbForce::XtbForce(XtbForce::Method method, double charge, int multiplicity, bool periodic,
                   const std::vector<int> &particleIndices, const std::vector<int> &atomicNumbers,
                   const std::vector<std::vector<XtbPointCharge> > &pointCharges,
                   const std::vector<int> &qmParticleIndices, double pcCutoff,
                   int boundaryUpdateFrequency) : method(method), charge(charge), multiplicity(multiplicity),
                                                          periodic(periodic), electrostaticEmbedding(true),
                                                          particleIndices(particleIndices),
                                                          atomicNumbers(atomicNumbers), pointCharges(pointCharges),
                                                          qmParticleIndices(qmParticleIndices), pcCutoff(pcCutoff),
                                                          boundaryUpdateFrequency(boundaryUpdateFrequency) {
}


XtbForce::Method XtbForce::getMethod() const {
    return method;
}

void XtbForce::setMethod(XtbForce::Method method) {
    this->method = method;
}

double XtbForce::getCharge() const {
    return charge;
}

void XtbForce::setCharge(double charge) {
    this->charge = charge;
}

int XtbForce::getMultiplicity() const {
    return multiplicity;
}

void XtbForce::setMultiplicity(int multiplicity) {
    this->multiplicity = multiplicity;
}

const vector<int> &XtbForce::getParticleIndices() const {
    return particleIndices;
}

void XtbForce::setParticleIndices(const vector<int> &indices) {
    particleIndices = indices;
}

const vector<int> &XtbForce::getAtomicNumbers() const {
    return atomicNumbers;
}

void XtbForce::setAtomicNumbers(const vector<int> &numbers) {
    atomicNumbers = numbers;
}

const std::vector<std::vector<XtbPointCharge> > &XtbForce::getPointCharges() const {
    if (!electrostaticEmbedding) throw OpenMMException("No Point Charges have been added");
    return pointCharges;
}

void XtbForce::setPointCharges(const std::vector<std::vector<XtbPointCharge> > &pc) {
    pointCharges = pc;
}

const std::vector<int> &XtbForce::getQMParticleIndices() const {
    if (!electrostaticEmbedding) throw OpenMMException("No Point Charges have been added");
    return qmParticleIndices;
}

void XtbForce::setQMParticleIndices(const std::vector<int> &qmIndices) {
    qmParticleIndices = qmIndices;
}

double XtbForce::getPointChargeCutoff() const {
    if (!electrostaticEmbedding) throw OpenMMException("No Point Charges have been added");
    return pcCutoff;
}

void XtbForce::setPointChargeCutoff(double cutoff) {
    pcCutoff = cutoff;
}

bool XtbForce::hasElectrostaticEmbedding() const {
    return electrostaticEmbedding;
}

int XtbForce::getBoundaryUpdateFrequency() const {
    if (!electrostaticEmbedding) throw OpenMMException("No Point Charges have been added");
    return boundaryUpdateFrequency;
}

void XtbForce::setBoundaryUpdateFrequency(int frequency) {
    boundaryUpdateFrequency = frequency;
}

bool XtbForce::usesPeriodicBoundaryConditions() const {
    return true; // periodic boundary conditions are used but not for xtb_singlepoint
//    return periodic;
}

void XtbForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    this->periodic = periodic;
}

ForceImpl *XtbForce::createImpl() const {
    return new XtbForceImpl(*this);
}
