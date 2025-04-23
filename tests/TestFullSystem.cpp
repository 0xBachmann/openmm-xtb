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

/**
 * This tests the Reference implementation of XtbForce.
 */

#include "XtbForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/LangevinMiddleIntegrator.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/reference/SimTKOpenMMRealType.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/CustomIntegrator.h"
#include <istream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <format>
#include <ranges>
#include <memory>

using namespace XtbPlugin;
using namespace OpenMM;
using namespace std;

std::vector<OpenMM::Vec3> extractPositionsFromPDB(const std::string &filename) {
    std::vector<OpenMM::Vec3> positions;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open PDB file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        // Check if this is an ATOM or HETATM line
        if (line.length() >= 54 && (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM")) {
            try {
                // PDB format: columns 31-38 = x, 39-46 = y, 47-54 = z (in Angstroms)
                double x = std::stod(line.substr(30, 8));
                double y = std::stod(line.substr(38, 8));
                double z = std::stod(line.substr(46, 8));

                // Convert from Angstroms to nanometers (PDB uses Angstroms, OpenMM uses nm)
                x *= 0.1;
                y *= 0.1;
                z *= 0.1;

                positions.push_back(OpenMM::Vec3(x, y, z));
            } catch (const std::exception &e) {
                std::cerr << "Warning: Failed to parse line: " << line << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                // Continue processing other lines
            }
        }
    }

    file.close();

    if (positions.empty()) {
        throw std::runtime_error("No valid positions found in PDB file: " + filename);
    }

    return positions;
}

CustomIntegrator qmedsIntegrator(size_t num_endstates) {
    const double dt = 0.0005;
    CustomIntegrator integrator(dt);

    integrator.addGlobalVariable("a", std::exp(-1 * dt));
    integrator.addGlobalVariable("b", std::sqrt(1 - std::exp(-2 * 1 * dt)));
    integrator.addGlobalVariable("kT", 1 / 0.40339545545103495);
    integrator.addGlobalVariable("beta", 0.40339545545103495);

    integrator.addPerDofVariable("x1", 0);
    integrator.addPerDofVariable("sigma", 0);
    integrator.beginIfBlock("sigma = 0");
    integrator.addComputePerDof("sigma", "sqrt(kT/m)");
    integrator.endBlock();

    integrator.addUpdateContextState();

    // Add global parameters for reporting;
    integrator.addGlobalVariable("VR", 1.0);
    // integrator.addGlobalVariable("WR", 1.0)
    for (size_t i = 0; i < num_endstates; ++i) {
        integrator.addGlobalVariable(std::format("V{}_vac", i), 1.0);
        integrator.addComputeGlobal(std::format("V{}_vac", i), std::format("energy{}", num_endstates + 1 + i));
        integrator.addGlobalVariable(std::format("V{}_ee", i), 1.0);
        integrator.addComputeGlobal(std::format("V{}_ee", i), std::format("energy{}", 2 * num_endstates + 1 + i));
        integrator.addGlobalVariable(std::format("V{}_vdw", i), 1.0);
        integrator.addComputeGlobal(std::format("V{}_vdw", i), std::format("energy{}", 1 + i));
        integrator.addGlobalVariable(std::format("V{}", i), 1.0);
        integrator.addComputeGlobal(std::format("V{}", i), std::format("V{}_ee - V{}_vac + V{}_vdw", i, i, i));

        // integrator.addGlobalVariable(std::format("V{i}", 1.0)
        // integrator.addComputeGlobal(f"V{i}", f"energy{i + 1}");
        integrator.addGlobalVariable(std::format("v_eff{}", i), 1.0);
        integrator.addComputeGlobal(std::format("v_eff{}", i), std::format("V{}-eoff{}", i, i));

    }

    if (num_endstates < 2) {
        integrator.addComputeGlobal("VR", "V0 - eoff0");
        integrator.addGlobalVariable("scal_0", 1.0);
    } else {
        integrator.addGlobalVariable("part0", 1.0);
        integrator.addComputeGlobal("part0", "-beta * s * v_eff0");
        integrator.addGlobalVariable("part1", 1.0);
        integrator.addComputeGlobal("part1", "-beta * s * v_eff1");
        integrator.addGlobalVariable("maxpart", 1.0);
        integrator.addComputeGlobal("maxpart", "max(part0, part1)");


        for (size_t i = 2; i < num_endstates; ++i) {
            integrator.addGlobalVariable(std::format("part{}", i), 1.0);
            integrator.addComputeGlobal(std::format("part{}", i), std::format("-beta * s * v_eff{}", i));
            integrator.addComputeGlobal("maxpart", std::format("max(maxpart, part{})", i));
        }

        integrator.addGlobalVariable("logsumexp", 1.0);
        integrator.addComputeGlobal("logsumexp", std::format("maxpart + log({});", [&]() {
            std::string sum = "exp(part0-maxpart)";
            for (size_t i = 1; i < num_endstates; ++i) {
                sum += std::format("+exp(part{}-maxpart)", i);
            }
            return sum;
        }()));
        integrator.addComputeGlobal("VR", "-1/(beta*s) * logsumexp");

        for (size_t i = 0; i < num_endstates; ++i) {
            integrator.addGlobalVariable(std::format("scal_{}", i), 1.0);
            integrator.addComputeGlobal(std::format("scal_{}", i), std::format("exp(part{}-logsumexp)", i));
        }
    }
    // B A O A (or LFMiddle);
    integrator.addComputePerDof("v", "v + dt*f0/m");
    for (size_t i = 1; i < num_endstates; ++i) {
        // integrator.addComputePerDof("v", std::format("v + dt * scal_{i - 1} * f{i}/m");
        integrator.addComputePerDof("v", std::format("v + dt * f{} / m", num_endstates + i));   // V_vac forces
        integrator.addComputePerDof("v", std::format("v + dt * scal_{} * f{} / m", i - 1,
                                                     2 * num_endstates + i));    // V_ee forces
        integrator.addComputePerDof("v", std::format("v - dt * scal_{} * f{} / m", i - 1,
                                                     num_endstates + i));    // -V_vac forces
        integrator.addComputePerDof("v", std::format("v + dt * scal_{} * f{} / m", i - 1, i));   // V_vdW forces
    }
    integrator.addConstrainVelocities();
    integrator.addComputePerDof("x", "x + 0.5*dt*v");
    integrator.addComputePerDof("v", "a*v + b*sigma*gaussian");
    integrator.addComputePerDof("x", "x + 0.5*dt*v");
    integrator.addComputePerDof("x1", "x");
    integrator.addConstrainPositions();
    integrator.addComputePerDof("v", "v + (x-x1)/dt");

    // QMEDSIntegrator.setRandomNumberSeed(42);
    integrator.setIntegrationForceGroups((1 << 0) | (1 << (3 * num_endstates + 1)));

    return integrator;
}

void testSetA(Platform &platform) {
    // Create a system representing a single water molecule.
    std::ifstream system_xml("input/seta_wat_qm.xml"); // TODO make this a path and make sure cmake knows where

    std::unique_ptr<System> system(OpenMM::XmlSerializer::deserialize<System>(system_xml));

    CustomIntegrator integrator = qmedsIntegrator(6);

    std::vector<OpenMM::Vec3> positions = extractPositionsFromPDB("input/seta_wat_qm.pdb");
    Context context(*system, integrator, platform);
    context.setPositions(positions);

    // Simulate it and make sure that the geometry remains reasonable.
    constexpr size_t num_steps = 100;
    for (std::size_t i: std::ranges::iota_view(0u, num_steps)) {
        std::cout << std::format("\r{}/{}", i+1, num_steps) << std::flush;
        integrator.step(1);
    }
    std::cout << "\n";
}

extern "C" void registerXtbSerializationProxies();

int main() {
    try {
        registerXtbSerializationProxies();
        Platform::loadPluginsFromDirectory(PLUGINS_DIR);
        printf("\n");
//        for (int i = 0; i < Platform::getNumPlatforms(); i++) {
//            Platform &platform = Platform::getPlatform(i);
//            printf("Testing %s\n", platform.getName().c_str());
//            testSetA(platform);
//        }
        printf("Testing CPU\n");
        testSetA(Platform::getPlatform("CPU"));
    }
    catch (const OpenMMException &e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
