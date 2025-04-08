%module openmmxtb

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "std_vector.i"

%{
#include "XtbForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%template(vectori) std::vector<int>;
%template(vectord) std::vector<double>;

namespace XtbPlugin {


class XtbForce : public OpenMM::Force {
public:
    enum Method {
        GFN1xTB = 0,
        GFN2xTB = 1,
        GFNFF = 2
    };
    XtbForce(Method method, double charge, int multiplicity, bool periodic, const std::vector<int>& particleIndices, const std::vector<int>& atomicNumbers);
    XtbForce(Method method, double charge, int multiplicity, bool periodic, const std::vector<int>& particleIndices, const std::vector<int>& atomicNumbers, const std::vector<int>& pcIndices, const std::vector<int>& pcNumbers, const std::vector<double>& pcCharges, const std::vector<int>& pcChargeGroups, double pcCutoff);
    Method getMethod() const;
    void setMethod(Method method);
    double getCharge() const;
    void setCharge(double charge);
    int getMultiplicity() const;
    void setMultiplicity(int multiplicity);
    const std::vector<int>& getParticleIndices() const;
    void setParticleIndices(const std::vector<int>& indices);
    const std::vector<int>& getAtomicNumbers() const;
    void setAtomicNumbers(const std::vector<int>& numbers);
    const std::vector<double>& getPointCharges() const;
    void setPointCharges(const std::vector<double>& charges);
    const std::vector<int>& getPointChargeIndices() const;
    void setPointChargeIndices(const std::vector<int>& indices);
    const std::vector<int>& getPointChargeNumbers() const;
    void setPointChargeNumbers(const std::vector<int>& numbers);
    const std::vector<int>& getChargeGroups() const;
    void setChargeGroups(const std::vector<int>& chargeGroups);
    double getPointChargeCutoff() const;
    void setPointChargeCutoff(double cutoff);
    bool hasElectrostaticEmbedding() const;
    bool usesPeriodicBoundaryConditions() const;
    void setUsesPeriodicBoundaryConditions(bool periodic);

    /*
     * Add methods for casting a Force to a XtbForce.
    */
    %extend {
        static XtbPlugin::XtbForce& cast(OpenMM::Force& force) {
            return dynamic_cast<XtbPlugin::XtbForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<XtbPlugin::XtbForce*>(&force) != NULL);
        }
    }
};

}
