/* -------------------------------------------------------------------------- *
 *                        Simbody(tm) Example: SIMBICON                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia, Jack Wang                                           *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "BipedSystem.h"
#include <Simbody.h>

#include <fstream>
#include <cassert>
#include <vector>

using namespace std;
using namespace SimTK;


// #define RIGID_CONTACT

// Normally SIMBICON has 4 states per gait cycle (i.e. 2 per leg), 2 state is
// simplified and not as realistic.
#define TWO_STATE
// TODO USE_GLOBAL_HIPROT
// Drop landing doesn't use controller so you can run with models other than
// the humanoid upon which the controller depends.
//#define DROP_LANDING

namespace { // file-scope symbols

const Vec3 UnitX(1.0, 0.0, 0.0);
const Vec3 UnitY(0.0, 1.0, 0.0);
const Vec3 UnitZ(0.0, 0.0, 1.0);

double clamp( double x, double max) {
	if (x > max)
		x = max;

	if (x < -max)
		x = -max;

	return x;
}
}

// TODO remove
const Real RealTimeFactor
    = 1.0; // try to run in real time

// TODO remove
const int NumActuators = 30;

// SIMBICON

#define EVENT_PERIOD 0.0005
//#define RIGID_CONTACT


// TODO
class SimbiconStateHandler;

class SIMBICON : public SimTK::Force::Custom::Implementation
{
public:

    /// The default arguments define the motion that the controller will try to
    /// execute. Only symmetrical 4-state motions are permitted. The first
    /// element of each Vec2 is for states 0 and 2. The second element is for
    /// states 1 and 3.
    ///
    ///   * deltaT: state dwell duration
    ///   * cd: position balance feedback coefficient
    ///   * cdLat: position balance feedback coefficient, lateral (for 3D gait)
    ///   * cv: velocity balance feedback coefficient
    ///   * cvLat: velocity balance feedback coefficient, lateral (for 3D gait)
    ///   * tor: torso target angle
    ///   * swh: swing-hip target angle
    ///   * swhLat: swing-hip target angle, lateral (for 3D gait)
    ///   * swk: swing-knee target angle
    ///   * swa: swing-ankle target angle
    ///   * stk: stance-knee target angle
    ///   * sta: stance-ankle target angle
    ///
    /// See the SIMBICON paper for information about these parameters.
    ///
	SIMBICON(Biped& biped,
            Real minSIMBICONStateDuration=0.1,
            Vec2 deltaT=Vec2(0.30, NaN),
            Vec2 cd=Vec2(0.5, 0.5),
            Vec2 cdLat=Vec2(0.5, 0.5),
            Vec2 cv=Vec2(0.2, 0.2),
            Vec2 cvLat=Vec2(0.2, 0.2),
            Vec2 tor=Vec2(0.0, 0.0),
            Vec2 swh=Vec2(0.5, -0.1),
            Vec2 swhLat=Vec2(0.0, 0.0),
            Vec2 swk=Vec2(-1.1, -0.05),
            Vec2 swa=Vec2(0.6, 0.15),
            Vec2 stk=Vec2(-0.05, -0.1),
            Vec2 sta=Vec2(0.0, 0.0)
            );

    // TODO
    ~SIMBICON() { m_global_angles_file.close(); }

    /// The possible states of the finite state machine.
    enum SIMBICONState {
        NOT_IN_STATE_MACHINE = -1, // initial state.
        STATE0 = 0, // Left stance.
        STATE1 = 1, // Right foot strike.
        STATE2 = 2, // Right stance.
        STATE3 = 3  // Left foot strike.
    };

    void realizeTopology(State& s) const OVERRIDE_11
    {
        // NOTE: the next two assignments discard the const qualifier of this
        // method. There is a very specific reason for this! Sherm explained
        // this to me as follows:
        //
        // 1. The state has a cache. Even when the state is const, it's okay to
        // write to the cache. This is because the cache is something that is
        // computed from a (const) state. So long as the state is const, we are
        // able to regenerate the cache appropriately.
        //
        // 2. We can think of the model-building process as having a state and
        // a cache. The state is the things that make up the model, such as the
        // MobilizedBody's. The cache in this case is quantities computed from
        // this model-buidling state, such as the number of mobilized bodies.
        //
        // 3. Therefore, quantities computed about the model-building state do
        // not need to be const.
        //
        // 4. The indicies in the state of my two AutoUpdateDiscreteVariable's
        // are part of the model-building cache, and so I can remove the const
        // qualifier.
        const_cast<SIMBICON*>(this)->m_simbiconStateIndex =
            m_forces.allocateAutoUpdateDiscreteVariable(s,
                    Stage::Acceleration,
                    new Value<SIMBICONState>(NOT_IN_STATE_MACHINE),
                    Stage::Dynamics);

        const_cast<SIMBICON*>(this)->m_simbiconStateCacheIndex =
            m_forces.getDiscreteVarUpdateIndex(s, m_simbiconStateIndex);
    }

    const SIMBICONState getSIMBICONState(const State& s) const
    {
        const AbstractValue& m = m_forces.getDiscreteVariable(
                s, m_simbiconStateIndex);
        return Value<SIMBICONState>::downcast(m).get();
    }
    Real getSIMBICONStateStartTime(const State& s) const
    {
        return m_forces.getDiscreteVarLastUpdateTime(s,
                m_simbiconStateIndex);
    }

    /// The control consists of a few components. This method calls other
    /// methods that take care of these individual components. TODO
    void calcForce(const SimTK::State&                state,
                   SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                   SimTK::Vector_<SimTK::Vec3>&       particleForces,
                   SimTK::Vector&                     mobilityForces) const
                   OVERRIDE_11;

    Real calcPotentialEnergy(const SimTK::State& state) const OVERRIDE_11
    { return 0; }

    /// Used to identify gains for different parts of the body.
    enum GainGroup {
        generic,
        neck,
        back,
        hip_flexion_adduction,
        hip_rotation,
        knee,
        arm_flexion_adduction,
        arm_rotation,
        ankle_flexion,
        ankle_inversion,
        toe
    };

    // TODO
	void updateGlobalAngles(const SimTK::State& s);

private:

    /// Determine if the SIMBICON state of the state machine has changed.
    void updateSIMBICONState(const State& s) const;

    void setSIMBICONState(const State& s, SIMBICONState simbiconState, bool LC, bool RC) const
    {
        m_forces.updDiscreteVarUpdateValue(s, m_simbiconStateIndex) =
            Value<SIMBICONState>(simbiconState);
        m_forces.markDiscreteVarUpdateValueRealized(s, m_simbiconStateIndex);

        std::cout << "t: " << s.getTime() << " SIMBICONState: " <<
            getSIMBICONState(s) << " LC: " << LC << " RC: " << RC << std::endl;

        // TODO
        /*
		if (simbiconState == STATE0 || simbiconState == STATE2) {
			// reset the swing thigh orientation when the swing leg changes
			for (int i = 0; i < 2; i++) {
				const_cast<SIMBICON*>(this)->_lastSWTAngle[i] = -100.0;
				const_cast<SIMBICON*>(this)->_curSWTAngle[i] = -100.0;
			}
        }
        */
    }

    /// Apply generalized force to a specific coord using a
    /// proportional-derivative control law that tracks thetad.
    /// The force is added into the appropriate place in `mobForces.`
    /// For convenience / inspection, the force is returned as well.
    // TODO
    Real coordPDControl(const State& s,
            const Biped::Coordinate coord, const GainGroup gainGroup,
            const Real thetad, Vector& mobilityForces) const
    {
        // Prepare quantities.
        Real kp = m_proportionalGains.at(gainGroup);
        Real kd = m_derivativeGains.at(gainGroup);
        Real q = m_biped.getQ(s, coord);
        Real u = m_biped.getU(s, coord);

        // PD control law:
        Real force = clamp(kp * (thetad - q) - kd * u, kp);

        // Update the proper entry in mobForces.
        m_biped.addInForce(coord, force, mobilityForces);

        return force;
    }
    // TODO
    Real coordPDControlControls(const State& s,
            const Biped::Coordinate coord, const GainGroup gainGroup,
            const Real thetad) const
    {
        // Prepare quantities.
        Real kp = m_proportionalGains.at(gainGroup);
        Real kd = m_derivativeGains.at(gainGroup);
        Real q = m_biped.getQ(s, coord);
        Real u = m_biped.getU(s, coord);

        // PD control law:
        Real force = clamp(kp * (thetad - q) - kd * u, kp);

        return force;
    }

    // TODO
    void calcGainsFromStrength(GainGroup group, double& kp, double& kd) const {
        kp = m_proportionalGains.at(group);
        kd = m_derivativeGains.at(group);
    }

    Biped& m_biped;
    const GeneralForceSubsystem& m_forces;

    const Real m_minSIMBICONStateDuration;

    /// SIMBICON parameters.
    const Vec2 m_deltaT;
    const Vec2 m_cd;
    const Vec2 m_cdLat;
    const Vec2 m_cv;
    const Vec2 m_cvLat;
    const Vec2 m_tor;
    const Vec2 m_swh;
    const Vec2 m_swhLat;
    const Vec2 m_swk;
    const Vec2 m_swa;
    const Vec2 m_stk;
    const Vec2 m_sta;

    /// Proportional (position) gains (kp).
    std::map<GainGroup, Real> m_proportionalGains;
    /// Derivative (speed) gains; mostly chosen for critical damping.
    std::map<GainGroup, Real> m_derivativeGains;

    // The SIMBICON controller itself has some state.
    // ----------------------------------------------
    // Written during the realizeTopology step.

    /// The SIMBICON state of the controller for the finite state machine
    /// (e.g., 0, 1, 2, or 3). We read the value of this variable using this
    //index.
    DiscreteVariableIndex m_simbiconStateIndex;

    /// We set the value of the SIMBICON state using this index.
    CacheEntryIndex m_simbiconStateCacheIndex;

    void computeControls(const State& s, Vector& controls, Vector& mobForces) const;

	void fillInHipJointControls( const SimTK::State& s,
		SimTK::Vector& controls ) const;

    /// Angles used for balance control. These angles describe the orientation
    /// of parts of the Biped relative to the sagittal and coronal planes.
    struct GlobalAngles
    {
        // TODO
        Real swingHipFlexion;
        Real swingHipAbduction;
        Real trunkExtension;
        Real trunkBending;
        void reset() {
            swingHipFlexion = -100;
            swingHipAbduction = -100;
            trunkExtension = -100;
            trunkBending = -100;
        }
    };

    /// Global angles calculated in the most recent call to updateGlobalAngles.
    GlobalAngles m_currentGlobalAngles;
    /// Global angles from the prevoius call to updateGlobalAngles.
    GlobalAngles m_previousGlobalAngles;

    // TODO
    friend SimbiconStateHandler;

    // TODO
    ofstream m_global_angles_file;

};

class SimbiconStateHandler : public SimTK::PeriodicEventHandler {
public:
    SimbiconStateHandler(Biped& m, SIMBICON& simctrl, SimTK::Real interval);

    void handleEvent(SimTK::State& s, SimTK::Real accuracy, bool& shouldTerminate) const;
private:
    Biped& _model;
    SIMBICON& _simctrl;
};
Real criticallyDampedDerivativeGain(Real proportionalGain) {
    return 2.0 * std::sqrt(proportionalGain);
}

SIMBICON::SIMBICON(Biped& biped,
        Real minSIMBICONStateDuration,
        Vec2 deltaT, Vec2 cd, Vec2 cdLat, Vec2 cv,
        Vec2 cvLat, Vec2 tor, Vec2 swh, Vec2 swhLat, Vec2 swk, Vec2 swa, Vec2
        stk, Vec2 sta)
    : m_biped(biped), m_forces(m_biped.getForceSubsystem()),
      m_minSIMBICONStateDuration(minSIMBICONStateDuration),
      m_deltaT(deltaT), m_cd(cd),
      m_cdLat(cdLat), m_cv(cv), m_cvLat(cvLat), m_tor(tor), m_swh(swh),
      m_swhLat(swhLat), m_swk(swk), m_swa(swa), m_stk(stk), m_sta(sta)

{
    // TODO
    m_global_angles_file.open("global_angles.txt");
    m_global_angles_file << "time swtangle0 swtangle1 trunkangle0 trunkangle1 swingHipFlexion swingHipAbduction TrunkExtension TrunkBending" << endl;
    m_proportionalGains[generic] = 300;
    m_proportionalGains[neck] = 100;
    m_proportionalGains[back] = 300;
    m_proportionalGains[hip_flexion_adduction] = 1000;
    m_proportionalGains[hip_rotation] = 300;
    m_proportionalGains[knee] = 300;
    m_proportionalGains[arm_flexion_adduction] = 300;
    m_proportionalGains[arm_rotation] = 300;
    m_proportionalGains[ankle_flexion] = 300;
    m_proportionalGains[ankle_inversion] = 30;
    m_proportionalGains[toe] = 30;

    m_derivativeGains[generic] = criticallyDampedDerivativeGain(
            m_proportionalGains[generic]);
    m_derivativeGains[neck] = criticallyDampedDerivativeGain(
            m_proportionalGains[neck]);
    m_derivativeGains[back] = criticallyDampedDerivativeGain(
            m_proportionalGains[back]);
    // Overdamped.
    m_derivativeGains[hip_flexion_adduction] = 100;
    m_derivativeGains[hip_rotation] = criticallyDampedDerivativeGain(
            m_proportionalGains[hip_rotation]);
    m_derivativeGains[knee] = criticallyDampedDerivativeGain(
            m_proportionalGains[knee]);
    m_derivativeGains[arm_flexion_adduction] = criticallyDampedDerivativeGain(
            m_proportionalGains[arm_flexion_adduction]);
    m_derivativeGains[arm_rotation] = criticallyDampedDerivativeGain(
            m_proportionalGains[arm_rotation]);
    m_derivativeGains[ankle_flexion] = criticallyDampedDerivativeGain(
            m_proportionalGains[ankle_flexion]);
    m_derivativeGains[ankle_inversion] = criticallyDampedDerivativeGain(
            m_proportionalGains[ankle_inversion]);
    m_derivativeGains[toe] = criticallyDampedDerivativeGain(
            m_proportionalGains[toe]);

        // TODO
    m_currentGlobalAngles.reset();
    m_previousGlobalAngles.reset();

}

void SIMBICON::calcForce(const State&         s,
                         Vector_<SpatialVec>& /*bodyForces*/,
                         Vector_<Vec3>&       /*particleForces*/,
                         Vector&              mobilityForces) const
{
    updateSIMBICONState(s);
    Vector controls(NumActuators);
    controls.setToZero();
    computeControls(s,controls, mobilityForces);

    for (int i=0; i<NumActuators; ++i)
    {
        if (controls[i] != 0)
        {
            mobilityForces[m_biped.getUIndex(Biped::Coordinate(i))] = controls[i];
        }
    }
}


void SIMBICON::fillInHipJointControls( const State& s, Vector& controls ) const {

    // Which leg is in stance?
    // =======================
	const SIMBICONState simbiconState = getSIMBICONState(s);

    Biped::Coordinate swing_hip_flexion;
    Biped::Coordinate swing_hip_adduction;
    Biped::Coordinate stance_hip_flexion;
    Biped::Coordinate stance_hip_adduction;
    const MobilizedBody* stanceFoot;
    if (simbiconState == STATE0 || simbiconState == STATE1)
    {
        // Left leg is in stance.
        swing_hip_flexion = Biped::hip_r_flexion;
        swing_hip_adduction = Biped::hip_r_adduction;
        stance_hip_flexion = Biped::hip_l_flexion;
        stance_hip_adduction = Biped::hip_l_adduction;
        stanceFoot = &m_biped.getBody(Biped::foot_l);
    }
    else if (simbiconState == STATE2 || simbiconState == STATE3)
    {
        // Right leg is in stance.
        swing_hip_flexion = Biped::hip_l_flexion;
        swing_hip_adduction = Biped::hip_l_adduction;
        stance_hip_flexion = Biped::hip_r_flexion;
        stance_hip_adduction = Biped::hip_r_adduction;
        stanceFoot = &m_biped.getBody(Biped::foot_r);
    }
    else SimTK_THROW(Exception::Base);

    // Get the correct SIMBICON parameters for the current SIMBICONState.
    // ==================================================================
    // simbiconState stateIdx
    // ------------- --------
    // 0             0
    // 1             1
    // 2             0
    // 3             1
    // TODO move stateIdx into a subfunction.
    #ifdef TWO_STATE
        const int stateIdx = 0;
    #else
        const int stateIdx = simbiconState % 2;
    #endif
    // Target angle for swing hip flexion.
	double swh = m_swh[stateIdx];
    double swhLat = m_swhLat[stateIdx];
    // Feedback gains for the global angles global angular rates.
	double cd = 0.2; // global tipping feedback m_cd[stateIdx] TODO
    double cdLat = 0.2; // TODO  m_cdLat[stateIdx];
	double cv = m_cv[stateIdx];
    double cvLat = m_cvLat[stateIdx];

    // Compute necessary translation-related quantities.
    // =================================================
    const SimbodyMatterSubsystem& matter = m_biped.getMatterSubsystem();

    // Whole-body mass center, expressed in ground.
	Vec3 massCenterLoc = matter.calcSystemMassCenterLocationInGround(s);
    // This is 'v' in Yin, 2007.
	Vec3 massCenterVel = matter.calcSystemMassCenterVelocityInGround(s);

    // Ankle, expressed in ground.
    const Transform& X = stanceFoot->getInboardFrame(s);
    Vec3 stanceAnkleLocInParent = X.p();
    Vec3 stanceAnkleLoc = stanceFoot->getParentMobilizedBody()
        .findStationLocationInGround(s, stanceAnkleLocInParent);

	Vec3 massCenterLocFromStanceAnkle = massCenterLoc - stanceAnkleLoc;

    // Normals to the sagittal and coronal planes.
    UnitVec3 sagittalNormal = m_biped.getNormalToSagittalPlane(s);
    UnitVec3 coronalNormal = m_biped.getNormalToCoronalPlane(s);

	double d_sag = dot(coronalNormal, massCenterLocFromStanceAnkle);
	double v_sag = dot(coronalNormal, massCenterVel);
	double d_cor = dot(sagittalNormal, massCenterLocFromStanceAnkle);
	double v_cor = dot(sagittalNormal, massCenterVel);

    // Compute necessary translation-related quantities.
    // =================================================
    // The angles are updated in updateGlobalAngles(); here we compute rates.
    // Check there's a valid value for the previous angles.  Otherwise just use
    // 0 for the angular rates.
    double swingHipFlexionRate = 0;
    double swingHipAbductionRate = 0;
    double trunkExtensionRate = 0;
    double trunkBendingRate = 0;

    const GlobalAngles& cur = m_currentGlobalAngles;
    const GlobalAngles& prev = m_previousGlobalAngles;

    if (prev.swingHipFlexion > -100) {
        swingHipFlexionRate =
            (cur.swingHipFlexion - prev.swingHipFlexion) / EVENT_PERIOD;
    }
    if (prev.swingHipAbduction > -100) {
        swingHipAbductionRate =
            (cur.swingHipAbduction - prev.swingHipAbduction) / EVENT_PERIOD;
    }
	if (prev.trunkExtension > -100) {
        trunkExtensionRate =
            (cur.trunkExtension - prev.trunkExtension) / EVENT_PERIOD;
	}
	if (prev.trunkBending > -100) {
        trunkBendingRate =
            (cur.trunkBending - prev.trunkBending) / EVENT_PERIOD;
	}

	// sign change is needed for one of the stance legs in the coronal plane
	double sign = 1;
    if (simbiconState == STATE0 || simbiconState == STATE1) // left stance
		sign = -1;

    // TODO
    double kp, kd; // position, derivative gains for hip flex/adduction
    calcGainsFromStrength(hip_flexion_adduction, kp, kd);

    double desiredSwingHipFlexionAngle = swh + (cd*d_sag + cv*v_sag);
    double swingHipFlexionTorque =
        + kp * (desiredSwingHipFlexionAngle - cur.swingHipFlexion)
        - kd * swingHipFlexionRate;
    controls[swing_hip_flexion] = clamp(swingHipFlexionTorque, kp);
    controls[swing_hip_adduction] = sign*(clamp(  kp*(swhLat + (cdLat*d_cor +
                        cvLat*v_cor) - cur.swingHipAbduction) -
                kd*swingHipAbductionRate, kp));
    // TODO abd/add

	// use stance hip to control the trunk
	controls[stance_hip_flexion] =  -kp*(0. - cur.trunkExtension) + kd*trunkExtensionRate;
	controls[stance_hip_flexion] -= controls[swing_hip_flexion];
	controls[stance_hip_flexion] = clamp(controls[stance_hip_flexion], kp);

	controls[stance_hip_adduction] = sign*(kp*(0. - cur.trunkBending) - kd*trunkBendingRate);
	controls[stance_hip_adduction] -= controls[swing_hip_adduction];
	controls[stance_hip_adduction] = clamp(controls[stance_hip_adduction], kp);

}

void SIMBICON::computeControls(const State& s, Vector& controls, Vector& mobForces) const
{
#ifdef DROP_LANDING
	for (int i = 0; i < controls.size(); i++)
			controls[i] = 0.0;
	return;
#else

    const SIMBICONState simbiconState = getSIMBICONState(s);

    // For most joints, track target theta of 0.0 degrees.
    // ===================================================
    coordPDControl(s, Biped::neck_extension, neck, 0.0, mobForces);
    coordPDControl(s, Biped::neck_bending, neck, 0.0, mobForces);
    coordPDControl(s, Biped::neck_rotation, neck, 0.0, mobForces);

    coordPDControl(s, Biped::back_tilt, back, 0.0, mobForces);
    coordPDControl(s, Biped::back_list, back, 0.0, mobForces);
    coordPDControl(s, Biped::back_rotation, back, 0.0, mobForces);

    coordPDControl(s, Biped::shoulder_r_flexion, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::shoulder_l_flexion, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::shoulder_r_adduction, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::shoulder_l_adduction, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::elbow_r_flexion, arm_flexion_adduction, 0.0, mobForces);
    coordPDControl(s, Biped::elbow_l_flexion, arm_flexion_adduction, 0.0, mobForces);

    coordPDControl(s, Biped::shoulder_r_rotation, arm_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::shoulder_l_rotation, arm_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::elbow_r_rotation, arm_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::elbow_l_rotation, arm_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::hip_r_rotation, hip_rotation, 0.0, mobForces);
    coordPDControl(s, Biped::hip_l_rotation, hip_rotation, 0.0, mobForces);

    // TODO change to hip_flexion_adduction gaingroup.
    coordPDControl(s, Biped::hip_r_flexion, generic, 0.0, mobForces);
    coordPDControl(s, Biped::hip_l_flexion, generic, 0.0, mobForces);
    coordPDControl(s, Biped::hip_r_adduction, generic, 0.0, mobForces);
    coordPDControl(s, Biped::hip_l_adduction, generic, 0.0, mobForces);

    coordPDControl(s, Biped::ankle_r_inversion, ankle_inversion, 0.0, mobForces);
    coordPDControl(s, Biped::ankle_l_inversion, ankle_inversion, 0.0, mobForces);

    coordPDControl(s, Biped::mtp_r_dorsiflexion, toe, 0.0, mobForces);
    coordPDControl(s, Biped::mtp_l_dorsiflexion, toe, 0.0, mobForces);

    // Deal with limbs whose target angle is affected by the state machine.
    // ====================================================================

    // Which leg is in stance?
    // -----------------------
    if (simbiconState != NOT_IN_STATE_MACHINE)
    {
        Biped::Coordinate swing_knee_extension;
        Biped::Coordinate swing_ankle_dorsiflexion;
        Biped::Coordinate stance_knee_extension;
        Biped::Coordinate stance_ankle_dorsiflexion;
        if (simbiconState == STATE0 || simbiconState == STATE1)
        {
            // Left leg is in stance.
            swing_knee_extension = Biped::knee_r_extension;
            swing_ankle_dorsiflexion = Biped::ankle_r_dorsiflexion;
            stance_knee_extension = Biped::knee_l_extension;
            stance_ankle_dorsiflexion = Biped::ankle_l_dorsiflexion;
        }
        else if (simbiconState == STATE2 || simbiconState == STATE3)
        {
            // Right leg is in stance.
            swing_knee_extension = Biped::knee_l_extension;
            swing_ankle_dorsiflexion = Biped::ankle_l_dorsiflexion;
            stance_knee_extension = Biped::knee_r_extension;
            stance_ankle_dorsiflexion = Biped::ankle_r_dorsiflexion;
        }


        // Apply swing/stance-dependent PD control to lower limb sagittal coords
        // ---------------------------------------------------------------------
        // simbiconState stateIdx
        // ------------- --------
        // 0             0
        // 1             1
        // 2             0
        // 3             1
        #ifdef TWO_STATE
            const int stateIdx = 0;
        #else
            const int stateIdx = simbiconState % 2;
        #endif

        coordPDControl(s, swing_knee_extension, knee, m_swk[stateIdx],
                mobForces);
        coordPDControl(s, stance_knee_extension, knee, m_stk[stateIdx],
                mobForces);

        coordPDControl(s, swing_ankle_dorsiflexion, ankle_flexion,
                m_swa[stateIdx], mobForces);
        coordPDControl(s, stance_ankle_dorsiflexion, ankle_flexion,
                m_sta[stateIdx], mobForces);
    }
    else
    {
        // TODO can remove later since i'm just doing this to identically match results.
        coordPDControl(s, Biped::knee_l_extension, knee, 0.0, mobForces);
        coordPDControl(s, Biped::knee_r_extension, knee, 0.0, mobForces);
        coordPDControl(s, Biped::ankle_l_dorsiflexion, ankle_flexion, 0.0, mobForces);
        coordPDControl(s, Biped::ankle_r_dorsiflexion, ankle_flexion, 0.0, mobForces);
    }

	if (simbiconState != NOT_IN_STATE_MACHINE) {
		fillInHipJointControls(s, controls);
    }

    // TODO can I remove this?
	Vec3 com = m_biped.getMatterSubsystem().calcSystemMassCenterLocationInGround(s);
	if (com[1] < 0.7) {
		for (int i = 0; i < controls.size(); i++)
			controls[i] = 0.0;
	}
	return;
#endif
}

Vec3 projectionOntoPlane(const Vec3& vecToProject,
                         const UnitVec3& normalToPlane)
{ return vecToProject - dot(vecToProject, normalToPlane) * normalToPlane; }

void SIMBICON::updateGlobalAngles(const State& s)
{
	const SIMBICONState simbiconState = getSIMBICONState(s);

    // Only compute angles if we have entered the state machine.
    if (simbiconState == NOT_IN_STATE_MACHINE) return;

    GlobalAngles& cur = m_currentGlobalAngles;

    // Save the previous values; we use these to approximate derivatives.
    m_previousGlobalAngles = cur;

    // Which leg is in stance, etc.?
    // =============================
    const MobilizedBody* stanceFoot;
    const MobilizedBody* swingThigh;
    if (simbiconState == STATE0 || simbiconState == STATE1)
    {
        // Left leg is in stance.
        stanceFoot = &m_biped.getBody(Biped::foot_l);
        swingThigh = &m_biped.getBody(Biped::thigh_r);
    }
    else if (simbiconState == STATE2 || simbiconState == STATE3)
    {
        // Right leg is in stance.
        stanceFoot = &m_biped.getBody(Biped::foot_r);
        swingThigh = &m_biped.getBody(Biped::thigh_l);
    }
    else SimTK_THROW(Exception::Base);

    // Preliminary quantities.
    // =======================
    const SimbodyMatterSubsystem& matter = m_biped.getMatterSubsystem();

    // Normals to the sagittal and coronal planes.
    UnitVec3 sagittalNormal = m_biped.getNormalToSagittalPlane(s);
    UnitVec3 coronalNormal = m_biped.getNormalToCoronalPlane(s);

    // UnitVec3(YAxis) gives a vector perpendicular to ground, directed up.
    UnitVec3 up(YAxis);

    // This vector is in sagittal plane, and is parallel to the ground.
    UnitVec3 sagittalHorizontal(cross(up, sagittalNormal));

    // A horizontal vector in the coronal plane, directed to the 'right'.
    UnitVec3 coronalHorizontal(-cross(up, coronalNormal));

    // Global hip angles.
    // ==================
    // This vector is along the axis of the thigh, directed proximally.
    UnitVec3 swingThighAxialDir(swingThigh->getBodyRotation(s).y());

    // Hip flexion.
    // ------------
    // Take the vector directed along the hip, project it into the sagittal
    // plane, and get the angle between that and the horizontal.

    // Projection of the swingThighAxialDir into the sagittal plane.
    UnitVec3 sagittalSwingThighAxialDir(
            projectionOntoPlane(swingThighAxialDir, sagittalNormal));

    // Angle between the thigh axis and the vertical, but done so in a way
    // that the angle is negative if the thigh axis is directed anteriorly.
    cur.swingHipFlexion =
        acos(dot(sagittalSwingThighAxialDir, sagittalHorizontal)) - Pi / 2;

    // Hip abduction.
    // --------------
    // Projection of the swingThighAxialDir into the coronal plane.
    UnitVec3 coronalSwingThighAxialDir(
            projectionOntoPlane(swingThighAxialDir, coronalNormal));

    // Angle between the thigh axis and the vertical, in the coronal plane.
    cur.swingHipAbduction =
        acos(dot(coronalSwingThighAxialDir, coronalHorizontal)) - Pi / 2;
    // TODO change sign for when left leg is in swing.

    // Global trunk angles (well, actually, pelvis angles).
    // ====================================================
	const MobilizedBody& pelvis = m_biped.getBody(Biped::pelvis);

    // Directed proximally.
    UnitVec3 pelvisAxialDir(pelvis.getBodyRotation(s).y());

    // Trunk extension (sagittal plane).
    // ---------------------------------
    // Take the vector directed up the pelvis, project it into the sagittal
    // plane, and get the angle between that and the horizontal.

    // Projection of the pelvisAxialDir into the sagittal plane.
    UnitVec3 sagittalPelvisAxialDir(
            projectionOntoPlane(pelvisAxialDir, sagittalNormal));

    // Angle between pelvis axis and the vertical, in the sagittal plane.
    // TODO move comments elsewhere
    cur.trunkExtension =
        acos(dot(sagittalPelvisAxialDir, sagittalHorizontal)) - Pi/2;

    // Trunk bending (coronal plane).
    // ------------------------------
    // Projection of the pelvisAxialDir into the coronal plane.
    UnitVec3 coronalPelvisAxialDir(
            projectionOntoPlane(pelvisAxialDir, coronalNormal));

    // Angle between the pelvis axis and the vertical, in the coronal plane.
    cur.trunkBending =
        acos(dot(coronalPelvisAxialDir, coronalHorizontal)) - Pi/2;
}

/* TODO
        // Translation-related quantities.
        // -------------------------------

        // Whole-body mass center, expressed in ground.
        const Vec3 massCenterLoc = matter.calcSystemMassCenterLocationInGround(s);
        // This is 'v' in Yin, 2007.
        const Vec3 massCenterVel = matter.calcSystemMassCenterVelocityInGround(s);

        // Ankle, expressed in ground.
        const Transform& X = stanceFoot->getInboardFrame(s);
        Vec3 stanceAnkleLocInParent = X.p();
        Vec3 stanceAnkleLoc = stanceFoot->getParentMobilizedBody()
            .findStationLocationInGround(s, stanceAnkleLocInParent);

        // Displacement of whole-body COM from ankle, expressed in ground.
        // This is called 'd' in Yin, 2007.
        Vec3 massCenterLocFromStanceAnkle = massCenterLoc - stanceAnkleLoc;

        // Normals to the sagittal and coronal planes.
        UnitVec3 sagittalNormal = m_biped.getNormalToSagittalPlane(s);
        UnitVec3 coronalNormal = m_biped.getNormalToCoronalPlane(s);

        // Projection of 'd' and 'v' into the forward / sagittal and lateral /
        // coronal planes.  We apply balance control separately for these two
        // planes.
        Real massCenterLocFromStanceAnkleForwardMeasure =
            dot(coronalNormal, massCenterLocFromStanceAnkle);
        Real massCenterVelFromStanceAnkleForwardMeasure =
            dot(coronalNormal, massCenterVel);
        Real massCenterLocFromStanceAnkleLateralMeasure =
            dot(sagittalNormal, massCenterLocFromStanceAnkle);
        Real massCenterVelFromStanceAnkleLateralMeasure =
            dot(sagittalNormal, massCenterVel);
            */

void SIMBICON::updateSIMBICONState(const State& s) const
{
    Real lForce, rForce;
    m_biped.findContactForces(s, lForce, rForce);
    const bool lContact = (lForce > 0);
    const bool rContact = (rForce > 0);


#ifndef DROP_LANDING
	const SIMBICONState simbiconState = getSIMBICONState(s);
	const Real duration = s.getTime() - getSIMBICONStateStartTime(s);

    // simbiconState stateIdx
    // ------------- --------
    // 0             0
    // 1             1
    // 2             0
    // 3             1
    const int stateIdx = simbiconState % 2;

    switch (simbiconState)
    {
        // Neither foot has contacted the ground yet.
        case NOT_IN_STATE_MACHINE:

            // Entering right stance.
            if (rContact) {
                setSIMBICONState(s, STATE2, lContact, rContact);
            }
            // Entering left stance.
            else if (lContact) {
                setSIMBICONState(s, STATE0, lContact, rContact);
            }
            break;

        // Left stance.
        case STATE0:

            // Stay in this state for \delta t seconds.
            if (duration > m_deltaT[stateIdx]) {
                setSIMBICONState(s, STATE1, lContact, rContact);
            }
            // Already entered right stance; skip STATE1.
            else if (rContact && duration > m_minSIMBICONStateDuration) {
                setSIMBICONState(s, STATE2, lContact, rContact);
            }
            break;

        // Right foot strike.
        case STATE1:

            // Stay in this state until the right foot makes contact.
            if (rContact && duration > m_minSIMBICONStateDuration) {
                setSIMBICONState(s, STATE2, lContact, rContact);
            }
            break;

        // Right stance.
        case STATE2:

            // Stay in this state for \delta t seconds.
            if (duration > m_deltaT[stateIdx]) {
                setSIMBICONState(s, STATE3, lContact, rContact);
            }
            // Already entered right stance; skip STATE3.
            else if (lContact && duration > m_minSIMBICONStateDuration) {
                setSIMBICONState(s, STATE0, lContact, rContact);
            }
            break;

        // Left foot strike.
        case STATE3:

            // Stay in this state until the left foot makes contact.
            if (lContact && duration > m_minSIMBICONStateDuration) {
                setSIMBICONState(s, STATE0, lContact, rContact);
            }
            break;
    }

#endif
}

SimbiconStateHandler::SimbiconStateHandler(Biped& model, SIMBICON& simctrl,
    Real interval)
        : PeriodicEventHandler(interval), _model(model), _simctrl(simctrl) {
	}

void SimbiconStateHandler::handleEvent(State& s, Real accuracy, bool& shouldTerminate) const
{
    shouldTerminate = false;

#ifndef DROP_LANDING
    _simctrl.updateGlobalAngles(s);
#endif
}




namespace {
// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo,
                     Real                   interval);
    void handleEvent(State& state, Real accuracy,
                     bool& shouldTerminate) const OVERRIDE_11;
private:
    Visualizer::InputSilo& m_silo;
};

// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input.
class OutputReporter : public PeriodicEventReporter {
public:
    OutputReporter(const Biped& biped,
                   Real            interval);
    void handleEvent(const State& state) const OVERRIDE_11;
private:
    const Biped& m_biped;
};

// Write interesting integrator info to stdout at end of simulation.
void dumpIntegratorStats(double startCPU, double startTime,
                         const Integrator& integ);
}

class ShowContact : public DecorationGenerator {
public:
    ShowContact(const Biped& unis, const SIMBICON* simctrl)
    :   m_unis(unis), m_simctrl(simctrl) {}

    void generateDecorations(const State&                state,
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        DecorativeText energy;
        energy.setIsScreenText(true);
        bool left, right;
        // TODO m_unis.findContactStatus(state, left, right);
        std::string str = "";
        /* TODO if (left) {
            str += "left on,";
        }
        if (right) {
            str += "right on";
        }
        */
        std::ostringstream oss;
        oss << "SIMBICONState: " <<
            m_simctrl->getSIMBICONState(state) << std::endl;
        energy.setText(oss.str());
        geometry.push_back(energy);
        /*
        for (int i=0; i < m_unis.getNumContactElements(); ++i) {
            const MyContactElement& contact = m_unis.getContactElement(i);
            const Vec3 loc = contact.whereToDisplay(state);
            if (!contact.isDisabled(state)) {
                geometry.push_back(DecorativeSphere(.025)
                    .setTransform(loc)
                    .setColor(Red).setOpacity(.25));
                String text = "LOCKED";
                if (contact.hasFrictionElement()) {
                    const MyFrictionElement& friction = contact.getFrictionElement();
                    text = friction.isSticking(state) ? "STICKING"
                                                      : "CONTACT";
                    m_unis.getMultibodySystem().realize(state, Stage::Acceleration);
                    friction.showFrictionForce(state, geometry);
                }
                geometry.push_back(DecorativeText(String(i)+"-"+text)
                    .setColor(White).setScale(.025)
                    .setTransform(loc+Vec3(0,.04,0)));
            } else {
                geometry.push_back(DecorativeText(String(i))
                    .setColor(White).setScale(.025)
                    .setTransform(loc+Vec3(0,.02,0)));
            }
        }
        */
    }
private:
    const Biped& m_unis;
    const SIMBICON* m_simctrl;
};

//==============================================================================
// MAIN FUNCTION
//==============================================================================
int main(int argc, char **argv)
{
// TODO
Vector orig;
std::string origAsString = "~[0.996609 0.0295878 0.064981 -0.0410475 1.58194 1.30901 -0.126126 -0.0226194 0.00668178 -0.00197447 0.0619336 -0.0497828 -0.0279272 0.00388359 0.00695986 0.00386386 0.000842296 0.00031449 -0.0130121 0.00666864 -0.00444311 -5.93963e-06 -0.000837548 0.205327 0.0229502 -0.0517189 -0.3854 0.170585 -0.035032 0.0773673 0.572483 0.0199089 0.0177721 -0.920657 0.524068 -0.0238833 0.235308 -0.563026 -1.93606 0.352619 0.9611 0.104494 -0.104818 0.0277442 -0.0119722 -0.00519416 -0.203116 0.272666 0.384122 -0.0707048 -0.0760685 -0.0571402 -0.024149 -0.0134003 0.116814 0.0146574 0.0571261 0.0004472 0.00915809 -2.17252 0.814301 1.33145 1.88792 -0.00255727 -0.0102517 0.508182 2.21343 -0.432101 0.0321825 -0.65675 -0.000677057 -0.0137553 -0.395188]";
istringstream iss(origAsString);
iss >> orig;
    try {

	double finalTime = 1000;

    Biped biped;
    SimTK::Visualizer                   viz(biped);
    SimTK::Visualizer::InputSilo* userInput = new Visualizer::InputSilo();
    viz.addInputListener(userInput);
    biped.addEventHandler
       (new UserInputHandler(*userInput, Real(0.1))); //check input every 100ms

    biped.addEventReporter(new OutputReporter(biped, .01));
    biped.addEventReporter(new Visualizer::Reporter(viz, RealTimeFactor/30));
    DecorativeText help("Any input to start; ESC to quit");
    help.setIsScreenText(true);
    viz.addDecoration(MobilizedBodyIndex(0),Vec3(0),help);
    biped.updMatterSubsystem().setShowDefaultGeometry(false);

    // Add the controller.
    SIMBICON* simctrl = new SIMBICON(biped);
    Force::Custom simbicon(biped.updForceSubsystem(), simctrl); // takes ownership

#ifndef RIGID_CONTACT
    biped.addEventHandler(new SimbiconStateHandler(biped,*simctrl,
                                                          EVENT_PERIOD));
#endif

    State s;
    biped.initialize(s);

    /* TODO
    printf("Act: u\n");
    for (int i=0; i < NumActuators; ++i) {
        printf("%2d: %d\n", i, int(biped.getUIndex(Biped::Coordinate(i))));
    }
    */

    //biped.toes_r.lockAt(s, .2); biped.toes_l.lockAt(s, .2);
    //biped.foot_r.lockAt(s, Vec2(.2,0)); biped.foot_l.lockAt(s, Vec2(.2,0));
    biped.realize(s, Stage::Instance);

    /* TODO
    printf("SIMBICON 3D:\n");
    printf("%d bodies, %d mobilities, -%d constraint equations -%d motions\n",
        biped.getMatterSubsystem().getNumBodies(), s.getNU(), s.getNMultipliers(),
        biped.getMatterSubsystem().getKnownUDotIndex(s).size());
        */

    biped.getBody(Biped::trunk).setQToFitTranslation(s, Vec3(0,1.5,0));
    biped.getBody(Biped::trunk).setUToFitLinearVelocity(s, Vec3(1,0,0));
    viz.report(s);

// TODO #ifdef RIGID_CONTACT
    viz.addDecorationGenerator(new ShowContact(biped, simctrl));
// TODO #endif

#ifndef RIGID_CONTACT
    // Simulate.
    //CPodesIntegrator integ(biped.system); integ.setOrderLimit(2); integ.setAccuracy(.01);

    //RungeKuttaMersonIntegrator integ(biped.system); integ.setAccuracy(1e-3);
    //RungeKutta2Integrator integ(biped.system); integ.setAccuracy(.1);
    SemiExplicitEuler2Integrator integ(biped); integ.setAccuracy(0.1);
    //SemiImplicitEulerIntegrator integ(biped.system, .002);
    //integ.setConstraintTolerance(.001);
    //integ.setMaximumStepSize(.005);
    TimeStepper ts(biped, integ);
#else
    SemiExplicitEulerTimeStepper ts(biped);
    ts.setImpulseSolverType(SemiExplicitEulerTimeStepper::PLUS);
    ts.setInducedImpactModel(SemiExplicitEulerTimeStepper::Simultaneous);
    //ts.setConstraintTol(1e-5);
    ts.setAccuracy(1e-4);
    ts.setMaxInducedImpactsPerStep(1000);
    //ts.setRestitutionModel(SemiExplicitEulerTimeStepper::Poisson);
#endif
    ts.initialize(s);
    viz.report(ts.getState());
    printf("Hit ENTER to simulate ... (ESC to quit)\n");
    userInput->waitForAnyUserInput(); userInput->clear();

    const double startCPU  = cpuTime(), startTime = realTime();

    try {
        ts.stepTo(1.5); // TODO Infinity); // RUN
        std::cout << "y=" << ts.getState().getY() << std::endl;
        std::cout << "normRMS: " << (orig - ts.getState().getY()).normRMS() << std::endl;

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        biped.realize(ts.getState());
        std::cout << "y=" << ts.getState().getY() << std::endl;
        std::cout << "ydot=" << ts.getState().getYDot() << std::endl;
        throw;
    }

    // dumpIntegratorStats(startCPU, startTime, integ);

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
	return 0;
}


//==============================================================================
//                           OUTPUT REPORTER
//==============================================================================
OutputReporter::OutputReporter(const Biped& biped,
                               Real            interval)
:   PeriodicEventReporter(interval), m_biped(biped) {}

void OutputReporter::handleEvent(const State& state) const {
    //printf("OutputReporter @t=%g:\n", state.getTime());
    //Real fLeft, fRight;
    //m_biped.findContactForces(state, fLeft, fRight);
    //printf("  forces left=%g, right=%g\n", fLeft, fRight);
}

//==============================================================================
//                           USER INPUT HANDLER
//==============================================================================
UserInputHandler::UserInputHandler(Visualizer::InputSilo& silo,
                                   Real                   interval)
:   PeriodicEventHandler(interval), m_silo(silo) {}

void UserInputHandler::handleEvent(State& state, Real accuracy,
                                   bool& shouldTerminate) const  {
    while (m_silo.isAnyUserInput()) {
        unsigned key, modifiers;
        while (m_silo.takeKeyHit(key,modifiers))
            if (key == Visualizer::InputListener::KeyEsc) {
                shouldTerminate = true;
                m_silo.clear();
                return;
            }
    }
}

//==============================================================================
//                        DUMP INTEGRATOR STATS
//==============================================================================
namespace {
void dumpIntegratorStats(double startCPU, double startTime,
                         const Integrator& integ) {
    std::cout << "DONE: Simulated " << integ.getTime() << " seconds in " <<
        realTime()-startTime << " elapsed s, CPU="<< cpuTime()-startCPU << "s\n";
    #ifdef ANIMATE
    printf("***CAUTION: CPU time not accurate when animation is enabled.\n");
    #endif

    const int evals = integ.getNumRealizations();
    std::cout << "\nUsed "  << integ.getNumStepsTaken() << " steps, avg step="
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms "
        << (1000*integ.getTime())/evals << "ms/eval\n";

    printf("Used Integrator %s at accuracy %g:\n",
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n",  integ.getNumStepsTaken(),
                                          integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n",     integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(),
                                          integ.getNumProjections());
    // Ratio of failed steps to successful ones is crude stiffness measure.
    printf("\nEstimated stiffness: %g\n",
       (double)integ.getNumErrorTestFailures()/integ.getNumStepsTaken());
}
}
