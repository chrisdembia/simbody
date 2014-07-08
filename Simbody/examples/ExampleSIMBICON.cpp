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


#include <cassert>
#include <vector>

using namespace std;
using namespace SimTK;

// #define RIGID_CONTACT

// TODO remove
const Real RealTimeFactor
    = 1; // try to run in real time

// TODO remove
const int NumActuators = 30;

// SIMBICON

#define STATE_UPD_STEPSIZE 0.0005
//#define RIGID_CONTACT

class SIMBICON : public SimTK::Force::Custom::Implementation
{
public:
	SIMBICON(Biped& biped);

    /// The possible states of the finite state machine.
    enum SIMBICONState {
        UNKNOWN = -1, // initial state at the beginning of a simulation.
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
                    Stage::Acceleration, new Value<SIMBICONState>(UNKNOWN),
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

	void computeSecondaryStateVals(const SimTK::State& s,
        SimTK::Real lForce, SimTK::Real rForce);

private:

    /// Determine if the SIMBICON state of the state machine has changed.
    void updateSIMBICONState(const State& s) const;

    void setSIMBICONState(const State& s, SIMBICONState simbiconState) const
    {
        m_forces.updDiscreteVarUpdateValue(s, m_simbiconStateIndex) =
            Value<SIMBICONState>(simbiconState);
        m_forces.markDiscreteVarUpdateValueRealized(s, m_simbiconStateIndex);

        // TODO
		if (simbiconState == STATE0 || simbiconState == STATE2) {
			// reset the swing thigh orientation when the swing leg changes
			for (int i = 0; i < 2; i++) {
				const_cast<SIMBICON*>(this)->_lastSWTAngle[i] = -100.0;
				const_cast<SIMBICON*>(this)->_curSWTAngle[i] = -100.0;
			}
        }
    }

    Biped& m_biped;
    const GeneralForceSubsystem& m_forces;

    // The SIMBICON controller itself has some state.
    // ----------------------------------------------
    // Written during the realizeTopology step.

    /// The SIMBICON state of the controller for the finite state machine
    /// (e.g., 0, 1, 2, or 3). We read the value of this variable using this
    //index.
    DiscreteVariableIndex m_simbiconStateIndex;

    // We set the value of the SIMBICON state using this index.
    CacheEntryIndex m_simbiconStateCacheIndex;


    void computeControls(const SimTK::State& s, SimTK::Vector& controls) const;

	void getSagCorNormals( const SimTK::State& s,
		SimTK::Vec3& sagN, SimTK::Vec3& corN ) const;
	void getUpVectorInGround( const SimTK::State& s, const SimTK::MobilizedBody& b,
		SimTK::Vec3& up ) const;
	void getFrontVectorInGround( const SimTK::State& s, const SimTK::MobilizedBody& b,
		SimTK::Vec3& up ) const;
	void fillInHipJointControls( const SimTK::State& s,
		SimTK::Vector& controls ) const;

	double _lastSWTAngle[2];
	double _curSWTAngle[2];
	double _lastTrunkAngle[2];
	double _curTrunkAngle[2];
	double _lastRFootAngle[2];
	double _curRFootAngle[2];
	double _lastLFootAngle[2];
	double _curLFootAngle[2];
	double _lastPelvisRotation;
    double _curPelvisRotation;

#ifndef RIGID_CONTACT
    double _lastRFootContactForce;
    double _lastLFootContactForce;
    double _curRFootContactForce;
    double _curLFootContactForce;
#endif

};

// Normally SIMBICON has 4 states per gait cycle (i.e. 2 per leg), 2 state is
// simplified and not as realistic.
#define TWO_STATE
// Torque angle to be flat relative to ground (not in original SIMBICON paper)
//#define USE_GLOBAL_ANKLE
// Aim hip towards a global orientation (prevents meandering) (Not implemented
// for rigid contact because it depends on checking force above some maximum.)
//#define USE_GLOBAL_HIPROT
// Drop landing doesn't use controller so you can run with models other than
// the humanoid upon which the controller depends.
//#define DROP_LANDING

namespace { // file-scope symbols

// Controller gains, given by "strength" in N-m/radian. Then kp=strength
// and kd=2*sqrt(strength) for critical damping.
#define USE_ORIG_GAINS
#ifdef USE_ORIG_GAINS
const Real DefaultStrength              = 300;
const Real NeckStrength                 = 100;
const Real BackStrength                 = 300;
const Real HipFlexionAdductionKp  = 1000;  // must be overdamped
const Real HipFlexionAdductionKd  = 100;
const Real HipRotationStrength          = 300;
const Real KneeStrength                 = 300;
const Real ArmFlexionAdductionStrength  = 300;
const Real ArmRotationStrength          = 300;
const Real AnkleFlexionStrength         = 300;
const Real AnkleInversionStrength       = 30;  //300
const Real ToeStrength                  = 30;
#else
const Real DefaultStrength              = 300;  // was 300,30
const Real NeckStrength                 = 20;  // was 100,10
const Real BackStrength                 = 200;  // was 300,30
const Real HipFlexionAdductionStrength  = 500;  // was 1000,100
const Real HipRotationStrength          = 100;  // was 300,30
const Real KneeStrength                 = 300;  // was 300,30
const Real ArmFlexionAdductionStrength  = 10;  // was 300,30
const Real ArmRotationStrength          = 5;  // was 300,30
const Real AnkleFlexionStrength         = 100;  // was 300,30
const Real AnkleInversionStrength       = 50;   // was 300,30
const Real ToeStrength                  = 25;   // was 30,3
#endif

const Vec3 UnitX(1.0, 0.0, 0.0);
const Vec3 UnitY(0.0, 1.0, 0.0);
const Vec3 UnitZ(0.0, 0.0, 1.0);

// Convert muscle strength into critically damped control gains.
void calcGainsFromStrength(double strength, double& kp, double& kd) {
    kp = strength;
    kd = 2*std::sqrt(strength);
}

double clamp( double x, double maxTorque ) {
	if (x > maxTorque)
		x = maxTorque;

	if (x < -maxTorque)
		x = -maxTorque;

	return x;
}
}

SIMBICON::SIMBICON(Biped& biped) : m_biped(biped),
    m_forces(m_biped.getForceSubsystem()) {
	for (int i = 0; i < 2; i++) {
		_lastSWTAngle[i] = -100.0;
		_curSWTAngle[i] = -100.0;
		_lastTrunkAngle[i] = -100.0;
		_curTrunkAngle[i] = -100.0;
		_lastRFootAngle[i] = -100.0;
		_curRFootAngle[i] = -100.0;
		_lastLFootAngle[i] = -100.0;
		_curLFootAngle[i] = -100.0;
	}
	_lastPelvisRotation = -100.0;
	_curPelvisRotation = -100.0;

}

void SIMBICON::calcForce(const State&         s,
                         Vector_<SpatialVec>& /*bodyForces*/,
                         Vector_<Vec3>&       /*particleForces*/,
                         Vector&              mobilityForces) const
{
    updateSIMBICONState(s);
    Vector controls(NumActuators);
    computeControls(s,controls);

    for (int i=0; i<NumActuators; ++i)
        mobilityForces[m_biped.getUIndex(Biped::Coordinate(i))] = controls[i];
}


// Project the pelvis z (right) and x (forward) directions onto the Ground
// (x-z) plane, and normalize. Projected z is the normal to the sagittal plane;
// projected x is the normal to the coronal plane.
void SIMBICON::getSagCorNormals(const State& s, Vec3& sagN, Vec3& corN ) const {
	const MobilizedBody& pelvis = m_biped.getBody(Biped::pelvis);

    sagN = Vec3(pelvis.getBodyRotation(s).z());
    corN = Vec3(pelvis.getBodyRotation(s).x());
	sagN[YAxis] = 0; // project to y==0
	sagN = sagN.normalize();
	corN[YAxis] = 0; // project to y==0
	corN = corN.normalize();
}

void SIMBICON::getUpVectorInGround( const State& s, const MobilizedBody& b,
	Vec3& up ) const {
	up = Vec3(b.getBodyRotation(s).y());
}

void SIMBICON::getFrontVectorInGround( const State& s, const MobilizedBody& b,
	Vec3& front ) const {
	front = Vec3(b.getBodyRotation(s).x());
}

void SIMBICON::fillInHipJointControls( const State& s, Vector& controls ) const {
	const SimbodyMatterSubsystem& matter = m_biped.getMatterSubsystem();
	const SIMBICONState simbiconState = getSIMBICONState(s);
	Vec3 sagN, corN;

    // Sagittal and coronal plane normals are right and front directions of
    // the pelvis, projected on the ground plane.
	getSagCorNormals(s, sagN, corN);

	int swh = Biped::hip_r_flexion;             // swing hip
	int sth = Biped::hip_l_flexion;             // stance hip
    MobilizedBody ankle = m_biped.getBody(Biped::foot_l); // stance ankle

    if (simbiconState == STATE2 || simbiconState == STATE3) { // right stance
		swh = Biped::hip_l_flexion;
		sth = Biped::hip_r_flexion;
        ankle = m_biped.getBody(Biped::foot_r);
	}
	int swhc = swh - 1; // coronal plane (hip adduction)
	int sthc = sth - 1; // stance hip adduction
	int sta = sth + 4;  // stance ankle dorsiflexion

    const Transform& X_PF = ankle.getInboardFrame(s);
    Vec3 ankleLocInParent = X_PF.p();
    Vec3 ankleLoc = ankle.getParentMobilizedBody()
                         .findStationLocationInGround(s,ankleLocInParent);

	Vec3 com = matter.calcSystemMassCenterLocationInGround(s);
	Vec3 d = com - ankleLoc;
	Vec3 v_com = matter.calcSystemMassCenterVelocityInGround(s);

	double d_sag = dot(corN, d);
	double v_sag = dot(corN, v_com);
	double d_cor = dot(sagN, d);
	double v_cor = dot(sagN, v_com);

	double thetad = 0.5;
#ifndef TWO_STATE
    if (simbiconState == STATE1 || simbiconState == STATE3)
		thetad = -0.1;
	}
#endif
    double kp, kd; // position, derivative gains for hip flex/adduction
    #ifdef USE_ORIG_GAINS
    kp = HipFlexionAdductionKp; kd = HipFlexionAdductionKd;
    #else
    calcGainsFromStrength(HipFlexionAdductionStrength, kp, kd);
    #endif
	double cd = 0.2; // global tipping feedback
	double cv = 0.2;

	double trunkAngleVelEst[2] = {0, 0};
	double SWTAngleVelEst[2] = {0, 0};
	double RFootAngleVelEst[2] = {0, 0};
	double LFootAngleVelEst[2] = {0, 0};
	double PelvisRotationVelEst = 0;
	if (_lastTrunkAngle[0] > -100) {
		// check there's a valid value for the _last*,
		// otherwise just use 0 for vel
		for (int i = 0; i < 2; i++) {
			trunkAngleVelEst[i] =
				(_curTrunkAngle[i] - _lastTrunkAngle[i])/STATE_UPD_STEPSIZE;
			RFootAngleVelEst[i] =
				(_curRFootAngle[i] - _lastRFootAngle[i])/STATE_UPD_STEPSIZE;
			LFootAngleVelEst[i] =
				(_curLFootAngle[i] - _lastLFootAngle[i])/STATE_UPD_STEPSIZE;
			if (_lastSWTAngle[i] > -100) {
				SWTAngleVelEst[i] =
					(_curSWTAngle[i] - _lastSWTAngle[i])/STATE_UPD_STEPSIZE;
			}
		}
		PelvisRotationVelEst =
			(_curPelvisRotation - _lastPelvisRotation)/STATE_UPD_STEPSIZE;
	}


	// sign change is needed for one of the stance legs in the coronal plane
	double sign = 1;
    if (simbiconState == STATE0 || simbiconState == STATE1) // left stance
		sign = -1;

	controls[swh] = clamp(  kp*(thetad + (cd*d_sag + cv*v_sag) - _curSWTAngle[0])
                          - kd*SWTAngleVelEst[0], kp);
	controls[swhc] = sign*(clamp(  kp*(0.0 + (cd*d_cor + cv*v_cor) - _curSWTAngle[1])
                                 - kd*SWTAngleVelEst[1], kp));

	// use stance hip to control the trunk
	controls[sth] =  -kp*(0. - _curTrunkAngle[0]) + kd*trunkAngleVelEst[0];
	controls[sth] -= controls[swh];
	controls[sth] = clamp(controls[sth], kp);

	controls[sthc] = sign*(kp*(0. - _curTrunkAngle[1]) - kd*trunkAngleVelEst[1]);
	controls[sthc] -= controls[swhc];
	controls[sthc] = clamp(controls[sthc], kp);

#ifdef USE_GLOBAL_ANKLE
    double kpaflex, kdaflex, kpainv, kdainv; // gains for ankle flex, inversion
    calcGainsFromStrength(AnkleFlexionStrength, kpaflex, kdaflex);
    calcGainsFromStrength(AnkleInversionStrength, kpainv, kdainv);
	controls[ankle_r_dorsiflexion] =
        clamp(  kpaflex*(0. - _curRFootAngle[0])
              - 0.*kdaflex*RFootAngleVelEst[0], kpaflex);
	controls[ankle_r_inversion] =
        clamp( -kpainv*(0. - _curRFootAngle[1])
              + kdainv*RFootAngleVelEst[1], kpainv);
	controls[ankle_l_dorsiflexion] =
        clamp(  kpaflex*(0. - _curLFootAngle[0])
              - 0.*kdaflex*LFootAngleVelEst[0], kpaflex);
	controls[ankle_l_inversion] =
        clamp(  kpainv*(0. - _curLFootAngle[1])
              - kdainv*LFootAngleVelEst[1], kpainv);
#endif

#ifdef USE_GLOBAL_HIPROT
    double kphrot, kdhrot; // gains for hip rotation
    calcGainsFromStrength(HipRotationStrength, kphrot, kdhrot);
	if (   (sth == hip_r_flexion && _curRFootContactForce > 100)
        || (sth == hip_l_flexion  && _curLFootContactForce > 100))
    {
		controls[sth+1] = clamp(sign*( -kphrot*(0. - _curPelvisRotation)
                                      + kdhrot*PelvisRotationVelEst ),
                                kphrot);
	}
#endif

}

void SIMBICON::computeControls(const State& s, Vector& controls) const
{
	const SIMBICONState simbiconState = getSIMBICONState(s);
#ifdef DROP_LANDING
	for (int i = 0; i < controls.size(); i++)
			controls[i] = 0.0;
	return;
#else
	int swh = Biped::hip_r_flexion;
	int sth = Biped::hip_l_flexion;

    if (simbiconState == STATE2 || simbiconState == STATE3) // right stance
    {
		swh = Biped::hip_l_flexion;
		sth = Biped::hip_r_flexion;
	}
	int swk = swh + 2; // swing knee
	int swa = swh + 4; // swing ankle
	int stk = sth + 2; // stance knee
	int sta = sth + 4; // stance ankle
	for (int i = 0; i < NumActuators; i++) {
        double kp, kd;          // position gain, derivative gain
        calcGainsFromStrength(DefaultStrength, kp, kd);
		double thetad = 0.0;                    // desired angle

		if (i == Biped::neck_extension || i == Biped::neck_bending || i == Biped::neck_rotation) {
            calcGainsFromStrength(NeckStrength, kp, kd);
		}
		else if (i == Biped::back_tilt || i == Biped::back_list || i == Biped::back_rotation) {
            calcGainsFromStrength(BackStrength, kp, kd);
		}
		else if (   i == Biped::shoulder_r_flexion   || i == Biped::shoulder_l_flexion
                 || i == Biped::shoulder_r_adduction || i == Biped::shoulder_l_adduction
                 || i == Biped::elbow_r_flexion      || i == Biped::elbow_l_flexion) {
            calcGainsFromStrength(ArmFlexionAdductionStrength, kp, kd);
		}
		else if (   i == Biped::shoulder_r_rotation || i == Biped::shoulder_l_rotation
                 || i == Biped::elbow_r_rotation    || i == Biped::elbow_l_rotation) {
            calcGainsFromStrength(ArmRotationStrength, kp, kd);
		}
		else if (i == Biped::hip_r_rotation || i == Biped::hip_l_rotation) {
            calcGainsFromStrength(HipRotationStrength, kp, kd);
		}
		else if (i == Biped::knee_r_extension || i == Biped::knee_l_extension) {
            calcGainsFromStrength(KneeStrength, kp, kd);
		}
        else if (i == Biped::ankle_r_dorsiflexion || i == Biped::ankle_l_dorsiflexion) {
            calcGainsFromStrength(AnkleFlexionStrength, kp, kd);
        }
        else if (i == Biped::ankle_r_inversion || i == Biped::ankle_l_inversion) {
            calcGainsFromStrength(AnkleInversionStrength, kp, kd);
        }
		else if (i == Biped::mtp_r_dorsiflexion || i == Biped::mtp_l_dorsiflexion) { // toe
            calcGainsFromStrength(ToeStrength, kp, kd);
		}

		if (simbiconState >= STATE0) {
			if (i == swk)
				thetad = -1.1;
			else if (i == swa)
				thetad = 0.6;
			else if (i == stk)
				thetad = -0.05;

            #ifndef TWO_STATE
            if (simbiconState == STATE1 || simbiconState ==  STATE3) {
				if (i == swk)
					thetad = -0.05;
				else if (i == swa)
					thetad = 0.15;
				else if (i == stk)
					thetad = -0.1;
			}
            #endif
		}

        const Real qi = s.getQ()[m_biped.getQIndex(Biped::Coordinate(i))];
        const Real ui = s.getU()[m_biped.getUIndex(Biped::Coordinate(i))];
		controls[i] = clamp(kp*(thetad - qi) - kd*ui, kp);
	}

	if (simbiconState >= STATE0) {
		fillInHipJointControls(s, controls);
    }

	Vec3 com = m_biped.getMatterSubsystem().calcSystemMassCenterLocationInGround(s);
	if (com[1] < 0.7) {
		for (int i = 0; i < controls.size(); i++)
			controls[i] = 0.0;
	}
	return;
#endif
}

void SIMBICON::
computeSecondaryStateVals(const State& s, Real lForce, Real rForce) {
	const MobilizedBody& pelvis = m_biped.getBody(Biped::pelvis);
	const SIMBICONState simbiconState = getSIMBICONState(s);

#ifndef RIGID_CONTACT
    _lastRFootContactForce = _curRFootContactForce;
    _lastLFootContactForce = _curLFootContactForce;
    _curRFootContactForce = rForce;
    _curLFootContactForce = lForce;
#endif

	Vec3 upThigh;
    if (simbiconState == STATE0 || simbiconState == STATE1)
		getUpVectorInGround(s, m_biped.getBody(Biped::thigh_r), upThigh);
    else if (simbiconState == STATE2 || simbiconState == STATE3)
		getUpVectorInGround(s, m_biped.getBody(Biped::thigh_l), upThigh);

	Vec3 upPelvis;
	getUpVectorInGround(s, pelvis, upPelvis);

	Vec3 sagN, corN;
	getSagCorNormals(s, sagN, corN);

	Vec3 XinSag = cross(UnitY, sagN);
	Vec3 ZinCor = -cross(UnitY, corN);

	// Store the current value, use these for velocity estimation
	for (int i = 0; i < 2; i++) {
		_lastSWTAngle[i] = _curSWTAngle[i];
		_lastTrunkAngle[i] = _curTrunkAngle[i];
		_lastLFootAngle[i] = _curLFootAngle[i];
		_lastRFootAngle[i] = _curRFootAngle[i];
	}
	_lastPelvisRotation = _curPelvisRotation;

	// Update trunk and swing thigh global orientations in the saggital and
	// coronal planes (by projecting the up vectors of the bodies in to the
	// planes and calculate angles)
	Vec3 projUpThigh = upThigh - dot(upThigh, sagN)*sagN;
	_curSWTAngle[0] =
		acos(dot(projUpThigh.normalize(), XinSag)) - Pi/2;

	Vec3 projUpPelvis = upPelvis - dot(upPelvis, sagN)*sagN;
	_curTrunkAngle[0] =
		acos(dot(projUpPelvis.normalize(), XinSag)) - Pi/2;

	Vec3 projUpPelvisCor = upPelvis - dot(upPelvis, corN)*corN;
	_curTrunkAngle[1] =
		acos(dot(projUpPelvisCor.normalize(), ZinCor)) - Pi/2;

	Vec3 projUpThighCor = upThigh - dot(upThigh, corN)*corN;
	_curSWTAngle[1] =
		acos(dot(projUpThighCor.normalize(), ZinCor)) - Pi/2;


#ifdef USE_GLOBAL_ANKLE
	const MobilizedBody& rFoot = m_biped.getBody(Biped::foot_r);
	const MobilizedBody& lFoot = m_biped.getBody(Biped::foot_l);

	Vec3 upRFoot, upLFoot;
	getUpVectorInGround(s, rFoot, upRFoot);
	getUpVectorInGround(s, lFoot, upLFoot);

	Vec3 projUpRFoot = upRFoot - dot(upRFoot, sagN)*sagN;
	_curRFootAngle[0] =
		acos(dot(projUpRFoot.normalize(), XinSag)) - Pi/2;

	Vec3 projUpLFoot = upLFoot - dot(upLFoot, sagN)*sagN;
	_curLFootAngle[0] =
		acos(dot(projUpLFoot.normalize(), XinSag)) - Pi/2;

	Vec3 projUpRFootCor = upRFoot - dot(upRFoot, corN)*corN;
	_curRFootAngle[1] =
		acos(dot(projUpRFootCor.normalize(), ZinCor)) - Pi/2;

	Vec3 projUpLFootCor = upLFoot - dot(upLFoot, corN)*corN;
	_curLFootAngle[1] =
		acos(dot(projUpLFootCor.normalize(), ZinCor)) - Pi/2;

#endif

#ifdef USE_GLOBAL_HIPROT
	Vec3 frontPelvis;
	getFrontVectorInGround(s, pelvis, frontPelvis);

	frontPelvis[1] = 0.0; // project to 2d
	_curPelvisRotation = acos(dot(frontPelvis.normalize(), Vec3(0.0, 0.0, 1.0))) - Pi/2;
#endif

}

void SIMBICON::updateSIMBICONState(const State& s) const
{
    // TODO remove this line?
    //m_biped.realize(s, Stage::Dynamics);

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
        case UNKNOWN:

            // Entering left stance.
            if (rContact) {
                setSIMBICONState(s, STATE0);
            }
            // Entering right stance.
            else if (lContact) {
                setSIMBICONState(s, STATE2);
            }
            break;

        // Left stance.
        case STATE0:

            // Stay in this state for \delta t seconds. TODO
            if (duration > 0.3) {
                setSIMBICONState(s, STATE1);
            }
            // Already entered right stance; skip STATE1.
            else if (rContact && duration > 0.1) { // TODO min duration.
                setSIMBICONState(s, STATE2);
            }
            break;

        // Right foot strike.
        case STATE1:

            // Stay in this state until the right foot makes contact.
            if (rContact && duration > 0.1) { // TODO min duration.
                setSIMBICONState(s, STATE2);
            }
            break;

        // Right stance.
        case STATE2:

            // Stay in this state for \delta t seconds. TODO
            if (duration > 0.3) {
                setSIMBICONState(s, STATE3);
            }
            // Already entered right stance; skip STATE3.
            else if (lContact && duration > 0.1) {
                setSIMBICONState(s, STATE0);
            }
            break;

        // Left foot strike.
        case STATE3:

            // Stay in this state until the left foot makes contact.
            if (lContact && duration > 0.1) {
                setSIMBICONState(s, STATE0);
            }
            break;
    }

    // TODO const cast.
    const_cast<SIMBICON*>(this)->computeSecondaryStateVals(s, lForce, rForce);

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
    ShowContact(const Biped& unis)
    :   m_unis(unis) {}

    void generateDecorations(const State&                state,
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        DecorativeText energy;
        energy.setIsScreenText(true);
        bool left, right;
        m_unis.findContactStatus(state, left, right);
        std::string str = "";
        if (left) {
            str += "left on,";
        }
        if (right) {
            str += "right on";
        }
        energy.setText(str);
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
};

//==============================================================================
// MAIN FUNCTION
//==============================================================================
int main(int argc, char **argv)
{
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

    State s;
    biped.initialize(s);

    printf("Act: u\n");
    for (int i=0; i < NumActuators; ++i) {
        printf("%2d: %d\n", i, int(biped.getUIndex(Biped::Coordinate(i))));
    }

    //biped.toes_r.lockAt(s, .2); biped.toes_l.lockAt(s, .2);
    //biped.foot_r.lockAt(s, Vec2(.2,0)); biped.foot_l.lockAt(s, Vec2(.2,0));
    biped.realize(s, Stage::Instance);

    printf("SIMBICON 3D:\n");
    printf("%d bodies, %d mobilities, -%d constraint equations -%d motions\n",
        biped.getMatterSubsystem().getNumBodies(), s.getNU(), s.getNMultipliers(),
        biped.getMatterSubsystem().getKnownUDotIndex(s).size());

    biped.getBody(Biped::trunk).setQToFitTranslation(s, Vec3(0,1.5,0));
    biped.getBody(Biped::trunk).setUToFitLinearVelocity(s, Vec3(1,0,0));
    viz.report(s);

#ifdef RIGID_CONTACT
    viz.addDecorationGenerator(new ShowContact(biped));
#endif

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

#ifdef RIGID_CONTACT
    int timeStep = 0;
    do {
        if (timeStep % 10 == 0) {
            viz.report(ts.getState());
        }
        ts.stepTo(ts.getState().getTime() + STATE_UPD_STEPSIZE);

        State& s = const_cast<State&>(ts.getState());

        biped.realize(s, Stage::Dynamics);

        bool lContact, rContact;
        biped.findContactStatus(s, lContact, rContact);


    #ifndef DROP_LANDING
        int curState = simctrl->getState();
        double duration = s.getTime() - simctrl->getStateStartTime();

        if (curState < 0) {
            if (rContact)
                simctrl->setState(2, s.getTime());
            else if (lContact)
                simctrl->setState(0, s.getTime());
        }
        else if (curState == 0) {
            if (duration > 0.3) {
                simctrl->setState(1, s.getTime());
    //        	std::cout << "Trans 0 -> 1" << std::endl;;
            }
            else if (rContact && duration > 0.1) {
                simctrl->setState(2, s.getTime());
      //      	std::cout << "Trans 0 -> 2" << std::endl;;
            }
        }
        else if (curState == 1) {
            if (rContact && duration > 0.1) {
                simctrl->setState(2, s.getTime());
        //    	std::cout << "Trans 1 -> 2" << std::endl;;
            }
        }
        else if (curState == 2) {
            if (duration > 0.3) {
                simctrl->setState(3, s.getTime());
          //  	std::cout << "Trans 2 -> 3" << std::endl;;
            }
            else if (lContact && duration > 0.1) {
                simctrl->setState(0, s.getTime());
           // 	std::cout << "Trans 2 -> 0" << std::endl;;
            }
        }
        else if (curState == 3) {
            if (lContact && duration > 0.1) {
                simctrl->setState(0, s.getTime());
           // 	std::cout << "Trans 3 -> 0" << std::endl;;
            }
        }

        simctrl->computeSecondaryStateVals(s, 0, 0);

#endif


        timeStep += 1;
    } while (ts.getTime() < 10);
#else
        ts.stepTo(Infinity); // RUN
#endif

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
