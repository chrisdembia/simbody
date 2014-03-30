/* -------------------------------------------------------------------------- *
 *                 Simbody(tm) Example: SIMBICON                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Jack Wang, Tim Dorn                                               *
 * Contributors: Chris Dembia                                                 *
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

/*

SIMBICON stands for Simple Biped Controller. It is described at
https://www.cs.ubc.ca/~van/papers/Simbicon.htm.


Brief description of the controller
===================================

The controller uses a "finite state machine" / "pose control graph": the biped
is always in some state, and all states have a target pose. The number of
states is typically 4.  The joints are torque-actuated, and the torques are
given by proportional-derivative (PD) control laws that minimize the error of a
given joint angle from that in the target pose:

    \tau = k_p (\theta_{des} - \theta) - k_d \dot{\theta}

There are some exceptions to this general control scheme. Namely, the target
angles for the torso and "swing-leg femur" are specified relative to the global
frame. These two angles are used to balance of the biped.

The parameters of the controller are used to define the motion that the biped
executes. The motions the controller can achieve consist of walking, running,
and variants on both (e.g, 'high-step walk'). The paramters are:

for each state:
    * state dwell duration
    * position balance feedback coefficient
    * velocity balance feedback coefficient
    * torso target angle
    * swing-hip target angle
    * swing-knee target angle
    * swing-ankle target angle
    * stance-knee target angle
    * stance-ankle target angle

See [1], particularly Table 1 of [1], for more information.


How to navigate this file
=========================
To see how X happens, look in method Y:

X                                        Y
---------------------------------------  ---------------------------------------
construction of the model                Biped::Biped()
control laws                             SIMBICON::calcForce()
logic to change state in state machine   SIMBICON::StateHandler::handleEvent()


Notes
=====

I (Chris Dembia) wrote this code pretty much by copying code that Michael
Sherman gave me. The code that Michael Sherman gave me was originally written
by Jack Wang and Tim Dorn.

Yin, KangKang, Kevin Loken, and Michiel van de Panne. "SIMBICON: Simple biped
locomotion control." ACM Transactions on Graphics (TOG). Vol. 26. No. 3. ACM,
2007.

*/

// Interval for updating the SIMBICON state (not same as SimTK State) of the
// biped.
#define SIMBICON_STATE_UPDATE_STEPSIZE 0.0005

#include "Simbody.h"

using namespace SimTK;

const Real D2R = Pi / 180.0;

//==============================================================================
// BIPED
//==============================================================================
/// The MultibodySystem that we will control with the SIMBICON controller
/// (though the system is independent of the controller, and coulc conceivably
/// be controlled by a different controller).
class Biped : public MultibodySystem {
public:
    Biped();

    const SimbodyMatterSubsystem& getMatterSubsystem() const {return m_matter;}
    SimbodyMatterSubsystem& updMatterSubsystem() {return m_matter;}

    const GeneralForceSubsystem& getForceSubsystem() const {return m_forces;}
    GeneralForceSubsystem& updForceSubsystem() {return m_forces;}

private:

    /// For convenience in building the model.
    /// Lower bound and upper bound should be in degrees.
    void addMobilityLinearStop(const MobilizedBody& mobod,
            unsigned int idx,
            Real lowerBound, Real upperBound,
            Real stopStiffness=1000 / D2R, Real stopDissipation=1.0)
    {
        Force::MobilityLinearStop(m_forces, mobod, MobilizerQIndex(idx),
                stopStiffness, stopDissipation,
                lowerBound * D2R, upperBound * D2R);
    }

    // Subsystems.
    SimbodyMatterSubsystem       m_matter;
    GeneralForceSubsystem        m_forces;
    ContactTrackerSubsystem      m_tracker;
    CompliantContactSubsystem    m_contact;

    // Mobilized bodies.
    MobilizedBody::Free          m_trunk;
    MobilizedBody::Gimbal        m_head;
    MobilizedBody::Gimbal        m_pelvis;
    MobilizedBody::Gimbal        m_upperarm_r;
    MobilizedBody::Universal     m_lowerarm_r;
    MobilizedBody::Weld          m_hand_r;
    MobilizedBody::Gimbal        m_upperarm_l;
    MobilizedBody::Universal     m_lowerarm_l;
    MobilizedBody::Weld          m_hand_l;
    MobilizedBody::Gimbal        m_thigh_r;
    MobilizedBody::Pin           m_shank_r;
    MobilizedBody::Universal     m_foot_r;
    MobilizedBody::Pin           m_toes_r;
    MobilizedBody::Gimbal        m_thigh_l;
    MobilizedBody::Pin           m_shank_l;
    MobilizedBody::Universal     m_foot_l;
    MobilizedBody::Pin           m_toes_l;
};

//==============================================================================
// SIMBICON
//==============================================================================
/// The actual controller that specifies the torques to apply to the Biped.
class SIMBICON : public Force::Custom::Implementation {
public:
    SIMBICON(Biped& biped)
    :   m_biped(biped)
    {
        m_biped.addEventHandler(
                new StateHandler(biped, *this,
                    SIMBICON_STATE_UPDATE_STEPSIZE));
    }
    void calcForce(const State&         state,
                   Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>&       particleForces,
                   Vector&              mobilityForces) const
                   OVERRIDE_11
    {
        // TODO
        // Which leg is in stance?
        // TODO double stance?
        // Pose graph control.
        // ===================
        // Apply PD control to most joints to track target pose.
        // 
        /*
        mobilityForces[NECK_EXTENSION] = calcPDControl(NECK, NECK_EXTENSION);
        mobilityForces[NECK_BENDING] = calcPDControl(NECK, NECK_BENDING);
        mobilityForces[NECK_ROTATION] = calcPDControl(NECK, NECK_ROTATION);
        mobilityForces[KNEE_STANCE_EXTENSION] = calcPDControl(KNEE, KNEE_STANCE_EXTENSION,
                    stk[m_simbicon_state]);
        }

        biped.addInForce("neck", calcPDControl(NECK,

        uIdx = m_biped.getUIndex("neck_extension");
        mobilityForces[uIdx] = calcPDControl(s, "neck", "neck_extension");

        mobilityForces[m_biped.getUIndex("neck_extension")] = calcPDControl(s,
            "neck", "neck_extension");


        addInPD(s, "neck", "neck_extension", 0.0, mobilityForces);
        addInPD(s, "neck", "neck_bending", 0.0, mobilityForces);
        addInPD(s, "neck", "neck_rotation", 0.0, mobilityForces);
        addInPD(s, "knee", stanceKneeRotation, stk[m_simbicon_state],
            mobilityForces);



            m_simbicon_state
        biped.addInForce("neck", calcPDControl(NECKmobilityForces);
                mobilityForces[uindex(uname)] += force;

                */

    }

    Real calcPotentialEnergy(const State& state) const OVERRIDE_11
    {
        return 0;
    }

    /// Updates the SIMBICON state of biped for the finite state machine.
    class StateHandler : public PeriodicEventHandler {
    public:
        StateHandler(Biped&       biped,
                     SIMBICON&    simbicon,
                     Real         interval)
        :   PeriodicEventHandler(interval),
            m_biped(biped), m_simbicon(simbicon) {}
        void handleEvent(State& s, Real accuracy, bool& shouldTerminate) const
        {
            // TODO
        }
    private:
        Biped& m_biped;
        SIMBICON& m_simbicon;
    };
private:
    Biped& m_biped;
/*
    Real calcPDControl(const State& s, string strengthGroup, string dof, Real thetad)
{
    kp = positionGains[strengthGroup];
    kd = speedGains[strengthGroup];
    q = m_biped.getQ(s, dof);
    u = m_biped.getU(s, dof);
    return clamp(kp * (thetad - q) - kd * u);
}

    void addInPDControl(const State& s, string strengthGroup, string dof, Real thetad)
{
    m_biped.addInForce(dof, calcPDControl(...), mobilityForces);
}
*/

};

//==============================================================================
// OUTPUT REPORTER
//==============================================================================
/// This is a periodic event handler that we use to print information to the
/// screen during the simulation.
class OutputReporter : public PeriodicEventReporter {
public:
    OutputReporter(const Biped& biped, Real interval)
    :   PeriodicEventReporter(interval), m_biped(biped) {}

    void handleEvent(const State& state) const OVERRIDE_11
    {
        // TODO
    }
private:
    const Biped& m_biped;
};

//==============================================================================
// MAIN FUNCTION
//==============================================================================
int main(int argc, char **argv)
{
    // Parameters.
    // -----------

    // Set this < 1 for slow motion, > 1 for speedup.
    // 1.0 indicates real time.
    const Real realTimeFactor = 1.0;

    // Create the system.
    Biped biped;

    // Periodically report state-related quantities to the screen.
    biped.addEventReporter(new OutputReporter(biped, 0.01));

    // Add controller to the force system. See Doxygen for Force::Custom.
    // The name of this Force is irrelevant.
    Force::Custom simbicon1(biped.updForceSubsystem(), new SIMBICON(biped));

    // Initialize the system.
    State state = biped.realizeTopology();
    biped.realizeModel(state);

    // TODO fillInActuatorMap

    biped.realize(state, Stage::Instance);

    // Start up the visualizer.
    Visualizer viz(biped);
    viz.setShowSimTime(true);
    viz.setShowFrameRate(true);

    // Update the visualization at a framerate of 1/interval.
    biped.addEventReporter(new Visualizer::Reporter(viz, realTimeFactor/30));

    // Show the biped.
    viz.report(state);

	return 0;
}


//==============================================================================
// BIPED
//==============================================================================
Biped::Biped()
    : m_matter(*this), m_forces(*this),
      m_tracker(*this), m_contact(*this, m_tracker)
{

    m_matter.setShowDefaultGeometry(false);

    //--------------------------------------------------------------------------
    //                          Constants, etc.
    //--------------------------------------------------------------------------

    // Properties of joint stops.
    const Real stopStiffness = 1000/D2R; // N-m/deg -> N-m/radian
    const Real stopDissipation = 1;

    // Contact parameters.
    // -------------------
    // slide -> stick velocity.
    const Real transitionVelocity = .03;

    // Friction coefficients.
    const Real mu_s = 5;
    const Real mu_d = 5;
    const Real mu_v = 0;

    // Contact spheres attached to feet.
    const Real contactSphereRadius = 1.5 * .01;

    // Contact material: rubber.
    const Real rubber_density = 1100.;  // kg/m^3
    const Real rubber_young   = 0.01e9; // pascals (N/m)
    const Real rubber_poisson = 0.5;    // ratio
    const Real rubber_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
    const Real rubber_dissipation = /*0.005*/1;

    const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,
                                   mu_s,mu_d,mu_v);

    // Miscellaneous.
    // --------------
    // Original OpenSim ellipsoid_center.vtp half dimensions.
    const Vec3 ectr(.03, .12, .03);
    // Original OpenSim sphere.vtp half dimension.
    const Real rad = .5;
    // Original OpenSim block.vtp half dimensions.
    const Vec3 blk(.05,.05,.05);

    // Some convenient rotation matrices for relabeling axes.
    const Rotation Rzmxmy(BodyRotationSequence,  Pi/2, XAxis,  Pi/2, ZAxis);
    const Rotation Rzxy(BodyRotationSequence,   -Pi/2, XAxis, -Pi/2, ZAxis);
    const Rotation Ryzx(BodyRotationSequence,    Pi/2, YAxis,  Pi/2, ZAxis);

    DecorativeSphere originMarker(0.04*rad);
    originMarker.setColor(Vec3(.1,.1,.1));


    //--------------------------------------------------------------------------
    //                          Body information
    //--------------------------------------------------------------------------
    // Trunk.
    // ------
    const Real trunkMass = 19.733716299458440;
    const Vec3 trunkCOM(0,0,0);
    Body trunkInfo(MassProperties(trunkMass, trunkCOM,
                                  Inertia(0.355666710204554,
                                          0.224533650416368,
                                          0.281526683481324,
                                          0.046737850487895,
                                          0.000655032101243,
                                         -0.000572934744554)));

    trunkInfo.addDecoration(Transform(Rotation(.01, ZAxis),
                                      Vec3(-0.025,-0.06,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.6,1.25,1.6)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Transform(Rotation(-.3, ZAxis),
                                      Vec3(-0.01,0.16,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.85,0.55,0.85)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Transform(Ryzx,Vec3(-0.025,0.09,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.43,1)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Vec3(0),originMarker);

    // Head.
    // -----
    const Real headMass = 4.340524803328146;
    const Vec3 headCOM(0.041488177946752, 0.085892319249215, 0);
    Body headInfo(MassProperties(headMass, headCOM,
                                 Inertia( 0.020576741740382,
                                          0.016008984554380,
                                          0.025555859085964,
                                         -0.004554656543977,
                                         -0.000103058383929,
                                          0.000186029116753)
                                 .shiftFromMassCenter(-headCOM,headMass)));
    headInfo.addDecoration(Vec3(0,.11,0),
        DecorativeEllipsoid(Vec3(0.2,0.22,0.2)/2.)
            .setColor(White).setOpacity(1));
    headInfo.addDecoration(Vec3(0),originMarker);

    // Pelvis.
    // -------
    const Real pelvisMass = 13.924855817213411;
    const Vec3 pelvisCOM(0.036907589663647, -0.142772863495411, 0);
    Body pelvisInfo(MassProperties(pelvisMass, pelvisCOM,
                                   Inertia( 0.172382614643800,
                                            0.137961114411544,
                                            0.128551359933154,
                                           -0.010239461806632,
                                            0.001027963710884,
                                           -0.003286514395970)
                                   .shiftFromMassCenter(-pelvisCOM,pelvisMass)));
    pelvisInfo.addDecoration(Transform(Ryzx,Vec3(0.02,-0.1,0.0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(3.3,1.5,2.4)))
            .setColor(White).setOpacity(1));
    pelvisInfo.addDecoration(Vec3(0),originMarker);

    // Upper arm.
    // ----------
    // COM z has opposite sign left to right.
    const Real upperarmMass = 2.070989783095760;
    const Vec3 upperarm_rCOM(0.003289136233947,-0.078058926824158,0.065606556342984);
    Body upperarm_rInfo(MassProperties(upperarmMass, upperarm_rCOM,
                      Inertia(0.014082101994221,
                              0.003604979502976,
                              0.013769380880709)
                      .shiftFromMassCenter(-upperarm_rCOM, upperarmMass)));
    upperarm_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,-0.47,XAxis,0.13,ZAxis),
                    Vec3(0.017,-0.12,0.06)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.2,1)))
            .setColor(White).setOpacity(1));
    upperarm_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 upperarm_lCOM(0.003289136233947,-0.078058926824158,-0.065606556342984);
    Body upperarm_lInfo(MassProperties(upperarmMass, upperarm_lCOM,
                      Inertia(0.014082101994221,
                              0.003604979502976,
                              0.013769380880709)
                      .shiftFromMassCenter(-upperarm_lCOM, upperarmMass)));
    upperarm_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,0.47,XAxis,0.13,ZAxis),
                    Vec3(0.017,-0.12,-0.06)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.2,1)))
            .setColor(White).setOpacity(1));
    upperarm_lInfo.addDecoration(Vec3(0),originMarker);

    // Lower arm.
    // ----------
    // Some signs differ left to right.
    const Real lowerarmMass = 1.106702647320712;
    const Vec3 lowerarm_rCOM(0.031656703591848,-0.089369993258598,0.017231110378866);
    Body lowerarm_rInfo(MassProperties(lowerarmMass, lowerarm_rCOM,
                            Inertia( 0.003846276658463,
                                     0.001704523106360,
                                     0.004819186789386,
                                     0.001553953681336,
                                    -0.000083971410109,
                                     0.000083971410109)
                            .shiftFromMassCenter(-lowerarm_rCOM,lowerarmMass)));
    lowerarm_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,-0.05,XAxis,0.45,ZAxis),
                    Vec3(0.053,-0.111,0.01)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.95,1.15,0.95)))
            .setColor(White).setOpacity(1));
    lowerarm_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 lowerarm_lCOM(0.031656703591848,-0.089369993258598,-0.017231110378866);
    Body lowerarm_lInfo(MassProperties(lowerarmMass, lowerarm_lCOM,
                            Inertia( 0.003846276658463,
                                     0.001704523106360,
                                     0.004819186789386,
                                     0.001553953681336,
                                     0.000083971410109,
                                    -0.000083971410109)
                            .shiftFromMassCenter(-lowerarm_lCOM,lowerarmMass)));
    lowerarm_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,0.05,XAxis,0.45,ZAxis),
                    Vec3(0.053,-0.111,-0.01)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.95,1.15,0.95)))
            .setColor(White).setOpacity(1));
    lowerarm_lInfo.addDecoration(Vec3(0),originMarker);

    // Hand.
    // -----
    // Some signs differ left to right.
    const Real handMass = 0.340742469583528;
    const Vec3 hand_rCOM(0.031681557027587,-0.041582042351409,-0.008872831097566);
    Body hand_rInfo(MassProperties(handMass, hand_rCOM,
                            Inertia( 0.000294382529694,
                                     0.000262531305170,
                                     0.000326233754218,
                                     0.000061772071805,
                                     0.000054050562829,
                                    -0.000036677167634)
                            .shiftFromMassCenter(-hand_rCOM,handMass)));
    hand_rInfo.addDecoration(Transform(
                    Rotation(0.8,ZAxis),
                    Vec3(0.025,-0.025,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.7,0.3,0.7)))
            .setColor(White).setOpacity(1));
    hand_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 hand_lCOM(0.031681557027587,-0.041582042351409, 0.008872831097566);
    Body hand_lInfo(MassProperties(handMass, hand_lCOM,
                            Inertia( 0.000294382529694,
                                     0.000262531305170,
                                     0.000326233754218,
                                     0.000061772071805,
                                    -0.000054050562829,
                                     0.000036677167634)
                            .shiftFromMassCenter(-hand_lCOM,handMass)));
    hand_lInfo.addDecoration(Transform(
                    Rotation(0.8,ZAxis),
                    Vec3(0.025,-0.025,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.7,0.3,0.7)))
            .setColor(White).setOpacity(1));
    hand_lInfo.addDecoration(Vec3(0),originMarker);

    // Thigh.
    // ------
    const Real thighMass = 8.082407914884000;
    const Vec3 thigh_rCOM(0,-0.178920728716523,0.001605747837523);
    Body thigh_rInfo(MassProperties(thighMass, thigh_rCOM,
                            Inertia( 0.116351777130644,
                                     0.030499980412887,
                                     0.122695077900275 )
                            .shiftFromMassCenter(-thigh_rCOM,thighMass)));
    thigh_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence, -.01,XAxis, -.012,ZAxis),
                    Vec3(-0.002,-0.205,0.003)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.3,1.9,1.3)))
            .setColor(White).setOpacity(1));
    thigh_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 thigh_lCOM(0,-0.178920728716523,-0.001605747837523);
    Body thigh_lInfo(MassProperties(thighMass, thigh_lCOM,
                            Inertia( 0.116351777130644,
                                     0.030499980412887,
                                     0.122695077900275 )
                            .shiftFromMassCenter(-thigh_lCOM,thighMass)));
    thigh_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence, .01,XAxis, -.012,ZAxis),
                    Vec3(-0.002,-0.205,-0.003)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.3,1.9,1.3)))
            .setColor(White).setOpacity(1));
    thigh_lInfo.addDecoration(Vec3(0),originMarker);

    // Shank.
    // ------
    const Real shankMass = 3.222323418816000;
    const Vec3 shank_rCOM(0,-0.182765070363067,0.005552190835500);
    Body shank_rInfo(MassProperties(shankMass, shank_rCOM,
                            Inertia( 0.043804477493817,
                                     0.004432595936874,
                                     0.044412873014564 )
                            .shiftFromMassCenter(-shank_rCOM,shankMass)));
    shank_rInfo.addDecoration(Transform(
                    Rotation(0.030,XAxis),
                    Vec3(0.0,-0.21,-0.005)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.27,1.8,1.27)))
            .setColor(White).setOpacity(1));
    shank_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 shank_lCOM(0,-0.182765070363067,-0.005552190835500);
    Body shank_lInfo(MassProperties(shankMass, shank_lCOM,
                            Inertia( 0.043804477493817,
                                     0.004432595936874,
                                     0.044412873014564 )
                            .shiftFromMassCenter(-shank_lCOM,shankMass)));
    shank_lInfo.addDecoration(Transform(
                    Rotation(-0.030,XAxis),
                    Vec3(0.0,-0.21,0.005)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.27,1.8,1.27)))
            .setColor(White).setOpacity(1));
    shank_lInfo.addDecoration(Vec3(0),originMarker);

    // Foot.
    // -----
    const Real footMass = 1.172905458264000;
    const Vec3 foot_rCOM(0.035606945567853,-0.051617802456029,-0.000574057583573);
    const Inertia foot_Ic(0.001313654113256, // central inertia
                          0.003659465029784,
                          0.003847129903106);
    const Real Ifac = 1; // for playing with foot inertia; should be 1
    Body foot_rInfo(MassProperties(footMass, foot_rCOM,
                      (Ifac*foot_Ic).shiftFromMassCenter(-foot_rCOM,footMass)));
    foot_rInfo.addDecoration(Vec3(0.052,-0.043,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(1.65,0.6,0.8)))
            .setColor(White).setOpacity(1));
    foot_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 foot_lCOM(0.035606945567853,-0.051617802456029,0.000574057583573);
    Body foot_lInfo(MassProperties(footMass, foot_lCOM,
                      (Ifac*foot_Ic).shiftFromMassCenter(-foot_lCOM,footMass)));
    foot_lInfo.addDecoration(Vec3(0.052,-0.043,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(1.65,0.6,0.8)))
            .setColor(White).setOpacity(1));
    foot_lInfo.addDecoration(Vec3(0),originMarker);

    // Toes.
    // -----
    const Real toesMass = 0.20;
    // This inertia looks artificially large.
    const Inertia toes_Ic(0.087132432150,0.0174264864299,0.0871324321496);
    const Vec3 toes_rCOM(0.023716003435794,-0.001184334594970,-0.002484544347746);
    Body toes_rInfo(MassProperties(toesMass, toes_rCOM,
                            toes_Ic.shiftFromMassCenter(-toes_rCOM,toesMass)));
    toes_rInfo.addDecoration(Vec3(0.03,0.014,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(0.6,0.3,0.8)))
            .setColor(White).setOpacity(1));
    toes_rInfo.addDecoration(Vec3(0),originMarker);

    const Vec3 toes_lCOM(0.023716003435794,-0.001184334594970,0.002484544347746);
    Body toes_lInfo(MassProperties(toesMass, toes_lCOM,
                            toes_Ic.shiftFromMassCenter(-toes_lCOM,toesMass)));
    toes_lInfo.addDecoration(Vec3(0.03,0.014,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(0.6,0.3,0.8)))
            .setColor(White).setOpacity(1));
    toes_lInfo.addDecoration(Vec3(0),originMarker);

    //--------------------------------------------------------------------------
    //                         Contact at feet and toes
    //--------------------------------------------------------------------------
    const Vec2 rlatmed(.04,-.04);
    const Vec2 heelball(-.02,.125);

    ContactSurface contactBall(ContactGeometry::Sphere(contactSphereRadius),
                               rubber);

    // Use this clique for contact surfaces on the humanoid that you don't want
    // to have bump into one another. They will still contact Ground.
    const ContactCliqueId clique = ContactSurface::createNewContactClique();

    contactBall.joinClique(clique);
    DecorativeSphere contactArt(contactSphereRadius);
    contactArt.setColor(Magenta);

    const Real footsphereht = -.07;
    for (int x=0; x <= 1; ++x)
        for (int z=0; z <= 1; ++z) {
            const Vec3 rctr(heelball[x], footsphereht,  rlatmed[z]);
            const Vec3 lctr(heelball[x], footsphereht, -rlatmed[z]);
            foot_rInfo.addContactSurface(rctr, contactBall);
            foot_rInfo.addDecoration(rctr, contactArt);
            foot_lInfo.addContactSurface(lctr, contactBall);
            foot_lInfo.addDecoration(lctr, contactArt);
        }

    // Balls just at toe tips.
    const Real toetip = .06;
    for (int z=0; z <= 1; ++z) {
        const Vec3 rctr(toetip, 0,  rlatmed[z]);
        const Vec3 lctr(toetip, 0, -rlatmed[z]);
        toes_rInfo.addContactSurface(rctr, contactBall);
        toes_rInfo.addDecoration(rctr, contactArt);
        toes_lInfo.addContactSurface(lctr, contactBall);
        toes_lInfo.addDecoration(lctr, contactArt);
    }

    //--------------------------------------------------------------------------
    //                         Mobilized Bodies
    //--------------------------------------------------------------------------
    // Trunk.
    // ------
    m_trunk = MobilizedBody::Free(
        m_matter.updGround(), Vec3(0),
        trunkInfo,          Vec3(0));

    // Neck angles are: q0=extension (about z), q1=bending (x), q2=rotation (y).
    m_head = MobilizedBody::Gimbal(
        m_trunk,    Transform(Rzxy, Vec3(0.010143822053248,0.222711680750785,0)),
        headInfo, Rzxy);
    addMobilityLinearStop(m_head, 0, -80, 50); // extension
    addMobilityLinearStop(m_head, 1, -60, 50); // bending
    addMobilityLinearStop(m_head, 2, -80, 80); // rotation

    // Back angles are: q0=tilt (about z), q1=list (x), q2=rotation (y).
    m_pelvis = MobilizedBody::Gimbal(
        m_trunk, Transform(Rzxy, Vec3(-0.019360589663647,-0.220484136504589,0)),
        pelvisInfo, Rzxy);
    addMobilityLinearStop(m_pelvis, 0, -5, 10); // tilt
    addMobilityLinearStop(m_pelvis, 1, -5, 5); // list
    addMobilityLinearStop(m_pelvis, 2, -15, 15); // rotation

    // Right arm.
    //-----------
    // Shoulder angles are q0=flexion, q1=adduction, q2=rotation
    m_upperarm_r = MobilizedBody::Gimbal(
        m_trunk, Transform(Rzxy,
                Vec3(-0.023921136233947,0.079313926824158,0.164710443657016)),
        upperarm_rInfo, Rzxy);
    addMobilityLinearStop(m_upperarm_r, 0, -80, 160); // flexion
    addMobilityLinearStop(m_upperarm_r, 1, -45, 45); // adduction
    addMobilityLinearStop(m_upperarm_r, 2, -20, 20); // rotation

    // Elbow angles are q0=flexion, q1=rotation
    m_lowerarm_r = MobilizedBody::Universal(
        m_upperarm_r, Transform(Rotation(-Pi/2,YAxis),
                Vec3(0.033488432642100,-0.239093933565560,0.117718445964118)),
        lowerarm_rInfo, Rotation(-Pi/2,YAxis));
    addMobilityLinearStop(m_lowerarm_r, 0, 0, 120); // flexion
    addMobilityLinearStop(m_lowerarm_r, 1, -90, 40); // rotation

    m_hand_r = MobilizedBody::Weld(
        m_lowerarm_r, Vec3(0.110610146564261,-0.232157950907188,0.014613941476432),
        hand_rInfo, Vec3(0));

    // Left arm.
    //----------
    m_upperarm_l = MobilizedBody::Gimbal(
        m_trunk, Transform(Rzmxmy,
                Vec3(-0.023921136233947,0.079313926824158,-0.164710443657016)),
        upperarm_lInfo, Rzmxmy);
    addMobilityLinearStop(m_upperarm_l, 0, -80, 160); // flexion
    addMobilityLinearStop(m_upperarm_l, 1, -45, 45); // adduction
    addMobilityLinearStop(m_upperarm_l, 2, -20, 20); // rotation

    m_lowerarm_l = MobilizedBody::Universal(
        m_upperarm_l, Transform(
                Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis),
                Vec3(0.033488432642100,-0.239093933565560,-0.117718445964118)),
        lowerarm_lInfo, Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis));
    addMobilityLinearStop(m_lowerarm_l, 0, 0, 120); // flexion
    addMobilityLinearStop(m_lowerarm_l, 1, -90, 40); // rotation

    MobilizedBody::Weld hand_l(
        m_lowerarm_l, Vec3(0.110610146564261,-0.232157950907188,-0.014613941476432),
        hand_lInfo, Vec3(0));


    // Right leg.
    //-----------
    // Hip angles are q0=flexion, q1=adduction, q2=rotation.
    m_thigh_r = MobilizedBody::Gimbal(
        m_pelvis, Transform(Rzxy,
                Vec3(0.029343644095793,-0.180413783750097,0.117592252162477)),
        thigh_rInfo, Rzxy);
    addMobilityLinearStop(m_thigh_r, 0, -60, 165); // flexion
    addMobilityLinearStop(m_thigh_r, 1, -20, 20); // adduction
    addMobilityLinearStop(m_thigh_r, 2, -120, 20); // rotation

    // Knee angle is q0=extension
    m_shank_r = MobilizedBody::Pin(
        m_thigh_r, Vec3(-0.005,-0.416780050422019,0.004172557002023),
        shank_rInfo, Vec3(0));
    addMobilityLinearStop(m_shank_r, 0, -165, 0); // extension

    // Ankle angles are q0=dorsiflexion, q1=inversion.
    m_foot_r = MobilizedBody::Universal(
        m_shank_r, Transform(
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,-0.011971751580927)),
        foot_rInfo,
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis));
    addMobilityLinearStop(m_foot_r, 0, -50, 30); // dorsiflexion
    addMobilityLinearStop(m_foot_r, 1, -2, 35); // inversion

    // Toe angle is q0=dorsiflexion
    m_toes_r = MobilizedBody::Pin(
        m_foot_r, Vec3(0.134331942132059,-0.071956467861059,-0.000000513235827),
        toes_rInfo, Vec3(0));
    addMobilityLinearStop(m_toes_r, 0, 0, 30); // dorsiflexion

    // Left leg.
    //----------
    m_thigh_l = MobilizedBody::Gimbal(
        m_pelvis, Transform(Rzmxmy,
                Vec3(0.029343644095793,-0.180413783750097,-0.117592252162477)),
        thigh_lInfo, Rzmxmy);
    addMobilityLinearStop(m_thigh_l, 0, -60, 165); // flexion
    addMobilityLinearStop(m_thigh_l, 1, -20, 20); // adduction
    addMobilityLinearStop(m_thigh_l, 2, -120, 20); // rotation

    m_shank_l = MobilizedBody::Pin(
        m_thigh_l, Vec3(-0.005,-0.416780050422019,-0.004172557002023),
        shank_lInfo, Vec3(0));
    addMobilityLinearStop(m_shank_l, 0, -165, 0); // extension

    m_foot_l = MobilizedBody::Universal(
        m_shank_l, Transform(
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,0.011971751580927)),
        foot_lInfo,
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis));
    addMobilityLinearStop(m_foot_l, 0, -50, 30); // dorsiflexion
    addMobilityLinearStop(m_foot_l, 1, -2, 35); // inversion

    m_toes_l = MobilizedBody::Pin(
        m_foot_l, Vec3(0.134331942132059,-0.071956467861059,0.000000513235827),
        toes_lInfo, Vec3(0));
    addMobilityLinearStop(m_toes_l, 0, 0, 30); // dorsiflexion
}









// TODO Make SIMBICONSTATEHANDLER not periodic, but triggered upon collisions?
// TODO make biped const within Simbicon?
// TODO hide all the visualization initializtion code.
// TODO make the body of main() really small.
// TODO summary of important pieces of code.
// TODO clean up model building
