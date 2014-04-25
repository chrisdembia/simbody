
#include <Simbody.h>


#include <cassert>
#include <vector>

using namespace std;
using namespace SimTK;

// This model should use Gimbal joints because it depends on qdot=u. But you
// can test what happens with Ball joints but picking them here.
//#define USE_BALLS
#ifdef USE_BALLS
    #define ORIENT Ball
#else
    #define ORIENT Gimbal
#endif

// #define RIGID_CONTACT

//==============================================================================
// Build the SIMBICON 3D Humanoid model
//==============================================================================

// These are the actuators in the order they were defined in the original
// model. These are indices into the controls Vector; the controls must then
// be mapped into the generalized forces Vector which is not ordered the same.
enum Actuator {
    neck_extension = 0,
    neck_bending = 1,
    neck_rotation = 2,
    back_tilt = 3,
    back_list = 4,
    back_rotation = 5,
    shoulder_r_flexion = 6,
    shoulder_r_adduction = 7,
    shoulder_r_rotation = 8,
    elbow_r_flexion = 9,
    elbow_r_rotation = 10,
    shoulder_l_flexion = 11,
    shoulder_l_adduction = 12,
    shoulder_l_rotation = 13,
    elbow_l_flexion = 14,
    elbow_l_rotation = 15,
    hip_r_adduction = 16,
    hip_r_flexion = 17,
    hip_r_rotation = 18,
    knee_r_extension = 19,
    ankle_r_inversion = 20,
    ankle_r_dorsiflexion = 21,
    mtp_r_dorsiflexion = 22,
    hip_l_adduction = 23,
    hip_l_flexion = 24,
    hip_l_rotation = 25,
    knee_l_extension = 26,
    ankle_l_inversion = 27,
    ankle_l_dorsiflexion = 28,
    mtp_l_dorsiflexion = 29,

    NumActuators = mtp_l_dorsiflexion + 1
};
class Humanoid {
public:
    Humanoid();

    // You must call this before using the model.
    void fillInActuatorMap(const SimTK::State& s);

    SimTK::QIndex getQIndex(Actuator act) const
    {   assert(!act2coord.empty()); return act2coord[act].first; }
    SimTK::UIndex getUIndex(Actuator act) const
    {   assert(!act2coord.empty()); return act2coord[act].second; }

    bool isRightFoot(const SimTK::MobilizedBody& mobod) const 
    {   return     mobod.isSameMobilizedBody(foot_r) 
                || mobod.isSameMobilizedBody(toes_r); }
    bool isLeftFoot(const SimTK::MobilizedBody& mobod) const 
    {   return     mobod.isSameMobilizedBody(foot_l) 
                || mobod.isSameMobilizedBody(toes_l); }

    // Return the number of active (force producing) contacts.
    int getNumContacts(const SimTK::State& s) const
    {   return contact.getNumContactForces(s); }

    // Set left,right to true if the corresponding foot is in contact.
    void findContactStatus(const SimTK::State& s, bool& left, bool& right) const;

    // Set fLeft, fRight to total force magnitude due to contact on that foot, 
    // zero if no contact.
    void findContactForces(const SimTK::State& s, 
                           SimTK::Real& fLeft, SimTK::Real& fRight) const;

    SimTK::MultibodySystem              system;
    SimTK::SimbodyMatterSubsystem       matter;
    SimTK::GeneralForceSubsystem        forces;
    SimTK::ContactTrackerSubsystem      tracker;
    SimTK::CompliantContactSubsystem    contact;
    SimTK::Visualizer                   viz;
    SimTK::Visualizer::InputSilo*       userInput; // just a ref; not owned here

    // Mobilized bodies.
    SimTK::MobilizedBody::Free          trunk;
    SimTK::MobilizedBody::ORIENT        head;
    SimTK::MobilizedBody::ORIENT        pelvis;
    SimTK::MobilizedBody::ORIENT        upperarm_r;
    SimTK::MobilizedBody::Universal     lowerarm_r;
    SimTK::MobilizedBody::Weld          hand_r;
    SimTK::MobilizedBody::ORIENT        upperarm_l;
    SimTK::MobilizedBody::Universal     lowerarm_l;
    SimTK::MobilizedBody::Weld          hand_l;
    SimTK::MobilizedBody::ORIENT        thigh_r;
    SimTK::MobilizedBody::Pin           shank_r;
    SimTK::MobilizedBody::Universal     foot_r;
    SimTK::MobilizedBody::Pin           toes_r;
    SimTK::MobilizedBody::ORIENT        thigh_l;
    SimTK::MobilizedBody::Pin           shank_l;
    SimTK::MobilizedBody::Universal     foot_l;
    SimTK::MobilizedBody::Pin           toes_l;
private:
    void constructSystem();
    void setUpVisualizer();

    static std::pair<SimTK::QIndex,SimTK::UIndex> getQUIx
       (const SimTK::State& s, const SimTK::MobilizedBody& mobod,int which)
    {   SimTK::QIndex q0=mobod.getFirstQIndex(s); 
        SimTK::UIndex u0=mobod.getFirstUIndex(s); 
        return std::make_pair(SimTK::QIndex(q0+which), SimTK::UIndex(u0+which)); 
    }

    // Use this to find indices to coordinate q and speed u associated with
    // each actuator.
    SimTK::Array_<std::pair<SimTK::QIndex,SimTK::UIndex>, int> act2coord;

#ifdef RIGID_CONTACT
    std::vector<SimTK::PointPlaneContact*> _leftContacts;
    std::vector<SimTK::PointPlaneContact*> _rightContacts;
#endif
};

#define ANIMATE

namespace {
const Real D2R = Pi/180; // convert degrees to radians

// Set this <1 for slow motion, >1 for speedup.
const Real RealTimeFactor 
    = 1; // try to run in real time
//  = 0.2; // 5X slower than real time

// Joint stops (coordinate limit forces).
const Real StopStiffness = 1000/D2R; // N-m/deg -> N-m/radian
const Real StopDissipation = 1;

const Real TransitionVelocity = .03; // slide->stick velocity
const Real mu_s = 5;       // Friction coefficients.
const Real mu_d = 5;
const Real mu_v = 0;

// Rubber for feet
#ifndef RIGID_CONTACT
const Real rubber_density = 1100.;  // kg/m^3
const Real rubber_young   = 0.01e9; // pascals (N/m)
const Real rubber_poisson = 0.5;    // ratio
const Real rubber_planestrain =
    ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
const Real rubber_dissipation = /*0.005*/1;

const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,
                               mu_s,mu_d,mu_v);

// Concrete for ground
const Real concrete_density = 2300.;  // kg/m^3
const Real concrete_young   = 25e9;  // pascals (N/m)
const Real concrete_poisson = 0.15;    // ratio
const Real concrete_planestrain =
    ContactMaterial::calcPlaneStrainStiffness(concrete_young,concrete_poisson);
const Real concrete_dissipation = 0.005;

const ContactMaterial concrete(concrete_planestrain,concrete_dissipation,
                               mu_s,mu_d,mu_v);

const Real ContactSphereRadius = 1.5*.01;
// Use this clique for contact surfaces on the humanoid that you don't want
// to have bump into one another. They will still contact Ground.
const ContactCliqueId ModelClique = ContactSurface::createNewContactClique();
#endif

}

Humanoid::Humanoid() 
:   system(), matter(system), forces(system), tracker(system),
    contact(system,tracker), viz(system), userInput(0)
{

    Force::Gravity(forces,matter,-YAxis,9.80660000);

#ifndef RIGID_CONTACT
    contact.setTransitionVelocity(TransitionVelocity);
    // Add the Ground contact geometry. Contact half space has -XAxis normal
    // (right hand wall) so we have to rotate.

    matter.updGround().updBody()
        .addContactSurface(Transform(Rotation(-Pi/2,ZAxis),Vec3(0)),
                    ContactSurface(ContactGeometry::HalfSpace(),concrete));
#endif

    constructSystem();
    setUpVisualizer();
}

void Humanoid::findContactStatus(const State& s, 
                                 bool& left, bool& right) const {

#ifdef RIGID_CONTACT
    left = false;
    right = false;
    for (int i = 1; i < 6; i++)
    {
        if (_leftContacts[i]) {
            left = left || _leftContacts[i]->isEnabled(s);
        }
        if (_rightContacts[i]) {
            right = right || _rightContacts[i]->isEnabled(s);
        }
    }
#else
    Real fLeft, fRight;
    findContactForces(s,fLeft,fRight);
    left = fLeft > 0; right = fRight > 0;
#endif

}

void Humanoid::findContactForces(const State& s, 
                                 Real& fLeft, Real& fRight) const {
#ifndef RIGID_CONTACT
    const int nContacts = contact.getNumContactForces(s);
    const ContactSnapshot& snapshot = tracker.getActiveContacts(s);
    //printf("  Humanoid::findContactForces(t=%g) %d contacts active\n",
    //    s.getTime(), nContacts);

    fLeft = fRight = 0;
    for (int i=0; i < nContacts; ++i) {
        const ContactForce& force = contact.getContactForce(s,i);
        const ContactId     id    = force.getContactId();
        assert(snapshot.hasContact(id));
        const Contact&      contact = snapshot.getContactById(id);
        const MobilizedBody& b1 = tracker.getMobilizedBody(contact.getSurface1());
        const MobilizedBody& b2 = tracker.getMobilizedBody(contact.getSurface2());
        const bool left = isLeftFoot(b1) || isLeftFoot(b2);
        const bool right = isRightFoot(b1) || isRightFoot(b2);
        //printf("    contact %d: bodies %d,%d (left=%d, right=%d)\n", (int)id,
        //    (int)b1.getMobilizedBodyIndex(), (int)b2.getMobilizedBodyIndex(),
        //    (int)left, (int)right);

        if (left)
            fLeft += force.getForceOnSurface2()[1].norm();
        else if (right)
            fRight += force.getForceOnSurface2()[1].norm();
    }
#endif
}


void Humanoid::constructSystem() {
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

    DecorativeSphere orgMarker(0.04*rad);
    orgMarker.setColor(Vec3(.1,.1,.1));


    //--------------------------------------------------------------------------
    //                          Body information 
    //--------------------------------------------------------------------------


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
    trunkInfo.addDecoration(Vec3(0),orgMarker);

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
    headInfo.addDecoration(Vec3(0),orgMarker);

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
    pelvisInfo.addDecoration(Vec3(0),orgMarker);

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
    upperarm_rInfo.addDecoration(Vec3(0),orgMarker);

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
    upperarm_lInfo.addDecoration(Vec3(0),orgMarker);

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
    lowerarm_rInfo.addDecoration(Vec3(0),orgMarker);

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
    lowerarm_lInfo.addDecoration(Vec3(0),orgMarker);

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
    hand_rInfo.addDecoration(Vec3(0),orgMarker);

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
    hand_lInfo.addDecoration(Vec3(0),orgMarker);

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
    thigh_rInfo.addDecoration(Vec3(0),orgMarker);

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
    thigh_lInfo.addDecoration(Vec3(0),orgMarker);


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
    shank_rInfo.addDecoration(Vec3(0),orgMarker);

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
    shank_lInfo.addDecoration(Vec3(0),orgMarker);

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
    foot_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 foot_lCOM(0.035606945567853,-0.051617802456029,0.000574057583573);
    Body foot_lInfo(MassProperties(footMass, foot_lCOM,
                      (Ifac*foot_Ic).shiftFromMassCenter(-foot_lCOM,footMass)));
    foot_lInfo.addDecoration(Vec3(0.052,-0.043,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(1.65,0.6,0.8)))
            .setColor(White).setOpacity(1));
    foot_lInfo.addDecoration(Vec3(0),orgMarker);


    const Real toesMass = 0.20;
    // This inertia looks artificially large.
    const Inertia toes_Ic(0.087132432150,0.0174264864299,0.0871324321496);
    const Vec3 toes_rCOM(0.023716003435794,-0.001184334594970,-0.002484544347746);
    Body toes_rInfo(MassProperties(toesMass, toes_rCOM,
                            toes_Ic.shiftFromMassCenter(-toes_rCOM,toesMass)));
    toes_rInfo.addDecoration(Vec3(0.03,0.014,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(0.6,0.3,0.8)))
            .setColor(White).setOpacity(1));
    toes_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 toes_lCOM(0.023716003435794,-0.001184334594970,0.002484544347746);
    Body toes_lInfo(MassProperties(toesMass, toes_lCOM,
                            toes_Ic.shiftFromMassCenter(-toes_lCOM,toesMass)));
    toes_lInfo.addDecoration(Vec3(0.03,0.014,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(0.6,0.3,0.8)))
            .setColor(White).setOpacity(1));
    toes_lInfo.addDecoration(Vec3(0),orgMarker);

    const Vec2 rlatmed(.04,-.04), heelball(-.02,.125);
#ifndef RIGID_CONTACT
    // Add contact spheres to feet and toes.
    ContactSurface contactBall(ContactGeometry::Sphere(ContactSphereRadius), 
                               rubber);
    contactBall.joinClique(ModelClique);
    DecorativeSphere contactArt(ContactSphereRadius);
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
#endif

    //--------------------------------------------------------------------------
    //                         Mobilized Bodies
    //--------------------------------------------------------------------------
    trunk = MobilizedBody::Free(matter.updGround(), Vec3(0),
                                trunkInfo,          Vec3(0));

    // Neck angles are: q0=extension (about z), q1=bending (x), q2=rotation (y).
    head = MobilizedBody::ORIENT(
        trunk,    Transform(Rzxy, Vec3(0.010143822053248,0.222711680750785,0)),
        headInfo, Rzxy);
    Force::MobilityLinearStop(forces, head, MobilizerQIndex(0), //extension
        StopStiffness,StopDissipation, -80*D2R, 50*D2R);
    Force::MobilityLinearStop(forces, head, MobilizerQIndex(1), //bending
        StopStiffness,StopDissipation, -60*D2R, 50*D2R);
    Force::MobilityLinearStop(forces, head, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -80*D2R, 80*D2R);

    // Back angles are: q0=tilt (about z), q1=list (x), q2=rotation (y).
    pelvis = MobilizedBody::ORIENT(
        trunk, Transform(Rzxy, Vec3(-0.019360589663647,-0.220484136504589,0)),
        pelvisInfo, Rzxy);
    Force::MobilityLinearStop(forces, pelvis, MobilizerQIndex(0), //tilt
        StopStiffness,StopDissipation, -5*D2R, 10*D2R);
    Force::MobilityLinearStop(forces, pelvis, MobilizerQIndex(1), //list
        StopStiffness,StopDissipation, -5*D2R, 5*D2R);
    Force::MobilityLinearStop(forces, pelvis, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -15*D2R, 15*D2R);

    // Right arm.
    //-----------
    // Shoulder angles are q0=flexion, q1=adduction, q2=rotation
    upperarm_r = MobilizedBody::ORIENT(
        trunk, Transform(Rzxy, 
                Vec3(-0.023921136233947,0.079313926824158,0.164710443657016)),
        upperarm_rInfo, Rzxy);
    Force::MobilityLinearStop(forces, upperarm_r, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, -80*D2R, 160*D2R);
    Force::MobilityLinearStop(forces, upperarm_r, MobilizerQIndex(1), //adduction
        StopStiffness,StopDissipation, -45*D2R, 45*D2R); // TODO: was -90:90
    Force::MobilityLinearStop(forces, upperarm_r, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -20*D2R, 20*D2R);

    // Elbow angles are q0=flexion, q1=rotation
    lowerarm_r = MobilizedBody::Universal(
        upperarm_r, Transform(Rotation(-Pi/2,YAxis),
                Vec3(0.033488432642100,-0.239093933565560,0.117718445964118)),
        lowerarm_rInfo, Rotation(-Pi/2,YAxis));
    Force::MobilityLinearStop(forces, lowerarm_r, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, 0*D2R, 120*D2R);
    Force::MobilityLinearStop(forces, lowerarm_r, MobilizerQIndex(1), //rotation
        StopStiffness,StopDissipation, -90*D2R, 40*D2R);

    hand_r = MobilizedBody::Weld(
        lowerarm_r, Vec3(0.110610146564261,-0.232157950907188,0.014613941476432),
        hand_rInfo, Vec3(0));

    // Left arm.
    //----------
    upperarm_l = MobilizedBody::ORIENT(
        trunk, Transform(Rzmxmy, 
                Vec3(-0.023921136233947,0.079313926824158,-0.164710443657016)),
        upperarm_lInfo, Rzmxmy);
    Force::MobilityLinearStop(forces, upperarm_l, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, -80*D2R, 160*D2R);
    Force::MobilityLinearStop(forces, upperarm_l, MobilizerQIndex(1), //adduction
        StopStiffness,StopDissipation, -45*D2R, 45*D2R); // TODO: was -90:90
    Force::MobilityLinearStop(forces, upperarm_l, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -20*D2R, 20*D2R);

    lowerarm_l = MobilizedBody::Universal(
        upperarm_l, Transform(
                Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis),
                Vec3(0.033488432642100,-0.239093933565560,-0.117718445964118)),
        lowerarm_lInfo, Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis));
    Force::MobilityLinearStop(forces, lowerarm_l, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, 0*D2R, 120*D2R);
    Force::MobilityLinearStop(forces, lowerarm_l, MobilizerQIndex(1), //rotation
        StopStiffness,StopDissipation, -90*D2R, 40*D2R);

    hand_l = MobilizedBody::Weld(
        lowerarm_l, Vec3(0.110610146564261,-0.232157950907188,-0.014613941476432),
        hand_lInfo, Vec3(0));


    // Right leg.
    //-----------
    // Hip angles are q0=flexion, q1=adduction, q2=rotation.
    thigh_r = MobilizedBody::ORIENT(
        pelvis, Transform(Rzxy, 
                Vec3(0.029343644095793,-0.180413783750097,0.117592252162477)),
        thigh_rInfo, Rzxy);
    Force::MobilityLinearStop(forces, thigh_r, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, -60*D2R, 165*D2R);
    Force::MobilityLinearStop(forces, thigh_r, MobilizerQIndex(1), //adduction
        StopStiffness,StopDissipation, -20*D2R, 20*D2R);
    Force::MobilityLinearStop(forces, thigh_r, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -120*D2R, 20*D2R);

    // Knee angle is q0=extension
    shank_r = MobilizedBody::Pin(
        thigh_r, Vec3(-0.005,-0.416780050422019,0.004172557002023),
        shank_rInfo, Vec3(0));
    Force::MobilityLinearStop(forces, shank_r, MobilizerQIndex(0), //extension
        StopStiffness,StopDissipation, -165*D2R, 0*D2R);

    // Ankle angles are q0=dorsiflexion, q1=inversion.
    foot_r = MobilizedBody::Universal(
        shank_r, Transform(
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,-0.011971751580927)),
        foot_rInfo, 
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis));
    Force::MobilityLinearStop(forces, foot_r, MobilizerQIndex(0), //dorsiflexion
        StopStiffness,StopDissipation, -50*D2R, 30*D2R);
    Force::MobilityLinearStop(forces, foot_r, MobilizerQIndex(1), //inversion
        StopStiffness,StopDissipation, -2*D2R, 35*D2R);

    // Toe angle is q0=dorsiflexion
    toes_r = MobilizedBody::Pin(
        foot_r, Vec3(0.134331942132059,-0.071956467861059,-0.000000513235827),
        toes_rInfo, Vec3(0));
    Force::MobilityLinearStop(forces, toes_r, MobilizerQIndex(0), //dorsiflexion
        StopStiffness,StopDissipation, 0*D2R, 30*D2R);

    // Left leg.
    //----------
    thigh_l = MobilizedBody::ORIENT(
        pelvis, Transform(Rzmxmy, 
                Vec3(0.029343644095793,-0.180413783750097,-0.117592252162477)),
        thigh_lInfo, Rzmxmy);
    Force::MobilityLinearStop(forces, thigh_l, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, -60*D2R, 165*D2R);
    Force::MobilityLinearStop(forces, thigh_l, MobilizerQIndex(1), //adduction
        StopStiffness,StopDissipation, -20*D2R, 20*D2R);
    Force::MobilityLinearStop(forces, thigh_l, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -120*D2R, 20*D2R);

    shank_l = MobilizedBody::Pin(
        thigh_l, Vec3(-0.005,-0.416780050422019,-0.004172557002023),
        shank_lInfo, Vec3(0));
    Force::MobilityLinearStop(forces, shank_l, MobilizerQIndex(0), //extension
        StopStiffness,StopDissipation, -165*D2R, 0*D2R);

    foot_l = MobilizedBody::Universal(
        shank_l, Transform(
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,0.011971751580927)),
        foot_lInfo, 
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis));
    Force::MobilityLinearStop(forces, foot_l, MobilizerQIndex(0), //dorsiflexion
        StopStiffness,StopDissipation, -50*D2R, 30*D2R);
    Force::MobilityLinearStop(forces, foot_l, MobilizerQIndex(1), //inversion
        StopStiffness,StopDissipation, -2*D2R, 35*D2R);

    toes_l = MobilizedBody::Pin (
        foot_l, Vec3(0.134331942132059,-0.071956467861059,0.000000513235827),
        toes_lInfo, Vec3(0));
    Force::MobilityLinearStop(forces, toes_l, MobilizerQIndex(0), //dorsiflexion
        StopStiffness,StopDissipation, 0*D2R, 30*D2R);

    //--------------------------------------------------------------------------
    //                         Rigid contact
    //--------------------------------------------------------------------------
#ifdef RIGID_CONTACT
    const Real CoefRest = 0;
    const Real mu_s = 0.8;
    const Real mu_d = 0.5;
    const Real mu_v = 0.0;
    const Real footsphereht = -.07;
    _leftContacts.resize(6);
    _rightContacts.resize(6);
    int i = 0;
    for (int x=0; x<=1; ++x) {
        for (int z=0; z<=1; ++z) {
            const Vec3 rctr(heelball[x], footsphereht,  rlatmed[z]);
            const Vec3 lctr(heelball[x], footsphereht, -rlatmed[z]);
            PointPlaneContact * ppcr = new PointPlaneContact(
                    matter.updGround(), YAxis, 0.,
                    foot_r, rctr, CoefRest, mu_s, mu_d, mu_v);
            PointPlaneContact * ppcl = new PointPlaneContact(
                    matter.updGround(), YAxis, 0.,
                    foot_l, lctr, CoefRest, mu_s, mu_d, mu_v);
            matter.adoptUnilateralContact(ppcr);
            matter.adoptUnilateralContact(ppcl);

            _rightContacts[i] = ppcr;
            _leftContacts[i] = ppcl;

            i += 1;
        }
    }
    i = 0;

    const Real toetip = .06;
    for (int z=0; z <= 1; ++z) {
        const Vec3 rctr(toetip, 0,  rlatmed[z]);
        const Vec3 lctr(toetip, 0, -rlatmed[z]);
        PointPlaneContact * ppcr = new PointPlaneContact(
                matter.updGround(), YAxis, 0.,
                toes_r, rctr, CoefRest, mu_s, mu_d, mu_v);
        PointPlaneContact * ppcl = new PointPlaneContact(
                matter.updGround(), YAxis, 0.,
                toes_l, lctr, CoefRest, mu_s, mu_d, mu_v);
        matter.adoptUnilateralContact(ppcr);
        matter.adoptUnilateralContact(ppcl);

        _rightContacts[i] = ppcr;
        _leftContacts[i] = ppcl;

        i += 1;
    }
#endif


}

void Humanoid::setUpVisualizer() {
    // Set up visualization and ask for a frame every 1/30 second.
    viz.setShowSimTime(true); viz.setShowFrameRate(true);
    userInput = new Visualizer::InputSilo();
    viz.addInputListener(userInput);   
    #ifdef ANIMATE
    system.addEventReporter(new Visualizer::Reporter(viz, RealTimeFactor/30));
    #endif
    DecorativeText help("Any input to start; ESC to quit");
    help.setIsScreenText(true);
    viz.addDecoration(MobilizedBodyIndex(0),Vec3(0),help);
    matter.setShowDefaultGeometry(false);
}

void Humanoid::fillInActuatorMap(const State& s) {
    act2coord.resize(NumActuators); 

    act2coord[neck_extension]= getQUIx(s,head, 0); 
    act2coord[neck_bending]  = getQUIx(s,head, 1); 
    act2coord[neck_rotation] = getQUIx(s,head, 2); 

    act2coord[back_tilt]     = getQUIx(s,pelvis, 0); 
    act2coord[back_list]     = getQUIx(s,pelvis, 1); 
    act2coord[back_rotation] = getQUIx(s,pelvis, 2); 

    act2coord[shoulder_r_flexion]   = getQUIx(s,upperarm_r, 0); 
    act2coord[shoulder_r_adduction] = getQUIx(s,upperarm_r, 1); 
    act2coord[shoulder_r_rotation]  = getQUIx(s,upperarm_r, 2);

    act2coord[elbow_r_flexion]  = getQUIx(s,lowerarm_r, 0); 
    act2coord[elbow_r_rotation] = getQUIx(s,lowerarm_r, 1); 

    act2coord[shoulder_l_flexion]   = getQUIx(s,upperarm_l, 0); 
    act2coord[shoulder_l_adduction] = getQUIx(s,upperarm_l, 1); 
    act2coord[shoulder_l_rotation]  = getQUIx(s,upperarm_l, 2);

    act2coord[elbow_l_flexion]  = getQUIx(s,lowerarm_l, 0); 
    act2coord[elbow_l_rotation] = getQUIx(s,lowerarm_l, 1); 

    act2coord[hip_r_flexion]    = getQUIx(s,thigh_r, 0); 
    act2coord[hip_r_adduction]  = getQUIx(s,thigh_r, 1); 
    act2coord[hip_r_rotation]   = getQUIx(s,thigh_r, 2); 

    act2coord[knee_r_extension] = getQUIx(s,shank_r, 0); 

    act2coord[ankle_r_dorsiflexion] = getQUIx(s,foot_r, 0); 
    act2coord[ankle_r_inversion]    = getQUIx(s,foot_r, 1); 

    act2coord[mtp_r_dorsiflexion] = getQUIx(s,toes_r, 0); 

    act2coord[hip_l_flexion]    = getQUIx(s,thigh_l, 0); 
    act2coord[hip_l_adduction]  = getQUIx(s,thigh_l, 1); 
    act2coord[hip_l_rotation]   = getQUIx(s,thigh_l, 2); 

    act2coord[knee_l_extension] = getQUIx(s,shank_l, 0); 

    act2coord[ankle_l_dorsiflexion] = getQUIx(s,foot_l, 0); 
    act2coord[ankle_l_inversion]    = getQUIx(s,foot_l, 1); 

    act2coord[mtp_l_dorsiflexion] = getQUIx(s,toes_l, 0); 
}


// SIMBICON

#define STATE_UPD_STEPSIZE 0.0005
//#define RIGID_CONTACT
 
class Simbicon : public SimTK::Force::Custom::Implementation
{
public:
	Simbicon(Humanoid& model); 
    void calcForce(const SimTK::State&                state, 
                   SimTK::Vector_<SimTK::SpatialVec>& bodyForces, 
                   SimTK::Vector_<SimTK::Vec3>&       particleForces, 
                   SimTK::Vector&                     mobilityForces) const 
                   OVERRIDE_11;

    SimTK::Real calcPotentialEnergy(const SimTK::State& state) const 
        OVERRIDE_11 {
        return 0;
    }

	int getState() const {
		return _state; 
	}
    void setState( int state, double time ) {
		_state = state; 
		_stateStartTime = time;
		if (state == 0 || state == 2) {
			// reset the swing thigh orientation when the swing leg changes
			for (int i = 0; i < 2; i++) {
				_lastSWTAngle[i] = -100.0;  
				_curSWTAngle[i] = -100.0;  
			}
        }
	}

	double getStateStartTime() const {
		return _stateStartTime; 
	}

	void computeSecondaryStateVals(const SimTK::State& s,
        SimTK::Real lForce, SimTK::Real rForce); 
	
private:
    void computeControls(const SimTK::State& s, SimTK::Vector& controls) const;

	void getSagCorNormals( const SimTK::State& s, 
		SimTK::Vec3& sagN, SimTK::Vec3& corN ) const; 
	void getUpVectorInGround( const SimTK::State& s, const SimTK::MobilizedBody& b, 
		SimTK::Vec3& up ) const; 
	void getFrontVectorInGround( const SimTK::State& s, const SimTK::MobilizedBody& b, 
		SimTK::Vec3& up ) const; 
	void fillInHipJointControls( const SimTK::State& s, 
		SimTK::Vector& controls ) const; 
	
    Humanoid& _model;
	int    _state;
	double _stateStartTime; 
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

class SimbiconStateHandler : public SimTK::PeriodicEventHandler {
public:
    SimbiconStateHandler(Humanoid& m, Simbicon& simctrl, SimTK::Real interval);

    void handleEvent(SimTK::State& s, SimTK::Real accuracy, bool& shouldTerminate) const;
private:
    Humanoid& _model;
    Simbicon& _simctrl;
};

// Normally Simbicon has 4 states per gait cycle (i.e. 2 per leg), 2 state is
// simplified and not as realistic.
#define TWO_STATE
// Torque angle to be flat relative to ground (not in original Simbicon paper)
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

Simbicon::Simbicon(Humanoid& model) : _model(model) {
	_state = -1;
	_stateStartTime = -1.0;
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

void Simbicon::calcForce(const State&         s, 
                         Vector_<SpatialVec>& /*bodyForces*/, 
                         Vector_<Vec3>&       /*particleForces*/, 
                         Vector&              mobilityForces) const 
{
    Vector controls(NumActuators);
    computeControls(s,controls);

    for (int i=0; i<NumActuators; ++i)
        mobilityForces[_model.getUIndex(Actuator(i))] = controls[i];
}


// Project the pelvis z (right) and x (forward) directions onto the Ground
// (x-z) plane, and normalize. Projected z is the normal to the sagittal plane;
// projected x is the normal to the coronal plane.
void Simbicon::getSagCorNormals(const State& s, Vec3& sagN, Vec3& corN ) const {
	const MobilizedBody& pelvis = _model.pelvis;
	
    sagN = Vec3(pelvis.getBodyRotation(s).z());
    corN = Vec3(pelvis.getBodyRotation(s).x());
	sagN[YAxis] = 0; // project to y==0 
	sagN = sagN.normalize();  
	corN[YAxis] = 0; // project to y==0
	corN = corN.normalize();  
}

void Simbicon::getUpVectorInGround( const State& s, const MobilizedBody& b, 
	Vec3& up ) const { 
	up = Vec3(b.getBodyRotation(s).y());   
}

void Simbicon::getFrontVectorInGround( const State& s, const MobilizedBody& b, 
	Vec3& front ) const { 
	front = Vec3(b.getBodyRotation(s).x());
}

void Simbicon::fillInHipJointControls( const State& s, Vector& controls ) const {
	const SimbodyMatterSubsystem& matter = _model.matter;
	Vec3 sagN, corN; 	
	
    // Sagittal and coronal plane normals are right and front directions of
    // the pelvis, projected on the ground plane.
	getSagCorNormals(s, sagN, corN); 	
	
	int swh = hip_r_flexion;             // swing hip
	int sth = hip_l_flexion;             // stance hip
    MobilizedBody ankle = _model.foot_l; // stance ankle
	
	if (_state == 2 || _state == 3) { // right stance
		swh = hip_l_flexion; 
		sth = hip_r_flexion; 
        ankle = _model.foot_r;
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
	if (_state == 1 || _state == 3) { 
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
	if (_state == 0 || _state == 1) // left stance	
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

void Simbicon::computeControls(const State& s, Vector& controls) const
{
#ifdef DROP_LANDING
	for (int i = 0; i < controls.size(); i++) 
			controls[i] = 0.0; 
	return; 
#else
	int swh = hip_r_flexion; 
	int sth = hip_l_flexion; 
	
	if (_state == 2 || _state == 3) { // right stance
		swh = hip_l_flexion; 
		sth = hip_r_flexion; 
	}
	int swk = swh + 2; // swing knee
	int swa = swh + 4; // swing ankle
	int stk = sth + 2; // stance knee
	int sta = sth + 4; // stance ankle 
	for (int i = 0; i < NumActuators; i++) {
        double kp, kd;          // position gain, derivative gain
        calcGainsFromStrength(DefaultStrength, kp, kd);
		double thetad = 0.0;                    // desired angle
		
		if (i == neck_extension || i == neck_bending || i == neck_rotation) {
            calcGainsFromStrength(NeckStrength, kp, kd);
		}
		else if (i == back_tilt || i == back_list || i == back_rotation) {
            calcGainsFromStrength(BackStrength, kp, kd);
		}
		else if (   i == shoulder_r_flexion   || i == shoulder_l_flexion 
                 || i == shoulder_r_adduction || i == shoulder_l_adduction
                 || i == elbow_r_flexion      || i == elbow_l_flexion) {
            calcGainsFromStrength(ArmFlexionAdductionStrength, kp, kd);
		}
		else if (   i == shoulder_r_rotation || i == shoulder_l_rotation 
                 || i == elbow_r_rotation    || i == elbow_l_rotation) {
            calcGainsFromStrength(ArmRotationStrength, kp, kd);
		}
		else if (i == hip_r_rotation || i == hip_l_rotation) {
            calcGainsFromStrength(HipRotationStrength, kp, kd);
		}
		else if (i == knee_r_extension || i == knee_l_extension) {
            calcGainsFromStrength(KneeStrength, kp, kd);
		}
        else if (i == ankle_r_dorsiflexion || i == ankle_l_dorsiflexion) {
            calcGainsFromStrength(AnkleFlexionStrength, kp, kd);
        }
        else if (i == ankle_r_inversion || i == ankle_l_inversion) {
            calcGainsFromStrength(AnkleInversionStrength, kp, kd);
        }
		else if (i == mtp_r_dorsiflexion || i == mtp_l_dorsiflexion) { // toe
            calcGainsFromStrength(ToeStrength, kp, kd);
		}

		if (_state >= 0) {
			if (i == swk) 
				thetad = -1.1;  
			else if (i == swa) 
				thetad = 0.6;  
			else if (i == stk) 
				thetad = -0.05;

            #ifndef TWO_STATE	
			if (_state == 1 || _state == 3) { 
				if (i == swk) 
					thetad = -0.05;  
				else if (i == swa) 
					thetad = 0.15;  
				else if (i == stk) 
					thetad = -0.1; 
			}
            #endif
		}

        const Real qi = s.getQ()[_model.getQIndex(Actuator(i))];
        const Real ui = s.getU()[_model.getUIndex(Actuator(i))];
		controls[i] = clamp(kp*(thetad - qi) - kd*ui, kp); 
	}
	
	if (_state >= 0) 
		fillInHipJointControls(s, controls); 

	Vec3 com = _model.matter.calcSystemMassCenterLocationInGround(s); 
	if (com[1] < 0.7) {
		for (int i = 0; i < controls.size(); i++) 
			controls[i] = 0.0; 
	}
	return;
#endif
}
	
void Simbicon::
computeSecondaryStateVals(const State& s, Real lForce, Real rForce) {
	const MobilizedBody& pelvis = _model.pelvis;

#ifndef RIGID_CONTACT
    _lastRFootContactForce = _curRFootContactForce;
    _lastLFootContactForce = _curLFootContactForce;
    _curRFootContactForce = rForce;
    _curLFootContactForce = lForce;
#endif

	Vec3 upThigh;
	if (_state == 0 || _state == 1)  
		getUpVectorInGround(s, _model.thigh_r, upThigh); 
	else if (_state == 2 || _state == 3) 
		getUpVectorInGround(s, _model.thigh_l, upThigh); 
	
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
	const MobilizedBody& rFoot = _model.foot_r;
	const MobilizedBody& lFoot = _model.foot_l;
	
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

SimbiconStateHandler::SimbiconStateHandler(Humanoid& model, Simbicon& simctrl,
    Real interval) 
        : PeriodicEventHandler(interval), _model(model), _simctrl(simctrl) {
	}

void SimbiconStateHandler::handleEvent(State& s, Real accuracy, bool& shouldTerminate) const
{
    shouldTerminate = false;
	Simbicon* simctrl = &_simctrl;

    _model.system.realize(s, Stage::Dynamics);

    Real lForce, rForce;
    _model.findContactForces(s, lForce, rForce);
    const bool lContact = (lForce > 0);
    const bool rContact = (rForce > 0);


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

    simctrl->computeSecondaryStateVals(s, lForce, rForce);

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
    OutputReporter(const Humanoid& model, 
                   Real            interval); 
    void handleEvent(const State& state) const OVERRIDE_11;
private:
    const Humanoid& _model;
};

// Write interesting integrator info to stdout at end of simulation.
void dumpIntegratorStats(double startCPU, double startTime,
                         const Integrator& integ);
}

class ShowContact : public DecorationGenerator {
public:
    ShowContact(const Humanoid& unis) 
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
    const Humanoid& m_unis;
};

//==============================================================================
// MAIN FUNCTION
//==============================================================================
int main(int argc, char **argv)
{
    try {

	double finalTime = 1000; 

    Humanoid model;
    model.system.addEventHandler
       (new UserInputHandler(*model.userInput, Real(0.1))); //check input every 100ms

    model.system.addEventReporter(new OutputReporter(model, .01));

    // Add the controller.
    Simbicon* simctrl = new Simbicon(model);
    Force::Custom simbicon(model.forces, simctrl); // takes ownership

#ifndef RIGID_CONTACT
    model.system.addEventHandler(new SimbiconStateHandler(model,*simctrl,
                                                          STATE_UPD_STEPSIZE));
#endif

    State s = model.system.realizeTopology();
    #ifdef USE_BALLS
        model.matter.setUseEulerAngles(s, true);
    #endif
    model.system.realizeModel(s);
    model.fillInActuatorMap(s);

    printf("Act: u\n");
    for (int i=0; i < NumActuators; ++i) {
        printf("%2d: %d\n", i, int(model.getUIndex(Actuator(i))));
    }

    //model.toes_r.lockAt(s, .2); model.toes_l.lockAt(s, .2);
    //model.foot_r.lockAt(s, Vec2(.2,0)); model.foot_l.lockAt(s, Vec2(.2,0));
    model.system.realize(s, Stage::Instance);

    printf("SIMBICON 3D:\n");
    printf("%d bodies, %d mobilities, -%d constraint equations -%d motions\n",
        model.matter.getNumBodies(), s.getNU(), s.getNMultipliers(), 
        model.matter.getKnownUDotIndex(s).size());

    model.trunk.setQToFitTranslation(s, Vec3(0,1.5,0));
    model.trunk.setUToFitLinearVelocity(s, Vec3(1,0,0));
    model.viz.report(s); 

#ifdef RIGID_CONTACT
    model.viz.addDecorationGenerator(new ShowContact(model));  
#endif

#ifndef RIGID_CONTACT
    // Simulate.
    //CPodesIntegrator integ(model.system); integ.setOrderLimit(2); integ.setAccuracy(.01);
    
    //RungeKuttaMersonIntegrator integ(model.system); integ.setAccuracy(1e-3);
    //RungeKutta2Integrator integ(model.system); integ.setAccuracy(.1);
    SemiExplicitEuler2Integrator integ(model.system); integ.setAccuracy(0.1);
    //SemiImplicitEulerIntegrator integ(model.system, .002);
    //integ.setConstraintTolerance(.001);
    //integ.setMaximumStepSize(.005);
    TimeStepper ts(model.system, integ);
#else
    SemiExplicitEulerTimeStepper ts(model.system);
    ts.setImpulseSolverType(SemiExplicitEulerTimeStepper::PLUS);
    ts.setInducedImpactModel(SemiExplicitEulerTimeStepper::Simultaneous);
    //ts.setConstraintTol(1e-5);
    ts.setAccuracy(1e-4);
    ts.setMaxInducedImpactsPerStep(1000);
    //ts.setRestitutionModel(SemiExplicitEulerTimeStepper::Poisson);
#endif
    ts.initialize(s);
    model.viz.report(ts.getState());
    printf("Hit ENTER to simulate ... (ESC to quit)\n");
    model.userInput->waitForAnyUserInput(); model.userInput->clear();

    const double startCPU  = cpuTime(), startTime = realTime();

    try {

#ifdef RIGID_CONTACT
    int timeStep = 0;
    do {
        if (timeStep % 10 == 0) {
            model.viz.report(ts.getState());
        }
        ts.stepTo(ts.getState().getTime() + STATE_UPD_STEPSIZE);

        State& s = const_cast<State&>(ts.getState());

        model.system.realize(s, Stage::Dynamics);

        bool lContact, rContact;
        model.findContactStatus(s, lContact, rContact);


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
        model.system.realize(ts.getState());
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
OutputReporter::OutputReporter(const Humanoid& model, 
                               Real            interval)
:   PeriodicEventReporter(interval), _model(model) {}

void OutputReporter::handleEvent(const State& state) const {
    //printf("OutputReporter @t=%g:\n", state.getTime());
    //Real fLeft, fRight;
    //_model.findContactForces(state, fLeft, fRight);
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
