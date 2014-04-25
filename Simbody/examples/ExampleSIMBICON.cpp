
#include <Simbody.h>


#include <cassert>
#include <vector>

using namespace std;
using namespace SimTK;

// #define RIGID_CONTACT

const Real D2R = Pi / 180.0;

// TODO remove
const Real RealTimeFactor 
    = 1; // try to run in real time

// TODO remove
const int NumActuators = 30;

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

    /// Realizes the topology and model, and uses the resulting initial state
    /// to perform further internal initialization.
    void initialize(State& state) {
        state = realizeTopology();
        realizeModel(state);

        // Initialization specific to Biped.
        fillInCoordinateMap(state);
    }

    /// Populates a map that allow users of Biped to access coordinates (Q) and
    /// speeds (U) with descriptive Enums (e.g., neck_extension) so the user
    /// (e.g., SIMBICON) needn't know how State indices correspond to
    /// coordinates.
    void fillInCoordinateMap(const State& s);

    void setTrunkOriginPosition(State& s, Vec3 posInGround) const {
        m_bodies.at(trunk).setQToFitTranslation(s, posInGround);
    }
    void setTrunkOriginVelocity(State& s, Vec3 velocityInGround) const {
        m_bodies.at(trunk).setUToFitLinearVelocity(s, velocityInGround);
    }

    /// param[out] fLeft vertical contact force on left foot.
    /// param[out] fRight vertical contact force on right foot.
    void findContactForces(const State& s, Real& fLeft, Real& fRight) const;

    /// param[out] left left foot is in contact with the ground.
    /// param[out] right right foot is in contact with the ground.
    void findContactStatus(const State& s, bool& left, bool& right) const;

    /// The names of the rigid bodies that make up the biped.
    /// Used for indexing into data structures.
    enum Segment {
        trunk, // called 'torso' in paper.
        head,
        pelvis,
        upperarm_r, lowerarm_r, hand_r, upperarm_l, lowerarm_l, hand_l,
        thigh_r, shank_r, foot_r, toes_r,
        thigh_l, shank_l, foot_l, toes_l
    };

    /// The names of the degrees of freedom of the biped.
    /// Used for indexing into data structures.
    enum Coordinate {
        neck_extension, neck_bending, neck_rotation,
        back_tilt, back_list, back_rotation,
        shoulder_r_flexion, shoulder_r_adduction, shoulder_r_rotation,
        elbow_r_flexion, elbow_r_rotation,
        shoulder_l_flexion, shoulder_l_adduction, shoulder_l_rotation,
        elbow_l_flexion, elbow_l_rotation,
        hip_r_adduction, hip_r_flexion, hip_r_rotation,
        knee_r_extension,
        ankle_r_inversion, ankle_r_dorsiflexion,
        mtp_r_dorsiflexion,
        hip_l_adduction, hip_l_flexion, hip_l_rotation,
        knee_l_extension,
        ankle_l_inversion, ankle_l_dorsiflexion,
        mtp_l_dorsiflexion,
    };

    /// Index, in the State vector, of this coordinate's value (Q).
    QIndex getQIndex(Coordinate coord) const {
        assert(!m_coordinates.empty());
        return m_coordinates.at(coord).first;
    }

    /// Coordinate coord's value (Q).
    Real getQ(const State& s, Coordinate coord) const {
        return s.getQ()[getQIndex(coord)];
    }

    /// Index, in the State vector, of this coordinate's speed (U).
    UIndex getUIndex(Coordinate coord) const {
        assert(!m_coordinates.empty());
        return m_coordinates.at(coord).second;
    }

    /// Coordinate coord's speed (U).
    Real getU(const State& s, Coordinate coord) const {
        return s.getU()[getUIndex(coord)];
    }

    /// Add in a generalized force into the proper place in `mobForces`,
    /// as determined by `coord`.
    void addInForce(const Coordinate coord, const Real force,
            Vector& mobForces) const {
        Real forceAccountingForSymmetry = force;
        // TODO
        //if (coord == hip_r_adduction) {
        //    // Negate so that we actually apply an adduction torque (not an
        //    // abduction torque).
        //    forceAccountingForSymmetry *= -1;
        //}
        mobForces[getUIndex(coord)] += forceAccountingForSymmetry;
    }

    /// The pelvis' z direction (which points to the biped's right) projected
    /// into the plane of the ground, and normalized into a unit vector.
    UnitVec3 getNormalToSagittalPlane(const State& s) const
    {
        Vec3 bodyZInGround = Vec3(m_bodies.at(pelvis).getBodyRotation(s).z());
        // Project onto a horizontal plane.
        bodyZInGround[YAxis] = 0;
        // Make into a unit vector.
        return UnitVec3(bodyZInGround);
    }

    /// The pelvis' x direction (which points in the direction of walking)
    /// projected into the plane of the ground, and normalized into a
    /// unit vector.
    UnitVec3 getNormalToCoronalPlane(const State& s) const
    {
        Vec3 bodyXInGround = Vec3(m_bodies.at(pelvis).getBodyRotation(s).x());
        // Project onto a horizontal plane.
        bodyXInGround[YAxis] = 0;
        // Make into a unit vector.
        return UnitVec3(bodyXInGround);
    }

    const MobilizedBody& getBody(Segment seg) const {
        return m_bodies.at(seg);
    }

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

    /// The indices, in the State vector, of a coordinate value (Q) and speed
    /// (U) in a specific MobilizedBody. Used in fillInCoordinateMap().
    static std::pair<QIndex, UIndex> getQandUIndices(
            const SimTK::State& s, const SimTK::MobilizedBody& mobod, int which)
    {
        SimTK::QIndex q0 = mobod.getFirstQIndex(s);
        SimTK::UIndex u0 = mobod.getFirstUIndex(s);
        return std::make_pair(QIndex(q0 + which), UIndex(u0 + which));
    }

    bool isRightFoot(const MobilizedBody& mobod) const {
        return     mobod.isSameMobilizedBody(m_bodies.at(foot_r))
                || mobod.isSameMobilizedBody(m_bodies.at(toes_r));
    }
    bool isLeftFoot(const MobilizedBody& mobod) const {
        return     mobod.isSameMobilizedBody(m_bodies.at(foot_l))
                || mobod.isSameMobilizedBody(m_bodies.at(toes_l));
    }

    // Subsystems.
    SimbodyMatterSubsystem       m_matter;
    GeneralForceSubsystem        m_forces;
    ContactTrackerSubsystem      m_tracker;
    CompliantContactSubsystem    m_contact;

    /// The MobilizedBody's that make up the Biped.
    std::map<Segment, MobilizedBody> m_bodies;

    /// The indices of each coordinate's Q and U in the State vector.
    std::map<Coordinate, std::pair<QIndex, UIndex> > m_coordinates;
};

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
    const Real stopStiffness = 1000 / D2R; // N-m/deg -> N-m/radian
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

    // Contact material: concrete.
    const Real concrete_density = 2300.;  // kg/m^3
    const Real concrete_young   = 25e9;  // pascals (N/m)
    const Real concrete_poisson = 0.15;    // ratio
    const Real concrete_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(concrete_young,concrete_poisson);
    const Real concrete_dissipation = 0.005;

    const ContactMaterial concrete(concrete_planestrain,concrete_dissipation,
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
    //                          Gravity and ground contact
    //--------------------------------------------------------------------------

    // Gravity.
    Force::Gravity(m_forces, m_matter, -YAxis, 9.8066);

    // Contact with the ground.
    m_matter.updGround().updBody().addContactSurface(
            Transform(Rotation(-0.5 * Pi, ZAxis), Vec3(0)),
            ContactSurface(ContactGeometry::HalfSpace(), concrete));


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
    m_bodies[trunk] = MobilizedBody::Free(
        m_matter.updGround(), Vec3(0),
        trunkInfo,          Vec3(0));

    // Neck angles are: q0=extension (about z), q1=bending (x), q2=rotation (y).
    m_bodies[head] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzxy, Vec3(0.010143822053248,0.222711680750785,0)),
        headInfo, Rzxy);
    addMobilityLinearStop(m_bodies[head], 0, -80, 50); // extension
    addMobilityLinearStop(m_bodies[head], 1, -60, 50); // bending
    addMobilityLinearStop(m_bodies[head], 2, -80, 80); // rotation

    // Back angles are: q0=tilt (about z), q1=list (x), q2=rotation (y).
    m_bodies[pelvis] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzxy, Vec3(-0.019360589663647,-0.220484136504589,0)),
        pelvisInfo, Rzxy);
    addMobilityLinearStop(m_bodies[pelvis], 0, -5, 10); // tilt
    addMobilityLinearStop(m_bodies[pelvis], 1, -5, 5); // list
    addMobilityLinearStop(m_bodies[pelvis], 2, -15, 15); // rotation

    // Right arm.
    //-----------
    // Shoulder angles are q0=flexion, q1=adduction, q2=rotation
    m_bodies[upperarm_r] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzxy,
                Vec3(-0.023921136233947,0.079313926824158,0.164710443657016)),
        upperarm_rInfo, Rzxy);
    addMobilityLinearStop(m_bodies[upperarm_r], 0, -80, 160); // flexion
    addMobilityLinearStop(m_bodies[upperarm_r], 1, -45, 45); // adduction
    addMobilityLinearStop(m_bodies[upperarm_r], 2, -20, 20); // rotation

    // Elbow angles are q0=flexion, q1=rotation
    m_bodies[lowerarm_r] = MobilizedBody::Universal(
        m_bodies[upperarm_r], Transform(Rotation(-Pi/2,YAxis),
                Vec3(0.033488432642100,-0.239093933565560,0.117718445964118)),
        lowerarm_rInfo, Rotation(-Pi/2,YAxis));
    addMobilityLinearStop(m_bodies[lowerarm_r], 0, 0, 120); // flexion
    addMobilityLinearStop(m_bodies[lowerarm_r], 1, -90, 40); // rotation

    m_bodies[hand_r] = MobilizedBody::Weld(
        m_bodies[lowerarm_r], Vec3(0.110610146564261,-0.232157950907188,0.014613941476432),
        hand_rInfo, Vec3(0));

    // Left arm.
    //----------
    m_bodies[upperarm_l] = MobilizedBody::Gimbal(
        m_bodies[trunk], Transform(Rzmxmy,
                Vec3(-0.023921136233947,0.079313926824158,-0.164710443657016)),
        upperarm_lInfo, Rzmxmy);
    addMobilityLinearStop(m_bodies[upperarm_l], 0, -80, 160); // flexion
    addMobilityLinearStop(m_bodies[upperarm_l], 1, -45, 45); // adduction
    addMobilityLinearStop(m_bodies[upperarm_l], 2, -20, 20); // rotation

    m_bodies[lowerarm_l] = MobilizedBody::Universal(
        m_bodies[upperarm_l], Transform(
                Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis),
                Vec3(0.033488432642100,-0.239093933565560,-0.117718445964118)),
        lowerarm_lInfo, Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis));
    addMobilityLinearStop(m_bodies[lowerarm_l], 0, 0, 120); // flexion
    addMobilityLinearStop(m_bodies[lowerarm_l], 1, -90, 40); // rotation

    m_bodies[hand_l] = MobilizedBody::Weld(
        m_bodies[lowerarm_l], Vec3(0.110610146564261,-0.232157950907188,-0.014613941476432),
        hand_lInfo, Vec3(0));


    // Right leg.
    //-----------
    // Hip angles are q0=flexion, q1=adduction, q2=rotation.
    m_bodies[thigh_r] = MobilizedBody::Gimbal(
        m_bodies[pelvis], Transform(Rzxy,
                Vec3(0.029343644095793,-0.180413783750097,0.117592252162477)),
        thigh_rInfo, Rzxy);
    addMobilityLinearStop(m_bodies[thigh_r], 0, -60, 165); // flexion
    addMobilityLinearStop(m_bodies[thigh_r], 1, -20, 20); // adduction
    addMobilityLinearStop(m_bodies[thigh_r], 2, -120, 20); // rotation

    // Knee angle is q0=extension
    m_bodies[shank_r] = MobilizedBody::Pin(
        m_bodies[thigh_r], Vec3(-0.005,-0.416780050422019,0.004172557002023),
        shank_rInfo, Vec3(0));
    addMobilityLinearStop(m_bodies[shank_r], 0, -165, 0); // extension

    // Ankle angles are q0=dorsiflexion, q1=inversion.
    m_bodies[foot_r] = MobilizedBody::Universal(
        m_bodies[shank_r], Transform(
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,-0.011971751580927)),
        foot_rInfo,
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis));
    addMobilityLinearStop(m_bodies[foot_r], 0, -50, 30); // dorsiflexion
    addMobilityLinearStop(m_bodies[foot_r], 1, -2, 35); // inversion

    // Toe angle is q0=dorsiflexion
    m_bodies[toes_r] = MobilizedBody::Pin(
        m_bodies[foot_r], Vec3(0.134331942132059,-0.071956467861059,-0.000000513235827),
        toes_rInfo, Vec3(0));
    addMobilityLinearStop(m_bodies[toes_r], 0, 0, 30); // dorsiflexion

    // Left leg.
    //----------
    m_bodies[thigh_l] = MobilizedBody::Gimbal(
        m_bodies[pelvis], Transform(Rzmxmy,
                Vec3(0.029343644095793,-0.180413783750097,-0.117592252162477)),
        thigh_lInfo, Rzmxmy);
    addMobilityLinearStop(m_bodies[thigh_l], 0, -60, 165); // flexion
    addMobilityLinearStop(m_bodies[thigh_l], 1, -20, 20); // adduction
    addMobilityLinearStop(m_bodies[thigh_l], 2, -120, 20); // rotation

    m_bodies[shank_l] = MobilizedBody::Pin(
        m_bodies[thigh_l], Vec3(-0.005,-0.416780050422019,-0.004172557002023),
        shank_lInfo, Vec3(0));
    addMobilityLinearStop(m_bodies[shank_l], 0, -165, 0); // extension

    m_bodies[foot_l] = MobilizedBody::Universal(
        m_bodies[shank_l], Transform(
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,0.011971751580927)),
        foot_lInfo,
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis));
    addMobilityLinearStop(m_bodies[foot_l], 0, -50, 30); // dorsiflexion
    addMobilityLinearStop(m_bodies[foot_l], 1, -2, 35); // inversion

    m_bodies[toes_l] = MobilizedBody::Pin(
        m_bodies[foot_l], Vec3(0.134331942132059,-0.071956467861059,0.000000513235827),
        toes_lInfo, Vec3(0));
    addMobilityLinearStop(m_bodies[toes_l], 0, 0, 30); // dorsiflexion

}

void Biped::fillInCoordinateMap(const State& s)
{
    // So all lines below fit within 80 columns:
    std::map<Coordinate, std::pair<QIndex, UIndex> >& coords = m_coordinates;

    // Upper body.
    coords[neck_extension] = getQandUIndices(s, m_bodies[head], 0);
    coords[neck_bending] = getQandUIndices(s, m_bodies[head], 1);
    coords[neck_rotation] = getQandUIndices(s, m_bodies[head], 2);

    coords[back_tilt]     = getQandUIndices(s, m_bodies[pelvis], 0);
    coords[back_list]     = getQandUIndices(s, m_bodies[pelvis], 1);
    coords[back_rotation] = getQandUIndices(s, m_bodies[pelvis], 2);

    // right.
    coords[shoulder_r_flexion]   = getQandUIndices(s, m_bodies[upperarm_r], 0);
    coords[shoulder_r_adduction] = getQandUIndices(s, m_bodies[upperarm_r], 1);
    coords[shoulder_r_rotation]  = getQandUIndices(s, m_bodies[upperarm_r], 2);

    coords[elbow_r_flexion]  = getQandUIndices(s, m_bodies[lowerarm_r], 0);
    coords[elbow_r_rotation] = getQandUIndices(s, m_bodies[lowerarm_r], 1);

    // left.
    coords[shoulder_l_flexion]   = getQandUIndices(s, m_bodies[upperarm_l], 0);
    coords[shoulder_l_adduction] = getQandUIndices(s, m_bodies[upperarm_l], 1);
    coords[shoulder_l_rotation]  = getQandUIndices(s, m_bodies[upperarm_l], 2);

    coords[elbow_l_flexion]  = getQandUIndices(s, m_bodies[lowerarm_l], 0);
    coords[elbow_l_rotation] = getQandUIndices(s, m_bodies[lowerarm_l], 1);

    // Lower body.
    // right.
    coords[hip_r_flexion]    = getQandUIndices(s, m_bodies[thigh_r], 0);
    coords[hip_r_adduction]  = getQandUIndices(s, m_bodies[thigh_r], 1);
    coords[hip_r_rotation]   = getQandUIndices(s, m_bodies[thigh_r], 2);

    coords[knee_r_extension] = getQandUIndices(s, m_bodies[shank_r], 0);

    coords[ankle_r_dorsiflexion] = getQandUIndices(s, m_bodies[foot_r], 0);
    coords[ankle_r_inversion]    = getQandUIndices(s, m_bodies[foot_r], 1);

    coords[mtp_r_dorsiflexion] = getQandUIndices(s, m_bodies[toes_r], 0);

    // left.
    coords[hip_l_flexion]    = getQandUIndices(s, m_bodies[thigh_l], 0);
    coords[hip_l_adduction]  = getQandUIndices(s, m_bodies[thigh_l], 1);
    coords[hip_l_rotation]   = getQandUIndices(s, m_bodies[thigh_l], 2);

    coords[knee_l_extension] = getQandUIndices(s, m_bodies[shank_l], 0);

    coords[ankle_l_dorsiflexion] = getQandUIndices(s, m_bodies[foot_l], 0);
    coords[ankle_l_inversion]    = getQandUIndices(s, m_bodies[foot_l], 1);

    coords[mtp_l_dorsiflexion] = getQandUIndices(s, m_bodies[toes_l], 0);
}

void Biped::findContactForces(const State& s, Real& fLeft, Real& fRight) const
    {
        const unsigned int nContacts = m_contact.getNumContactForces(s);
        const ContactSnapshot& snapshot = m_tracker.getActiveContacts(s);

        fLeft = 0;
        fRight = 0;
        for (unsigned int i = 0; i < nContacts; ++i)
        {
            const ContactForce& force = m_contact.getContactForce(s, i);
            const ContactId id = force.getContactId();
            assert(snapshot.hasContact(id));
            const Contact& contact = snapshot.getContactById(id);
            const MobilizedBody& b1 = m_tracker.getMobilizedBody(contact.getSurface1());
            const MobilizedBody& b2 = m_tracker.getMobilizedBody(contact.getSurface2());
            const bool left = isLeftFoot(b1) || isLeftFoot(b2);
            const bool right = isRightFoot(b1) || isRightFoot(b2);
            if (left)
            {
                fLeft += force.getForceOnSurface2()[1].norm();
            }
            else if (right)
            {
                fRight += force.getForceOnSurface2()[1].norm();
            }
        }
    }

void Biped::findContactStatus(const State& s, bool& left, bool& right) const
{
    Real fLeft;
    Real fRight;
    findContactForces(s, fLeft, fRight);
    left = fLeft > 0;
    right = fRight > 0;
}


// SIMBICON

#define STATE_UPD_STEPSIZE 0.0005
//#define RIGID_CONTACT
 
class Simbicon : public SimTK::Force::Custom::Implementation
{
public:
	Simbicon(Biped& model); 
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
	
    Biped& _model;
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
    SimbiconStateHandler(Biped& m, Simbicon& simctrl, SimTK::Real interval);

    void handleEvent(SimTK::State& s, SimTK::Real accuracy, bool& shouldTerminate) const;
private:
    Biped& _model;
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

Simbicon::Simbicon(Biped& model) : _model(model) {
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
        mobilityForces[_model.getUIndex(Biped::Coordinate(i))] = controls[i];
}


// Project the pelvis z (right) and x (forward) directions onto the Ground
// (x-z) plane, and normalize. Projected z is the normal to the sagittal plane;
// projected x is the normal to the coronal plane.
void Simbicon::getSagCorNormals(const State& s, Vec3& sagN, Vec3& corN ) const {
	const MobilizedBody& pelvis = _model.getBody(Biped::pelvis);
	
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
	const SimbodyMatterSubsystem& matter = _model.getMatterSubsystem();
	Vec3 sagN, corN; 	
	
    // Sagittal and coronal plane normals are right and front directions of
    // the pelvis, projected on the ground plane.
	getSagCorNormals(s, sagN, corN); 	
	
	int swh = Biped::hip_r_flexion;             // swing hip
	int sth = Biped::hip_l_flexion;             // stance hip
    MobilizedBody ankle = _model.getBody(Biped::foot_l); // stance ankle
	
	if (_state == 2 || _state == 3) { // right stance
		swh = Biped::hip_l_flexion; 
		sth = Biped::hip_r_flexion; 
        ankle = _model.getBody(Biped::foot_r);
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
	int swh = Biped::hip_r_flexion; 
	int sth = Biped::hip_l_flexion; 
	
	if (_state == 2 || _state == 3) { // right stance
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

        const Real qi = s.getQ()[_model.getQIndex(Biped::Coordinate(i))];
        const Real ui = s.getU()[_model.getUIndex(Biped::Coordinate(i))];
		controls[i] = clamp(kp*(thetad - qi) - kd*ui, kp); 
	}
	
	if (_state >= 0) 
		fillInHipJointControls(s, controls); 

	Vec3 com = _model.getMatterSubsystem().calcSystemMassCenterLocationInGround(s); 
	if (com[1] < 0.7) {
		for (int i = 0; i < controls.size(); i++) 
			controls[i] = 0.0; 
	}
	return;
#endif
}
	
void Simbicon::
computeSecondaryStateVals(const State& s, Real lForce, Real rForce) {
	const MobilizedBody& pelvis = _model.getBody(Biped::pelvis);

#ifndef RIGID_CONTACT
    _lastRFootContactForce = _curRFootContactForce;
    _lastLFootContactForce = _curLFootContactForce;
    _curRFootContactForce = rForce;
    _curLFootContactForce = lForce;
#endif

	Vec3 upThigh;
	if (_state == 0 || _state == 1)  
		getUpVectorInGround(s, _model.getBody(Biped::thigh_r), upThigh); 
	else if (_state == 2 || _state == 3) 
		getUpVectorInGround(s, _model.getBody(Biped::thigh_l), upThigh); 
	
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
	const MobilizedBody& rFoot = _model.getBody(Biped::foot_r);
	const MobilizedBody& lFoot = _model.getBody(Biped::foot_l);
	
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

SimbiconStateHandler::SimbiconStateHandler(Biped& model, Simbicon& simctrl,
    Real interval) 
        : PeriodicEventHandler(interval), _model(model), _simctrl(simctrl) {
	}

void SimbiconStateHandler::handleEvent(State& s, Real accuracy, bool& shouldTerminate) const
{
    shouldTerminate = false;
	Simbicon* simctrl = &_simctrl;

    _model.realize(s, Stage::Dynamics);

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
    OutputReporter(const Biped& model, 
                   Real            interval); 
    void handleEvent(const State& state) const OVERRIDE_11;
private:
    const Biped& _model;
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

    Biped model;
    SimTK::Visualizer                   viz(model);
    SimTK::Visualizer::InputSilo* userInput = new Visualizer::InputSilo();
    viz.addInputListener(userInput);
    model.addEventHandler
       (new UserInputHandler(*userInput, Real(0.1))); //check input every 100ms

    model.addEventReporter(new OutputReporter(model, .01));
    model.addEventReporter(new Visualizer::Reporter(viz, RealTimeFactor/30));
    DecorativeText help("Any input to start; ESC to quit");
    help.setIsScreenText(true);
    viz.addDecoration(MobilizedBodyIndex(0),Vec3(0),help);
    model.updMatterSubsystem().setShowDefaultGeometry(false);

    // Add the controller.
    Simbicon* simctrl = new Simbicon(model);
    Force::Custom simbicon(model.updForceSubsystem(), simctrl); // takes ownership

#ifndef RIGID_CONTACT
    model.addEventHandler(new SimbiconStateHandler(model,*simctrl,
                                                          STATE_UPD_STEPSIZE));
#endif

    State s;
    model.initialize(s);

    printf("Act: u\n");
    for (int i=0; i < NumActuators; ++i) {
        printf("%2d: %d\n", i, int(model.getUIndex(Biped::Coordinate(i))));
    }

    //model.toes_r.lockAt(s, .2); model.toes_l.lockAt(s, .2);
    //model.foot_r.lockAt(s, Vec2(.2,0)); model.foot_l.lockAt(s, Vec2(.2,0));
    model.realize(s, Stage::Instance);

    printf("SIMBICON 3D:\n");
    printf("%d bodies, %d mobilities, -%d constraint equations -%d motions\n",
        model.getMatterSubsystem().getNumBodies(), s.getNU(), s.getNMultipliers(), 
        model.getMatterSubsystem().getKnownUDotIndex(s).size());

    model.getBody(Biped::trunk).setQToFitTranslation(s, Vec3(0,1.5,0));
    model.getBody(Biped::trunk).setUToFitLinearVelocity(s, Vec3(1,0,0));
    viz.report(s); 

#ifdef RIGID_CONTACT
    viz.addDecorationGenerator(new ShowContact(model));  
#endif

#ifndef RIGID_CONTACT
    // Simulate.
    //CPodesIntegrator integ(model.system); integ.setOrderLimit(2); integ.setAccuracy(.01);
    
    //RungeKuttaMersonIntegrator integ(model.system); integ.setAccuracy(1e-3);
    //RungeKutta2Integrator integ(model.system); integ.setAccuracy(.1);
    SemiExplicitEuler2Integrator integ(model); integ.setAccuracy(0.1);
    //SemiImplicitEulerIntegrator integ(model.system, .002);
    //integ.setConstraintTolerance(.001);
    //integ.setMaximumStepSize(.005);
    TimeStepper ts(model, integ);
#else
    SemiExplicitEulerTimeStepper ts(model);
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

        model.realize(s, Stage::Dynamics);

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
        model.realize(ts.getState());
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
OutputReporter::OutputReporter(const Biped& model, 
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
