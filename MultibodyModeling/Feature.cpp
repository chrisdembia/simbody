/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * Implementations of API-level multibody modeling objects for Simbody.
 */

#include "SimbodyCommon.h"
#include "Feature.h"
#include "FeatureRep.h"

#include <string>
#include <iostream> 
#include <sstream>


namespace simtk {

    // FEATURE //

Feature::Feature(const Feature& f) : rep(0) {
    if (f.rep) 
        f.rep->cloneWithoutParentOrExternalPlacements(*this);
}
Feature& Feature::operator=(const Feature& f) {
    if (this != &f) {
        // This will blow up if rep doesn't have a handle -- we shouldn't
        // be pointing to it in that case!
        if (rep && (&rep->getMyHandle() == this)) delete rep; 
        rep = 0;
        if (f.rep) 
            f.rep->cloneWithoutParentOrExternalPlacements(*this);
    }
    return *this;
}
Feature::~Feature() {
    // This will blow up if rep doesn't have a handle -- we shouldn't
    // be pointing to it in that case!
    if (rep && (&rep->getMyHandle() == this)) delete rep; 
    rep = 0;
}
static String featureHasNoRep(const Feature& f) {
    std::ostringstream s;
    s << "<FEATURE AT 0x" << &f << " WITH NULL REP>";
    return String(s.str());
}
String Feature::getName() const {
    return rep ? String(rep->getName()) : featureHasNoRep(*this);
}
String Feature::getFullName() const {
    return rep ? String(rep->getFullName()) : featureHasNoRep(*this);
}
String Feature::getFeatureTypeName() const {
    return rep ? String(rep->getFeatureTypeName()) 
               : featureHasNoRep(*this);
}

bool Feature::hasPlacement() const {
    return getRep().hasPlacement();
}
void Feature::place(const Placement& p) {
    try {
        updRep().place(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Feature::place", exc.getMessage());
    }
}

const Placement& Feature::getPlacement() const {
    assert(hasPlacement());
    return getRep().getPlacement();
}

String Feature::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "Feature ";
    if (!rep) {
        s << featureHasNoRep(*this);
        return String(s.str());
    }

    const FeatureRep& f = *rep;
    s << f.getFeatureTypeName() << " " << f.getFullName() << ": ";
    s << (f.hasPlacement() ? f.getPlacement().toString(linePrefix)
                           : String("NO PLACEMENT"));

    const size_t nSubfeatures     = f.getNSubfeatures();
    const size_t nPlacement       = f.getNPlacementExpressions();
    const size_t nPlacementValues = f.getNPlacementValues();
    const std::string nextIndent  = linePrefix + "    ";

    if (nSubfeatures) {
        s << std::endl << linePrefix << "  Subfeatures (" << nSubfeatures << "):";
        for (size_t i=0; i < nSubfeatures; ++i)
            s  << std::endl << nextIndent << f.getSubfeature(i).toString(nextIndent);
    }
    if (nPlacement) {
        s << std::endl << linePrefix << "  Placement Expressions (" << nPlacement << "):";
        for (size_t i=0; i < nPlacement; ++i)
            s  << std::endl << nextIndent << f.getPlacementExpression(i).toString(nextIndent);
    }
    if (nPlacementValues) {
        s << std::endl << linePrefix << "  Placement Values (" << nPlacementValues << "):";
        for (size_t i=0; i < nPlacementValues; ++i)
            s  << std::endl << nextIndent << f.getPlacementValue(i).toString(nextIndent);
    }
    return s.str();
}

std::ostream& operator<<(std::ostream& o, const Feature& f) {
    return o << f.toString() << std::endl;
}
bool Feature::hasParentFeature() const {
    return getRep().hasParentFeature();
}
int Feature::getIndexInParent() const {
    assert(getRep().hasParentFeature());
    return getRep().getIndexInParent();
}
const Feature& Feature::getParentFeature() const {
    assert(getRep().hasParentFeature());
    return getRep().getParentFeature();
}
bool Feature::dependsOn(const Feature& f) const {
    if (isSameFeature(f)) return true;
    return hasRep() && getRep().dependsOn(f);
}

int Feature::getNSubfeatures() const
  { return getRep().getNSubfeatures(); }
const Feature& Feature::getSubfeature(int i) const
  { return getRep().getSubfeature(i); }
Feature& Feature::updSubfeature(int i)
  { return updRep().updSubfeature(i); }

// getXXX() methods
const Feature& Feature::getSubfeature(const String& n) const 
  { return getRep().getSubfeature(n); }
const RealParameter& Feature::getRealParameter(const String& n) const
  { return RealParameter::downcast(getSubfeature(n)); }
const StationParameter& Feature::getStationParameter(const String& n) const
  { return StationParameter::downcast(getSubfeature(n)); }
const RealMeasure& Feature::getRealMeasure(const String& n) const
  { return RealMeasure::downcast(getSubfeature(n)); }
const StationMeasure& Feature::getStationMeasure(const String& n) const
  { return StationMeasure::downcast(getSubfeature(n)); }
const Station& Feature::getStation(const String& n) const
  { return Station::downcast(getSubfeature(n)); }
const Direction& Feature::getDirection(const String& n) const
  { return Direction::downcast(getSubfeature(n)); }
const Orientation& Feature::getOrientation(const String& n) const
  { return Orientation::downcast(getSubfeature(n)); }
const Frame& Feature::getFrame(const String& n) const
  { return Frame::downcast(getSubfeature(n)); }

// updXXX() methods
Feature& Feature::updSubfeature(const String& n)
  { return updRep().updSubfeature(n); }
RealParameter& Feature::updRealParameter(const String& n)
  { return RealParameter::downcast(updSubfeature(n)); }
StationParameter& Feature::updStationParameter(const String& n)
  { return StationParameter::downcast(updSubfeature(n)); }
RealMeasure& Feature::updRealMeasure(const String& n)
  { return RealMeasure::downcast(updSubfeature(n)); }
StationMeasure& Feature::updStationMeasure(const String& n)
  { return StationMeasure::downcast(updSubfeature(n)); }
Station& Feature::updStation(const String& n)
  { return Station::downcast(updSubfeature(n)); }
Direction& Feature::updDirection(const String& n)
  { return Direction::downcast(updSubfeature(n)); }
Orientation& Feature::updOrientation(const String& n)
  { return Orientation::downcast(updSubfeature(n)); }
Frame& Feature::updFrame(const String& n)
  { return Frame::downcast(updSubfeature(n)); }

// addXXX() methods
RealParameter& Feature::addRealParameter(const String& n, const Placement& p) {
    RealParameter& rp = RealParameter::downcast(updRep().addSubfeatureLike(RealParameter(n), n));
    if (p.hasRep()) rp.place(p);
    return rp;
}
RealMeasure& Feature::addRealMeasure(const String& n, const Placement& p) {
    RealMeasure& rm = RealMeasure::downcast(updRep().addSubfeatureLike(RealMeasure(n), n));
    if (p.hasRep()) rm.place(p);
    return rm;
}
StationParameter& Feature::addStationParameter(const String& n, const Placement& p) {
    StationParameter& sp = StationParameter::downcast(updRep().addSubfeatureLike(StationParameter(n), n));
    if (p.hasRep()) sp.place(p);
    return sp;
}
StationMeasure& Feature::addStationMeasure(const String& n, const Placement& p) {
    StationMeasure& sm = StationMeasure::downcast(updRep().addSubfeatureLike(StationMeasure(n), n));
    if (p.hasRep()) sm.place(p);
    return sm;
}
Station& Feature::addStation(const String& n, const Placement& p) {
    Station& s = Station::downcast(updRep().addSubfeatureLike(Station(n), n));
    if (p.hasRep()) s.place(p);
    return s;
}
Direction& Feature::addDirection(const String& n, const Placement& p) {
    Direction& d = Direction::downcast(updRep().addSubfeatureLike(Direction(n), n));
    if (p.hasRep()) d.place(p);
    return d;
}
Orientation& Feature::addOrientation(const String& n, const Placement& p) {
    Orientation& o = Orientation::downcast(updRep().addSubfeatureLike(Orientation(n), n));
    if (p.hasRep()) o.place(p);
    return o;
}
Frame& Feature::addFrame(const String& n, const Placement& p) {
    Frame& f = Frame::downcast(updRep().addSubfeatureLike(Frame(n), n));
    if (p.hasRep()) f.place(p);
    return f;
}

Feature& Feature::addSubfeatureLike(const Feature& f, const String& n, const Placement& p) {
    Feature& fnew = updRep().addSubfeatureLike(f,n);
    if (p.hasRep()) fnew.place(p);
    return fnew;
}

    // REAL PARAMETER //

RealParameter::RealParameter(const String& nm) : RealMeasure() {
    rep = new RealParameterRep(*this, std::string(nm)); 
    rep->initializeStandardSubfeatures();
}
RealParameter::RealParameter(const RealParameter& src) : RealMeasure(src) { }
RealParameter& RealParameter::operator=(const RealParameter& src)
  { RealMeasure::operator=(src); return *this; }
RealParameter::~RealParameter() { }

/*static*/ bool             
RealParameter::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return RealParameterRep::isA(f.getRep());
}
/*static*/ const RealParameter& 
RealParameter::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const RealParameter&>(f);
}

/*static*/ RealParameter&       
RealParameter::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<RealParameter&>(f);
}

    // VEC3 PARAMETER //

Vec3Parameter::Vec3Parameter(const String& nm) : Vec3Measure() {
    rep = new Vec3ParameterRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Vec3Parameter::Vec3Parameter(const Vec3Parameter& src) : Vec3Measure(src) { }
Vec3Parameter& Vec3Parameter::operator=(const Vec3Parameter& src)
  { Vec3Measure::operator=(src); return *this; }
Vec3Parameter::~Vec3Parameter() { }

/*static*/ bool             
Vec3Parameter::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return Vec3ParameterRep::isA(f.getRep());
}
/*static*/ const Vec3Parameter& 
Vec3Parameter::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Vec3Parameter&>(f);
}

/*static*/ Vec3Parameter&       
Vec3Parameter::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Vec3Parameter&>(f);
}

    // STATION PARAMETER //
StationParameter::StationParameter(const String& nm) : StationMeasure() {
    rep = new StationParameterRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
StationParameter::StationParameter(const StationParameter& src) : StationMeasure(src) { }
StationParameter& StationParameter::operator=(const StationParameter& src)
  { StationMeasure::operator=(src); return *this; }
StationParameter::~StationParameter() { }

/*static*/ bool             
StationParameter::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return StationParameterRep::isA(f.getRep());
}
/*static*/ const StationParameter& 
StationParameter::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const StationParameter&>(f);
}

/*static*/ StationParameter&       
StationParameter::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<StationParameter&>(f);
}

    // REAL MEASURE //
RealMeasure::RealMeasure(const String& nm) : Feature() {
    rep = new RealMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
RealMeasure::RealMeasure(const RealMeasure& src) : Feature(src) { }
RealMeasure& RealMeasure::operator=(const RealMeasure& src)
  { Feature::operator=(src); return *this; }
RealMeasure::~RealMeasure() { }

/*static*/ bool             
RealMeasure::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return RealMeasureRep::isA(f.getRep());
}
/*static*/ const RealMeasure& 
RealMeasure::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const RealMeasure&>(f);
}

/*static*/ RealMeasure&       
RealMeasure::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<RealMeasure&>(f);
}

    // VEC3 MEASURE //
Vec3Measure::Vec3Measure(const String& nm) : Feature() {
    rep = new Vec3MeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Vec3Measure::Vec3Measure(const Vec3Measure& src) : Feature(src) { }
Vec3Measure& Vec3Measure::operator=(const Vec3Measure& src)
  { Feature::operator=(src); return *this; }
Vec3Measure::~Vec3Measure() { }

/*static*/ bool             
Vec3Measure::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return Vec3MeasureRep::isA(f.getRep());
}
/*static*/ const Vec3Measure& 
Vec3Measure::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Vec3Measure&>(f);
}

/*static*/ Vec3Measure&       
Vec3Measure::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Vec3Measure&>(f);
}

    // STATION MEASURE //
StationMeasure::StationMeasure(const String& nm) : Feature() {
    rep = new StationMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
StationMeasure::StationMeasure(const StationMeasure& src) : Feature(src) { }
StationMeasure& StationMeasure::operator=(const StationMeasure& src)
  { Feature::operator=(src); return *this; }
StationMeasure::~StationMeasure() { }

/*static*/ bool             
StationMeasure::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return StationMeasureRep::isA(f.getRep());
}
/*static*/ const StationMeasure& 
StationMeasure::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const StationMeasure&>(f);
}

/*static*/ StationMeasure&       
StationMeasure::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<StationMeasure&>(f);
}

    // STATION //
Station::Station(const String& nm) : Feature() {
    rep = new StationRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Station::Station(const Station& src) : Feature(src) { }
Station& Station::operator=(const Station& src)
  { Feature::operator=(src); return *this; }
Station::~Station() { }

/*static*/ bool             
Station::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return StationRep::isA(f.getRep());
}
/*static*/ const Station& 
Station::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Station&>(f);
}

/*static*/ Station&       
Station::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Station&>(f);
}

    // DIRECTION MEASURE //
DirectionMeasure::DirectionMeasure(const String& nm) : Feature() {
    rep = new DirectionMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
DirectionMeasure::DirectionMeasure(const DirectionMeasure& src) : Feature(src) { }
DirectionMeasure& DirectionMeasure::operator=(const DirectionMeasure& src)
  { Feature::operator=(src); return *this; }
DirectionMeasure::~DirectionMeasure() { }

/*static*/ bool             
DirectionMeasure::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return DirectionMeasureRep::isA(f.getRep());
}
/*static*/ const DirectionMeasure& 
DirectionMeasure::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const DirectionMeasure&>(f);
}

/*static*/ DirectionMeasure&       
DirectionMeasure::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<DirectionMeasure&>(f);
}

    // DIRECTION //
Direction::Direction(const String& nm) : Feature() {
    rep = new DirectionRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Direction::Direction(const Direction& src) : Feature(src) { }
Direction& Direction::operator=(const Direction& src)
  { Feature::operator=(src); return *this; }
Direction::~Direction() { }

/*static*/ bool             
Direction::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return DirectionRep::isA(f.getRep());
}
/*static*/ const Direction& 
Direction::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Direction&>(f);
}

/*static*/ Direction&       
Direction::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Direction&>(f);
}

    // ORIENTATION MEASURE //
OrientationMeasure::OrientationMeasure(const String& nm) : Feature() {
    rep = new OrientationMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
OrientationMeasure::OrientationMeasure(const OrientationMeasure& src) : Feature(src) { }
OrientationMeasure& OrientationMeasure::operator=(const OrientationMeasure& src)
  { Feature::operator=(src); return *this; }
OrientationMeasure::~OrientationMeasure() { }

/*static*/ bool             
OrientationMeasure::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return OrientationMeasureRep::isA(f.getRep());
}
/*static*/ const OrientationMeasure& 
OrientationMeasure::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const OrientationMeasure&>(f);
}

/*static*/ OrientationMeasure&       
OrientationMeasure::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<OrientationMeasure&>(f);
}

    // ORIENTATION //
Orientation::Orientation(const String& nm) : Feature() {
    rep = new OrientationRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Orientation::Orientation(const Orientation& src) : Feature(src) { }
Orientation& Orientation::operator=(const Orientation& src)
  { Feature::operator=(src); return *this; }
Orientation::~Orientation() { }

const Direction& 
Orientation::getAxis(int i) const {
    return OrientationRep::downcast(getRep()).getAxis(i);
}

/*static*/ bool             
Orientation::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return OrientationRep::isA(f.getRep());
}
/*static*/ const Orientation& 
Orientation::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Orientation&>(f);
}

/*static*/ Orientation&       
Orientation::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Orientation&>(f);
}

    // FRAME //
Frame::Frame(const String& nm) : Feature() {
    rep = new FrameRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Frame::Frame(const Frame& src) : Feature(src) { }
Frame& Frame::operator=(const Frame& src)
  { Feature::operator=(src); return *this; }
Frame::~Frame() { }

const Orientation& Frame::getOrientation() const {
    return FrameRep::downcast(getRep()).getOrientation();
}
const Station& Frame::getOrigin() const {
    return FrameRep::downcast(getRep()).getOrigin();
}

/*static*/ bool             
Frame::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return FrameRep::isA(f.getRep());
}
/*static*/ const Frame& 
Frame::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Frame&>(f);
}

/*static*/ Frame&       
Frame::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Frame&>(f);
}

} // namespace simtk



static int caseInsensitiveCompare(const std::string& key, const std::string& test) {
    const size_t minlen = std::min(key.size(), test.size());
    for (size_t i=0; i < minlen; ++i) {
        const int k = tolower(key[i]), t = tolower(test[i]);
        if (k < t) return -1;
        else if (k > t) return 1;
    }
    // caution -- size() is unsigned, don't get clever here
    if (key.size() > minlen) return 1;
    else if (test.size() > minlen) return -1;
    return 0;
}
