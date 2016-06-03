/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-16 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
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

#include "SimTKcommon.h"

using namespace SimTK;

// This is an approximation to an equality operator (State does not have an
// equality operator).
bool operator==(const SimTK::Vector& vA, const SimTK::Vector& vB) {
    if (vA.size() != vB.size()) return false;
    for (int i = 0; i < vA.size(); ++i) {
        if (vA[i] != vB[i]) return false;
    }
    return true;
}
void testSomewhatEqual(const State& sA, const State& sB) {
    SimTK_TEST(sA.isConsistent(sB));
    SimTK_TEST(   (isNaN(sA.getTime()) && isNaN(sB.getTime()))
               ||  sA.getTime() == sB.getTime());
    SimTK_TEST(sA.getQ() == sB.getQ());
    SimTK_TEST(sA.getU() == sB.getU());
    SimTK_TEST(sA.getZ() == sB.getZ());
    SimTK_TEST(sA.getUWeights() == sB.getUWeights());
    SimTK_TEST(sA.getZWeights() == sB.getZWeights());
    SimTK_TEST(sA.getQErrWeights() == sB.getQErrWeights());
    SimTK_TEST(sA.getUErrWeights() == sB.getUErrWeights());
}

void testEmptyState() {

    // Test a roundtrip for an empty state.
    {
        State s0;
        auto xml0 = s0.toXmlElement();
        State s1;
        s1.fromXmlElement(xml0);
        testSomewhatEqual(s0, s1);
    }

    // Test naming the xml element.
    {
        State s0;
        auto xml0 = s0.toXmlElement("nork9");
        State s1;
        s1.fromXmlElement(xml0);
        testSomewhatEqual(s0, s1);
    }
    {
        State s0;
        auto xml0 = s0.toXmlElement("nork9");
        State s1;
        s1.fromXmlElement(xml0, "nork9");
        testSomewhatEqual(s0, s1);
    }
    {
        State s0;
        auto xml0 = s0.toXmlElement("nork9");
        State s1;
        SimTK_TEST_MUST_THROW(s1.fromXmlElement(xml0, "nork8"));
    }
    {
        State s0;
        auto xml0 = s0.toXmlElement();
        State s1;
        SimTK_TEST_MUST_THROW(s1.fromXmlElement(xml0, "nork8"));
    }
}

// test giving a bogus Xml element.

// test filling up a state that already had stuff in it (different size, etc.)
//
// tamper with the Xml format version to catch the version exception.
// q
// u
// z
// qerr
// uerr
// udoterr
// multiple subsystems
// different stages of realization?
// test all the measures

void testStateXml() {
}

int main() {
    SimTK_START_TEST("TestStateXml");
        SimTK_SUBTEST(testEmptyState);
        SimTK_SUBTEST(testStateXml);
    SimTK_END_TEST();
}
