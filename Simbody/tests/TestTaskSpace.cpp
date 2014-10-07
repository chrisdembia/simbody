/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
 * Contributors: Michael Sherman                                              *
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

#include "Simbody.h"

using namespace SimTK;
using namespace std;

class DoublePendulumOnCart : public MultibodySystem {
public:
    DoublePendulumOnCart() : m_matter(*this), m_forces(*this) {
        const Real m1 = 1.0;
        const Real m2 = 1.0;
        const Real m3 = 1.0;
        const Real L2 = 1.0;
        const Real L3 = 1.0;

        Force::Gravity(m_forces, m_matter, -YAxis, 9.8);

        Body::Rigid cartInfo(MassProperties(m1, 0, Inertia(0)));
        Body::Rigid link2Info(MassProperties(m2, 0, Inertia(0)));
        Body::Rigid link3Info(MassProperties(m3, 0, Inertia(0)));

        Transform cartTransform(Vec3(0)); // TODO Rotation(-0.5 * Pi, YAxis));
        MobilizedBody::Slider cart(
                m_matter.Ground(), cartTransform,
                cartInfo, cartTransform);

        Transform pinTransform(Vec3(0));
        MobilizedBody::Pin link2(cart, pinTransform, link2Info, pinTransform);
        MobilizedBody::Pin link3(link2, pinTransform, link3Info, pinTransform);

    }

private:
    SimbodyMatterSubystem m_matter;
    GeneralForceSubsystem m_forces;

    const Real m1;
    const Real m2;
    const Real m3;
    const Real L2;
    const Real L3;
};

void testTaskSpace() {
}

int main() {
    try {
        cout << "*** TEST TASKSPACE ***\n\n"; testTaskSpace();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
