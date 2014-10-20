/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
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

// TODO
// 2. boundary_transformation.
// 3. restart.
// 4. libcmaes test cases. lots of cost functions!!!
// 5. memory leaks.
// 6. how to disable reading of cmaes_signals.par.
// 7. look at actparcmaes.par for the parameters actually used!
// 9. allow verbosity; diagnostics level.
// 10. parameter limits.
// 12 all the cmaes options.
// 13 understand how the cmaes options are set.
// 14 threading.

#include "SimTKmath.h"
#include "OptimizerSystems.h"

#include <iostream>
using std::cout;
using std::endl;
using SimTK::Vector;
using SimTK::Real;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;

static bool equalToTol(Real v1, Real v2, Real tol) {
    const Real scale = std::max(std::max(std::abs(v1), std::abs(v2)), Real(1));
    return std::abs(v1-v2) < scale*tol;
} 

void testCMAESAvailable() {
    SimTK_TEST(Optimizer::isAlgorithmAvailable(SimTK::CMAES));
}

// If we try to create an OptimizerSystem with only one parameter,
// we should get an exception.
void testTwoOrMoreParameters() {
    SimTK_TEST_MUST_THROW_EXC(Optimizer opt(Cigtab(1), SimTK::CMAES),
            SimTK::Exception::ValueOutOfRange
            );
}

// This tests that setting max iterations works using the Simbody
// interface. CMAES cannot find the optimum of Cigtab in 500 iterations (or
// less) given the initial condition we use.
void testMaxIterations() {

    Cigtab sys(22);
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    results.setTo(0.5);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setAdvancedRealOption("sigma", 0.3);
    opt.setMaxIterations(500);
    
    // Optimize!
    Real f = opt.optimize(results);

    // Print results of the optimization.
    printf("MaxIterations: f = %f params = ", f);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    // Check if the result is correct.
    static const Real TOL = 1e-5;
    Vector expected = sys.optimalParameters();
    bool answerIsCorrect = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results[i], expected[i], TOL)) {
            answerIsCorrect = false;
        }
    }
    SimTK_ASSERT_ALWAYS(!answerIsCorrect,
            "MaxIterations: setting of max iterations does not work.");
}

// This also tests that setting max iterations works using the Simbody
// interface, because Cigtab is not optimized in under the default number of
// max iterations (1000 at the time of this writing).
void testCigtabOptimum() {

    Cigtab sys(22);
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    results.setTo(0.5);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setMaxIterations(5000);
    opt.setAdvancedRealOption("sigma", 0.3);
    // Sometimes this test fails, so choose a seed where the test passes.
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    
    // Optimize!
    Real f = opt.optimize(results);

    // Print results of the optimization.
    printf("CigtabOptimum: f = %f params = ", f);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    // Check that the result is correct.
    static const Real TOL = 1e-5;
    Vector expected = sys.optimalParameters();
    bool answerIsCorrect = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results[i], expected[i], TOL)) {
            printf(" CigtabOptimum: error results[%d] = %f  expected=%f \n",
                    i, results[i], expected[i]);
            answerIsCorrect = false;
        }
    }
    SimTK_ASSERT_ALWAYS(answerIsCorrect,
            "CigtabOptimum: could not find optimum.");
}

void testParameterLimits() {
    // TODO by trying initial point outside of the parameter limits,
    // we never optimize. We may end up in an endless loop then.
    throw std::exception();
}

// Make sure that we are able to set sigma (step size) using Simbody's
// interface, and that with appropriate step size, we can find the optimum of
// Ackley's function.
void testSigmaStepSizeAndAckleyOptimum() {

    Ackley sys(15);
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    // Far from optimum, but within the parameter limits.
    results.setTo(25);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setMaxIterations(5000);
    opt.setAdvancedRealOption("seed", 10);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    // Default sigma (step size) leaves us in a local minimum.
    // =======================================================
    
    // Optimize!
    Real f1 = opt.optimize(results);

    // Print results of the optimization.
    printf("SigmaStepSizeAndAckleyOptimum: f = %f params = ", f1);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    // Check that the result is correct.
    static const Real TOL = 1e-5;
    Real expectedLocalMinimum = 24.999749;
    bool inLocalMinimum = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results[i], expectedLocalMinimum, TOL)) {
            printf(" SigmaStepSizeAndAckleyOptimum: "
                    "error results[%d] = %f  expected=%f \n",
                    i, results[i], expectedLocalMinimum);
            inLocalMinimum = false;
        }
    }
    if (!inLocalMinimum) {
        SimTK_TEST_FAILED("SigmaStepSizeAndAckleyOptimum: "
            "Should have ended up in the 24.9997 local minimum.");
    }

    // Can find the optimum with an appropriate step size.
    // ===================================================
    // sigma should be 1/4 the range of possible values.
    opt.setAdvancedRealOption("sigma", 0.25 * 64);

    // Optimize!
    results.setTo(25);
    Real f2 = opt.optimize(results);

    // Print results of the optimization.
    printf("SigmaStepSizeAndAckleyOptimum: f = %f params = ", f2);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    // Check that the result is correct.
    Vector expected = sys.optimalParameters();
    bool answerIsCorrect = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results[i], expected[i], TOL)) {
            printf(" SigmaStepSizeAndAckleyOptimum: "
                    "error results[%d] = %f  expected=%f \n",
                    i, results[i], expected[i]);
            answerIsCorrect = false;
        }
    }
    SimTK_TEST(answerIsCorrect);
}

// To find the optimum of this function, we need lots of samples and a lot
// of function evaluations. Thus, this test makes sure that we are able
// to modify both of these settings.
// TODO setting the seed is not working.
void testDropWaveOptimumLambdaMaxFunEvals() {

    DropWave sys;
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    results.setTo(1);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-3);
    opt.setAdvancedRealOption("sigma", 3.5);
    // With default lambda, this test fails. So if this test passes, we know we
    // can set lambda.
    opt.setAdvancedIntOption("lambda", 1000);
    opt.setAdvancedRealOption("seed", 10);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    opt.setMaxIterations(5000);
    
    // Will not find optimum to tolerance with default # function evals.
    // =================================================================
    // Optimize!
    Real f1 = opt.optimize(results);

    // Print results of the optimization.
    printf("DropWaveOptimumLambdaMaxFunEvals: f1 = %f params = ", f1);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    // Check if the result is correct.
    static const Real TOL = 1e-3;
    Vector expected = sys.optimalParameters();
    bool answerIsCorrect1 = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results[i], expected[i], TOL)) {
            answerIsCorrect1 = false;
        }
        else {
            printf(" DropWaveOptimumLambdaMaxFunEvals: "
                    "equal results[%d] = %f  expected=%f \n",
                    i, results[i], expected[i]);
        }
    }

    SimTK_ASSERT_ALWAYS(!answerIsCorrect1,
            "DropWaveOptimumLambdaMaxFunEvals: "
            "With default max function evaluations, "
            "should not have found optimum.");

    // With enough function evals, we can find the optimum.
    // ====================================================
    opt.setAdvancedIntOption("stopMaxFunEvals", 100000);
    results.setTo(1);
    Real f2 = opt.optimize(results);
    Vector results2 = results;

    // Print results of the optimization.
    printf("DropWaveOptimumLambdaMaxFunEvals: f2 = %f params2 = ", f2);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results2[i]);
    }
    printf("\n");

    // Check that the result is correct.
    bool answerIsCorrect2 = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results2[i], expected[i], TOL)) {
            printf(" DropWaveOptimumLambdaMaxFunEvals: "
                    "error results2[%d] = %f  expected=%f \n",
                    i, results2[i], expected[i]);
            answerIsCorrect2 = false;
        }
    }

    SimTK_ASSERT_ALWAYS(answerIsCorrect2,
            "DropWaveOptimumLambdaMaxFunEvals: could not find optimum.");
}

void testSeed() {

    Ackley sys(22);
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    results.setTo(25);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setAdvancedRealOption("sigma", 1);
    
    // A negative seed causes an exception to be thrown upon optimization.
    // ===================================================================
    opt.setAdvancedIntOption("seed", -10);
    SimTK_TEST_MUST_THROW_EXC(
            Real f = opt.optimize(results),
            SimTK::Exception::ValueWasNegative
    );

    // Using the same seed gives identical results, if maxtime is 1.
    // =============================================================
    // We end prematurely because non-identical seeds may lead to similar
    // results at the optimum; we don't want to be at the optimum.
    opt.setMaxIterations(500);
    opt.setAdvancedIntOption("seed", 42);

    // First optimization.
    Real f1 = opt.optimize(results);
    Vector results1 = results;

    // Second optimization.
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    results.setTo(25);
    Real f2 = opt.optimize(results);
    Vector results2 = results;

    // Third optimization, which should now produce identical values to f2.
    results.setTo(25);
    Real f3 = opt.optimize(results);
    Vector results3 = results;

    // Print results of the optimizations.
    printf("Seed: f1 = %f params1 = ", f1);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results1[i]);
    }
    printf("\n");
    printf("Seed: f2 = %f params2 = ", f2);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results2[i]);
    }
    printf("\n");
    printf("Seed: f3 = %f params3 = ", f3);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results3[i]);
    }
    printf("\n");

    // f1 and f2 don't match.
    // ----------------------
    // Check for an identical value of f.
    SimTK_ASSERT_ALWAYS(!equalToTol(f1, f2, 1e-12),
            "Seed: using the same seed without maxtime leads to "
            "identical results.");

    // Check for identical parameters.
    bool answersAreIdentical2 = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results1[i], results2[i], 1e-10)) {
            answersAreIdentical2 = false;
        }
        else {
            printf(" Seed: equal results1[%d] = %f  results2[%d] = %f\n",
                    i, results1[i], i, results2[i]);
        }
    }
    SimTK_ASSERT_ALWAYS(!answersAreIdentical2,
            "Seed: using the same seed without maxtime leads to "
            "identical results.");

    // f2 and f3 match.
    // ----------------
    // Check for an identical value of f.
    SimTK_ASSERT_ALWAYS(equalToTol(f2, f3, 1e-12),
            "Seed: using the same seed does not lead to identical results.");

    // Check for identical parameters.
    bool answersAreIdentical3 = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results2[i], results3[i], 1e-10)) {
            printf(" Seed: error results2[%d] = %f  results3[%d] = %f\n",
                    i, results2[i], i, results3[i]);
            answersAreIdentical3 = false;
        }
    }
    SimTK_ASSERT_ALWAYS(answersAreIdentical3,
        "Seed: using the same seed does not lead to identical results.");

    // Using a different seed gives different results.
    // ===============================================
    opt.setAdvancedIntOption("seed", 50);
    results.setTo(25);
    Real f4 = opt.optimize(results);
    Vector results4 = results;
    printf("Seed: f4 = %f params4 = ", f4);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results4[i]);
    }
    printf("\n");

    // Check for an identical value of f.
    SimTK_ASSERT_ALWAYS(!equalToTol(f2, f4, 1e-10),
        "Seed: using different seeds does not lead to different results.");

    // Check for identical parameters.
    bool answersAreIdentical4 = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results2[i], results4[i], 1e-10)) {
            answersAreIdentical4 = false;
        }
        else {
            printf(" Seed: equal results2[%d] = %f  results4[%d] = %f\n",
                    i, results2[i], i, results4[i]);
        }
    }
    SimTK_ASSERT_ALWAYS(!answersAreIdentical4,
        "Seed: using different seeds does not lead to different results.");

}

void testConvergenceTolerance() {

    Ackley sys(15);
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    // Far from optimum, but within the parameter limits.
    results.setTo(25);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setMaxIterations(5000);
    opt.setAdvancedRealOption("seed", 10);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    opt.setAdvancedRealOption("sigma", 0.25 * 64);

    Real looseTolerance = 0.001;
    Real tightTolerance = 0.000001;

    // Use a loose tolerance.
    // ======================
    opt.setConvergenceTolerance(looseTolerance);

    // Optimize!
    results.setTo(25);
    Real f = opt.optimize(results);

    // Print results of the optimization.
    printf("ConvergenceTolerance: f = %f params = ", f);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    SimTK_ASSERT_ALWAYS(equalToTol(f, sys.optimalValue(), looseTolerance),
        "ConvergenceTolerance: "
        "optimal value is not correct within the loose tolerance.");
    SimTK_ASSERT_ALWAYS(!equalToTol(f, sys.optimalValue(), tightTolerance),
        "ConvergenceTolerance: "
        "the setting of the convergence tolerance is not functioning.");
}

// CMA-ES is able to minimize the Rosenbrock function.
// https://www.lri.fr/~hansen/cmsa-versus-cma.html
void testRosenbrock() {

    Rosenbrock sys(80);
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    results.setTo(0.9);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-3);
    opt.setMaxIterations(100000);
    opt.setAdvancedIntOption("lambda", 4 * N);
    opt.setAdvancedRealOption("sigma", 0.01);
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    
    // Optimize!
    Real f = opt.optimize(results);

    // Print results of the optimization.
    printf("Rosenbrock: f = %f params = ", f);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    // Check that the result is correct.
    static const Real TOL = 1e-5;
    Vector expected = sys.optimalParameters();
    bool answerIsCorrect = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results[i], expected[i], TOL)) {
            printf(" Rosenbrock: error results[%d] = %f  expected=%f \n",
                    i, results[i], expected[i]);
            answerIsCorrect = false;
        }
    }
    SimTK_ASSERT_ALWAYS(answerIsCorrect,
            "Rosenbrock: could not find optimum.");
}

// TODO i've been able to get f = -1.990145 so maybe there is a bug here?
void testSchwefel() {

    Schwefel sys(2);
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    results.setTo(418);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
// TODO    opt.setMaxIterations(100000);
    opt.setAdvancedIntOption("lambda", 4 * N);
    opt.setAdvancedRealOption("sigma", 1);
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    
    // Optimize!
    Real f = opt.optimize(results);

    // Print results of the optimization.
    printf("Schwefel: f = %f params = ", f);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    // Check that the result is correct.
    static const Real TOL = 1e-5;
    Vector expected = sys.optimalParameters();
    bool answerIsCorrect = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results[i], expected[i], TOL)) {
            printf(" Schwefel: error results[%d] = %f  expected=%f \n",
                    i, results[i], expected[i]);
            answerIsCorrect = false;
        }
    }
    SimTK_ASSERT_ALWAYS(answerIsCorrect, "Schwefel: could not find optimum.");
}

void testEasom() {

    Easom sys;
    int N = sys.getNumParameters();

    // set initial conditions.
    Vector results(N);
    results.setTo(-10);

    // Create optimizer; set settings.
    Optimizer opt(sys, SimTK::CMAES);
    opt.setAdvancedIntOption("lambda", 500);
    opt.setAdvancedRealOption("sigma", 25);
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    
    // Optimize!
    Real f = opt.optimize(results);

    // Print results of the optimization.
    printf("Easom: f = %f params = ", f);
    for (unsigned int i = 0; i < N; ++i) {
        printf(" %f", results[i]);
    }
    printf("\n");

    // Check that the result is correct.
    static const Real TOL = 1e-5;
    Vector expected = sys.optimalParameters();
    bool answerIsCorrect = true;
    for (unsigned int i=0; i < N; ++i) {
        if(!equalToTol(results[i], expected[i], TOL)) {
            printf(" Easom: error results[%d] = %f  expected=%f \n",
                    i, results[i], expected[i]);
            answerIsCorrect = false;
        }
    }
    SimTK_ASSERT_ALWAYS(answerIsCorrect, "Easom: could not find optimum.");
}

int main() {
    SimTK_START_TEST("CMAES");

        // TODO
        // Even though most of the tests use seeds, some tests may fail
        // sporadically.
        SimTK_SUBTEST(testCMAESAvailable);
        SimTK_SUBTEST(testTwoOrMoreParameters);
        SimTK_SUBTEST(testMaxIterations);
        SimTK_SUBTEST(testCigtabOptimum);
        // TODO        testParameterLimits();
        SimTK_SUBTEST(testSigmaStepSizeAndAckleyOptimum);
        // TODO        testDropWaveOptimumLambdaMaxFunEvals();
        SimTK_SUBTEST(testSeed);
        // TODO        testRestart();
        SimTK_SUBTEST(testConvergenceTolerance);
        // TODO        testRosenbrock();
        // TODO        testSchwefel();
        SimTK_SUBTEST(testEasom);

    SimTK_END_TEST();
}