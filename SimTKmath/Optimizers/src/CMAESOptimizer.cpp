/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
 * Contributors: Michael Sherman, Nikolaus Hansen, Jack Wang                  *
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

#include "CMAESOptimizer.h"

#if SimTK_SIMMATH_MPI
    #include <mpi.h>
    #define SimTK_COMM_WORLD MPI_COMM_WORLD
    #define SimTK_CMAES_USE_MPI(use_mpi) (SimTK_SIMMATH_MPI && use_mpi)
#endif

#include <bitset>

namespace SimTK {

#define SimTK_CMAES_PRINT(diag, cmds) \
do { \
if (std::bitset<2>(diag).test(0)) { cmds; } \
} \
while(false)

#define SimTK_CMAES_FILE(diag, cmds) \
do { \
if (std::bitset<2>(diag).test(1)) { cmds; } \
} \
while(false)

CMAESOptimizer::CMAESOptimizer(const OptimizerSystem& sys) : OptimizerRep(sys)
{
    SimTK_VALUECHECK_ALWAYS(2, sys.getNumParameters(), INT_MAX, "nParameters",
            "CMAESOptimizer");
}

Optimizer::OptimizerRep* CMAESOptimizer::clone() const {
    return new CMAESOptimizer(*this);
}

Real CMAESOptimizer::optimize(SimTK::Vector& results)
{
    // Initialize parallelism, if requested.
    std::string parallel;
    bool useMPI = false;
    if (getAdvancedStrOption("parallel", parallel)) {
        if (parallel == "mpi") {
            useMPI = true;
        }
    }
 
//    if (SimTK_CMAES_USE_MPI(use_mpi)) {
//    // TODO don't use comM_world.
//    // TODO use MPI_Scatter.
//    // TODO test on an actual cluster first.
//    // TODO test on windows.
//    // TODO example that uses MPI.
//    // TODO http://www.lam-mpi.org/tutorials/one-step/ezstart.php
//    // TODO libraries need private communicators (COMM_CREATE_GROUP).
//    // TODO http://stackoverflow.com/questions/13867809/how-are-mpi-scatter-and-mpi-gather-used-from-c
//        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//        printf("My rank: %d.\n", myRank);
//        MPI_Finalized(&finalized);
//        if (!finalized) MPI_Finalize();
//        initialized = false;
//        finalized = false;
//    }

    // If using MPI, figure out whether we're the master or a worker.
    #if SimTK_CMAES_USE_MPI
        if (useMPI) {
            
            // Initialize MPI if it hasn't been initialized yet.
            int initialized;
            MPI_Initialized(&initialized);
            if (!initialized) {
                MPI_Init(NULL, NULL);
            }

            // Are we a worker?
            int myRank = 0;
            MPI_Comm_rank(SimTK_COMM_WORLD, &myRank);
            if (myRank == 0) {
                return master(results);
            }
            else {
                mpi_worker();
                return 0;
            }
        }
    #endif
    
    return master(results);
}

Real CMAESOptimizer::master(SimTK::Vector& results)
{ 

    const OptimizerSystem& sys = getOptimizerSystem();
    int n = sys.getNumParameters();

    // Initialize objective function value and cmaes data structure.
    // =============================================================
    Real f; 
    cmaes_t evo;
    
    // Initialize parallelism, if requested.
    std::string parallel;
    std::auto_ptr<ParallelExecutor> executor;
    bool use_mpi = false;
    if (getAdvancedStrOption("parallel", parallel)) {

        // Number of parallel processes/threads.
        int number_of_threads = ParallelExecutor::getNumProcessors();
        getAdvancedIntOption("number_of_threads", number_of_threads);

        // Multithreading.
        if (parallel == "multithreading") {
            executor.reset(new ParallelExecutor(number_of_threads));
        }

    }

    // Check that the initial point is feasible.
    // =========================================
    checkInitialPointIsFeasible(results);
    
    // Initialize cmaes.
    // =================
    double* funvals = init(evo, results);
    SimTK_CMAES_PRINT(diagnosticsLevel, printf("%s\n", cmaes_SayHello(&evo)));
    
    // Optimize.
    // =========
    while (!cmaes_TestForTermination(&evo)) {

        // Sample a population.
        // ====================
        double*const* pop = cmaes_SamplePopulation(&evo);

        // Resample to keep population within limits.
        // ==========================================
        resampleToObeyLimits(evo, pop);

        // Evaluate the objective function on the samples.
        // ===============================================
        evaluateObjectiveFunctionOnPopulation(evo, pop, funvals, executor);
        
        // Update the distribution (mean, covariance, etc.).
        // =================================================
        cmaes_UpdateDistribution(&evo, funvals);
    }

    // Wrap up.
    // ========
    SimTK_CMAES_PRINT(diagnosticsLevel,
            printf("Stop:\n%s\n", cmaes_TestForTermination(&evo)));

    // Update results and objective function value.
    const double* xbestever = cmaes_GetPtr(&evo, "xbestever");
    for (int i = 0; i < n; i++) {
        results[i] = xbestever[i]; 
    }
    f = cmaes_Get(&evo, "fbestever");

    SimTK_CMAES_FILE(diagnosticsLevel,
            cmaes_WriteToFile(&evo, "all", "allcmaes.dat"));

    // Free memory.
    cmaes_exit(&evo);
    
    return f;  
}

void CMAESOptimizer::checkInitialPointIsFeasible(const Vector& x) const {

    const OptimizerSystem& sys = getOptimizerSystem();

    if( sys.getHasLimits() ) {
        Real *lower, *upper;
        sys.getParameterLimits( &lower, &upper );
        for (int i = 0; i < sys.getNumParameters(); i++) {
            SimTK_APIARGCHECK4_ALWAYS(
                    lower[i] <= x[i] && x[i] <= upper[i],
                    "CMAESOptimizer", "checkInitialPointIsFeasible",
                    "Initial guess results[%d] = %f "
                    "is not within limits [%f, %f].",
                    i, x[i], lower[i], upper[i]);
        }
    }
}

double* CMAESOptimizer::init(cmaes_t& evo, SimTK::Vector& results) const
{
    const OptimizerSystem& sys = getOptimizerSystem();
    int n = sys.getNumParameters();

    // Prepare to call cmaes_init_para.
    // ================================

    // lambda
    // ------
    int lambda = 0;
    getAdvancedIntOption("lambda", lambda);
    
    // sigma
    // -----
    double sigma = 0;
    double* stddev = NULL;
    Vector sigmaArray;
    if (getAdvancedRealOption("sigma", sigma)) {
        sigmaArray.resize(n);
        for (int i = 0; i < n; i++) {
            sigmaArray[i] = sigma;
        }
        stddev = &sigmaArray[0];
    }

    // seed
    // ----
    int seed = 0;
    if (getAdvancedIntOption("seed", seed)) {
        SimTK_VALUECHECK_NONNEG_ALWAYS(seed, "seed",
                "CMAESOptimizer::init");
    }

    // input parameter filename
    // ------------------------
    std::string input_parameter_filename = "none";
    SimTK_CMAES_FILE(diagnosticsLevel,
            input_parameter_filename = "writeonly";);

    // Call cmaes_init_para.
    // =====================
    // Here, we specify the subset of options that can be passed to
    // cmaes_init_para.
    cmaes_init_para(&evo,
            n,                 // dimension
            &results[0],       // xstart
            stddev,            // stddev
            seed,              // seed
            lambda,            // lambda
            input_parameter_filename.c_str() // input_parameter_filename
            ); 

    // Set settings that are usually read in from cmaes_initials.par.
    // ==============================================================
    process_readpara_settings(evo);

    // Once we've updated settings in readpara_t, finalize the initialization.
    return cmaes_init_final(&evo);
}

void CMAESOptimizer::process_readpara_settings(cmaes_t& evo) const
{
    // Termination criteria
    // ====================

    // stopMaxIter
    // -----------
    // maxIterations is a protected member variable of OptimizerRep.
    evo.sp.stopMaxIter = maxIterations;

    // stopTolFun
    // ----------
    // convergenceTolerance is a protected member variable of OptimizerRep.
    evo.sp.stopTolFun = convergenceTolerance;

    // stopMaxFunEvals
    // ---------------
    int stopMaxFunEvals;
    if (getAdvancedIntOption("stopMaxFunEvals", stopMaxFunEvals)) {
        SimTK_VALUECHECK_NONNEG_ALWAYS(stopMaxFunEvals, "stopMaxFunEvals",
                "CMAESOptimizer::process_readpara_settings");
        evo.sp.stopMaxFunEvals = stopMaxFunEvals;
    }

    // stopFitness
    // -----------
    double stopFitness;
    if (getAdvancedRealOption("stopFitness", stopFitness)) {
        evo.sp.stStopFitness.flg = 1;
        evo.sp.stStopFitness.val = stopFitness;
    }

    // stopTolFunHist
    // --------------
    double stopTolFunHist;
    if (getAdvancedRealOption("stopTolFunHist", stopTolFunHist)) {
        evo.sp.stopTolFunHist = stopTolFunHist;
    }

    // stopTolX
    // --------
    double stopTolX;
    if (getAdvancedRealOption("stopTolX", stopTolX)) {
        evo.sp.stopTolX = stopTolX;
    }

    // stopTolXFactor
    // --------------
    double stopTolUpXFactor;
    if (getAdvancedRealOption("stopTolUpXFactor", stopTolUpXFactor)) {
        evo.sp.stopTolUpXFactor = stopTolUpXFactor;
    }

    // maxtime
    // =======
    double maxtime;
    if (getAdvancedRealOption("maxTimeFractionForEigendecomposition", maxtime))
    {
        evo.sp.updateCmode.maxtime = maxtime;
    }
}

void CMAESOptimizer::resampleToObeyLimits(cmaes_t& evo, double*const* pop)
{
    const OptimizerSystem& sys = getOptimizerSystem();
    if( sys.getHasLimits() ) {

        Real *lower, *upper;
        sys.getParameterLimits( &lower, &upper );

        for (int i = 0; i < cmaes_Get(&evo, "lambda"); i++) {
            bool feasible = false; 
            while (!feasible) {
                feasible = true; 
                for (int j = 0; j < sys.getNumParameters(); j++) {
                    if (pop[i][j] < lower[j] || pop[i][j] > upper[j]) {
                        feasible = false; 
                        pop = cmaes_ReSampleSingle(&evo, i); 
                        break; 
                    }
                }
            }
        }
    }
}

void CMAESOptimizer::evaluateObjectiveFunctionOnPopulation(
        cmaes_t& evo, double*const* pop, double* funvals,
        const std::auto_ptr<ParallelExecutor>& executor)
{
    const OptimizerSystem& sys = getOptimizerSystem();

    // Execute in parallel.
    if (executor.get()) {
        Task task(*this, sys.getNumParameters(), pop, funvals);
        executor->execute(task, cmaes_Get(&evo, "lambda"));
    }
    // Execute normally.
    else {
        for (int i = 0; i < cmaes_Get(&evo, "lambda"); i++) {
            objectiveFuncWrapper(sys.getNumParameters(),
                    pop[i], true, &funvals[i], this);
        }
    }
}

void CMAESOptimizer::slave() {
    #if SimTK_CMAES_USE_MPI
        // To get the die tag.
        MPI_Status status;

        while (true) {
            MPI_Recv
        }

    #endif
}

#undef SimTK_CMAES_PRINT
#undef SimTK_CMAES_FILE
#undef SimTK_CMAES_USE_MPI

} // namespace SimTK
