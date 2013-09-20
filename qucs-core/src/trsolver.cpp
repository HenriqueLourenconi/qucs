/*
 * trsolver.cpp - transient solver class implementation
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2009 Stefan Jahn <stefan@lkcc.org>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * $Id$
 *
 */

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <iostream>

#include "compat.h"
#include "object.h"
#include "logging.h"
#include "complex.h"
#include "circuit.h"
#include "sweep.h"
#include "net.h"
#include "netdefs.h"
#include "analysis.h"
#include "nasolver.h"
#include "history.h"
#include "trsolver.h"
#include "transient.h"
#include "exception.h"
#include "exceptionstack.h"

#define STEPDEBUG   0 // set to zero for release
#define BREAKPOINTS 0 // exact breakpoint calculation

#define dState 0 // delta T state

// Macro for the n-th state of the solution vector history.
#define SOL(state) (solution[getIndex (state)])
#define UPD(state) (update[getIndex (state)])
#define RHS(state) (rhs[getIndex (state)])

using namespace transient;

// Constructor creates an unnamed instance of the trsolver class.
trsolver::trsolver ()
    : nasolver<nr_double_t> (), states<nr_double_t> ()
{
    swp = NULL;
    type = ANALYSIS_TRANSIENT;
    setDescription ("transient");
    for (int i = 0; i < 8; i++) solution[i] = NULL;
    for (int i = 0; i < 8; i++) update[i] = NULL;
    for (int i = 0; i < 8; i++) rhs[i] = NULL;
    tHistory = NULL;
    relaxTSR = false;
    initialDC = true;
    corrType = INTEGRATOR_UNKNOWN;
    eqns_c = new eqnsys<nr_complex_t> (); // just in case
    MA_c = NULL;
    mx_c = NULL;
    dmx_c = NULL;
    dmxsum_c = NULL;
    mz_c = NULL;
    mxprev_c = NULL;
    dmxsumprev_c = NULL;
}

// Constructor creates a named instance of the trsolver class.
trsolver::trsolver (char * n)
    : nasolver<nr_double_t> (n), states<nr_double_t> ()
{
    swp = NULL;
    type = ANALYSIS_TRANSIENT;
    setDescription ("transient");
    for (int i = 0; i < 8; i++) solution[i] = NULL;
    for (int i = 0; i < 8; i++) update[i] = NULL;
    for (int i = 0; i < 8; i++) rhs[i] = NULL;
    tHistory = NULL;
    relaxTSR = false;
    initialDC = true;
    corrType = INTEGRATOR_UNKNOWN;
    eqns_c = new eqnsys<nr_complex_t> (); // just in case
    MA_c = NULL;
    mx_c = NULL;
    dmx_c = NULL;
    dmxsum_c = NULL;
    mz_c = NULL;
    mxprev_c = NULL;
    dmxsumprev_c = NULL;
}

// Destructor deletes the trsolver class object.
trsolver::~trsolver ()
{
    if (swp) delete swp;
    for (int i = 0; i < 8; i++)
    {
        if (solution[i] != NULL)
            delete solution[i];
        if (update[i] != NULL)
            delete update[i];
	if (rhs[i] != NULL)
	    delete rhs[i];
    }

    if (MA_c != NULL)
	delete MA_c;
    if (mx_c != NULL)
	delete mx_c;
    if (dmx_c != NULL)
	delete dmx_c;
    if (dmxsum_c != NULL)
	delete dmxsum_c;
    if (mz_c != NULL)
	delete mz_c;
    if (mxprev_c != NULL)
	delete mxprev_c;
    if (dmxsumprev_c != NULL)
	delete dmxsumprev_c;
    if (eqns_c != NULL)
	delete eqns_c;

    if (tHistory) delete tHistory;
}

/* The copy constructor creates a new instance of the trsolver class
   based on the given trsolver object. */
trsolver::trsolver (trsolver & o)
    : nasolver<nr_double_t> (o), states<nr_double_t> (o)
{
    swp = o.swp ? new sweep (*o.swp) : NULL;
    for (int i = 0; i < 8; i++) solution[i] = NULL;
    for (int i = 0; i < 8; i++) update[i] = NULL;
    for (int i = 0; i < 8; i++) rhs[i] = NULL;
    tHistory = o.tHistory ? new history (*o.tHistory) : NULL;
    relaxTSR = o.relaxTSR;
    initialDC = o.initialDC;
    eqns_c = o.eqns_c;
    mxprev_c = o.mxprev_c;
    dmxsumprev_c = o.dmxsumprev_c;
}

// This function creates the time sweep if necessary.
void trsolver::initSteps (void)
{
    if (swp != NULL) delete swp;
    swp = createSweep ("time");
}

// Performs the initial DC analysis.
int trsolver::dcAnalysis (void)
{
    int error = 0;

    // First calculate a initial state using the non-linear DC analysis.
    setDescription ("initial DC");
    initDC ();
    setCalculation ((calculate_func_t) &calcDC);
    solve_pre ();
    applyNodeset ();

    fixpoint = 0;

    // Run the DC solver once.
    try_running ()
    {
        told = 1;
        error = solve_nonlinear ();
    }
    // Appropriate exception handling.
    catch_exception ()
    {
    case EXCEPTION_NO_CONVERGENCE:
        pop_exception ();
	convHelper = CONV_SourceStepping;
	//        convHelper = CONV_GMinStepping;
        logprint (LOG_ERROR, "WARNING: %s: %s analysis failed, using source"
                  " stepping fallback\n", getName (), getDescription ());
	for (int i = 0; i < x->getSize (); i++)
	    x->set (i, 0);
        applyNodeset ();
        restart ();
	//        applyNodeset ();
        error = solve_nonlinear ();
        break;
    default:
        // Otherwise return.
        estack.print ();
        error++;
        break;
    }

    // Save the DC solution.
    storeSolution ();

    // Cleanup nodal analysis solver.
    solve_post ();

    // Really failed to find initial DC solution?
    if (error)
    {
        logprint (LOG_ERROR, "ERROR: %s: %s analysis failed\n",
                  getName (), getDescription ());
    }
    return error;
}

/* This is the transient netlist solver.  It prepares the circuit list
   for each requested time and solves it then. */
int trsolver::solve (void)
{
    nr_double_t time, saveCurrent;
    int error = 0, convError = 0;
    char * solver = getPropertyString ("Solver");
    relaxTSR = !strcmp (getPropertyString ("relaxTSR"), "yes") ? true : false;
    initialDC = !strcmp (getPropertyString ("initialDC"), "yes") ? true : false;

    runs++;
    saveCurrent = current = 0;
    stepDelta = -1;
    converged = 0;
    fixpoint = 0;
    statRejected = statSteps = statIterations = statConvergence = 0;

    // Choose a solver.
    if (!strcmp (solver, "CroutLU"))
        eqnAlgo = ALGO_LU_DECOMPOSITION;
    else if (!strcmp (solver, "DoolittleLU"))
        eqnAlgo = ALGO_LU_DECOMPOSITION_DOOLITTLE;
    else if (!strcmp (solver, "HouseholderQR"))
        eqnAlgo = ALGO_QR_DECOMPOSITION;
    else if (!strcmp (solver, "HouseholderLQ"))
        eqnAlgo = ALGO_QR_DECOMPOSITION_LS;
    else if (!strcmp (solver, "GolubSVD"))
        eqnAlgo = ALGO_SV_DECOMPOSITION;

    predType = INTEGRATOR_UNKNOWN;
    corrType = INTEGRATOR_UNKNOWN;

    // Perform initial DC analysis.
    if (initialDC)
    {
        error = dcAnalysis ();
        if (error)
            return -1;
    }

    convHelper = CONV_None;
    fixpoint = 1;

    // Initialize transient analysis.
    setDescription ("transient");
    initTR ();
    setCalculation ((calculate_func_t) &calcTR);
    solve_pre ();

    // Create time sweep if necessary.
    initSteps ();
    swp->reset ();

    // Recall the DC solution.
    recallSolution ();

    // Apply the nodesets and adjust previous solutions.
    applyNodeset (false);
    fillSolution (x);

    // Tell integrators to be initialized.
    setMode (MODE_INIT);

    int running = 0;
    rejected = 0;
    delta /= 10;
    fillState (dState, delta);
    adjustOrder (1);

    current = delta;

    // Start to sweep through time.
    for (int i = 0; i < swp->getSize (); i++)
    {
        time = swp->next ();
        if (progress) logprogressbar (i, swp->getSize (), 40);

#if DEBUG && 0
        logprint (LOG_STATUS, "NOTIFY: %s: solving netlist for t = %e\n",
                  getName (), (double) time);
#endif

        do // while (saveCurrent < time), i.e. until a requested breakpoint is hit
        {
	integrate:
#if STEPDEBUG
            if (delta == deltaMin)
            {
                // the integrator step size has become smaller than the
                // specified allowed minimum, Qucs is unable to solve the circuit
                // while meeting the tolerance conditions
                logprint (LOG_ERROR,
                          "WARNING: %s: minimum delta h = %.3e at t = %.3e\n",
                          getName (), (double) delta, (double) current);
            }
#endif
            // updates the integrator coefficients, and updates the array of prev
            // 8 deltas with the new delta for this step
            updateCoefficients (delta);

            // Run predictor to get a start value for the solution vector for
            // the successive iterative corrector process
            error += predictor ();

            // restart Newton iteration
            if (rejected)
            {
                restart ();      // restart non-linear devices
                rejected = 0;
            }

            // Run corrector process with appropriate exception handling.
            // The corrector iterates through the solutions of the integration
            // process until a certain error tolerance has been reached.
            try_running () // #defined as:    do {
            {
                error += corrector ();
            }
            catch_exception () // #defined as:   } while (0); if (estack.top ()) switch (estack.top()->getCode ())
            {
            case EXCEPTION_NO_CONVERGENCE:
                pop_exception ();

                // step back from the current time value to the previous time
                if (current > 0) current -= delta;
                // Reduce step-size (by half) if failed to converge.
                delta /= 2;
                if (delta <= deltaMin)
                {
                    // but do not reduce the step size below a specified minimum
                    delta = deltaMin;
                    // instead reduce the order of the integration
                    adjustOrder (1);
                }
                // step forward to the new current time value
                if (current > 0) current += delta;

                // Update statistics.
                statRejected++;
                statConvergence++;
                rejected++; // mark the previous step size choice as rejected
                converged = 0;
                error = 0;

                // Start using damped Newton-Raphson.
                convHelper = CONV_SteepestDescent;
                convError = 2;
#if DEBUG
                logprint (LOG_ERROR, "WARNING: delta rejected at t = %.3e, h = %.3e "
                          "(no convergence)\n", (double) saveCurrent, (double) delta);
#endif
                break;
            default:
                // Otherwise return.
                estack.print ();
                error++;
                break;
            }
            // return if any errors occured other than convergence failure
            if (error) return -1;

            // if the step was rejected, the solution loop is restarted here
            if (rejected) goto integrate;

            // check whether Jacobian matrix is still non-singular
            if (!A->isFinite ())
            {
                logprint (LOG_ERROR, "ERROR: %s: Jacobian singular at t = %.3e, "
                          "aborting %s analysis\n", getName (), (double) current,
                          getDescription ());
                return -1;
            }

            // Update statistics and no more damped Newton-Raphson.
            statIterations += iterations;
            if (--convError < 0) convHelper = 0;

            // Now advance in time or not...
            if (running > 1)
            {
                adjustDelta (time);
                adjustOrder ();
            }
            else
            {
                fillStates ();
                nextStates ();
                rejected = 0;
            }

#if STEPDEBUG
	    logprint (LOG_STATUS, "DEBUG: order: %d, %d, %d\n", predOrder, corrOrder,
		      activeStates);
#endif

            saveCurrent = current;
            current += delta;
            running++;
            converged++;

            // Tell integrators to be running.
            setMode (MODE_NONE);

            // Initialize or update history.
            if (running > 1)
            {
                updateHistory (saveCurrent);
            }
            else
            {
                initHistory (saveCurrent);
            }
        }
        while (saveCurrent < time); // Hit a requested time point?

        // Save results.
#if STEPDEBUG
        logprint (LOG_STATUS, "DEBUG: save point at t = %.3e, h = %.3e\n",
                  (double) saveCurrent, (double) delta);
#endif

#if BREAKPOINTS
        saveAllResults (saveCurrent);
#else
        saveAllResults (time);
#endif
    } // for (int i = 0; i < swp->getSize (); i++)

    solve_post ();
    if (progress) logprogressclear (40);
    logprint (LOG_STATUS, "NOTIFY: %s: average time-step %g, %d rejections\n",
              getName (), (double) (saveCurrent / statSteps), statRejected);
    logprint (LOG_STATUS, "NOTIFY: %s: average NR-iterations %g, "
              "%d non-convergences\n", getName (),
              (double) statIterations / statSteps, statConvergence);

    // cleanup
    deinitTR ();
    return 0;
}

// The function initializes the history.
void trsolver::initHistory (nr_double_t t)
{
    // initialize time vector
    tHistory = new history ();
    tHistory->append (t);
    tHistory->self ();
    // initialize circuit histories
    nr_double_t age = 0.0;
    circuit * root = subnet->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
    {
        if (c->hasHistory ())
        {
            c->applyHistory (tHistory);
            saveHistory (c);
            if (c->getHistoryAge () > age)
            {
                age = c->getHistoryAge ();
            }
        }
    }
    // set maximum required age for all circuits
    tHistory->setAge (age);
}

/* The following function updates the histories for the circuits which
   requested them. */
void trsolver::updateHistory (nr_double_t t)
{
    if (t > tHistory->last ())
    {
        // update time vector
        tHistory->append (t);
        // update circuit histories
        circuit * root = subnet->getRoot ();
        for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
        {
            if (c->hasHistory ()) saveHistory (c);
        }
        tHistory->drop ();
    }
}

// Stores node voltages and branch currents in the given circuits history.
void trsolver::saveHistory (circuit * c)
{

    int N = countNodes ();
    int r, i, s = c->getSize ();

    for (i = 0; i < s; i++)
    {
        // save node voltages
        r = findAssignedNode (c, i);
        if (r < 0)
            // the node was not found, append a zero to the history
            // matching this index
            c->appendHistory (i, 0.0);
        else
            // the node was found, append the voltage value to
            // that node's history
            c->appendHistory (i, x->get (r));
    }

    for (i = 0; i < c->getVoltageSources (); i++)
    {
        // save branch currents
        r = c->getVoltageSource () + i;
        c->appendHistory (i + s, x->get (r + N));
    }

}

/* This function predicts a start value for the solution vector for
   the successive iterative corrector process. */
int trsolver::predictor (void)
{
    int error = 0;

    assert (predType == INTEGRATOR_GEAR);

    switch (predType)
    {
    case INTEGRATOR_GEAR: // explicit GEAR
        predictGear ();
        break;
    case INTEGRATOR_ADAMSBASHFORD: // ADAMS-BASHFORD
        predictBashford ();
        break;
    case INTEGRATOR_EULER: // FORWARD EULER
        predictEuler ();
        break;
    default:
        *x = *SOL (1);  // This is too a simple predictor...
        break;
    }

    saveSolution ();
    *SOL (0) = *x;

    if (corrType == INTEGRATOR_RADAU5)
    {
	int n = countNodes () + countVoltageSources ();
	for (int k = 0; k < n; k++)
	{
	    (*dmxsum)(k) = 0;
	    (*dmxsum_c)(k) = 0;
	}
	*x = *SOL (1);
	saveSolution ();

	//int n = countNodes () + countVoltageSources ();
	//for (int k = 0; k < n; k++)
	//{
	//    nr_double_t dif = (*x)(k) - (*SOL (1))(k);
	//    (*dmxsum)(k) = 0;
	//    (*dmxsum_c)(k) = 0;
	//    for (int j = 0; j < 3; j++)
	//    {
	//	(*dmxsum)(k) += real (radau_Pi[0][j]) * dif;
	//	(*dmxsum_c)(k) += radau_Pi[1][j] * dif;
	//    }
	//}
    }
    else
	*dmxsum = *x - *SOL (1);

    return error;
}

// Stores the given vector into all the solution vectors.
void trsolver::fillSolution (tvector<nr_double_t> * s)
{
    for (int i = 0; i < 8; i++)
    {
        *SOL (i) = *s;
    }
}

/* The function predicts the successive solution vector using the
   explicit Adams-Bashford integration formula. */
void trsolver::predictBashford (void)
{
    int N = countNodes ();
    int M = countVoltageSources ();
    nr_double_t xn, dd, hn;

    // go through each solution
    for (int r = 0; r < N + M; r++)
    {
        xn = predCoeff[0] * SOL(1)->get (r); // a0 coefficient
        for (int o = 1; o <= predOrder; o++)
        {
            hn = getState (dState, o);         // previous time-step
            // divided differences
            dd = (SOL(o)->get (r) - SOL(o + 1)->get (r)) / hn;
            xn += predCoeff[o] * dd;           // b0, b1, ... coefficients
        }
        x->set (r, xn);                      // save prediction
    }
}

/* The function predicts the successive solution vector using the
   explicit forward Euler integration formula.  Actually this is
   Adams-Bashford order 1. */
void trsolver::predictEuler (void)
{
    int N = countNodes ();
    int M = countVoltageSources ();
    nr_double_t xn, dd, hn;

    for (int r = 0; r < N + M; r++)
    {
      //      printf("r=%i\n", r);
        hn = getState (dState, 1);
        dd = (SOL(1)->get (r) - SOL(2)->get (r)) / hn;
//	printf("(%g-%g)/%g\n", SOL(1)->get (r), SOL(2)->get (r), hn);
//	printf("%g vs. %g\n", dd, CHFL(1)->get(r));

        xn = predCoeff[0] * SOL(1)->get (r);
        xn += predCoeff[1] * dd;
        x->set (r, xn);
    }
}

void trsolver::calcMatrices (void)
{
    switch (corrType) {
    case INTEGRATOR_GEAR:
	calcGear ();
	break;
    case INTEGRATOR_TRAPEZOIDAL:
	calcBilinear ();
	break;
    case INTEGRATOR_EULER:
	calcEuler ();
	break;
    case INTEGRATOR_ADAMSMOULTON:
	calcMoulton ();
	break;
    case INTEGRATOR_RADAU5:
	calcRadau5 ();
	break;
    default:
        nasolver::calcMatrices ();
    }
}

void trsolver::calcMAMulti (void)
{
    calculate ();
    createMatrix ();

    if (updateMatrix)
    {
	*MA = *F;
	*MA *= corrCoeff[0];
        *MA += *A;

	scaleMatrix (true);
    }
}

void trsolver::scaleMatrix (bool ldlt)
{
    int N = countNodes ();
    int M = countVoltageSources ();
    int n = N + M;

    if (L) delete (L);
    if (ldlt)
    {
	L = new tmatrix<nr_double_t> ();

	getLDLt (*L);

	// Perform the operations on MA
	for (int a = 0; a < N + M; a++)
        {
	    for (int i = a+1; i < N + M; i++)
	    {
		nr_double_t l = (*L)(i, a);
		if (l == 0) continue;
		for (int j = 0; j < N + M; j++)
		    (*MA)(i, j) -= (*MA)(a, j) * l;
	    }
	    for (int i = a+1; i < N + M; i++)
	    {
		nr_double_t l = (*L)(i, a);
		if (l == 0) continue;
		for (int j = 0; j < N + M; j++)
		    (*MA)(j, i) -= (*MA)(j, a) * l;
	    }
	}
    }
    
    if (RS != NULL) delete RS;
    RS = new tvector<nr_double_t> (n);
    if (CS != NULL) delete CS;
    CS = new tvector<nr_double_t> (n);    

    // First, attempt diagonal dominance and normalize
    for (int i = N; i < N + M; i++)
    {
    	nr_double_t d;
    
    	d = fabs ((*MA)(i, i));
    
    	nr_double_t Bmax = 0, Cmax = 0;
    
    	// Check the B part
    	for (int j = 0; j < N; j++)
    	    if (fabs ((*MA)(j, i)) > Bmax)
    		Bmax = fabs ((*MA)(j, i));
    
    	// Check the C part
    	for (int j = 0; j < N; j++)
    	    if (fabs ((*MA)(i, j)) > Cmax)
    		Cmax = fabs ((*MA)(i, j));
    
    	if (d > Bmax && d > Cmax)
    	{
    	    // Just scale the row
    	    RS->set (i, 1 / d);
    	    CS->set (i, 1);
    	}
    	else if (d <= Bmax && d <= Cmax)
    	{
    	    // Ditto
    	    RS->set (i, 1 / Cmax);
    	    CS->set (i, 1);
    	}
    	else if (Bmax > Cmax)
    	{
    	    // Scale the row, normalize the row and the column
	    if (Cmax * Bmax > d)
		RS->set (i, 1 / Cmax);
	    else
		RS->set (i, Bmax / d);
    	    CS->set (i, 1 / Bmax);
    	}
    	else
    	{
    	    // Scale the column, normalize the row
    	    CS->set (i, Cmax / d);
    	    RS->set (i, 1 / Cmax);
    	}
    }
    
    // Apply to the matrix
    for (int i = 0; i < N; i++)
    	for (int j = N; j < N + M; j++)
    	    (*MA)(i, j) *= (*CS)(j);
    for (int i = N; i < N + M; i++)
    	for (int j = 0; j < N; j++)
    	    (*MA)(i, j) *= (*RS)(i);
    for (int i = N; i < N + M; i++)
    	for (int j = N; j < N + M; j++)
    	    (*MA)(i, j) *= (*RS)(i) * (*CS)(j);

    // Scale the first N rows
    for (int i = 0; i < N; i++)
    {
	CS->set (i, 1);
	if (!ldlt)
	{
	    RS->set (i, 1);
	    continue;
	}

	nr_double_t rmax = 0;
	for (int j = 0; j < N + M; j++)
	{
	    if (fabs ((*MA)(i, j)) > rmax)
		rmax = fabs ((*MA)(i, j));
	}
	assert (rmax > 0);
	RS->set (i, 1 / rmax);
    }

    // Apply to the matrix
    if (ldlt)
	for (int i = 0; i < N; i++)
	    for (int j = 0; j < N + M; j++)
		(*MA)(i, j) *= (*RS)(i);
}

// Do an LDLt decomposition on F
void trsolver::getLDLt (tmatrix<nr_double_t> &L)
{
    int N = countNodes ();
    int M = countVoltageSources ();
    tmatrix <nr_double_t> D = *F;

    L = teye<nr_double_t> (N + M);
    for (int a = 0; a < N-1; a++)
    {
	nr_double_t div = D(a, a);
	if (div == 0)
	    continue;
	else
	    div = 1 / div;
	for (int i = a+1; i < N; i++)
	{
	    nr_double_t d = D(i, a);
	    if (d == 0) continue;
	    for (int j = a+1; j <= i; j++)
		D(i, j) -= d * D(j, a) * div;
	    L(i, a) = d * div;
	}
    }
}

void trsolver::calcEuler (void)
{
    calcMAMulti ();
    *mz = *z - *A * *mx;

    //fprintf (stderr, "A, x, z:\n");
    //A->print (1);
    //x->print (1);
    //z->print (1);
    //mz->print (1);

    *mz += corrCoeff[1] * (*F * *dmxsum);
}

void trsolver::calcBilinear (void)
{
    calcMAMulti ();
    *mz = *z - *A * *mx + *RHS(1);
    *mz += corrCoeff[1] * (*F * *dmxsum);
}

void trsolver::calcGear (void)
{
    nr_double_t dc = 0;

    calcMAMulti ();
    *mz = *z - *A * *mx;
    for (int i = corrOrder; i > 1; i--) {
	dc += corrCoeff[i];
	*mz += dc * (*F * *UPD(i-1));
    }
    dc += corrCoeff[1];
    *mz += dc * (*F * *dmxsum);

    //if (updateMatrix)
    //	MA->print(1);
    //mx->print (1);
    //mz->print (1);
}

void trsolver::calcMoulton (void)
{
    calcMAMulti ();
    *mz = *z - *A * *mx;
    *mz += corrCoeff[1] * (*F * *dmxsum);
    for (int i = 1; i < corrOrder; i++) {
	*mz -= corrCoeff[i+1] * *RHS(i);
    }
}

void trsolver::calcRadau5 (void)
{
    //const nr_double_t radau_A[3][3] =
    //	{{(88-7*sqrt(6))/360, (296-169*sqrt(6))/1800, (-2+3*sqrt(6))/225},
    //	 {(296+169*sqrt(6))/1800, (88+7*sqrt(6))/360, (-2-3*sqrt(6))/225},
    //	 {(16-sqrt(6))/36, (16+sqrt(6))/36, 1./9}};
    //const nr_double_t radau_b[3] =
    //	{(16-sqrt(6))/36, (16+sqrt(6))/36, 1./9};

    nr_double_t newCurrent = current;
    nr_double_t saveCurrent = current - delta;
    int n = countNodes () + countVoltageSources ();
    nr_double_t dr = 1. / delta;

    if (updateMatrix)
    {
	current = saveCurrent;
	// Is this necessary?
	for (int k = 0; k < n; k++)
	    x->set(k, SOL (1)->get(k));
	saveSolution ();
	calculate ();
	createMatrix ();

	// Compute the real and complex Jacobians for the RadauIIA method
	*MA = radau_A_1 * *A;
	*MA += dr * *F;

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                (*MA_c)(i, j) = radau_A_2 * (*A)(i, j) + dr * (*F)(i, j);

	scaleMatrix (false);
	for (int i = 0; i < n; i++)
	    for (int j = 0; j < n; j++)
		(*MA_c)(i, j) *= (*RS)(i) * (*CS)(j);
	//logprint (LOG_STATUS, "DEBUG: condition (complex): %g\n", condition (*MA_c));
    }
    
    //tmatrix<nr_complex_t> checkident = tmatrix<nr_complex_t> (3);
    //
    //for (int i = 0; i < 3; i++)
    //	for (int j = 0; j < 3; j++)
    //	    for (int k = 0; k < 3; k++)
    //		checkident(i, j) += radau_Pi[i][k] * radau_P[k][j];
    //checkident.print ();

    // Transform mx back
    tvector<nr_double_t> mx_r[3];
    for (int i = 0; i < 3; i++)
    {
	mx_r[i] = tvector<nr_double_t> (n);
	for (int k = 0; k < n; k++)
	{
	    mx_r[i](k) = real (radau_P[i][0]) * (*mx)(k)
		+ real (radau_P[i][1] * (*mx_c)(k))
		+ real (radau_P[i][2] * conj ((*mx_c)(k)));
	    //fprintf (stderr, "real: %g\n", mx_r[i](k));
	    //fprintf (stderr, "imaginary: %g\n",
	    //	     imag (radau_P[i][1] * (*mx_c)(k))
	    //	     + imag (radau_P[i][2] * conj ((*mx_c)(k))));
	}
    }

    // Evaluate the circuit
    tvector<nr_double_t> fsave[3];
    for (int i = 0; i < 3; i++)
    {
	current = saveCurrent + radau_c[i] * delta;
	*x = mx_r[i];

	saveSolution ();
	calculate ();
	createMatrix ();

	fsave[i] = *z - *A * *x;
	//fprintf (stderr, "A, x, z:\n");
	//A->print (1);
	//x->print (1);
	//z->print (1);
	//fsave[i].print (1);
    }

    // Combine
    for (int k = 0; k < n; k++)
    {
	nr_double_t mzk = 0;
	nr_complex_t mzk_c = 0;

	// Evaluate capacitors
	for (int l = 0; l < n; l++)
	{
	    mzk -= F->get(k, l) * dmxsum->get(l) * dr;
	    mzk_c -= F->get(k, l) * dmxsum_c->get(l) * dr;
	}

	// Add evaluated resistors etc.
	for (int j = 0; j < 3; j++)
	{
	    mzk += radau_A_1 * real (radau_Pi[0][j]) * fsave[j](k);
	    mzk_c += radau_A_2 * radau_Pi[1][j] * fsave[j](k);
	}	    

	mz->set (k, mzk);
	mz_c->set (k, mzk_c);
    }

    //mz->print (1);

    //fprintf (stderr, "complex system:\n");
    //if (updateMatrix)
    //	MA_c->print(0);
    //mz_c->print (0);

    //if (updateMatrix)
    //	MA->print(1);
    //mx_r[2].print (1);
    //fsave[2].print (1);

    // We should already be there, but still:
    current = newCurrent;
}

void trsolver::setInitX (void)
{
    int n = countNodes () + countVoltageSources ();

    switch (corrType)
    {
    case INTEGRATOR_RADAU5:
	for (int k = 0; k < n; k++)
	{
	    (*mx)(k) = 0;
	    (*mx_c)(k) = 0;
	    for (int j = 0; j < 3; j++)
	    {
		(*mx)(k) += real (radau_Pi[0][j]) * (*x)(k);
		(*mx_c)(k) += radau_Pi[1][j] * (*x)(k);
	    }
	}
	break;
    default:
	nasolver::setInitX ();
    }
}

void trsolver::setInit (void)
{
    nasolver::setInit ();
    if (corrType == INTEGRATOR_RADAU5)
    {
	if (mxprev_c != NULL) delete mxprev_c;
	mxprev_c = new tvector<nr_complex_t> (*mx_c);
	if (dmxsumprev_c != NULL) delete dmxsumprev_c;
	dmxsumprev_c = new tvector<nr_complex_t> (*dmxsum_c);
    }
}

void trsolver::extractSol (void)
{
    int n = countNodes () + countVoltageSources ();

    switch (corrType)
    {
    case INTEGRATOR_RADAU5:
	for (int k = 0; k < n; k++)
	{
	    (*x)(k) = real (radau_P[2][0]) * (*mx)(k)
		+ real (radau_P[2][1] * (*mx_c)(k))
		+ real (radau_P[2][2] * conj ((*mx_c)(k)));
	    (*dx)(k) = real (radau_P[2][0]) * (*dmxsum)(k)
		+ real (radau_P[2][1] * (*dmxsum_c)(k))
		+ real (radau_P[2][2] * conj ((*dmxsum_c)(k)));
	}
	break;
    default:
	nasolver::extractSol ();
    }

    //x->print (1);
    //dx->print (1);
    //if (corrType != INTEGRATOR_UNKNOWN)
    //{
    //	tvector<nr_double_t> dx1 = *x - *SOL (1);
    //	dx1.print (1);
    //}
}

int trsolver::getSysSize (void)
{
    int n = countVoltageSources () + countNodes ();

    return n;
}

void trsolver::solve_pre (void)
{
    nasolver::solve_pre ();

    if (corrType == INTEGRATOR_RADAU5)
    {
	int n = getSysSize ();

	if (MA_c != NULL) delete MA_c;
	MA_c = new tmatrix<nr_complex_t> (n);
	if (mx_c != NULL) delete mx_c;
	mx_c = new tvector<nr_complex_t> (n);
	if (dmx_c != NULL) delete dmx_c;
	dmx_c = new tvector<nr_complex_t> (n);
	if (dmxsum_c != NULL) delete dmxsum_c;
	dmxsum_c = new tvector<nr_complex_t> (n);
	if (mz_c != NULL) delete mz_c;
	mz_c = new tvector<nr_complex_t> (n);
    }
}

void trsolver::update_mx ()
{
    nasolver::update_mx ();

    if (corrType == INTEGRATOR_RADAU5)
    {
	if (mxprev_c != NULL)
	    *mx_c = *mxprev_c + *dmx_c;
	else
	    *mx_c += *dmx_c;

	if (dmxsumprev_c != NULL)
	    *dmxsum_c = *dmxsumprev_c + *dmx_c;
	else
	    *dmxsum_c += *dmx_c;

	//fprintf (stderr, "update:\n");
	//dmx_c->print ();
	//mx_c->print ();
    }
}

void trsolver::solveEquation ()
{
    nasolver::solveEquation ();

    if (corrType == INTEGRATOR_RADAU5)
    {
	int n = getSysSize ();

	if (RS)
	    for (int k = 0; k < n; k++)
		(*mz_c)(k) *= (*RS)(k);

	eqns_c->setAlgo (eqnAlgo);
	eqns_c->passEquationSys (updateMatrix ? MA_c : NULL, dmx_c, mz_c);
	eqns_c->solve ();

	if (CS)
	    for (int k = 0; k < n; k++)
	    (*dmx_c)(k) *= (*CS)(k);
    }
}

int trsolver::checkConvergence (void)
{
    int retval = nasolver::checkConvergence ();

    if (retval <= 0 || corrType != INTEGRATOR_RADAU5)
	return retval;

    int N = countNodes ();
    int M = countVoltageSources ();
    nr_double_t v_abs, v_rel, i_abs, i_rel;
    int r;

#if STEPDEBUG
    logprint (LOG_STATUS, "DEBUG: covariant norm (complex): %g\n", norm (*mx_c));
#endif

    for (r = 0; r < N; r++)
    {
        v_abs = fabs (dmx_c->get (r));
        v_rel = fabs (mx_c->get (r));
        if (told * v_abs >= vntol + reltol * v_rel) goto noconv;
    }

    for (r = 0; r < M; r++)
    {
        i_abs = fabs (dmx_c->get (r + N));
        i_rel = fabs (mx_c->get (r + N));
        if (told * i_abs >= abstol + reltol * i_rel) goto noconv;
    }
    return 1;
    
 noconv:
    return 0;
}

void trsolver::savePreviousIteration (void)
{
    nasolver::savePreviousIteration ();
    if (corrType == INTEGRATOR_RADAU5)
    {
	*mxprev_c = *mx_c;

	//if (dmxprev_c != NULL)
	//	*dmxprev_c = *dmx_c;
	//else
	//	dmxprev_c = new tvector<nr_complex_t> (*dmx_c);

	*dmxsumprev_c = *dmxsum_c;
    }
}

/* The function predicts the successive solution vector using the
   explicit Gear integration formula. */
void trsolver::predictGear (void)
{
    int N = countNodes ();
    int M = countVoltageSources ();
    nr_double_t xn;

    // go through each solution
    for (int r = 0; r < N + M; r++)
    {
        xn = 0;
        for (int o = 0; o <= predOrder; o++)
        {
            // a0, a1, ... coefficients
            xn += predCoeff[o] * SOL(o + 1)->get (r);
        }
        x->set (r, xn); // save prediction
    }
}

/* The function iterates through the solutions of the integration
   process until a certain error tolerance has been reached. */
int trsolver::corrector (void)
{
    int error = 0;
    told = delta;
    error += solve_nonlinear ();
    return error;
}

// The function advances one more time-step.
void trsolver::nextStates (void)
{
    circuit * root = subnet->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
    {
        // for each circuit get the next state
        c->nextState ();
    }

    *SOL (0) = *x;              // save current solution
    *UPD (0) = *dx;             // save update
    *RHS (0) = *z - *A * *x;    // save resistive current

    nextState ();
    statSteps++;
}

/* This function stores the current state of each circuit into all
   other states as well.  It is useful for higher order integration
   methods in order to initialize the states after the initial
   transient solution. */
void trsolver::fillStates (void)
{
    circuit * root = subnet->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
    {
        for (int s = 0; s < c->getStates (); s++)
            c->fillState (s, c->getState (s));
    }
}

// The function modifies the circuit lists integrator mode.
void trsolver::setMode (int state)
{
    circuit * root = subnet->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
        c->setMode (state);
}

// The function passes the time delta array to the circuit list.
void trsolver::setDelta (void)
{
    circuit * root = subnet->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
        c->setDelta (deltas);
}

/* This function tries to adapt the current time-step according to the
   global truncation error. */
void trsolver::adjustDelta (nr_double_t t)
{
    deltaOld = delta;
    delta = checkDelta ();
    if (delta > deltaMax) delta = deltaMax;
    if (delta < deltaMin) delta = deltaMin;

    // delta correction in order to hit exact breakpoint
    int good = 0;
    if (!relaxTSR)   // relaxed step raster?
    {
        if (!statConvergence || converged > 64)   /* Is this a good guess? */
        {
            // check next breakpoint
            if (stepDelta > 0.0)
            {
                // restore last valid delta
                delta = stepDelta;
                stepDelta = -1.0;
            }
            else
            {
                if (delta > (t - current) && t > current)
                {
                    // save last valid delta and set exact step
                    stepDelta = deltaOld;
                    delta = t - current;
                    good = 1;
                }
                else
                {
                    stepDelta = -1.0;
                }
            }
            if (delta > deltaMax) delta = deltaMax;
            if (delta < deltaMin) delta = deltaMin;
        }
    }

    // usual delta correction
    if (delta > 0.9 * deltaOld || good)   // accept current delta
    {
        nextStates ();
        rejected = 0;
#if STEPDEBUG
        logprint (LOG_STATUS,
                  "DEBUG: delta accepted at t = %.3e, h = %.3e\n",
                  (double) current, (double) delta);
#endif
    }
    else if (deltaOld > delta)   // reject current delta
    {
        rejected++;
        statRejected++;
#if STEPDEBUG
        logprint (LOG_STATUS,
                  "DEBUG: delta rejected at t = %.3e, h = %.3e\n",
                  (double) current, (double) delta);
#endif
        if (current > 0) current -= deltaOld;
    }
    else
    {
        nextStates ();
        rejected = 0;
    }
}

/* The function can be used to increase the current order of the
   integration method or to reduce it. */
void trsolver::adjustOrder (int reduce)
{
    if (corrOrder < corrMaxOrder || predOrder < predMaxOrder)
    {
        if (reduce)
        {
            corrOrder = 1;
        }
        else if (!rejected && corrOrder < corrMaxOrder)
        {
            corrOrder++;
        }

        // adjust type and order of corrector and predictor
        corrType = correctorType (CMethod, corrOrder);
	predOrder = corrOrder;
        predType = predictorType (corrType, predOrder);

        // apply new corrector method and order to each circuit
        circuit * root = subnet->getRoot ();
        for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
        {
            c->setOrder (corrOrder);
            setIntegrationMethod (c, corrType);
        }
    }
}

/* Goes through the list of circuit objects and runs its calcDC()
   function. */
void trsolver::calcDC (trsolver * self)
{
    circuit * root = self->getNet()->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
    {
        c->calcDC ();
    }
}

/* Goes through the list of circuit objects and runs its calcTR()
   function. */
void trsolver::calcTR (trsolver * self)
{
    circuit * root = self->getNet()->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
    {
        c->calcTR (self->current);
    }
}

/* Goes through the list of non-linear circuit objects and runs its
   restartDC() function. */
void trsolver::restart (void)
{
    circuit * root = subnet->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
    {
        if (c->isNonLinear ()) c->restartDC ();
    }
}

/* Goes through the list of circuit objects and runs its initDC()
   function. */
void trsolver::initDC (void)
{
    circuit * root = subnet->getRoot ();
    for (circuit * c = root; c != NULL; c = (circuit *) c->getNext ())
    {
        c->initDC ();
    }
}

/* Goes through the list of circuit objects and runs its initTR()
   function. */
void trsolver::initTR (void)
{
    char * IMethod = getPropertyString ("IntegrationMethod");
    nr_double_t start = getPropertyDouble ("Start");
    nr_double_t stop = getPropertyDouble ("Stop");
    nr_double_t points = getPropertyDouble ("Points");

    activeStates = 1;

    // fetch corrector integration method and determine predicor method
    corrMaxOrder = getPropertyInteger ("Order");
    corrType = CMethod = correctorType (IMethod, corrMaxOrder);
    predMaxOrder = corrMaxOrder;
    predOrder = min (predMaxOrder, max (1, activeStates - 1));
    predType = PMethod = predictorType (CMethod, predOrder);
    corrOrder = corrMaxOrder;

    // initialize step values
    delta = getPropertyDouble ("InitialStep");
    deltaMin = getPropertyDouble ("MinStep");
    deltaMax = getPropertyDouble ("MaxStep");
    if (deltaMax == 0.0)
        deltaMax = MIN ((stop - start) / (points - 1), stop / 200);
    if (deltaMin == 0.0)
        deltaMin = NR_TINY * 10 * deltaMax;
    if (delta == 0.0)
        delta = MIN (stop / 200, deltaMax) / 10;
    if (delta < deltaMin) delta = deltaMin;
    if (delta > deltaMax) delta = deltaMax;

    // initialize step history
    setStates (1);
    initStates ();
    // initialise the history of states, setting them all to 'delta'
    fillState (dState, delta);

    // copy the initialised states to the 'deltas' array
    saveState (dState, deltas);
    // copy the deltas to all the circuits
    setDelta ();
    // set the initial corrector and predictor coefficients
    calcCorrectorCoeff (corrType, corrOrder, corrCoeff, deltas);
    calcPredictorCoeff (predType, predOrder, predCoeff, deltas);

    // initialize history of solution vectors (solutions)
    for (int i = 0; i < 8; i++)
    {
        // solution contains the last sets of node voltages and branch
        // currents at each of the last 8 'deltas'.
        solution[i] = new tvector<nr_double_t>;
        update[i] = new tvector<nr_double_t>;
        rhs[i] = new tvector<nr_double_t>;
    }

    // tell circuits about the transient analysis
    circuit *c, * root = subnet->getRoot ();
    for (c = root; c != NULL; c = (circuit *) c->getNext ())
        initCircuitTR (c);
    // also initialize created circuits
    for (c = root; c != NULL; c = (circuit *) c->getPrev ())
        initCircuitTR (c);
}

// This function cleans up some memory used by the transient analysis.
void trsolver::deinitTR (void)
{
    // cleanup solutions
    for (int i = 0; i < 8; i++)
    {
        delete solution[i];
        solution[i] = NULL;
        delete update[i];
        update[i] = NULL;
        delete rhs[i];
        rhs[i] = NULL;
    }
    // cleanup history
    if (tHistory)
    {
        delete tHistory;
        tHistory = NULL;
    }
}

// The function initialize a single circuit.
void trsolver::initCircuitTR (circuit * c)
{
    c->initTR ();
    c->initStates ();
    c->setCoefficients (corrCoeff);
    c->setOrder (corrOrder);
    setIntegrationMethod (c, corrType);
}

/* This function saves the results of a single solve() functionality
   (for the given timestamp) into the output dataset. */
void trsolver::saveAllResults (nr_double_t time)
{
    ::vector * t;
    // add current frequency to the dependency of the output dataset
    if ((t = data->findDependency ("time")) == NULL)
    {
      t = new ::vector ("time");
        data->addDependency (t);
    }
    if (runs == 1) t->add (time);
    saveResults ("Vt", "It", 0, t);
}

/* This function is meant to adapt the current time-step the transient
   analysis advanced.  For the computation of the new time-step the
   truncation error depending on the integration method is used. */
nr_double_t trsolver::checkDelta (void)
{
    nr_double_t LTEreltol = getPropertyDouble ("LTEreltol");
    nr_double_t LTEabstol = getPropertyDouble ("LTEabstol");
    nr_double_t LTEfactor = getPropertyDouble ("LTEfactor");
    nr_double_t dif, rel, tol, lte, q, n = NR_MAX;
    int N = countNodes ();
    int M = countVoltageSources ();

    // cec = corrector error constant
    nr_double_t cec = getCorrectorError (corrType, corrOrder);
    // pec = predictor error constant
    nr_double_t pec = getPredictorError (predType, predOrder);

    // go through each solution
    for (int r = 0; r < N + M; r++)
    {

        // skip real voltage sources
        if (r >= N)
        {
            if (findVoltageSource(r - N)->isVSource ())
                continue;
        }

        dif = x->get (r) - SOL(0)->get (r);
        if (finite (dif) && dif != 0)
        {
            // use Milne' estimate for the local truncation error
            rel = MAX (fabs (x->get (r)), fabs (SOL(0)->get (r)));
            tol = LTEreltol * rel + LTEabstol;
            lte = LTEfactor * (cec / (pec - cec)) * dif;
            q =  delta * exp (log (fabs (tol / lte)) / (corrOrder + 1));
            n = MIN (n, q);
        }
    }
#if STEPDEBUG
    logprint (LOG_STATUS, "DEBUG: delta according to local truncation "
              "error h = %.3e\n", (double) n);
#endif
    delta = MIN ((n > 1.9 * delta) ? 2 * delta : delta, n);
    return delta;
}

// The function updates the integration coefficients.
void trsolver::updateCoefficients (nr_double_t delta)
{
    setState (dState, delta);
    saveState (dState, deltas);
    calcCorrectorCoeff (corrType, corrOrder, corrCoeff, deltas);
    calcPredictorCoeff (predType, predOrder, predCoeff, deltas);
}

// properties
PROP_REQ [] =
{
    {
        "Type", PROP_STR, { PROP_NO_VAL, "lin" },
        PROP_RNG_STR2 ("lin", "log")
    },
    { "Start", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
    { "Stop", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
    { "Points", PROP_INT, { 10, PROP_NO_STR }, PROP_MIN_VAL (2) },
    PROP_NO_PROP
};
PROP_OPT [] =
{
    {
        "IntegrationMethod", PROP_STR, { PROP_NO_VAL, "Trapezoidal" },
        PROP_RNG_STR5 ("Euler", "Trapezoidal", "Gear", "AdamsMoulton", "Radau5")
    },
    { "Order", PROP_INT, { 2, PROP_NO_STR }, PROP_RNGII (1, 6) },
    { "InitialStep", PROP_REAL, { 1e-9, PROP_NO_STR }, PROP_POS_RANGE },
    { "MinStep", PROP_REAL, { 1e-16, PROP_NO_STR }, PROP_POS_RANGE },
    { "MaxStep", PROP_REAL, { 0, PROP_NO_STR }, PROP_POS_RANGE },
    { "MaxIter", PROP_INT, { 150, PROP_NO_STR }, PROP_RNGII (2, 10000) },
    { "abstol", PROP_REAL, { 1e-12, PROP_NO_STR }, PROP_RNG_X01I },
    { "vntol", PROP_REAL, { 1e-6, PROP_NO_STR }, PROP_RNG_X01I },
    { "reltol", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_RNG_X01I },
    { "LTEabstol", PROP_REAL, { 1e-6, PROP_NO_STR }, PROP_RNG_X01I },
    { "LTEreltol", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_RNG_X01I },
    { "LTEfactor", PROP_REAL, { 1, PROP_NO_STR }, PROP_RNGII (1, 16) },
    { "Temp", PROP_REAL, { 26.85, PROP_NO_STR }, PROP_MIN_VAL (K) },
    { "Solver", PROP_STR, { PROP_NO_VAL, "CroutLU" }, PROP_RNG_SOL },
    { "relaxTSR", PROP_STR, { PROP_NO_VAL, "no" }, PROP_RNG_YESNO },
    { "initialDC", PROP_STR, { PROP_NO_VAL, "yes" }, PROP_RNG_YESNO },
    PROP_NO_PROP
};
struct define_t trsolver::anadef =
    { "TR", 0, PROP_ACTION, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };
