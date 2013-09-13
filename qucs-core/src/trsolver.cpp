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
            if (rejected) continue;

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
	for (int i = 0; i < 3; i++)
	    for (int k = 0; k < n; k++)
		dmxsum->set (i*n + k, x->get (k) - SOL (1)->get (k));
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

	//fprintf (stderr, "performing SVD...\n");
	//eqnsys<nr_double_t> F_eqns;
	////tmatrix<nr_double_t> U = *MA;
	//tmatrix<nr_double_t> U = *F;
	//if (Ut) delete (Ut);
	//if (Vt) delete (Vt);
	//Vt = new tmatrix<nr_double_t> ();
	//F_eqns.get_svd (&U, NULL, Vt);
	//
	//Ut = new tmatrix<nr_double_t> (U);
	//Ut->transpose ();
	//tmatrix<nr_double_t> V = *Vt;
	//V.transpose ();
	//MA->print (true);
	//*MA = *Ut * *MA * V;
	//F->print (true);
	//(*Ut * *F * V).print (true);
	//normalizeRows ();

	
	//fprintf (stderr, "normalizing...\n");
	//MA->print (true);
	scaleMatrix ();
    }
}

//void trsolver::normalizeRows (void)
//{
//    int n = getSysSize ();
//    if (Nv != NULL) delete Nv;
//    Nv = new tvector<nr_double_t> (n);
//
//    for (int i = 0; i < n; i++)
//    {
//	nr_double_t gmax = 0.0;
//	for (int j = 0; j < n; j++)
//	{
//	    nr_double_t g = fabs (MA->get (i, j));
//	    if (g > gmax)
//		gmax = g;
//	}
//	assert (gmax > 0.0);
//	nr_double_t mult = 1.0/gmax;
//	Nv->set (i, mult);
//	for (int j = 0; j < n; j++)
//	{
//	    nr_double_t g = MA->get (i, j);
//	    MA->set (i, j, g * mult);
//	}
//    }
//}

void trsolver::scaleMatrix (void)
{
    int N = countNodes ();
    int M = countVoltageSources ();
    int n = N + M;

    if (L) delete (L);
    L = new tmatrix<nr_double_t> ();

    //fprintf (stderr, "doing LDLt\n");
    getLDLt (*L);
    //tmatrix<nr_double_t> Lt = *L;
    //Lt.transpose ();
    //tmatrix<nr_double_t> MA2 = *MA;
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
    
    //fprintf (stderr, "scaling...\n");
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
	//RS->set (i, 1);
	//continue;

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
    for (int i = 0; i < N; i++)
    	for (int j = 0; j < N + M; j++)
    	    (*MA)(i, j) *= (*RS)(i);

    //MA->print (1);
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
    const nr_double_t radau_A[3][3] =
	{{(88-7*sqrt(6))/360, (296-169*sqrt(6))/1800, (-2+3*sqrt(6))/225},
	 {(296+169*sqrt(6))/1800, (88+7*sqrt(6))/360, (-2-3*sqrt(6))/225},
	 {(16-sqrt(6))/36, (16+sqrt(6))/36, 1./9}};
    //const nr_double_t radau_b[3] =
    //	{(16-sqrt(6))/36, (16+sqrt(6))/36, 1./9};
    const nr_double_t radau_c[3] =
	{(4-sqrt(6))/10, (4+sqrt(6))/10, 1};

    nr_double_t saveCurrent = current - delta;
    int n = countNodes () + countVoltageSources ();

    if (updateMatrix)
    {
	current = saveCurrent;
	// Is this necessary?
	for (int k = 0; k < n; k++)
	    x->set(k, SOL (0)->get(k));
	calculate ();
	createMatrix ();

	// Compute the (naive) Jacobian for the RadauIIA method
	for (int i = 0; i < 3; i++)
	    for (int j = 0; j < 3; j++)
		for (int k = 0; k < n; k++)
		    for (int l = 0; l < n; l++)
		    {
			nr_double_t g = radau_A[i][j] * A->get(k, l);

			if (i == j)
			    g += F->get(k, l) / delta;

			MA->set(i*n + k, j*n + l, g);
		    }
	//MA->print (1);
    }
    
    // Evaluate the circuit
    tvector<nr_double_t> fsave[3];

    for (int i = 0; i < 3; i++)
    {
	current = saveCurrent + radau_c[i] * delta;
	for (int k = 0; k < n; k++)
	    x->set(k, mx->get(i*n + k));
	saveSolution ();
	calculate ();
	createMatrix ();

	fsave[i] = *z - *A * *x;

	for (int k = 0; k < n; k++)
	{
	    nr_double_t mzk = 0;
	    for (int l = 0; l < n; l++)
		mzk -= F->get(k, l) * dmxsum->get(i*n + l) / delta;
	    mz->set(i*n + k, mzk);
	}
    }

    // Combine
    for (int i = 0; i < 3; i++)
	for (int k = 0; k < n; k++)
        {
	    nr_double_t mzk = mz->get(i*n + k);
	    for (int j = 0; j < 3; j++)
		mzk += radau_A[i][j] * fsave[j].get(k);
	    
	    mz->set(i*n + k, mzk);
	}

    // We should already be there:
    // current = saveCurrent + delta;
    assert (current - delta - saveCurrent < 1e-16);
}

void trsolver::setInitX (void)
{
    int n = countNodes () + countVoltageSources ();

    switch (corrType)
    {
    case INTEGRATOR_RADAU5:
	for (int i = 0; i < 3; i++)
	    for (int k = 0; k < n; k++)
		mx->set (i*n + k, x->get (k));
	break;
    default:
	nasolver::setInitX ();
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
	    x->set (k, mx->get(2*n + k));
	    dx->set (k, dmxsum->get(2*n + k));
	}
	break;
    default:
	nasolver::extractSol ();
    }
}

int trsolver::getSysSize (void)
{
    int n = countVoltageSources () + countNodes ();

    if (corrType == INTEGRATOR_RADAU5)
	return n * 3;
    else
	return n;
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
