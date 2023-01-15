/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2021 Mathieu Fortin for Rouge Epicea.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.stats.estimators.mcmc;

import repicea.math.Matrix;
import repicea.stats.distributions.GaussianDistribution;
import repicea.stats.estimators.AbstractEstimator.EstimatorCompatibleModel;

/**
 * Ensure the compatibility with the Metropolis-Hastings algorithm.
 * @author Mathieu Fortin - November 2021
 */
public interface MetropolisHastingsCompatibleModel extends EstimatorCompatibleModel {

	/**
	 * Return the log-likelihood of the parameters. <br>
	 * <br>
	 * The model implementation is handled by the class implementing this interface. In the
	 * context of mixed-effects model, the vector of parameters (the argument parms) must 
	 * also include the random effects.
	 * 
	 * @param parms the model parameters (a Matrix instance)
	 * @return the log-likelihood of the parameters.
	 */
	public double getLogLikelihood(Matrix parms);

	/**
	 * Return the number of subjects. <br>
	 * <br>
	 * If the model is a mixed-effects model, the number of subjects must match the number of random effects.
	 * Otherwise, it should be the number of observations. 
	 * 
	 * @return an integer
	 */
	public int getNbSubjects();
	
	public double getLikelihoodOfThisSubject(Matrix parms, int subjectId);
	
	public GaussianDistribution getStartingParmEst(double coefVar);
	
}
