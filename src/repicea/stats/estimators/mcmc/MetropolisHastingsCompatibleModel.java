/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2021-24 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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
	 * Return the log-likelihood of the parameters. <p>
	 * The model implementation is handled by the class implementing this interface. In the
	 * context of mixed-effects model, the vector of parameters (the argument parms) must 
	 * also include the random effects.<p>
	 * By default, this method 
	 * <ol>
	 * <li> iterates on the subjects, 
	 * <li> computes their likelihood, 
	 * <li> converts it to log-likelihood and
	 * <li> returns the sum of the subject log-likelihoods
	 * <ol>
	 * @param parms the model parameters (a Matrix instance)
	 * @return the log-likelihood of the parameters.
	 */
	public default double getLogLikelihood(Matrix parms) {
		double llk = 0;
		for (int i = 0; i < getNbSubjects(); i++) {
			llk += Math.log(getLikelihoodOfThisSubject(parms, i));
		}
		return llk;
	}

	
	/**
	 * Return the number of subjects. <p>
	 * If the model is a mixed-effects model, the number of subjects must match the number of random effects.
	 * Otherwise, it should be the number of observations. 
	 * 
	 * @return an integer
	 */
	public int getNbSubjects();
	
	/**
	 * Provide the likelihood of a particular subject.<p>
	 * In the mixed model implementation, the subject is the highest hierarchical level (e.g. the plot).
	 * @param parms the parameter estimates
	 * @param subjectId an integer that stands for the id of the subject
	 * @return the likelihood (a double)
	 */
	public double getLikelihoodOfThisSubject(Matrix parms, int subjectId);
	
	/**
	 * Return the sampler.<p>
	 * The sampler contains the starting values of the parameter estimates plus
	 * their variances. The variances are often assumed to be the square of the product 
	 * of the parameter estimate by the coefVar argument. 
	 * @param coefVar a multiplicative factor for the variance. 
	 * @return a GaussianDistribution instance that will act as the sampler in the Metropolis-Hastings algorithm
	 */
	public GaussianDistribution getStartingParmEst(double coefVar);
	
	
	/**
	 * Set the prior distributions of the different fixed and random parameters.<p>
	 * If random effects appear in the model, their standard deviation should be modeled under the assumption
	 * it is uniformly distributed (Gelman 2006). 
	 * These distributions are set using the {@link MetropolisHastingsPriorHandler#addFixedEffectDistribution(repicea.stats.distributions.ContinuousDistribution, int)} 
	 * and the {@link MetropolisHastingsPriorHandler#addRandomEffectStandardDeviation(GaussianDistribution, repicea.stats.distributions.ContinuousDistribution, int)} public methods.
	 * @param handler a MetropolisHastingsPriorHandler instance
	 * @see <a href=https://doi.org/10.1214/06-BA117A> Gelman 2006 </a>
	 */
	public void setPriorDistributions(MetropolisHastingsPriorHandler handler);
	
}
