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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.distributions.ContinuousDistribution;
import repicea.stats.distributions.GaussianDistribution;

/** 
 * A class to handle prior distributions.
 * @author Mathieu Fortin - November 2021
 */
public class MetropolisHastingsPriorHandler {
	
	private final Map<ContinuousDistribution, Integer> distributions;
	private final Map<GaussianDistribution, ContinuousDistribution> randomEffectDistributions;
	private final List<GaussianDistribution> randomEffectList;
	private int nbElements;
	private List<Integer> commodityArray; // a commodity array for random effect density calculation

	/**
	 * Package constructor.
	 */
	MetropolisHastingsPriorHandler() {
		distributions = new LinkedHashMap<ContinuousDistribution, Integer>();
		randomEffectDistributions = new HashMap<GaussianDistribution, ContinuousDistribution>();
		randomEffectList = new ArrayList<GaussianDistribution>();
	}

	private List<Integer> getCommodityArray() {
		if (commodityArray == null) {
			commodityArray = new ArrayList<Integer>();
		}
		return commodityArray;
	}
	
	/**
	 * Provide a realization of the parameters (fixed and random).
	 * @return a Matrix instance
	 */
	Matrix getRandomRealization() {
		Matrix realizedParameters = new Matrix(nbElements, 1);
		for (ContinuousDistribution d : distributions.keySet()) {
			updateRandomEffectVariance(d, realizedParameters);
			Matrix thisR = d.getRandomRealization();
			int index = distributions.get(d);
			realizedParameters.setElements(turnIntIntoArray(index), thisR);
		}
		return realizedParameters;
	}

	/**
	 * Update the variance of the random effects on the fly.
	 * @param d a ContinuousDistribution instance that stands for the distribution of a particular random effect
	 * @param realizedParameters the realized parameters generated from one realization of the Metropolis-Hastings algorithm.
	 */
	private void updateRandomEffectVariance(ContinuousDistribution d, Matrix realizedParameters) {
		if (randomEffectDistributions.containsKey(d)) {	// it is a random effect. So we must update its variance
			ContinuousDistribution stdPrior = randomEffectDistributions.get(d);
			int index = distributions.get(stdPrior);	
			Matrix realizedRandomEffectVariance = realizedParameters.getSubMatrix(index, index, 0, 0).elementWisePower(2d); // power 2 because we are dealing with the standard deviation
			((GaussianDistribution) d).setVariance(SymmetricMatrix.convertToSymmetricIfPossible(realizedRandomEffectVariance));
		}
	}

	/**
	 * Return the log probability density of the parameters (only fixed) with respect to the priors.
	 * @param realizedParameters 
	 * @return a double
	 */
	double getLogProbabilityDensity(Matrix realizedParameters) {
		double logProb = 0;
		
		for (ContinuousDistribution d : distributions.keySet()) {
			if (!randomEffectDistributions.containsKey(d)) {	// we do not consider the random effects in the probability density of the prior
				int index = distributions.get(d);
				double thisProb = d.getProbabilityDensity(realizedParameters.getSubMatrix(turnIntIntoArray(index), null));
				if (thisProb == 0d) {
					return Double.NEGATIVE_INFINITY;
				}
				logProb += Math.log(thisProb);
			}
		}
		return logProb;
	}

	double getLogProbabilityDensityOfRandomEffects(Matrix realizedParameters) {
		double logProb = 0;
		for (int i = 0; i < randomEffectList.size(); i++) {
			double thisProb = getProbabilityDensityOfThisRandomEffect(realizedParameters, i);
			if (thisProb == 0d) {
				return Double.NEGATIVE_INFINITY;
			} else {
				logProb += Math.log(thisProb);
			}
		}
		return logProb;
	}
	
	double getProbabilityDensityOfThisRandomEffect(Matrix realizedParameters, int i) {
		if (randomEffectList.isEmpty()) {
			return 1d;		// no random effect then prob = 1
		} else {
			GaussianDistribution d = randomEffectList.get(i);
			updateRandomEffectVariance(d, realizedParameters);
			int index = distributions.get(d);
			return d.getProbabilityDensity(realizedParameters.getSubMatrix(turnIntIntoArray(index), null));
		}
	}

	
	private List<Integer> turnIntIntoArray(int i) {
		List<Integer> array = getCommodityArray();
		array.clear();
		array.add(i);
		return array;
	}
	
	/**
	 * Add a prior distribution for a fixed effect parameter.
	 * @param dist a ContinuousDistribution instance
	 * @param index the index of the parameter associated with this distribution
	 */
	public void addFixedEffectDistribution(ContinuousDistribution dist, int index) {
//		List<Integer> ind = Arrays.asList(indices);
		distributions.put(dist, index);
//		nbElements += ind.size();
		nbElements++;
	}

	/**
	 * Add a prior distribution for a random effect standard deviation.<p>
	 * 
	 * When assuming the existence of a random effect, it is good practice to assume
	 * the standard deviation is uniformly distributed and not the variance (Gelman 2006). 
	 * @param dist a GaussianDistribution instance that stands for the distribution of a particular random effect.
	 * @param stdPrior a ContinuousDistribution instance that stands for the prior distribution of the random effect standard deviation.
	 * @param index the index of the parameter associated with this random effect
	 * @see <a href=https://doi.org/10.1214/06-BA117A> Gelman 2006 </a>
	 */
	public void addRandomEffectStandardDeviation(GaussianDistribution dist, ContinuousDistribution stdPrior, int index) {
		addFixedEffectDistribution(dist, index);
		randomEffectDistributions.put(dist, stdPrior);
		randomEffectList.add(dist);
	}

	/**
	 * Check if any prior distribution has been set.
	 * @return true if there is NO prior distribution in the handler yet.
	 */
	public boolean isEmpty() {
		return distributions.isEmpty();
	}
	
	/**
	 * Clear the prior distributions.
	 */
	public void clear() {
		nbElements = 0;
		distributions.clear();
		randomEffectDistributions.clear();
		randomEffectList.clear();
	}

}
