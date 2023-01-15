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
	
	private final Map<ContinuousDistribution, List<Integer>> distributions;
	private final Map<GaussianDistribution, ContinuousDistribution> randomEffectDistributions;
	private final List<GaussianDistribution> randomEffectList;
	private int nbElements;

	MetropolisHastingsPriorHandler() {
		distributions = new LinkedHashMap<ContinuousDistribution, List<Integer>>();
		randomEffectDistributions = new HashMap<GaussianDistribution, ContinuousDistribution>();
		randomEffectList = new ArrayList<GaussianDistribution>();
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
			List<Integer> indices = distributions.get(d);
			realizedParameters.setElements(indices, thisR);
		}
		return realizedParameters;
	}

	/**
	 * Update the variance of the random effects on the fly.
	 * @param d the distribution
	 * @param realizedParameters the realized parameters
	 */
	private void updateRandomEffectVariance(ContinuousDistribution d, Matrix realizedParameters) {
		if (randomEffectDistributions.containsKey(d)) {	// it is a random effect. So we must update its variance
			ContinuousDistribution varianceDist = randomEffectDistributions.get(d);
			int index = distributions.get(varianceDist).get(0);	// TODO FP MF2021-11-01 here we assume that there is only one index 
			Matrix realizedRandomEffectVariance = realizedParameters.getSubMatrix(index, index, 0, 0);
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
				List<Integer> indices = distributions.get(d);
				double thisProb = d.getProbabilityDensity(realizedParameters.getSubMatrix(indices, null));
				if (thisProb == 0d) {
					return Double.NEGATIVE_INFINITY;
				}
				logProb += Math.log(thisProb);
			}
//			updateRandomEffectVarianceIfNeedsBe(d, m);
//			List<Integer> indices = distributions.get(d);
//			double thisProb = d.getProbabilityDensity(m.getSubMatrix(indices, null));
//			if (thisProb == 0d) {
//				return 0d;
//			}
//			logProb += Math.log(thisProb);
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
//		for (GaussianDistribution d : randomEffectList) {
//			updateRandomEffectVariance(d, realizedParameters);
//			List<Integer> indices = distributions.get(d);
//			double thisProb = d.getProbabilityDensity(realizedParameters.getSubMatrix(indices, null));
//			if (thisProb == 0d) {
//				return Double.NEGATIVE_INFINITY;
//			}
//			logProb += Math.log(thisProb);
//		}
		return logProb;
	}
	
	private double getProbabilityDensityOfThisRandomEffect(Matrix realizedParameters, int i) {
		if (randomEffectList.isEmpty()) {
			return 1d;		// no random effect then prob = 1
		} else {
			GaussianDistribution d = randomEffectList.get(i);
			updateRandomEffectVariance(d, realizedParameters);
			List<Integer> indices = distributions.get(d);
			return d.getProbabilityDensity(realizedParameters.getSubMatrix(indices, null));
		}
		
	}
	
	

	public void addFixedEffectDistribution(ContinuousDistribution dist, Integer... indices) {
		List<Integer> ind = Arrays.asList(indices);
		distributions.put(dist, ind);
		nbElements += ind.size();
	}

	public void addRandomEffectVariance(GaussianDistribution dist, ContinuousDistribution variancePrior, Integer... indices) {
		addFixedEffectDistribution(dist, indices);
		randomEffectDistributions.put(dist, variancePrior);
		randomEffectList.add(dist);
	}

}
