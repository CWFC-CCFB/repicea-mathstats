/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
 * Copyright (C) 2025 His Majesty the King in right of Canada
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
package repicea.stats.sampling;

/**
 * An implementation of the estimator of the mean and the total for finite populations.
 * Population units are assumed to have the same size. The variance corrected 
 * is corrected by the (1 - n/N) factor.
 * @author Mathieu Fortin - April 2016, January 2025
 */
@SuppressWarnings("serial")
public class FinitePopulationEstimate extends PopulationMeanEstimate implements FinitePopulationPointEstimate {

	final double populationSize;
	
	/**
	 * Constructor with population size.
	 * @param populationSize the number of units in the population.
	 */
	public FinitePopulationEstimate(double populationSize) {
		super();
		this.populationSize = populationSize;
	}

	@Override
	protected final void recalculate() {
		mean = computeMeanInternally();
		double finitePopulationCorrectionFactor = 1d - getSampleSize()/getPopulationSize();
		variance = computeVarianceInternally().scalarMultiply(finitePopulationCorrectionFactor);
		getDistribution().setMean(mean);		// the mean and variance and not tied to the the distribution
		getDistribution().setVariance(variance);	// consequently, they have to be specified before drawing the random deviates
		needsToBeRecalculated = false;
	}

	@Override
	public double getPopulationSize() {return populationSize;}

	@Override
	FinitePopulationEstimate getEmptyEstimate() {
		return new FinitePopulationEstimate(getPopulationSize());
	}

}
