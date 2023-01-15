/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.estimates;

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.distributions.EmpiricalDistribution;
import repicea.stats.sampling.PopulationUnitWithEqualInclusionProbability;

/**
 * This class implements an estimator of the mean of the population, which is the 
 * sample mean and the sample variance corrected by n/(n-1) being the estimator of the
 * variance.
 * @author Mathieu Fortin - April 2016
 */
@SuppressWarnings("serial")
public class PopulationMeanEstimate extends PointEstimate<PopulationUnitWithEqualInclusionProbability> {
		
	private final EmpiricalDistribution sample;

	/**
	 * Basic constructor without population size.
	 */
	public PopulationMeanEstimate() {
		super();
		sample = new EmpiricalDistribution();
	}

	/**
	 * Constructor with population size.
	 * @param populationSize the number of units in the population.
	 */
	public PopulationMeanEstimate(double populationSize) {
		super(populationSize);
		sample = new EmpiricalDistribution();
	}
	

	@Override
	protected Matrix getMeanFromDistribution() {
		return sample.getMean();
	}
	
	private int getSampleSize() {return getObservations().size();}
	
	@Override
	protected SymmetricMatrix getVarianceFromDistribution() {
		double smallAreaCorrectionFactor = 1d;
		if (isPopulationSizeKnown()) {
			smallAreaCorrectionFactor = 1d - getSampleSize()/getPopulationSize();
		}
		return sample.getVariance().scalarMultiply(1d/getSampleSize() * smallAreaCorrectionFactor);
	}

	@Override
	public Matrix getRandomDeviate() {
		getDistribution().setMean(getMean());		// the mean and variance and not tied to the the distribution
		getDistribution().setVariance(getVariance());	// consequently, they have to be specified before drawing the random deviates
		return super.getRandomDeviate();
	}

	@Override
	public void addObservation(PopulationUnitWithEqualInclusionProbability obs) {
		super.addObservation(obs);
		sample.addRealization(obs.getData());
	}

	@Override
	protected PopulationMeanEstimate add(PointEstimate<?> pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			PopulationMeanEstimate newEstimate = new PopulationMeanEstimate();
			PopulationMeanEstimate meanEstimate = (PopulationMeanEstimate) pointEstimate;
			for (String sampleId : getSampleIds()) {
				PopulationUnitWithEqualInclusionProbability thisUnit = getObservations().get(sampleId);
				PopulationUnitWithEqualInclusionProbability thatUnit = meanEstimate.getObservations().get(sampleId);
				PopulationUnitWithEqualInclusionProbability newUnit = new PopulationUnitWithEqualInclusionProbability(sampleId, thisUnit.getData().add(thatUnit.getData()));
				newEstimate.addObservation(newUnit);
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected PopulationMeanEstimate subtract(PointEstimate<?> pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			PopulationMeanEstimate newEstimate = new PopulationMeanEstimate();
			PopulationMeanEstimate meanEstimate = (PopulationMeanEstimate) pointEstimate;
			for (String sampleId : getSampleIds()) {
				PopulationUnitWithEqualInclusionProbability thisUnit = getObservations().get(sampleId);
				PopulationUnitWithEqualInclusionProbability thatUnit = meanEstimate.getObservations().get(sampleId);
				PopulationUnitWithEqualInclusionProbability newUnit = new PopulationUnitWithEqualInclusionProbability(sampleId, thisUnit.getData().subtract(thatUnit.getData())); 
				newEstimate.addObservation(newUnit);
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected PopulationMeanEstimate multiply(double scalar) {
		PopulationMeanEstimate newEstimate = new PopulationMeanEstimate();
		for (String sampleId : getSampleIds()) {
			PopulationUnitWithEqualInclusionProbability thisUnit = getObservations().get(sampleId);
			PopulationUnitWithEqualInclusionProbability newUnit = new PopulationUnitWithEqualInclusionProbability(sampleId, thisUnit.getData().scalarMultiply(scalar)); 
			newEstimate.addObservation(newUnit);
		}
		return newEstimate;
	}

}
