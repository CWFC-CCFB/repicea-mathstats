/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import repicea.math.Matrix;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.stats.sampling.SamplingUtility;

/**
 * A class of static methods for quantile calculation or estimation. 
 * @author Mathieu Fortin - August 2022
 */
public class QuantileUtility {
	
	/**
	 * Return the quantile of a distribution estimated from a sample. <br>
	 * <br>
	 * The quantile is calculated following the Definition 8 found in <a href=https://doi.org/10.1080/00031305.1996.10473566>
	 * Hyndman, R. J. and Fan, Y. 1996. Sample quantiles in statistical packages. The American Statistician
	 * 50(4): 361-365. </a> <br>
	 * <br>
	 * The weighted version implements the technique shown on the wikipedia entry for <a href=https://en.wikipedia.org/wiki/Percentile#The_weighted_percentile_method>
	 *  Percentile </a>.
	 * 
	 * @param sample the sample of the distribution
	 * @param p the probability of the quantile (between 0 and 1)
	 * @param weights an optional list of integers representing the weighting (must be positive)
	 * @return the estimated quantile of the distribution
	 */
	public static double getQuantileEstimateFromSample(List<Double> sample, double p, List<Double> weights) {
		if (weights == null)
			return getInternalUnweightedQuantileEstimationFromSample(sample, p, true);
		else 
			return getInternalWeightedQuantileEstimationFromSample(sample, p, weights, true);
	}

	/**
	 * Internal estimation for unweighted quantiles.
	 * @param sample the sample of the distribution
	 * @param p the probability of the quantile (between 0 and 1)
	 * @param weights an optional list of integers representing the weighting (must be positive)
	 * @param boolean performChecks checks whether the input are correct
	 * @return a double
	 */
	private static double getInternalUnweightedQuantileEstimationFromSample(List<Double> sample, double p, boolean performChecks) {
		if (performChecks) 
			checkInputBeforeQuantileEstimation(sample, p, null);
		int nbObs = sample.size(); // default value
		
		List<Double> copyList = new ArrayList<Double>(nbObs);
		copyList.addAll(sample);
		Collections.sort(copyList);
		
		double N = copyList.size();
		double h = (N + 1d/3) * p + 1d/3;
		int h_floor = (int) Math.floor(h);
		int h_ceiling = (int) Math.ceil(h);
		double x_floor = copyList.get(h_floor - 1);
		double q = x_floor + (h - h_floor) * (copyList.get(h_ceiling - 1) - x_floor);
		return q;
	}
	
	
	private static void checkInputBeforeQuantileEstimation(List<Double> sample, double p, List<Double> weights) {
		if (p < 0d || p > 1d)
			throw new InvalidParameterException("The p argument must range from 0 to 1!");
		if (sample == null || sample.isEmpty()) {
			throw new InvalidParameterException("The sample argument should be a non empty list of doubles!");
		}
		if (weights != null) {
			if (weights.size() != sample.size()) {
				throw new InvalidParameterException("If not null, the weights argument should be a list of the same size as sample!");
			}
			if (weights.stream().anyMatch(n -> n <= 0)) {
				throw new InvalidParameterException("If not null, the weights argument must contain strictly positive values (i.e. > 0)!");
			}
		}
	}
	
	
	/**
	 * Return an estimated quantile as well as it variability. <br>
	 * <br>
	 * The variability is assessed through bootstrap.
	 * 
	 * @param sample a list of doubles
	 * @param p the quantile probability
	 * @param weights a list of doubles 
	 * @param nReal the number of bootstrap realizations 
	 * @return a MonteCarloEstimate
	 * @throws InterruptedException 
	 */
	public static MonteCarloEstimate getQuantileEstimateFromSample(List<Double> sample, 
			double p, 
			List<Double> weights, 
			int nReal) {
		checkInputBeforeQuantileEstimation(sample, p, weights);
		if (nReal <= 0) {
			throw new InvalidParameterException("The nReal argument should be a strictly positive integer (i.e. > 0)!");
		}
		int n = sample.size();
//		List<Integer> indices = new ArrayList<Integer>(n);
//		for (int i = 0; i < n; i++)
//			indices.add(i);
		
		List<Double> bootstrapSample = new ArrayList<Double>(n);
		List<Double> bootstrapWeights = new ArrayList<Double>(n);
		List<Integer> selectedIndices = new ArrayList<Integer>(n);
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		for (int i = 0; i < nReal; i++) {
//			List<Integer> selectedIndices = SamplingUtility.getSample(indices, indices.size(), true);
			SamplingUtility.getSampleIndex(n, n, true, selectedIndices);
			bootstrapSample.clear();
			bootstrapWeights.clear();
			for (Integer index : selectedIndices) {
				bootstrapSample.add(sample.get(index));
				if (weights != null) {
					bootstrapWeights.add(weights.get(index));
				}
			}
			double quantile = weights == null ? 
					getInternalUnweightedQuantileEstimationFromSample(bootstrapSample, p, false) :	// no checks needed they've been done at the beginning of the method
						getInternalWeightedQuantileEstimationFromSample(bootstrapSample, p, weights, false);
			estimate.addRealization(new Matrix(1,1,quantile,0));
		}
		return estimate;
	}

	
	/**
	 * Return the quantile of a distribution calculated from the population. <br>
	 * <br>
	 * 
	 * @param population the population
	 * @param p the probability of the quantile (between 0 and 1)
	 * @return the calculated quantile of the distribution
	 */
	public static double getQuantileFromPopulation(List<Double> population, double p) {
		if (p < 0d || p > 1d)
			throw new InvalidParameterException("The p argument must range from 0 to 1!");
		if (population == null || population.isEmpty()) {
			throw new InvalidParameterException("The population argument should be a non empty list of doubles!");
		}
		List<Double> copyList = new ArrayList<Double>();
		copyList.addAll(population);
		Collections.sort(copyList);
		double N = copyList.size();
		int h = (int) Math.round(p*N);
		double q = copyList.get(h - 1);
		return q;
	}

	
	
	


	static class WeightedSampleUnit implements Comparable<WeightedSampleUnit> {

		final double value;
		final double weight;
		double p_n;

		WeightedSampleUnit(double value, double weight) {
			this.value = value;
			this.weight = weight;
		}
		
		@Override
		public int compareTo(WeightedSampleUnit o) {
			if (value < o.value) return -1;
			else if (value > o.value) return 1;
			else return 0;
		}
		
	}
	
	/**
	 * Internal estimation for weighted quantiles.
	 * @param sample the sample of the distribution
	 * @param p the probability of the quantile (between 0 and 1)
	 * @param weights an optional list of integers representing the weighting (must be positive)
	 * @param boolean performChecks checks whether the input are correct
	 * @return a double
	 */
	protected static double getInternalWeightedQuantileEstimationFromSample(List<Double> sample, double p, List<Double> weights, boolean performChecks) {
		if (performChecks) 
			checkInputBeforeQuantileEstimation(sample, p, weights);
		
		List<WeightedSampleUnit> copyList = new ArrayList<WeightedSampleUnit>(sample.size());
		for (int i = 0; i < sample.size(); i++) {
			WeightedSampleUnit unit = weights != null ?
					new WeightedSampleUnit(sample.get(i), weights.get(i)) :
						new WeightedSampleUnit(sample.get(i), 1d);
			copyList.add(unit);
		}
		setRanks(copyList);
		int floor_x = findLowerBoundAmongRanks(p, copyList);
		if (floor_x == -1) {	// smaller than the first value
			return copyList.get(0).value;
		} else if (floor_x == copyList.size() - 1) {
			return copyList.get(floor_x).value; 	// larger than the last value
		} else {
			double value = copyList.get(floor_x).value;
			double p_k = copyList.get(floor_x).p_n;
			if (p > p_k) {
				double p_kPlusOne = copyList.get(floor_x + 1).p_n;
				double ratio = (p - p_k) / (p_kPlusOne - p_k);
				value += ratio * (copyList.get(floor_x + 1).value - value);
			}
			return value;
		}
	}
	
	protected static void setRanks(List<WeightedSampleUnit> units) {
		Collections.sort(units);
		double S_N = 0d;
		for (WeightedSampleUnit u : units) {
			S_N += u.weight;
		}
		double S_n = 0d;
		for (WeightedSampleUnit u : units) {
			S_n += u.weight;
			u.p_n = (S_n - 1d/3 * u.weight) / (S_N + 1d/3 * u.weight) ;
		}
	}
	

	protected static int findLowerBoundAmongRanks(double p, List<WeightedSampleUnit> units) {
		return findLowerBoundAmongRanks(p, units, 0, units.size() - 1);
	}
	
	protected static int findLowerBoundAmongRanks(double p, List<WeightedSampleUnit> units, int lowerBound, int upperBound) {
		if (p < units.get(lowerBound).p_n) {
			return -1;
		} 
//		System.out.println("Lower = " + lowerBound + "; Upper = " + upperBound);
		if (p < units.get(upperBound).p_n) {
			if (p >= units.get(lowerBound).p_n && p < units.get(lowerBound + 1).p_n) {
				return lowerBound;
			} else {
				int midPoint = ((Number) ((upperBound + lowerBound) * .5)).intValue();
				if (p < units.get(midPoint).p_n) {
					return findLowerBoundAmongRanks(p, units, lowerBound, midPoint);
				} else {
					return findLowerBoundAmongRanks(p, units, midPoint, upperBound);
				}
			}
		} else {
			return upperBound;	// return the upperBound
		}
	}
	
	

}
