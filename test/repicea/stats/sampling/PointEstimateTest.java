/*
 * This file is part of the repicea-mathstats library.
 *
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

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.junit.FixMethodOrder;
import org.junit.Test;
import org.junit.runners.MethodSorters;

import repicea.math.Matrix;

@FixMethodOrder(MethodSorters.NAME_ASCENDING)
public class PointEstimateTest {

	@Test
	public void test01TotalAndVarianceForFinitePopulationMeanPointEstimate() {
		List<Double> sample = new ArrayList<Double>();
		sample.add(2d);
		sample.add(4d);
		sample.add(2d);
		sample.add(5d);
		sample.add(7d);
		sample.add(1d);
		sample.add(5d);
		sample.add(4d);
		sample.add(7d);
		FinitePopulationEstimate estimate = new FinitePopulationEstimate(1000);
		Matrix obs;
		for (int i = 0; i < sample.size(); i++) {
			double value = sample.get(i);
			obs = new Matrix(1,1);
			obs.setValueAt(0, 0, value);
			estimate.addObservation(obs, i + "");
		}
		Matrix total = estimate.getTotal();
		Assert.assertEquals("Testing the estimate of the total", 4111.11111111111, total.getValueAt(0, 0), 1E-8);
		Matrix totalVariance = estimate.getVarianceOfTotal();
		Assert.assertEquals("Testing the variance of the total", 507734.5679012318, totalVariance.getValueAt(0, 0), 1E-8);
	}
	

	
	private static void fillStratum(StratifiedPopulationEstimate estimate, String stratumName, List<Double> values) {
		Matrix obs;
		for (int i = 0; i < values.size(); i++) {
			double value = values.get(i);
			obs = new Matrix(1,1);
			obs.setValueAt(0, 0, value);
			estimate.addObservation(obs, i + "", stratumName);
		}
	}

	private static FinitePopulationEstimate producePopulationMeanEstimate(List<Double> values, double populationSize) {
		FinitePopulationEstimate meanEstimate = new FinitePopulationEstimate(populationSize);
		Matrix obs;
		for (int i = 0; i < values.size(); i++) {
			double value = values.get(i);
			obs = new Matrix(1,1);
			obs.setValueAt(0, 0, value);
			meanEstimate.addObservation(obs, i + "");
		}
		return meanEstimate;
	}

	
	@Test
	public void test02TotalAndVarianceStratifiedPopulationTotalEstimate() {
		StratifiedPopulationEstimate estimate = new StratifiedPopulationEstimate(Arrays.asList(new String[] {"Stratum1", "Stratum2", "Stratum3"}), 
				Arrays.asList(new Double[] {100d, 200d, 300d}));

		List<Double> sampleStratum1 = new ArrayList<Double>();
		sampleStratum1.add(2d);
		sampleStratum1.add(4d);
		sampleStratum1.add(2d);
		sampleStratum1.add(5d);
		fillStratum(estimate, "Stratum1", sampleStratum1);
		FinitePopulationEstimate meanEstimate1 = producePopulationMeanEstimate(sampleStratum1, 100d);
		
		List<Double> sampleStratum2 = new ArrayList<Double>();
		sampleStratum2.add(0d);
		sampleStratum2.add(6d);
		sampleStratum2.add(2d);
		sampleStratum2.add(5d);
		fillStratum(estimate, "Stratum2", sampleStratum2);
		FinitePopulationEstimate meanEstimate2 = producePopulationMeanEstimate(sampleStratum2, 200d);

		List<Double> sampleStratum3 = new ArrayList<Double>();
		sampleStratum3.add(2d);
		sampleStratum3.add(10d);
		sampleStratum3.add(2d);
		fillStratum(estimate, "Stratum3", sampleStratum3);
		FinitePopulationEstimate meanEstimate3 = producePopulationMeanEstimate(sampleStratum3, 300d);

		double popSize = estimate.getPopulationSize();
		Assert.assertEquals("Testing total population size", 600, popSize, 1E-8);
		Matrix actualTotal = estimate.getTotal();
		Matrix expectedTotal = meanEstimate1.getMean().scalarMultiply(meanEstimate1.getPopulationSize())
				.add(meanEstimate2.getMean().scalarMultiply(meanEstimate2.getPopulationSize()))
				.add(meanEstimate3.getMean().scalarMultiply(meanEstimate3.getPopulationSize()));
				
		Assert.assertEquals("Testing the estimate of the total", 
				expectedTotal.getValueAt(0, 0), 
				actualTotal.getValueAt(0, 0), 
				1E-8);
		Matrix actualTotalVariance = estimate.getVarianceOfTotal();
		Matrix expectedTotalVariance = meanEstimate1.getVariance().scalarMultiply(meanEstimate1.getPopulationSize() * meanEstimate1.getPopulationSize())
				.add(meanEstimate2.getVariance().scalarMultiply(meanEstimate2.getPopulationSize() * meanEstimate2.getPopulationSize()))
				.add(meanEstimate3.getVariance().scalarMultiply(meanEstimate3.getPopulationSize() * meanEstimate3.getPopulationSize()));

		Assert.assertEquals("Testing the variance of the total", 
				expectedTotalVariance.getValueAt(0,0), 
				actualTotalVariance.getValueAt(0, 0), 
				1E-8);
	}
	
	@Test
	public void test03FailConstructStratifiedPopulationTotalEstimate() {
		try {
			new StratifiedPopulationEstimate(Arrays.asList(new String[] {"Stratum1", "Stratum2", "Stratum2"}), 
					Arrays.asList(new Double[] {100d, 200d, 300d}));
			Assert.fail("Should have thrown an InvalidParameterException");
		} catch (InvalidParameterException e) {
			System.err.println(e.getMessage());
			System.out.println("This error was expected");
		}

	}

}
