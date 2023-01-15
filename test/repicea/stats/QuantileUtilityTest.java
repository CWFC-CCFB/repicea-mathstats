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

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.stats.QuantileUtility.WeightedSampleUnit;
import repicea.stats.estimates.MonteCarloEstimate;

public class QuantileUtilityTest {

	@Test
	public void quantileTest() {
		List<Double> myList = new ArrayList<Double>();
		for (double i = 0; i < 100;  i++) {
			myList.add(i);
		}
		double q = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, null);
		Assert.assertEquals("Comparing quantile", 29.4333333333, q, 1E-8);
	}
	
	@Test
	public void quantileTest2() {
		List<Double> myList = new ArrayList<Double>();
		for (double i = 0; i < 200;  i++) {
			myList.add(i);
		}
		double q = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, null);
		Assert.assertEquals("Comparing quantile", 59.4333333333, q, 1E-8);
	}

	
	@Test
	public void quantileTest3() {
		List<Double> myList = new ArrayList<Double>();
		for (double i = 0; i < 200;  i = i + 2) {
			myList.add(i);
		}
		double q = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, null);
		Assert.assertEquals("Comparing quantile", 58.8666666666, q, 1E-8);
	}

	
	@Test
	public void quantileTest4() {
		List<Double> myList = new ArrayList<Double>();
		for (double i = 0; i < 100;  i++) {
			myList.add(i * i);
		}
		double q = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, null);
		Assert.assertEquals("Comparing quantile", 866.5666666666, q, 1E-8);
	}
	
	@Test
	public void quantileTest5() {
		List<Double> myList = new ArrayList<Double>();
		for (double i = 0; i < 200;  i++) {
			myList.add(i * i);
		}
		double q = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, null);
		Assert.assertEquals("Comparing quantile", 3532.5666666666666, q, 1E-8);
	}

	
	@Test
	public void quantileTest6() {
		List<Double> myList = new ArrayList<Double>();
		for (double i = 0; i < 200;  i = i + 2) {
			myList.add(i * i);
		}
		double q = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, null);
		Assert.assertEquals("Comparing quantile", 3466.2666666666, q, 1E-8);
	}

	


	@Test
	public void findingIndexTest1() {
		List<WeightedSampleUnit> myList = new ArrayList<WeightedSampleUnit>();
		for (double i = 0; i < 200;  i = i + 2) {
			myList.add(new WeightedSampleUnit(i * i, i == 0 ? 20d : 1d));
		}
		QuantileUtility.setRanks(myList);
		int index = QuantileUtility.findLowerBoundAmongRanks(0.2, myList);
		Assert.assertEquals("Comparing index", 4, index);

		index = QuantileUtility.findLowerBoundAmongRanks(1.0, myList);
		Assert.assertEquals("Comparing index", 99, index);

		index = QuantileUtility.findLowerBoundAmongRanks(0, myList);
		Assert.assertEquals("Comparing index", -1, index);
	}

	
	@Test
	public void weightedQuantileTest1() {
		List<Double> myList = new ArrayList<Double>();
		List<Double> myWeights = new ArrayList<Double>();
		for (double i = 0; i < 200;  i = i + 2) {
			myList.add(i * i);
			if (i==0) {
				myWeights.add(20d);
			} else {
				myWeights.add(1d);
			}
		}
		double q = QuantileUtility.getInternalWeightedQuantileEstimationFromSample(myList, 0.3, myWeights, true);
		Assert.assertEquals("Comparing quantile", 1041.6, q, 1E-8);
	}
	
	@Test
	public void weightedQuantileTest2() {
		List<Double> myList = new ArrayList<Double>();
		List<Double> myWeights = new ArrayList<Double>();
		for (double i = 0; i < 200;  i = i + 2) {
			myList.add(i * i);
			if (i==0) {
				myWeights.add(20d);
			} else {
				myWeights.add(1d);
			}
		}
		double q = QuantileUtility.getInternalWeightedQuantileEstimationFromSample(myList, 0, myWeights, true);
		Assert.assertEquals("Comparing quantile", myList.get(0), q, 1E-8);
	}

	@Test
	public void weightedQuantileTest3() {
		List<Double> myList = new ArrayList<Double>();
		List<Double> myWeights = new ArrayList<Double>();
		for (double i = 0; i < 200;  i = i + 2) {
			myList.add(i * i);
			if (i==0) {
				myWeights.add(20d);
			} else {
				myWeights.add(1d);
			}
		}
		double q = QuantileUtility.getInternalWeightedQuantileEstimationFromSample(myList, 1.0, myWeights, true);
		Assert.assertEquals("Comparing quantile", myList.get(myList.size() - 1), q, 1E-8);
	}

	
	@Test
	public void weightedQuantileTest4() {
		List<Double> myList = new ArrayList<Double>();
		List<Double> myWeights = new ArrayList<Double>();
		for (double i = 0; i < 200;  i = i + 2) {
			myList.add(i * i);
			myWeights.add(1d);
		}
		double q1 = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, myWeights);
		double q2 = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, null);
		Assert.assertEquals("Comparing quantile", q1, q2, 1E-8);
	}

	@Test
	public void weightedQuantileEstimateTest1() {
		List<Double> myList = new ArrayList<Double>();
		List<Double> myWeights = new ArrayList<Double>();
		for (double i = 0; i < 100;  i = i + 2) {
			myList.add(i * i);
			myWeights.add(1d);
		}
		MonteCarloEstimate est = QuantileUtility.getQuantileEstimateFromSample(myList, 0.3, myWeights, 1000);	// first call is to make sure the class is loaded
		double mean = est.getMean().getValueAt(0, 0);
		double variance = est.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing mean", 895.6358666666663, mean, 35);
		Assert.assertEquals("Testing variance", 137033.9766902726, variance, 3E4);
	}

}
