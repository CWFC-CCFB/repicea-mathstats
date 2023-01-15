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
package repicea.stats.functions;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import repicea.stats.integral.TrapezoidalRule;


public class EmpiricalDistributionProbabilityDensityFunctionTest {

	@Test
	public void simpleUnweightedValuesTest() {
		List<Double> values = new ArrayList<Double>();
		for (double i = 1; i <= 100; i++) {
			values.add(i);
		}
		EmpiricalDistributionProbabilityDensityFunction edf = new EmpiricalDistributionProbabilityDensityFunction(values, null);
		edf.setVariableValue(0, 50);
		double density = edf.getValue();
		Assert.assertEquals("Comparing densities", 0.010553531401233376, density, 1E-8);
	}

	@Test
	public void simpleUnweightedValuesExceedingLowerBoundTest() {
		List<Double> values = new ArrayList<Double>();
		for (double i = 1; i <= 100; i++) {
			values.add(i);
		}
		EmpiricalDistributionProbabilityDensityFunction edf = new EmpiricalDistributionProbabilityDensityFunction(values, null);
		edf.setVariableValue(0, 0);
		double density = edf.getValue();
		Assert.assertEquals("Comparing densities", 0d, density, 1E-8);
	}

	@Test
	public void simpleUnweightedValuesExactLowerBoundTest() {
		List<Double> values = new ArrayList<Double>();
		for (double i = 1; i <= 100; i++) {
			values.add(i);
		}
		EmpiricalDistributionProbabilityDensityFunction edf = new EmpiricalDistributionProbabilityDensityFunction(values, null);
		edf.setVariableValue(0, 1);
		double density = edf.getValue();
		Assert.assertEquals("Comparing densities", 0.009594119455666702, density, 1E-8);
	}

	@Test
	public void simpleUnweightedValuesExceedingUpperBoundTest() {
		List<Double> values = new ArrayList<Double>();
		for (double i = 1; i <= 100; i++) {
			values.add(i);
		}
		EmpiricalDistributionProbabilityDensityFunction edf = new EmpiricalDistributionProbabilityDensityFunction(values, null);
		edf.setVariableValue(0, 101);
		double density = edf.getValue();
		Assert.assertEquals("Comparing densities", 0d, density, 1E-8);
	}

	@Test
	public void simpleUnweightedValuesExactUpperBoundTest() {
		List<Double> values = new ArrayList<Double>();
		for (double i = 1; i <= 100; i++) {
			values.add(i);
		}
		EmpiricalDistributionProbabilityDensityFunction edf = new EmpiricalDistributionProbabilityDensityFunction(values, null);
		edf.setVariableValue(0, 100);
		double density = edf.getValue();
		Assert.assertEquals("Comparing densities", 0.011512943346800067, density, 1E-8);
	}

	@Test
	public void simpleUnweightedValuesTestForCumulativeProbability() {
		TrapezoidalRule tr = new TrapezoidalRule(0.1);
		tr.setLowerBound(-5);
		tr.setUpperBound(105);
		List<Double> values = new ArrayList<Double>();
		for (double i = 1; i <= 100; i++) {
			values.add(i);
		}
		EmpiricalDistributionProbabilityDensityFunction edf = new EmpiricalDistributionProbabilityDensityFunction(values, null);
		double cumulProb = tr.getIntegralApproximation(edf, 0, false);	// false : integrate over the variable
		Assert.assertEquals("Comparing densities", 0.9979340827157241, cumulProb, 1E-8);
	}

	// TODO implement tests with weights MF 2022-07-28
}
