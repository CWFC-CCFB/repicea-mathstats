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
package repicea.math;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.functions.GaussianFunction;
import repicea.stats.distributions.utility.GaussianUtility;

public class ProductFunctionWrapperTest {
	
	@Test
	public void simpleTestWithTwoIndependentFunctions() {
		InternalMathematicalFunctionWrapper wrapper1 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{0,1}),
				Arrays.asList(new Integer[]{0}));
		InternalMathematicalFunctionWrapper wrapper2 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{2,3}),
				Arrays.asList(new Integer[]{1}));
		ProductFunctionWrapper pfw = new ProductFunctionWrapper(wrapper1, wrapper2);
		Matrix parms = new Matrix(4,1);
		parms.setValueAt(0, 0, 0);
		parms.setValueAt(1, 0, 1);
		parms.setValueAt(2, 0, 2);
		parms.setValueAt(3, 0, 3);
		pfw.setParameters(parms);

		Matrix m = new Matrix(1,2);
		m.setValueAt(0, 1, 2);
		
		pfw.setVariables(m);
		double observedValue = pfw.getValue();
		double expectedValue = GaussianUtility.getProbabilityDensity(0) * GaussianUtility.getProbabilityDensity(2, 2, 3);
		Assert.assertEquals("Comparing density product", expectedValue, observedValue, 1E-8);
		
		m.setValueAt(0, 0, -1);
		m.setValueAt(0, 1, 3);
		pfw.setVariables(m);
		observedValue = pfw.getValue();
		expectedValue = GaussianUtility.getProbabilityDensity(-1) * GaussianUtility.getProbabilityDensity(3, 2, 3);
		Assert.assertEquals("Comparing density product", expectedValue, observedValue, 1E-8);
	}
	
	
	@Test
	public void simpleTestWithTwoFunctionsDependingOnTheSameVariable() {
		InternalMathematicalFunctionWrapper wrapper1 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{0,1}),
				Arrays.asList(new Integer[]{0}));
		InternalMathematicalFunctionWrapper wrapper2 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{2,3}),
				Arrays.asList(new Integer[]{0}));
		ProductFunctionWrapper pfw = new ProductFunctionWrapper(wrapper1, wrapper2);
		
		Matrix parms = new Matrix(4,1);
		parms.setValueAt(0, 0, 0);
		parms.setValueAt(1, 0, 1);
		parms.setValueAt(2, 0, 2);
		parms.setValueAt(3, 0, 3);
		pfw.setParameters(parms);
		Matrix m = new Matrix(1,1);
		m.setValueAt(0, 0, 1);
		pfw.setVariables(m);
		
		double observedValue = pfw.getValue();
		double expectedValue = GaussianUtility.getProbabilityDensity(1) * GaussianUtility.getProbabilityDensity(1, 2, 3);
		Assert.assertEquals("Comparing density product", expectedValue, observedValue, 1E-8);
		
		m.setValueAt(0, 0, 2);
		pfw.setVariables(m);
		observedValue = pfw.getValue();
		expectedValue = GaussianUtility.getProbabilityDensity(2) * GaussianUtility.getProbabilityDensity(2, 2, 3);
		Assert.assertEquals("Comparing density product", expectedValue, observedValue, 1E-8);
	}

	@Test
	public void simpleTestWithTwoFunctionsDependingOnTheSameParameter() {
		InternalMathematicalFunctionWrapper wrapper1 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{0,1}),
				Arrays.asList(new Integer[]{0}));
		InternalMathematicalFunctionWrapper wrapper2 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{0,2}),
				Arrays.asList(new Integer[]{1}));
		ProductFunctionWrapper pfw = new ProductFunctionWrapper(wrapper1, wrapper2);
		
		Matrix parms = new Matrix(3,1);
		parms.setValueAt(0, 0, 0);
		parms.setValueAt(1, 0, 1);
		parms.setValueAt(2, 0, 2);
		pfw.setParameters(parms);

		Matrix m = new Matrix(1,2);
		m.setValueAt(0, 0, 1);
		m.setValueAt(0, 1, 3);
		pfw.setVariables(m);
		double observedValue = pfw.getValue();
		double expectedValue = GaussianUtility.getProbabilityDensity(1) * GaussianUtility.getProbabilityDensity(3, 0, 2);
		Assert.assertEquals("Comparing density product", expectedValue, observedValue, 1E-8);
		
		pfw.setParameterValue(0, 1);
		observedValue = pfw.getValue();
		expectedValue = GaussianUtility.getProbabilityDensity(1, 1, 1) * GaussianUtility.getProbabilityDensity(3, 1, 2);
		Assert.assertEquals("Comparing density product", expectedValue, observedValue, 1E-8);
	}

	@Test
	public void simpleGradientTestWithTwoIndependentFunctions() {
		InternalMathematicalFunctionWrapper wrapper1 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{0,1}),
				Arrays.asList(new Integer[]{0}));
		InternalMathematicalFunctionWrapper wrapper2 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{2,3}),
				Arrays.asList(new Integer[]{1}));
		ProductFunctionWrapper pfw = new ProductFunctionWrapper(wrapper1, wrapper2);
		Matrix parms = new Matrix(4,1);
		parms.setValueAt(0, 0, 0);
		parms.setValueAt(1, 0, 1);
		parms.setValueAt(2, 0, 2);
		parms.setValueAt(3, 0, 3);
		pfw.setParameters(parms);

		Matrix m = new Matrix(1,2);
		m.setValueAt(0, 1, 1);
		
		pfw.setVariables(m);
		GaussianFunction gf1 = new GaussianFunction();
		gf1.setVariableValue(0, 0);
		GaussianFunction gf2 = new GaussianFunction(2,3);
		gf2.setVariableValue(0, 1);
		
		Matrix expectedGradient = gf1.getGradient().scalarMultiply(gf2.getValue()).matrixStack(gf2.getGradient().scalarMultiply(gf1.getValue()), true);
		Matrix observedGradient = pfw.getGradient();
		boolean areDifferent = observedGradient.subtract(expectedGradient).getAbsoluteValue().anyElementLargerThan(1E-8);
		Assert.assertTrue("Comparing gradients", !areDifferent);
	}

	@Test
	public void simpleGradientTestWithThreeIndependentFunctions() {
		InternalMathematicalFunctionWrapper wrapper1 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{0,1}),
				Arrays.asList(new Integer[]{0}));
		InternalMathematicalFunctionWrapper wrapper2 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{2,3}),
				Arrays.asList(new Integer[]{1}));
		InternalMathematicalFunctionWrapper wrapper3 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{4,5}),
				Arrays.asList(new Integer[]{2}));

		ProductFunctionWrapper pfw = new ProductFunctionWrapper(wrapper1, wrapper2, wrapper3);
		Matrix parms = new Matrix(6,1);
		parms.setValueAt(0, 0, 0);
		parms.setValueAt(1, 0, 1);
		parms.setValueAt(2, 0, 2);
		parms.setValueAt(3, 0, 3);
		parms.setValueAt(4, 0, -1);
		parms.setValueAt(5, 0, 2);
		pfw.setParameters(parms);

		Matrix m = new Matrix(1,3);
		m.setValueAt(0, 1, 1);
		m.setValueAt(0, 2, 0);
		pfw.setVariables(m);
		
		GaussianFunction gf1 = new GaussianFunction();
		gf1.setVariableValue(0, 0);
		GaussianFunction gf2 = new GaussianFunction(2,3);
		gf2.setVariableValue(0, 1);
		GaussianFunction gf3 = new GaussianFunction(-1,2);
		gf3.setVariableValue(0, 0);
		
		Matrix expectedGradient = gf1.getGradient().scalarMultiply(gf2.getValue()*gf3.getValue()).matrixStack(
				gf2.getGradient().scalarMultiply(gf1.getValue()*gf3.getValue()), true).matrixStack(
						gf3.getGradient().scalarMultiply(gf1.getValue()*gf2.getValue()), true);
		
		Matrix observedGradient = pfw.getGradient();
		boolean areDifferent = observedGradient.subtract(expectedGradient).getAbsoluteValue().anyElementLargerThan(1E-8);
		Assert.assertTrue("Comparing gradients", !areDifferent);
	}

	@Test
	public void simpleHessianTestWithTwoIndependentFunctions() {
		InternalMathematicalFunctionWrapper wrapper1 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{0,1}),
				Arrays.asList(new Integer[]{0}));
		InternalMathematicalFunctionWrapper wrapper2 = new InternalMathematicalFunctionWrapper(new GaussianFunction(), 
				Arrays.asList(new Integer[]{2,3}),
				Arrays.asList(new Integer[]{1}));
		ProductFunctionWrapper pfw = new ProductFunctionWrapper(wrapper1, wrapper2);
		Matrix parms = new Matrix(4,1);
		parms.setValueAt(0, 0, 0);
		parms.setValueAt(1, 0, 1);
		parms.setValueAt(2, 0, 2);
		parms.setValueAt(3, 0, 3);
		pfw.setParameters(parms);

		Matrix m = new Matrix(1,2);
		m.setValueAt(0, 1, 1);
		
		pfw.setVariables(m);
		GaussianFunction gf1 = new GaussianFunction();
		gf1.setVariableValue(0, 0);
		GaussianFunction gf2 = new GaussianFunction(2,3);
		gf2.setVariableValue(0, 1);
		
		Matrix h11 = gf1.getHessian().scalarMultiply(gf2.getValue());
		Matrix h22 = gf2.getHessian().scalarMultiply(gf1.getValue());
		Matrix h12 = gf1.getGradient().matrixStack(new Matrix(2,1), true).multiply(
				new Matrix(2,1).matrixStack(gf2.getGradient(), true).transpose());
		Matrix expectedHessian = h11.matrixDiagBlock(h22).add(h12).add(h12.transpose());
		Matrix observedHessian = pfw.getHessian();
		boolean areDifferent = observedHessian.subtract(expectedHessian).getAbsoluteValue().anyElementLargerThan(1E-8);
		Assert.assertTrue("Comparing gradients", !areDifferent);
	}

}
