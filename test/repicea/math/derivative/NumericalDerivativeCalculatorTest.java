/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, CFS
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
package repicea.math.derivative;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.derivative.NumericalDerivativeCalculator.Technique;
import repicea.stats.model.glm.LinkFunction;
import repicea.stats.model.glm.LinkFunction.Type;

public class NumericalDerivativeCalculatorTest {

	@Test
	public void firstDerivativeWithNewtonQuotient1() {
		LinkFunction lf = new LinkFunction(Type.Logit);
		lf.setVariableValue(0, 1);
		lf.setParameterValue(0, 0);
		double originalValue = lf.getValue();
		double refValue = lf.getGradient().getValueAt(0, 0);
		double obsValue = NumericalDerivativeCalculator.computeNumericalFirstDerivative(Technique.NewtonDifferenceQuotient, lf, 0, true);
		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
	}
	
	@Test
	public void firstDerivativeWithNewtonQuotient2() {
		LinkFunction lf = new LinkFunction(Type.Logit);
		lf.setVariableValue(0, 1);
		lf.setParameterValue(0, 1);
		double originalValue = lf.getValue();
		double refValue = lf.getGradient().getValueAt(0, 0);
		double obsValue = NumericalDerivativeCalculator.computeNumericalFirstDerivative(Technique.NewtonDifferenceQuotient, lf, 0, true);
		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
	}

	@Test
	public void firstDerivativeWithNewtonQuotient3() {
		LinkFunction lf = new LinkFunction(Type.Logit);
		lf.setVariableValue(0, 1);
		lf.setParameterValue(0, -1);
		double originalValue = lf.getValue();
		double refValue = lf.getGradient().getValueAt(0, 0);
		double obsValue = NumericalDerivativeCalculator.computeNumericalFirstDerivative(Technique.NewtonDifferenceQuotient, lf, 0, true);
		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
	}
	
	@Test
	public void firstDerivativeWithNewtonQuotient4() {
		LinkFunction lf = new LinkFunction(Type.CLogLog);
		lf.setVariableValue(0, 1);
		lf.setParameterValue(0, -1);
		double originalValue = lf.getValue();
		double refValue = lf.getGradient().getValueAt(0, 0);
		double obsValue = NumericalDerivativeCalculator.computeNumericalFirstDerivative(Technique.NewtonDifferenceQuotient, lf, 0, true);
		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
	}
	
	@Test
	public void firstDerivativeWithSymmetricQuotient1() {
		LinkFunction lf = new LinkFunction(Type.Logit);
		lf.setVariableValue(0, 1);
		lf.setParameterValue(0, 0);
		double originalValue = lf.getValue();
		double refValue = lf.getGradient().getValueAt(0, 0);
		double obsValue = NumericalDerivativeCalculator.computeNumericalFirstDerivative(Technique.SymmetricDifferenceQuotient, lf, 0, true);
		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
	}
	
	@Test
	public void firstDerivativeWithSymmetricQuotient2() {
		LinkFunction lf = new LinkFunction(Type.Logit);
		lf.setVariableValue(0, 1);
		lf.setParameterValue(0, 1);
		double originalValue = lf.getValue();
		double refValue = lf.getGradient().getValueAt(0, 0);
		double obsValue = NumericalDerivativeCalculator.computeNumericalFirstDerivative(Technique.SymmetricDifferenceQuotient, lf, 0, true);
		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
	}

	@Test
	public void firstDerivativeWithSymmetricQuotient3() {
		LinkFunction lf = new LinkFunction(Type.Logit);
		lf.setVariableValue(0, 1);
		lf.setParameterValue(0, -1);
		double originalValue = lf.getValue();
		double refValue = lf.getGradient().getValueAt(0, 0);
		double obsValue = NumericalDerivativeCalculator.computeNumericalFirstDerivative(Technique.SymmetricDifferenceQuotient, lf, 0, true);
		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
	}
	
	@Test
	public void firstDerivativeWithSymmetricQuotient4() {
		LinkFunction lf = new LinkFunction(Type.CLogLog);
		lf.setVariableValue(0, 1);
		lf.setParameterValue(0, -1);
		double originalValue = lf.getValue();
		double refValue = lf.getGradient().getValueAt(0, 0);
		double obsValue = NumericalDerivativeCalculator.computeNumericalFirstDerivative(Technique.SymmetricDifferenceQuotient, lf, 0, true);
		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
	}

	
	
	
//	@Test
//	public void secondtDerivativeWithNewtonQuotient1() {
//		LinkFunction lf = new LinkFunction(Type.Logit);
//		lf.setVariableValue(0, 1);
//		lf.setParameterValue(0, 0);
//		double originalValue = lf.getValue();
//		double refValue = lf.getHessian().getValueAt(0, 0);
//		double obsValue = NumericalDerivativeCalculator.computeNumericalSecondDerivative(lf, 0, true);
//		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
//		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
//		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
//	}
//	
//	@Test
//	public void secondtDerivativeWithNewtonQuotient2() {
//		LinkFunction lf = new LinkFunction(Type.Logit);
//		lf.setVariableValue(0, 1);
//		lf.setParameterValue(0, 1);
//		double originalValue = lf.getValue();
//		double refValue = lf.getHessian().getValueAt(0, 0);
//		double obsValue = NumericalDerivativeCalculator.computeNumericalSecondDerivative(lf, 0, true);
//		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
//		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
//		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
//	}
//
//	@Test
//	public void secondtDerivativeWithNewtonQuotient3() {
//		LinkFunction lf = new LinkFunction(Type.Logit);
//		lf.setVariableValue(0, 1);
//		lf.setParameterValue(0, -1);
//		double originalValue = lf.getValue();
//		double refValue = lf.getHessian().getValueAt(0, 0);
//		double obsValue = NumericalDerivativeCalculator.computeNumericalSecondDerivative(lf, 0, true);
//		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
//		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
//		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
//	}
//	
//	@Test
//	public void secondtDerivativeWithNewtonQuotient4() {
//		LinkFunction lf = new LinkFunction(Type.CLogLog);
//		lf.setVariableValue(0, 1);
//		lf.setParameterValue(0, -1);
//		double originalValue = lf.getValue();
//		double refValue = lf.getHessian().getValueAt(0, 0);
//		double obsValue = NumericalDerivativeCalculator.computeNumericalSecondDerivative(lf, 0, true);
//		System.out.println("Observed is " + obsValue + "; Expected was " + refValue);
//		Assert.assertEquals("Testing numerical derivative value", refValue, obsValue, 1E-8);
//		Assert.assertEquals("Testing original value", originalValue, lf.getValue(), 1E-8);
//	}

	
}
