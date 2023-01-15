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

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import repicea.math.AbstractMathematicalFunction;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

@SuppressWarnings("serial")
public class EmpiricalDistributionProbabilityDensityFunction extends AbstractMathematicalFunction {

	private static class ValueWeight implements Comparable<ValueWeight> {

		final Double x;
		double w;
		double F;
		double f;
		double f_prime;

		private ValueWeight(double v, double w) {
			this.x = v;
			this.w = w;
		}
		
		@Override
		public int compareTo(ValueWeight o) {
			return x.compareTo(o.x);
		}
		
	}
	
	private final double K_alpha_0_05 = 1.358;
	
	private final List<ValueWeight> valueWeights;
	private final double lambda;
	private final double lowerBound;
	private final double upperBound;
	
	public EmpiricalDistributionProbabilityDensityFunction(List<Double> values, List<Double> weights){
		valueWeights = new ArrayList<ValueWeight>();
		if (weights != null && weights.size() != values.size()) {
			throw new InvalidParameterException("The weights argument should be a list of the same size as that of the values argument!");
		}
		List<ValueWeight> tmpValueWeights = new ArrayList<ValueWeight>();
		for (int i = 0; i < values.size(); i++) {
			double w = (weights != null) ? weights.get(i) : 1d; 
			tmpValueWeights.add(new ValueWeight(values.get(i), w));
		}
		Collections.sort(tmpValueWeights);
		for (int i = 0; i < tmpValueWeights.size(); i++) {
			ValueWeight lastValueWeight = valueWeights.isEmpty() ? null : valueWeights.get(valueWeights.size() - 1);
			ValueWeight thisValueWeight = tmpValueWeights.get(i);
			if (lastValueWeight != null && Math.abs(lastValueWeight.x - thisValueWeight.x) < 1E-8) {
				lastValueWeight.w += thisValueWeight.w;
			} else {
				valueWeights.add(thisValueWeight);
			}
		}
		double S_N = 0d;
		for (ValueWeight vw : valueWeights) {
			S_N += vw.w;
		}
		double S_n = 0d;
		double max_f_prime = Double.NEGATIVE_INFINITY;
		double max_delta = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < valueWeights.size(); i++) {
			ValueWeight previousVw  = i == 0 ? null : valueWeights.get(i - 1);
			ValueWeight vw  = valueWeights.get(i);
			
			S_n += vw.w;
			vw.F = S_n / S_N;
			vw.f = previousVw == null ? 0 : (vw.F - previousVw.F) / (vw.x - previousVw.x);
			vw.f_prime = previousVw == null ? 0 : (vw.f - previousVw.f) / (vw.x - previousVw.x);
			double delta_x = previousVw == null ? Double.NEGATIVE_INFINITY : (vw.x - previousVw.x);
			if (vw.f_prime > max_f_prime)
				max_f_prime = vw.f_prime; 
			if (delta_x > max_delta)
				max_delta = delta_x;
		}
		double lambdaTmp = Math.sqrt(2 * K_alpha_0_05 / max_f_prime) * Math.pow(S_N, -0.25);
		lambda = lambdaTmp < max_delta ? max_delta : lambdaTmp;
		lowerBound = valueWeights.get(0).x;
		upperBound = valueWeights.get(valueWeights.size() - 1).x;
		setVariableValue(0, valueWeights.get(0).x);
	}
	
	private double findFEqualToOrSmallerThan(double x) {
		int lower = 0;
		int upper = valueWeights.size() - 1;
		if (x < valueWeights.get(0).x) {
			return 0d;
		} else if (x >= valueWeights.get(upper).x) {
			return 1d;
		} else {
			return findBetweenTheseBounds(x, lower, upper);
		}
	}
	
	private double findBetweenTheseBounds(double x, int lower, int upper) {
		if (x >= valueWeights.get(lower).x && x < valueWeights.get(lower + 1).x) {
			return valueWeights.get(lower).F;
		} else if (x >= valueWeights.get(upper - 1).x && x < valueWeights.get(upper).x) {
			return valueWeights.get(upper - 1).F;
		} else {
			int middle = (int) Math.floor(0.5 * (upper + lower));
			return x < valueWeights.get(middle).x ? 
					findBetweenTheseBounds(x, lower, middle) :
						findBetweenTheseBounds(x, middle, upper);
		}
	}

	
	@Override
	public Double getValue() {
		double x = getVariableValue(0);
		double F_upper, F_lower, g_x;
		if (x < lowerBound + lambda) {
			if (x < lowerBound) {
				g_x = 0d;
			} else {
	 			F_upper = findFEqualToOrSmallerThan(x + lambda);
				F_lower = findFEqualToOrSmallerThan(lowerBound);
				g_x = (F_upper - F_lower) / (x + lambda - lowerBound);
			}
		} else if (x > upperBound - lambda) {
			if (x > upperBound) {
				g_x = 0;
			} else {
				F_upper = findFEqualToOrSmallerThan(upperBound);
				F_lower = findFEqualToOrSmallerThan(x - lambda);
				g_x = (F_upper - F_lower) / (upperBound + lambda - x);
			}
		} else {
			F_upper = findFEqualToOrSmallerThan(x + lambda);
			F_lower = findFEqualToOrSmallerThan(x - lambda);
			g_x = (F_upper - F_lower) / (2 * lambda);
		}
		return g_x;
	}

	@Override
	public Matrix getGradient() {
		return null;	// this function has no parameter
	}

	@Override
	public SymmetricMatrix getHessian() {
		return null;	// this function has no parameter
	}

	@Override
	public void setParameterValue(int parameterIndex, double parameterValue) {
		throw new UnsupportedOperationException("The EmpiricalDistributionFunction class has no parameter on its own!");
	}

	@Override
	public void setVariableValue(int variableIndex, double variableValue) {
		if (variableIndex > 0) {
			throw new UnsupportedOperationException("The EmpiricalDistributionFunction class has a single variable namely the x observation!");
		} else {
			super.setVariableValue(variableIndex, variableValue);
		}
	}

	/**
	 * Return the lower bound of the empirical distribution.
	 * @return a double
	 */
	public double getLowerBoundForX() {return lowerBound;}
	
	/**
	 * Return the lower bound of the empirical distribution.
	 * @return a double
	 */
	public double getUpperBoundForX() {return upperBound;}
	
}