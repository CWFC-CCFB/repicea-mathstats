/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.integral;

import java.security.InvalidParameterException;
import java.util.List;

import repicea.math.EvaluableFunction;
import repicea.math.Matrix;

/**
 * This class implements the Composite Simpson's rule.
 * @author Mathieu Fortin - July 2012
 */
@SuppressWarnings("serial")
public final class CompositeSimpsonRule extends AbstractNumericalIntegrationMethod implements UnidimensionalIntegralApproximation<EvaluableFunction<Double>>,
																								UnidimensionalIntegralApproximationForMatrix<EvaluableFunction<Matrix>> {

	private static final double EPSILON = 1E-12;
	private final int numberOfSubintervals;


	/**
	 * Constructor.
	 * @param numberOfSubintervals the number of sub intervals which must be even
	 */
	public CompositeSimpsonRule(int numberOfSubintervals) {
		if ((numberOfSubintervals % 2) != 0) {
			throw new InvalidParameterException("The number of subintervals must be even!");
		} else {
			this.numberOfSubintervals = numberOfSubintervals;
		}
	}

	/**
	 * Constructor 2.
	 * @param points a List of evenly spaced points.
	 */
	public CompositeSimpsonRule(List<Double> points) {
		if (points.size() % 2 != 1) {		// checks if the number of subsections is even
			throw new InvalidParameterException("The number of subintervals must be even!");
		}
		double previousDiff = 0;
		double diff;
		for (int i = 1; i < points.size(); i++) {
			diff = points.get(i) - points.get(i - 1);
			if (i == 1) {
				previousDiff = diff;
			}
			if (Math.abs(1 - diff / previousDiff) > EPSILON) {		// checks if the sections are evenly spaced
				throw new InvalidParameterException("The points are not evenly spaced!");
			}
		}
		setXValuesFromListOfPoints(points);
		this.numberOfSubintervals = points.size() - 1;

	}

	@Override
	public List<Double> getWeights() {
		if (weights.isEmpty()) {
			for (int i = 0; i < getXValues().size(); i++) {
				if (i == 0 || i == getXValues().size() - 1) {
					weights.add(1d);
				} else if ((i%2) == 1) {
					weights.add(4d);
				} else {
					weights.add(2d);
				}
			}
		}
		return weights;
	}

	private void setXValues() {
		if (Double.isNaN(getLowerBound()) || Double.isNaN(getUpperBound())) {
			return;
		} else {
			double difference = getUpperBound() - getLowerBound();
			double slope = difference / numberOfSubintervals;
			double height;
			int i = 0;
			while (i <= numberOfSubintervals) {
				height = getLowerBound() + i * slope;
				xValues.add(height);
				i++;
			}
		}
	}



	@Override
	public List<Double> getXValues() {
		if (xValues.isEmpty()) {
			setXValues();
		}
		return xValues;
	}

	@Override
	public List<Double> getRescalingFactors() {
		if (rescalingFactors.isEmpty()) {
			double rf = (getUpperBound() - getLowerBound()) / numberOfSubintervals / 3;
			for (int i = 0; i < getXValues().size(); i++) {
				rescalingFactors.add(rf);
			}
		}
		return rescalingFactors;
	}

	@Override
	public double getIntegralApproximation(EvaluableFunction<Double> functionToEvaluate, 
			int index,
			boolean isParameter) {

		double originalValue;
		if (isParameter) {
			originalValue = functionToEvaluate.getParameterValue(index);
		} else {
			originalValue = functionToEvaluate.getVariableValue(index);
		}

		double sum = 0d;
		double point;
		for (int i = 0; i < getXValues().size(); i++) {
			point = getXValues().get(i);
			if (isParameter) {
				functionToEvaluate.setParameterValue(index, point);
			} else {
				functionToEvaluate.setVariableValue(index, point);
			}
			sum += functionToEvaluate.getValue() * getWeights().get(i) * getRescalingFactors().get(i);
		}
		
		if (isParameter) {
			functionToEvaluate.setParameterValue(index, originalValue);
		} else {
			functionToEvaluate.setVariableValue(index, originalValue);
		}

		return sum;
	}

	@Override
	public Matrix getIntegralApproximationForMatrixFunction(EvaluableFunction<Matrix> functionToEvaluate, 
											int index,
											boolean isParameter) {
		double originalValue;
		if (isParameter) {
			originalValue = functionToEvaluate.getParameterValue(index);
		} else {
			originalValue = functionToEvaluate.getVariableValue(index);
		}
		
		Matrix sum = null;
		double point;
		for (int i = 0; i < getXValues().size(); i++) {
			point = getXValues().get(i);
			if (isParameter) {
				functionToEvaluate.setParameterValue(index, point);
			} else {
				functionToEvaluate.setVariableValue(index, point);
			}
			Matrix value = functionToEvaluate.getValue().scalarMultiply(getWeights().get(i) * getRescalingFactors().get(i));
			sum = i == 0 ? value : sum.add(value);
		}
		
		if (isParameter) {
			functionToEvaluate.setParameterValue(index, originalValue);
		} else {
			functionToEvaluate.setVariableValue(index, originalValue);
		}
		
		return sum;
	}
	
}
