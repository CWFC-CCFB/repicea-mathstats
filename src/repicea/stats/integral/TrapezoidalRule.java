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
package repicea.stats.integral;

import java.security.InvalidParameterException;
import java.util.List;

import repicea.math.EvaluableFunction;
import repicea.math.Matrix;
import repicea.math.utility.MatrixUtility;

/**
 * The TrapezoidalRule class implements the trapezoidal rule integration method.<p>
 * The resolution parameter in the constructor sets the distance between the x values.
 * The class is adapted to resolution that does not perfectly match the range of the integral.
 * @author Mathieu Fortin - July 2012
 */
@SuppressWarnings("serial")
public class TrapezoidalRule extends AbstractNumericalIntegrationMethod implements UnidimensionalIntegralApproximation<EvaluableFunction<Double>>,
																					UnidimensionalIntegralApproximationForMatrix<EvaluableFunction<Matrix>> {

	private double resolution;
	
	/**
	 * Constructor.
	 * @param resolution the distance between the x points (must be larger than 0)
	 */
	public TrapezoidalRule(double resolution) {
		if (resolution <= 0) {
			throw new InvalidParameterException("The resolution must be larger than 0!");
		}
		this.resolution = resolution;
	}
	
	
	/**
	 * This constructor is set through a list of double instances that represents the point along the
	 * integral. Any replicate is eliminated in the process.
	 * @param points a List of Double instances
	 */
	public TrapezoidalRule(List<Double> points) {
		setXValuesFromListOfPoints(points);
	}

	/**
	 * Constructor for derived classes.
	 */
	protected TrapezoidalRule() {
		super();
	}

	private void setXValues() {
		if (Double.isNaN(getLowerBound()) || Double.isNaN(getUpperBound())) {
			return;
		} else {
			double difference = getUpperBound() - getLowerBound();
			double numberOfSectionsWithDigits = difference / resolution;
			int numberOfSections = (int) (numberOfSectionsWithDigits);
			if (numberOfSectionsWithDigits > numberOfSections) {
				numberOfSections++;			// add one section for the top section
			}
			numberOfSections++;
			double height;
			int i = 0;
			while (i < numberOfSections) {
				height = getLowerBound() + i * resolution;
				if (height >= getUpperBound()) {
					height = getUpperBound();
				}
				xValues.add(height);
				i++;
			}
		}
	}

	/*
	 * Calculates the weights based on the distance between the sections. 
	 * The distances does not need to be even. 
	 * (non-Javadoc)
	 * @see repicea.stats.integral.NumericalIntegrationMethod#getWeights()
	 */
	@Override
	public List<Double> getWeights() {
		if (weights.isEmpty()) {
			List<Double> xValues = getXValues();
			if (xValues != null) {
				double diff;
				double lag = 0d;
				for (int i = 0; i < xValues.size(); i++) {
					if (i == xValues.size() - 1) {
						diff = 0;
					} else {
						diff = (xValues.get(i + 1) - xValues.get(i)) * .5;
					}
					weights.add(diff + lag);
					lag = diff;
				}
			}
		}
		return weights;
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
			List<Double> xValues = getXValues();
			if (!xValues.isEmpty()) {
				for (int i = 0; i < xValues.size(); i++) {
					rescalingFactors.add(1d);
				}
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
			if (i == 0) {
				sum = value;
			} else {
				sum = sum.add(value);
			}
		}
		
		if (isParameter) {
			functionToEvaluate.setParameterValue(index, originalValue);
		} else {
			functionToEvaluate.setVariableValue(index, originalValue);
		}
		
		return sum;
	}

//	public static void main(String[] args) {
//		TrapezoidalRule tr = new TrapezoidalRule(3);
//		tr.setLowerBound(0);
//		tr.setUpperBound(20);
//		List<Double> xValues = tr.getXValues();
//		List<Double> weights = tr.getWeights();
//		List<Double> mf = tr.getRescalingFactors();
//		int u = 0;
//	}
	
}
