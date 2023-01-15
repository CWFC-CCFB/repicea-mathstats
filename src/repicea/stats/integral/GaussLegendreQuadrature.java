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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import repicea.math.EvaluableFunction;
import repicea.math.Matrix;

/**
 * The GaussLegendreQuadrature class implements a numerical integration method based
 * on Legendre polynomials. The current implementation is based on 2, 3, 4 or 5 points,
 * which is set in the constructor.
 * @author Mathieu Fortin - July 2012
 */
@SuppressWarnings("serial")
public class GaussLegendreQuadrature extends AbstractGaussQuadrature implements UnidimensionalIntegralApproximation<EvaluableFunction<Double>>,
																				UnidimensionalIntegralApproximationForMatrix<EvaluableFunction<Matrix>> {

	private static Map<NumberOfPoints, Set<QuadratureNode>> NODE_MAP = new HashMap<NumberOfPoints, Set<QuadratureNode>>();
	static {
		Set<QuadratureNode> nodes = new HashSet<QuadratureNode>();
		nodes.add(new QuadratureNode(Math.sqrt(3d) / 3d, 1d));
		NODE_MAP.put(NumberOfPoints.N2, nodes);
		
		nodes = new HashSet<QuadratureNode>();
		nodes.add(new QuadratureNode(0d, 8d / 9d));
		nodes.add(new QuadratureNode(Math.sqrt(15d) / 5d, 5d / 9d));
		NODE_MAP.put(NumberOfPoints.N3, nodes);
		
		nodes = new HashSet<QuadratureNode>();
		nodes.add(new QuadratureNode(Math.sqrt(525d - 70d * Math.sqrt(30)) / 35d, (18d + Math.sqrt(30d)) / 36d));
		nodes.add(new QuadratureNode(Math.sqrt(525d + 70d * Math.sqrt(30)) / 35d, (18d - Math.sqrt(30d)) / 36d));
		NODE_MAP.put(NumberOfPoints.N4, nodes);

		nodes = new HashSet<QuadratureNode>();
		nodes.add(new QuadratureNode(0d, 128d / 225d));
		nodes.add(new QuadratureNode(Math.sqrt(245d - 14d * Math.sqrt(70)) / 21d, (322d + 13 * Math.sqrt(70d)) / 900d));
		nodes.add(new QuadratureNode(Math.sqrt(245d + 14d * Math.sqrt(70)) / 21d, (322d - 13 * Math.sqrt(70d)) / 900d));
		NODE_MAP.put(NumberOfPoints.N5, nodes);
	}
	
	private NumberOfPoints numberOfPoints;
	
	/**
	 * Constructor.
	 * @param numberOfPoints a NumberOfPoints enum variable (either NumberOfPoints.N5, NumberOfPoints.N10, or NumberOfPoints.N15) 
	 */
	public GaussLegendreQuadrature(NumberOfPoints numberOfPoints) {
		if (!NODE_MAP.containsKey(numberOfPoints)) {
			throw new InvalidParameterException("The Gauss-Legendre quadrature with this number of points is not implemented!");
		}
		this.numberOfPoints = numberOfPoints;
		setLowerBound(-1);
		setUpperBound(1);
	}
	
	@Override
	public List<Double> getWeights() {
		if (weights == null) {
			getXValues();
		}
		return weights;
	}



	@Override
	public List<Double> getXValues() {
		if (xValues.isEmpty()) {
			weights.clear();
			List<QuadratureNode> orderedNodes = getOrderedNodes(GaussLegendreQuadrature.NODE_MAP.get(numberOfPoints));
			double intercept = (getLowerBound() + getUpperBound()) * .5;
			double slope = (getUpperBound() - getLowerBound()) * .5;
			for (QuadratureNode node : orderedNodes) {
				xValues.add(node.getValue() * slope + intercept);
				weights.add(node.getWeight());
			}
		}
		return xValues;
	}



	@Override
	public List<Double> getRescalingFactors() {
		if (rescalingFactors.isEmpty()) {
			List<Double> xValues = getXValues();
			double rescaling = (getUpperBound() - getLowerBound()) * .5;
			for (int i = 0; i < xValues.size(); i++) {
				rescalingFactors.add(rescaling);
			}
		}
		return rescalingFactors;
	}

	@Override
	public double getIntegralApproximation(EvaluableFunction<Double> functionToEvaluate, int index,
			boolean isParameter) {
		
		double originalValue;
		if (isParameter) {
			originalValue = functionToEvaluate.getParameterValue(index);
		} else {
			originalValue = functionToEvaluate.getVariableValue(index);
		}

		double sum = 0;
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
