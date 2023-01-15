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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import repicea.serial.xml.XmlSerializerChangeMonitor;
import repicea.stats.integral.GaussHermiteQuadrature.GaussHermiteQuadratureCompatibleFunction;

@SuppressWarnings("serial")
public abstract class AbstractGaussHermiteQuadrature extends AbstractGaussQuadrature {

	static {
		XmlSerializerChangeMonitor.registerClassNameChange("repicea.stats.integral.GaussQuadrature$NumberOfPoints", "repicea.stats.integral.AbstractGaussQuadrature$NumberOfPoints");
	}

	protected static Map<NumberOfPoints, Set<QuadratureNode>> NODE_MAP = new HashMap<NumberOfPoints, Set<QuadratureNode>>();
	static {
		Set<QuadratureNode> nodes = new HashSet<QuadratureNode>();
		nodes.add(new QuadratureNode(0d, 0.945308720482942));
		nodes.add(new QuadratureNode(0.958572464613819, 0.393619323152241));
		nodes.add(new QuadratureNode(0.202018287045609E1, 0.199532420590459E-1));
		NODE_MAP.put(NumberOfPoints.N5, nodes);
		
		nodes = new HashSet<QuadratureNode>();
		nodes.add(new QuadratureNode(0.342901327223704609, 0.610862633735325799));
		nodes.add(new QuadratureNode(0.103661082978951365E1, 0.240138611082314686));
		nodes.add(new QuadratureNode(0.175668364929988177E1, 0.338743944554810631E-1));
		nodes.add(new QuadratureNode(0.253273167423278980E1, 0.134364574678123269E-2));
		nodes.add(new QuadratureNode(0.343615911883773760E1, 0.764043285523262063E-5));
		NODE_MAP.put(NumberOfPoints.N10, nodes);
		
		nodes = new HashSet<QuadratureNode>();
		nodes.add(new QuadratureNode(0d, 0.564100308726417532853));
		nodes.add(new QuadratureNode(0.565069583255575748526, 0.412028687498898627026));
		nodes.add(new QuadratureNode(0.113611558521092066632E1, 0.158488915795935746884));
		nodes.add(new QuadratureNode(0.171999257518648893242E1, 0.307800338725460822287E-1));
		nodes.add(new QuadratureNode(0.232573248617385774545E1, 0.277806884291277589608E-2));
		nodes.add(new QuadratureNode(0.296716692790560324849E1, 0.100004441232499868127E-3));
		nodes.add(new QuadratureNode(0.366995037340445253473E1, 0.105911554771106663578E-5));
		nodes.add(new QuadratureNode(0.449999070730939155366E1, 0.152247580425351702016E-8));
		NODE_MAP.put(NumberOfPoints.N15, nodes);
	}

	
	protected final NumberOfPoints numberOfPoints;
	
	/**
	 * Constructor.
	 * @param numberOfPoints a NumberOfPoints enum variable (either NumberOfPoints.N5, NumberOfPoints.N10, or NumberOfPoints.N15) 
	 */
	protected AbstractGaussHermiteQuadrature(NumberOfPoints numberOfPoints) {
		if (!NODE_MAP.containsKey(numberOfPoints)) {
			throw new InvalidParameterException("The Gauss-Hermite quadrature with this number of points is not implemented!");
		}
		this.numberOfPoints = numberOfPoints;
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
			List<QuadratureNode> orderedNodes = getOrderedNodes(AbstractGaussHermiteQuadrature.NODE_MAP.get(numberOfPoints));
			for (QuadratureNode node : orderedNodes) {
				xValues.add(node.getValue());
				weights.add(node.getWeight());
			}
		}
		return xValues;
	}



	@Override
	public List<Double> getRescalingFactors() {
		if (rescalingFactors.isEmpty()) {
			List<Double> xValues = getXValues();
			for (int i = 0; i < xValues.size(); i++) {
				rescalingFactors.add(1d);
			}
		}
		return rescalingFactors;
	}

	/**
	 * This method makes it possible to integrate an AbstractStatisticalExpression through Gauss-Hermite quadrature. 
	 * @param functionToEvaluate an EvaluableFunction instance that returns Double 
	 * @param index the index of the parameter over which the integration is made
	 * @param isParameter a boolean to indicate that indices refer to parameters. If false, it is assumed that the
	 * indices refer to variables.
	 * @param standardDeviation the standard deviation of this variable
	 * @return the approximation of the integral
	 */
	protected abstract double getOneDimensionIntegral(GaussHermiteQuadratureCompatibleFunction<Double> functionToEvaluate,
			List<Integer> indices, 
			boolean isParameter,
			int startingIndex);
	
	/**
	 * This method returns the value of a multi-dimension integral
	 * @param functionToEvaluate an EvaluableFunction instance that returns Double 
	 * @param indices the indices of the parameters over which the integration is made
	 * @param isParameter a boolean to indicate that indices refer to parameters. If false, it is assumed that the
	 * indices refer to variables.
	 * @param lowerCholeskyTriangle the lower triangle of the Cholesky factorization of the variance-covariance matrix
	 * @return the approximation of the integral
	 */
	protected double getMultiDimensionIntegral(GaussHermiteQuadratureCompatibleFunction<Double> functionToEvaluate,
												List<Integer> indices, 
												boolean isParameter,
												int startingIndex) {
		if (startingIndex == indices.size() - 1) {
			return getOneDimensionIntegral(functionToEvaluate, indices, isParameter, startingIndex);
		}
		Integer thisIndex = indices.get(startingIndex);
		double originalValue = isParameter ? functionToEvaluate.getParameterValue(thisIndex) : functionToEvaluate.getVariableValue(thisIndex);
		double sum = 0;
		double value;
		for (int i = 0; i < getXValues().size(); i++) {
			double thisValueOnTheOriginalScale = functionToEvaluate.convertFromGaussToOriginal(getXValues().get(i), originalValue, startingIndex, startingIndex); // the first variance in the vcov matrix
			if (isParameter) {
				functionToEvaluate.setParameterValue(thisIndex, thisValueOnTheOriginalScale);
			} else {
				functionToEvaluate.setVariableValue(thisIndex, thisValueOnTheOriginalScale);
			}
			Map<Integer, Double> originalValues = new HashMap<Integer, Double>();
			for (int ii = startingIndex + 1; ii < indices.size(); ii++) {
				int index = indices.get(ii);
				double original = isParameter ? functionToEvaluate.getParameterValue(index) : functionToEvaluate.getVariableValue(index);
				originalValues.put(ii, original);
				double valueOnTheOriginalScale = functionToEvaluate.convertFromGaussToOriginal(getXValues().get(i), original, ii, startingIndex);
				if (isParameter) {
					functionToEvaluate.setParameterValue(index, valueOnTheOriginalScale);
				} else {
					functionToEvaluate.setVariableValue(index, valueOnTheOriginalScale);
				}

			}
			int newStartingIndex = startingIndex + 1;
			value = getMultiDimensionIntegral(functionToEvaluate,
					indices, 
					isParameter,
					newStartingIndex) * getWeights().get(i);
			sum += value;
			for (int ii = startingIndex + 1; ii < indices.size(); ii++) {
				if (isParameter) {
					functionToEvaluate.setParameterValue(indices.get(ii), originalValues.get(ii));
				} else {
					functionToEvaluate.setVariableValue(indices.get(ii), originalValues.get(ii));
				}
			}
		}
		if (isParameter) {
			functionToEvaluate.setParameterValue(thisIndex, originalValue);
		} else {
			functionToEvaluate.setVariableValue(thisIndex, originalValue);
		}
		return sum;
	}

}
