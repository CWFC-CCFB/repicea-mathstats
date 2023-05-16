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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The NumericalIntegrationMethod class is the basic class for all numerical integration, such as
 * Gauss quadrature, trapezoidal rule, Simpson's composite rule, etc...
 * @author Mathieu Fortin - October 2012
 */
@SuppressWarnings("serial")
public abstract class AbstractNumericalIntegrationMethod implements Serializable {
	
	private double upperBound = Double.NaN;
	private double lowerBound = Double.NaN;
	protected final List<Double> xValues;
	protected final List<Double> weights;
	protected final List<Double> rescalingFactors;

	protected AbstractNumericalIntegrationMethod() {
		xValues = new ArrayList<Double>();
		weights = new ArrayList<Double>();
		rescalingFactors = new ArrayList<Double>();
	}
	
	/**
	 * This method sets the upper and lower bounds from a list of Double instances. It also sets the
	 * xValues member and sorts it in ascending order. NOTE: replicates are omitted in the xValues 
	 * member. If points = {1d, 1d, 2d, 3d}, xValues will be set to {1d, 2d, 3d}.
	 * @param points a List of Double instances
	 */
	protected void setXValuesFromListOfPoints(List<Double> points) {
		setUpperAndLowerBoundsFromListOfDouble(points);
		xValues.clear();
		for (Double point : points) {
			if (!xValues.contains(point)) {
				xValues.add(point);
			}
		}
		Collections.sort(xValues);
	}
	
	
	/**
	 * This method finds the minimum and the maximum values in a list of Double instances. The 
	 * lower and the upper bounds are then set according to these values.
	 * @param points a List od Double instance
	 */
	private void setUpperAndLowerBoundsFromListOfDouble(List<Double> points) {
		double lowerLimit = Double.NaN;
		double upperLimit = Double.NaN; 
		for (Double point : points) {
			if (Double.isNaN(upperLimit)) {
				upperLimit = point;
			}
			if (Double.isNaN(lowerLimit)) {
				lowerLimit = point;
			}
			if (point < lowerLimit) {
				lowerLimit = point;
			}
			if (point > upperLimit) {
				upperLimit = point;
			}
		}
		setLowerBound(lowerLimit);
		setUpperBound(upperLimit);
	}
	
	/**
	 * Set the upper bound of the numerical integration. <p>
	 * Setting the upper bound to Double.NaN (default value) means plus infinity. IMPORTANT:
	 * changing the upper bound resets the arrays that contain the x values, the weights,
	 * and the rescaling factors.
	 * @param upperBound the value of the upper bound
	 */
	public void setUpperBound(double upperBound) {
		this.upperBound = upperBound; 
		resetArrays();
	}
	
	/**
	 * This method sets the lower bound of the numerical integration. Setting the 
	 * lower bound to Double.NaN (default value) means minus infinity. IMPORTANT:
	 * changing the upper bound resets the arrays that contain the x values, the weights,
	 * and the rescaling factors.
	 * @param lowerBound a double
	 */
	public void setLowerBound(double lowerBound) {
		this.lowerBound = lowerBound;
		resetArrays();
	}

	private void resetArrays() {
		xValues.clear();
		weights.clear();
		rescalingFactors.clear();
	}
	
	
	/**
	 * This method returns the weights associated to the numerical integration.
	 * @return a List of Double instances
	 */
	public abstract List<Double> getWeights();

	/**
	 * This method returns the x values for the numerical integration.
	 * @return a List of Double instances
	 */
	public abstract List<Double> getXValues();
	
	/**
	 * This method returns the rescaling factor for the numerical integration.
	 * @return a List of Double instances
	 */
	public abstract List<Double> getRescalingFactors();

	/**
	 * This method returns the lower bound of the numerical integration.
	 * @return a double or Double.NaN in case of minus infinity
	 */
	public double getLowerBound() {return lowerBound;}
	
	/**
	 * This method returns the upper bound of the numerical integration.
	 * @return a double or Double.NaN in case of plus infinity
	 */
	public double getUpperBound() {return upperBound;}
	
}
