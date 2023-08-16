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
package repicea.math.integral;

/**
 * This class represents the implementation of a Gaussian Quadrature Node. Every instance contains
 * a value and a weight.
 * @author Mathieu Fortin - July 2012
 */
class QuadratureNode implements Comparable<Object> {

	private double value;
	private double weight;

	/**
	 * Protected Constructor.
	 * @param value the value of the node
	 * @param weight the weight of the node
	 */
	protected QuadratureNode(double value, double weight) {
		this.value = value;
		this.weight = weight;
	}

	/**
	 * This method returns the value of the node.
	 * @return a double
	 */
	public double getValue() {
		return value;
	}


	/**
	 * This method returns the weight of the node.
	 * @return a double
	 */
	public double getWeight() {
		return weight;
	}

	@Override
	public int compareTo(Object arg0) {
		QuadratureNode otherNode = (QuadratureNode) arg0;
		if (this.getValue() < otherNode.getValue()) {
			return -1;
		} else if (this.getValue() == otherNode.getValue()) {
			return 0;
		} else {
			return 1;
		}
	}
}
