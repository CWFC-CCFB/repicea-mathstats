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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * The GaussQuadrature class is the basic class for all the Gaussian quadrature
 * methods, such as the Gauss-Hermite and the Gauss-Legendre quadrature.
 * @author Mathieu Fortin- October 2012
 */
@SuppressWarnings("serial")
public abstract class AbstractGaussQuadrature extends AbstractNumericalIntegrationMethod {
	
	/**
	 * The enum NumberOfPoints defines the number of quadrature points.
	 * Theoretically, the larger the number of points, the more accurate
	 * the approximation is.
	 * @author Mathieu Fortin - October 2012
	 */
	public static enum NumberOfPoints {N2, N3, N4, N5, N10, N15}

	/**
	 * This method returns the correct number of nodes in the correct order.
	 * @param nodes the original set of QuadratureNode instances
	 * @return a List of QuadratureNode instances
	 */
	protected List<QuadratureNode> getOrderedNodes(Set<QuadratureNode> nodes) {
		List<QuadratureNode> orderedNodes = new ArrayList<QuadratureNode>();
		for (QuadratureNode node : nodes) {
			if (node.getValue() != 0) {
				orderedNodes.add(new QuadratureNode(-node.getValue(), node.getWeight()));
			} 
			orderedNodes.add(node);
		}
		Collections.sort(orderedNodes);
		return orderedNodes;
	}
}
