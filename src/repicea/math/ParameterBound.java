/*
 * This file is part of the repicea library.
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
package repicea.math;

import java.io.Serializable;

/**
 * This class defines the bound for a particular parameter.
 * @author Mathieu Fortin - October 2011
 */
@SuppressWarnings("serial")
public class ParameterBound implements Serializable, Cloneable {
	
	private Double lowerBound;
	private Double upperBound;
	
	/**
	 * Constructor. <br>
	 * <br>
	 * The value of the bound should be set to null if there is no bound.
	 * 
	 * @param lowerBound a Double instance
	 * @param upperBound a Double instance
	 */
	public ParameterBound(Double lowerBound, Double upperBound) {
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}
	
	/**
	 * Check and eventually set the parameter value to its bound.
	 * @param parameter the suggested parameter value
	 * @return return either the suggested value or the bound
	 */
	protected double validateParameter(double parameter) {
		if (lowerBound != null && parameter < lowerBound) {
			return lowerBound;
		}
		if (upperBound != null && parameter > upperBound) {
			return upperBound;
		}
		return parameter;
	}

	/**
	 * Check if this value is within the parameter bounds if any.
	 * @param parameter the suggested value
	 * @return a boolean
	 */
	protected boolean isParameterValueValid(double parameter) {
		if (lowerBound != null && parameter < lowerBound) {
			return false;
		}
		if (upperBound != null && parameter > upperBound) {
			return false;
		}
		return true;
	}

	@Override
	public ParameterBound clone() {
		return new ParameterBound(lowerBound, upperBound);
	}
	
}
