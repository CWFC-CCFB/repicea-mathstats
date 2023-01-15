/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2016 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.estimates;

import repicea.math.Matrix;
import repicea.stats.distributions.BasicBound;

public class ConfidenceInterval {

	@SuppressWarnings("serial")
	class CIBound extends BasicBound {

		CIBound(boolean isUpperBound) {
			super(isUpperBound);
		}
		
		@Override
		protected void setBoundValue(Matrix value) {
			super.setBoundValue(value);
		}
		
		@Override
		public Matrix getBoundValue() {
			return super.getBoundValue();
		}
	}
	
	private final CIBound lowerBound;
	private final CIBound upperBound;
	private final double probabilityLevel;
	
	protected ConfidenceInterval(Matrix lowerBoundValue, Matrix upperBoundValue, double probabilityLevel) {
		lowerBound = new CIBound(false);
		upperBound = new CIBound(true);
		lowerBound.setBoundValue(lowerBoundValue);
		upperBound.setBoundValue(upperBoundValue);
		this.probabilityLevel = probabilityLevel;
	}

	/**
	 * This method returns the lower bound of the interval.
	 * @return a Matrix instance
	 */
	public Matrix getLowerLimit() {return lowerBound.getBoundValue();}
	
	/**
	 * This method returns the upper bound of the interval.
	 * @return a Matrix instance
	 */
	public Matrix getUpperLimit() {return upperBound.getBoundValue();}
	
	/**
	 * This method returns the probability level of the interval.
	 * @return a double
	 */
	public double getProbabilityLevel() {return probabilityLevel;}

	
	/**
	 * This method returns true if one of the bound of the confidence intervals contains a NaN
	 * @return a boolean
	 */
	public boolean isThereAnyNaN() {
		return getLowerLimit().doesContainAnyNaN() || getUpperLimit().doesContainAnyNaN();
	}
	
}
