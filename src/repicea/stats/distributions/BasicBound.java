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
package repicea.stats.distributions;

import java.io.Serializable;

import repicea.math.Matrix;

/**
 * The basic class for bounded distributions.
 * @author Mathieu Fortin - April 2016
 *
 */
@SuppressWarnings("serial")
public class BasicBound implements Serializable {

	private final boolean isUpperBound;
	private Matrix value;

	
	protected BasicBound(boolean isUpperBound) {
		this.isUpperBound = isUpperBound;
	}
	
	protected void setBoundValue(Matrix value) {
		this.value = value;
	}

	public Matrix getBoundValue() {return value;}
	
	protected boolean isUpperBound() {return isUpperBound;}
}
