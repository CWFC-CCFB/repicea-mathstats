/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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

import repicea.stats.StatisticalUtility;
import repicea.stats.distributions.GaussianErrorTermList.IndexableErrorTerm;

@SuppressWarnings("serial")
public class GaussianErrorTerm implements Serializable, Comparable<GaussianErrorTerm> {

	protected final int distanceIndex;
	protected Double value;
	protected final Double normalizedValue;
	
	public GaussianErrorTerm(IndexableErrorTerm caller) {
		this(caller, StatisticalUtility.getRandom().nextGaussian());
	}
	
	public GaussianErrorTerm(IndexableErrorTerm caller, double normalizedValue) {
		this.distanceIndex = ((IndexableErrorTerm) caller).getErrorTermIndex();
		this.normalizedValue = normalizedValue;
	}

	@Override
	public int compareTo(GaussianErrorTerm errorTerm) {
		if (this.distanceIndex < errorTerm.distanceIndex) {
			return -1;
		} else if (this.distanceIndex == errorTerm.distanceIndex) {
			return 0;
		} else {
			return 1;
		}
	}
	
	
	
}

