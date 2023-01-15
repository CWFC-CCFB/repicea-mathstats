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
package repicea.stats.estimates;

import repicea.math.Matrix;


/**
 * The CorrelationEstimate class represents a Pearson's correlation
 * estimate.
 * @author Mathieu Fortin - November 2012
 */
public class CorrelationEstimate extends SimpleEstimate {
	
	private static final long serialVersionUID = 20110912L;
	private int sampleSize;
	
	/**
	 * Constructor.
	 * @param mean the mean correlation
	 * @param sampleSize the sample size
	 */
	public CorrelationEstimate(double mean, int sampleSize) {
		super();
		Matrix meanMat = new Matrix(1,1);
		meanMat.setValueAt(0, 0, mean);
		setMean(meanMat);
		this.sampleSize = sampleSize;
	}
	
	/**
	 * This method returns the t value under the assumption that the true correlation is null.
	 * @return a double
	 */
	public double getStudentT() {
		double mean = getMean().getValueAt(0, 0);
		return mean * Math.sqrt((double) (sampleSize - 2)/(1 - mean * mean));
	}
	
	/**
	 * This method return the sample size.
	 * @return an integer
	 */
	public int getSampleSize() {return sampleSize;}

}
