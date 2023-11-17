/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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
import java.util.ArrayList;
import java.util.List;

import repicea.math.AbstractMatrix;
import repicea.math.Matrix;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;

/**
 * The AbstractEmpirical class is an abstract class that supports the EmpiricalDistribution and 
 * ComplexEmpiricalDistribution classes.
 * 
 * @author Mathieu Fortin - November 2023
 *
 * @param <P> Matrix or ComplexMatrix
 */
@SuppressWarnings({ "rawtypes", "serial" })
public abstract class AbstractEmpiricalDistribution<P extends AbstractMatrix> implements Distribution<P>, Serializable {

	protected final List<P> observations;

	/**
	 * Constructor.
	 */
	public AbstractEmpiricalDistribution() {
		observations = new ArrayList<P>();
	}
	
	/**
	 * This method returns the number of observations in this nonparametric distribution.
	 * @return an integer
	 */
	public int getNumberOfRealizations() {return observations.size();}
	
	/**
	 * This method sets a given observation of the nonparametric distribution.
	 * @param value the value of the observation
	 */
	public void addRealization(P value) {observations.add(value);}
	
	/**
	 * This method returns the array that contains all the observations of this distribution.
	 * @return an array of Matrix instances
	 */
	public List<P> getRealizations() {return observations;}

	@Override
	public boolean isParametric() {
		return false;
	}

	@Override
	public final boolean isMultivariate() {
		if (observations != null && observations.size() > 0) {
			return observations.get(0) instanceof Matrix && observations.get(0).m_iRows > 1;
		} else {
			return false;
		}
	}

	@Override
	public final Type getType() {return Type.NONPARAMETRIC;}


//	@Override
//	public double getQuantile(double... values) {
//		if (observationsgetM)
//		// TODO to be implemented
//		return -1;
//	}

	@Override
	public P getRandomRealization() {
		int observationIndex = (int) (StatisticalUtility.getRandom().nextDouble() * getNumberOfRealizations());
		return getRealizations().get(observationIndex);
	}

}
